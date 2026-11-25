#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>
//[[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

namespace {

bool encode_flat_membership_cpp(SEXP mb,
                                std::vector<int>& codes,
                                std::vector<int>& sizes,
                                int& n_clusters) {
    const R_xlen_t n_xlen = Rf_xlength(mb);
    if (n_xlen <= 0 || n_xlen > INT_MAX) {
        return false;
    }
    const int n = static_cast<int>(n_xlen);
    codes.resize(n);
    sizes.clear();

    const int type = TYPEOF(mb);
    if (type == INTSXP) {
        const int* ptr = INTEGER(mb);
        int min_v = ptr[0];
        int max_v = ptr[0];
        for (int i = 0; i < n; ++i) {
            const int v = ptr[i];
            if (v == NA_INTEGER) {
                return false;
            }
            if (v < min_v) min_v = v;
            if (v > max_v) max_v = v;
        }
        const int64_t span = static_cast<int64_t>(max_v) - static_cast<int64_t>(min_v) + 1;
        const int64_t dense_limit = std::min<int64_t>(std::max<int64_t>(2LL * n, 1024LL), 1048576LL);
        if (span > 0 && span <= dense_limit) {
            std::vector<int> remap(static_cast<size_t>(span), -1);
            int k = 0;
            for (int i = 0; i < n; ++i) {
                const size_t idx = static_cast<size_t>(static_cast<int64_t>(ptr[i]) - static_cast<int64_t>(min_v));
                int id = remap[idx];
                if (id < 0) {
                    id = k++;
                    remap[idx] = id;
                    sizes.push_back(1);
                } else {
                    sizes[id]++;
                }
                codes[i] = id;
            }
            n_clusters = k;
            return true;
        }
        std::unordered_map<int, int> remap;
        remap.reserve(std::min<int>(n, 1024));
        int k = 0;
        for (int i = 0; i < n; ++i) {
            const int v = ptr[i];
            auto it = remap.find(v);
            if (it == remap.end()) {
                const int id = k++;
                remap.emplace(v, id);
                sizes.push_back(1);
                codes[i] = id;
            } else {
                sizes[it->second]++;
                codes[i] = it->second;
            }
        }
        n_clusters = k;
        return true;
    }

    if (type == REALSXP) {
        const double* ptr = REAL(mb);
        bool small_int = true;
        int min_v = INT_MAX;
        int max_v = INT_MIN;
        const int64_t dense_limit = std::min<int64_t>(std::max<int64_t>(2LL * static_cast<int64_t>(n), 1024LL), 1048576LL);
        for (int i = 0; i < n; ++i) {
            const double v = ptr[i];
            if (ISNA(v) || ISNAN(v)) {
                return false;
            }
            if (small_int) {
                if (!R_FINITE(v) || v < -1048576.0 || v > 1048576.0 || v != std::trunc(v)) {
                    small_int = false;
                } else {
                    const int iv = static_cast<int>(v);
                    if (iv < min_v) min_v = iv;
                    if (iv > max_v) max_v = iv;
                    if (static_cast<int64_t>(max_v) - static_cast<int64_t>(min_v) + 1 > dense_limit) {
                        small_int = false;
                    }
                }
            }
        }
        if (small_int) {
            const size_t span = static_cast<size_t>(static_cast<int64_t>(max_v) - static_cast<int64_t>(min_v) + 1);
            std::vector<int> remap(span, -1);
            int k = 0;
            for (int i = 0; i < n; ++i) {
                const size_t idx = static_cast<size_t>(static_cast<int>(ptr[i]) - min_v);
                int id = remap[idx];
                if (id < 0) {
                    id = k++;
                    remap[idx] = id;
                    sizes.push_back(1);
                } else {
                    sizes[id]++;
                }
                codes[i] = id;
            }
            n_clusters = k;
            return true;
        }
        std::unordered_map<double, int> remap;
        remap.reserve(std::min<int>(n, 1024));
        int k = 0;
        for (int i = 0; i < n; ++i) {
            double v = ptr[i];
            if (v == 0.0) v = 0.0; // normalize -0.0 to +0.0 to match R's ==
            auto it = remap.find(v);
            if (it == remap.end()) {
                const int id = k++;
                remap.emplace(v, id);
                sizes.push_back(1);
                codes[i] = id;
            } else {
                sizes[it->second]++;
                codes[i] = it->second;
            }
        }
        n_clusters = k;
        return true;
    }

    if (type == STRSXP) {
        std::unordered_map<std::string, int> remap;
        remap.reserve(std::min<int>(n, 1024));
        int k = 0;
        const void* vmax = vmaxget();
        for (int i = 0; i < n; ++i) {
            SEXP s = STRING_ELT(mb, i);
            if (s == NA_STRING || Rf_getCharCE(s) == CE_BYTES) {
                vmaxset(vmax);
                return false;
            }
            std::string key(Rf_translateCharUTF8(s));
            vmaxset(vmax);
            auto it = remap.find(key);
            if (it == remap.end()) {
                const int id = k++;
                remap.emplace(std::move(key), id);
                sizes.push_back(1);
                codes[i] = id;
            } else {
                sizes[it->second]++;
                codes[i] = it->second;
            }
        }
        n_clusters = k;
        return true;
    }

    return false;
}

inline double pair_cluster_ecs(int intersect_size, int size1, int size2, bool use_r_l1_formula) {
    if (!use_r_l1_formula) {
        return static_cast<double>(intersect_size) / std::max(size1, size2);
    }
    const double csize1 = static_cast<double>(size1);
    const double csize2 = static_cast<double>(size2);
    const double inter = static_cast<double>(intersect_size);
    double temp_ecs = (1.0 / csize1 - 1.0 / csize2) * inter;
    if (temp_ecs < 0.0) {
        temp_ecs = -temp_ecs;
    }
    return 1.0 - 0.5 * (temp_ecs + (csize1 - inter) / csize1 + (csize2 - inter) / csize2);
}

// Compute pair ECS (and optionally weighted accumulation and R-exact compensated mean)
// for two pre-encoded partitions. Uses a dense contingency table when K1*K2 is bounded,
// and an O(n) hash map when K1*K2 is large (e.g. singleton-heavy partitions).
double score_encoded_pair(const int* c1, const int* s1, int K1,
                          const int* c2, const int* s2, int K2,
                          int n,
                          bool use_r_l1_formula,
                          bool accumulate_ecc,
                          double wi, double wj, double norm_factor,
                          double* cons,
                          bool compute_mean,
                          double* out_elscores,
                          std::vector<int>& cont,
                          std::vector<double>& unique_ecs,
                          std::vector<double>& unique_term,
                          std::vector<int>& elem_pair_ids) {
    const size_t table_sz = static_cast<size_t>(K1) * static_cast<size_t>(K2);
    const size_t dense_cap = std::max<size_t>(4 * static_cast<size_t>(n), 65536);

    if (table_sz <= dense_cap) {
        cont.assign(table_sz, 0);
        for (int k = 0; k < n; ++k) {
            cont[static_cast<size_t>(c1[k]) * K2 + c2[k]]++;
        }

        unique_ecs.resize(table_sz);
        if (accumulate_ecc) {
            unique_term.resize(table_sz);
        }
        for (int a = 0; a < K1; ++a) {
            const int sa = s1[a];
            const size_t row_off = static_cast<size_t>(a) * K2;
            for (int b = 0; b < K2; ++b) {
                const int cnt = cont[row_off + b];
                if (cnt > 0) {
                    const double val = pair_cluster_ecs(cnt, sa, s2[b], use_r_l1_formula);
                    unique_ecs[row_off + b] = val;
                    if (accumulate_ecc) {
                        // Exact left-to-right R precedence: current_ecs * weights[i] / norm_factor * weights[j]
                        unique_term[row_off + b] = ((val * wi) / norm_factor) * wj;
                    }
                }
            }
        }

        if (out_elscores != nullptr) {
            for (int k = 0; k < n; ++k) {
                out_elscores[k] = unique_ecs[static_cast<size_t>(c1[k]) * K2 + c2[k]];
            }
        }
        if (compute_mean) {
            // Exact two-pass long-double compensated mean matching R's mean.default()
            long double s = 0.0L;
            if (accumulate_ecc) {
                for (int k = 0; k < n; ++k) {
                    const size_t idx = static_cast<size_t>(c1[k]) * K2 + c2[k];
                    cons[k] += unique_term[idx];
                    s += unique_ecs[idx];
                }
            } else {
                for (int k = 0; k < n; ++k) {
                    s += unique_ecs[static_cast<size_t>(c1[k]) * K2 + c2[k]];
                }
            }
            s /= n;
            if (R_FINITE(static_cast<double>(s))) {
                long double t = 0.0L;
                for (int k = 0; k < n; ++k) {
                    t += (unique_ecs[static_cast<size_t>(c1[k]) * K2 + c2[k]] - s);
                }
                s += t / n;
            }
            return static_cast<double>(s);
        }
        if (accumulate_ecc) {
            for (int k = 0; k < n; ++k) {
                cons[k] += unique_term[static_cast<size_t>(c1[k]) * K2 + c2[k]];
            }
        }
        return 0.0;
    }

    // Sparse O(n) path for high-K pairs
    elem_pair_ids.resize(n);
    std::unordered_map<uint64_t, int> pair_map;
    pair_map.reserve(std::min<int>(n, 4096));
    std::vector<int> pair_counts;
    std::vector<int> pair_a;
    std::vector<int> pair_b;

    int n_pairs = 0;
    for (int k = 0; k < n; ++k) {
        const uint64_t key = (static_cast<uint64_t>(static_cast<uint32_t>(c1[k])) << 32) |
                             static_cast<uint64_t>(static_cast<uint32_t>(c2[k]));
        auto it = pair_map.find(key);
        if (it == pair_map.end()) {
            const int pid = n_pairs++;
            pair_map.emplace(key, pid);
            pair_counts.push_back(1);
            pair_a.push_back(c1[k]);
            pair_b.push_back(c2[k]);
            elem_pair_ids[k] = pid;
        } else {
            pair_counts[it->second]++;
            elem_pair_ids[k] = it->second;
        }
    }

    unique_ecs.resize(n_pairs);
    if (accumulate_ecc) {
        unique_term.resize(n_pairs);
    }
    for (int p = 0; p < n_pairs; ++p) {
        const double val = pair_cluster_ecs(pair_counts[p], s1[pair_a[p]], s2[pair_b[p]], use_r_l1_formula);
        unique_ecs[p] = val;
        if (accumulate_ecc) {
            unique_term[p] = ((val * wi) / norm_factor) * wj;
        }
    }

    if (out_elscores != nullptr) {
        for (int k = 0; k < n; ++k) {
            out_elscores[k] = unique_ecs[elem_pair_ids[k]];
        }
    }
    if (compute_mean) {
        long double s = 0.0L;
        if (accumulate_ecc) {
            for (int k = 0; k < n; ++k) {
                const int pid = elem_pair_ids[k];
                cons[k] += unique_term[pid];
                s += unique_ecs[pid];
            }
        } else {
            for (int k = 0; k < n; ++k) {
                s += unique_ecs[elem_pair_ids[k]];
            }
        }
        s /= n;
        if (R_FINITE(static_cast<double>(s))) {
            long double t = 0.0L;
            for (int k = 0; k < n; ++k) {
                t += (unique_ecs[elem_pair_ids[k]] - s);
            }
            s += t / n;
        }
        return static_cast<double>(s);
    }
    if (accumulate_ecc) {
        for (int k = 0; k < n; ++k) {
            cons[k] += unique_term[elem_pair_ids[k]];
        }
    }
    return 0.0;
}

SEXP weighted_ecc_flat_impl(List mb_list,
                            NumericVector weights,
                            bool calculate_sim_matrix,
                            bool use_r_l1_formula) {
    const int U = mb_list.size();
    if (U < 2) {
        stop("At least two clusterings are required to calculate the consistency.");
    }
    const R_xlen_t n_xlen = Rf_xlength(mb_list[0]);
    if (n_xlen <= 0 || n_xlen > INT_MAX) {
        return R_NilValue;
    }
    const int n = static_cast<int>(n_xlen);

    // Encode each partition once in O(U * n)
    std::vector<std::vector<int>> codes(U);
    std::vector<std::vector<int>> sizes(U);
    std::vector<int> K(U, 0);
    for (int i = 0; i < U; ++i) {
        SEXP mb = mb_list[i];
        if (Rf_xlength(mb) != n_xlen) {
            return R_NilValue;
        }
        if (!encode_flat_membership_cpp(mb, codes[i], sizes[i], K[i])) {
            return R_NilValue;
        }
    }

    // Match R's sum(weights) long-double accumulation and norm_factor
    long double total_w_ld = 0.0L;
    for (int i = 0; i < U; ++i) {
        total_w_ld += weights[i];
    }
    const double total_weights = static_cast<double>(total_w_ld);
    const double norm_factor = (total_weights * (total_weights - 1.0)) / 2.0;

    NumericVector consistency(n);
    double* cons = REAL(consistency);
    NumericMatrix sim_mat;
    if (calculate_sim_matrix) {
        sim_mat = NumericMatrix(U, U);
    }

    std::vector<int> cont;
    std::vector<double> unique_ecs;
    std::vector<double> unique_term;
    std::vector<int> elem_pair_ids;

    for (int i = 0; i < U - 1; ++i) {
        checkUserInterrupt();
        const double wi = weights[i];
        // Exact R left-to-right precedence: weights[i] / norm_factor * (weights[i] - 1) / 2
        const double diag_i = ((wi / norm_factor) * (wi - 1.0)) / 2.0;
        for (int k = 0; k < n; ++k) {
            cons[k] += diag_i;
        }

        const int* c1 = codes[i].data();
        const int* s1 = sizes[i].data();
        const int K1 = K[i];

        for (int j = i + 1; j < U; ++j) {
            const double wj = weights[j];
            const double pair_mean = score_encoded_pair(
                c1, s1, K1,
                codes[j].data(), sizes[j].data(), K[j],
                n, use_r_l1_formula,
                true, wi, wj, norm_factor, cons,
                calculate_sim_matrix, nullptr,
                cont, unique_ecs, unique_term, elem_pair_ids
            );
            if (calculate_sim_matrix) {
                sim_mat(i, j) = pair_mean;
                sim_mat(j, i) = pair_mean;
            }
        }

        if (calculate_sim_matrix) {
            sim_mat(i, i) = 1.0;
        }
    }

    const double w_last = weights[U - 1];
    const double diag_last = ((w_last / norm_factor) * (w_last - 1.0)) / 2.0;
    for (int k = 0; k < n; ++k) {
        cons[k] += diag_last;
    }

    if (!calculate_sim_matrix) {
        return consistency;
    }

    sim_mat(U - 1, U - 1) = 1.0;
    return List::create(Named("ecc") = consistency, Named("ecs_matrix") = sim_mat);
}

SEXP sim_matrix_flat_impl(List mb_list, bool use_r_l1_formula) {
    const int U = mb_list.size();
    if (U == 1) {
        NumericMatrix m(1, 1);
        m(0, 0) = 1.0;
        return m;
    }
    const R_xlen_t n_xlen = Rf_xlength(mb_list[0]);
    if (n_xlen <= 0 || n_xlen > INT_MAX) {
        return R_NilValue;
    }
    const int n = static_cast<int>(n_xlen);

    std::vector<std::vector<int>> codes(U);
    std::vector<std::vector<int>> sizes(U);
    std::vector<int> K(U, 0);
    for (int i = 0; i < U; ++i) {
        SEXP mb = mb_list[i];
        if (Rf_xlength(mb) != n_xlen) {
            return R_NilValue;
        }
        if (!encode_flat_membership_cpp(mb, codes[i], sizes[i], K[i])) {
            return R_NilValue;
        }
    }

    NumericMatrix sim_mat(U, U);
    std::vector<int> cont;
    std::vector<double> unique_ecs;
    std::vector<double> unique_term;
    std::vector<int> elem_pair_ids;

    for (int i = 0; i < U - 1; ++i) {
        checkUserInterrupt();
        sim_mat(i, i) = 1.0;
        const int* c1 = codes[i].data();
        const int* s1 = sizes[i].data();
        const int K1 = K[i];
        for (int j = i + 1; j < U; ++j) {
            const double pair_mean = score_encoded_pair(
                c1, s1, K1,
                codes[j].data(), sizes[j].data(), K[j],
                n, use_r_l1_formula,
                false, 0.0, 0.0, 1.0, nullptr,
                true, nullptr,
                cont, unique_ecs, unique_term, elem_pair_ids
            );
            sim_mat(i, j) = pair_mean;
            sim_mat(j, i) = pair_mean;
        }
    }
    sim_mat(U - 1, U - 1) = 1.0;
    return sim_mat;
}

} // namespace

// [[Rcpp::export]]
SEXP weighted_ecc_flat_cpp(List mb_list, NumericVector weights, bool calculate_sim_matrix) {
    return weighted_ecc_flat_impl(mb_list, weights, calculate_sim_matrix, false);
}

// [[Rcpp::export(rng = false)]]
SEXP weighted_ecc_flat_char_cpp(List mb_list, NumericVector weights, bool calculate_sim_matrix) {
    return weighted_ecc_flat_impl(mb_list, weights, calculate_sim_matrix, true);
}

// [[Rcpp::export]]
SEXP sim_matrix_flat_cpp(List mb_list) {
    return sim_matrix_flat_impl(mb_list, false);
}

// [[Rcpp::export(rng = false)]]
SEXP sim_matrix_flat_char_cpp(List mb_list) {
    return sim_matrix_flat_impl(mb_list, true);
}

// [[Rcpp::export(rng = false)]]
SEXP corrected_l1_mb_cpp(SEXP mb1, SEXP mb2) {
    const R_xlen_t n_xlen = Rf_xlength(mb1);
    if (n_xlen != Rf_xlength(mb2)) {
        return R_NilValue;
    }
    if (n_xlen == 0) {
        return NumericVector(0);
    }
    if (n_xlen > INT_MAX) {
        return R_NilValue;
    }
    const int n = static_cast<int>(n_xlen);
    std::vector<int> c1, s1, c2, s2;
    int K1 = 0, K2 = 0;
    if (!encode_flat_membership_cpp(mb1, c1, s1, K1) ||
        !encode_flat_membership_cpp(mb2, c2, s2, K2)) {
        return R_NilValue;
    }
    NumericVector out(n);
    std::vector<int> cont;
    std::vector<double> unique_ecs;
    std::vector<double> unique_term;
    std::vector<int> elem_pair_ids;
    score_encoded_pair(
        c1.data(), s1.data(), K1,
        c2.data(), s2.data(), K2,
        n, true,
        false, 0.0, 0.0, 1.0, nullptr,
        false, REAL(out),
        cont, unique_ecs, unique_term, elem_pair_ids
    );
    return out;
}

// [[Rcpp::export]]
IntegerMatrix myContTable(IntegerVector &a, IntegerVector &b, int minim_mb1, int minim_mb2) {
	int maxim_mb1 = max(a), maxim_mb2 = max(b);
    int nClustersA = maxim_mb1 - minim_mb1 + 1;
	int nClustersB = maxim_mb2 - minim_mb2 + 1, n = a.size();

    IntegerMatrix result(nClustersA, nClustersB);

    for(int i = 0; i < n; i++) {
        result(a(i) - minim_mb1, b(i) - minim_mb2)++;
    }

	return(result);
}

// [[Rcpp::export]]
NumericVector disjointECS(IntegerVector mb1, IntegerVector mb2) {
	const int n = mb1.size();
	if (n == 0) {
		return NumericVector(0);
	}
	int minim_mb1 = mb1[0], maxim_mb1 = mb1[0];
	int minim_mb2 = mb2[0], maxim_mb2 = mb2[0];
	for (int i = 1; i < n; ++i) {
		const int v1 = mb1[i], v2 = mb2[i];
		if (v1 < minim_mb1) minim_mb1 = v1;
		if (v1 > maxim_mb1) maxim_mb1 = v1;
		if (v2 < minim_mb2) minim_mb2 = v2;
		if (v2 > maxim_mb2) maxim_mb2 = v2;
	}
	const int64_t span1 = static_cast<int64_t>(maxim_mb1) - static_cast<int64_t>(minim_mb1) + 1;
	const int64_t span2 = static_cast<int64_t>(maxim_mb2) - static_cast<int64_t>(minim_mb2) + 1;
	const int64_t dense_cap = std::max<int64_t>(4LL * static_cast<int64_t>(n), 65536LL);

	if (span1 > 0 && span2 > 0 &&
	    static_cast<double>(span1) * static_cast<double>(span2) <= static_cast<double>(dense_cap)) {
		int cluster1, cluster2;
		IntegerMatrix contTable = myContTable(mb1, mb2, minim_mb1, minim_mb2);
		int nClusters1 = contTable.nrow(), nClusters2 = contTable.ncol();

		NumericMatrix uniqueECSvals(nClusters1, nClusters2);
		NumericVector ecsScore(n);
		IntegerVector clustersSizes1 = rowSums(contTable), clustersSizes2 = colSums(contTable);

		for (int i = 0; i < nClusters1; i++) {
			for (int j = nClusters2 - 1; j >= 0; j--) {
				uniqueECSvals(i, j) = (double) contTable(i, j) / std::max(clustersSizes1(i), clustersSizes2(j));
			}
		}

		for (int i = 0; i < n; i++) {
			cluster1 = mb1(i) - minim_mb1;
			cluster2 = mb2(i) - minim_mb2;
			ecsScore(i) = uniqueECSvals(cluster1, cluster2);
		}

		return ecsScore;
	}

	// Sparse O(n) path for high-cardinality pairs (avoids O(K1 * K2) allocation).
	std::vector<int> c1, s1, c2, s2;
	int K1 = 0, K2 = 0;
	encode_flat_membership_cpp(mb1, c1, s1, K1);
	encode_flat_membership_cpp(mb2, c2, s2, K2);
	NumericVector out(n);
	std::vector<int> cont;
	std::vector<double> unique_ecs;
	std::vector<double> unique_term;
	std::vector<int> elem_pair_ids;
	score_encoded_pair(
		c1.data(), s1.data(), K1,
		c2.data(), s2.data(), K2,
		n, false,
		false, 0.0, 0.0, 1.0, nullptr,
		false, REAL(out),
		cont, unique_ecs, unique_term, elem_pair_ids
	);
	return out;
}

// [[Rcpp::export]]
double disjointECSaverage(IntegerVector mb1, IntegerVector mb2) {
	// NOTE Eigen doesn't seem to improve performance
	const int n = mb1.size();
	if (n == 0) {
		return 0.0;
	}
	int minim_mb1 = mb1[0], maxim_mb1 = mb1[0];
	int minim_mb2 = mb2[0], maxim_mb2 = mb2[0];
	for (int i = 1; i < n; ++i) {
		const int v1 = mb1[i], v2 = mb2[i];
		if (v1 < minim_mb1) minim_mb1 = v1;
		if (v1 > maxim_mb1) maxim_mb1 = v1;
		if (v2 < minim_mb2) minim_mb2 = v2;
		if (v2 > maxim_mb2) maxim_mb2 = v2;
	}
	const int64_t span1 = static_cast<int64_t>(maxim_mb1) - static_cast<int64_t>(minim_mb1) + 1;
	const int64_t span2 = static_cast<int64_t>(maxim_mb2) - static_cast<int64_t>(minim_mb2) + 1;
	const int64_t dense_cap = std::max<int64_t>(4LL * static_cast<int64_t>(n), 65536LL);

	if (span1 > 0 && span2 > 0 &&
	    static_cast<double>(span1) * static_cast<double>(span2) <= static_cast<double>(dense_cap)) {
		IntegerMatrix contTable = myContTable(mb1, mb2, minim_mb1, minim_mb2);
		int nClusters1 = contTable.nrow(), nClusters2 = contTable.ncol();

		double avgECS = 0, intersectSize;
		IntegerVector clustersSizes1 = rowSums(contTable), clustersSizes2 = colSums(contTable);

		for (int i = 0; i < nClusters1; i++) {
			for (int j = nClusters2 - 1; j >= 0; j--) {
				intersectSize = contTable(i, j);
				int maxClusterSize = std::max(clustersSizes1(i), clustersSizes2(j));
				if (maxClusterSize > 0) {
					// Promote before multiplying: large cell counts exceed integer range.
					avgECS += intersectSize * intersectSize /
						(static_cast<double>(n) * maxClusterSize);
				}
			}
		}

		return avgECS;
	}

	// Sparse O(n log n) path for high-cardinality pairs: sort observed pairs by
	// (label1 ASC, label2 DESC) to match the dense kernel's exact summation order.
	std::unordered_map<int, int> size1;
	std::unordered_map<int, int> size2;
	std::unordered_map<uint64_t, int> pair_counts;
	size1.reserve(std::min<int>(n, 4096));
	size2.reserve(std::min<int>(n, 4096));
	pair_counts.reserve(std::min<int>(n, 4096));

	for (int i = 0; i < n; ++i) {
		const int a = mb1[i];
		const int b = mb2[i];
		size1[a]++;
		size2[b]++;
		const uint64_t key = (static_cast<uint64_t>(static_cast<uint32_t>(a)) << 32) |
		                     static_cast<uint64_t>(static_cast<uint32_t>(b));
		pair_counts[key]++;
	}

	struct SparsePair {
		int a;
		int b;
		int cnt;
	};
	std::vector<SparsePair> pairs;
	pairs.reserve(pair_counts.size());
	for (const auto& kv : pair_counts) {
		const int a = static_cast<int>(static_cast<uint32_t>(kv.first >> 32));
		const int b = static_cast<int>(static_cast<uint32_t>(kv.first & 0xFFFFFFFFu));
		pairs.push_back(SparsePair{a, b, kv.second});
	}
	std::sort(pairs.begin(), pairs.end(), [](const SparsePair& x, const SparsePair& y) {
		if (x.a != y.a) return x.a < y.a;
		return x.b > y.b;
	});

	double avgECS = 0.0;
	const double n_d = static_cast<double>(n);
	for (const SparsePair& sp : pairs) {
		const int maxClusterSize = std::max(size1[sp.a], size2[sp.b]);
		const double intersectSize = static_cast<double>(sp.cnt);
		avgECS += intersectSize * intersectSize / (n_d * maxClusterSize);
	}
	return avgECS;
}
