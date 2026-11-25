#include <Rcpp.h>
#include <RcppEigen.h>
#include <algorithm>
#include <functional>
#include <vector>
//[[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

namespace {

struct DisjointSet {
    std::vector<int> parent;
    std::vector<int> rank;
    int components;

    explicit DisjointSet(int n = 0) {
        reset(n);
    }

    void reset(int n) {
        parent.resize(n);
        rank.assign(n, 0);
        for (int i = 0; i < n; ++i) {
            parent[i] = i;
        }
        components = n;
    }

    inline int find(int x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    }

    inline void unite(int a, int b) {
        int ra = find(a);
        int rb = find(b);
        if (ra == rb) {
            return;
        }
        if (rank[ra] < rank[rb]) {
            parent[ra] = rb;
        } else if (rank[ra] > rank[rb]) {
            parent[rb] = ra;
        } else {
            parent[rb] = ra;
            rank[ra]++;
        }
        components--;
    }
};

struct PruneEdge {
    int u;
    int v;
    int c;
};

int run_prune_binary_search(int n_neigh,
                            int target_n_conn_comps,
                            const std::function<int(int)>& count_components_for_idx) {
    // Mirrors R's 1-based binary search over possible_values = (0:n_neigh) / (2*n_neigh - (0:n_neigh)).
    int start_n = 1;
    int stop_n = n_neigh + 1;
    int middle = 1;
    while (start_n <= stop_n) {
        middle = (stop_n + start_n) / 2;
        const int current_n_conn_comps = count_components_for_idx(middle - 1);
        if (current_n_conn_comps > target_n_conn_comps) {
            stop_n = middle;
            if (start_n == stop_n) {
                start_n = middle - 1;
                if (start_n < 1) {
                    stop("No valid pruning threshold preserves the graph's connected components.");
                }
            }
        } else {
            start_n = middle + 1;
            if (start_n == stop_n) {
                break;
            }
        }
    }
    return middle - 1;
}

} // namespace

// [[Rcpp::export(rng = false)]]
Eigen::SparseMatrix<double> pruneSNN(Eigen::SparseMatrix<double> snnMatrix, double prune = 0) {
    if (prune <= 0) {
        return snnMatrix;
    }

    for (int i = 0; i < snnMatrix.outerSize(); i++) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(snnMatrix, i); it; ++it) {
            if (it.value() <= prune) {
                it.valueRef() = 0;
            }
        }
    }

    snnMatrix.prune(0.0); // actually remove pruned values

    return snnMatrix;
}

// TODO check if k is lower than the max possible value
// [[Rcpp::export(rng = false)]]
Eigen::SparseMatrix<double> computeSNN(Eigen::SparseMatrix<double> &nnMatrix, int k, double prune = 0) {
    Eigen::SparseMatrix<double> SNN = nnMatrix * (nnMatrix.transpose());

    for (int i = 0; i < SNN.outerSize(); i++) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(SNN, i); it; ++it) {
            it.valueRef() = it.value() / (k + (k - it.value()));
            if (prune <= 0) {
                continue;
            }

            if (it.value() <= prune) {
                it.valueRef() = 0;
            }
        }
    }

    if (prune > 0) {
        SNN.prune(0.0); // actually remove pruned values
    }

    return SNN;
}

//' Computes the NN adjacency matrix given the neighbours
//'
//' @param nnRanked A matrix with the lists of the nearest neighbours for each point
//' @param k The number of neighbours to consider. Defaults to `-1`, which means all neighbours.
//' @param start The index of the first neighbour to consider. Defaults to `0`.
//' @param prune The threshold to prune the SNN matrix. If -1, the function will only return the NN matrix. Defaults to `0`.
//'
//' @export
//'
//' @return A list with the NN and SNN adjacency matrices.
// [[Rcpp::export(rng = false)]]
List getNNmatrix(Eigen::MatrixXd nnRanked, int k = -1, int start = 0, double prune = 0) {
    int nRows = nnRanked.rows(), nCols = nnRanked.cols();

    if (k == -1 || k > nCols) {
        k = nCols;
    }

    if (start + k > nCols) {
        k = nCols - start;
    }

    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(static_cast<size_t>(nRows) * static_cast<size_t>(k));

    for (int j = start; j < start + k; j++) {
        for (int i = 0; i < nRows; i++) {
            tripletList.push_back(Eigen::Triplet<double>(i, nnRanked(i, j) - 1, 1));
        }
    }

    Eigen::SparseMatrix<double> NN(nRows, nRows);
    NN.setFromTriplets(tripletList.begin(), tripletList.end());

    if (prune < 0) {
        return List::create(Named("nn") = NN);
    }

    return List::create(Named("nn") = NN, Named("snn") = computeSNN(NN, k, prune));
}

// Compute the highest pruning parameter for the SNN graph and its pruned
// adjacency matrix directly on CSC sparse storage, avoiding igraph conversion
// and unpruned SNN materialization for standard binary NN matrices.
// [[Rcpp::export(rng = false)]]
List highest_prune_param_cpp(const Eigen::MappedSparseMatrix<double> nn_matrix,
                             int n_neigh,
                             bool return_matrix = true,
                             bool use_nn_target_comps = false) {
    const int n = nn_matrix.rows();
    if (n != nn_matrix.cols()) {
        stop("nn_matrix must be square.");
    }
    if (n_neigh < 1 || n_neigh > INT_MAX / 2) {
        stop("n_neigh must be an integer between 1 and INT_MAX / 2.");
    }

    const int* c_outer = nn_matrix.outerIndexPtr();
    const int* c_inner = nn_matrix.innerIndexPtr();
    const double* c_val = nn_matrix.valuePtr();
    const int nnz = nn_matrix.nonZeros();

    // Check whether nn_matrix is binary with row degrees <= n_neigh.
    bool is_standard_binary = true;
    std::vector<int> row_deg(n, 0);
    for (int p = 0; p < nnz; ++p) {
        if (c_val[p] != 1.0) {
            is_standard_binary = false;
            break;
        }
        const int r = c_inner[p];
        row_deg[r]++;
        if (row_deg[r] > n_neigh) {
            is_standard_binary = false;
            break;
        }
    }

    if (is_standard_binary && n > 0) {
        // Build CSR view of nn_matrix in O(nnz).
        std::vector<int> row_ptr(n + 1, 0);
        for (int i = 0; i < n; ++i) {
            row_ptr[i + 1] = row_ptr[i] + row_deg[i];
        }
        std::vector<int> row_col(nnz);
        std::vector<int> write_pos = row_ptr;
        DisjointSet uf_nn(use_nn_target_comps ? n : 0);
        for (int col = 0; col < n; ++col) {
            for (int p = c_outer[col]; p < c_outer[col + 1]; ++p) {
                const int r = c_inner[p];
                row_col[write_pos[r]++] = col;
                if (use_nn_target_comps) {
                    uf_nn.unite(r, col);
                }
            }
        }

        // Single upper-triangle pass computing SNN overlaps c = (NN * NN^T)_{i, v}.
        // Edges with c == 1 only contribute to target_n_conn_comps (and to the
        // returned matrix if prune_val == 0), so only c >= 2 edges are stored.
        DisjointSet uf_snn(n);
        std::vector<int> diag_c(n, 0);
        std::vector<PruneEdge> edges_ge2;
        edges_ge2.reserve(static_cast<size_t>(nnz) * 2);

        std::vector<int> counts(n, 0);
        std::vector<int> touched;
        touched.reserve(std::min<int>(n, 1024));

        for (int i = 0; i < n; ++i) {
            if ((i & 1023) == 0) {
                checkUserInterrupt();
            }
            touched.clear();
            const int r_start = row_ptr[i];
            const int r_end = row_ptr[i + 1];
            for (int rp = r_start; rp < r_end; ++rp) {
                const int u = row_col[rp];
                const int c_start = c_outer[u];
                // CSC inner indices are sorted ascending: scan backwards until v < i.
                for (int cp = c_outer[u + 1] - 1; cp >= c_start; --cp) {
                    const int v = c_inner[cp];
                    if (v < i) {
                        break;
                    }
                    if (counts[v] == 0) {
                        touched.push_back(v);
                    }
                    counts[v]++;
                }
            }
            for (int v : touched) {
                const int c = counts[v];
                counts[v] = 0;
                if (v == i) {
                    diag_c[i] = c;
                } else {
                    uf_snn.unite(i, v);
                    if (c >= 2) {
                        edges_ge2.push_back(PruneEdge{i, v, c});
                    }
                }
            }
        }

        const int target_n_conn_comps = use_nn_target_comps ? uf_nn.components : uf_snn.components;
        DisjointSet uf_cur(n);
        const int best_c = run_prune_binary_search(
            n_neigh,
            target_n_conn_comps,
            [&](int c_thresh) -> int {
                if (c_thresh == 0) {
                    return uf_snn.components;
                }
                uf_cur.reset(n);
                for (const PruneEdge& e : edges_ge2) {
                    if (e.c > c_thresh) {
                        uf_cur.unite(e.u, e.v);
                    }
                }
                return uf_cur.components;
            }
        );

        const double prune_val = static_cast<double>(best_c) / (2.0 * static_cast<double>(n_neigh) - static_cast<double>(best_c));
        if (!return_matrix) {
            return List::create(Named("prune_value") = prune_val);
        }

        if (best_c >= 1) {
            // Every surviving off-diagonal edge (c > best_c >= 1) has c >= 2 and
            // is already in edges_ge2. Build exact symmetric CSC directly.
            std::vector<int> col_nnz(n, 0);
            for (int i = 0; i < n; ++i) {
                if (diag_c[i] > best_c) {
                    col_nnz[i]++;
                }
            }
            for (const PruneEdge& e : edges_ge2) {
                if (e.c > best_c) {
                    col_nnz[e.u]++;
                    col_nnz[e.v]++;
                }
            }

            size_t total_nnz = 0;
            for (int i = 0; i < n; ++i) {
                total_nnz += static_cast<size_t>(col_nnz[i]);
            }

            Eigen::SparseMatrix<double> SNN(n, n);
            SNN.resize(n, n);
            SNN.resizeNonZeros(static_cast<Eigen::Index>(total_nnz));
            int* out_p = SNN.outerIndexPtr();
            int* out_i = SNN.innerIndexPtr();
            double* out_x = SNN.valuePtr();

            out_p[0] = 0;
            for (int i = 0; i < n; ++i) {
                out_p[i + 1] = out_p[i] + col_nnz[i];
            }
            std::vector<int> col_pos(out_p, out_p + n);

            // Precompute exact Jaccard weight table for c = 0..n_neigh matching computeSNN:
            // c / (k + (k - c)).
            std::vector<double> w_table(n_neigh + 1, 0.0);
            for (int c = 1; c <= n_neigh; ++c) {
                w_table[c] = static_cast<double>(c) / (static_cast<double>(n_neigh) + static_cast<double>(n_neigh - c));
            }

            for (int i = 0; i < n; ++i) {
                if (diag_c[i] > best_c) {
                    const int pos = col_pos[i]++;
                    out_i[pos] = i;
                    out_x[pos] = w_table[diag_c[i]];
                }
            }
            for (const PruneEdge& e : edges_ge2) {
                if (e.c > best_c) {
                    const double w = w_table[e.c];
                    const int pos_u = col_pos[e.u]++;
                    out_i[pos_u] = e.v;
                    out_x[pos_u] = w;
                    const int pos_v = col_pos[e.v]++;
                    out_i[pos_v] = e.u;
                    out_x[pos_v] = w;
                }
            }

            // Sort each column's row indices (typically ~20-60 entries) for strict CSC invariant.
            std::vector<std::pair<int, double>> col_buf;
            for (int col = 0; col < n; ++col) {
                const int start = out_p[col];
                const int end = out_p[col + 1];
                bool sorted = true;
                for (int p = start + 1; p < end; ++p) {
                    if (out_i[p] < out_i[p - 1]) {
                        sorted = false;
                        break;
                    }
                }
                if (!sorted) {
                    const int len = end - start;
                    col_buf.resize(len);
                    for (int k = 0; k < len; ++k) {
                        col_buf[k] = std::make_pair(out_i[start + k], out_x[start + k]);
                    }
                    std::sort(col_buf.begin(), col_buf.end(),
                              [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
                                  return a.first < b.first;
                              });
                    for (int k = 0; k < len; ++k) {
                        out_i[start + k] = col_buf[k].first;
                        out_x[start + k] = col_buf[k].second;
                    }
                }
            }

            return List::create(
                Named("prune_value") = prune_val,
                Named("adj_matrix") = SNN
            );
        }

        // best_c == 0 (prune_val == 0): return full unpruned SNN matrix.
        Eigen::SparseMatrix<double> nn_copy(nn_matrix);
        return List::create(
            Named("prune_value") = prune_val,
            Named("adj_matrix") = computeSNN(nn_copy, n_neigh, 0.0)
        );
    }

    // General fallback for non-binary or high-degree nn_matrix inputs.
    Eigen::SparseMatrix<double> nn_copy(nn_matrix);
    Eigen::SparseMatrix<double> SNN = computeSNN(nn_copy, n_neigh, 0.0);

    DisjointSet uf_target(n);
    if (use_nn_target_comps) {
        for (int col = 0; col < n; ++col) {
            for (int p = c_outer[col]; p < c_outer[col + 1]; ++p) {
                if (c_val[p] != 0.0) {
                    uf_target.unite(c_inner[p], col);
                }
            }
        }
    } else {
        for (int col = 0; col < SNN.outerSize(); ++col) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(SNN, col); it; ++it) {
                if (it.value() != 0.0) {
                    uf_target.unite(it.row(), col);
                }
            }
        }
    }
    const int target_n_conn_comps = uf_target.components;

    DisjointSet uf_cur(n);
    const int best_c = run_prune_binary_search(
        n_neigh,
        target_n_conn_comps,
        [&](int c_thresh) -> int {
            const double thr = static_cast<double>(c_thresh) / (2.0 * static_cast<double>(n_neigh) - static_cast<double>(c_thresh));
            uf_cur.reset(n);
            for (int col = 0; col < SNN.outerSize(); ++col) {
                for (Eigen::SparseMatrix<double>::InnerIterator it(SNN, col); it; ++it) {
                    if (it.value() > thr) {
                        uf_cur.unite(it.row(), col);
                    }
                }
            }
            return uf_cur.components;
        }
    );

    const double prune_val = static_cast<double>(best_c) / (2.0 * static_cast<double>(n_neigh) - static_cast<double>(best_c));
    if (!return_matrix) {
        return List::create(Named("prune_value") = prune_val);
    }
    if (prune_val > 0.0) {
        for (int col = 0; col < SNN.outerSize(); ++col) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(SNN, col); it; ++it) {
                if (it.value() <= prune_val) {
                    it.valueRef() = 0.0;
                }
            }
        }
        SNN.prune(0.0);
    }
    return List::create(
        Named("prune_value") = prune_val,
        Named("adj_matrix") = SNN
    );
}

// Compute weak connected-component counts incrementally from a 1-based kNN
// index matrix across a sorted sequence of neighbour counts.
// [[Rcpp::export(rng = false)]]
NumericVector nn_conn_comps_from_idx_cpp(SEXP nn_idx, IntegerVector n_neigh_sequence) {
    if (!Rf_isMatrix(nn_idx)) {
        stop("nn_idx must be a matrix.");
    }
    const int n = Rf_nrows(nn_idx);
    const int max_cols = Rf_ncols(nn_idx);
    const int n_steps = n_neigh_sequence.size();
    NumericVector out(n_steps);

    DisjointSet uf(n);
    int prev_k = 0;
    const int type = TYPEOF(nn_idx);
    const int* iptr = (type == INTSXP) ? INTEGER(nn_idx) : nullptr;
    const double* dptr = (type == REALSXP) ? REAL(nn_idx) : nullptr;
    if (iptr == nullptr && dptr == nullptr) {
        stop("nn_idx must be an integer or numeric matrix.");
    }

    for (int s = 0; s < n_steps; ++s) {
        int cur_k = n_neigh_sequence[s];
        if (cur_k < 0) {
            cur_k = 0;
        }
        if (cur_k > max_cols) {
            cur_k = max_cols;
        }
        if (cur_k < prev_k) {
            uf.reset(n);
            prev_k = 0;
        }
        for (int col = prev_k; col < cur_k; ++col) {
            const size_t col_off = static_cast<size_t>(col) * static_cast<size_t>(n);
            if (iptr != nullptr) {
                for (int row = 0; row < n; ++row) {
                    const int raw = iptr[col_off + row];
                    if (raw != NA_INTEGER && raw >= 1 && raw <= n) {
                        uf.unite(row, raw - 1);
                    }
                }
            } else {
                for (int row = 0; row < n; ++row) {
                    const double raw = dptr[col_off + row];
                    if (R_FINITE(raw) && raw >= 1.0 && raw <= static_cast<double>(n)) {
                        const int nbr = static_cast<int>(raw) - 1;
                        if (nbr >= 0 && nbr < n) {
                            uf.unite(row, nbr);
                        }
                    }
                }
            }
        }
        out[s] = static_cast<double>(uf.components);
        prev_k = cur_k;
    }

    return out;
}
