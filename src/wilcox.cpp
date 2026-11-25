
#include <Rcpp.h>
#include <map>
#include <unordered_map>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <vector>

using namespace Rcpp;
using namespace std::chrono;

// [[Rcpp::export]]
NumericVector wilcox_test(IntegerMatrix rank_values, int n1, int max_rank) {
    if (n1 < 0 || n1 > rank_values.ncol()) {
        stop("n1 must be between 0 and ncol(rank_values).");
    }
    if (max_rank < 0) {
        stop("max_rank must be non-negative.");
    }
    double n = rank_values.ncol(), m = rank_values.nrow();
    double n1_d = static_cast<double>(n1);
    double n2 = n - n1_d;
    long long value, n_unique;
    double mu = (n1_d * n2) / 2.0, sigma2 = (n1_d * n2 * (n + 1.0)) / 12.0;
    double U = n1_d * n2 + n1_d * (n1_d + 1.0) / 2.0 - mu, factor_n = n * (n + 1.0) * (n - 1.0);
    double freq_val, sum_r1, adjustment;
    double calculated_rank, start, less_val, temp;

    NumericVector zuppertail(m), zlowertail(m); 

    std::vector<int> unique_rankings(static_cast<size_t>(n) + 1, 0);
    std::vector<double> new_ranks(static_cast<size_t>(n) + 1, 0.0);
    std::vector<int> mappings(static_cast<size_t>(max_rank) + 1, 0);
    std::vector<int> freq(static_cast<size_t>(max_rank) + 1, 0);

    for (int j = 0; j < m; j++) {
        n_unique = 0;
        std::fill(freq.begin(), freq.end(), 0);
        for (int i = 0; i < n; i++) {
            value = rank_values(j, i);
            if (value == NA_INTEGER || value < 0 || value > max_rank) {
                stop("rank_values entries must be integers in [0, max_rank].");
            }
            freq[value]++;

            if (freq[value] == 1) {
                unique_rankings[n_unique] = value;
                n_unique++;
            }
        }

        std::sort(unique_rankings.begin(), unique_rankings.begin() + n_unique);

        start = 1;
        for (int i = 0; i < n_unique; i++) {
            freq_val = freq[unique_rankings[i]];
            if (freq_val == 1) {
                new_ranks[i] = start;
                start++;
            } else {
                calculated_rank = ((start+freq_val)*(start+freq_val-1) / 2 - (start-1) * start / 2) / 
                freq_val;
                new_ranks[i] = calculated_rank;
                start += freq_val;
            }

            mappings[unique_rankings[i]] = i;
        }

        sum_r1 = 0;
        for (int i = 0; i < n1; i++) {
            sum_r1 += new_ranks[mappings[rank_values(j, i)]];
        }

        if (n_unique != n) {
            adjustment = 0; 

            for (int i = 0; i < n_unique; i++) {
                freq_val = freq[unique_rankings[i]];
                adjustment += freq_val * (freq_val - 1) * (freq_val + 1);
            }

            adjustment /= factor_n;

            sigma2 = sigma2 * (1 - adjustment);
        }

        zlowertail(j) = (U - sum_r1 + 0.5) / sqrt(sigma2); 
        zuppertail(j) = (U - sum_r1 - 0.5) / sqrt(sigma2); 

        if (n_unique != n) {
            sigma2 = (n1_d * n2 * (n + 1.0)) / 12.0;
        }

    }

    NumericVector pvalue_less = pt(zuppertail, R_PosInf, false);
    NumericVector pvalue_greater = pt(zlowertail, R_PosInf, true);

    for (int i = 0; i < m; i++) {
        less_val = pvalue_less(i), temp = pvalue_greater(i);
        
        if (temp < less_val) {
            less_val = temp;
        }

        less_val *= 2;

        if (less_val > 1) {
            less_val = 1;
        }

        pvalue_less[i] = less_val;
    }


    return(pvalue_less);
}
