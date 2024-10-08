
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
    double n = rank_values.ncol(), m = rank_values.nrow();
    double n2 = n - n1;
    long long value, n_unique;
    double mu = (n1 * n2) / 2, sigma2 = (n1 * n2 * (n+1)) / 12;
    double U = n1 * n2 + n1 * (n1 + 1) / 2 - mu, factor_n = n * (n+1) * (n-1);
    double freq_val, sum_r1, adjustment;
    double calculated_rank, start, less_val, temp;

    NumericVector zuppertail(m), zlowertail(m); 

    int* unique_rankings = new int[(int) n+1];
    for (int i = 0; i < n+1; i++) {
        unique_rankings[i] = 0;
    }
    double* new_ranks = new double[(int) n+1];
    for (int i = 0; i < n+1; i++) {
        new_ranks[i] = 0;
    }
    int* mappings = new int[max_rank+1];
    for (int i = 0; i < max_rank+1; i++) {
        mappings[i] = 0;
    }

    int* freq = new int[max_rank+1];


    for (int j = 0; j < m; j++) {
        int max = 0;

        n_unique = 0;
        for (int i = 0; i < max_rank+1; i++) {
            freq[i] = 0;
        }
        for(int i = 0; i < n; i++) {
            value = rank_values(j, i);
            freq[value]++;
            
            if (freq[value] == 1) {
                unique_rankings[n_unique] = value;
                n_unique++;
            }
        }

        std::sort(unique_rankings, unique_rankings+n_unique);

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
            sigma2 = (double) (n1 * n2 * (n+1)) / 12;
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
