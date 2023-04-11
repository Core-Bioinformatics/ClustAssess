
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
    // int duration_actual_rankings = 0, duration_initial_rankings = 0, duration_sum = 0, duration_adjustment = 0;
    double freq_val, sum_r1, adjustment;
    double calculated_rank, start, less_val, temp;

    NumericVector zuppertail(m), zlowertail(m); // result(m);

    int unique_rankings[(int) n+1] = {0};
    double new_ranks[(int) n+1] = {0};
    int mappings[max_rank+1] = {0};


    for (int j = 0; j < m; j++) {
        // auto start_time = high_resolution_clock::now();
        // int freq[max_rank+1] = {0};
        // std::cout << initial_rankings.size() << '\n';
        int max = 0;

        n_unique = 0;
        int freq[max_rank+1] = {0};
        for(int i = 0; i < n; i++) {
            value = rank_values(j, i);
            freq[value]++;
            
            if (freq[value] == 1) {
                // unique_rankings.push_back(value);
                unique_rankings[n_unique] = value;
                n_unique++;
            }
        }

        std::sort(unique_rankings, unique_rankings+n_unique);

        // auto stop_time = high_resolution_clock::now();
        // auto duration_time = duration_cast<microseconds>(stop_time - start_time);
        
        // duration_initial_rankings += duration_time.count();
        // start_time = high_resolution_clock::now();
        start = 1;
        for (int i = 0; i < n_unique; i++) {
            freq_val = freq[unique_rankings[i]];
            if (freq_val == 1) {
                // new_ranks[i] = 1;
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

        // stop_time = high_resolution_clock::now();
        // duration_time = duration_cast<microseconds>(stop_time - start_time);
        // duration_actual_rankings += duration_time.count();

        sum_r1 = 0;
        // start_time = high_resolution_clock::now();
        for (int i = 0; i < n1; i++) {
            sum_r1 += new_ranks[mappings[rank_values(j, i)]];

            // if (j == 3497) {
            //     std::cout << rank_values(j, i) << ' ' << mappings[rank_values(j,i)] << ' ' << new_ranks[mappings[rank_values(j, i)]] << '\n';
            // }
        }


        // if (j == 0) {
        //     for(int i = 0; i < n_unique; i++) {
        //         std::cout << unique_rankings[i] << ' ' << new_ranks[i] << ' ' << mappings[unique_rankings[i]] << '\n';
        //     }
        // }

        // stop_time = high_resolution_clock::now();
        // duration_time = duration_cast<microseconds>(stop_time - start_time);
        // duration_sum += duration_time.count();


        // start_time = high_resolution_clock::now();
        if (n_unique != n) {
            adjustment = 0; 

            for (int i = 0; i < n_unique; i++) {
                freq_val = freq[unique_rankings[i]];
                // if (j == 3497) {
                //     std::cout << i << ' ' << freq_val << ' ' << factor_n << '\n';
                // }
                adjustment += freq_val * (freq_val - 1) * (freq_val + 1);
            }

            adjustment /= factor_n;

            sigma2 = sigma2 * (1 - adjustment);
        }
        // stop_time = high_resolution_clock::now();
        // duration_time = duration_cast<microseconds>(stop_time - start_time);
        // duration_adjustment += duration_time.count();


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
        // pvalue_less[i] = zlowertail(i);
    }


    return(pvalue_less);
}

//  NumericVector wilcox_test(IntegerMatrix rank_values, int n1, int max_rank) {
//     std::unordered_map<int, int> initial_rankings;
//     std::unordered_map<double, int> actual_rankings;
//     std::unordered_map<int, double> ranking_mappings;

//     int n = rank_values.ncol(), m = rank_values.nrow();
//     int n2 = n - n1;
//     double mu = n1 * n2 / 2, sigma2;
//     double U = n1 * n2 + n1 * (n1 + 1) / 2, factor_n = n * (n+1) * (n-1);
//     int duration_actual_rankings = 0, duration_initial_rankings = 0, duration_sum = 0, duration_adjustemt = 0,freq_val;
//     double calculated_rank, start;

//     NumericVector zuppertail(m), zlowertail(m), result(m);

//     for (int j = 0; j < m; j++) {
//         sigma2 = n1 * n2 * (n+1) / 12;
//         auto start_time = high_resolution_clock::now();

//         for(int i = 0; i < n; i++) {
//             initial_rankings[rank_values(j,i)]++;

//         }

//         // std::cout << max << ' ';
//         // std::cout << "OUT\n";
//         auto stop_time = high_resolution_clock::now();
//         auto duration_time = duration_cast<microseconds>(stop_time - start_time);
        
//         duration_initial_rankings += duration_time.count();

        
//         start_time = high_resolution_clock::now();
//         start = 1;
//         for (auto p : initial_rankings) {
//             if (p.second == 1) {
//                 actual_rankings[start] = 1;
//                 ranking_mappings[p.first] = start;
//                 if (j == 0) {
//                     std::cout << start << ' ';
//                 }
//                 start++;
//             } else {
//                 calculated_rank = ((start+p.second)*(start+p.second-1) / 2 - (start-1) * start / 2) / p.second;
//                 actual_rankings[calculated_rank] = p.second;
//                 ranking_mappings[p.first] = calculated_rank;
//                 start += p.second;
//                 if (j == 0) {
//                     std::cout << calculated_rank << ' ';
//                 }
//             }
//         }
//         stop_time = high_resolution_clock::now();
//         duration_time = duration_cast<microseconds>(stop_time - start_time);

//         duration_actual_rankings += duration_time.count();

//         double sum_r1 = 0;
//         start_time = high_resolution_clock::now();
//         for (int i = 0; i < n1; i++) {
//             sum_r1 += ranking_mappings[rank_values(j, i)];
//         }
//         stop_time = high_resolution_clock::now();
//         duration_time = duration_cast<microseconds>(stop_time - start_time);
//         duration_sum += duration_time.count();

//         // std::cout << "sum_r1" << ' ' << sum_r1 << '\n';
        
//         // if (j < 5) {
//         //     std::cout << sum_r1 << ' ';
//         // }

//         start_time = high_resolution_clock::now();
//         if (actual_rankings.size() != n) {
//             double adjustment = 0; 

//             for (auto p : actual_rankings) {
//                 adjustment += p.second * (p.second - 1) * (p.second + 1) / factor_n;
//             }

//             // std::cout << "TIE " << adjustment << ' ' << sigma2 << '\n';

//             sigma2 = sigma2 * (1 - adjustment);
//         }
//         stop_time = high_resolution_clock::now();
//         duration_time = duration_cast<microseconds>(stop_time - start_time);
//         duration_adjustemt += duration_time.count();

//         zlowertail(j) = (U - sum_r1 ) / sqrt(sigma2); 
//         zuppertail(j) = (U - sum_r1 ) / sqrt(sigma2); 
//         initial_rankings.clear();
//         actual_rankings.clear();
//         ranking_mappings.clear();
//     }

//     std::cout << "Initial rankings " << duration_initial_rankings << '\n';
//     std::cout << "Actual rankings " << duration_actual_rankings << '\n';
//     std::cout << "Sum r1 " << duration_sum << '\n';
//     std::cout << "Adjustments " << duration_adjustemt << '\n';

//     auto start_time = high_resolution_clock::now();
//     NumericVector pvalue_less = pt(zuppertail, R_PosInf, false);
//     NumericVector pvalue_greater = pt(zlowertail, R_PosInf, true);
//     auto stop_time = high_resolution_clock::now();
//     auto duration = duration_cast<microseconds>(stop_time - start_time);
//     std::cout << "T-student distr " << duration.count() << '\n';

//     // std::cout << zuppertail(0) << ' ' << zlowertail(0) << '\n';
//     // std::cout << pvalue_less(0) << ' ' <<pvalue_greater(0) << '\n';

//     start_time = high_resolution_clock::now();
//     for (int i = 0; i < m; i++) {
//         // double temp = pvalue_less(i);
//         // if (pvalue_greater(i) < temp) {
//         //     temp = pvalue_greater(i);
//         // }
//         // temp *= 2;

//         // if (temp > 1) {
//         //     temp = 1;
//         // }

//         // result(i) = temp;
//     }
//     stop_time = high_resolution_clock::now();
//     duration = duration_cast<microseconds>(stop_time - start_time);
//     std::cout << "Calculate 2*.." << duration.count() << '\n';

//     return(result);
// }