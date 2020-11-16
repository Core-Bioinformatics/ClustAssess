#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix update_connectivity_cpp(NumericMatrix connectivity,
                         IntegerVector sampling_indices,
                         IntegerVector cluster_assignments){
  int n_samples = sampling_indices.size();
  for (int i_cell=0; i_cell<n_samples; i_cell++){
    int cell_1 = sampling_indices[i_cell] - 1;
    for (int j_cell=(i_cell+1); j_cell<n_samples; j_cell++){
      if (j_cell>n_samples){
        break;
      }
      //Rcout << sampling_indices[i_cell];
      //Rcout << sampling_indices[j_cell];
      //Rcout << cluster_assignments[i_cell];
      //Rcout << cluster_assignments[j_cell];
      int cell_2 = sampling_indices[j_cell] - 1;
      if (cluster_assignments[i_cell] == cluster_assignments[j_cell]){
        connectivity(cell_1, cell_2)++;
        connectivity(cell_2, cell_1)++;
      }
    }
    checkUserInterrupt();
  }
  return connectivity;
}
