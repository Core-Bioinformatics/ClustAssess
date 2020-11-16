#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double calculate_pac_cpp(IntegerMatrix indicator,
                                NumericMatrix connectivity,
                                double lower_lim,
                                double upper_lim){

  int n_samples = indicator.nrow();
  double fraction_inside_interval = 0;

  for (int i_cell=0; i_cell<(n_samples-1); i_cell++){
    for (int j_cell=(i_cell+1); j_cell<n_samples; j_cell++){

      // make sure we don't divide by zero
      if (indicator(i_cell, j_cell) != 0){
        double cell_fraction = connectivity(i_cell, j_cell) /
          indicator(i_cell, j_cell);

        if (cell_fraction > lower_lim && cell_fraction < upper_lim){
          fraction_inside_interval++;
        }
      }
    }
    checkUserInterrupt();
  }
  fraction_inside_interval = fraction_inside_interval /
    (n_samples * (n_samples-1) / 2);
  return fraction_inside_interval;
}
