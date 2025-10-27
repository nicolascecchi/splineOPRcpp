#include <Rcpp.h>
#include "splineOP.h"

// [[Rcpp::export]]
Rcpp::IntegerVector dp_matrix_Rcpp(std::vector<double> data,
                                   double beta,
                                   int S = 10,
                                   int nb_initSpeed = 5,
                                   double data_var = 1.0)
{
  // Call your pure C++ backend
  std::vector<int> result = dp_matrix(data, beta, S, nb_initSpeed, data_var);

  // Return it to R
  return Rcpp::wrap(result);
}
