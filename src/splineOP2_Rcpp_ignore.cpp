#if (NEVER_DEFINED_FLAG)
#include <Rcpp.h>
#include "splineOP.h"
#include "Matrix.h"
#include "Rcpp_Wrap_Specialization.h"

Rcpp::NumericMatrix dp_matrix_Rcpp(std::vector<double> data,
                                   double beta,
                                   int S = 10,
                                   int nb_initSpeed = 5,
                                   double data_var = 1.0)
{
  // Call your pure C++ backend
  spop::Matrix<double> result = dp_matrix(data, beta, S, nb_initSpeed, data_var);

  // Return it to R
  return Rcpp::wrap(result);
}

#endif
