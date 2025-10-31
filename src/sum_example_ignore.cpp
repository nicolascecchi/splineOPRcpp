
#include <Rcpp.h>
using namespace Rcpp;

//
//' Compute the sum of a numeric vector (C++ implementation)
//'
//' This function computes the sum of a numeric vector using a standard C++
//' vector structure. It provides performance comparable to the built-in
//' \code{sum()} function in R and serves as an example of integrating
//' native C++ code with R using Rcpp.
//'
//' @param v A numeric vector.
//'
//' @return A single numeric value: the sum of all elements in \code{v}.
//'
//' @examples
//' x <- runif(10)
//' cpp_sum2(x)
//'
//' @export
// [[Rcpp::export]]
double cpp_sum2(NumericVector v)
{
  double total = 0.0;
  int n = v.size();
  for (int i = 0; i < n; ++i)
  {
    total += v[i]; // accesses R memory directly
  }
  return total;
}




