#include <Rcpp.h>
using namespace Rcpp;

//' Square Each Element of a States Matrix
//'
//' @description
//' Given a matrix of states of size S x N, this function returns a new matrix
//' where each element is squared. Additional arguments `data` and `beta` are
//' included for future extensions but are not used in the squaring.
//'
//' @param data Numeric vector of length N (currently unused).
//' @param beta Numeric penalty value (currently unused).
//' @param states Numeric matrix of size S x N.
//'
//' @return Numeric matrix of size S x N where each element is squared.
//'
//' @examples
//' library(Rcpp)
//' # Create example states matrix
//' states <- matrix(1:6, nrow = 2, ncol = 3)
//' data <- c(10, 20, 30)
//' beta <- 0.5
//' # Call the function
//' square_states(data, beta, states)
//'
//' @export
// [[Rcpp::export]]
NumericMatrix square_states(NumericVector data, double beta, NumericMatrix states)
{
  int S = states.nrow();
  int N = states.ncol();

  NumericMatrix result(S, N);
  double current_MIN;
  double candidate;

  /// INITIAL COLUMN
  for(int i = 0; i < S; i++)
  {
    result(i,0) = 100000;
  }


  ///
  /// UPDATE RULE
  ///
  ///
  for(int t = 1; t < N; t++)
  {
    for(int j = 0; j < S; j++)
    {
      /////// ///////
      /////// ///////
      /// compute Q state k, time l
      /////// ///////
      /////// ///////

      current_MIN = 100000;
      for(int s = 0; s < t; s++) /// time s smaller than t
      {
        for(int i = 0; i < S; i++) // for all states
        {
          /// find candidate
          candidate = (t-s)*(S-i);
          if(candidate < current_MIN)
          {
            /// UPDATE
            current_MIN = candidate;
          }
        }
      }

      result(j, t) = current_MIN;
      // to do add best_s, best_i (matrix B)

      /////// ///////
      /////// ///////
    }
  }

  return result;
}
