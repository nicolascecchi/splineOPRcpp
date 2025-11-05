#include <cmath> // Required for std::pow
#include <cassert>
#include <Rcpp.h>
#include <RcppEigen.h>

#include "Faulhaber.h"
#include "QuadraticCost.h"


// Constructor: precompute cumulative sums depending on y
QuadraticCost::QuadraticCost(Eigen::MatrixXd data)
  : y(data.rows(), data.cols()),
    nobs{data.cols()},
    ndims{data.rows()},
    // Folowing are 0-initialized by default
    cumsum_y(ndims,nobs+1),
    cumsum_y2(ndims,nobs+1),
    cumsum_yL1(ndims,nobs+1),
    cumsum_yL2(ndims,nobs+1){
  // Fill cumsum one dimension at a time
  y = data;
  for (size_t j = 0; j < ndims; ++j)  // Loop over dimensions remove cast?
  {
    for (size_t i = 0; i < nobs; ++i) // start at 1 since first is dummy at value=0 
    {
      cumsum_y(j,i+1)    = cumsum_y(j,i)    + y(j,i);
      cumsum_y2(j,i+1)   = cumsum_y2(j,i)   + y(j,i) * y(j,i);
      cumsum_yL1(j,i+1)  = cumsum_yL1(j,i)  + y(j,i) * i;
      cumsum_yL2(j,i+1)  = cumsum_yL2(j,i)  + y(j,i) * i * i;
    }
  }
}

// --- Compute cost C_{s:t}(p_s, p_t, v_t) ---
double QuadraticCost::interval_cost( size_t s
                                   , size_t t
                                   , const Eigen::VectorXd& p_s
                                   , const Eigen::VectorXd& p_t
                                   , const Eigen::VectorXd& v_s)
{
  //assert(t > s && s >= 0 && t <= this->nobs);
  int n = t - s;
  double mse = 0.0;

  // Coefficients of the quadratic p(x) = a(x - x_s)^2 + b(x - x_s) + c
  double L = static_cast<double>(n);
  Eigen::ArrayXd a = 2./std::pow(L,2) * (p_t - p_s - v_s*L);
  Eigen::ArrayXd b = v_s.eval();
  Eigen::ArrayXd c = p_s.eval();
  //// Retrieve y-based sums from cumulative arrays.
  //// These are effectively cusum[t-1] - cusum[s-1] 
  //// when we think our mathematical cost function
  Eigen::ArrayXd sum_y    = (cumsum_y.col(t)    - cumsum_y.col(s)).array();
  Eigen::ArrayXd sum_y2   = (cumsum_y2.col(t)   - cumsum_y2.col(s)).array();
  Eigen::ArrayXd sum_yL1  = ((cumsum_yL1.col(t)  - cumsum_yL1.col(s)).array() - (s * sum_y)).array();
  Eigen::ArrayXd aux = (cumsum_yL1.col(t)  - cumsum_yL1.col(s)).array();
  Eigen::ArrayXd sum_yL2  = (cumsum_yL2.col(t)  - cumsum_yL2.col(s)).array()-2*s*aux+std::pow(s,2)*sum_y;
  //// Compute L-based sums via Faulhaber
  double sum_L1 = Faulhaber(n-1,1); //#S1(n-1);
  double sum_L2 = Faulhaber(n-1,2); //#S2(n-1);
  double sum_L3 = Faulhaber(n-1,3); //#S3(n-1);
  double sum_L4 = Faulhaber(n-1,4); //#S4(n-1);
  //// Expanded quadratic cost
  Eigen::ArrayXd dimensionCosts = Eigen::ArrayXd::Zero(y.rows());
  dimensionCosts += a * a * sum_L4 / 4.;
  dimensionCosts += a * b * sum_L3;
  dimensionCosts += (a * c + b * b) * sum_L2;
  dimensionCosts += 2.0 * b * c * sum_L1;
  dimensionCosts += c * c * n;
  dimensionCosts -=  a * sum_yL2;
  dimensionCosts -= 2.0 * b * sum_yL1;
  dimensionCosts -= 2.0 * c * sum_y;
  dimensionCosts += sum_y2;
  mse = dimensionCosts.sum()/n;
  return mse;
}


