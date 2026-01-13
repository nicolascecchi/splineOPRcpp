#include <cmath> // Required for std::pow
#include <cassert>
#include <Eigen/Dense>
//#include <Rcpp.h>
//#include <RcppEigen.h>

#include "Faulhaber.h"
#include "QuadraticCost.h"


// Constructor: precompute cumulative sums depending on y
QuadraticCost::QuadraticCost(Eigen::MatrixXd data)
  : y(data.rows(), data.cols()),
    nobs{static_cast<int>(data.cols())},
    ndims{static_cast<int>(data.rows())},
    // Folowing are 0-initialized by default
    cumsum_y{Eigen::MatrixXd::Zero(ndims, nobs + 1)},
    cumsum_y2{Eigen::MatrixXd::Zero(ndims, nobs + 1)},
    cumsum_yL1{Eigen::MatrixXd::Zero(ndims, nobs + 1)},
    cumsum_yL2{Eigen::MatrixXd::Zero(ndims, nobs + 1)},
    subsum_y{Eigen::VectorXd::Zero(ndims)},
    subsum_y2{Eigen::VectorXd::Zero(ndims)},
    subsum_yL1{Eigen::VectorXd::Zero(ndims)},
    subsum_yL2{Eigen::VectorXd::Zero(ndims)}
  {
  // Fill cumsum one dimension at a time
  y = data;
  a = Eigen::ArrayXd::Zero(ndims);
  b = Eigen::ArrayXd::Zero(ndims);
  c = Eigen::ArrayXd::Zero(ndims);

  for (size_t j = 0; j < ndims; ++j)  // Loop over dimensions remove cast?
  {
    for (size_t i = 0; i < nobs; ++i) // start at 1 since first is dummy at value=0 
    {
      cumsum_y(j,   i + 1) = cumsum_y(j, i)  + y(j, i);
      cumsum_y2(j,  i + 1) = cumsum_y2(j,i)  + y(j, i) * y(j, i);
      cumsum_yL1(j, i + 1) = cumsum_yL1(j,i) + y(j, i) * i;
      cumsum_yL2(j, i + 1) = cumsum_yL2(j,i) + y(j, i) * i * i;
    }
  }
}
///////////////////////////////////////////////////////////////////////////////////////////////////
// --- Compute cost C_{s:t}(p_s, p_t, v_t) ---
double QuadraticCost::interval_cost( size_t s
                                   , size_t t
                                   , const Eigen::VectorXd& p_s
                                   , const Eigen::VectorXd& p_t
                                   , const Eigen::VectorXd& v_s)
{
  //assert(t > s && s >= 0 && t <= this->nobs);
  int n = t - s;
  // Compute sums of integers with Faulhaber 
  double sum_L1 = Faulhaber(n - 1, 1); // #S1(n-1);
  double sum_L2 = Faulhaber(n - 1, 2); // #S2(n-1);
  double sum_L3 = Faulhaber(n - 1, 3); // #S3(n-1);
  double sum_L4 = Faulhaber(n - 1, 4); // #S4(n-1);
  
  // Coefficients of the quadratic p(x) = a(x - x_s)^2 + b(x - x_s) + c
  double L = static_cast<double>(n);
  a = 2./std::pow(L,2) * (p_t - p_s - v_s*L);
  b = v_s.eval();
  c = p_s.eval();

  //// Retrieve y-based sums from cumulative arrays.
  //// These are effectively cusum[t-1] - cusum[s-1] 
  //// when we think our mathematical cost function
  
  //// PRIORITE 1 : REDUIRE A UN SEUL CAS AVEC UNE COLONNE EN PLUS ???
  if (s==0)
  {
    subsum_y = cumsum_y.col(t - 1);
    subsum_y2 = cumsum_y2.col(t - 1);
    subsum_yL1 = cumsum_yL1.col(t - 1);
    aux = cumsum_yL1.col(t - 1);
    subsum_yL2  = cumsum_yL2.col(t - 1);  
  }
  else
  {
    subsum_y = cumsum_y.col(t - 1)    - cumsum_y.col(s - 1);
    subsum_y2 = cumsum_y2.col(t - 1)   - cumsum_y2.col(s - 1);
    subsum_yL1 = (cumsum_yL1.col(t - 1)  - cumsum_yL1.col(s - 1)) - (s * subsum_y);
    aux = cumsum_yL1.col(t - 1)  - cumsum_yL1.col(s - 1);
    subsum_yL2  = (cumsum_yL2.col(t - 1)  - cumsum_yL2.col(s - 1)) - 2*s*aux+std::pow(s,2)*subsum_y;  
  }
  
  //// Expanded quadratic cost
  // We use type ArrayXd because that allows for the behaviour
  // of element-wise multiplication (instead of dot product with VectorXd)
  
  dimensionCosts = a * a * sum_L4 / 4.;
  dimensionCosts += a * b * sum_L3;
  dimensionCosts += (a * c + b * b) * sum_L2;
  dimensionCosts += 2.0 * b * c * sum_L1;
  dimensionCosts += c * c * n;
  dimensionCosts -=  a * subsum_yL2.array();
  dimensionCosts -= 2.0 * b * subsum_yL1.array();
  dimensionCosts -= 2.0 * c * subsum_y.array();
  dimensionCosts += subsum_y2.array();
  
  double cost = dimensionCosts.sum();
  return cost;
}


