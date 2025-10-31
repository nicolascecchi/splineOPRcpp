#include <cmath> // Required for std::pow
#include "QuadraticCost.h"
#include <cassert>
#include <Rcpp.h>
#include "Faulhaber.h"

// Constructor: precompute cumulative sums depending on y
QuadraticCost::QuadraticCost(std::vector<double> data)
  : y(data),
    N(static_cast<int>(data.size())),
    cumsum_y(N+1, 0.0),
    cumsum_y2(N+1, 0.0),
    cumsum_yL1(N+1, 0.0),
    cumsum_yL2(N+1, 0.0)
{
  double yi;
  double I;
  for (int i = 0; i < N; ++i)
  {
    yi = y[i];
    I = i; // equispaced, L_i = i

    cumsum_y[i+1]    = cumsum_y[i]    + yi;
    cumsum_y2[i+1]   = cumsum_y2[i]   + yi * yi;
    cumsum_yL1[i+1]  = cumsum_yL1[i]  + yi * I;
    cumsum_yL2[i+1]  = cumsum_yL2[i]  + yi * I * I;
  }
}

// --- Compute cost C_{s:t}(p_s, p_t, v_t) ---
double QuadraticCost::quadratic_cost_interval(int s, int t, double p_s, double p_t, double v_s) const
{
  assert(t > s && s >= 0 && t <= N);
  int n = t - s;

  // Coefficients of the quadratic p(x) = a(x - x_s)^2 + b(x - x_s) + c
  double L = static_cast<double>(n);
  double a = 2./std::pow(L,2) * (p_t - p_s - v_s*L);
  double b = v_s;
  double c = p_s;

  // Retrieve y-based sums from cumulative arrays.
  // These are effectively cusum[t-1] - cusum[s-1] 
  // when we think our mathematical cost function
  double sum_y    = cumsum_y[t]    - cumsum_y[s];
  double sum_y2   = cumsum_y2[t]   - cumsum_y2[s];
  double sum_yL1  = (cumsum_yL1[t]  - cumsum_yL1[s]) - (s * sum_y);
  double aux = (cumsum_yL1[t]  - cumsum_yL1[s]);
  double sum_yL2  = (cumsum_yL2[t]  - cumsum_yL2[s])-2*s*aux+std::pow(s,2)*sum_y;

  // Compute L-based sums via Faulhaber
  double sum_L1 = Faulhaber(n-1,1); //#S1(n-1);
  double sum_L2 = Faulhaber(n-1,2); //#S2(n-1);
  double sum_L3 = Faulhaber(n-1,3); //#S3(n-1);
  double sum_L4 = Faulhaber(n-1,4); //#S4(n-1);

  // Expanded quadratic cost
  double cost = 0.0;
  cost += a * a * sum_L4 / 4.;
  cost += a * b * sum_L3;
  cost += (a * c + b * b) * sum_L2;
  cost += 2.0 * b * c * sum_L1;
  cost += c * c * n;
  cost -=  a * sum_yL2;
  cost -= 2.0 * b * sum_yL1;
  cost -= 2.0 * c * sum_y;
  cost += sum_y2;

  //Rcpp::Rcout << "n " << n << " " << "\n";
  //Rcpp::Rcout << "L " << L << " " << "\n";
  //Rcpp::Rcout << "a " << a << " " << "\n";
  //Rcpp::Rcout << "b " << b << " " << "\n";
  //Rcpp::Rcout << "c " << c << " " << "\n";
  //Rcpp::Rcout << "sum_y " << sum_y << " " << "\n";
  //Rcpp::Rcout << "sum_y2 " << sum_y2 << " " << "\n";
  //Rcpp::Rcout << "sum_yL1 " << sum_yL1 << " " << "\n";
  //Rcpp::Rcout << "sum_yL2 " << sum_yL2 << " " << "\n";
  //Rcpp::Rcout << "sum_L1 " << sum_L1 << " " << "\n";
  //Rcpp::Rcout << "sum_L2 " << sum_L2 << " " << "\n";
  //Rcpp::Rcout << "sum_L3 " << sum_L3 << " " << "\n";
  //Rcpp::Rcout << "sum_L4 " << sum_L4 << " " << "\n";

  return cost;
}

