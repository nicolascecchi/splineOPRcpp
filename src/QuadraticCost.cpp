#include "QuadraticCost.h"
#include <cassert>

// Constructor: precompute cumulative sums depending on y
QuadraticCost::QuadraticCost(const std::vector<double>& data)
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

// --- Faulhaber formulas for sums of integer powers ---
double QuadraticCost::S1(int n)
{
  return 0.5 * n * (n - 1);
}

double QuadraticCost::S2(int n)
{
  return n * (n - 1) * (2 * n - 1) / 6.0;
}

double QuadraticCost::S3(int n)
{
  double x = 0.5 * n * (n - 1);
  return x * x;
}

double QuadraticCost::S4(int n)
{
  double nn = static_cast<double>(n);
  return (nn - 1) * nn * (2 * nn - 1) * (3 * nn * nn - 3 * nn - 1) / 30.0;
}

// --- Compute cost C_{s:t}(p_s, p_t, v_t) ---
double QuadraticCost::quadratic_cost_interval(int s, int t, double p_s, double p_t, double v_t) const
{
  assert(t > s && s >= 0 && t <= N);
  int n = t - s;

  // Coefficients of the quadratic p(x) = a(x - x_s)^2 + b(x - x_s) + c
  double L = static_cast<double>(n);
  double a = (v_t * L - 2.0 * (p_t - p_s)) / (L * L);
  double b = v_t - 2.0 * a * L;
  double c = p_s;

  // Retrieve y-based sums from cumulative arrays
  double sum_y    = cumsum_y[t]    - cumsum_y[s];
  double sum_y2   = cumsum_y2[t]   - cumsum_y2[s];
  double sum_yL1  = cumsum_yL1[t]  - cumsum_yL1[s];
  double sum_yL2  = cumsum_yL2[t]  - cumsum_yL2[s];

  // Compute L-based sums via Faulhaber
  double sum_L1 = S1(n);
  double sum_L2 = S2(n);
  double sum_L3 = S3(n);
  double sum_L4 = S4(n);

  // Expanded quadratic cost
  double cost = 0.0;
  cost += a * a * sum_L4;
  cost += 2.0 * a * b * sum_L3;
  cost += (2.0 * a * c + b * b) * sum_L2;
  cost += 2.0 * b * c * sum_L1;
  cost += c * c * n;
  cost -= 2.0 * a * sum_yL2;
  cost -= 2.0 * b * sum_yL1;
  cost -= 2.0 * c * sum_y;
  cost += sum_y2;

  return cost;
}

