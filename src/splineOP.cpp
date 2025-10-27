#include <vector>
#include <limits>
#include "Matrix.h"
#include "QuadraticCost.h"

// -----------------------------------------------------------------------------
// Dynamic Programming with Quadratic Cost Integration
// -----------------------------------------------------------------------------
std::vector<int> dp_matrix(const std::vector<double>& data,
                           double beta,
                           int S = 10,
                           int nb_initSpeed = 5)
{
  int N = data.size();
  QuadraticCost qc(data); // Precompute all cumulative sums ONCE

  Matrix value(S, N, std::numeric_limits<double>::infinity());
  Matrix argmin_i(S, N, -1);
  Matrix argmin_s(S, N, -1);

  // Initialization (first column)
  for (int j = 0; j < S; j++)
  {
    value(j, 0) = 0.0;
  }

  // ---------------------------------------------------------------------------
  // Dynamic programming loop
  // ---------------------------------------------------------------------------
  for (int t = 1; t < N; t++)
  {
    for (int j = 0; j < S; j++)
    { // current state
      double current_MIN = std::numeric_limits<double>::infinity();
      int best_i = -1;
      int best_s = -1;

      for (int s = 0; s < t; s++)
      {   // previous time
        for (int i = 0; i < S; i++)
        { // previous state
          /// TODO TODO TODO TODO TODO TODO TODO
          // Example parameters (you can model how p_s, p_t, v_t depend on states)
          double p_s = static_cast<double>(i);
          double p_t = static_cast<double>(j);
          double v_t = static_cast<double>(j - i); // or use any other state→velocity rule

          // Cost for interval [s, t)
          double interval_cost = qc.quadratic_cost_interval(s, t, p_s, p_t, v_t);

          // Regularization term (β penalty on transitions)
          double candidate = interval_cost + beta + value(i, s);

          if (candidate < current_MIN)
          {
            current_MIN = candidate;
            best_i = i;
            best_s = s;
          }
        }
      }

      value(j, t) = current_MIN;
      argmin_i(j, t) = best_i;
      argmin_s(j, t) = best_s;
    }
  }

  // ---------------------------------------------------------------------------
  // Backtracking: recover optimal path
  // ---------------------------------------------------------------------------
  std::vector<int> optimal_path(N, -1);

  // Find best final state
  double min_final = std::numeric_limits<double>::infinity();
  int best_final = -1;
  for (int j = 0; j < S; j++) {
    if (value(j, N - 1) < min_final) {
      min_final = value(j, N - 1);
      best_final = j;
    }
  }
  optimal_path[N - 1] = best_final;

  // Backtrack in time
  for (int t = N - 1; t > 0; t--)
  {
    int j = optimal_path[t];
    if (j < 0) break;
    optimal_path[t - 1] = argmin_i(j, t);
  }

  return optimal_path;
}

