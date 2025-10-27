#include <iostream>
#include <vector>
#include <limits>
#include <random>      // for random number generation
#include "Matrix.h"
#include "QuadraticCost.h"


std::vector<int> dp_matrix(const std::vector<double>& data,
                           double beta,
                           int S = 10,
                           int nb_initSpeed = 5,
                           double data_var = 1.0)
{
  int N = data.size();
  QuadraticCost qc(data); // Precompute all cumulative sums once

  Matrix<double> speeds(S, N, 0.0);

  //////////////////////////////////////////////////////////////////////////////
  // ---------------------------------------------------------------------------
  // Initialize Gaussian random matrix: S x N with iid N(data[i], data_var)
  // ---------------------------------------------------------------------------
  // rand_init: S x N
  Matrix<double> states(S, N, 0.0);

  std::random_device rd;
  std::mt19937 gen(rd());

  for (int t = 0; t < N; t++)
  {
    // first row = data[t]
    states(0, t) = data[t];

    // remaining rows: Gaussian centered on data[j]
    std::normal_distribution<double> normal_dist(data[t], std::sqrt(data_var));
    for (int j = 1; j < S; j++)
    {
      states(j, t) = normal_dist(gen);
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // ---------------------------------------------------------------------------
  // Allocate DP matrices
  // ---------------------------------------------------------------------------
  Matrix<double> value(S, N, std::numeric_limits<double>::infinity());
  Matrix<int> argmin_i(S, N, -1);
  Matrix<int> argmin_s(S, N, -1);

  // Initialization (first column)
  for (int j = 0; j < S; j++)
  {
    value(j, 0) = 0;  // random initial cost for 0 data point
  }

  //////////////////////////////////////////////////////////////////////////////
  // ---------------------------------------------------------------------------
  // Dynamic programming loop
  // ---------------------------------------------------------------------------
  for (int t = 1; t < N; t++)
  {
    for (int j = 0; j < S; j++)
    { // current state
      double current_MIN = std::numeric_limits<double>::infinity();
      double best_speed;
      int best_i = -1;
      int best_s = -1;

      for (int s = 0; s < t; s++)
      {   // previous time
        for (int i = 0; i < S; i++)
        { // previous state
          // Example mapping of state indices to spline parameters
          double p_s = states(i, s);
          double p_t = states(j, t);
          double v_t = 2*(p_t - p_s)/(t - s) - speeds(i, s); // simple slope rule

          // Quadratic cost for interval [s, t)
          double interval_cost = qc.quadratic_cost_interval(s, t, p_s, p_t, v_t);

          // Candidate cost (DP recurrence)
          double candidate = value(i, s) + interval_cost + beta;

          if (candidate < current_MIN)
          {
            current_MIN = candidate;
            best_speed = v_t;
            best_i = i;
            best_s = s;
          }
        }
      }

      value(j, t) = current_MIN;
      speeds(j, t) = best_speed;
      argmin_i(j, t) = best_i;
      argmin_s(j, t) = best_s;
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // ---------------------------------------------------------------------------
  // Backtracking: recover change points
  // ---------------------------------------------------------------------------
  std::vector<int> change_points;
  change_points.push_back(N - 1);  // last time index is always a change point

  // Find best final state
  double min_final = std::numeric_limits<double>::infinity();
  int best_final_state = -1;
  for (int j = 0; j < S; j++)
  {
    if (value(j, N - 1) < min_final)
    {
      min_final = value(j, N - 1);
      best_final_state = j;
    }
  }

  int j = best_final_state;
  int t = N - 1;

  // Backtrack using argmin_s and argmin_i
  while (true)
  {
    int s_prev = argmin_s(j, t);
    int i_prev = argmin_i(j, t);

    if (s_prev < 0)
      break;  // reached the beginning or invalid index

    change_points.push_back(s_prev);  // record change boundary

    t = s_prev;
    j = i_prev;
  }

  // Reverse the order to chronological (0 → N)
  std::reverse(change_points.begin(), change_points.end());

  return change_points;
}


