#include <Rcpp.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <limits>
#include <cmath>
#include <random>      // for random number generation
#include "Matrix.h"
#include "QuadraticCost.h"


// Check if the number is denormalized

spop::Matrix<double> dp_matrix(const std::vector<double>& data,
                           double beta,
                           int S = 10,
                           int nb_initSpeed = 1,
                           double data_var = 1.0)
{
  int N = data.size();
  QuadraticCost qc(data); // Precompute all cumulative sums once

  spop::Matrix<double> speeds(S, N, 0.0);
  spop::Matrix<double> initspeeds(nb_initSpeed,1,-0.017485295);


  //////////////////////////////////////////////////////////////////////////////
  // ---------------------------------------------------------------------------
  // Initialize Gaussian random spop::Matrix: S x N with iid N(data[i], data_var)
  // ---------------------------------------------------------------------------
  // rand_init: S x N
  spop::Matrix<double> states(S, N, 0.0);

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
  spop::Matrix<double> value(S, N, std::numeric_limits<double>::infinity());
  spop::Matrix<int> argmin_i(S, N, -1);
  spop::Matrix<int> argmin_s(S, N, -1);

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
      //Rcpp::Rcout << "Enter cost at time: " << t << "\n";
      double current_MIN = std::numeric_limits<double>::infinity();
      double best_speed;
      int best_i = -1;
      int best_s = -1;

      for (int s = 0; s < t; s++){   // previous time
        //Rcpp::Rcout << "Enter loop of previous times"<< "\n";
        //Rcpp::Rcout << "======================================================="<< "\n";
        for (int i = 0; i < S; i++){ // previous state
          // Fix start and end position in time and space
          double p_s = states(i, s);
          double p_t = states(j, t);
          if (s == 0){
            //Rcpp::Rcout << "Enter loop of NO previous changes: s=0"<< "\n";
            for (int j = 0; j < nb_initSpeed; j++){ // init speed loop
              double v_t = 2*(p_t - p_s)/(t - s) - initspeeds(j,0); // simple slope rule
              // Quadratic cost for interval [s, t)
              double interval_cost = qc.quadratic_cost_interval(s, t, p_s, p_t, v_t);
              // Candidate cost (DP recurrence)
              double candidate = interval_cost;
              //Rcpp::Rcout << "s: " << s <<" t: "<<t << " ps: "<<p_s << " pt: "<< p_t <<" v_t: "<<v_t<<"\n";
              //Rcpp::Rcout << "Cost of no change: "<< candidate;
              
              if (candidate < current_MIN){
                current_MIN = candidate;
                best_speed = v_t;
                best_i = i;
                best_s = s;
              }
            }
          }
          else{
           //Rcpp::Rcout << "Enter loop of WITH previous changes: s!=0"<< "\n";
            //Rcpp::Rcout << "Previous time: " << s << " current time: "<< t<<"\n";
            // simple slope rule
            double v_t = 2*(p_t - p_s)/(t - s) - speeds(i, s); 
            
            // Quadratic cost for interval [s, t)
            //Rcpp::Rcout << "s: " << s <<" t: "<<t << " ps: "<<p_s << " pt: "<< p_t <<" v_t: "<<v_t<<"\n";
            double interval_cost = qc.quadratic_cost_interval(s, t, p_s, p_t, v_t);
            
            // Candidate cost (DP recurrence)
            double candidate = value(i, s) + interval_cost + beta;
            //Rcpp::Rcout << "vt = " << v_t << " int cost: " << interval_cost << "candidate " << candidate<< "\n";
            if (candidate < current_MIN){
              current_MIN = candidate;
              best_speed = v_t;
              best_i = i;
              best_s = s;
            }
          }
        }
      }
      if (std::fpclassify(current_MIN) == FP_SUBNORMAL) {
        current_MIN = 0.0;
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
  //Rcpp::Rcout << " ======================================= \n";
  //Rcpp::Rcout << " ======================================= \n";
  //Rcpp::Rcout <<"Change points: " << change_points;
  //Rcpp::Rcout << " ======================================= \n";
  //Rcpp::Rcout << " ======================================= \n";
  return value;
}


