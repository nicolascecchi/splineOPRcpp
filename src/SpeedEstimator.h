// SpeedEstimator.h
#ifndef SPEED_ESTIMATOR_H
#define SPEED_ESTIMATOR_H

#include <vector>

// Declaration of the function to be called by other C++ files
double get_approx_speed(const std::vector<double>& data, int N);
std::vector<double> EsimateSpeeds(std::vector<double>& data, std::vector<int>sizes);

#endif // SPEED_ESTIMATOR_H
