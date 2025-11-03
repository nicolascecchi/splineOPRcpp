// SpeedEstimator.h
#ifndef SPEED_ESTIMATOR_H
#define SPEED_ESTIMATOR_H

#include <vector>
#include <RcppEigen.h>

// Declaration of the function to be called by other C++ files
double get_approx_speed_1d(const Eigen::VectorXd& data_row, size_t N);
Eigen::MatrixXd EstimateSpeeds(Eigen::MatrixXd data, const std::vector<int> sizes);   

#endif // SPEED_ESTIMATOR_H
