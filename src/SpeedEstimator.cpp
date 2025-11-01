#include <Rcpp.h>
#include <RcppEigen.h> 
#include <vector>
#include <stdexcept>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * @brief 
 * 
 * for the first N points of the data and returns the coefficient 'b' (speed).
 */
double get_approx_speed(const std::vector<double>& data, int N) {
    
    // 1. Input Validation
    if (N < 3 || N > data.size()) {
        throw std::runtime_error("N must be at least 3 and not exceed the data size.");
    }
    
    // 2. Setup Matrices for LSQ (Y = X * Beta)
    // X is the design matrix: [ x^2, x, 1 ]
    MatrixXd X(N, 3);
    // Y is the response vector: [ y ]
    VectorXd Y(N);
    
    for (int i = 0; i < N; ++i) {
        double x = (double)i;
        double x_sq = x * x;
        
        // Design Matrix X (column-major order in Eigen is fine)
        X(i, 0) = x_sq; // x^2
        X(i, 1) = x;    // x
        X(i, 2) = 1.0;  // 1
        
        // Response Vector Y (the data value)
        Y(i) = data[i];
    }
    
    // 3. Solve the Least Squares System: Beta = (X' * X)^-1 * X' * Y
    // Use Eigen's QR decomposition for a stable solution
    VectorXd Beta = X.householderQr().solve(Y);
    
    // Beta vector holds the coefficients: [ a, b, c ]
    
    // 4. Return the coefficient 'b' (the linear term, approximation of the first derivative/speed)
    // The 'b' coefficient is at index 1.
    return Beta(1);
}

std::vector<double> EsimateSpeeds(std::vector<double>& data, std::vector<int>sizes){
    std::vector<double> speeds = std::vector<double>(sizes.size(), 0.);
    for (size_t j = 0; j < sizes.size(); ++j){
        speeds[j] = get_approx_speed(data, sizes[j]);
    }
    return speeds;
}
