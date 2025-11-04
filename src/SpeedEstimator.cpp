#include <Rcpp.h>
#include <RcppEigen.h> 
#include <vector>
#include <stdexcept>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * @brief Computes the least squares quadratic approximation for one dimension.
 * @param data_row The single dimension data (N x 1 vector) passed by CONST REFERENCE.
 */
double get_approx_speed_1d(const Eigen::VectorXd& data_row, size_t N) {
    
    // ... (Input Validation and Setup Matrices remain the same) ...
    
    // 1. Input Validation
    if (N < 3 || N > static_cast<size_t>(data_row.size())) {
        
        throw std::runtime_error("N (endpoint) must be at least 3 and not exceed the data size.");
    }
    
    // 2. Setup Matrices for LSQ (Y = X * Beta)
    Eigen::MatrixXd X(N, 3);
    Eigen::VectorXd Y(N);
    
    for (int i = 0; i < N; ++i) {
        double x = static_cast<double>(i);
        
        // Design Matrix X (independent of data)
        X(i, 0) = x * x; 
        X(i, 1) = x;     
        X(i, 2) = 1.0;   
        
        // Response Vector Y (uses the referenced data_row)
        Y(i) = data_row(i);
    }
    
    // 3. Solve the Least Squares System: Beta = (X' * X)^-1 * X' * Y
    Eigen::VectorXd Beta = X.householderQr().solve(Y);
    
    // 4. Return the coefficient 'b'
    return Beta(1);
}


/**
 * @brief Calculates the initial speed for each dimension using multiple window sizes.
 * @param data The input data matrix (Rows = Dimensions, Columns = Observations).
 * @param sizes A vector of window lengths (N) to use for the LSQ fit.
 * @return Eigen::MatrixXd A matrix (Dimensions x Sizes) containing the estimated speeds.
 */
Eigen::MatrixXd EstimateSpeeds(Eigen::MatrixXd data, const std::vector<int> sizes){    
    Eigen::Index ndims = data.rows(); // size of output speed vector
    Eigen::Index nsamples_to_estimate = sizes.size(); // quantity of speed estimations
    Eigen::MatrixXd initSpeeds(ndims, nsamples_to_estimate); // Place holder for the different initSpeeds
    
    // Loop through each dimension (row)
    for (Eigen::Index j = 0; j < ndims; ++j) 
    {
        // get current row
        const Eigen::VectorXd data_row = data.row(j).transpose(); 
        
        // Loop through each window size 
        for (Eigen::Index k = 0; k < nsamples_to_estimate; ++k)
        {
            size_t N_window = static_cast<size_t>(sizes[k]);
            
            // Pass the materialized data_row by CONST REFERENCE
            initSpeeds(j, k) = get_approx_speed_1d(data_row, N_window);
        }
    }
    
    return initSpeeds;
}
