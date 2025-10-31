#ifndef MATRIX_EXPOSURE_H
#define MATRIX_EXPOSURE_H

#include <Rcpp.h>
#include "Matrix.h"

// In SplineOP_Rcpp.cpp
namespace Rcpp {
    
    // Rcpp::wrap for spop::Matrix<double> (C++ -> Rcpp::NumericMatrix)
    template<>
    inline SEXP wrap(const spop::Matrix<double>& cpp_mat) {
        
        size_t rows = cpp_mat.nrows();
        size_t cols = cpp_mat.ncols();
        
        Rcpp::NumericMatrix r_mat(rows, cols);
        
        // Data Copy: Requires cpp_mat.get_data() to return a const double*
        std::copy(cpp_mat.get_data(), 
                  cpp_mat.get_data() + (rows * cols), 
                  r_mat.begin());

        return r_mat;
    }
    
    // Rcpp::wrap for spop::Matrix<int> (C++ -> Rcpp::IntegerMatrix)
    template<>
    inline SEXP wrap(const spop::Matrix<int>& cpp_mat) {
        
        size_t rows = cpp_mat.nrows();
        size_t cols = cpp_mat.ncols();
        
        Rcpp::IntegerMatrix r_mat(rows, cols);
        
        // Data Copy: Requires cpp_mat.get_data() to return a const int*
        std::copy(cpp_mat.get_data(), 
                  cpp_mat.get_data() + (rows * cols), 
                  r_mat.begin());
        
        return r_mat;
    }
    // Rcpp::as for spop::Matrix<int> (R -> C++)
    template<>
    inline spop::Matrix<int> as(SEXP r_mat_sexp) {
        
        Rcpp::IntegerMatrix r_mat(r_mat_sexp);
        
        spop::Matrix<int> cpp_mat(r_mat.nrow(), r_mat.ncol(), 0); 
        
        // Data Copy: Requires cpp_mat.get_data_pointer() to return a writable int*
        std::copy(r_mat.begin(), r_mat.end(), cpp_mat.get_data_pointer());
        
        return cpp_mat;
    }

}
#endif
