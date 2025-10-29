#ifndef RCPP_WRAP_SPECIALIZATION_H
#define RCPP_WRAP_SPECIALIZATION_H

#include <Rcpp.h>
#include "splineOP.h" // Includes your spop::Matrix<T> definition

// Make sure the access operator is defined for your custom matrix.
// Assuming your spop::Matrix<T> class has a public member function:
// T operator()(size_t i, size_t j) const; 

namespace Rcpp {

template <typename T>
inline SEXP wrap(const spop::Matrix<T>& matrix) {
    
    // 1. *** CORRECTION ***
    //    Get the RTYPE integer constant (e.g., REALSXP for double) 
    //    using the correct Rcpp trait 'r_sexptype_traits'
    //    and define RType (e.g., Rcpp::Matrix<REALSXP>) in one step.
    typedef typename Rcpp::Matrix<Rcpp::traits::r_sexptype_traits<T>::rtype> RType;
    
    // Get dimensions using the custom getters
    size_t rows = matrix.nrows();
    size_t cols = matrix.ncols();

    // 3. Create the Rcpp Matrix object
    // This should now work as 'RType' is correctly defined.
    RType wrapped_matrix((int)rows, (int)cols); 

    // 4. Manually copy data 
    for (size_t j = 0; j < cols; ++j) {
        for (size_t i = 0; i < rows; ++i) {
            // *** Using the parenthesis access operator matrix(i, j) ***
            wrapped_matrix(i, j) = matrix(i, j); 
        }
    }

    // 5. Return the Rcpp Matrix object (as SEXP)
    return wrapped_matrix;
}

} // end namespace Rcpp

#endif // RCPP_WRAP_SPECIALIZATION_H

