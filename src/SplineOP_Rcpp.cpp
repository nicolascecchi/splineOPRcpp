#include <RcppCommon.h>
#include <Rcpp.h>
#include "Matrix.h"
#include "SplineOP.h"

//====================================================================
// This is what exposes your C++ class to R.
//====================================================================

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
    
// ... continue with Rcpp::as specializations below ...

RCPP_MODULE(splineop) {
    // This line exposes the C++ class "myApp::Accumulator"
    // as an R class named "Accumulator"
    Rcpp::class_<SplineOP>("SplineOP")
    // 1. Expose the constructor
        .constructor<
                std::vector<double>  // data (Must match const reference)
                ,size_t                        // nstates
                ,size_t                        // nspeeds
                ,double                     // data_var
                ,int>()                       // seed
    // 2. Expose the methods
    // Format: .method("R_name", &C++_class::C++_method, "Optional docstring")
    // all of these require the exposition of Matrix class to R     
        .property("get_changepoints", &SplineOP::get_changepoints) 
        .property("get_speeds", &SplineOP::get_speeds) 
        //.property("get_costs", &SplineOP::get_costs) 
        //.property("get_initspeeds", &SplineOP::get_initspeeds) 
        //.property("get_states", &SplineOP::get_states) 
        //.property("get_argmin_i", &SplineOP::get_argmin_i) 
        //.property("get_argmin_s", &SplineOP::get_argmin_s) 

        .method("predict", &SplineOP::predict, 
            "Predicts changepoints with given penalty");
        
    // EXPOSE QUADRATIC COST CLASS TO R

    Rcpp::class_<QuadraticCost>("QuadraticCost")
        .constructor<std::vector<double> >("data")
        // Format: .method("R_name", &C++_class::C++_method, "Optional docstring")
        .property("get_cumsum_y", &QuadraticCost::get_cumsum_y) 
        .property("get_cumsum_y2", &QuadraticCost::get_cumsum_y2)
        .property("get_cumsum_yL1", &QuadraticCost::get_cumsum_yL1)
        .property("get_cumsum_yL2", &QuadraticCost::get_cumsum_yL2)
        // .property("get_", &QuadraticCost::get_)
        
        .method("segmentcost", &QuadraticCost::quadratic_cost_interval, 
            "Adds a number to the accumulator");
    
    //Rcpp::class_<spop::Matrix<double>>("MatrixDouble")
    //    .constructor<size_t, size_t, double>() 
    //    //expose other methods
    //    ; 
    //Rcpp::class_<spop::Matrix<int>>("MatrixInt")
    //    .constructor<size_t, size_t, int>()
    //    //expose other methods
    //    ;
            
    }

