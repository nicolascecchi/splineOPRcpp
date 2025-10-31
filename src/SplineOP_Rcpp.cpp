#include <RcppCommon.h>
#include <Rcpp.h>
#include "SplineOP.h"
//#include "Matrix.h"

//====================================================================
// This is what exposes your C++ class to R.
//====================================================================
//RCPP_EXPOSED_AS(spop::Matrix<double>)
//RCPP_EXPOSED_AS(spop::Matrix<int>)

//RCPP_EXPOSED_CLASS(std::vector<double>);

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
        //.property("get_speeds", &SplineOP::get_speeds) 
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

