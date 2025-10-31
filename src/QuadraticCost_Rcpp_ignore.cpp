#include <Rcpp.h>
#include "QuadraticCost.h"
#include "Matrix.h"


//====================================================================
// This is what exposes your C++ class to R.
//====================================================================
RCPP_MODULE(quadraticcost_module) {
    // This line exposes the C++ class "myApp::Accumulator"
    // as an R class named "Accumulator"
    Rcpp::class_<QuadraticCost>("QuadraticCost")

    // 1. Expose the constructor
        .constructor<std::vector<double>>("data")

    // 2. Expose the methods
    // Format: .method("R_name", &C++_class::C++_method, "Optional docstring")
        .property("get_cumsum_y", &QuadraticCost::get_cumsum_y) 
        .property("get_cumsum_y2", &QuadraticCost::get_cumsum_y2)
        .property("get_cumsum_yL1", &QuadraticCost::get_cumsum_yL1)
        .property("get_cumsum_yL2", &QuadraticCost::get_cumsum_yL2)
       // .property("get_", &QuadraticCost::get_)
        .property("get_a", &QuadraticCost::get_a)
        .property("get_b", &QuadraticCost::get_b)
        .property("get_c", &QuadraticCost::get_c)
        .property("get_L", &QuadraticCost::get_L)
        
        .property("get_sum_y", &QuadraticCost::get_sum_y)
        .property("get_sum_y2", &QuadraticCost::get_sum_y2)
        .property("get_sum_yL1", &QuadraticCost::get_sum_yL1)
        .property("get_sum_yL2", &QuadraticCost::get_sum_yL2)

        .method("segmentcost", &QuadraticCost::quadratic_cost_interval, 
            "Adds a number to the accumulator");
            
    }

