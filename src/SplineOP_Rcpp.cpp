#include <Rcpp.h>
#include "SplineOP.h"
#include "Matrix.h"

template<> 
struct Rcpp::traits::is_wrapable<spop::Matrix<double>> : public Rcpp::traits::true_type {};

template<> 
struct Rcpp::traits::is_wrapable<spop::Matrix<int>> : public Rcpp::traits::true_type {};

//====================================================================
// This is what exposes your C++ class to R.
//====================================================================
RCPP_MODULE(splineop) {
    // This line exposes the C++ class "myApp::Accumulator"
    // as an R class named "Accumulator"
    Rcpp::class_<SplineOP>("SplineOP")
    // 1. Expose the constructor
        .constructor<
                const std::vector<double>&  // data (Must match const reference)
                ,size_t                        // nstates
                ,size_t                        // nspeeds
                ,double                     // data_var
                ,int>()                       // seed
    // 2. Expose the methods
    // Format: .method("R_name", &C++_class::C++_method, "Optional docstring")
        .property("get_changepoints", &SplineOP::get_changepoints) 
        .property("get_speeds", &SplineOP::get_speeds) 
        .property("get_costs", &SplineOP::get_costs) 
        .property("get_initspeeds", &SplineOP::get_initspeeds) 
        .property("get_states", &SplineOP::get_states) 
        .property("get_argmin_i", &SplineOP::get_argmin_i) 
        .property("get_argmin_s", &SplineOP::get_argmin_s) 

        .method("predict", &SplineOP::predict, 
            "Predicts changepoints with given penalty");
        
    // EXPOSE QUADRATIC COST CLASS TO R

    Rcpp::class_<QuadraticCost>("QuadraticCost")
        .constructor<std::vector<double>>("data")
        // Format: .method("R_name", &C++_class::C++_method, "Optional docstring")
        .property("get_cumsum_y", &QuadraticCost::get_cumsum_y) 
        .property("get_cumsum_y2", &QuadraticCost::get_cumsum_y2)
        .property("get_cumsum_yL1", &QuadraticCost::get_cumsum_yL1)
        .property("get_cumsum_yL2", &QuadraticCost::get_cumsum_yL2)
        // .property("get_", &QuadraticCost::get_)
        // probably to drop these getters 
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
    
    Rcpp::class_<spop::Matrix<double>>("MatrixDouble")
        .constructor<size_t, size_t, double>() 
        //expose other methods
        ; 
        
    Rcpp::class_<spop::Matrix<int>>("MatrixInt")
        .constructor<size_t, size_t, int>()
        //expose other methods
        ;
            
    }

