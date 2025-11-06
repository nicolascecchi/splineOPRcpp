#include <RcppCommon.h>
#include <Rcpp.h>
#include "SplineOP.h"
#include "QuadraticCost.h"
//====================================================================
// This is what exposes your C++ class to R.
//====================================================================

RCPP_EXPOSED_CLASS(SplineOP)    
// ... continue with Rcpp::as specializations below ...
RCPP_MODULE(splineop) {
    // This line exposes the C++ class "myApp::Accumulator"
    // as an R class named "Accumulator"
    Rcpp::class_<SplineOP>("SplineOP")
    // 1. Expose the constructor
        .constructor<
                Eigen::MatrixXd   // data (Must match const reference)
                ,size_t                        // nstates
                ,size_t                        // nspeeds
                ,std::vector<int>           // initial speed size for estimator
                ,double                     // data_var
                ,int>()                       // seed
    // 2. Expose the methods
    // Format: .method("R_name", &C++_class::C++_method, "Optional docstring")
    // all of these require the exposition of Matrix class to R     
        .property("get_changepoints", &SplineOP::get_changepoints) 
        .property("get_speeds", &SplineOP::get_speeds) 
        .property("get_costs", &SplineOP::get_costs) 
        .property("get_initSpeeds", &SplineOP::get_initSpeeds) 
        .property("get_states", &SplineOP::get_states) 
        .property("get_argmin_i", &SplineOP::get_argmin_i) 
        .property("get_argmin_s", &SplineOP::get_argmin_s)

        .property("get_cumsum_y", &SplineOP::get_cumsum_y)
        .property("get_cumsum_y2", &SplineOP::get_cumsum_y2)
        .property("get_cumsum_yL1", &SplineOP::get_cumsum_yL1)
        .property("get_cumsum_yL2", &SplineOP::get_cumsum_yL2)
        .property("get_sum_y", &SplineOP::get_sum_y)
        .property("get_sum_y2", &SplineOP::get_sum_y2)
        .property("get_sum_yL1", &SplineOP::get_sum_yL1)
        .property("get_sum_yL2", &SplineOP::get_sum_yL2)

        .method("set_qc", &SplineOP::set_qc)
        .method("get_segment_cost", &SplineOP::get_segment_cost) 
        .method("predict", &SplineOP::predict, "Predicts changepoints with given penalty");
        
    // EXPOSE QUADRATIC COST CLASS TO R

    Rcpp::class_<QuadraticCost>("QuadraticCost")
        .constructor<Eigen::MatrixXd>("data")
        // Format: .method("R_name", &C++_class::C++_method, "Optional docstring")
        .property("get_cumsum_y", &QuadraticCost::get_cumsum_y) 
        .property("get_cumsum_y2", &QuadraticCost::get_cumsum_y2)
        .property("get_cumsum_yL1", &QuadraticCost::get_cumsum_yL1)
        .property("get_cumsum_yL2", &QuadraticCost::get_cumsum_yL2)
        .property("get_nobs", &QuadraticCost::get_nobs)
        .property("get_ndims", &QuadraticCost::get_ndims)
        .property("get_data", &QuadraticCost::get_data)
        // .property("get_", &QuadraticCost::get_)
        
        .method("interval_cost", &QuadraticCost::interval_cost, 
            "Computes the cost of the intervanl, given border parameters.");
    
            
    }
