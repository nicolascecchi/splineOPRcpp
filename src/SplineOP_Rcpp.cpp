#include <RcppCommon.h>
#include <Rcpp.h>
#include <RcppEigen.h>
#include "SplineOP.h"
#include "SplineOP_constrained.h"
#include "QuadraticCost.h"
//====================================================================
// This is what exposes your C++ class to R.
//====================================================================

RCPP_EXPOSED_CLASS(SplineOP)  
RCPP_EXPOSED_CLASS(SplineOP_constrained)  
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
        .method("predict", &SplineOP::predict, "Predicts changepoints with given penalty.")
        .method("pruningv1", &SplineOP::pruningv1, "Predicts changepoints with given penalty and pruning.")
        .method("pruningv2", &SplineOP::pruningv2, "Predicts changepoints with given penalty and pruning.")
        .method("get_pruning_costs", &SplineOP::get_pruning_costs, "Should be useless at the end, all +inf.")
        .method("get_non_pruned_times", &SplineOP::get_non_pruned_times, "Non pruned times at the end.")
        
        .method("set_initSpeeds", &SplineOP::set_initSpeeds)
        .method("set_states", &SplineOP::set_states);
        
    // Expose constrained version of SplineOP 
    Rcpp::class_<SplineOP_constrained>("SplineOP_constrained")
    // 1. Expose the constructor
        .constructor<
                Eigen::MatrixXd // data (Must match const reference)
                ,size_t // nstates
                ,std::vector<int> // initial speed size for estimator
                ,double  // data_var
                ,int // K
                ,int>() // seed
        .property("get_changepoints", &SplineOP_constrained::get_changepoints) 
        //.property("get_speeds", &SplineOP_constrained::get_speeds) 
        //.property("get_costs", &SplineOP_constrained::get_costs) 
        .property("get_initSpeeds", &SplineOP_constrained::get_initSpeeds) 
        .property("get_states", &SplineOP_constrained::get_states) 
        //.property("get_argmin_i", &SplineOP_constrained::get_argmin_i) 
        //.property("get_argmin_s", &SplineOP_constrained::get_argmin_s)

        .property("get_cumsum_y", &SplineOP_constrained::get_cumsum_y)
        .property("get_cumsum_y2", &SplineOP_constrained::get_cumsum_y2)
        .property("get_cumsum_yL1", &SplineOP_constrained::get_cumsum_yL1)
        .property("get_cumsum_yL2", &SplineOP_constrained::get_cumsum_yL2)
        .property("get_sum_y", &SplineOP_constrained::get_sum_y)
        .property("get_sum_y2", &SplineOP_constrained::get_sum_y2)
        .property("get_sum_yL1", &SplineOP_constrained::get_sum_yL1)
        .property("get_sum_yL2", &SplineOP_constrained::get_sum_yL2)

        .method("set_qc", &SplineOP_constrained::set_qc)
        .method("get_segment_cost", &SplineOP_constrained::get_segment_cost) 
        .method("predict", &SplineOP_constrained::predict, "Predicts changepoints with K changepoints.")
    
        .method("set_initSpeeds", &SplineOP_constrained::set_initSpeeds)
        .method("set_states", &SplineOP_constrained::set_states);
    

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
