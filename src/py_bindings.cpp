#include <pybind11/pybind11.h>
// Include for automatic conversion of Eigen matrices to NumPy arrays
#include <pybind11/eigen.h> 
// Include for automatic conversion of standard C++ containers (like std::vector) to Python lists
#include <pybind11/stl.h>   

// Include your C++ headers
#include "SplineOP.h"
#include "SplineOP_constrained.h"
#include "QuadraticCost.h"
#include "Faulhaber.h"
namespace py = pybind11;

// Define the Python module named 'splineop_py'
// The macro takes the module name and a variable 'm' for the module object
PYBIND11_MODULE(splineop_cpp, m) {
    m.doc() = "Python bindings for SplineOP, a change point detection library using optimal partitioning.";
    
    // Expose Faulhaber function for debugging and compare with Python version.
    m.def("faulhaber", &Faulhaber, 
          "Calculates the sum of the first n integers raised to the power of deg.",
          py::arg("n"), py::arg("deg"));
    // =================================================================
    // 1. Expose SplineOP (Unconstrained Optimal Partitioning)
    // =================================================================
    py::class_<SplineOP>(m, "SplineOP")
        // Constructor: .constructor<...>() becomes .def(py::init<...>())
        .def(py::init<
             Eigen::MatrixXd,   // data
             size_t,            // nstates
             std::vector<int>,  // initial speed size for estimator
             double,            // data_var
             int>(),            // seed
             py::arg("data"),
             py::arg("nstates"),
             py::arg("initial_speed_size"),
             py::arg("data_var"),
             py::arg("seed")
             )

        // Read-Only Properties (Getters)
        // .property("R_name", &Class::get_method) becomes .def_property_readonly("py_name", &Class::get_method)
        .def_property_readonly("changepoints", &SplineOP::get_changepoints) 
        .def_property_readonly("speeds", &SplineOP::get_speeds) 
        .def_property_readonly("costs", &SplineOP::get_costs) 
        .def_property_readonly("init_speeds", &SplineOP::get_initSpeeds) 
        .def_property_readonly("states", &SplineOP::get_states) 
        .def_property_readonly("argmin_i", &SplineOP::get_argmin_i) 
        .def_property_readonly("argmin_s", &SplineOP::get_argmin_s)

        // Cumulative Sum Getters
        .def_property_readonly("cumsum_y", &SplineOP::get_cumsum_y)
        .def_property_readonly("cumsum_y2", &SplineOP::get_cumsum_y2)
        .def_property_readonly("cumsum_yL1", &SplineOP::get_cumsum_yL1)
        .def_property_readonly("cumsum_yL2", &SplineOP::get_cumsum_yL2)
        .def_property_readonly("sum_y", &SplineOP::get_sum_y)
        .def_property_readonly("sum_y2", &SplineOP::get_sum_y2)
        .def_property_readonly("sum_yL1", &SplineOP::get_sum_yL1)
        .def_property_readonly("sum_yL2", &SplineOP::get_sum_yL2)

        // Methods
        .def("set_qc", &SplineOP::set_qc)
        .def("get_segment_cost", &SplineOP::get_segment_cost) 
        .def("predict", &SplineOP::predict, "Predicts changepoints with given penalty.")
        .def("pruningv1", &SplineOP::pruningv1, "Predicts changepoints with given penalty and pruning.")
        .def("pruningv2", &SplineOP::pruningv2, "Predicts changepoints with given penalty and pruning.")
        .def("get_pruning_costs", &SplineOP::get_pruning_costs, "Should be useless at the end, all +inf.")
        .def("get_non_pruned_times", &SplineOP::get_non_pruned_times, "Non pruned times at the end.")
        .def("set_init_speeds", &SplineOP::set_initSpeeds)
        .def("set_states", &SplineOP::set_states);
        
    // =================================================================
    // 2. Expose SplineOP_constrained (Fixed-K Optimal Partitioning)
    // =================================================================
    py::class_<SplineOP_constrained>(m, "SplineOP_constrained")
        .def(py::init<
             Eigen::MatrixXd, // data
             size_t, // nstates
             std::vector<int>, // initial speed size for estimator
             double,  // data_var
             int, // K
             int>(), // seed
             py::arg("data"),
             py::arg("nstates"),
             py::arg("initial_speed_size"),
             py::arg("data_var"),
             py::arg("K"),
             py::arg("seed")
             )
        .def_property_readonly("changepoints", &SplineOP_constrained::get_changepoints) 
        .def_property_readonly("init_speeds", &SplineOP_constrained::get_initSpeeds) 
        .def_property_readonly("states", &SplineOP_constrained::get_states) 

        // Cumulative Sum Getters
        .def_property_readonly("cumsum_y", &SplineOP_constrained::get_cumsum_y)
        .def_property_readonly("cumsum_y2", &SplineOP_constrained::get_cumsum_y2)
        .def_property_readonly("cumsum_yL1", &SplineOP_constrained::get_cumsum_yL1)
        .def_property_readonly("cumsum_yL2", &SplineOP_constrained::get_cumsum_yL2)
        .def_property_readonly("sum_y", &SplineOP_constrained::get_sum_y)
        .def_property_readonly("sum_y2", &SplineOP_constrained::get_sum_y2)
        .def_property_readonly("sum_yL1", &SplineOP_constrained::get_sum_yL1)
        .def_property_readonly("sum_yL2", &SplineOP_constrained::get_sum_yL2)

        // Methods
        .def("set_qc", &SplineOP_constrained::set_qc)
        .def("get_segment_cost", &SplineOP_constrained::get_segment_cost) 
        .def("predict", &SplineOP_constrained::predict, "Predicts changepoints with K changepoints.")
        .def("backtrack_changes", &SplineOP_constrained::backtrack_changes, "Backtracks changepoints from a model with computed costs.")
        .def("get_costs", &SplineOP_constrained::get_costs, "Gets the slice K of the get_costs tensor.")
        .def("get_argmin_i", &SplineOP_constrained::get_argmin_i, "Gets the slice K of the get_argmin_i tensor.")
        .def("get_argmin_s", &SplineOP_constrained::get_argmin_s, "Gets the slice K of the get_argmin_s tensor.")
        .def("set_init_speeds", &SplineOP_constrained::set_initSpeeds)
        .def("set_states", &SplineOP_constrained::set_states);
    

    // =================================================================
    // 3. Expose QuadraticCost Class
    // =================================================================

    py::class_<QuadraticCost>(m, "QuadraticCost")
        .def(py::init<Eigen::MatrixXd>(), py::arg("data"))
        
        // Properties
        .def_property_readonly("cumsum_y", &QuadraticCost::get_cumsum_y) 
        .def_property_readonly("cumsum_y2", &QuadraticCost::get_cumsum_y2)
        .def_property_readonly("cumsum_yL1", &QuadraticCost::get_cumsum_yL1)
        .def_property_readonly("cumsum_yL2", &QuadraticCost::get_cumsum_yL2)
        .def_property_readonly("nobs", &QuadraticCost::get_nobs)
        .def_property_readonly("ndims", &QuadraticCost::get_ndims)
        .def_property_readonly("data", &QuadraticCost::get_data)
        
        // Method
        .def("interval_cost", &QuadraticCost::interval_cost, 
            "Computes the cost of the interval, given border parameters.");
            
}
