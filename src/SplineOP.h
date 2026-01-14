#ifndef SPLINEOP_H
#define SPLINEOP_H

//#include <Rcpp.h> // For binding with R
#include <iostream> // WHY?
#include <algorithm> // USE IT BUT CHEKC WHY
#include <vector>   // Store changepoints, and other (?)
#include <limits>   // Use it to have the infinity for the CURRENT_MIN cost
#include <cmath>   // use it in Faulhaber inside the cost 
#include <random>      // for random number generation
//#include <RcppEigen.h> // for storing values in DP loo
#include <Eigen/Dense>
#include <optional> // To avoid initializing matrices at construction

#include "QuadraticCost.h" // for computing the cost 


class SplineOP
{
    //public methods and attributes
    public:
        //constructor
        explicit SplineOP(Eigen::MatrixXd data
                      ,int nstates
                      ,std::vector<int> sp
                      ,double data_var
                      ,int seed);    
        
        void set_qc(Eigen::MatrixXd& data);
    
        
        //void set_speeds(const Eigen::MatrixXdXd& speeds);
        void predict(double beta); // predicts with a given penalty
        
        void pruningv1(double beta, double margin); // predicts with a given penalty
        void pruningv2(double beta); // predicts with a given penalty


        //// Getters
        std::vector<int> get_changepoints() const {return changepoints;}
        std::vector<Eigen::MatrixXd> get_speeds() const {return speeds;}
        Eigen::MatrixXd get_costs() const {return costs;}
        Eigen::MatrixXd get_initSpeeds() const {return initSpeeds;}
        std::vector<Eigen::MatrixXd> get_states() const {return states;}
        Eigen::MatrixXi get_argmin_i() const {return argmin_i;}
        Eigen::MatrixXi get_argmin_s() const {return argmin_s;}
        double get_segment_cost(int s, int t, Eigen::VectorXd p_s, Eigen::VectorXd p_t, Eigen::VectorXd v_s); // compute a cost with given parameters
        Eigen::MatrixXd get_pruning_costs() const {return pruning_costs;}
        std::vector<std::vector<int>> get_non_pruned_times() const {return times_for_states;}
        //setters 
        void set_states(std::vector<Eigen::MatrixXd> new_states);// {states = new_states;}
        void set_initSpeeds(Eigen::MatrixXd new_initSpeeds);// {initSpeeds = new_initSpeeds;}


        // Cost getters
        Eigen::MatrixXd get_cumsum_y()   const {return qc.get_cumsum_y();}
        Eigen::MatrixXd get_cumsum_y2()  const {return qc.get_cumsum_y2();}
        Eigen::MatrixXd get_cumsum_yL1() const {return qc.get_cumsum_yL1();}
        Eigen::MatrixXd get_cumsum_yL2() const {return qc.get_cumsum_yL2();}
        Eigen::VectorXd get_sum_y()   const {return qc.get_sum_y();}
        Eigen::VectorXd get_sum_y2()  const {return qc.get_sum_y2();}
        Eigen::VectorXd get_sum_yL1() const {return qc.get_sum_yL1();}
        Eigen::VectorXd get_sum_yL2() const {return qc.get_sum_yL2();}

 
    // private methods and attributes
    private:
        int nobs; // Size of the data --> Change to npoints for being more explicit?
        int ndims;
        std::vector<int> sp;
        int nstates;
        int nspeeds;
        //QuadraticCost qc; // Cost object to compute intervals. Needs DATA to precompute stuff. 
        
        std::vector<Eigen::MatrixXd> speeds; // best speed holder
        Eigen::MatrixXd costs;  // matrix of costs
        Eigen::MatrixXd pruning_costs;  // matrix of costs
        std::vector<std::vector<int>> times_for_states;
        Eigen::MatrixXd initSpeeds; //  to check if matrix or vector ; set of initial speeds;
        std::vector<Eigen::MatrixXd> states; // sets of state for each time 
        
        Eigen::MatrixXi argmin_i; // Store INDEX of best previous state
        Eigen::MatrixXi argmin_s; // Store best previous time
        
        QuadraticCost qc;
        
        std::vector<int> changepoints;
        std::vector<Eigen::MatrixXd> generate_states(
                                            int nstates,
                                            Eigen::MatrixXd data, // Input is now MatrixXd
                                            double data_var, 
                                            int seed);
        Eigen::MatrixXd generate_matrix_of_noise( //MON
                                            std::mt19937& gen, 
                                            double std_dev, 
                                            int rows, 
                                            int cols);
        void backtrack_changes();
        void prunev1(int t, double margin);
        void prunev2(int t);
};

#endif
