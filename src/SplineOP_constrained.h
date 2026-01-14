#ifndef SPLINEOPCONSTRAINED_H
#define SPLINEOPCONSTRAINED_H

#include <iostream> // WHY?
#include <algorithm> // USE IT BUT CHEKC WHY
#include <vector>   // Store changepoints, and other (?)
#include <limits>   // Use it to have the infinity for the CURRENT_MIN cost
#include <cmath>   // use it in Faulhaber inside the cost 
#include <random>      // for random number generation
#include <optional> // To avoid initializing matrices at construction
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor> // 3D Arrays

#include "QuadraticCost.h" // for computing the cost 


class SplineOP_constrained
{
        //public methods and attributes
    public:
    //constructor
        explicit SplineOP_constrained(Eigen::MatrixXd data
                      ,int nstates
                      ,std::vector<int> sp
                      ,double data_var
                      ,int K
                      ,int seed);
        
        //void set_speeds(const Eigen::MatrixXdXd& speeds);
        void predict(int K); // predicts with a given penalty
        void backtrack_changes(int K);
        
        //setters 
        void set_qc(Eigen::MatrixXd& data);
        void set_states(std::vector<Eigen::MatrixXd> new_states);// {states = new_states;}
        void set_initSpeeds(Eigen::MatrixXd new_initSpeeds);// {initSpeeds = new_initSpeeds;}

        //// Getters
        std::vector<int> get_changepoints() const {return changepoints;}
        std::vector<Eigen::MatrixXd> get_states() const {return states;}
        double get_segment_cost(int s, int t, Eigen::VectorXd p_s, Eigen::VectorXd p_t, Eigen::VectorXd v_s); // compute a cost with given parameters
        Eigen::MatrixXd get_initSpeeds() const {return initSpeeds;}
       
        // More complex getters
        //Eigen::MatrixXd get_speeds(int k) const {return speeds;}
        Eigen::MatrixXd get_costs(int k) const;
        Eigen::MatrixXi get_argmin_i(int k) const;
        Eigen::MatrixXi get_argmin_s(int k) const;        
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
        int K;
        std::vector<int> sp;
        int nstates;
        int nspeeds;
        //QuadraticCost qc; // Cost object to compute intervals. Needs DATA to precompute stuff. 
        
        Eigen::Tensor<double, 4> speeds; // best speed holder
        Eigen::Tensor<double, 3> costs;  // K+1 matrices of costs
        Eigen::MatrixXd initSpeeds; //  to check if matrix or vector ; set of initial speeds;
        std::vector<Eigen::MatrixXd> states; // sets of state for each time 
        
        Eigen::Tensor<int,3> argmin_i; // Store INDEX of best previous state
        Eigen::Tensor<int,3> argmin_s; // Store best previous time
        
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
        
        void compute_best_with_change(int K, int t, int j);
        void compute_best_without_change( Eigen::VectorXd& p_t, double& current_MIN,Eigen::VectorXd& best_out_speed,int& best_i,int& best_s, int& t);

};

#endif
