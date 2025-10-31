#ifndef SPLINEOP_H
#define SPLINEOP_H

#include <Rcpp.h> // For binding with R
#include <iostream> // WHY?
#include <algorithm> // USE IT BUT CHEKC WHY
#include <vector>   // Store changepoints, and other (?)
#include <limits>   // Use it to have the infinity for the CURRENT_MIN cost
#include <cmath>   // use it in Faulhaber inside the cost 
#include <random>      // for random number generation
#include "Matrix.h" // for storing values in DP loo
#include "QuadraticCost.h" // for computing the cost 
#include <optional> // To avoid initializing matrices at construction

class SplineOP
{
    // private methods and attributes
    private:
        int nobs; // Size of the data --> Change to npoints for being more explicit?
        size_t nstates;
        size_t nspeeds;
        //QuadraticCost qc; // Cost object to compute intervals. Needs DATA to precompute stuff. 

        spop::Matrix<double> speeds; // best speed holder
        spop::Matrix<double> costs;  // matrix of costs
        spop::Matrix<double> initspeeds; //  to check if matrix or vector ; set of initial speeds;
        spop::Matrix<double> states; // sets of state for each time 
        
        spop::Matrix<int> argmin_i; // Store INDEX of best previous state
        spop::Matrix<int> argmin_s; // Store best previous time
        
        QuadraticCost qc;
        
        std::vector<int> changepoints;
        spop::Matrix<double> generate_states(size_t nstates, const std::vector<double>& data, double data_var, int seed); // State generator
        void backtrack_changes();
    //public methods and attributes
    public:
        //void set_speeds(const spop::Matrix<double>& speeds);
        void predict(double beta); // predicts with a given penalty

        //// Getters
        std::vector<int> get_changepoints() const {return changepoints;}
        spop::Matrix<double> get_speeds() const {return speeds;}
        spop::Matrix<double> get_costs() const {return costs;}
        spop::Matrix<double> get_initspeeds() const {return initspeeds;}
        spop::Matrix<double> get_states() const {return states;}
        spop::Matrix<int> get_argmin_i() const {return argmin_i;}
        spop::Matrix<int> get_argmin_s() const {return argmin_s;}
    //constructor
    explicit SplineOP(std::vector<double>  data
                      ,size_t nstates
                      ,size_t nspeeds
                      ,double data_var
                      ,int seed);

};

#endif
