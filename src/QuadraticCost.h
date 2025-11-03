#ifndef QUADRATIC_COST_H
#define QUADRATIC_COST_H

#include <vector>
#include <RcppEigen.h>
//#include <RcppEigen.h> 

// Class for fast O(1) evaluation of quadratic spline cost on equispaced intervals
class QuadraticCost
{
  private:
    Eigen::MatrixXd y;
    Eigen::Index nobs;
    Eigen::Index ndims;
    // Precomputed cumulative sums depending on 
    Eigen::MatrixXd cumsum_y;
    Eigen::MatrixXd cumsum_y2;
    Eigen::MatrixXd cumsum_yL1;
    Eigen::MatrixXd cumsum_yL2;
    // Polynomial coefficients and interval length
    Eigen::VectorXd sum_y;
    Eigen::VectorXd sum_y2;
    Eigen::VectorXd sum_yL1;
    Eigen::VectorXd sum_yL2;

  public:
    // getters
    Eigen::MatrixXd get_cumsum_y()   const {return cumsum_y;}
    Eigen::MatrixXd get_cumsum_y2()  const {return cumsum_y2;}
    Eigen::MatrixXd get_cumsum_yL1() const {return cumsum_yL1;}
    Eigen::MatrixXd get_cumsum_yL2() const {return cumsum_yL2;}
    Eigen::VectorXd get_sum_y()   const {return sum_y;}
    Eigen::VectorXd get_sum_y2()  const {return sum_y2;}
    Eigen::VectorXd get_sum_yL1() const {return sum_yL1;}
    Eigen::VectorXd get_sum_yL2() const {return sum_yL2;}
    Eigen::Index get_nobs() {return nobs;}
    Eigen::Index get_ndims() {return ndims;}
    Eigen::MatrixXd get_data() {return y;}
    // Constructor
    //explicit QuadraticCost(Eigen::VectorXd data);
    explicit QuadraticCost(Eigen::MatrixXd data);

    // Compute quadratic cost on interval [s, t) for given (p_s, p_t, v_t)
    //double interval_cost(int s, int t, double p_s, double p_t, double v_t) const;
    double interval_cost(Eigen::Index s
      , Eigen::Index t
      , Eigen::VectorXd p_s
      , Eigen::VectorXd p_t
      , Eigen::VectorXd v_s);
};

#endif // QUADRATIC_COST_H
