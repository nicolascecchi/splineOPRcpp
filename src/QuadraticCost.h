#ifndef QUADRATIC_COST_H
#define QUADRATIC_COST_H

#include <vector>

// Class for fast O(1) evaluation of quadratic spline cost on equispaced intervals
class QuadraticCost
{
  private:
    std::vector<double> y;
    int N;
    // Precomputed cumulative sums depending on 
    std::vector<double> cumsum_y;
    std::vector<double> cumsum_y2;
    std::vector<double> cumsum_yL1;
    std::vector<double> cumsum_yL2;
    // Polynomial coefficients and interval length
    //double L;
    //double a;
    //double b;
    //double c;
    double sum_y;
    double sum_y2;
    double sum_yL1;
    double sum_yL2;

  public:
    // getters
    std::vector<double> get_cumsum_y()   const {return cumsum_y;}
    std::vector<double> get_cumsum_y2()  const {return cumsum_y2;}
    std::vector<double> get_cumsum_yL1() const {return cumsum_yL1;}
    std::vector<double> get_cumsum_yL2() const {return cumsum_yL2;}
    //double get_L() const {return L;}
    //double get_a() const {return a;}
    //double get_b() const {return b;}
    //double get_c() const {return c;}
    double get_sum_y()   const {return sum_y;}
    double get_sum_y2()  const {return sum_y2;}
    double get_sum_yL1() const {return sum_yL1;}
    double get_sum_yL2() const {return sum_yL2;}
    
    // Constructor
    explicit QuadraticCost(std::vector<double> data);

    // Compute quadratic cost on interval [s, t) for given (p_s, p_t, v_t)
    double interval_cost(int s, int t, double p_s, double p_t, double v_t) const;
};

#endif // QUADRATIC_COST_H
