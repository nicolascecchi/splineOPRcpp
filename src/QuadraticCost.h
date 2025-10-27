#ifndef QUADRATIC_COST_H
#define QUADRATIC_COST_H

#include <vector>

// Class for fast O(1) evaluation of quadratic spline cost on equispaced intervals
class QuadraticCost
{
  private:
    std::vector<double> y;
    int N;

    // Precomputed cumulative sums depending on y
    std::vector<double> cumsum_y;
    std::vector<double> cumsum_y2;
    std::vector<double> cumsum_yL1;
    std::vector<double> cumsum_yL2;

    // Internal helper: Faulhaber sums for equispaced points (0..n-1)
    static double S1(int n);
    static double S2(int n);
    static double S3(int n);
    static double S4(int n);

  public:
    // Constructor
    explicit QuadraticCost(const std::vector<double>& data);

    // Compute quadratic cost on interval [s, t) for given (p_s, p_t, v_t)
    double quadratic_cost_interval(int s, int t, double p_s, double p_t, double v_t) const;
};

#endif // QUADRATIC_COST_H
