#ifndef DP_MATRIX_H
#define DP_MATRIX_H

#include <vector>

// Pure C++ dynamic programming function
std::vector<int> dp_matrix(const std::vector<double>& data,
                           double beta,
                           int S = 10,
                           int nb_initSpeed = 5,
                           double data_var = 1.0);

#endif // DP_MATRIX_H
