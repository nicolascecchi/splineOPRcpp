#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <cstddef>

class Matrix
{
  size_t rows;
  size_t cols;
  std::vector<double> data;

public:
  Matrix(size_t r, size_t c, double init = 0.0);

  double& operator()(size_t i, size_t j);
  const double& operator()(size_t i, size_t j) const;

  size_t nrows() const;
  size_t ncols() const;
};

#endif
