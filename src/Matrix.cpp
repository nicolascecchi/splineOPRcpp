#include "Matrix.h"



Matrix::Matrix(size_t r, size_t c, double init) : rows(r), cols(c), data(r * c, init) {}



double& Matrix::operator()(size_t i, size_t j)
{
  return data[i * cols + j];
}


const double& Matrix::operator()(size_t i, size_t j) const
{
  return data[i * cols + j];
}


size_t Matrix::nrows() const
{
  return rows;
}


size_t Matrix::ncols() const
{
  return cols;
}
