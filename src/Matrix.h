#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <cstddef>
#include <stdexcept>

namespace spop{
template <typename T>
class Matrix
{
  private:
    size_t rows;
    size_t cols;
    std::vector<T> data;
    //std::vector<int> change_points;

  public:
    // Constructors
    // 1. Default constructor
    Matrix(); 
    //2. The good constructor, the one we like
    Matrix(size_t r, size_t c, T init = T());
    // 3. Copy Constructor (Deep copy: Needed if you pass Matrix objects around)
    Matrix(const Matrix& other);

    // Element access
    T& operator()(size_t i, size_t j);
    const T& operator()(size_t i, size_t j) const;

    // Resize method
    //void resize(size_t nrows, size_t ncols, T initial_value);
    //std::vector<int> get_changes() const {return change_points;}

    // Dimensions
    size_t nrows() const;
    size_t ncols() const;
};

// Implementation (must be in header for templates)

// Constructor
template <typename T>
Matrix<T>::Matrix(size_t r, size_t c, T init) : rows(r), cols(c), data(r * c, init) {}
template <typename T>
Matrix<T>::Matrix() : rows(0), cols(0), data(0) {}
template <typename T>
Matrix<T>::Matrix(const Matrix<T>& other)
    : rows(other.rows),
      cols(other.cols),
      data(other.data) {} // body can be empty, everything is already handled

// Element access
template <typename T>
T& Matrix<T>::operator()(size_t i, size_t j)
{
  //if (i >= rows || j >= cols)
  //  throw std::out_of_range("Matrix index out of range");
  return data[i * cols + j];
}

template <typename T>
const T& Matrix<T>::operator()(size_t i, size_t j) const
{
  //if (i >= rows || j >= cols)
  //  throw std::out_of_range("Matrix index out of range");
  return data[i * cols + j];
}

// Dimensions
template <typename T>
size_t Matrix<T>::nrows() const { return rows; }

template <typename T>
size_t Matrix<T>::ncols() const { return cols; }

}

#endif // MATRIX_H
