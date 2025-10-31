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
    std::vector<T> data_; // call data_ to avoid collision with std::vector<T>.data method
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
    

    // Get Dimensions
    size_t nrows() const;
    size_t ncols() const;
    //std::vector<T> get_data() const;

    // For Rcpp::wrap (C++ -> R): Returns a CONST pointer for reading
    // For Rcpp::wrap (C++ -> R): Returns a CONST pointer for reading
    const T* get_data() const {return data_.data();}
    // std::vector::data() returns a pointer to the vector's underlying array
    T* get_data_pointer() {return data_.data();}
    // This allows writing data back into the matrix's internal vector
    
    
};

// Implementation (must be in header for templates)

// Constructor
template <typename T>
Matrix<T>::Matrix(size_t r, size_t c, T init) : rows(r), cols(c), data_(r * c, init) {}
template <typename T>
Matrix<T>::Matrix() : rows(0), cols(0), data_(0) {}
template <typename T>
Matrix<T>::Matrix(const Matrix<T>& other)
    : rows(other.rows),
      cols(other.cols),
      data_(other.data_) {} // body can be empty, everything is already handled

// Element access
// Column-major access to interact easily with R
template <typename T>
T& Matrix<T>::operator()(size_t i, size_t j) { return data_[j * rows + i]; }

template <typename T>
const T& Matrix<T>::operator()(size_t i, size_t j) const { return data_[j * rows + i];}

// Getters to member variables
template <typename T>
size_t Matrix<T>::nrows() const { return rows; }

template <typename T>
size_t Matrix<T>::ncols() const { return cols; }

//template <typename T>
//std::vector<T> Matrix<T>::get_data() const { return data; }


}




#endif 
