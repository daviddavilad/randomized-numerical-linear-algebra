#pragma once

#include <cstddef>
#include <vector>

namespace rnla {

class Matrix {
 public:
  Matrix(int rows, int cols)
      : rows_(rows), cols_(cols), data_(static_cast<std::size_t>(rows) * cols) {}

  int rows() const { return rows_; }
  int cols() const { return cols_; }
  int ld() const { return rows_; }

  double& operator()(int i, int j) {
    return data_[static_cast<std::size_t>(j) * rows_ + i];
  }
  double operator()(int i, int j) const {
    return data_[static_cast<std::size_t>(j) * rows_ + i]; 
  }

  double* data() { return data_.data(); }
  const double* data() const { return data_.data(); }

 private:
  int rows_;
  int cols_;
  std::vector<double> data_;
};

// C = op(A) * op(B), where op(X) is X or X^T.
Matrix matmul(const Matrix& A, const Matrix& B, bool transA = false, bool transB = false);

}  // namespace rnla