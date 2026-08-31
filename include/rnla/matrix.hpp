#pragma once

#include <cstddef>
#include <vector>

namespace rnla {

class Matrix {
 public:
  Matrix() = default;
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
  int rows_ = 0;
  int cols_ = 0;
  std::vector<double> data_;
};

// C = op(A) * op(B), where op(X) is X or X^T.
Matrix matmul(const Matrix& A, const Matrix& B, bool transA = false, bool transB = false);

// Copies src into dst starting at column col_offset.
void set_block(Matrix& dst, const Matrix& src, int col_offset);

}  // namespace rnla