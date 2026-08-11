#include <cstdio>

#include "rnla/matrix.hpp"

int main() {
  rnla::Matrix A(2, 3);
  rnla::Matrix B(3, 2);

  for (int j = 0; j < A.cols(); ++j)
    for (int i = 0; i < A.rows(); ++i) A(i, j) = i + 10 * j;

  for (int j = 0; j < B.cols(); ++j)
    for (int i = 0; i < B.rows(); ++i) B(i, j) = i + j;

  rnla::Matrix C = rnla::matmul(A, B);

  for (int i = 0; i < C.rows(); ++i) {
    for (int j = 0; j < C.cols(); ++j) std::printf("%g ", C(i, j));
    std::printf("\n");
  }

  rnla::Matrix D = rnla::matmul(A, A, false, true);  // A * A^T, 2x2

  for (int i = 0; i < D.rows(); ++i) {
    for (int j = 0; j < D.cols(); ++j) std::printf("%g ", D(i, j));
    std::printf("\n");
  }

  return 0;
}