#include <cstdio>

#include "rnla/matrix.hpp"

int main() {
  rnla::Matrix A(3, 2);

  // Fill A(i,j) = i + 10*j
  for (int i = 0; i < A.rows(); ++i) {
    for (int j = 0; j < A.cols(); ++j) {
      A(i, j) = i + 10 * j;
    }
  }

  // Print A as a 3x2 grid
  for (int i = 0; i < A.rows(); ++i) {
    for (int j = 0; j < A.cols(); ++j) {
      std::printf("%g ", A(i, j));
    }
    std::printf("\n");
  }

  std::printf("raw storage order: ");

  // Print A.data()[0] .. A.data()[5] on one line
  const int n = A.rows() * A.cols();
  for (int k = 0; k < n; ++k) {
    std::printf("%g ", A.data()[k]);
  }
  std::printf("\n");

  return 0;
}