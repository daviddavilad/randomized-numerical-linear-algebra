#include <cstdio>
#include <vector>

#include "rnla/matrix.hpp"
#include "rnla/linalg.hpp"

int main() {
  rnla::Matrix B(2, 3);
  for (int j = 0; j < B.cols(); ++j)
    for (int i = 0; i < B.rows(); ++i) B(i, j) = i + 10 * j;

  rnla::Matrix D = rnla::matmul(B, B, false, true);   // 2x2, symmetric PSD

  std::vector<double> s = rnla::singular_values(B);  // singular values of B
  std::vector<double> t = rnla::singular_values(D);  // singular values of D

  for (double x : s) std::printf("%g ", x);
  std::printf("\n");
  for (double x : t) std::printf("%g ", x);
  std::printf("\n");

  double traceD = 0.0;
  for (int i = 0; i < D.rows(); ++i) traceD += D(i, i);
  std::printf("trace(D) = %g\n", traceD);
  std::printf("sum(s^2) = %g\n", s[0] * s[0] + s[1] * s[1]);
  std::printf("sum(t) = %g\n", t[0] + t[1]);

  for (std::size_t i = 0; i < t.size(); ++i) {
    double expected = s[i] * s[i];
    std::printf("t[%zu] = %g, expected = %g\n", i, t[i], expected);
  }

  return 0;
}