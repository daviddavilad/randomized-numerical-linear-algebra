#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "rnla/linalg.hpp"
#include "rnla/matrix.hpp"

namespace {

int failures = 0;

void check(bool ok, const char* what) {
  std::printf("  [%s] %s\n", ok ? "PASS" : "FAIL", what);
  if (!ok) ++failures;
}

bool close(double a, double b, double tol = 1e-12) {
  return std::abs(a - b) <= tol * std::max(1.0, std::abs(b));
}

void test_diagonal_singular_values() {
  rnla::Matrix A(3, 3);
  A(0, 0) = 1.0;
  A(1, 1) = 2.0;
  A(2, 2) = 3.0;

  std::vector<double> s = rnla::singular_values(A);

  check(s.size() == 3, "size of singular values");
  check(close(s[0], 3.0), "largest singular value");
  check(close(s[1], 2.0), "middle singular value");
  check(close(s[2], 1.0), "smallest singular value");
}

void test_gram_identity() {
  rnla::Matrix B(2, 3);
  B(0, 0) = 1.0;
  B(0, 1) = 2.0;
  B(0, 2) = 3.0;
  B(1, 0) = 4.0;
  B(1, 1) = 5.0;
  B(1, 2) = 6.0;

  rnla::Matrix D = rnla::matmul(B, B, false, true);   // D = B * B^T

  std::vector<double> s = rnla::singular_values(B);   // singular values of B
  std::vector<double> t = rnla::singular_values(D);   // singular values of D

  check(s.size() == t.size(), "size of singular values");
  for (std::size_t i = 0; i < s.size(); ++i) {
    check(close(t[i], s[i] * s[i]), "singular value squared");
  }
}

}  // namespace

int main() {
  std::printf("rnla tests\n");
  test_diagonal_singular_values();
  test_gram_identity();
  std::printf("%s (%d failures)\n", failures ? "FAILED" : "OK", failures);
  return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}