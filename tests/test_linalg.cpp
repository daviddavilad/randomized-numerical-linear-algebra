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

rnla::Matrix make_hilbert(int m, int n) {
  rnla::Matrix A(m, n);
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < m; ++i) A(i, j) = 1.0 / (i + j + 1.0);
  return A;
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

void test_svd_reconstruction() {
  rnla::Matrix A = make_hilbert(4, 3);   // Hilbert-like, full rank

  const int k = 3;   // = min(m, n), so this should be exact
  rnla::TruncatedSVD svd = rnla::truncated_svd(A, k);

  check(svd.U.rows() == 4 && svd.U.cols() == k, "U has shape m x k");
  check(svd.Vt.rows() == k && svd.Vt.cols() == 3, "Vt has shape k x n");
  check(svd.s.size() == static_cast<std::size_t>(k), "s has k entries");

  // Scale the columns of a copy of U by s, then multiply by Vt
  rnla::Matrix Uscaled = svd.U;
  for (int j = 0; j < k; ++j)
    for (int i = 0; i < Uscaled.rows(); ++i)
      Uscaled(i, j) *= svd.s[j];

  rnla::Matrix reconstruction = rnla::matmul(Uscaled, svd.Vt);

  // Check every entry against A
  for (int j = 0; j < A.cols(); ++j)
    for (int i = 0; i < A.rows(); ++i)
      check(close(reconstruction(i, j), A(i, j), 1e-10), "reconstruction matches original");
}

void test_svd_orthonormality() {
  // Same A, k = 2 this time
  rnla::Matrix A = make_hilbert(4, 3);   // Hilbert-like, full rank
  const int k = 2;
  rnla::TruncatedSVD svd = rnla::truncated_svd(A, k);

  // G = U^T U should be the 2x2 identity
  rnla::Matrix G = rnla::matmul(svd.U, svd.U, true, false);
  check(G.rows() == k && G.cols() == k, "G has shape k x k");
  for (int j = 0; j < G.cols(); ++j)
    for (int i = 0; i < G.rows(); ++i) {
      double expected = (i == j) ? 1.0 : 0.0;
      check(close(G(i, j), expected), "U^T U is identity");
    }
}

}  // namespace

int main() {
  std::printf("rnla tests\n");
  test_diagonal_singular_values();
  test_gram_identity();
  test_svd_reconstruction();
  test_svd_orthonormality();
  std::printf("%s (%d failures)\n", failures ? "FAILED" : "OK", failures);
  return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}