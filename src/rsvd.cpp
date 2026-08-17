#include "rnla/rsvd.hpp"

#include <algorithm>
#include <stdexcept>

#include "rnla/random.hpp"

namespace rnla {

TruncatedSVD randomized_svd(const Matrix& A, int k, int p, std::uint64_t seed) {
  const int m = A.rows();
  const int n = A.cols();

  if (k < 0 || k > std::min(m, n))
    throw std::invalid_argument("k must be in [0, min(m, n)]");
  if (p < 0) throw std::invalid_argument("p must be non-negative");

  const int l = std::min(k + p, std::min(m, n));

  // Omega = gaussian(n, l, seed) (n x l)
  Matrix Omega = gaussian(n, l, seed);

  // Y = A * Omega (m x l)
  Matrix Y = matmul(A, Omega);

  // Q = orth(Y) (m x l)
  Matrix Q = orth(Y);

  // B = Q^T * A (l x n)
  Matrix B = matmul(Q, A, true, false);
  
  // TruncatedSVD svd_B = truncated_svd(B, k)
  TruncatedSVD svd_B = truncated_svd(B, k);

  TruncatedSVD out;
  // out.U = Q * svd_B.U (m x k)
  out.U = matmul(Q, svd_B.U);
  out.s = svd_B.s;
  out.Vt = svd_B.Vt;
  return out;
}

}  // namespace rnla