#include "rnla/rsvd.hpp"

#include <algorithm>
#include <stdexcept>

#include "rnla/random.hpp"

namespace rnla {

TruncatedSVD randomized_svd(const Matrix& A, int k, int p, int q, std::uint64_t seed) {
  const int m = A.rows();
  const int n = A.cols();

  if (k < 0 || k > std::min(m, n))
    throw std::invalid_argument("k must be in [0, min(m, n)]");
  if (p < 0) throw std::invalid_argument("p must be non-negative");

  if (q < 0) throw std::invalid_argument("q must be non-negative");

  const int l = std::min(k + p, std::min(m, n));

  // Omega = gaussian(n, l, seed) (n x l)
  Matrix Omega = gaussian(n, l, seed);

  // Y = A * Omega (m x l)
  Matrix Y = matmul(A, Omega);

  // Naive power iteration: Y <- (A A^T)^q A Omega.
  // Numerically unstable — see step 3.
  // Power iteration with re-orthonormalization after every application of
  // A^T and A (HMT Algorithm 4.4). Without this, the columns of Y span
  // directions differing by (sigma_1/sigma_l)^(2q), which exceeds double
  // precision once that quantity passes 1/eps.
  for (int t = 0; t < q; ++t) {
    Y = orth(matmul(A, Y, true, false));   // n x l
    Y = orth(matmul(A, Y, false, false));  // m x l
  }

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