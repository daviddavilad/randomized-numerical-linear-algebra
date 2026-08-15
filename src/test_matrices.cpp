#include "rnla/test_matrices.hpp"

#include <algorithm>
#include <cmath>

#include "rnla/linalg.hpp"
#include "rnla/random.hpp"

namespace rnla {

TestMatrix make_test_matrix(int m, int n, Spectrum kind, double alpha, std::uint64_t seed) {
  const int r = std::min(m, n);

  std::vector<double> sigma(r);
  // Fill sigma according to kind
  switch (kind) {
    case Spectrum::Polynomial:
      for (int i = 0; i < r; ++i) {
        sigma[i] = std::pow(i + 1.0, -alpha);
      }
      break;
    case Spectrum::Exponential:
      for (int i = 0; i < r; ++i) {
        sigma[i] = std::exp(-alpha * i);
      }
      break;
  }

  // U = orth(gaussian(m, r, seed))
  Matrix U = orth(gaussian(m, r, seed));

  // V = orth(gaussian(n, r, seed + 1))
  Matrix V = orth(gaussian(n, r, seed + 1));

  // Scale column j of U by sigma[j]
  for (int j = 0; j < r; ++j) {
    for (int i = 0; i < m; ++i) {
      U(i, j) *= sigma[j];
    }
  }

  TestMatrix out;
  out.A = matmul(U, V, false, true);
  out.sigma = std::move(sigma);
  return out;
}

}  // namespace rnla