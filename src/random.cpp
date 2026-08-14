#include "rnla/random.hpp"

#include <random>

namespace rnla {

Matrix gaussian(int m, int n, std::uint64_t seed) {
  Matrix G(m, n);

  // Create the generator, seeded
  std::mt19937_64 gen(seed);

  // Create the standard normal distribution
  std::normal_distribution<double> dist(0.0, 1.0);

  // Fill G, j outer and i inner
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < m; ++i) {
      G(i, j) = dist(gen);
    }
  }

  return G;
}

}  // namespace rnla