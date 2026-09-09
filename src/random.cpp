#include "rnla/random.hpp"

#include <random>
#include <algorithm>
#include <stdexcept>
#include <vector>

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

Matrix sparse_sign(int m, int n, int zeta, std::uint64_t seed) {
  if (zeta < 1 || zeta > m)
    throw std::invalid_argument("zeta must be in [1, m]");

  Matrix S(m, n);

  std::mt19937_64 gen(seed);
  std::uniform_int_distribution<int> row_dist(0, m - 1);
  std::uniform_int_distribution<int> sign_dist(0, 1);

  const double scale = 1.0 / std::sqrt(static_cast<double>(zeta));

  for (int j = 0; j < n; ++j) {
    // zeta distinct rows, by rejection. With zeta << m collisions are rare,
    // so the linear scan costs little.
    std::vector<int> rows;
    rows.reserve(static_cast<std::size_t>(zeta));
    while (rows.size() < static_cast<std::size_t>(zeta)) {
      const int r = row_dist(gen);
      if (std::find(rows.begin(), rows.end(), r) == rows.end()) rows.push_back(r);
    }

    for (int i : rows) {
      const double sign = sign_dist(gen) == 0 ? -1.0 : 1.0;
      S(i, j) = sign * scale;
    }
  }

  return S;
}

Matrix make_sketch(const Sketch& spec, int m, int n, std::uint64_t seed) {
  switch (spec.kind) {
    case Sketch::Gaussian:
      return gaussian(m, n, seed);
    case Sketch::SparseSign:
      return sparse_sign(m, n, spec.zeta, seed);
  }
  throw std::invalid_argument("unknown sketch kind");
}

}  // namespace rnla