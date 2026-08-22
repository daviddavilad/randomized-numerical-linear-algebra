#include <algorithm>
#include <cstdio>
#include <numeric>
#include <vector>

#include "rnla/linalg.hpp"
#include "rnla/matrix.hpp"
#include "rnla/rsvd.hpp"
#include "rnla/test_matrices.hpp"

namespace {

struct Case {
  rnla::Spectrum kind;
  const char* label;
  double alpha;
};

}  // namespace

int main() {
  const int m = 300, n = 200, k = 20, p = 10;
  const int n_seeds = 10;

  const Case cases[] = {
      {rnla::Spectrum::Exponential, "exp", 1.0},
      {rnla::Spectrum::Exponential, "exp", 0.5},
      {rnla::Spectrum::Exponential, "exp", 0.3},
      {rnla::Spectrum::Exponential, "exp", 0.2},
      {rnla::Spectrum::Exponential, "exp", 0.1},
      {rnla::Spectrum::Polynomial, "poly", 3.0},
      {rnla::Spectrum::Polynomial, "poly", 2.0},
      {rnla::Spectrum::Polynomial, "poly", 1.5},
      {rnla::Spectrum::Polynomial, "poly", 1.0},
      {rnla::Spectrum::Polynomial, "poly", 0.5},
      {rnla::Spectrum::Flat,       "flat", 0.0},
  };

  std::printf("%-6s %6s %8s %10s %12s %12s %12s\n", "family", "alpha", "gap", "opt_rel", "q=0", "q=1", "q=2");

  for (const Case& c : cases) {
    auto tm = rnla::make_test_matrix(m, n, c.kind, c.alpha, 7);

    const double gap = tm.sigma[k] / tm.sigma[k - 1];
    const double optimal = rnla::eckart_young_fro(tm.sigma, k);
    const double opt_rel = optimal / rnla::norm_fro(tm.A);

    double means[3];
    for (int q = 0; q <= 2; ++q) {
      std::vector<double> ratios;
      for (int seed = 0; seed < n_seeds; ++seed) {
        auto svd = rnla::randomized_svd(tm.A, k, p, q, seed);
        ratios.push_back(rnla::reconstruction_error(tm.A, svd) / optimal);
      }
      means[q] = std::accumulate(ratios.begin(), ratios.end(), 0.0) / ratios.size();
    }

    std::printf("%-6s %6.1f %8.4f %10.2e %12.9f %12.9f %12.9f\n", c.label, c.alpha, gap, opt_rel, means[0], means[1], means[2]);
  }
  return 0;
}