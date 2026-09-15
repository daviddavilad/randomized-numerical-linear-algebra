#include <cstdint>
#include <cstdio>
#include <iterator>
#include <vector>
#include <cmath>

#include "rnla/linalg.hpp"
#include "rnla/matrix.hpp"
#include "rnla/random.hpp"
#include "rnla/rsvd.hpp"
#include "rnla/test_matrices.hpp"

namespace {

struct Stat { double mean, se; };

Stat summarize(const std::vector<double>& v) {
  const double n = static_cast<double>(v.size());
  double mean = 0.0;
  for (double x : v) mean += x;
  mean /= n;
  double var = 0.0;
  for (double x : v) var += (x - mean) * (x - mean);
  var /= (n - 1.0);           // Bessel: sample variance, not population
  return {mean, std::sqrt(var / n)};
}

}  // namespace

int main() {
  const int m = 300, n = 200, k = 20, p = 10;
  const int n_seeds = 200;

  struct Case { rnla::Spectrum kind; const char* label; double alpha; };
  const Case cases[] = {
      {rnla::Spectrum::Exponential, "exp", 0.3},
      {rnla::Spectrum::Exponential, "exp", 0.1},
      {rnla::Spectrum::Polynomial, "poly", 2.0},
      {rnla::Spectrum::Polynomial, "poly", 1.0},
  };
  const int zetas[] = {1, 2, 4, 8, 16};
  const std::size_t n_zetas = std::size(zetas);

  std::printf("%-6s %6s %8s %10s", "family", "alpha", "gap", "gaussian");
  for (int z : zetas) {
    char buf[16];
    std::snprintf(buf, sizeof buf, "d(zeta=%d)", z);
    std::printf(" %17s", buf);
  }
  std::printf("\n");

  for (const Case& c : cases) {
    // Spectrum-derived quantities are the same for every rotation, so take
    // them from one representative matrix.
    auto tm0 = rnla::make_test_matrix(m, n, c.kind, c.alpha, 0);
    const double gap = tm0.sigma[k] / tm0.sigma[k - 1];
    const double opt = rnla::eckart_young_fro(tm0.sigma, k);

    std::vector<std::vector<double>> samples(1 + n_zetas);

    for (std::uint64_t seed = 0; seed < static_cast<std::uint64_t>(n_seeds);
         ++seed) {
      // New Haar rotation per seed: U and V change, the spectrum does not.
      auto tm = rnla::make_test_matrix(m, n, c.kind, c.alpha, seed);

      auto svd_g = rnla::randomized_svd(tm.A, k, p, 0, seed);
      samples[0].push_back(rnla::reconstruction_error(tm.A, svd_g) / opt);

      for (std::size_t i = 0; i < n_zetas; ++i) {
        rnla::Sketch spec;
        spec.kind = rnla::Sketch::SparseSign;
        spec.zeta = zetas[i];
        auto svd_s = rnla::randomized_svd(tm.A, k, p, 0, seed, spec);
        samples[i + 1].push_back(rnla::reconstruction_error(tm.A, svd_s) / opt);
      }
    }

    std::printf("%-6s %6.1f %8.4f %10.6f", c.label, c.alpha, gap, summarize(samples[0]).mean);

    for (std::size_t i = 0; i < n_zetas; ++i) {
      std::vector<double> d(samples[i + 1].size());
      for (std::size_t j = 0; j < d.size(); ++j)
        d[j] = samples[i + 1][j] - samples[0][j];
      const Stat st = summarize(d);
      std::printf(" %9.6f±%.6f", st.mean, st.se);
    }
    std::printf("\n");
  }
  return 0;
}