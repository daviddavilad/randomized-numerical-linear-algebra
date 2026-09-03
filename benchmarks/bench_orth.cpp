#include <cmath>
#include <cstdio>
#include <exception>

#include "rnla/linalg.hpp"
#include "rnla/matrix.hpp"
#include "rnla/test_matrices.hpp"

int main() {
  const int m = 200, n = 30;
  const double kappas[] = {1e2, 1e4, 1e6, 1e7, 1e8, 1e10, 1e12};

    
  std::printf("%10s %12s %12s %14s %14s\n", "kappa", "kappa_act", "kappa^2", "householder", "choleskyqr");

  for (double kappa : kappas) {
    // exponential spectrum over r = min(m,n) values has
    //   kappa = sigma_0 / sigma_{r-1} = e^{alpha(r-1)}
    // so alpha = ln(kappa) / (r-1)
    const double alpha = std::log(kappa) / (n - 1);
    auto tm = rnla::make_test_matrix(m, n, rnla::Spectrum::Exponential, alpha, 7);
    const double kappa_actual = tm.sigma.front() / tm.sigma.back();

    const double eh = rnla::orthogonality_error(rnla::orth(tm.A));

    // cholesky_qr in a try/catch — it throws past the threshold.
    //       print the error on success, "FAILED" on catch.
    double ec = 0.0;
    try {
      ec = rnla::orthogonality_error(rnla::cholesky_qr(tm.A));
      std::printf("%10.2e %12.2e %12.2e %14.2e %14.2e\n", kappa, kappa_actual, kappa * kappa, eh, ec);
    } catch (const std::exception& e) {
      std::printf("%10.2e %12.2e %12.2e %14.2e   FAILED: %s\n", kappa, kappa_actual, kappa * kappa, eh, e.what());
    }
  }
  return 0;
}