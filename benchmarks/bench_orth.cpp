#include <cmath>
#include <cstdio>
#include <exception>

#include "rnla/linalg.hpp"
#include "rnla/matrix.hpp"
#include "rnla/test_matrices.hpp"
#include "rnla/timing.hpp"

int main() {
  const int m = 200, n = 30;
  const double kappas[] = {1e2, 1e4, 1e6, 1e7, 1e8, 1e10, 1e12};

  std::printf("%10s %14s %14s %12s %14s %10s %10s %10s\n", "kappa", "householder", "cholqr", "cond(Q1)", "cholqr2", "t_hh(us)", "t_c(us)", "t_c2(us)");
  for (double kappa : kappas) {
    // exponential spectrum over r = min(m,n) values has
    //   kappa = sigma_0 / sigma_{r-1} = e^{alpha(r-1)}
    // so alpha = ln(kappa) / (r-1)
    const double alpha = std::log(kappa) / (n - 1);
    auto tm = rnla::make_test_matrix(m, n, rnla::Spectrum::Exponential, alpha, 7);
    // Sanity: the constructed spectrum really has the requested conditioning.
    const double kappa_actual = tm.sigma.front() / tm.sigma.back();
    if (std::abs(kappa_actual - kappa) > 1e-6 * kappa)
      std::printf("  WARNING: kappa mismatch %.3e vs %.3e\n", kappa_actual, kappa);

    const double eh = rnla::orthogonality_error(rnla::orth(tm.A));

    // cholesky_qr in a try/catch — it throws past the threshold.
    //       print the error on success, "FAILED" on catch.
    double ec = 0.0;
    try {
      rnla::Matrix Q1 = rnla::cholesky_qr(tm.A);
      ec = rnla::orthogonality_error(Q1);
      auto s1 = rnla::singular_values(Q1);
      const double kappa_q1 = s1.front() / s1.back();
      const double ec2 = rnla::orthogonality_error(rnla::cholesky_qr(Q1));
      const double th = rnla::time_median([&] {
        volatile double s = rnla::norm_fro(rnla::orth(tm.A));
        (void)s;
      });
      const double tc = rnla::time_median([&] {
        volatile double s = rnla::norm_fro(rnla::cholesky_qr(tm.A));
        (void)s;
      });
      const double tc2 = rnla::time_median([&] {
        volatile double s = rnla::norm_fro(rnla::cholesky_qr2(tm.A));
        (void)s;
      });
      std::printf("%10.2e %14.2e %14.2e %12.2e %14.2e %10.1f %10.1f %10.1f\n", kappa, eh, ec, kappa_q1, ec2, th * 1e6, tc * 1e6, tc2 * 1e6);
    } catch (const std::exception& e) {
      std::printf("%10.2e %14.2e   FAILED: %s\n", kappa, eh, e.what());
    }
  }
  return 0;
}