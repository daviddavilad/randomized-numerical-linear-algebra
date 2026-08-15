#pragma once

#include <cstdint>
#include <vector>

#include "rnla/matrix.hpp"

namespace rnla {

enum class Spectrum {
  // sigma_i = (i+1)^(-alpha). Heavy-tailed: the HARD case for randomized
  // methods, and where power iteration and block Krylov earn their keep.
  Polynomial,
  // sigma_i = exp(-alpha * i). Fast decay: the easy case.
  Exponential,
};

// A = U * diag(sigma) * V^T with random orthonormal U and V, so the singular
// values are exactly the ones requested. Knowing sigma in closed form is what
// makes the Eckart-Young optimum computable.
struct TestMatrix {
  Matrix A;
  std::vector<double> sigma;  // descending
};

TestMatrix make_test_matrix(int m, int n, Spectrum kind, double alpha, std::uint64_t seed);

}  // namespace rnla