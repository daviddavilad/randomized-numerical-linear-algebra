#pragma once

#include <cstdint>

#include "rnla/matrix.hpp"

namespace rnla {

// m x n matrix of i.i.d. standard normal entries.
// The seed is required: every experiment in this project must be
// exactly reproducible from its inputs.
Matrix gaussian(int m, int n, std::uint64_t seed);

// m x n sparse sign matrix: each column has exactly zeta nonzeros, each
// +/- 1/sqrt(zeta), at uniformly random rows without replacement.
// Returned dense for now — the cost advantage needs a sparse apply to realize.
Matrix sparse_sign(int m, int n, int zeta, std::uint64_t seed);

}  // namespace rnla