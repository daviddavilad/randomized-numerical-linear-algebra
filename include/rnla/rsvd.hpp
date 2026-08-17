#pragma once

#include <cstdint>

#include "rnla/linalg.hpp"
#include "rnla/matrix.hpp"

namespace rnla {

// Randomized rank-k SVD (Halko-Martinsson-Tropp, Algorithm 4.1 + 5.1).
// p is the oversampling parameter: the sketch uses k + p columns and the
// result is truncated back to k.
TruncatedSVD randomized_svd(const Matrix& A, int k, int p, std::uint64_t seed);

}  // namespace rnla