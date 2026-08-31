#pragma once

#include <cstdint>

#include "rnla/linalg.hpp"
#include "rnla/matrix.hpp"

namespace rnla {

// Randomized rank-k SVD (Halko-Martinsson-Tropp, Algorithms 4.1 / 4.3 / 5.1).
// p: oversampling — the sketch uses k + p columns, truncated back to k.
// q: power iterations — the sketch becomes (A A^T)^q A Omega, which raises
//    each singular value to the power 2q+1 and so widens the relative gaps.
//    Each iteration costs two additional passes over A.
TruncatedSVD randomized_svd(const Matrix& A, int k, int p, int q, std::uint64_t seed);

// Randomized block Krylov SVD (Musco & Musco 2015).
// Builds the Krylov block [A*Omega, (A A^T)A*Omega, ..., (A A^T)^q A*Omega],
// an m x (q+1)l matrix, and uses its range instead of only the last iterate.
// Same number of passes over A as power iteration at the same q; more memory
// and a larger orthogonalization and small SVD.
TruncatedSVD randomized_svd_krylov(const Matrix& A, int k, int p, int q, std::uint64_t seed);

}  // namespace rnla