#pragma once

#include <cstdint>

#include "rnla/matrix.hpp"

namespace rnla {

// m x n matrix of i.i.d. standard normal entries.
// The seed is required: every experiment in this project must be
// exactly reproducible from its inputs.
Matrix gaussian(int m, int n, std::uint64_t seed);

}  // namespace rnla