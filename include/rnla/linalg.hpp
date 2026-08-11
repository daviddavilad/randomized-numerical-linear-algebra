#pragma once

#include <vector>

#include "rnla/matrix.hpp"

namespace rnla {

// All singular values of A, in descending order.
std::vector<double> singular_values(const Matrix& A);

}  // namespace rnla