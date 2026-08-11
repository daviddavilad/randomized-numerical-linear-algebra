#include "rnla/linalg.hpp"

#include <algorithm>
#include <stdexcept>
#include <vector>
#include <string>

#include "rnla/lapack.hpp"

namespace rnla {

std::vector<double> singular_values(const Matrix& A) {
  Matrix work_copy = A;

  const int m = A.rows();
  const int n = A.cols();
  const int lda = A.ld();
  const int mn = std::min(m, n);

  const char jobz = 'N';
  int ldu = 1, ldvt = 1;

  std::vector<double> s(mn);
  std::vector<int> iwork(8 * mn);
  int info = 0;

  // --- workspace query -------------------------------------------------
  double wkopt = 0.0;
  int lwork = -1;
  // Call dgesdd_ with &wkopt as work and lwork = -1
  dgesdd_(&jobz, &m, &n, work_copy.data(), &lda, s.data(), nullptr, &ldu, nullptr, &ldvt, &wkopt, &lwork, iwork.data(), &info);

  lwork = static_cast<int>(wkopt);
  std::vector<double> work(std::max(lwork, 1));

  // --- real call -------------------------------------------------------
  // Call dgesdd_ again, with work.data() and the real lwork
  dgesdd_(&jobz, &m, &n, work_copy.data(), &lda, s.data(), nullptr, &ldu, nullptr, &ldvt, work.data(), &lwork, iwork.data(), &info);

  if (info != 0)
    throw std::runtime_error("dgesdd failed, info = " + std::to_string(info));

  return s;
}

}  // namespace rnla