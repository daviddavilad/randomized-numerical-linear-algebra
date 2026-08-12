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

TruncatedSVD truncated_svd(const Matrix& A, int k) {
  Matrix work_copy = A;

  const int m = A.rows();
  const int n = A.cols();
  const int lda = A.ld();
  const int mn = std::min(m, n);

  if (k < 0 || k > mn)
    throw std::invalid_argument("k must be in [0, min(m, n)]");

  const char jobz = 'S';

  Matrix Ufull(m, mn);
  Matrix Vtfull(mn, n);
  const int ldu = Ufull.ld();
  const int ldvt = Vtfull.ld();

  std::vector<double> s(mn);
  std::vector<int> iwork(8 * mn);
  int info = 0;

  double wkopt = 0.0;
  int lwork = -1;
  // Workspace query
  dgesdd_(&jobz, &m, &n, work_copy.data(), &lda, s.data(), Ufull.data(), &ldu, Vtfull.data(), &ldvt, &wkopt, &lwork, iwork.data(), &info);

  lwork = static_cast<int>(wkopt);
  std::vector<double> work(std::max(lwork, 1));

  // Real call
  dgesdd_(&jobz, &m, &n, work_copy.data(), &lda, s.data(), Ufull.data(), &ldu, Vtfull.data(), &ldvt, work.data(), &lwork, iwork.data(), &info);

  if (info != 0)
    throw std::runtime_error("dgesdd failed, info = " + std::to_string(info));

  TruncatedSVD out;
  out.U = Matrix(m, k);
  out.Vt = Matrix(k, n);
  out.s.assign(s.begin(), s.begin() + k);

  // Copy first k columns of Ufull into out.U
  for (int j = 0; j < k; ++j)
    for (int i = 0; i < m; ++i) out.U(i, j) = Ufull(i, j);

  // Copy first k rows of Vtfull into out.Vt
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < k; ++i) out.Vt(i, j) = Vtfull(i, j);

  return out;
}

}  // namespace rnla