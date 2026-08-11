#pragma once

extern "C" {

// C := alpha * op(A) * op(B) + beta * C, op(X) = X or X^T.
void dgemm_(
    const char* transa,
    const char* transb,
    const int* m,
    const int* n,
    const int* k,
    const double* alpha,
    const double* a,
    const int* lda,
    const double* b,
    const int* ldb,
    const double* beta,
    double* c,
    const int* ldc
);

}