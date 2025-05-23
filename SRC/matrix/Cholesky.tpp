//===----------------------------------------------------------------------===//
//
//                                   xara
//
//===----------------------------------------------------------------------===//
//                              https://xara.so
//===----------------------------------------------------------------------===//
#pragma once
#include <array>
#include <cmath>
#include <cstddef>


template<std::size_t N>
class Cholesky {
  public:
    // Attempt to factor A (row-major size N*N) on construction.
    // If A is not SPD, then ok will be false.
    template <typename MatType>
    constexpr explicit 
    Cholesky(const MatType& A) noexcept
    : ok{true}, L{}
    {
        // Perform in-place Cholesky into L
        for (std::size_t j = 0; j < N; ++j) {
            // diagonal
            double sum_sq = 0;
            for (std::size_t k = 0; k < j; ++k)
                sum_sq += L[j*N + k] * L[j*N + k];
            double diag = A(j, j) - sum_sq;
            if (diag <= 0.0) {
                ok = false; return;
            }
            L[j*N + j] = std::sqrt(diag);

            // sub-diagonal
            for (std::size_t i = j + 1; i < N; ++i) {
                double sum_pr = 0;
                for (std::size_t k = 0; k < j; ++k)
                    sum_pr += L[i*N + k] * L[j*N + k];
                L[i*N + j] = (A(i, j) - sum_pr) / L[j*N + j];
            }
        }
    }

    template<typename MatType>
    int 
    invert(MatType& B) const noexcept
    {
      if (!ok) 
        return -1;

      double e[N], x[N];
      // Solve N systems A·x = e_j to build A⁻¹ one column at a time
      for (std::size_t j = 0; j < N; ++j) {
        // form unit vector e_j
        for (std::size_t i = 0; i < N; ++i) 
            e[i] = (i == j ? 1.0 : 0.0);

        if (solve(e, x) != 0) 
            return -1;

        // write column j of the inverse
        for (std::size_t i = 0; i < N; ++i) 
            B(i, j) = x[i];
      }
      return 0;
    }

    // Solve Ax = b by forward/back substitution through L.  
    // Returns 0 on success; -1 if factorization failed.
    int
    solve(const double* b, double* x) const noexcept
    {
        if (!ok)
          return -1;

        // forward: Ly = b
        double y[N];
        for (std::size_t i = 0; i < N; ++i) {
            double sum = 0;
            for (std::size_t k = 0; k < i; ++k)
                sum += L[i*N + k] * y[k];
            y[i] = (b[i] - sum) / L[i*N + i];
        }

        // backward: Lᵀx = y
        for (std::size_t i = N; i-- > 0; ) {
            double sum = 0;
            for (std::size_t k = i + 1; k < N; ++k)
                sum += L[k*N + i] * x[k];
            x[i] = (y[i] - sum) / L[i*N + i];
        }

        return 0;
    }


    constexpr bool is_ok() const noexcept {
        return ok; 
    }

private:
    bool               ok;
    std::array<double, N*N> L;
};

