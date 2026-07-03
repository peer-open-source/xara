//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, OpenSees/Xara Developers
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
//
// June 2026
//
#pragma once

#include <cmath>
#include <cstring>
#include <type_traits>

#ifndef RESTRICT
#  if defined(_MSC_VER)
#    define RESTRICT __restrict
#  elif defined(__GNUC__) || defined(__clang__)
#    define RESTRICT __restrict__
#  else
#    define RESTRICT
#  endif
#endif

#ifndef SMALL_SOLVE_INLINE
#  if defined(_MSC_VER)
#    define SMALL_SOLVE_INLINE __forceinline
#  elif defined(__GNUC__) || defined(__clang__)
#    define SMALL_SOLVE_INLINE inline __attribute__((always_inline))
#  else
#    define SMALL_SOLVE_INLINE inline
#  endif
#endif


namespace small_solve_detail {

template<int I, int End, class F>
SMALL_SOLVE_INLINE void static_for(F&& f) noexcept {
  if constexpr (I < End) {
    f(std::integral_constant<int, I>{});
    static_for<I + 1, End>(f);
  }
}

template<int Begin, int End, class F>
SMALL_SOLVE_INLINE void static_for_rev(F&& f) noexcept {
  if constexpr (Begin < End) {
    f(std::integral_constant<int, End - 1>{});
    static_for_rev<Begin, End - 1>(f);
  }
}

SMALL_SOLVE_INLINE double abs_d(double x) noexcept {
#if defined(__GNUC__) || defined(__clang__)
  return __builtin_fabs(x);
#else
  return std::fabs(x);
#endif
}


// Swap rows R0 and R1, but only on active columns K..N-1.
// Earlier columns are not needed by this fused one-RHS solve.
template<int R0, int R1, int K, int N>
SMALL_SOLVE_INLINE void swap_active_rows(
    double* RESTRICT W,
    double* RESTRICT x) noexcept
{
  if constexpr (R0 != R1) {
    static_for<K, N>([&](auto Jc) noexcept {
        constexpr int j = decltype(Jc)::value;

        const double t = W[R0 + j * N];
        W[R0 + j * N] = W[R1 + j * N];
        W[R1 + j * N] = t;
    });

    const double tx = x[R0];
    x[R0] = x[R1];
    x[R1] = tx;
  }
}


// Convert dynamic pivot p into compile-time row-swap indices.
template<int K, int P, int N>
SMALL_SOLVE_INLINE void dispatch_active_row_swap(
    int p,
    double* RESTRICT W,
    double* RESTRICT x) noexcept
{
  if constexpr (P < N) {
    if (p == P) {
      swap_active_rows<K, P, K, N>(W, x);
    } else {
      dispatch_active_row_swap<K, P + 1, N>(p, W, x);
    }
  }
}


// LU factorization with partial pivoting.
// W is packed column-major, W[i + j*N].
// x is the RHS workspace.
//
// Diagonal entries W[k + k*N] are replaced by their reciprocals.
template<int K, int N>
SMALL_SOLVE_INLINE int factor_step(
    double* RESTRICT W,
    double* RESTRICT x) noexcept
{
  if constexpr (K >= N) {
    return 0;
  }
  else {
    int p = K;
    double amax = abs_d(W[K + K * N]);

    static_for<K + 1, N>([&](auto Ic) noexcept {
      constexpr int i = decltype(Ic)::value;

      const double v = abs_d(W[i + K * N]);
      const bool take = v > amax;

      amax = take ? v : amax;
      p = take ? i : p;
    });

    if (amax == 0.0) [[unlikely]] {
      return K + 1;
    }

    dispatch_active_row_swap<K, K, N>(p, W, x);

    const double inv_pivot = 1.0 / W[K + K * N];
    W[K + K * N] = inv_pivot;

    const double xk = x[K];

    // Row-oriented update.
    // The multiplier lik is not stored because this one-RHS fused solve
    // never reuses the L factor.
    static_for<K + 1, N>([&](auto Ic) noexcept {
        constexpr int i = decltype(Ic)::value;

      const double lik = W[i + K * N] * inv_pivot;

      x[i] -= lik * xk;

      static_for<K + 1, N>([&](auto Jc) noexcept {
        constexpr int j = decltype(Jc)::value;
        W[i + j * N] -= lik * W[K + j * N];
      });
    });

    return factor_step<K + 1, N>(W, x);
  }
}


template<int N>
SMALL_SOLVE_INLINE void back_substitute(
    double* RESTRICT W,
    double* RESTRICT x) noexcept
{
  static_for_rev<0, N>([&](auto Kc) noexcept {
    constexpr int k = decltype(Kc)::value;

    const double xk = x[k] * W[k + k * N];
    x[k] = xk;

    static_for<0, k>([&](auto Ic) noexcept {
        constexpr int i = decltype(Ic)::value;
        x[i] -= W[i + k * N] * xk;
    });
  });
}


template<int N>
SMALL_SOLVE_INLINE int solve_generic(
    const double* RESTRICT A,
    const double* RESTRICT b,
    double* RESTRICT x) noexcept
{
  double W[N * N];

  std::memcpy(W, A, sizeof(W));
  std::memcpy(x, b, sizeof(double) * N);

  const int info = factor_step<0, N>(W, x);

  if (info != 0) [[unlikely]]
    return -info;


  back_substitute<N>(W, x);
  return 0;
}


SMALL_SOLVE_INLINE int 
solve_1(
    const double* RESTRICT A,
    const double* RESTRICT b,
    double* RESTRICT x) noexcept
{
  if (A[0] == 0.0) [[unlikely]]
    return -1;

  x[0] = b[0] / A[0];
  return 0;
}


SMALL_SOLVE_INLINE int 
solve_2(
    const double* RESTRICT A,
    const double* RESTRICT b,
    double* RESTRICT x) noexcept
{
  double w00 = A[0];
  double w10 = A[1];

  double w01 = A[2];
  double w11 = A[3];

  double y0 = b[0];
  double y1 = b[1];

  double amax = abs_d(w00);
  const double v10 = abs_d(w10);

  if (v10 > amax) {
    amax = v10;

    const double t0 = w00;
    w00 = w10;
    w10 = t0;

    const double t1 = w01;
    w01 = w11;
    w11 = t1;

    const double tb = y0;
    y0 = y1;
    y1 = tb;
  }

  if (amax == 0.0) [[unlikely]]
    return -1;


  const double inv0 = 1.0 / w00;
  const double l10 = w10 * inv0;

  y1 -= l10 * y0;
  w11 -= l10 * w01;

  if (w11 == 0.0) [[unlikely]]
    return -2;

  const double inv1 = 1.0 / w11;

  const double z1 = y1 * inv1;
  const double z0 = (y0 - w01 * z1) * inv0;

  x[0] = z0;
  x[1] = z1;

  return 0;
}


SMALL_SOLVE_INLINE int 
solve_3(
    const double* RESTRICT A,
    const double* RESTRICT b,
    double* RESTRICT x) noexcept
{
  double w00 = A[0];
  double w10 = A[1];
  double w20 = A[2];

  double w01 = A[3];
  double w11 = A[4];
  double w21 = A[5];

  double w02 = A[6];
  double w12 = A[7];
  double w22 = A[8];

  double y0 = b[0];
  double y1 = b[1];
  double y2 = b[2];

  int p0 = 0;
  double amax = abs_d(w00);

  {
    const double v = abs_d(w10);
    const bool take = v > amax;
    amax = take ? v : amax;
    p0 = take ? 1 : p0;
  }

  {
    const double v = abs_d(w20);
    const bool take = v > amax;
    amax = take ? v : amax;
    p0 = take ? 2 : p0;
  }

  if (amax == 0.0) [[unlikely]]
    return -1;


  if (p0 == 1) {
    const double t00 = w00;
    w00 = w10;
    w10 = t00;

    const double t01 = w01;
    w01 = w11;
    w11 = t01;

    const double t02 = w02;
    w02 = w12;
    w12 = t02;

    const double tb = y0;
    y0 = y1;
    y1 = tb;
  } 
  else if (p0 == 2) {
    const double t00 = w00;
    w00 = w20;
    w20 = t00;

    const double t01 = w01;
    w01 = w21;
    w21 = t01;

    const double t02 = w02;
    w02 = w22;
    w22 = t02;

    const double tb = y0;
    y0 = y2;
    y2 = tb;
  }

  const double inv0 = 1.0 / w00;

  const double l10 = w10 * inv0;
  y1 -= l10 * y0;
  w11 -= l10 * w01;
  w12 -= l10 * w02;

  const double l20 = w20 * inv0;
  y2 -= l20 * y0;
  w21 -= l20 * w01;
  w22 -= l20 * w02;

  double amax1 = abs_d(w11);
  const double v21 = abs_d(w21);
  const bool take21 = v21 > amax1;

  amax1 = take21 ? v21 : amax1;

  if (amax1 == 0.0) [[unlikely]]
    return -2;


  if (take21) {
    const double t11 = w11;
    w11 = w21;
    w21 = t11;

    const double t12 = w12;
    w12 = w22;
    w22 = t12;

    const double tb = y1;
    y1 = y2;
    y2 = tb;
  }

  const double inv1 = 1.0 / w11;

  const double l21 = w21 * inv1;
  y2 -= l21 * y1;
  w22 -= l21 * w12;

  if (w22 == 0.0) [[unlikely]] {
    return -3;
  }

  const double inv2 = 1.0 / w22;

  const double z2 = y2 * inv2;
  const double z1 = (y1 - w12 * z2) * inv1;
  const double z0 = (y0 - w01 * z1 - w02 * z2) * inv0;

  x[0] = z0;
  x[1] = z1;
  x[2] = z2;

  return 0;
}

} // namespace small_solve_detail


//
// Solve A x = b for column-major A.
//
// Storage:
//     A[i + j * N] is entry (i,j).
//
// NOTE:
//     A, b, and x must not overlap.
//
template<int N>
SMALL_SOLVE_INLINE 
int small_solve(
    const double* RESTRICT A,
    const double* RESTRICT b,
    double* RESTRICT x) noexcept
{
  static_assert(N > 0, "N must be positive");

  if constexpr (N == 1) {
    return small_solve_detail::solve_1(A, b, x);
  } else if constexpr (N == 2) {
    return small_solve_detail::solve_2(A, b, x);
  } else if constexpr (N == 3) {
    return small_solve_detail::solve_3(A, b, x);
  } else {
    return small_solve_detail::solve_generic<N>(A, b, x);
  }
}
