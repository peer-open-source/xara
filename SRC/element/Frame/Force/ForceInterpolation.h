//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, Claudio M. Perez
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
#include <FrameSection.h>
#include <MatrixND.h>
#include <Vector3D.h>

#if defined(__GNUC__) || defined(__clang__)
#  define XARA_FORCE_INTERP_INLINE [[gnu::always_inline]] inline
#elif defined(_MSC_VER)
#  define XARA_FORCE_INTERP_INLINE __forceinline
#else
#  define XARA_FORCE_INTERP_INLINE inline
#endif


#if defined(__GNUC__) || defined(__clang__)
  #define X_RESTRICT __restrict__
#elif defined(_MSC_VER)
  #define X_RESTRICT __restrict
#else
  #define X_RESTRICT
#endif

template <int nsr, int nwm, int NBV, int NDF, const FrameStressLayout& scheme>
class ForceInterpolation {
public:
  enum : int {
    imy =   3, //  4
    imz =   1, //  5
    iwx =   6, //
    //
    jnx =   0, //  6
    jmx =   5, //  9
    jmy =   4, // 10
    jmz =   2, // 11
    jwx =   7,
  };

  enum : int {
    NNW = 6,
    NLF = 2*NDF
  };

  static_assert(NDF >= 6 + nwm,
                "ForceInterpolation: NDF must contain 6 frame DOFs plus nwm warping DOFs");
  static_assert(NBV >= 6 + 2*nwm,
                "ForceInterpolation: NBV must contain 6 basic resultants plus 2*nwm warping resultants");

  MatrixND<2*NDF,NBV>
  reshape_matrix()
  {
    MatrixND<2*NDF, NBV> Tb {};
    Tb(0*NDF+0, jnx) = -1.0;
    Tb(0*NDF+3, jmx) = -1.0;
    Tb(0*NDF+4, imy) =  1.0;
    Tb(0*NDF+5, imz) =  1.0;
    Tb(1*NDF+0, jnx) =  1.0;
    Tb(1*NDF+3, jmx) =  1.0;
    Tb(1*NDF+4, jmy) =  1.0;
    Tb(1*NDF+5, jmz) =  1.0;
    for (int i=0; i<nwm; i++) {
      Tb(0*NDF+6+i, NNW+i*2)   =  1.0;
      Tb(1*NDF+6+i, NNW+i*2+1) =  1.0;
    }
    return Tb;
  }

  XARA_FORCE_INTERP_INLINE static void
  expand(const VectorND<NBV>& q_pres, VectorND<2*NDF>& pl) 
  noexcept
  {
    const double* X_RESTRICT q = &q_pres[0];
    double*       X_RESTRICT p = &pl[0];

    pl.zero();

    p[0]       = -q[jnx];
    p[NDF]     =  q[jnx];
    p[3]       = -q[jmx];
    p[NDF+3]   =  q[jmx];
    p[4]       =  q[imy];
    p[NDF+4]   =  q[jmy];
    p[5]       =  q[imz];
    p[NDF+5]   =  q[jmz];

    expand_vector_warp<0>(q, p);
  }

  XARA_FORCE_INTERP_INLINE static void
  expand(const MatrixND<NBV,NBV>& kb, MatrixND<2*NDF,2*NDF>& kl) 
  noexcept
  {
    const double* X_RESTRICT k = kb.data();
    double*       X_RESTRICT K = kl.data();

    kl.zero();

    expand_matrix_column(K, k, 0,       jnx, -1.0);
    expand_matrix_column(K, k, NDF,     jnx,  1.0);
    expand_matrix_column(K, k, 3,       jmx, -1.0);
    expand_matrix_column(K, k, NDF+3,   jmx,  1.0);
    expand_matrix_column(K, k, 4,       imy,  1.0);
    expand_matrix_column(K, k, NDF+4,   jmy,  1.0);
    expand_matrix_column(K, k, 5,       imz,  1.0);
    expand_matrix_column(K, k, NDF+5,   jmz,  1.0);

    expand_matrix_warp_columns<0>(K, k);
  }


  //
  // Interpolation
  //
  XARA_FORCE_INTERP_INLINE const MatrixND<nsr,NBV>&
  b(double xL, double L) noexcept
  {
    B.zero();

    const double jsx = (L == 0.0) ? 0.0 : 1.0/L;
    const double xL1 = xL - 1.0;
    b_impl<0>(B.data(), xL, xL1, jsx);
    return B;
  }

  XARA_FORCE_INTERP_INLINE static void
  interpolate(const VectorND<NBV>& q, double xL, double jsx,
              VectorND<nsr>& s) noexcept
  {
    const double* X_RESTRICT qp = &q[0];
    double*       X_RESTRICT sp = &s[0];
    const double xL1 = xL - 1.0;

    interpolate_impl<0>(qp, xL, xL1, jsx, sp);
  }


  //
  // Integration
  //
  XARA_FORCE_INTERP_INLINE static void
  integrate(const VectorND<nsr>& e, double xL, double weight, double jsx,
            VectorND<NBV>& v)
  noexcept
  {
    // do v += B(x)' * e(x) * weight with B = b(x)
    const double* X_RESTRICT ep = &e[0];
    double*       X_RESTRICT vp = &v[0];
    const double xL1 = xL - 1.0;

    integrate_vector_impl<0>(ep, xL, xL1, weight, jsx, vp);
  }


  XARA_FORCE_INTERP_INLINE static void
  integrate(const MatrixND<nsr,nsr>& Fs,
            double xL, double wtL, double jsx,
            MatrixND<NBV,NBV>& F) 
  noexcept
  {
    const double* X_RESTRICT f = Fs.data();
    double*       X_RESTRICT g = F.data();
    const double xL1 = xL - 1.0;
    static_assert(scheme_is_unique<0,1>(), 
                  "stress components in scheme must be unique.");

    if constexpr (scheme_is_unique<0,1>()) {
      integrate_matrix_output_impl<0,0>(f, xL, xL1, wtL, jsx, g);
    }
    else {
      integrate_matrix_scatter_impl<0,0>(f, xL, xL1, wtL, jsx, g);
    }
  }

  //
  // Sensitivity
  //
  const MatrixND<nsr,NBV> &
  getDB(double dxLdh, double d1oLdh)
  {
    B.zero();
    for (int i = 0; i < nsr; i++) {
      switch (scheme[i]) {
      case FrameStress::N:
        B(i, jnx) = 0.0;
        break;
      case FrameStress::Vy:
        B(i, imz) = -d1oLdh;
        B(i, jmz) = -d1oLdh;
        break;
      case FrameStress::Vz:
        B(i, imy) =  d1oLdh;
        B(i, jmy) =  d1oLdh;
        break;
      case FrameStress::T:
        B(i, jmx) = 0.0;
        break;
      case FrameStress::My:
        B(i, imy) = dxLdh;
        B(i, jmy) = dxLdh;
        break;
      case FrameStress::Mz:
        B(i, imz) = dxLdh;
        B(i, jmz) = dxLdh;
        break;
      case FrameStress::Bimoment:
        B(i, iwx) = dxLdh;
        B(i, jwx) = dxLdh;
        break;
      case FrameStress::Bishear:
        B(i, iwx) = d1oLdh;
        B(i, jwx) = d1oLdh;
        break;
      }
    }
    return B;
  }

private:

  template <int w>
  XARA_FORCE_INTERP_INLINE static void
  expand_vector_warp(const double* X_RESTRICT q, double* X_RESTRICT p) noexcept
  {
    if constexpr (w < nwm) {
      p[6+w]       = q[NNW + 2*w];
      p[NDF+6+w]   = q[NNW + 2*w + 1];
      expand_vector_warp<w+1>(q, p);
    }
  }

  template <int w>
  XARA_FORCE_INTERP_INLINE static void
  expand_matrix_column_warp_rows(const double* X_RESTRICT kc,
                                 double* X_RESTRICT Kc,
                                 double col_sign) noexcept
  {
    if constexpr (w < nwm) {
      Kc[6+w]       = col_sign * kc[NNW + 2*w];
      Kc[NDF+6+w]   = col_sign * kc[NNW + 2*w + 1];
      expand_matrix_column_warp_rows<w+1>(kc, Kc, col_sign);
    }
  }

  XARA_FORCE_INTERP_INLINE static void
  expand_matrix_column(double* X_RESTRICT K, 
                       const double* X_RESTRICT k,
                       int out_col, int src_col, 
                       double col_sign) noexcept
  {
    const double* X_RESTRICT kc = k + src_col*NBV;
    double*       X_RESTRICT Kc = K + out_col*NLF;

    Kc[0]       = -col_sign * kc[jnx];
    Kc[NDF]     =  col_sign * kc[jnx];
    Kc[3]       = -col_sign * kc[jmx];
    Kc[NDF+3]   =  col_sign * kc[jmx];
    Kc[4]       =  col_sign * kc[imy];
    Kc[NDF+4]   =  col_sign * kc[jmy];
    Kc[5]       =  col_sign * kc[imz];
    Kc[NDF+5]   =  col_sign * kc[jmz];

    expand_matrix_column_warp_rows<0>(kc, Kc, col_sign);
  }

  template <int w>
  XARA_FORCE_INTERP_INLINE static void
  expand_matrix_warp_columns(double* X_RESTRICT K,
                             const double* X_RESTRICT k) noexcept
  {
    if constexpr (w < nwm) {
      expand_matrix_column(K, k, 6+w,       NNW + 2*w,     1.0);
      expand_matrix_column(K, k, NDF+6+w,   NNW + 2*w + 1, 1.0);
      expand_matrix_warp_columns<w+1>(K, k);
    }
  }

  template <int S>
  XARA_FORCE_INTERP_INLINE static void
  b_case(double* X_RESTRICT Bd, int row,
         double xL, double xL1, double jsx) noexcept
  {
    if constexpr (S == FrameStress::N) {
      Bd[jnx*nsr + row] = 1.0;
    }
    else if constexpr (S == FrameStress::Vy) {
      Bd[imz*nsr + row] = -jsx;
      Bd[jmz*nsr + row] = -jsx;
    }
    else if constexpr (S == FrameStress::Vz) {
      Bd[imy*nsr + row] =  jsx;
      Bd[jmy*nsr + row] =  jsx;
    }
    else if constexpr (S == FrameStress::T) {
      Bd[jmx*nsr + row] = 1.0;
    }
    else if constexpr (S == FrameStress::My) {
      Bd[imy*nsr + row] = xL1;
      Bd[jmy*nsr + row] = xL;
    }
    else if constexpr (S == FrameStress::Mz) {
      Bd[imz*nsr + row] = xL1;
      Bd[jmz*nsr + row] = xL;
    }
    else if constexpr (S == FrameStress::Bimoment) {
      Bd[iwx*nsr + row] = xL1;
      Bd[jwx*nsr + row] = xL;
    }
    else if constexpr (S == FrameStress::Bishear) {
      Bd[iwx*nsr + row] = jsx;
      Bd[jwx*nsr + row] = jsx;
    }
  }

  template <int i>
  XARA_FORCE_INTERP_INLINE static void
  b_impl(double* X_RESTRICT Bd, double xL, double xL1, double jsx) noexcept
  {
    if constexpr (i < nsr) {
      b_case<scheme[i]>(Bd, i, xL, xL1, jsx);
      b_impl<i+1>(Bd, xL, xL1, jsx);
    }
  }

  template <int S>
  XARA_FORCE_INTERP_INLINE static double
  interpolate_case(const double* X_RESTRICT q,
                   double xL, double xL1, double jsx) noexcept
  {
    if constexpr (S == FrameStress::N)
      return q[jnx];

    else if constexpr (S == FrameStress::Vy)
      return -jsx * (q[imz] + q[jmz]);

    else if constexpr (S == FrameStress::Vz)
      return  jsx * (q[imy] + q[jmy]);

    else if constexpr (S == FrameStress::T)
      return q[jmx];

    else if constexpr (S == FrameStress::My)
      return xL1*q[imy] + xL*q[jmy];

    else if constexpr (S == FrameStress::Mz)
      return xL1*q[imz] + xL*q[jmz];

    else if constexpr (S == FrameStress::Bimoment)
      return xL1*q[iwx] + xL*q[jwx];

    else if constexpr (S == FrameStress::Bishear)
      return jsx * (q[iwx] + q[jwx]);

    else
      return 0.0;

  }

  template <int i>
  XARA_FORCE_INTERP_INLINE static void
  interpolate_impl(const double* X_RESTRICT q, double xL, double xL1, double jsx,
                   double* X_RESTRICT s) noexcept
  {
    if constexpr (i < nsr) {
      s[i] += interpolate_case<scheme[i]>(q, xL, xL1, jsx);
      interpolate_impl<i+1>(q, xL, xL1, jsx, s);
    }
  }

  template <int S>
  XARA_FORCE_INTERP_INLINE static void
  integrate_vector_case(double* X_RESTRICT v, double ew,
                        double xL, double xL1, double jsx) noexcept
  {
    if constexpr (S == FrameStress::N) {
      v[jnx] += ew;
    }
    else if constexpr (S == FrameStress::Vy) {
      const double a = -jsx*ew;
      v[imz] += a;
      v[jmz] += a;
    }
    else if constexpr (S == FrameStress::Vz) {
      const double a = jsx*ew;
      v[imy] += a;
      v[jmy] += a;
    }
    else if constexpr (S == FrameStress::T) {
      v[jmx] += ew;
    }
    else if constexpr (S == FrameStress::My) {
      v[imy] += xL1*ew;
      v[jmy] += xL *ew;
    }
    else if constexpr (S == FrameStress::Mz) {
      v[imz] += xL1*ew;
      v[jmz] += xL *ew;
    }
    else if constexpr (S == FrameStress::Bimoment) {
      v[iwx] += xL1*ew;
      v[jwx] += xL *ew;
    }
    else if constexpr (S == FrameStress::Bishear) {
      const double a = jsx*ew;
      v[iwx] += a;
      v[jwx] += a;
    }
  }

  template <int i>
  XARA_FORCE_INTERP_INLINE static void
  integrate_vector_impl(const double* X_RESTRICT e,
                        double xL, double xL1, double weight, double jsx,
                        double* X_RESTRICT v) noexcept
  {
    if constexpr (i < nsr) {
      integrate_vector_case<scheme[i]>(v, e[i]*weight, xL, xL1, jsx);
      integrate_vector_impl<i+1>(e, xL, xL1, weight, jsx, v);
    }
  }

  template <int i, int j>
  XARA_FORCE_INTERP_INLINE static consteval bool
  scheme_is_unique() noexcept
  {
    if constexpr (i >= nsr)
      return true;

    else if constexpr (j >= nsr)
      return scheme_is_unique<i+1, i+2>();

    else if constexpr (scheme[i] == scheme[j])
      return false;

    else
      return scheme_is_unique<i, j+1>();
  }

  template <int S, int i = 0>
  XARA_FORCE_INTERP_INLINE static consteval int
  stress_index() noexcept
  {
    if constexpr (i >= nsr)
      return -1;
  
    else if constexpr (scheme[i] == S)
      return i;

    else
      return stress_index<S, i+1>();
  }

  template <int P>
  XARA_FORCE_INTERP_INLINE static consteval int
  basic_term_count() noexcept
  {
    if constexpr (P == jnx || P == jmx)
      return 1;

    else if constexpr (P == imy || P == jmy ||
                       P == imz || P == jmz ||
                       P == iwx || P == jwx)
      return 2;

    else
      return 0;
  }

  template <int P, int term>
  XARA_FORCE_INTERP_INLINE static consteval FrameStress
  basic_term_stress() noexcept
  {
    if constexpr (P == jnx) {
      return FrameStress::N;
    }
    else if constexpr (P == jmx) {
      return FrameStress::T;
    }
    else if constexpr (P == imy || P == jmy) {
      if constexpr (term == 0)
        return FrameStress::Vz;
      else
        return FrameStress::My;
    }
    else if constexpr (P == imz || P == jmz) {
      if constexpr (term == 0)
        return FrameStress::Vy;
      else
        return FrameStress::Mz;
    }
    else if constexpr (P == iwx || P == jwx) {
      if constexpr (term == 0)
        return FrameStress::Bimoment;
      else
        return FrameStress::Bishear;
    }
    else {
      return FrameStress::N;
    }
  }

  template <int P, int term>
  XARA_FORCE_INTERP_INLINE static double
  basic_term_coeff(double xL, double xL1, double jsx) noexcept
  {
    if constexpr (P == jnx || P == jmx) {
      return 1.0;
    }
    else if constexpr (P == imy || P == jmy) {
      if constexpr (term == 0)
        return jsx;
      else if constexpr (P == imy)
        return xL1;
      else
        return xL;
    }
    else if constexpr (P == imz || P == jmz) {
      if constexpr (term == 0)
        return -jsx;
      else if constexpr (P == imz)
        return xL1;
      else
        return xL;
    }
    else if constexpr (P == iwx || P == jwx) {
      if constexpr (term == 0) {
        if constexpr (P == iwx)
          return xL1;
        else
          return xL;
      }
      else {
        return jsx;
      }
    }
    else {
      return 0.0;
    }
  }

  template <int P, int Q, int ip, int iq>
  XARA_FORCE_INTERP_INLINE static void
  add_basic_matrix_term(double& sum, const double* X_RESTRICT Fs,
                        double xL, double xL1, double jsx) noexcept
  {
    if constexpr (ip < basic_term_count<P>() && iq < basic_term_count<Q>()) {
      constexpr FrameStress Sp = basic_term_stress<P, ip>();
      constexpr FrameStress Sq = basic_term_stress<Q, iq>();
      constexpr int row = stress_index<Sp>();
      constexpr int col = stress_index<Sq>();

      if constexpr (row >= 0 && col >= 0) {
        const double cp = basic_term_coeff<P, ip>(xL, xL1, jsx);
        const double cq = basic_term_coeff<Q, iq>(xL, xL1, jsx);
        sum += cp * Fs[col*nsr + row] * cq;
      }
    }
  }

  template <int P, int Q>
  XARA_FORCE_INTERP_INLINE static void
  add_basic_matrix_entry(const double* X_RESTRICT Fs,
                         double xL, double xL1, double wtL, double jsx,
                         double* X_RESTRICT F)
  noexcept
  {
    if constexpr (basic_term_count<P>() != 0 && basic_term_count<Q>() != 0) {
      double sum = 0.0;

      add_basic_matrix_term<P,Q,0,0>(sum, Fs, xL, xL1, jsx);
      add_basic_matrix_term<P,Q,0,1>(sum, Fs, xL, xL1, jsx);
      add_basic_matrix_term<P,Q,1,0>(sum, Fs, xL, xL1, jsx);
      add_basic_matrix_term<P,Q,1,1>(sum, Fs, xL, xL1, jsx);

      F[Q*NBV + P] += wtL * sum;
    }
  }

  template <int P, int Q>
  XARA_FORCE_INTERP_INLINE static void
  integrate_matrix_output_impl(const double* X_RESTRICT Fs,
                               double xL, double xL1, double wtL, double jsx,
                               double* X_RESTRICT F) noexcept
  {
    if constexpr (Q < NBV) {
      if constexpr (P < NBV) {
        add_basic_matrix_entry<P,Q>(Fs, xL, xL1, wtL, jsx, F);
        integrate_matrix_output_impl<P+1,Q>(Fs, xL, xL1, wtL, jsx, F);
      }
      else {
        integrate_matrix_output_impl<0,Q+1>(Fs, xL, xL1, wtL, jsx, F);
      }
    }
  }


  template <int SRight>
  XARA_FORCE_INTERP_INLINE static void
  add_matrix_right_scatter(double* X_RESTRICT F, int left, double left_val,
                           double xL, double xL1, double jsx)
  noexcept
  {
    if constexpr (SRight == FrameStress::N) {
      F[jnx*NBV + left] += left_val;
    }
    else if constexpr (SRight == FrameStress::Vy) {
      const double a = -jsx*left_val;
      F[imz*NBV + left] += a;
      F[jmz*NBV + left] += a;
    }
    else if constexpr (SRight == FrameStress::Vz) {
      const double a = jsx*left_val;
      F[imy*NBV + left] += a;
      F[jmy*NBV + left] += a;
    }
    else if constexpr (SRight == FrameStress::T) {
      F[jmx*NBV + left] += left_val;
    }
    else if constexpr (SRight == FrameStress::My) {
      F[imy*NBV + left] += xL1*left_val;
      F[jmy*NBV + left] += xL *left_val;
    }
    else if constexpr (SRight == FrameStress::Mz) {
      F[imz*NBV + left] += xL1*left_val;
      F[jmz*NBV + left] += xL *left_val;
    }
    else if constexpr (SRight == FrameStress::Bimoment) {
      F[iwx*NBV + left] += xL1*left_val;
      F[jwx*NBV + left] += xL *left_val;
    }
    else if constexpr (SRight == FrameStress::Bishear) {
      const double a = jsx*left_val;
      F[iwx*NBV + left] += a;
      F[jwx*NBV + left] += a;
    }
  }

  template <int SLeft, int SRight>
  XARA_FORCE_INTERP_INLINE static void
  add_matrix_pair_scatter(double* X_RESTRICT F, double a,
                          double xL, double xL1, double jsx)
  noexcept
  {
    if constexpr (SLeft == FrameStress::N) {
      add_matrix_right_scatter<SRight>(F, jnx, a, xL, xL1, jsx);
    }
    else if constexpr (SLeft == FrameStress::Vy) {
      const double b = -jsx*a;
      add_matrix_right_scatter<SRight>(F, imz, b, xL, xL1, jsx);
      add_matrix_right_scatter<SRight>(F, jmz, b, xL, xL1, jsx);
    }
    else if constexpr (SLeft == FrameStress::Vz) {
      const double b = jsx*a;
      add_matrix_right_scatter<SRight>(F, imy, b, xL, xL1, jsx);
      add_matrix_right_scatter<SRight>(F, jmy, b, xL, xL1, jsx);
    }
    else if constexpr (SLeft == FrameStress::T) {
      add_matrix_right_scatter<SRight>(F, jmx, a, xL, xL1, jsx);
    }
    else if constexpr (SLeft == FrameStress::My) {
      add_matrix_right_scatter<SRight>(F, imy, xL1*a, xL, xL1, jsx);
      add_matrix_right_scatter<SRight>(F, jmy, xL *a, xL, xL1, jsx);
    }
    else if constexpr (SLeft == FrameStress::Mz) {
      add_matrix_right_scatter<SRight>(F, imz, xL1*a, xL, xL1, jsx);
      add_matrix_right_scatter<SRight>(F, jmz, xL *a, xL, xL1, jsx);
    }
    else if constexpr (SLeft == FrameStress::Bimoment) {
      add_matrix_right_scatter<SRight>(F, iwx, xL1*a, xL, xL1, jsx);
      add_matrix_right_scatter<SRight>(F, jwx, xL *a, xL, xL1, jsx);
    }
    else if constexpr (SLeft == FrameStress::Bishear) {
      const double b = jsx*a;
      add_matrix_right_scatter<SRight>(F, iwx, b, xL, xL1, jsx);
      add_matrix_right_scatter<SRight>(F, jwx, b, xL, xL1, jsx);
    }
  }

  template <int i, int j>
  XARA_FORCE_INTERP_INLINE static void
  integrate_matrix_scatter_impl(const double* X_RESTRICT Fs,
                                double xL, double xL1, double wtL, double jsx,
                                double* X_RESTRICT F) noexcept
  {
    if constexpr (i < nsr) {
      if constexpr (j < nsr) {
        const double a = Fs[j*nsr + i] * wtL;
        add_matrix_pair_scatter<scheme[i], scheme[j]>(F, a, xL, xL1, jsx);
        integrate_matrix_scatter_impl<i, j+1>(Fs, xL, xL1, wtL, jsx, F);
      }
      else {
        integrate_matrix_scatter_impl<i+1, 0>(Fs, xL, xL1, wtL, jsx, F);
      }
    }
  }

private:
  MatrixND<nsr,NBV> B;
};

#undef X_RESTRICT
#undef XARA_FORCE_INTERP_INLINE
