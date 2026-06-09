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
// Claudio M. Perez
//
#include <stdio.h>
#include <stdlib.h>
#include "GaussLegendre1D.hpp"

template<template<int, int> class Gauss> void
quadrature(int n, const double **const x, const double **const w)
{
  switch (n) {
    case  1: *x = Gauss<1, 1>::pts; *w = Gauss<1, 1>::wts; break;
    case  2: *x = Gauss<1, 2>::pts; *w = Gauss<1, 2>::wts; break;
    case  3: *x = Gauss<1, 3>::pts; *w = Gauss<1, 3>::wts; break;
    case  4: *x = Gauss<1, 4>::pts; *w = Gauss<1, 4>::wts; break;
    case  5: *x = Gauss<1, 5>::pts; *w = Gauss<1, 5>::wts; break;
    case  6: *x = Gauss<1, 6>::pts; *w = Gauss<1, 6>::wts; break;
    case  7: *x = Gauss<1, 7>::pts; *w = Gauss<1, 7>::wts; break;
    case  8: *x = Gauss<1, 8>::pts; *w = Gauss<1, 8>::wts; break;
    case  9: *x = Gauss<1, 9>::pts; *w = Gauss<1, 9>::wts; break;
    case 10: *x = Gauss<1,10>::pts; *w = Gauss<1,10>::wts; break;
    case 11: *x = Gauss<1,11>::pts; *w = Gauss<1,11>::wts; break;
    case 12: *x = Gauss<1,12>::pts; *w = Gauss<1,12>::wts; break;
    case 13: *x = Gauss<1,13>::pts; *w = Gauss<1,13>::wts; break;
    case 14: *x = Gauss<1,14>::pts; *w = Gauss<1,14>::wts; break;
    case 15: *x = Gauss<1,15>::pts; *w = Gauss<1,15>::wts; break;
    case 16: *x = Gauss<1,16>::pts; *w = Gauss<1,16>::wts; break;
    case 17: *x = Gauss<1,17>::pts; *w = Gauss<1,17>::wts; break;
    case 18: *x = Gauss<1,18>::pts; *w = Gauss<1,18>::wts; break;
    case 19: *x = Gauss<1,19>::pts; *w = Gauss<1,19>::wts; break;
    case 20: *x = Gauss<1,20>::pts; *w = Gauss<1,20>::wts; break;
  }
}

template<template<int, int> class Gauss, int ndm> void
quadrature(int n, const double **const x[ndm], const double **const w)
{
  switch (n) {
    case  1: *x = Gauss<ndm, 1>::pts; *w = Gauss<ndm, 1>::wts; break;
    case  2: *x = Gauss<ndm, 2>::pts; *w = Gauss<ndm, 2>::wts; break;
    case  3: *x = Gauss<ndm, 3>::pts; *w = Gauss<ndm, 3>::wts; break;
    case  4: *x = Gauss<ndm, 4>::pts; *w = Gauss<ndm, 4>::wts; break;
    case  5: *x = Gauss<ndm, 5>::pts; *w = Gauss<ndm, 5>::wts; break;
    case  6: *x = Gauss<ndm, 6>::pts; *w = Gauss<ndm, 6>::wts; break;
    case  7: *x = Gauss<ndm, 7>::pts; *w = Gauss<ndm, 7>::wts; break;
    case  8: *x = Gauss<ndm, 8>::pts; *w = Gauss<ndm, 8>::wts; break;
    case  9: *x = Gauss<ndm, 9>::pts; *w = Gauss<ndm, 9>::wts; break;
    case 10: *x = Gauss<ndm,10>::pts; *w = Gauss<ndm,10>::wts; break;
    case 11: *x = Gauss<ndm,11>::pts; *w = Gauss<ndm,11>::wts; break;
    case 12: *x = Gauss<ndm,12>::pts; *w = Gauss<ndm,12>::wts; break;
    case 13: *x = Gauss<ndm,13>::pts; *w = Gauss<ndm,13>::wts; break;
    case 14: *x = Gauss<ndm,14>::pts; *w = Gauss<ndm,14>::wts; break;
    case 15: *x = Gauss<ndm,15>::pts; *w = Gauss<ndm,15>::wts; break;
    case 16: *x = Gauss<ndm,16>::pts; *w = Gauss<ndm,16>::wts; break;
    case 17: *x = Gauss<ndm,17>::pts; *w = Gauss<ndm,17>::wts; break;
    case 18: *x = Gauss<ndm,18>::pts; *w = Gauss<ndm,18>::wts; break;
    case 19: *x = Gauss<ndm,19>::pts; *w = Gauss<ndm,19>::wts; break;
    case 20: *x = Gauss<ndm,20>::pts; *w = Gauss<ndm,20>::wts; break;
  }
}

int
main(int argc, char** argv)
{
  const double *x, *w;
  int n = atoi(argv[1]);

  quadrature<GaussLegendre>(n, &x, &w);

  for (int i=0; i<n; i++)
    printf("%lf\t%lf\n", x[i], w[i]);

}

