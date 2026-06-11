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
#ifndef LegendreFixedQuadrilateral_H
#define LegendreFixedQuadrilateral_H

template <int n> struct LegendreFixedQuadrilateral;

template <>
struct LegendreFixedQuadrilateral<4> {
  constexpr static int nip = 4;

  constexpr static double pts[4][2] = {
             {-0.577350269189626, -0.577350269189626},
             { 0.577350269189626, -0.577350269189626},
             { 0.577350269189626,  0.577350269189626},
             {-0.577350269189626,  0.577350269189626}};
  
  constexpr static double wts[4] = {1.0, 1.0, 1.0, 1.0};
};


template <>
struct LegendreFixedQuadrilateral<9> {
  constexpr static int    nip = 9;

  constexpr static double pts[9][2] = {
             { -0.7745966692414834, -0.7745966692414834},
             {  0.7745966692414834, -0.7745966692414834},
             {  0.7745966692414834,  0.7745966692414834},
             { -0.7745966692414834,  0.7745966692414834},
             {  0.0,                -0.7745966692414834},
             {  0.7745966692414834,  0.0               },
             {  0.0,                 0.7745966692414834},
             { -0.7745966692414834,  0.0               },
             {  0.0,                 0.0               }};

  constexpr static double wts[9] = {
               0.30864197530864196,
               0.30864197530864196,
               0.30864197530864196,
               0.30864197530864196,
               0.49382716049382713,
               0.49382716049382713,
               0.49382716049382713,
               0.49382716049382713,
               0.7901234567901234 };
};

#endif // LegendreFixedQuadrilateral_H
