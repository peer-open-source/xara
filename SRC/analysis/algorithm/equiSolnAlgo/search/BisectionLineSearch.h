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
// Description: This file contains the class definition for BisectinLineSearch.
// This performs the search for U(i+1) = U(i) + eta * deltaU(i) by using the 
// bisection method to find the best solution.
//
//                eta(j+1) = eta(l) + eta(u)
//                           ---------------
//                                2.0
//
// where     s(j) = U(i+1,j) ^ R(U(i+1, j))
//
//  and      U(i+1,j) = U(i) + eta(j)*deltaU(i)
//
// note eta(u) and eta(l) must bracket the solution, i.e. s(u)*s(l)<0,
//
//      if s(eta(j+1))*s(l) < 0 { eta(u) = eta(j+1) and s(u) = s(eta(j+1))
//      if s(eta(j+1))*s(u) < 0 { eta(l) = eta(j+1) and s(l) = s(eta(j+1))
//      if s(eta(j+1))*s(u) == 0  SOLN FOUND.
//
// Written: cmp
// Adapted from BisectionLineSearch.cpp by fmk, dated 11/01
//
#pragma once

#include <LineSearch.h>
#include <Vector.h>

class BisectionLineSearch: public LineSearch
{
  public:
    BisectionLineSearch(double tolerance = 0.8, 
                        int maxIter = 10, 
                        double minEta = 0.1, 
                        double maxEta = 10.0, 
                        int flag = 1);

    ~BisectionLineSearch();

    int newStep(const Vector &) override;
    int search(double s0, 
               double s1, 
               const Vector& dU,
               Vector& G,
               Vector& Xs,
               IncrementalResidual &) override;

    void Print(OPS_Stream &s, int flag) override;    

  private:
    double tolerance;
    int    maxIter;
    double minEta;
    double maxEta;
    int    printFlag;
};
