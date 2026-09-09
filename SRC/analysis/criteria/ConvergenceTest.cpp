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
// Purpose: This file contains the class definition for ConvergenceTest,
// which is an abstract class. Objects of concrete subclasses can be used
// to test the convergence of an algorithm.
//
// Written: fmk
// Date: 09/98
//
#include <ConvergenceTest.h>
#include <string>
#include <stdlib.h>
#include <iomanip>
#include <sstream>
 

ConvergenceTest::ConvergenceTest(int clasTag)
:MovableObject(clasTag)
{

}

ConvergenceTest::~ConvergenceTest()
{

}

std::string
ConvergenceTest::pad(double x)
{
  std::ostringstream oss;
  oss << std::setw(11) << x;
  return oss.str();
}


std::string
ConvergenceTest::pad(int i)
{
  const int n = 5; // i < 100 ? 3 : 5;
  std::ostringstream oss;
  oss << std::setw(n) << i;
  return oss.str();
}
