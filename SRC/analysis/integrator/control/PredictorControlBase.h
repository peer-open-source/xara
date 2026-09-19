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
// PredictorControlBase selects the load-factor increment at the beginning of
// a load step.
//
#pragma once

class LinearSOE;
class AnalysisModel;
class Vector;

class PredictorControlBase {
public:
  enum class Type {
    IterationTarget,   // -j
    StiffnessParameter
  };
  virtual int domainChanged(LinearSOE&,  AnalysisModel&) { return 0; }

  // Select the increment for a new step and initialize any trial history for
  // that step.
  virtual int predict(LinearSOE &, double &dLam) = 0;

  virtual int update(const Vector &dX) = 0;

  virtual void commit() {}
  virtual void revert() {}
};
