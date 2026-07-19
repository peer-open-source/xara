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
#pragma once

#include <cmath>
#include <LinearSOE.h>

struct TrackSign {
  enum class Type {
    Determinant,          // cumulative Clarke-Hancock determinant rule
    DeterminantOpenSees,  // reproduce the OpenSees determinant multiplier
    DeltaLambdaStep       // sign of previous accumulated load step
  };

  explicit TrackSign(Type type)
    : type(type),
      lastSign(1),
      lastDetSign(1),
      last_det(0.0),
      last_dlambda(0.0),
      haveLastDet(false),
      haveLastStep(false)
  {}

  double newStep(LinearSOE& soe, double dlambda) {
    // Return a multiplier m in {-1,+1}.

    const int inputSign = signNonZero(dlambda);
    int multiplier = 1;

    switch (type) {
    case Type::DeltaLambdaStep: {
      // OpenSees SIGN_LAST_STEP logic. 
      // The next predictor follows the sign of the previous accumulated
      // load-step increment.
      const int targetSign =
          haveLastStep && last_dlambda != 0.0
        ? signNonZero(last_dlambda)
        : inputSign;

      multiplier = targetSign * inputSign;
      lastSign = targetSign;
      break;
    }

    case Type::Determinant: {
      // Clarke-Hancock rule:
      // keep the previous predictor sign unless det(K_T) changes sign.
      //
      // Equivalently:
      //
      //     s_i = s_{i-1} sign(det K_i) sign(det K_{i-1})
      //
      // but with a separate initialization for the first step.
      const double det = soe.getDeterminant();
      const int detSign = determinantSign(det);

      if (!haveLastDet) {
        // No previous determinant exists.  Do not invent a reversal.
        lastSign = inputSign;
      } else if (detSign != lastDetSign) {
        lastSign = -lastSign;
      }

      last_det = det;
      lastDetSign = detSign;
      haveLastDet = true;

      multiplier = lastSign * inputSign;
      break;
    }

    case Type::DeterminantOpenSees: {
      // This reproduces the OpenSees code:
      //
      //     dLambda *= sign(det_i) * signLastDeterminant;
      //     signLastDeterminant = sign(det_i);
      //
      // This is not cumulative.  It reverses only on the step where the
      // determinant sign changes, then flips back on the following step if
      // the determinant sign remains changed.
      const double det = soe.getDeterminant();
      const int detSign = openSeesSign(det);
      const int prevDetSign = haveLastDet ? lastDetSign : 1;

      multiplier = detSign * prevDetSign;

      last_det = det;
      lastDetSign = detSign;
      haveLastDet = true;

      lastSign = signNonZero(dlambda * multiplier);
      break;
    }
    }

    last_dlambda = dlambda * multiplier;
    haveLastStep = true;

    return static_cast<double>(multiplier);
  }

  double update(LinearSOE& soe, double ddlambda) {
    (void)soe;

    // Do not change the sign of the MinUnbalDispNorm iteration correction.
    // The iteration strategy already computes the signed ddlambda.
    //
    // This method only records the accumulated load increment for the
    // current load step, so that DeltaLambdaStep can use it at the next
    // newStep().
    constexpr int multiplier = 1;

    last_dlambda += ddlambda * multiplier;
    haveLastStep = true;

    if (type == Type::DeltaLambdaStep && last_dlambda != 0.0)
      lastSign = signNonZero(last_dlambda);

    return static_cast<double>(multiplier);
  }

  double lastStepDeltaLambda() const {
    return last_dlambda;
  }

  double lastDeterminant() const {
    return last_det;
  }

  int currentPredictorSign() const {
    return lastSign;
  }

private:
  static int signNonZero(double x) {
    return x < 0.0 ? -1 : 1;
  }

  static int openSeesSign(double det) {
    // Matches OpenSees: det < 0 gives -1, otherwise +1.
    return det < 0.0 ? -1 : 1;
  }

  int determinantSign(double det) const {
    // For the cumulative determinant rule, exactly zero is not treated as
    // a sign change by itself.  The sign changes when the determinant moves
    // to the opposite side of zero.
    if (det < 0.0) return -1;
    if (det > 0.0) return 1;
    return haveLastDet ? lastDetSign : 1;
  }

private:
  const Type type;

  // For Type::Determinant, this is the cumulative predictor sign.
  // For Type::DeltaLambdaStep, it records the sign of the last accumulated load step.
  int lastSign;

  // Sign and value of the determinant at the last newStep().
  int lastDetSign;
  double last_det;

  // Accumulated signed load increment in the current or most recent step.
  double last_dlambda;

  bool haveLastDet;
  bool haveLastStep;
};
