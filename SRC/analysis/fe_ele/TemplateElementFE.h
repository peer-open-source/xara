//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2026, OpenSees/Xara Developers
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause.
//
//===----------------------------------------------------------------------===//
//
// Description: Fixed-size implementation of the ordinary Element FE adapter.
//
//===----------------------------------------------------------------------===//
#pragma once

#include <FE_Element.h>
#include <MatrixND.h>
#include <VectorND.h>

class AnalysisModel;
class Element;
class Integrator;
class Matrix;
class Vector;


template <int N>
class TemplateElementFE : public FE_Element
{
public:
  static_assert(N > 0, "TemplateElementFE requires at least one DOF");

  TemplateElementFE(int tag, Element &element);
  ~TemplateElementFE() override = default;

  TemplateElementFE(const TemplateElementFE &) = delete;
  TemplateElementFE(TemplateElementFE &&) = delete;
  TemplateElementFE &operator=(const TemplateElementFE &) = delete;
  TemplateElementFE &operator=(TemplateElementFE &&) = delete;

  const char *getClassName() const override { return "TemplateElementFE"; }

  int setID(AnalysisModel &) override;
  const ID &getID() const override;

  const Matrix &getTangent(Integrator *) override;
  const Vector &getResidual(Integrator *) override;

  void zeroTangent() override;
  void addKtToTang(double fact = 1.0) override;
  void addKiToTang(double fact = 1.0) override;
  void addCtoTang(double fact = 1.0) override;
  void addMtoTang(double fact = 1.0) override;
  void addKpToTang(double fact = 1.0, int numP = 0) override;
  int storePreviousK(int numP) override;

  void zeroResidual() override;
  void addRtoResidual(double fact = 1.0) override;
  void addRIncInertiaToResidual(double fact) override;

  const Vector &getTangForce(const Vector &x, double fact = 1.0) override;
  const Vector &getK_Force(const Vector &x, double fact = 1.0) override;
  const Vector &getKi_Force(const Vector &x, double fact = 1.0) override;
  const Vector &getC_Force(const Vector &x, double fact = 1.0) override;
  const Vector &getM_Force(const Vector &x, double fact = 1.0) override;
  void addM_Force(const Vector &accel, double fact = 1.0) override;
  void addD_Force(const Vector &vel, double fact = 1.0) override;
  void addK_Force(const Vector &disp, double fact = 1.0) override;

  virtual int updateElement();

  Integrator *getLastIntegrator() override;
  const Vector &getLastResponse() override;
  Element *getElement() override;

  void Print(OPS_Stream &, int) override {}

  void addResistingForceSensitivity(int gradNumber, double fact = 1.0) override;
  void addM_ForceSensitivity(int gradNumber, const Vector &vect,
                             double fact = 1.0) override;
  void addD_ForceSensitivity(int gradNumber, const Vector &vect,
                             double fact = 1.0) override;
  int commitSensitivity(int gradNum, int numGrads) override;

protected:
  void addLocalM_Force(const Vector &accel, double fact = 1.0);
  void addLocalD_Force(const Vector &vel, double fact = 1.0);
  void addLocalM_ForceSensitivity(int gradNumber, const Vector &accel,
                                  double fact = 1.0);
  void addLocalD_ForceSensitivity(int gradNumber, const Vector &vel,
                                  double fact = 1.0);

  ID myID;

private:
  static Matrix &tangentView();
  static Vector &residualView();
  static Vector &scratchView();
  void gatherResponse(OpenSees::VectorND<N> &target,
                      const Vector &response) const;

  Element &myEle;
  Integrator *theIntegrator;

  inline static OpenSees::MatrixND<N, N> tangentStorage{};
  inline static OpenSees::VectorND<N> residualStorage{};
  inline static OpenSees::VectorND<N> scratchStorage{};
};

#include "TemplateElementFE.tpp"

