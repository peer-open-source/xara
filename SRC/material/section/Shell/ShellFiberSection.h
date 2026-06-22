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

#include <vector>
#include <tuple>

#include <Vector.h>
#include <Matrix.h>
#include <ID.h>
#include <NDMaterial.h>
#include <SectionForceDeformation.h>

class Channel;
class FEM_ObjectBroker;
class Information;
class OPS_Stream;
class Response;

class ShellFiberSection : public SectionForceDeformation {
public:
  using Fiber = std::tuple<double, double, NDMaterial*>; // z, weight, material
  using FiberVector = std::vector<Fiber>;

  ShellFiberSection();

  ShellFiberSection(int tag,
                    const FiberVector& fibers);

  // Compatibility constructor for the old layered-shell representation.
  ShellFiberSection(int tag,
                    int nLayers,
                    double* thickness,
                    NDMaterial** materials);

  ShellFiberSection(const ShellFiberSection& other);
  ShellFiberSection& operator=(const ShellFiberSection&) = delete;

  ~ShellFiberSection();

  const char* getClassType() const
  {
    return "ShellFiberSection";
  }

  SectionForceDeformation* getCopy();

  double getRho();

  int getOrder() const;
  const ID& getType();

  int commitState();
  int revertToLastCommit();
  int revertToStart();

  int setTrialSectionDeformation(const Vector& strain_from_element);

  const Vector& getSectionDeformation();
  const Vector& getStressResultant();
  const Matrix& getSectionTangent();

  const Matrix& getInitialTangent()
  {
    return this->getSectionTangent();
  }

  Response* setResponse(const char** argv, int argc, OPS_Stream& s);
  int getResponse(int responseID, Information& info);

  void Print(OPS_Stream& s, int flag);

  int sendSelf(int commitTag, Channel&);
  int recvSelf(int commitTag, Channel&, FEM_ObjectBroker&);

private:
  void clearFibers();
  void copyFiber(double z, double w, NDMaterial* material);
  int getNumFibers() const;

private:
  FiberVector fibers;  // owned PlateFiber material copies
  double h;            // total section thickness, equal to sum of weights

  Vector strainResultant;

  static Vector stressResultant;
  static Matrix tangent;
  static ID array;
};
