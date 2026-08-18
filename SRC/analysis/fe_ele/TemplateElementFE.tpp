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
#ifndef TemplateElementFE_tpp
#define TemplateElementFE_tpp

#include <analysis/model/AnalysisModel.h>
#include <DOF_Group.h>
#include <Domain.h>
#include <Element.h>
#include <Integrator.h>
#include <Logging.h>
#include <Matrix.h>
#include <Node.h>
#include <Vector.h>

#include <cassert>

template <int N>
TemplateElementFE<N>::TemplateElementFE(int tag, Element &element)
  : FE_Element(tag, element.getExternalNodes().Size(), element.getNumDOF())
  , myID(N)
  , myEle(element)
  , theIntegrator(nullptr)
{
  assert(myEle.getNumDOF() == N);
  assert(!myEle.isSubdomain());

  Domain *domain = myEle.getDomain();
  assert(domain != nullptr);

  const int numGroups = myEle.getNumExternalNodes();
  const ID &nodes = myEle.getExternalNodes();
  for (int i = 0; i < numGroups; ++i) {
    Node *node = domain->getNode(nodes(i));
    assert(node != nullptr);

    DOF_Group *dofGroup = node->getDOF_GroupPtr();
    assert(dofGroup != nullptr);
    myDOF_Groups(i) = dofGroup->getTag();
  }
}

template <int N>
Matrix &
TemplateElementFE<N>::tangentView()
{
  // Matrix is a non-owning view: its destructor must not release the
  // MatrixND storage backing it.
  static Matrix view(tangentStorage.data(), N, N);
  return view;
}

template <int N>
Vector &
TemplateElementFE<N>::residualView()
{
  // Vector is a non-owning view of the class-wide VectorND workspace.
  static Vector view(residualStorage.values, N);
  return view;
}

template <int N>
Vector &
TemplateElementFE<N>::scratchView()
{
  // Vector is a non-owning view of the class-wide VectorND workspace.
  static Vector view(scratchStorage.values, N);
  return view;
}

template <int N>
const ID &
TemplateElementFE<N>::getID() const
{
  return myID;
}

template <int N>
int
TemplateElementFE<N>::setID(AnalysisModel &model)
{
  int current = 0;

  const int numGroups = myDOF_Groups.Size();
  for (int i = 0; i < numGroups; ++i) {
    DOF_Group *dofGroup = model.getDOF_GroupPtr(myDOF_Groups(i));
    assert(dofGroup != nullptr);

    const ID &dofID = dofGroup->getID();
    for (int j = 0; j < dofID.Size(); ++j) {
      assert(current < N);
      myID(current++) = dofID(j);
    }
  }
  return 0;
}

template <int N>
int
TemplateElementFE<N>::updateElement()
{
  return myEle.update();
}

template <int N>
const Matrix &
TemplateElementFE<N>::getTangent(Integrator *newIntegrator)
{
  theIntegrator = newIntegrator;

  if (newIntegrator != nullptr)
    newIntegrator->formEleTangent(this);

  return tangentView();
}

template <int N>
void
TemplateElementFE<N>::zeroTangent()
{
  tangentStorage.Zero();
}

template <int N>
void
TemplateElementFE<N>::addKtToTang(double fact)
{
  if (fact != 0.0)
    tangentStorage.addMatrix(myEle.getTangentStiff(), fact);
}

template <int N>
void
TemplateElementFE<N>::addCtoTang(double fact)
{
  if (fact != 0.0)
    tangentStorage.addMatrix(myEle.getDamp(), fact);
}

template <int N>
void
TemplateElementFE<N>::addMtoTang(double fact)
{
  if (fact != 0.0)
    tangentStorage.addMatrix(myEle.getMass(), fact);
}

template <int N>
void
TemplateElementFE<N>::addKiToTang(double fact)
{
  if (fact != 0.0)
    tangentStorage.addMatrix(myEle.getInitialStiff(), fact);
}

template <int N>
void
TemplateElementFE<N>::addKpToTang(double fact, int numP)
{
  if (fact == 0.0)
    return;

  const Matrix *previous = myEle.getPreviousK(numP);
  if (previous != nullptr)
    tangentStorage.addMatrix(*previous, fact);
}

template <int N>
int
TemplateElementFE<N>::storePreviousK(int numP)
{
  return myEle.storePreviousK(numP);
}

template <int N>
const Vector &
TemplateElementFE<N>::getResidual(Integrator *newIntegrator)
{
  theIntegrator = newIntegrator;

  if (newIntegrator != nullptr)
    newIntegrator->formEleResidual(this);

  return residualView();
}

template <int N>
void
TemplateElementFE<N>::zeroResidual()
{
  residualStorage.zero();
}

template <int N>
void
TemplateElementFE<N>::addRtoResidual(double fact)
{
  if (fact != 0.0)
    residualStorage.addVector(1.0, myEle.getResistingForce(), -fact);
}

template <int N>
void
TemplateElementFE<N>::addRIncInertiaToResidual(double fact)
{
  if (fact != 0.0)
    residualStorage.addVector(1.0, myEle.getResistingForceIncInertia(), -fact);
}

template <int N>
void
TemplateElementFE<N>::gatherResponse(OpenSees::VectorND<N> &target,
                                     const Vector &response) const
{
  for (int i = 0; i < N; ++i) {
    const int dof = myID(i);
    target(i) = dof >= 0 ? response(dof) : 0.0;
  }
}

template <int N>
const Vector &
TemplateElementFE<N>::getTangForce(const Vector &disp, double fact)
{
  residualStorage.zero();
  if (fact == 0.0)
    return residualView();

  OpenSees::VectorND<N> temporary{};
  gatherResponse(temporary, disp);
  theIntegrator->formEleTangent(this);
  residualStorage.addMatrixVector(tangentStorage, temporary, fact);

  return residualView();
}

template <int N>
const Vector &
TemplateElementFE<N>::getK_Force(const Vector &disp, double fact)
{
  residualStorage.zero();
  if (fact == 0.0)
    return residualView();

  OpenSees::VectorND<N> temporary{};
  gatherResponse(temporary, disp);
  Vector temporaryView(temporary.values, N);
  residualStorage.addMatrixVector(1.0, myEle.getTangentStiff(), temporaryView, fact);

  return residualView();
}

template <int N>
const Vector &
TemplateElementFE<N>::getKi_Force(const Vector &disp, double fact)
{
  residualStorage.zero();
  if (fact == 0.0)
    return residualView();

  OpenSees::VectorND<N> temporary{};
  gatherResponse(temporary, disp);
  Vector temporaryView(temporary.values, N);
  residualStorage.addMatrixVector(1.0, myEle.getInitialStiff(), temporaryView, fact);

  return residualView();
}

template <int N>
const Vector &
TemplateElementFE<N>::getM_Force(const Vector &accel, double fact)
{
  residualStorage.zero();
  if (fact == 0.0)
    return residualView();

  gatherResponse(scratchStorage, accel);
  residualStorage.addMatrixVector(1.0, myEle.getMass(), scratchView(), fact);

  return residualView();
}

template <int N>
const Vector &
TemplateElementFE<N>::getC_Force(const Vector &veloc, double fact)
{
  residualStorage.zero();
  if (fact == 0.0)
    return residualView();

  gatherResponse(scratchStorage, veloc);
  residualStorage.addMatrixVector(1.0, myEle.getDamp(), scratchView(), fact);

  return residualView();
}

template <int N>
Integrator *
TemplateElementFE<N>::getLastIntegrator()
{
  return theIntegrator;
}

template <int N>
const Vector &
TemplateElementFE<N>::getLastResponse()
{
  if (theIntegrator != nullptr) {
    if (theIntegrator->getLastResponse(residualView(), myID) < 0) {
      opserr << "WARNING TemplateElementFE::getLastResponse()";
      opserr << " - the Integrator had problems with getLastResponse()\n";
    }
  }
  else {
    residualStorage.zero();
    opserr << "WARNING TemplateElementFE::getLastResponse()";
    opserr << " No Integrator yet passed\n";
  }

  return residualView();
}

template <int N>
void
TemplateElementFE<N>::addM_Force(const Vector &accel, double fact)
{
  if (fact == 0.0)
    return;

  gatherResponse(scratchStorage, accel);
  residualStorage.addMatrixVector(1.0, myEle.getMass(), scratchView(), fact);
}

template <int N>
void
TemplateElementFE<N>::addD_Force(const Vector &veloc, double fact)
{
  if (fact == 0.0)
    return;

  OpenSees::VectorND<N> temporary{};
  gatherResponse(temporary, veloc);
  Vector temporaryView(temporary.values, N);
  residualStorage.addMatrixVector(1.0, myEle.getDamp(), temporaryView, fact);
}

template <int N>
void
TemplateElementFE<N>::addK_Force(const Vector &disp, double fact)
{
  if (fact == 0.0)
    return;

  OpenSees::VectorND<N> temporary{};
  gatherResponse(temporary, disp);
  Vector temporaryView(temporary.values, N);
  residualStorage.addMatrixVector(1.0, myEle.getTangentStiff(), temporaryView, fact);
}

template <int N>
void
TemplateElementFE<N>::addLocalM_Force(const Vector &accel, double fact)
{
  if (fact != 0.0)
    residualStorage.addMatrixVector(1.0, myEle.getMass(), accel, fact);
}

template <int N>
void
TemplateElementFE<N>::addLocalD_Force(const Vector &veloc, double fact)
{
  if (fact != 0.0)
    residualStorage.addMatrixVector(1.0, myEle.getDamp(), veloc, fact);
}

template <int N>
Element *
TemplateElementFE<N>::getElement()
{
  return &myEle;
}

template <int N>
void
TemplateElementFE<N>::addResistingForceSensitivity(int gradNumber, double fact)
{
  residualStorage.addVector(1.0,
                            myEle.getResistingForceSensitivity(gradNumber),
                            -fact);
}

template <int N>
void
TemplateElementFE<N>::addM_ForceSensitivity(int gradNumber,
                                             const Vector &vect,
                                             double fact)
{
  OpenSees::VectorND<N> temporary{};
  gatherResponse(temporary, vect);
  Vector temporaryView(temporary.values, N);
  residualStorage.addMatrixVector(1.0, myEle.getMassSensitivity(gradNumber),
                                  temporaryView, fact);
}

template <int N>
void
TemplateElementFE<N>::addD_ForceSensitivity(int gradNumber,
                                             const Vector &vect,
                                             double fact)
{
  if (fact == 0.0)
    return;

  OpenSees::VectorND<N> temporary{};
  gatherResponse(temporary, vect);
  Vector temporaryView(temporary.values, N);
  residualStorage.addMatrixVector(1.0, myEle.getDampSensitivity(gradNumber),
                                  temporaryView, fact);
}

template <int N>
void
TemplateElementFE<N>::addLocalD_ForceSensitivity(int gradNumber,
                                                  const Vector &veloc,
                                                  double fact)
{
  if (fact != 0.0)
    residualStorage.addMatrixVector(1.0,
                                    myEle.getDampSensitivity(gradNumber),
                                    veloc, fact);
}

template <int N>
void
TemplateElementFE<N>::addLocalM_ForceSensitivity(int gradNumber,
                                                  const Vector &accel,
                                                  double fact)
{
  if (fact != 0.0)
    residualStorage.addMatrixVector(1.0,
                                    myEle.getMassSensitivity(gradNumber),
                                    accel, fact);
}

template <int N>
int
TemplateElementFE<N>::commitSensitivity(int gradNum, int numGrads)
{
  myEle.commitSensitivity(gradNum, numGrads);
  return 0;
}

#endif
