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
#include <array>

#include <Domain.h>
#include <Element.h>
#include <VectorND.h>
#include <MatrixND.h>
#include <ID.h>
#include <State.h>

class Node;
class Domain;
class Response;
class Rotation;

using namespace OpenSees;


template <int nen, int ndm, int ndf>
class FiniteElement : public Element {
public:
    FiniteElement(int tag, int classtag)
      : Element(tag, classtag),
        connectedExternalNodes(nen),
        parameterID(0),
        e_state(State::None)
    {
      for (int i=0; i<nen; i++)
        theNodes[i] = nullptr;
    }

    FiniteElement(int tag, int classtag, 
                  std::array<int, nen>& nodes, 
                  int mass_flag=0)
      : Element(tag, classtag),
        connectedExternalNodes(nen),
        parameterID(0),
        e_state(State::None)
    {
      for (int i=0; i<nen; i++) {
        connectedExternalNodes(i) = nodes[i];
        theNodes[i] = nullptr;
      }
    }


    // For Element
    const ID& getExternalNodes() final {
      return connectedExternalNodes;
    }

    Node **getNodePtrs() final {return theNodes.data();}
    int  getNumExternalNodes() const final {return nen;}
    int  getNumDOF() final {return nen*ndf;}

    void setDomain(Domain *) final;

    const Vector &
    getResistingForceIncInertia() override
    {
      static VectorND<nen*ndf> P_{0.0};
      static Vector P(P_);

      //
      // 1) static residual
      //
      P = this->getResistingForce(); 
      
      //
      // 2) Rayleigh damping
      //
      // add the damping forces if rayleigh damping
      if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
        P.addVector(1.0, this->getRayleighDampingForces(), 1.0);

      // if (total_mass == 0.0)
      //   return P;

      //
      // 3) Inertial forces
      //
      {
        VectorND<nen*ndf> accel{};
        for (int i=0; i<nen; i++) {
          const Vector& trialAccel = theNodes[i]->getTrialAccel();
          for (int j=0; j<6; j++) {
            accel[i*ndf+j] = trialAccel(j);
          }
        }
        P.addMatrixVector(1.0, this->getMass(), accel, 1.0);
      }
      return P;
    }

    int
    setParameter(const char **argv, int argc, Parameter &param) override
    {
      if (argc < 1)
        return -1;

      // don't do anything if MaterialStageParameter calls this element
      if (strcmp(argv[0], "updateMaterialStage") == 0) {
        return -1;
      }

      return -1;
    }

    int
    updateParameter(int paramID, Information &info) override
    {
      return -1;
    }

    int
    activateParameter(int passedParameterID) override
    {
      parameterID = passedParameterID; 
      return 0;
    }


protected:

#ifdef FEFT
    // Implemented by children
    virtual MatrixND<ndf*nen,ndf*nen> getTangent(State state, int rate) =0;
    virtual VectorND<ndf*nen>         getForce(State state, int rate) =0;
#endif

// TODO: Rename setNodes to setReference
    virtual int setNodes() = 0;

    // Supplied for children
    inline int 
    setState(State state) {

      if ((e_state & state) == state)
        return 0;

      else {
        int status = -1;
        switch (state) {
          case State::Init:
            if ((status = this->setNodes()) == 0)
              e_state |= State::Init;
            break;

          case State::Pres:
            if ((status = this->update()) == 0)
              e_state |= State::Pres;
            break;

          default:
            break;
        }
        return status;
      }
    }

protected:
  //
  std::array<Node*, nen> theNodes;

  ID  connectedExternalNodes;

  int  parameterID;

private:
  State  e_state;
  static MatrixND<ndf*nen,ndf*nen> D; // Damping matrix
  static MatrixND<ndf*nen,ndf*nen> M; // Mass matrix
  static MatrixND<ndf*nen,ndf*nen> K; // Stiffness matrix
};

#include "FiniteElement.tpp"