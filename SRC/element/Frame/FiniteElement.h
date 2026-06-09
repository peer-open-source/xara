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
        p_iner(nen*ndf),
        parameterID(0),
        e_state(State::None),
        cMass(1)
    {
      for (int i=0; i<nen; i++)
        theNodes[i] = nullptr;
    }

    FiniteElement(int tag, int classtag, std::array<int, nen>& nodes, int mass_flag)
      : Element(tag, classtag),
        connectedExternalNodes(nen),
        p_iner(nen*ndf),
        parameterID(0),
        e_state(State::None),
        cMass(mass_flag)
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
    void zeroLoad() override {
    // TODO: need to reconcile with BasicFrame3d::zeroLoad()
      p_iner.Zero();
    }

    // addInertia(a)
    // addDamping()
    virtual VectorND<nen*ndf>
    getInertia(VectorND<nen*ndf>& accel) {
      VectorND<nen*ndf> zero{};
      return zero;
    }

    virtual int 
    getLumpedInertia(VectorND<nen*ndf>& m) {
      return -1;
    }

    virtual int addResidual(VectorND<nen*ndf>& R, double c, int flag) {
      return 0;
    }

    virtual int addTangent(MatrixND<nen*ndf,nen*ndf>&k, double c, int flag) {
      if (flag == 1) {
        k.addMatrix(this->getMass(), c);
      }
      if (flag == 2) {
        k.addMatrix(this->getTangentStiff(), c);
      }
      if (flag == 3) {
        k.addMatrix(this->getInitialStiff(), c);
      }
      return 0;
    }

    void
    setDomain(Domain *theDomain) final
    {
      if (theDomain == nullptr) {
        for (int i=0; i<nen; i++)
          theNodes[i] = nullptr;
        return;
      }

      for (int i=0; i<nen; i++) {
        theNodes[i] = theDomain->getNode(connectedExternalNodes(i));
        if (theNodes[i] == nullptr) {
          opserr << "FiniteElement::setDomain  tag: " << this->getTag() << " -- Node " 
                 << connectedExternalNodes(i) << " does not exist\n";
          return;
        }

        if (theNodes[i]->getNumberDOF() != ndf) {
          opserr << "FiniteElement::setDomain  tag: " << this->getTag() << " -- Node " << connectedExternalNodes(i) 
                  << " has incorrect number of DOF\n";
          opserr << " " << theNodes[i]->getNumberDOF() << " should be " << ndf << endln;
          return;
        }
      }

      if (theDomain != nullptr)
        this->Element::link(*theDomain);

      if (this->setState(State::Init) != 0)
        return;

//    if (this->setState(State::Pres) != 0)
//      return;
    }


    const Vector &
    getResistingForceIncInertia()
    {
      static VectorND<nen*ndf> P_{0.0};
      static Vector P(P_);
      P = this->getResistingForce(); 
      
      // add the damping forces if rayleigh damping
      if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
        P.addVector(1.0, this->getRayleighDampingForces(), 1.0);

      // if (total_mass == 0.0)
      //   return P;

      // add inertia forces from element mass

      // if (cMass == 0)  {
      //   // take advantage of lumped mass matrix
      //   double m = total_mass/double(nen);
      //   for (int i=0; i<nen; i++) {
      //     const Vector& accel = theNodes[i]->getTrialAccel();
      //     for (int j=0; j<3; j++) 
      //       P[i*ndf+j] += m * accel(j);
      //   }
      // } else  
      {
        // TODO!!!! update for nen>2, ndf != 6
        // use matrix-vector mult against consistent mass matrix
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

    //
    std::array<Node*, nen> theNodes;

    ID  connectedExternalNodes;    

    Vector p_iner;

    int  parameterID;

private:
    State  e_state;
    int    cMass;
    static MatrixND<ndf*nen,ndf*nen> D; // Damping matrix

};

