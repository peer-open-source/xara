#pragma once
#include <array>

#include <Domain.h>
#include <Element.h>
#include <VectorND.h>
#include <MatrixND.h>
#include <ID.h>

class Node;
class Domain;
class Response;
class Rotation;

using namespace OpenSees;

#include <State.h>

template <int nen, int ndm, int ndf, int mass_flag=0>
class FiniteElement : public Element {
public:
     FiniteElement(int tag, int classtag)
       : Element(tag, classtag),
         connectedExternalNodes(nen),
         e_state(State::None),
         p_iner(nen*ndf),
         parameterID(0)
    {
         for (int i=0; i<nen; i++)
           theNodes[i] = nullptr;
     }

     FiniteElement(int tag, int classtag, std::array<int, nen>& nodes)
       : Element(tag, classtag),
         connectedExternalNodes(nen),
         e_state(State::None),
         p_iner(nen*ndf),
         parameterID(0)
    {
         for (int i=0; i<nen; i++) {
           connectedExternalNodes(i) = nodes[i];
           theNodes[i] = nullptr;
         }
     }


    // For Element
    virtual const ID& getExternalNodes() final {
      return connectedExternalNodes;
    }
    virtual Node **getNodePtrs() final {return theNodes.data();}
    virtual int  getNumExternalNodes() const final {return nen;}
    virtual int  getNumDOF() final {return nen*ndf;}
    virtual void zeroLoad() {
    // TODO: need to reconcile with BasicFrame3d::zeroLoad()
      p_iner.Zero();
    }

    void
    setDomain(Domain *theDomain) final
    {
        if (theDomain == nullptr) {
         for (int i=0; i<nen; i++) {
           theNodes[i] = nullptr;
         }
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

        this->DomainComponent::setDomain(theDomain);

        if (this->setState(State::Init) != 0)
          return;

//      if (this->setState(State::Pres) != 0)
//        return;
    }

        
    int
    addInertiaLoadToUnbalance(const Vector &accel)
    {
      if (total_mass == 0.0)
        return 0;

      // add ( - fact * M R * accel ) to unbalance
      if (cMass == 0) {
        // take advantage of lumped mass matrix
        double m = 0.5*total_mass;
        for (int i=0; i<nen; i++) {
          const Vector& Raccel = theNodes[i]->getRV(accel);
          for (int j=0; i<3; i++) {
            p_iner[i*ndf+j] -= m * Raccel(j);
          }
        }
      }

      else  {
        // TODO: Move this to FiniteElement::getAcceleration() ?

        constexpr static int nrv = 6;

        // use matrix vector multip. for consistent mass matrix
        static VectorND<nen*ndf> Raccel;
        for (int i=0; i<nen; i++) {
          // get R * accel from the nodes
          const Vector& rv = theNodes[i]->getRV(accel);
          if (nrv != rv.Size()) {
            opserr << "addInertiaLoadToUnbalance matrix and vector sizes are incompatible\n";
            return -1;
          }
          for (int j=0; i<nrv; i++)  {
            Raccel[i*ndf+j] = rv[j];
          }
        }
        p_iner.addMatrixVector(1.0, this->getMass(), Raccel, -1.0);
      }

      return 0;
    }

    const Vector &getResistingForceIncInertia()
    {
      // TODO!!!! update for nen>2, ndf != 6
      static VectorND<nen*ndf> P_{0.0};
      static Vector P(P_);
      P = this->getResistingForce(); 
      
      // add the damping forces if rayleigh damping
      if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
        P.addVector(1.0, this->getRayleighDampingForces(), 1.0);


      if (total_mass == 0.0)
        return P;

      // add inertia forces from element mass
      const Vector &accel1 = theNodes[0]->getTrialAccel();
      const Vector &accel2 = theNodes[1]->getTrialAccel();    
        
      if (cMass == 0)  {
        // take advantage of lumped mass matrix
        double m = 0.5*total_mass;

        P(0) += m * accel1(0);
        P(1) += m * accel1(1);
        P(2) += m * accel1(2);

        P(6) += m * accel2(0);
        P(7) += m * accel2(1);
        P(8) += m * accel2(2);

      } else  {
        // use matrix vector multip. for consistent mass matrix
        static Vector accel(12);
        for (int i=0; i<6; i++)  {
          accel(i)   = accel1(i);
          accel(i+6) = accel2(i);
        }
        P.addMatrixVector(1.0, this->getMass(), accel, 1.0);
      }
      
      return P;
    }
    int setParameter(const char **argv, int argc, Parameter &param)
    {
      if (argc < 1)
        return -1;

      // don't do anything if MaterialStageParameter calls this element
      if (strcmp(argv[0],"updateMaterialStage") == 0) {
        return -1;
      }

      return -1;
    }

    int updateParameter(int paramID, Information &info)
    {
      return -1;
    }

    int activateParameter(int passedParameterID)
    {
      parameterID = passedParameterID; 
      return 0;
    }

#if 0
    const Matrix &getMassSensitivity(int gradNumber)
    {
      // From DispBeamColumn
      static MatrixND<12,12> M_;
      static Matrix M(M_);
      M.Zero();

      if (total_mass == 0.0 || parameterID != 1)
        return M;
      
      double L = theCoordTransf->getInitialLength();
      if (cMass == 0)  {
        // lumped mass matrix
        //double m = 0.5*rho*L;
        double m = 0.5*L;
        M(0,0) = M(1,1) = M(2,2) = M(6,6) = M(7,7) = M(8,8) = m;
      }

      else {
        // consistent mass matrix
        static Matrix ml(12,12);
        //double m = rho*L/420.0;
        double m = L/420.0;
        ml(0,0) = ml(6,6) = m*140.0;
        ml(0,6) = ml(6,0) = m*70.0;
        //ml(3,3) = ml(9,9) = m*(Jx/A)*140.0;  // TODO: CURRENTLY NO TORSIONAL MASS 
        //ml(3,9) = ml(9,3) = m*(Jx/A)*70.0;   // TODO: CURRENTLY NO TORSIONAL MASS
        
        ml(2, 2) = ml( 8, 8) =  m*156.0;
        ml(2, 8) = ml( 8, 2) =  m*54.0;
        ml(4, 4) = ml(10,10) =  m*4.0*L*L;
        ml(4,10) = ml(10, 4) = -m*3.0*L*L;
        ml(2, 4) = ml( 4, 2) = -m*22.0*L;
        ml(8,10) = ml(10, 8) = -ml(2,4);
        ml(2,10) = ml(10, 2) =  m*13.0*L;
        ml(4, 8) = ml( 8, 4) = -ml(2,10);
        
        ml(1, 1) = ml(7,7) = m*156.0;
        ml(1, 7) = ml(7,1) = m*54.0;
        ml(5, 5) = ml(11,11) = m*4.0*L*L;
        ml(5,11) = ml(11,5) = -m*3.0*L*L;
        ml(1, 5) = ml(5,1) = m*22.0*L;
        ml(7,11) = ml(11,7) = -ml(1,5);
        ml(1,11) = ml(11,1) = -m*13.0*L;
        ml(5, 7) = ml(7,5) = -ml(1,11);
        
        // transform local mass matrix to global system
        M = theCoordTransf->getGlobalMatrixFromLocal(ml);
      }
      
      return M;
    }
#endif
protected:

#ifdef FEFT
    // Implemented by children
    virtual MatrixND<ndf*nen,ndf*nen> getTangent(State state, int rate) =0;
    virtual VectorND<ndf*nen>         getForce(State state, int rate) =0;

    // Supplied to children
    const VectorND<ndf>& getNodeUnknowns(int tag, int rate);
    const VectorND<ndm>& getNodePosition(int tag, State state);
    const Rotation&      getNodeRotation(int tag, State state);
    const VectorND<ndm>& getNodeVelocity(int tag);
    const VectorND<ndm>& getNodeLocation(int tag, State state);
#endif

// TODO: Rename setNodes to setReference
    virtual int setNodes() = 0;

    // Supplied for children
    inline int setState(State state) {

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

   int parameterID;

private:
    State  e_state;
    int cMass;
    double rho;

    double total_mass,
            twist_mass;

};

