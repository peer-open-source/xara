/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
//
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Purpose: This file contains the class interface for Node.
// A Node provides the abstraction for a node in the FEM.
// Nodes have original position, trial displacement, velocity and 
// acceleration and committed displacement, velocity and acceleration 
// (the last committed() trial quantities).
//
#ifndef Node_h
#define Node_h

#ifdef VIRTUAL
#undef VIRTUAL
#endif

#define VIRTUAL virtual
#include <TaggedObject.h>
#include <MovableObject.h>
// TODO: Remove include of NodeData
#include "NodeData.h"
#include <Vector.h>
#include <Matrix3D.h>
#include <Rotations.h>

class Vector;
class Matrix;
class Channel;
class DOF_Group;
class NodalThermalAction; //L.Jiang [ SIF ]
class Domain;
class Element;
namespace OpenSees {
  struct Versor;
}
using OpenSees::Versor;
using OpenSees::Matrix3D;

class Node :
  public TaggedObject,
  public MovableObject
{
  public:
    enum class Field {
      R1, R2, R3, SE2, SE3, None
    } field;

    // constructors
    Node(int classTag);
    Node(int tag, int classTag);
    Node(int tag, int ndof, double Crd1);
    Node(int tag, int ndof, double Crd1, double Crd2);
    Node(int tag, int ndof, double Crd1, double Crd2, double Crd3, Rotations::Parameters rotationType = Rotations::Parameters::None);
    Node(const Node &theCopy, bool copyMass = true);
    
    // destructor
    VIRTUAL ~Node();

    // public methods dealing with the DOF at the node
    VIRTUAL int  getNumberDOF() const;
    VIRTUAL void setDOF_GroupPtr(DOF_Group *);
    VIRTUAL DOF_Group *getDOF_GroupPtr();

    // public methods for obtaining the nodal coordinates
    VIRTUAL const Vector &getCrds() const;
    VIRTUAL void  setCrds(double Crd1);
    VIRTUAL void  setCrds(double Crd1, double Crd2);
    VIRTUAL void  setCrds(double Crd1, double Crd2, double Crd3);
    VIRTUAL void  setCrds(const Vector &);

    //
    // State
    //
    // public methods dealing with the committed state of the node
    virtual int commitState() noexcept;
    virtual int revertToLastCommit();
    virtual int revertToStart();

    //
    // Response
    //
    // public methods for obtaining committed and trial 
    // response quantities of the node
    VIRTUAL const Vector &getDisp() noexcept;
    VIRTUAL const Vector &getIncrDisp() noexcept;
    VIRTUAL const Vector &getIncrDeltaDisp() noexcept;
    VIRTUAL const Vector &getTrialDisp() noexcept {return *trialDisp;}
    VIRTUAL       Versor  getTrialRotation();
    VIRTUAL       Matrix3D  getRotationTangent(Rotations::Parameters);
    VIRTUAL       Versor  getCommitRotation();

    VIRTUAL const Vector &getVel();
    VIRTUAL const Vector &getAccel();
    VIRTUAL const Vector &getTrialVel();
    VIRTUAL const Vector &getTrialAccel();

    // public methods for updating the trial response quantities
    virtual int setTrialDisp  (double value, int dof) noexcept;
    virtual int setTrialDisp  (const Vector &) noexcept;
    virtual int incrTrialDisp (const Vector &) noexcept;
    VIRTUAL int setTrialVel   (const Vector &);
    VIRTUAL int setTrialAccel (const Vector &);
    VIRTUAL int incrTrialVel  (const Vector &);
    VIRTUAL int incrTrialAccel(const Vector &);

    // Dynamics
    VIRTUAL const Matrix &getMass();
    VIRTUAL const Matrix &getDamp();
    VIRTUAL int setMass(const Matrix &);
    VIRTUAL int addPositionInertia(double value);
    VIRTUAL const Vector &getRV(const Vector &V);
    VIRTUAL int setRayleighDampingFactor(double alphaM);

    // Eigen vectors
    VIRTUAL int setNumEigenvectors(int numVectorsToStore);
    VIRTUAL int setEigenvector(int mode, const Vector &eigenVector);
    VIRTUAL const Matrix &getEigenvectors(void);

    // Thermodynamics
    VIRTUAL NodalThermalAction* getNodalThermalActionPtr();
    VIRTUAL void setNodalThermalActionPtr(NodalThermalAction*);
    
    //
    // Load information
    //
    VIRTUAL void zeroUnbalancedLoad() noexcept;
    VIRTUAL int addUnbalancedLoad(const Vector &load, double fact = 1.0);
    // VIRTUAL int addInertiaLoadToUnbalance(const Vector &accel, double fact = 1.0);
    VIRTUAL int addResidualInertia(int dof, double accel);
    VIRTUAL int addResidual(const Vector &, double scale);
    VIRTUAL const Vector &getUnbalancedLoad() noexcept;
    VIRTUAL const Vector &getUnbalancedLoadIncInertia();


    // Misc. Responses
    VIRTUAL const Vector &getReaction();
    VIRTUAL int   addReactionForce(const Vector &, double factor);
    VIRTUAL int   resetReactionForce(int flag);

    VIRTUAL const Vector *getResponse(NodeData);
    int fillResponse(NodeData responseType, Vector& result, int offset=0);
    
    //
    // Parallel
    //
    VIRTUAL int sendSelf(int commitTag, Channel &);
    VIRTUAL int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &);

    //
    // Misc
    //
    VIRTUAL void Print(OPS_Stream &s, int flag=0);


    // Sensitivity
    int addInertiaLoadSensitivityToUnbalance(const Vector &accel, 
                                             double fact = 1.0, 
                                             bool tag=false);
    Matrix getMassSensitivity();
    VIRTUAL const Matrix &getDampSensitivity();
    int    getCrdsSensitivity();
    int    saveDispSensitivity(const Vector &v, int gradNum, int numGrads);
    int    saveVelSensitivity(const Vector &vdot, int gradNum, int numGrads);
    int    saveAccelSensitivity(const Vector &vdot, int gradNum, int numGrads);
    double getDispSensitivity(int dof, int gradNum);
    double getVelSensitivity(int dof, int gradNum);
    double getAccSensitivity(int dof, int gradNum);
    int    setParameter(const char **argv, int argc, Parameter &);
    int    updateParameter(int parameterID, Information &);
    int    activateParameter(int parameterID);

#if 1
    //
    // Display
    //
    VIRTUAL int getDisplayCrds(Vector &results, double fact, int displayMode=0);
#endif

    // Reset global matrices
    static int resetGlobalMatrices();


    Domain *getDomain() {return theDomain;}
    void setDomain(Domain *model) {theDomain = model;}

  protected:
    double *vel, *accel;              // double arrays holding the vel and accel values
    Vector *trialDisp, *commitDisp, *incrDisp, *incrDeltaDisp;
    Vector *commitVel, *commitAccel;  // committed quantities
    Vector *trialVel,  *trialAccel;   // trial quantities
    Vector *unbalLoad;                // unbalanced load
    Versor *rotation;
    Rotations::Parameters rotationType;

  private:
    double *disp;

    // State
    Vector *unbalLoadWithInertia;
    Matrix *theEigenvectors;

#if 1
    Domain* theDomain;
#endif

    // Private global state
    int setGlobalMatrices();

    static Matrix **theMatrices;
    static int numMatrices;
    int index;


    // priavte methods used to create the Vector objects 
    // for the committed and trial response quantaties.
    virtual int createDisp();
    int createVel();
    int createAccel();

    // private data associated with each node object
    int numberDOF;                    // number of dof at Node
    DOF_Group *theDOF_GroupPtr;       // pointer to associated DOF_Group
    NodalThermalAction *theNodalThermalActionPtr; //Added by Liming Jiang for pointer to nodalThermalAction, [SIF]

    Vector xyz;                       // coordinates in ref. configuration
    double coord_data[3];
    

    Matrix *R;                          // nodal participation matrix
    Matrix *mass;                       // pointer to mass matrix
    double alphaM;                      // rayleigh damping factor 
    double position_inertia;
    enum class MassType {
      Full,
      Classical, 
      None
    } mass_type = MassType::None;



    int dbTag1, dbTag2, dbTag3, dbTag4; // needed for database

    // AddingSensitivity:BEGIN /////////////////////////////////////////
    Matrix *dispSensitivity;
    Matrix *velSensitivity;
    Matrix *accSensitivity;
    int parameterID;

    Vector *reaction;
};

#undef VIRTUAL
#endif

