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
// Written: MHS
// Created: Feb 2001
//
// Description: This file contains the class definition for TimoshenkoBeamColumn2d.
// The element displacement field gives rise to constant axial strain, constant shearing and
// linear curvature.

#ifndef TimoshenkoBeamColumn2d_h
#define TimoshenkoBeamColumn2d_h

#ifndef _bool_h
#include <stdbool.h>
#endif

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <BeamIntegration.h>

class Node;
class SectionForceDeformation;
class CrdTransf;
class Response;

class TimoshenkoBeamColumn2d : public Element
{
  public:
    TimoshenkoBeamColumn2d(int tag, int nd1, int nd2,
		     int numSections, SectionForceDeformation **s,
		     BeamIntegration &bi, CrdTransf &coordTransf,
		     double rho = 0.0);
    TimoshenkoBeamColumn2d();
    ~TimoshenkoBeamColumn2d();

    const char *getClassType() const {return "TimoshenkoBeamColumn2d";};
    static constexpr const char* class_name = "TimoshenkoBeamColumn2d";

    int getNumExternalNodes() const;
    const ID &getExternalNodes();
    Node **getNodePtrs();

    int getNumDOF();
    void setDomain(Domain *theDomain);

    // public methods to set the state of the element    
    int commitState();
    int revertToLastCommit();
    int revertToStart();

    // public methods to obtain stiffness, mass, damping and residual information    
    int update();
    const Matrix &getTangentStiff();
    const Matrix &getInitialStiff();
    const Matrix &getMass();

    void zeroLoad();
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    const Vector &getResistingForce();
    const Vector &getResistingForceIncInertia();            

    // public methods for element output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker 
		  &theBroker);
    void Print(OPS_Stream &s, int flag =0);

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInfo);

    // AddingSensitivity:BEGIN //////////////////////////////////////////
    int setParameter(const char **argv, int argc, Parameter &param);
    int            updateParameter(int parameterID, Information &info);
    int            activateParameter(int parameterID);
    const Vector & getResistingForceSensitivity(int gradNumber);
    const Matrix & getKiSensitivity(int gradNumber);
    const Matrix & getMassSensitivity(int gradNumber);
    int            commitSensitivity(int gradNumber, int numGrads);
    // AddingSensitivity:END ///////////////////////////////////////////

  protected:
    
  private:
    const Matrix &getInitialBasicStiff();

    int numSections;
    SectionForceDeformation **theSections; // pointer to the ND material objects
    CrdTransf *crdTransf;        // pointer to coordinate tranformation object 

    BeamIntegration *beamInt;

    ID connectedExternalNodes; // Tags of quad nodes

    Node *theNodes[2];

    static Matrix K;		// Element stiffness, damping, and mass Matrix
    static Vector P;		// Element resisting force vector

    Vector Q;		        // Applied nodal loads
    Vector q;		        // Basic force
    double q0[3];       // Fixed end forces in basic system
    double p0[3];       // Reactions in basic system

    double rho;			// Mass density per unit length

    enum {maxNumSections = 20};

    static double workArea[];

    // AddingSensitivity:BEGIN //////////////////////////////////////////
    int parameterID;
    // AddingSensitivity:END ///////////////////////////////////////////
};

#endif

