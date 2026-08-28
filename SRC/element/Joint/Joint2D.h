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
// Written: Arash & GGD
// Created: 03/02
// Revision: Arash

// Joint2D.h: interface for the Joint2d class.
//
//////////////////////////////////////////////////////////////////////

#ifndef Joint2D_h
#define Joint2D_h

#include <vector>
#include <Matrix.h>
#include <Vector.h>
#include <Element.h>
#include <ID.h>
#include <Domain.h>

class Node;
class UniaxialMaterial;
class Response;
class DamageModel;

class Joint2D : public Element
{
public:
  Joint2D();

  Joint2D(int tag, int nd1, int nd2, int nd3, int nd4, int IntNodeTag,
    UniaxialMaterial* springModels[],
    Domain* theDomain,
    int LrgDisp);

  Joint2D(int tag, int nd1, int nd2, int nd3, int nd4, int IntNodeTag,
    UniaxialMaterial* springModels[],
    Domain* theDomain,
    int LrgDisp,
    DamageModel* damageModels[]);

  ~Joint2D();

  const char* getClassType() const { return "Joint2D"; }

  // methods dealing with domain
  int	getNumExternalNodes() const;
  const	ID& getExternalNodes();
  Node** getNodePtrs();
  int	getNumDOF();

  void	setDomain(Domain* theDomain);
  bool	isSubdomain() { return false; };

  // methods dealing with committed state and update
  int update();
  int commitState();
  int revertToLastCommit();
  int revertToStart();

  // methods to return the current linearized stiffness,
  // damping and mass matrices
  const	Matrix& getTangentStiff();
  const Matrix& getInitialStiff();
  const	Matrix& getDamp();
  const	Matrix& getMass();

  // methods for returning and applying loads
  //virtual Vector &getUVLoadVector(double q1, double q2);
  void	zeroLoad();
  int addLoad(ElementalLoad* theLoad, double loadFactor);
  int addInertiaLoadToUnbalance(const Vector& accel);

  const	Vector& getResistingForce();
  const	Vector& getResistingForceIncInertia();

  // method for obtaining information specific to an element
  Response* setResponse(const char** argv, int argc, OPS_Stream& s);
  int getResponse(int responseID, Information& eleInformation);
  int sendSelf(int commitTag, Channel& theChannel);
  int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);
  void Print(OPS_Stream& s, int flag);

  // AddingSensitivity:BEGIN //////////////////////////////////////////
  int	  addInertiaLoadSensitivityToUnbalance(const Vector& accel, bool tag);
  int   setParameter(const char** argv, int argc, Parameter& param);
  const Vector& getResistingForceSensitivity(int gradNumber);
  const Matrix& getKiSensitivity(int gradNumber);
  const Matrix& getMassSensitivity(int gradNumber);
  int   commitSensitivity(int gradNumber, int numGrads);
  // AddingSensitivity:END ///////////////////////////////////////////

protected:
  int 	addMP_Joint(Domain* theDomain, int RnodeID, int CnodeID, int MainDOF, int FixedEnd, int LrgDispFlag);

private:
  UniaxialMaterial* theSprings[5];
  DamageModel* theDamages[5];
  ID		ExternalNodes;
  ID		InternalConstraints;
  int       fixedEnd[5];
  Node* theNodes[5];
  Domain* TheDomain;
  int		numDof, nodeDbTag, dofDbTag;
  static	Matrix K;
  static	Vector V;

  // AddingSensitivity:BEGIN //////////////////////////////////////////
  int parameterID;
  Vector* theLoadSens;
  // AddingSensitivity:END ///////////////////////////////////////////
};

#endif
