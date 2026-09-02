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

// $Revision: 1.5 $
// $Date: 2010-04-23 22:53:56 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/joint/Joint3D.h,v $

// Written: Arash Altoontash, Gregory Deierlein
// Created: 04/03
// Revision: Arash

// Joint3D.h: interface for the Joint3D class.
//
//////////////////////////////////////////////////////////////////////

#ifndef Joint3D_h
#define Joint3D_h

#include <Matrix.h>
#include <Vector.h>
#include <Element.h>
#include <ID.h>
#include <Domain.h>

class Node;
class UniaxialMaterial;
class Response;

class Joint3D : public Element
{
public:
  Joint3D();

  Joint3D(int tag, int nd1, int nd2, int nd3, int nd4, int nd5, int nd6, int IntNodeTag,
    UniaxialMaterial* springModels[],
    Domain* theDomain,
    int LrgDisp);


  ~Joint3D();

  const char* getClassType() const { return "Joint3D"; }

  // methods dealing with domain
  int	getNumExternalNodes() const;
  const	ID& getExternalNodes();
  Node** getNodePtrs();
  int	getNumDOF();

  void	setDomain(Domain* theDomain);
  bool	isSubdomain() { return false; };

  // methods dealing with committed state and update
  int update() override;
  int commitState();
  int revertToLastCommit();
  int revertToStart();

  // methods to return the current linearized stiffness,
  // damping and mass matrices
  const	Matrix& getTangentStiff() override;
  const Matrix& getInitialStiff() override;
  const	Matrix& getDamp() override;
  const	Matrix& getMass() override;

  // methods for returning and applying loads
  //virtual Vector &getUVLoadVector(double q1, double q2);
  void	zeroLoad();
  int addLoad(ElementalLoad* theLoad, double loadFactor);
  int addInertiaLoadToUnbalance(const Vector& accel);

  const	Vector& getResistingForce(void);
  const	Vector& getResistingForceIncInertia(void);

  // method for obtaining information specific to an element
  Response* setResponse(const char** argv, int argc, OPS_Stream& s);
  int getResponse(int responseID, Information& eleInformation);
  int sendSelf(int commitTag, Channel& theChannel);
  int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);
  void Print(OPS_Stream& s, int flag);


protected:
  int addMP_Joint(Domain* theDomain, int RetNodeID, int ConNodeID,
    int RotNodeID, int Rdof, int DspNodeID, int Ddof,
    int LrgDispFlag);

private:
  UniaxialMaterial* theSprings[3];
  ID		ExternalNodes;
  ID		InternalConstraints;
  Node* theNodes[7];
  Domain* TheDomain;
  int		numDof, nodeDbTag, dofDbTag;
  static	Matrix K;
  static	Vector V;
};

#endif
