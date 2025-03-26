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
// Description: This file contains the class definition for 
// CrdTransf.h. CrdTransf provides the abstraction of a frame 
// coordinate transformation. It is an abstract base class.
//
// Written: Remo Magalhaes de Souza (rmsouza@ce.berkeley.edu)
// Created: 04/2000
// Revision: A
//
#ifndef CrdTransf_h
#define CrdTransf_h

#include <MovableObject.h>
#include <TaggedObject.h>

class Vector;
class ID;
class Matrix;
class Node;
class Response;

class CrdTransf: public TaggedObject, public MovableObject
{
public:
    CrdTransf(int tag, int classTag);
    CrdTransf();
    virtual ~CrdTransf();

    virtual CrdTransf *getCopy2d() {return nullptr;};
    virtual CrdTransf *getCopy3d() {return nullptr;};

    virtual int getLocalAxes(Vector &xAxis, Vector &yAxis, Vector &zAxis);
    virtual int getRigidOffsets(Vector &offsets);
  
    virtual int    initialize(Node *node1Pointer, Node *node2Pointer) = 0;
    virtual int    update() = 0;
    virtual double getInitialLength() = 0;
    virtual double getDeformedLength() = 0;
    
    virtual int commitState() = 0;
    virtual int revertToLastCommit() = 0;
    virtual int revertToStart() = 0;

    virtual const Vector &getBasicTrialDisp() = 0;
    virtual const Vector &getBasicIncrDisp() = 0;
    virtual const Vector &getBasicIncrDeltaDisp() = 0;
    virtual const Vector &getBasicTrialVel() = 0;
    // virtual const Vector &getBasicTrialAccel() = 0;

    virtual const Vector &getGlobalResistingForce(const Vector &basicForce, const Vector &uniformLoad) = 0;
    virtual const Matrix &getGlobalStiffMatrix(const Matrix &basicStiff, const Vector &basicForce) = 0;
    virtual const Matrix &getInitialGlobalStiffMatrix(const Matrix &basicStiff) = 0;

    // method used to rotate consistent mass matrix
    virtual const Matrix &getGlobalMatrixFromLocal(const Matrix &local) = 0;


    // method for obtaining information specific to a coordinate transformation
    virtual Response *setResponse(const char **argv, int argc, 
                                  OPS_Stream &theHandler);
    virtual int getResponse(int responseID, Information &eleInformation);

    // methods used in post-processing only
    virtual const Vector &getPointGlobalCoordFromLocal(const Vector &localCoords) = 0;
    virtual const Vector &getPointGlobalDisplFromBasic(double xi, const Vector &basicDisps) = 0;
    virtual const Vector &getPointLocalDisplFromBasic(double xi, const Vector &basicDisps) = 0;

    // Sensitivity
    virtual const Vector &getBasicDisplTotalGrad(int grad);
    virtual const Vector &getBasicDisplFixedGrad();
    virtual const Vector &getGlobalResistingForceShapeSensitivity(const Vector &pb, const Vector &p0, int gradNumber);
    virtual bool   isShapeSensitivity() {return false;}
    virtual double getLengthGrad() {return 0.0;}
    virtual double getd1overLdh() {return 0.0;}
    //
protected:
    
private:
};

// some additional methods related to prototypes created for copy constructors
#if !defined(OPS_USE_RUNTIME)
extern bool       OPS_addCrdTransf(CrdTransf *newComponent);
extern CrdTransf *OPS_getCrdTransf(int tag);
#endif
extern bool       OPS_removeCrdTransf(int tag);
extern void       OPS_clearAllCrdTransf();
extern void       OPS_printCrdTransf(OPS_Stream &s, int flag=0);
extern ID       OPS_getAllCrdTransfTags();

#endif
