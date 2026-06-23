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
                                                                        
// $Revision: 1.4 $
// $Date: 2003-02-14 23:01:05 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/beam2d/beam2d02.h,v $
                                                                        
                                                                        
// File: ~/model/element/beam2d02.h
//
// Written: fmk 11/95
// Revised:
//
// Purpose: This file contains the class definition for beam2d02.
// beam2d02 is a linear elastic plane frame member.

#ifndef beam2d02_h
#define beam2d02_h

#include <Element.h>
#include <Node.h>
#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
class CrdTransf2d;

class beam2d02 : public Element
{
  public:
    beam2d02();    
    beam2d02(int tag, double A, double E, double I, int Nd1, int Nd2,
	     CrdTransf2d &theTrans, double rho = 0.0);
    ~beam2d02();

    int getNumExternalNodes() const;
    const ID &getExternalNodes();
    Node **getNodePtrs();

    int getNumDOF();
    void setDomain(Domain *);
    
    int commitState();
    int revertToLastCommit(void);        
    int revertToStart(void);            
    
    const Matrix &getTangentStiff();
    const Matrix &getInitialStiff();
    const Matrix &getMass();    

    void zeroLoad();	
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    const Vector &getResistingForce();
    const Vector &getResistingForceIncInertia();            
    
    int sendSelf(int commitTag, Channel &);
    int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &);

    void Print(OPS_Stream &s, int flag);
    
  private:
    const Matrix &getStiff();    
    
    // private attributes - a copy for each object of the class        
    double A,E,I,M;
    double L,sn,cs;
    ID  connectedExternalNodes;    
    Node *theNodes[2];

    Matrix Kd; // the stiffness matrix
    Matrix m; // the mass matrix	
    Matrix d; // the damping matrix

    Vector q;
    Vector rForce;
    Vector load;
    Matrix trans;  // for the transformation matrix
//    int isStiffFormed;
    CrdTransf2d *theCoordTrans;
};

#endif

