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
// $Date: 2003-02-14 23:01:24 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/BeamFiberMaterial2dPS.h,v $

// Written: MHS
// Created: Aug 2001
//
// Description: This file contains the class definition of BeamFiberMaterial2dPS.
// The BeamFiberMaterial2dPS class is a wrapper class that performs static
// condensation on a three-dimensional material model to give the 11 and 12
// stress components which can then be integrated over an area to model a
// shear flexible 2D beam.
//
#include <stdlib.h> 
#include <math.h> 

#include <Vector.h>
#include <Matrix.h>
#include <ID.h> 
#include <NDMaterial.h>

class BeamFiberMaterial2dPS: public NDMaterial {

  public:
    BeamFiberMaterial2dPS(int tag, NDMaterial &theMat);
    BeamFiberMaterial2dPS();
    virtual ~BeamFiberMaterial2dPS();

    int setTrialStrain( const Vector &strainFromElement);
    const Vector& getStrain();
    const Vector& getStress();
    const Matrix& getTangent();
    const Matrix& getInitialTangent();

    double getRho();

    int commitState();
    int revertToLastCommit();
    int revertToStart();

    NDMaterial *getCopy();
    NDMaterial *getCopy(const char *type);
    const char *getType() const;
    int getOrder() const; 

    void Print(OPS_Stream &s, int flag);

    int sendSelf(int commitTag, Channel &);
    int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &);

    int setParameter(const char **argv, int argc, Parameter &param);

    const Vector& getStressSensitivity(int gradIndex, bool conditional);
    int commitSensitivity(const Vector &depsdh, int gradIndex, int numGrads);

  private:
    double Tstrain22;
    double Cstrain22;

    NDMaterial *theMaterial;

    Vector strain;

    static Vector stress;
    static Matrix tangent;
};





