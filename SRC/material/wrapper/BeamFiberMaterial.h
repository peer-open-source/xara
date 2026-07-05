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
// Description: This file contains the class definition of BeamFiberMaterial.
// The BeamFiberMaterial class is a wrapper class that performs static
// condensation on a three-dimensional material model to give the 11, 12, and 13
// stress components which can then be integrated over an area to model a
// shear flexible 3D beam.
//
// Written: MHS
// Created: Aug 2001
//

#include <Vector.h>
#include <Matrix.h>
#include <VectorND.h>
#include <MatrixND.h>
#include <Matrix3D.h>
#include <ID.h> 
#include <NDMaterial.h>


namespace OpenSees {

class BeamFiberMaterial: public NDMaterial {

  public:
    BeamFiberMaterial(int tag, NDMaterial &);
    BeamFiberMaterial();
    virtual ~BeamFiberMaterial();
    const char* getClassType() const override {return "BeamFiberMaterial";}

    int setTrialStrain( const Vector &);
    const Vector& getStrain();
    const Vector& getStress();
    const Matrix& getTangent();
    const Matrix& getInitialTangent();
    bool threadSafe() const override {
      if (theMaterial != nullptr)
        return theMaterial->threadSafe();
      else
        return false;
    }

    double getRho();

    int commitState();
    int revertToLastCommit();
    int revertToStart();

    NDMaterial *getCopy();
    NDMaterial *getCopy(const char *type);
    const char *getType() const;
    int getOrder() const;

    Response *setResponse(const char **argv, int argc, OPS_Stream &) override;
    int getResponse(int responseID, Information &) override;

    void Print(OPS_Stream &s, int flag);

    int sendSelf(int commitTag, Channel &);
    int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &);

    int setParameter(const char **argv, int argc, Parameter &param);

    const Vector& getStressSensitivity(int gradIndex, bool conditional);

  private:
    double Tstrain22;
    double Tstrain33;
    double Tgamma23;
    double Cstrain22;
    double Cstrain33;
    double Cgamma23;
    // struct {
    //   double strain22, strain33, gamma23;
    // } past, pres;

    NDMaterial *theMaterial;

    VectorND<3>   strain;
    VectorND<3>   stress;
    MatrixND<3,3> tangent;

    Vector rvec;
    Matrix rmat;
};
} // namespace OpenSees
