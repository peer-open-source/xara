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
**                                                                    **
** ****************************************************************** */

// $Revision: 1.23 $
// $Date: 2010-09-13 21:29:28 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/NDMaterial.h,v $


#ifndef NDMaterial_h
#define NDMaterial_h
//
// Description: This file contains the class definition for NDMaterial.
// NDMaterial is an abstract base class and thus no objects of it's type
// can be instantiated. It has pure virtual functions which must be
// implemented in it's derived classes.
//
//
// Written: MHS
// Created: Feb 2000
// Revision: A
//
#include <cassert>
#include <Domain.h>
#include <TaggedObject.h>
#include <MovableObject.h>
#include <MaterialBuilder.h>
class Matrix;
class ID;
class Vector;
class Information;
class Response;
class Domain;

using namespace OpenSees;

class NDMaterial : public MaterialBuilder, public MovableObject
{
  public:
    NDMaterial(int tag, int classTag);
    NDMaterial();
    virtual ~NDMaterial();

    int setDomain(Domain *domain) {m_domain = domain; return 0;}

    virtual NDMaterial *getCopy() = 0;
    virtual NDMaterial *getCopy(const char *type) override;

    // methods to set state and retrieve state
    virtual double getRho() { return 0.0; }

    virtual int setTrialStrain(const Vector &) =0;
    virtual int setTrialStrainIncr(const Vector &);
    virtual const Matrix &getTangent();
    virtual const Matrix &getInitialTangent() {return this->getTangent();}

    //Added by L.Jiang, [SIF]
    virtual double getThermalTangentAndElongation(double &TempT, double &, double &);
    virtual double setThermalTangentAndElongation(double &TempT, double &, double &);
    virtual const Vector& getTempAndElong();
    //Added by L.Jiang, [SIF]

    virtual const Vector &getStress();
    virtual const Vector &getStrain();
    virtual bool threadSafe() const {return false;}

    // State management
    // Return 0 on success.
    virtual int commitState() = 0;
    virtual int revertToLastCommit() = 0;
    virtual int revertToStart() = 0;


    virtual const char *getType() const = 0;
    virtual int getOrder() const {return 0;};  //??

    virtual Response *setResponse(const char **argv, int argc, OPS_Stream &);
    virtual int getResponse(int responseID, Information &);

    void Print(OPS_Stream &s, int flag) override;

    // Sensitivity
    virtual const Vector & getStressSensitivity         (int gradIndex, bool conditional);
    virtual const Vector & getStrainSensitivity         (int gradIndex);
    virtual const Matrix & getTangentSensitivity        (int gradIndex);
    virtual const Matrix & getInitialTangentSensitivity (int gradIndex);
    virtual const Matrix & getDampTangentSensitivity    (int gradIndex);
    virtual double         getRhoSensitivity            (int gradIndex);
    virtual int            commitSensitivity            (const Vector & strainGradient, int gradIndex, int numGrads);

  

  protected:
    int getCharacteristicLength(double& l) const { l = c_length; return c_length_set ? 0 : -1;}
    int getTimeStep(double& dt) const { 
      if (m_domain) {
        dt = m_domain->getDT();
        return 0;
      } else {
        return -1;
      }
    }

  private:
    static Matrix errMatrix;
    static Vector errVector;
    double c_length = 0.0;     // characteristic length for nonlocal damage models
    bool c_length_set = false; // flag to indicate if characteristic length has been set

    Domain *m_domain = nullptr; // pointer to the domain in which the material is used
};

#endif
