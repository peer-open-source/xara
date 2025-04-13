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

// $Revision: 1.14 $
// $Date: 2008-08-26 16:47:42 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/FiberSection3dThermal.h,v $

// Written: fmk
// Created: 04/01
//
// Description: This file contains the class definition for
// FiberSection3dThermal.h. FiberSection3dThermal provides the abstraction of a
// 3d beam section discretized by fibers. The section stiffness and
// stress resultants are obtained by summing fiber contributions.
// Modified for SIF modelling by Jian Jiang,Liming Jiang [http://openseesforfire.github.io]
// Corrected by Giovanni Rinaldin, 2024

#ifndef FiberSection3dThermal_h
#define FiberSection3dThermal_h

#include <vector>
#include <FrameSection.h>
#include <Vector.h>
#include <VectorND.h>
#include <Matrix.h>

class UniaxialMaterial;
// class Fiber;
class Response;
class SectionIntegration;

namespace OpenSees {

class FiberSection3dThermal : public FrameSection
{
  public:
    FiberSection3dThermal();
//  FiberSection3dThermal(int tag, int numFibers, Fiber **fibers, 
//  				UniaxialMaterial &torsion, bool compCentroid=true);
    FiberSection3dThermal(int tag, int numFibers, UniaxialMaterial &torsion, bool compCentroid=true);
    ~FiberSection3dThermal();

    const char *getClassType() const {return "FiberSection3dThermal";};

    int   setTrialSectionDeformation(const Vector &deforms);
    const Vector &getSectionDeformation();

	const Vector &getTemperatureStress(const Vector& dataMixed); //JJadd to get Ft=EA*Elongation//

    const Vector &getStressResultant();
    const Matrix &getSectionTangent();
    const Matrix &getInitialTangent();

    int   commitState();
    int   revertToLastCommit();
    int   revertToStart();
 
    FrameSection *getFrameCopy();
    const ID &getType ();
    int getOrder () const;
    
    int sendSelf(int cTag, Channel &theChannel);
    int recvSelf(int cTag, Channel &theChannel,
		 FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag = 0);
	    
    Response *setResponse(const char **argv, int argc, 
			  OPS_Stream &s);
    int getResponse(int responseID, Information &info);

//  int addFiber(Fiber &theFiber);
    int addFiber(UniaxialMaterial &theMat, const double area, const double y, const double z);

    // AddingSensitivity:BEGIN //////////////////////////////////////////
    int setParameter(const char **argv, int argc, Parameter &param);

    const Vector & getStressResultantSensitivity(int gradIndex, bool conditional);
    const Matrix & getSectionTangentSensitivity(int gradIndex);
    int   commitSensitivity(const Vector& sectionDeformationGradient, int gradIndex, int numGrads);

    const Vector & getSectionDeformationSensitivity(int gradIndex);
    // AddingSensitivity:END ///////////////////////////////////////////

	double determineFiberTemperature(const Vector& , double , double);
    const Vector& getThermalElong();
  protected:

  private:
    int numFibers, sizeFibers;                   // number of fibers in the section
    UniaxialMaterial **theMaterials; // array of pointers to materials
    struct FiberData {
      double area,
             y,
             z,
             temp,
             temp_max;
      UniaxialMaterial *material;
    };
    std::vector<FiberData> fibers;
//  double   *matData;               // data for the materials [yloc, zloc, area]
    double   kData[16];               // data for ks matrix
    double   sData[4];               // data for s vector

    double QzBar, QyBar, ABar;
    double yBar;       // Section centroid
    double zBar;
    bool computeCentroid;

    static ID code;

    Vector e;          // trial section deformations
    Vector eCommit;    // committed section deformations
    Vector *s;         // section resisting forces  (axial force, bending moment)
    Matrix *ks;        // section stiffness

    UniaxialMaterial *theTorsion;
    Vector sT;  // JZ  section resisting forces, caused by the temperature
    VectorND<4> AverageThermalElong;
    //double  *TemperatureTangent; // JZ  the E of E*A*alpha*DeltaT

    // AddingSensitivity:BEGIN //////////////////////////////////////////
    int parameterID;
    Matrix *SHVs;
    // AddingSensitivity:END ///////////////////////////////////////////
};

} // namespace
#endif
