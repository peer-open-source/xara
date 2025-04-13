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
// Created: 04/01
//
// Description: This file contains the class definition for 
// FiberSection3d.h. FiberSection3d provides the abstraction of a 
// 3d beam section discretized by fibers. The section stiffness and
// stress resultants are obtained by summing fiber contributions.

#ifndef FiberSection3d_h
#define FiberSection3d_h

#include <FrameSection.h>
#include <Vector.h>
#include <Matrix.h>
#include <VectorND.h>
#include <memory>

class Response;
class UniaxialMaterial;

class FiberSection3d : public FrameSection
{
  public:
    FiberSection3d(); 
    FiberSection3d(int tag, int numFibers, UniaxialMaterial &torsion, bool compCentroid=true);
#if 0
    FiberSection3d(int tag, int numFibers, Fiber **fibers, 
		   UniaxialMaterial &torsion, bool compCentroid=true);
    FiberSection3d(int tag, int numFibers, UniaxialMaterial **mats,
		   SectionIntegration &si, UniaxialMaterial &torsion, bool compCentroid=true);
#endif
    ~FiberSection3d();

    const char *getClassType() const {return "FiberSection3d";};

    int   setTrialSectionDeformation(const Vector &deforms);
    const Vector &getSectionDeformation();

    int   getIntegral(Field field, State state, double& value) const;
    const Vector &getStressResultant();
    const Matrix &getSectionTangent();
    const Matrix &getInitialTangent();

    int   commitState();
    int   revertToLastCommit();    
    int   revertToStart();
 
    FrameSection *getFrameCopy();
    const ID &getType();
    int getOrder () const;
 
    int sendSelf(int cTag, Channel &theChannel);
    int recvSelf(int cTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag = 0);
	    
    Response *setResponse(const char **argv, int argc, 
			  OPS_Stream &s);
    int getResponse(int responseID, Information &info);

    int addFiber(UniaxialMaterial &theMat, const double area, const double y, const double z);

    // AddingSensitivity:BEGIN //////////////////////////////////////////
    int setParameter(const char **argv, int argc, Parameter &param);

    const Vector & getStressResultantSensitivity(int gradIndex, bool conditional);
    const Matrix & getSectionTangentSensitivity(int gradIndex);
    int   commitSensitivity(const Vector& sectionDeformationGradient, int gradIndex, int numGrads);

    const Vector & getSectionDeformationSensitivity(int gradIndex);
    // AddingSensitivity:END ///////////////////////////////////////////

    //by SAJalali
    double getEnergy() const;


  protected:

  private:
    int numFibers, sizeFibers;         // number of fibers in the section
    UniaxialMaterial **theMaterials;   // array of pointers to materials
    std::shared_ptr<double[]> matData; // data for the materials [yloc, zloc, and area]
    double   kData[16];                // data for ks matrix 

    double QzBar, QyBar, Abar;
    double yBar;                       // Section centroid
    double zBar;
    bool computeCentroid;

    static ID code;

    Vector  e;         // trial section deformations 
    Vector  s;         // section resisting forces  (axial force, bending moment)
    Matrix  ks;        // section stiffness

    OpenSees::VectorND<4> eData, sData;
    UniaxialMaterial *theTorsion;
    void *pool;        // thread pool
};

#endif
