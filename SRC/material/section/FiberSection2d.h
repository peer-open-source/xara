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
// FiberSection2d.h. FiberSection2d provides the abstraction of a 
// 2d beam section discretized by fibers. The section stiffness and
// stress resultants are obtained by summing fiber contributions.
//
// Written: fmk
// Created: 04/01
//
#ifndef FiberSection2d_h
#define FiberSection2d_h

#include <FrameSection.h>
#include <Vector.h>
#include <Matrix.h>
#include <memory>
#include <vector>

#define SHARE_FIBERS
class UniaxialMaterial;
class Response;

class FiberSection2d : public FrameSection
{
  public:
    FiberSection2d(); 
    FiberSection2d(int tag, int numFibers, bool compCentroid);
  private:
    FiberSection2d(const FiberSection2d&);
  public:
    ~FiberSection2d();
    FiberSection2d& operator=(const FiberSection2d&) = delete;
    FiberSection2d(FiberSection2d&&) = delete;
    FiberSection2d& operator=(FiberSection2d&&) = delete;

    const char *getClassType() const {return "FiberSection2d";}

    int   setTrialSectionDeformation(const Vector &deforms); 
    const Vector &getSectionDeformation();

    const Vector &getStressResultant();
    const Matrix &getSectionTangent();
    const Matrix &getInitialTangent();

    int   commitState();
    int   revertToLastCommit();    
    int   revertToStart();
 
    FrameSection *getFrameCopy();
    const ID &getType();
    int getOrder () const;
    
    int sendSelf(int cTag, Channel &);
    int recvSelf(int cTag, Channel &, FEM_ObjectBroker &);
    void Print(OPS_Stream &s, int flag);

    Response *setResponse(const char **argv, int argc,  OPS_Stream &s);
    int getResponse(int responseID, Information &);

    int addFiber(UniaxialMaterial &theMat, const double area, const double yLoc);

    // Sensitivity
    int setParameter(const char **argv, int argc, Parameter &);
    const Vector& getStressResultantSensitivity(int gradIndex, bool conditional);
    const Vector& getSectionDeformationSensitivity(int gradIndex);
    const Matrix& getInitialTangentSensitivity(int gradIndex);
    int commitSensitivity(const Vector& sectionDeformationGradient,  int gradIndex, int numGrads);


    double getEnergy() const; // by SAJalali
    int   getIntegral(Field field, State state, double& value) const final;

  protected:
    
  private:
    struct FiberData {
      double area;
      double y;
    };
    std::shared_ptr<std::vector<FiberData>> fibers;
    std::vector<UniaxialMaterial*> theMaterials; // array of pointers to the materials of the fibers
    double   kData[4];                 // data for ks matrix 
    double   sData[2];                 // data for s vector

    
    double QzBar, ABar, yBar;       // Section centroid
    bool computeCentroid;

    static ID code;

    Vector  e;         // trial section deformations 
    Vector  s;         // section resisting forces  (axial force, bending moment)
    Matrix  ks;        // section stiffness

    // Sensitivity
    Vector dedh; // MHS hack
};

#endif
