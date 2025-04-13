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
// SectionAggregator.  SectionAggregator decorates an MP
// section (couple bending and axial) with an uncoupled shear
// relation.
//
// Written: MHS
// Created: Jun 2000
// Revision: A
//
// What: "@(#) SectionAggregator.h, revA"

#ifndef SectionAggregator_h
#define SectionAggregator_h

#include <FrameSection.h>
#include <UniaxialMaterial.h>

#include <Vector.h>
#include <Matrix.h>

class SectionAggregator : public FrameSection
{
  public:
    SectionAggregator(); 

    SectionAggregator(int tag, FrameSection &theSection,
                      int numAdditions, UniaxialMaterial **theAdditions,
                      const ID &code); 
    SectionAggregator(int tag, int numAdditions,
                      UniaxialMaterial **theAdditions, const ID &code); 
    SectionAggregator(int tag, FrameSection &thesection,
                      UniaxialMaterial &theAddition, int c);

    ~SectionAggregator();

    const char *getClassType() const {return "SectionAggregator";};

    int   setTrialSectionDeformation(const Vector &deforms); 
    const Vector &getSectionDeformation();

    const Vector &getStressResultant();
    const Matrix &getSectionTangent();
    const Matrix &getInitialTangent();
    const Matrix &getSectionFlexibility();
    const Matrix &getInitialFlexibility();

    int   commitState();
    int   revertToLastCommit();    
    int   revertToStart();
 
    FrameSection *getFrameCopy();
    const ID &getType();
    int getOrder() const;

    int sendSelf(int cTag, Channel &theChannel);
    int recvSelf(int cTag, Channel &theChannel,
                 FEM_ObjectBroker &theBroker);

    int getIntegral(Field field, State state, double& value) const override final;

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &info);

    void Print(OPS_Stream &s, int flag =0);

    int getVariable(const char *, Information &);

    // AddingSensitivity:BEGIN //////////////////////////////////////////
    int setParameter(const char **argv, int argc, Parameter &param);
    const Vector & getStressResultantSensitivity(int gradIndex, bool conditional);
    const Vector & getSectionDeformationSensitivity(int gradIndex);
    const Matrix & getSectionTangentSensitivity(int gradIndex);
    const Matrix & getInitialTangentSensitivity(int gradIndex);
    int   commitSensitivity(const Vector& sectionDeformationGradient, int gradIndex, int numGrads);

    const Vector &getdedh(); // MHS hack
    // AddingSensitivity:END ///////////////////////////////////////////

    FrameSection* getSection() {return theSection;}

  protected:
    
  private:
    
    FrameSection *theSection;
    UniaxialMaterial **theAdditions;

    ID *matCodes;
    int numMats;
    
    Vector *e;    // Storage for section deformations
    Vector *s;    // Storage for stress resultants
    Matrix *ks;   // Storage for section stiffness
    Matrix *fs;   // Storage for section flexibility
    ID     *theCode;     // Storage for section type information
   
    int otherDbTag;

    static double workArea[];
    static int codeArea[];

// AddingSensitivity:BEGIN //////////////////////////////////////////
    Vector dedh; // MHS hack
// AddingSensitivity:END ///////////////////////////////////////////

};

#endif
