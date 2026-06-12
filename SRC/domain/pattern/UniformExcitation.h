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
// Written: fmk 11/98
// Revised:
//
// Purpose: This file contains the class definition for UniformExcitation.
// UniformExcitation is a concrete class. It sets the R for a single
// ground motion acting on a structure.

#ifndef UniformExcitation_h
#define UniformExcitation_h

#include <LoadPattern.h>

class UniformExcitation : public LoadPattern
{
  public:
    UniformExcitation();  
    UniformExcitation(GroundMotion *, 
                      int ndm,
                      int dof, 
                      int tag, 
                      double vel0 = 0.0, 
                      double fact = 1.0);  
    ~UniformExcitation();
    int getDirection() {return theDof;}

    const char* getClassType() const override {return "UniformExcitation";}

    void setDomain(Domain *) override;
    void applyLoad(double time) override;
    int applyResidual(AnalysisModel &, LinearSOE &, double) override;

    bool addSP_Constraint(SP_Constraint *) override {return false;}

    void Print(OPS_Stream &s, int flag);

    int sendSelf(int tag, Channel &) override;
    int recvSelf(int tag, Channel &, FEM_ObjectBroker &) override;

    virtual int setParameter(const char **argv, int argc, Parameter &);
    virtual int  updateParameter(int parameterID, Information &);
    virtual int  activateParameter(int parameterID);
    void applyLoadSensitivity(double time);
    
    const GroundMotion *getGroundMotion();

    
 private:
    GroundMotion *theMotion; // the ground motion
    int theDof;      // the dof corrseponding to the ground motion
    double vel0;     // the initial velocity, should be neg of ug dot(0)
    double fact;
    int NDM;
    double currentTime;
    int parameterID;
};

#endif
