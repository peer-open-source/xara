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
// Description: This file contains the class definition for AnalysisModel.
// AnalysisModel is a container class. This class is responsible for holding
// and providing access to the FE_Element and DOF_Group objects that the 
// ConstraintHandler creates. It is also responsible for updating the 
// response quantities at the DOF_Groups and for triggering methods 
// in the associated Domain.
//
// Written: fmk 
// Created: 9/96
// Revision: A
//
#ifndef AnalysisModel_h
#define AnalysisModel_h

#include <MovableObject.h>
#define VIRTUAL

class TaggedObjectStorage;
class Domain;
class FE_EleIter;
class DOF_GrpIter;
class Graph;
class FE_Element;
class DOF_Group;
class Vector;
class FEM_ObjectBroker;
class ConstraintHandler;

class AnalysisModel: public MovableObject
{
  public:
    AnalysisModel();    
    AnalysisModel(int classTag);
    AnalysisModel(TaggedObjectStorage &, TaggedObjectStorage &);

    VIRTUAL ~AnalysisModel();

    // methods to populate/depopulate the AnalysisModel
    VIRTUAL bool addFE_Element(FE_Element *theFE_Ele);
    VIRTUAL bool addDOF_Group(DOF_Group *theDOF_Grp); // called by Handler
    VIRTUAL void clearAll();
    VIRTUAL void clearDOFGraph();                 // called by Numberer and Analysis
    VIRTUAL void clearDOFGroupGraph();

    // methods to access the FE_Elements and DOF_Groups and their numbers
    VIRTUAL int getNumDOF_Groups(void) const;		
    VIRTUAL DOF_Group *getDOF_GroupPtr(int tag);	
    VIRTUAL FE_EleIter &getFEs();
    VIRTUAL DOF_GrpIter &getDOFs();

    // method to access the connectivity for SysOfEqn to size itself
    VIRTUAL void   setNumEqn(int) ;	
    VIRTUAL int    getNumEqn() const ; 
    VIRTUAL Graph &getDOFGraph();
    VIRTUAL Graph &getDOFGroupGraph();
    
    // methods to update the response quantities at the DOF_Groups,
    // which in turn set the new nodal trial response quantities.
    VIRTUAL void setResponse(const Vector &disp, 
			     const Vector &vel, 
			     const Vector &accel);
    VIRTUAL void setDisp(const Vector &disp);
    VIRTUAL void setVel(const Vector &vel);
    VIRTUAL void setAccel(const Vector &vel);            

    VIRTUAL void incrDisp(const Vector &disp);    
    VIRTUAL void incrVel(const Vector &vel);

    int getState(Vector&, Vector&, Vector&, int flag); // cmp

//  VIRTUAL void incrAccel(const Vector &vel);            

    // methods added to store the eigenvalues and vectors in the domain
    VIRTUAL void setNumEigenvectors(int numEigenvectors);
    VIRTUAL void setEigenvector(int mode, const Vector &);
    VIRTUAL void setEigenvalues(const Vector &);    
    VIRTUAL const Vector &getEigenvalues(void);    
    const Vector *getModalDampingFactors(void);
    bool inclModalDampingMatrix(void);
    
    // methods which trigger operations in the Domain
    VIRTUAL void setLinks(Domain &theDomain, ConstraintHandler &theHandler);
	
    VIRTUAL void   applyLoadDomain(double newTime);
    VIRTUAL int    updateDomain(void);
    VIRTUAL int    updateDomain(double newTime, double dT);
    VIRTUAL int    analysisStep(double dT =0.0);
    VIRTUAL int    eigenAnalysis(int numMode, bool generalized, bool findSmallest);
    VIRTUAL int    commitDomain(void);
    VIRTUAL int    revertDomainToLastCommit(void);
    VIRTUAL double getCurrentDomainTime(void);
    VIRTUAL void   setCurrentDomainTime(double newTime);    
    VIRTUAL void   setRayleighDampingFactors(double alphaM, double betaK, double betaKi, double betaKc);    
    
    // Parallel
    VIRTUAL int sendSelf(int commitTag, Channel &theChannel);
    VIRTUAL int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);

    Domain *getDomainPtr(void) const;

  protected:

    
  private:
    Domain *myDomain;
    ConstraintHandler *myHandler;

    Graph *myDOFGraph;
    Graph *myGroupGraph;    
    
    int numFE_Ele;             // number of FE_Elements objects added
    int numDOF_Grp;            // number of DOF_Group objects added
    int numEqn;                // numEqn set by the ConstraintHandler typically

    TaggedObjectStorage  *theFEs;
    TaggedObjectStorage  *theDOFs;
    
    FE_EleIter    *theFEiter;     
    DOF_GrpIter   *theDOFiter;    
};

#endif
