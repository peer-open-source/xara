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
    ~AnalysisModel();

    void setLinks(Domain &, ConstraintHandler &);

    // methods to populate/depopulate the AnalysisModel
    bool addFE_Element(FE_Element *theFE_Ele);
    bool addDOF_Group(DOF_Group *theDOF_Grp); // called by Handler
    void clearAll();
    void clearDOFGraph();                 // called by Numberer and Analysis
    void clearDOFGroupGraph();

    // methods to access the FE_Elements and DOF_Groups and their numbers
    int getNumDOF_Groups() const;		
    DOF_Group *getDOF_GroupPtr(int tag);	
    FE_EleIter &getFEs();
    DOF_GrpIter &getDOFs();

    // Access the connectivity for SysOfEqn to size itself
    void   setNumEqn(int) ;	
    int    getNumEqn() const; 
    Graph &getDOFGraph();
    Graph &getDOFGroupGraph();
    
    // Update the response quantities at the DOF_Groups,
    // which in turn set the new nodal trial response quantities.
    void setResponse(const Vector &disp, const Vector &vel, const Vector &accel);
    void setDisp(const Vector &disp);
    void setVel(const Vector &vel);
    void setAccel(const Vector &vel);            

    void incrDisp(const Vector &disp);    
    void incrVel(const Vector &vel);
//  void incrAccel(const Vector &vel);

    int getState(Vector&, Vector&, Vector&, int flag); // cmp
      


    // Store the eigenvalues and vectors in the domain
    void setNumEigenvectors(int numEigenvectors);
    void setEigenvector(int mode, const Vector &);
    void setEigenvalues(const Vector &);
    const Vector *getModalDampingFactors();
    bool inclModalDampingMatrix();

    
    // Useful methods
    void   applyLoadDomain(double newTime);
    int    updateDomain();
    int    updateDomain(double newTime, double dT);

    // Simple wrappers for the Domain methods; remove these!
    int    commitDomain();
    int    analysisStep(double dT =0.0);
    double getCurrentDomainTime();
    void   setCurrentDomainTime(double newTime);
    void   setRayleighDampingFactors(double alphaM, double betaK, double betaKi, double betaKc);
    const Vector &getEigenvalues();
    #if 0
    int    eigenAnalysis(int numMode, bool generalized, bool findSmallest);
    int    revertDomainToLastCommit();
    #endif

    // Parallel
    int sendSelf(int commitTag, Channel &);
    int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &);

    Domain *getDomainPtr() const;
    
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
