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
// Description: This file contains the class definition for Subdomain.
// Subdomain is a container class. The class is responsible for holding
// and providing access to the Elements, Nodes, LoadCases, SP_Constraints 
// and MP_Constraints that have been added to the subdomain.
//
// Written: fmk 
// Created:  11/96
// Revision: A
// Revision: B 03/98 - revised to allow parallel model generation
//
#ifndef Subdomain_h
#define Subdomain_h

#include <Domain.h>
#include <Element.h>

class Node;
class ID;
class TaggedObjectStorage;
class DomainDecompositionAnalysis;
//class PartitionedModelBuilder;
class EquiSolnAlgo;
class IncrementalIntegrator;
class LinearSOE;
class EigenSOE;
class ConvergenceTest;
class FE_Element;

#include "SubdomainNodIter.h"

class Subdomain: public Element, public Domain
{
  public:
    Subdomain(int tag);

    Subdomain(int tag, 
              TaggedObjectStorage &theInternalNodeStorage,
              TaggedObjectStorage &theExternalNodeStorage,
              TaggedObjectStorage &theElementsStorage,
              TaggedObjectStorage &theLoadPatternsStorage,	      
              TaggedObjectStorage &theMPsStorage,
              TaggedObjectStorage &theSPsStorage);
    
    virtual  ~Subdomain();    

    // method added for parallel domain generation
    //    virtual int buildSubdomain(int numSubdomains, 
    //			       PartitionedModelBuilder &theBuilder); 

    // Domain methods which must be rewritten
    virtual void clearAll();
    virtual bool addNode(Node *);	
    virtual Node *removeNode(int tag);        
    virtual NodeIter &getNodes();    
    virtual Node *getNode(int tag);            
    virtual Node **getNodePtrs();            

    virtual bool hasNode(int tag);
    virtual bool hasElement(int tag);

    virtual int getNumNodes() const;    
    virtual int commit();
    virtual int revertToLastCommit();    
    virtual int revertToStart();        
    virtual int update();
    virtual int update(double newTime, double dT);

    virtual  int barrierCheckIN() {return 0;};
    virtual  int barrierCheckOUT(int) {return 0;};
   
    virtual  void Print(OPS_Stream &s, int flag =0);
    virtual void Print(OPS_Stream &s, ID *nodeTags, ID *eleTags, int flag =0);
    
    // Domain type methods unique to a Subdomain
    virtual NodeIter &getInternalNodeIter();
    virtual NodeIter &getExternalNodeIter();
    virtual bool addExternalNode(Node *);

    virtual void wipeAnalysis();
    virtual void setDomainDecompAnalysis(DomainDecompositionAnalysis &);
    virtual int setAnalysisAlgorithm(EquiSolnAlgo &);
    virtual int setAnalysisIntegrator(IncrementalIntegrator &);
    virtual int setAnalysisLinearSOE(LinearSOE &);
    virtual int setAnalysisEigenSOE(EigenSOE &);
    virtual int setAnalysisConvergenceTest(ConvergenceTest &);
    virtual int invokeChangeOnAnalysis();

    // Element methods which must be written
    virtual int getNumExternalNodes() const;    
    virtual const ID &getExternalNodes();
    virtual int getNumDOF();

    virtual int commitState();    
    
    virtual const Matrix &getTangentStiff();
    virtual const Matrix &getInitialStiff();    
    virtual const Matrix &getDamp();    
    virtual const Matrix &getMass();    

    virtual void  zeroLoad();
    virtual int addLoad(ElementalLoad *theLoad, double loadFactor);
    virtual int addInertiaLoadToUnbalance(const Vector &accel);

    virtual const Vector &getResistingForce();    
    virtual const Vector &getResistingForceIncInertia();        
    virtual bool isSubdomain();    
    virtual int setRayleighDampingFactors(double alphaM, 
					  double betaK, 
					  double betaK0, 
					  double betaKc);

//  virtual  int  updateParameter(int tag, int value);
//  virtual  int  updateParameter(int tag, double value);    

    // Element type methods unique to a subdomain
    virtual int computeTang();
    virtual int computeResidual();
    virtual const Matrix &getTang();    

    void setFE_ElementPtr(FE_Element *theFE_Ele);
    virtual const Vector &getLastExternalSysResponse();
    virtual int computeNodalResponse();    
    virtual int analysisStep(double deltaT);
    virtual int eigenAnalysis(int numMode, bool generalized, bool findSmallest);
    virtual bool doesIndependentAnalysis();

    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    virtual double getCost();
    virtual int addResistingForceToNodalReaction(int inclInertia);
    
  protected:    
    virtual int buildMap();
    bool mapBuilt;
    ID *map;
    Vector *mappedVect;
    Matrix *mappedMatrix;


    FE_Element *getFE_ElementPtr();
    TaggedObjectStorage  *internalNodes;
    TaggedObjectStorage  *externalNodes;    

    DomainDecompositionAnalysis *getDDAnalysis();

  private:
    double realCost;
    double cpuCost;
    int pageCost;
    DomainDecompositionAnalysis *theAnalysis;
    ID *extNodes;
    FE_Element *theFEele;
    
    //    TaggedObjectStorage  *realExternalNodes;        

    SingleDomNodIter   *internalNodeIter;
    SingleDomNodIter   *externalNodeIter;    
    SubdomainNodIter   *theNodIter;

    //    PartitionedModelBuilder *thePartitionedModelBuilder;
    static Matrix badResult;
};

#endif


