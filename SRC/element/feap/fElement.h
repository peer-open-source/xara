/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
//
// Description: This file contains the class definition for fElement. fElement
// is a wrapper used to call fortran element subroutine. It is an abstract class.
//
// File: ~/element/fortran/fElement.h
// 
// Written: fmk 
// Created: 03/99
// Revision: A
//
// What: "@(#) fElement.h, revA"
//
#ifndef fElement_h
#define fElement_h

#include <Element.h>
#include <Matrix.h>

class Node;
class Channel;
class UniaxialMaterial;

class fElement : public Element
{
  public:
    // Constructors
    fElement(int tag, 
             int classTag,
             int eleType,
             int sizeD, int nen,
             int ndm, int ndf,
             int nh1, int nh3);

    fElement(int tag, 
             int classTag,
             int eleType,
             int sizeD, int nen,
             int ndm, int ndf, 
             int iow);
    
    fElement(int classTag);    
    
    // Destructor
    virtual ~fElement();

    // Public methods for element operations
    virtual int getNumExternalNodes() const;
    virtual const ID &getExternalNodes();
    Node **getNodePtrs();

    virtual int getNumDOF();        
    virtual void setDomain(Domain *theDomain);
    
    virtual int commitState();
    virtual int revertToLastCommit();        
    virtual int revertToStart();        
    virtual int update();
    
    virtual const Matrix &getTangentStiff();
    const Matrix &getInitialStiff();
    virtual const Matrix &getDamp();    
    virtual const Matrix &getMass();    

    virtual void zeroLoad();        
    virtual int addLoad(ElementalLoad *theLoad, double loadFactor);
    virtual int addInertiaLoadToUnbalance(const Vector &accel);

    virtual const Vector &getResistingForce();
    virtual const Vector &getResistingForceIncInertia();            

    // public methods for output
    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, 
                         FEM_ObjectBroker &theBroker);
    
    virtual void Print(OPS_Stream &s, int flag =0);    

  protected:
    // Protected methods 
    virtual int invokefRoutine(int ior, int iow, double *ctan, int isw);
    virtual int readyfRoutine(bool incInertia);
    virtual int invokefInit(int isw, int iow); 

    // protected data
    Vector *data;
    ID *connectedNodes;

  private:
    // private attributes - a copy for each object of the class      
    double *h;            // the elements local h array
    Node **theNodes;       // pointer to the elements nodes
    double *u;            // to hold u^(k-1)_(n+1) as nodes don't save this
    double *d;            // data stored for each element
    int eleType;          // defines which elmtnn subroutine to invoke
    int ndf;              // number of dof at nodes  nst = nen x ndf
    int nen;              // number of element nodes
    int ndm;              // dimension of mesh
    int nh1, nh3;         // the size of nh1(nh2) and nh3 data in h
    int nrCount;          // needed for Prof. Fillipou's Elmt05

    Vector *theLoad;    // vector to hold the applied load P
    Matrix *Ki;

    // static data - single copy for all objects of the class
    static Matrix **fElementM;   // class wide matrices - use s array
    static Vector **fElementV;   // class wide vectors - use r array
    static double *s;  // element tangent (nst x nst)
    static double *r;  // element residual (ndf x nen) == (nst)
    static double *ul; // nodal responses (ndf x nen x 5) == (nst X 5)
    static double *xl; // nodal coordinates (ndf x nen) == (nst)
    static double *tl; // nodal temp (nen)
    static int    *ix; // nodal tags (nen)
    static int numfElements;
};

#endif




