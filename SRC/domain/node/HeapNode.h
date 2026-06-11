//
#ifndef HeapNode_h
#define HeapNode_h


#include <Node.h>

// TODO: Remove include of NodeData
#include "NodeData.h"

class Vector;
class Matrix;
class Channel;
class DOF_Group;
class NodalThermalAction; //L.Jiang [ SIF ]
class Domain;
class Element;

class HeapNode : public Node
{
  public:

    // constructors
    HeapNode(int classTag);
    HeapNode(int tag, int classTag);
    HeapNode(int tag, int ndof, double Crd1);
    HeapNode(int tag, int ndof, double Crd1, double Crd2);
    HeapNode(int tag, int ndof, double Crd1, double Crd2, double Crd3);
    HeapNode(const Node &theCopy, bool copyMass = true);
    
    // destructor
    virtual ~HeapNode();

    // public methods dealing with the DOF at the node
    virtual int  getNumberDOF() const;

    // public methods for obtaining the nodal coordinates
    virtual const Vector &getCrds() const;
    virtual void  setCrds(const Vector &);

    //
    // State
    //
    // public methods dealing with the committed state of the node
    virtual int commitState() noexcept;
    virtual int revertToLastCommit();
    virtual int revertToStart();

    //
    // Response
    //
    // public methods for obtaining committed and trial 
    // response quantities of the node
    virtual const Vector &getDisp() noexcept;
    virtual const Vector &getIncrDisp() noexcept;
    virtual const Vector &getIncrDeltaDisp() noexcept;
    virtual const Vector &getTrialDisp() noexcept;

    virtual const Vector &getVel();
    virtual const Vector &getAccel();
    virtual const Vector &getTrialVel();
    virtual const Vector &getTrialAccel();
    virtual const Vector &getTotalAccel();
//  virtual const Vector &getRigidAccel();

    // public methods for updating the trial response quantities
    virtual int setTrialDisp  (double value, int dof) noexcept final;
    virtual int setTrialDisp  (const Vector &) noexcept final;
    virtual int incrTrialDisp (const Vector &) final;
    virtual int setTrialVel   (const Vector &) ;
    virtual int setTrialAccel (const Vector &) ;
    virtual int incrTrialVel  (const Vector &) ;
    virtual int incrTrialAccel(const Vector &) ;

    // Dynamics
    virtual const Matrix &getMass() ;
    virtual const Matrix &getDamp() final;
    virtual int setMass(const Matrix &theMass) final;
    virtual int setNumColR(int numCol) final;
    virtual int setR(int row, int col, double Value) final;
    virtual const Vector &getRV(const Vector &V) final;
    virtual int setRayleighDampingFactor(double alphaM) final;

    // Eigen vectors
    virtual int setNumEigenvectors(int numVectorsToStore);
    virtual int setEigenvector(int mode, const Vector &eigenVector);
    virtual const Matrix &getEigenvectors();
 
    //
    // Load information
    //
    virtual void zeroUnbalancedLoad() noexcept;
    virtual int addRigidAccleration(const Vector& accel, double fact);
    virtual int addUnbalancedLoad(const Vector &load, double fact = 1.0);
    virtual int addInertiaLoadToUnbalance(const Vector &accel, double fact = 1.0);
    virtual const Vector &getUnbalancedLoad();
    virtual const Vector &getUnbalancedLoadIncInertia();


    // Misc. Responses
    virtual const Vector &getReaction();
    virtual int   addReactionForce(const Vector &, double factor);
    virtual int   resetReactionForce(int flag);

    virtual const Vector *getResponse(NodeData);
    int fillResponse(NodeData responseType, Vector& result, int offset=0);
    
    //
    // Parallel
    //
    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, 
                         FEM_ObjectBroker &theBroker);

    //
    // Misc
    //
    virtual void Print(OPS_Stream &s, int flag = 0);


    // AddingSensitivity:BEGIN /////////////////////////////////////////
    int addInertiaLoadSensitivityToUnbalance(const Vector &accel, 
                                             double fact = 1.0, 
                                             bool tag=false);
    Matrix getMassSensitivity();
    virtual const Matrix &getDampSensitivity();
    int    getCrdsSensitivity();
    int    saveDispSensitivity(const Vector &v, int gradNum, int numGrads);
    int    saveVelSensitivity(const Vector &vdot, int gradNum, int numGrads);
    int    saveAccelSensitivity(const Vector &vdot, int gradNum, int numGrads);
    double getDispSensitivity(int dof, int gradNum);
    double getVelSensitivity(int dof, int gradNum);
    double getAccSensitivity(int dof, int gradNum);
    int    setParameter(const char **argv, int argc, Parameter &param);
    int    updateParameter(int parameterID, Information &info);
    int    activateParameter(int parameterID);
    // AddingSensitivity:END ///////////////////////////////////////////


  private:
    double *disp;
    double *vel, *accel;               // double arrays holding the vel and accel values
    Vector *trialDisp, *commitDisp, *incrDisp, *incrDeltaDisp;
    Vector *commitVel, *commitAccel;   // committed quantities
    Vector *trialVel,   *trialAccel;   // trial quantities
    Vector *totalAccel, *rigidAccel;               //
    Vector *unbalLoad;                // unbalanced load


#if 1
    Domain* theDomain;
#endif

    // Private global state
    int setGlobalMatrices();

    static Matrix **theMatrices;
    static int numMatrices;
    static Matrix **theVectors;
    static int numVectors;
    int index;


    // priavte methods used to create the Vector objects 
    // for the committed and trial response quantaties.
    virtual int createDisp() final;
    virtual int createVel() final;
    virtual int createAccel() final;

    // private data associated with each node object
    int numberDOF;                    // number of dof at Node

    Vector *Crd;                      // original nodal coords

    

    Matrix *R;                          // nodal participation matrix
    Matrix *mass;                       // pointer to mass matrix
    Vector *unbalLoadWithInertia;
    double alphaM;                      // rayleigh damping factor 
    Matrix *theEigenvectors;


    int dbTag1, dbTag2, dbTag3, dbTag4; // needed for database

    // AddingSensitivity:BEGIN /////////////////////////////////////////
    Matrix *dispSensitivity;
    Matrix *velSensitivity;
    Matrix *accSensitivity;
    int parameterID;
    // AddingSensitivity:END ///////////////////////////////////////////

    Vector *reaction;
};

#endif

