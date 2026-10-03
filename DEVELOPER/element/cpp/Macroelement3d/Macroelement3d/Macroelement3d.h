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

/*
Source: /DEVELOPER/element/cpp/Macroelement3d/Macroelement3d/Macroelement3d.h
Written by Francesco Vanin (francesco.vanin@epfl.ch)
Ecole Polytechnique Federale de Lausanne, Switzerland,
Earthquake Engineering and Structural Dynamics laboratory, 2019

Reference: Vanin F., Penna A., Beyer K.;"A three dimensional macro-element
for modelling of the in-plane and out-of-plane response of masonry walls",
submitted to Earthquake Engineering and Structural Dynamics (2019)

Last edit: 27 Feb 2019
*/

#ifndef Macroelement3d_h
#define Macroelement3d_h

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <NDMaterial.h>
#include <UniaxialMaterial.h>
#include <vector>
#include <string>

class Node;
class SectionForceDeformation;
class NDMaterial;
class Response;

class Macroelement3d : public Element
{
  public:
    Macroelement3d(int tag, int nd1, int nd2, int ndE, SectionForceDeformation *sI, SectionForceDeformation *sE, SectionForceDeformation *sJ, NDMaterial* shearModel, NDMaterial* shearModelOOP, 
		double h, double E_, 
		Vector driftF, Vector driftF_ALR, Vector driftS, Vector driftS_ALR, double Ltfc, double alphaAC_HC, double betaShearSpanF, double betaShearSpanS, double failureFactorF, double failureFactorS,
		Vector axis, Vector oop, Vector intLength, Vector intLengthMasses, Vector massDir, 
		int PDelta=0, double rho=0.0, int cm=0.0, double isGable=false, int dampingModel=0);
    Macroelement3d();
    ~Macroelement3d();

    const char *getClassType(void) const {return "Macroelement3d";};

    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);

    int getNumDOF(void);
    void setDomain(Domain *theDomain);

    // public methods to set the state of the element    
    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);

    // public methods to obtain stiffness, mass, damping and residual information    
    int update(void);
	Vector getBasicDisplacement(Vector uI, Vector uJ, Vector uE);
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);
    const Matrix &getMass(void);

    void zeroLoad();
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);    

    // public methods for element output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag =0);

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInfo);

	//drift model
	void driftModel(double currentDriftF, double currentDriftS, double axialLoadRatio, double H0overL=1.0);

	// damping methods
	int setRayleighDampingFactors(double alphaM, double betaK, double betaK0, double betaKc);
	const Matrix &getDamp(void);
	const Matrix &getSecantStiff(void);

  protected:
	  const Vector &getRayleighDampingForces(void);
    
  private:
    // builds initial basic stiffness. Zero stiffness if the element failed for excessive drift demand 
    const Matrix &getInitialBasicStiff(void);

	// incremental compatibility matrix (derivative of P-delta compatibility equations)
	// its transpose is used as incremental equilibrium matrix
    const Matrix &getIncrementalCompatibilityMatrix(bool flagIncremental=true);

	// explicit matrix operations to transform from local to global accounting for node offset(s)
	int trasformMatrixToGlobal(Matrix& A);                                           

    const int numSections;                      // number of sections = 3
	const int numShearModels;                   // number of shear models = 2

    SectionForceDeformation **theSections;      // pointers to sectional models
	NDMaterial **theShearModel;                 // pointers to shear models

    ID connectedExternalNodes;                  // Tags of elements nodes (i,j,e)
    Node *theNodes[3];                          // pointers to element nodes (i,j,e)

    static Matrix K;                            // stores element stiffness, damping, or mass Matrix.
    Matrix GammaC;                              // Compatibility matrix, including or not P-Delta effects
    static Vector P;                            // Element resisting force vector
    Matrix Tgl, Tgl6;                           // Rotation 3d matrix from global to local coordinates. Full version and 6x6 version to rotate only one node dofs
    double R[3][3];                            

    Vector Q;                                   // Applied nodal loadsalphaNC_SD
    Vector q;                                   // Basic force
    Vector uBasic, uBasicCommitted;             // Basic displcaments

    double q0[12];                              // Element forces producing a p-delta effect
    double p0[12];                              // Reactions in basic system of distributed loads

    double rho;                                 // Mass density per unit length
    int cMass;                                  // consistent mass flag
	Vector massglobalDir;                       // participation factors for masses (in the global system)
    int PDelta;                                 // flag for p-delta formulation (default: false)
    double deltaW1, deltaV1, deltaW3, deltaV3;  // relative displacements used to define the compatibility matrix with P-Delta transformation

	int parameterID;
    enum {maxNumSections = 20}; // not used 

    static double workArea[];  // not used here

    double* nodeIInitialDisp;                    // initial displacements of the end nodes 
    double* nodeJInitialDisp; 

    double* nodeIOffset;                         // node offsets of the end nodes
    double* nodeJOffset;

    Vector xAxis, yAxis, zAxis;                  // local orientation axes expressed in global coordinates
	Vector intLength, intLengthMasses;           // integration lengths for sectional responses and lumped masses
    double L;                                    // half-element length
	double E;                                    // elastic modulus

	// drift model parameters
	double driftF, driftS;                       // current flexural and shear drift                         
	double Ltfc;                                 // product L*t*fc, maximum axial load (used for axial load ratio, only for the drift model)
	double betaShearSpanF, betaShearSpanS;       // exponent beta applied to the shear span ratio in the drift model
	Vector limitDriftF, limitDriftS;             // given points of the drift/axial load ratio relationship
	Vector limitDriftF_ALR, limitDriftS_ALR;     // corresponding axial load ratios
	double alphaAC_HC;                           // factor defining axial collapse drift as a function of lateral load collapse drift

	bool failedF, failedS;                       // flags for loss of lateral force capacity (current and committed) 
	bool failedFcommitted, failedScommitted;
	bool collapsedF, collapsedS;                 // flag for loss of axial force capacity (current and committed)
	bool collapsedFcommitted, collapsedScommitted;

	double failureFactorF,failureFactorS;        // factor defining the loss of force capacity
	
	bool isGable;                                // flag: gable element (used for mass matrix) 
	double wx, wy, wz;                           // applied distributed element loads

	double ductilityDemand, ductilityDemandCommitted, ductilityDemandPeak, ductilityDemandCycle; // ductility demands for secant stiffness - damping model
	bool triggerDuctility;
	Vector CommittedDeformation[3];    // committed section deformations
	int dampingModel;
               

};

#endif

