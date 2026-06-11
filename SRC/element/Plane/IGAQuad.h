/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
** ****************************************************************** */
//
// Description: This file contains the class definition for IGAQuad.
//

#ifndef IGAQuad_h
#define IGAQuad_h

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
//typedef struct {
//    double J2;
//    Vector WXg;
//    Vector Xi;
//    Matrix* N;
//} resultDersBasisFunsAtGPS;
//typedef struct {
//    Vector R0;
//    Matrix R1;
//} resultRationalize;
class Node;
class NDMaterial;
class Response;

class IGAQuad : public Element
{
public:
    IGAQuad(int tag, int numCPs1, int ex, int nex, int ey, int ney, int* obfs1, int ndof1,
        int* CPIds, Vector KnotVect_x1, Vector KnotVect_y1, int* numMults,
        NDMaterial& m, const char* type,
        double t, double pressure = 0.0,
        double rho = 0.0,
        double b1 = 0.0, double b2 = 0.0);
    IGAQuad();
    ~IGAQuad();

    const char* getClassType() const { return "IGAQuad"; };

    int getNumExternalNodes() const;
    const ID& getExternalNodes();
    Node** getNodePtrs();

    int getNumDOF();
    void setDomain(Domain* theDomain);

    // public methods to set the state of the element    
    int commitState();
    int revertToLastCommit();
    int revertToStart();
    int update();

    // public methods to obtain stiffness, mass, damping and residual information    
    const Matrix& getTangentStiff();
    const Matrix& getInitialStiff();
    const Matrix& getMass();

    void zeroLoad();
    int addLoad(ElementalLoad* theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector& accel);

    const Vector& getResistingForce();
    const Vector& getResistingForceIncInertia();

    // public methods for element output
    int sendSelf(int commitTag, Channel& theChannel);
    int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker
        & theBroker);

    void Print(OPS_Stream& s, int flag = 0);

    Response* setResponse(const char** argv, int argc,
        OPS_Stream& s);

    int getResponse(int responseID, Information& eleInformation);

    int setParameter(const char** argv, int argc, Parameter& param);
    int updateParameter(int parameterID, Information& info);

private:
    // private attributes - a copy for each object of the class

    NDMaterial** theMaterial; // pointer to the ND material objects

    ID connectedExternalNodes; // Tags of quad nodes

    Node* theNodes[9];

    static double matrixData[324];  // array data for matrix
    static Matrix K;		            // Element stiffness, damping, and mass Matrix
    static Vector P;		            // Element resisting force vector
    Vector Q;		                    // Applied nodal loads
    double b[2];		                // Body forces

    double appliedB[2]; // Body forces applied with load pattern, C.McGann, U.Washington
    int applyLoad;      // flag for body force in load, C.McGann, U.Washington

    Vector pressureLoad;	    // Pressure load at nodes

    double thickness;	        // Element thickness
    double pressure;	        // Normal surface traction (pressure) over entire element
                              // Note: positive for outward normal
    void setPressureLoadAtNodes();

    double rho;
    static double shp[3][9];	// Stores shape functions and derivatives (overwritten)
    static double pts[9][2];	// Stores quadrature points
    static double wts[9];		  // Stores quadrature weights
                              //

    // IGA

    //int numKnots;         // number of knots inside IGA element (to be refined later)
    Matrix* N0n;  // stores B-Spline Basis Function and 1st order 
    Matrix* N0nx; // stores shape function and derivatives in both direction
    Matrix* N0ny;
    Matrix* Ki;
    Matrix SpanIdxList; // the list of knot span index in x/y directions 
    Matrix refinedKnotVect;// the refined knot vector in x/y directions
    int numCPs;
    int ndof;
    int obfs[2];
    Vector KnotVect_x;
    Vector KnotVect_y;
    ID eleIdInfo;// stores ex, nex, ey, ney

    Vector shapeFunction(int qx, int qy, ID& eleIdInfo, Vector& KnotVect_x, Vector& KnotVect_y);

};

#endif
