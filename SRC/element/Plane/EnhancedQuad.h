//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
#include <array>
#include <ID.h> 
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <NDMaterial.h>
#include <quadrature/Plane/LegendreFixedQuadrilateral.h>

namespace OpenSees {template<int n, int m, typename T> struct MatrixND;};

class EnhancedQuad : public Element,
                     protected LegendreFixedQuadrilateral<4>
{
  public:

    EnhancedQuad(int tag, 
                  std::array<int,4>& nodes,
                  NDMaterial &theMaterial, 
                  double thickness);

    EnhancedQuad();

    ~EnhancedQuad();

    const char *getClassType() const {return "EnhancedQuad";}

    void setDomain( Domain *theDomain );

    int getNumExternalNodes( ) const ;
    const ID &getExternalNodes( ) ;
    Node **getNodePtrs();

    int getNumDOF( ) ;

    // methods dealing with state updates
    int commitState( ) ;
    int revertToLastCommit( ) ;
    int revertToStart( ) ;
    int update();
	
    //return stiffness matrix 
    const Matrix &getTangentStiff();
    const Matrix &getInitialStiff();
    const Matrix &getMass();

    // zero the load -- what load?
    void zeroLoad();
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);


    const Vector &getResistingForce();    
    const Vector &getResistingForceIncInertia();

    // public methods for element output
    Response *setResponse(const char **argv, int argc, OPS_Stream &s);

    int getResponse(int responseID, Information &eleInformation);
    int sendSelf (int commitTag, Channel &);
    int recvSelf (int commitTag, Channel &, FEM_ObjectBroker &);

    void Print( OPS_Stream &s, int flag );

  private:
    constexpr static int NEN = 4;
    static constexpr int ndm = 2 ;
    static constexpr int NDF = 2 ; 
    static constexpr int nstress = 3 ; 
    static constexpr int numberNodes = 4 ;
    static constexpr int nip = 4 ;
    static constexpr int nEnhanced = 4 ; 
    static constexpr int nModes = 2 ;
    static constexpr int numberDOF = 8 ;
    static constexpr int nShape = 3 ;
    static constexpr int massIndex = nShape - 1 ;

    //static data
    static Matrix stiff ;
    static Vector resid ;
    static Matrix mass ;
    static Matrix damping ;

    // //stress data
    // static double stressData[][4] ;

    // //tangent data 
    // static double tangentData[][3][4] ;

    ID connectedExternalNodes ;  // node tags
    Node *theNodes[NEN] ;        // pointers to nodes

    Vector alpha ; // enhanced strain parameters

    double thickness ;

    //material information
    NDMaterial *materialPointers[4] ; //pointers to four materials
					  
    //local nodal coordinates, two coordinates for each of four nodes
    //    static double xl[2][4] ; 
    static double xl[][4] ; 

    //compute local coordinates and basis
    void computeBasis( ) ;
        
    //
    int formResidAndTangent( int tang_flag ) ;

    void formInertiaTerms( int tangFlag ) ;

    void computeB(int node, const double shp[3][4], 
                          OpenSees::MatrixND<3,2,double> &B);

    //shape function routine for four node quads
    void shape2d( double ss, double tt, 
		  const double x[2][4], 
		  double shp[3][4], 
		  double &xsj ) ;

    Vector *load;
    Matrix *Ki;
} ; 


