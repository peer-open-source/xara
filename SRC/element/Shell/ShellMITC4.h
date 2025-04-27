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
// Original implementation: Ed "C++" Love
// Reimplementation: Leopoldo Tesser, Diego A. Talledo, V�ronique Le Corvec
//
// Bathe MITC 4 four node shell element with membrane and drill
// Ref: Dvorkin,Bathe, A continuum mechanics based four node shell
//      element for general nonlinear analysis,
//      Eng.Comput.,1,77-88,1984
//
#ifndef ShellMITC4_h
#define ShellMITC4_h

#include <ID.h> 
#include <array>
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <quadrature/GaussLegendre2D.hpp>

class SectionForceDeformation;

namespace OpenSees {template<int n, int m, typename T> struct MatrixND;};

class ShellMITC4 : public    Element,
                   protected GaussLegendre<2,4>
{
 public:
  
    // null constructor
    ShellMITC4( );
    
    // full constructor
    ShellMITC4( int tag, 
                int node1,
                int node2,
                int node3,
                int node4,
                SectionForceDeformation &theMaterial,
                bool updateBasis=false) ;
    
    virtual ~ShellMITC4( ) ;

    const char* getClassType() const final {
      return "ShellMITC4";
    }

    void setDomain( Domain *theDomain ) ;
    
    // get the number of external nodes
    int getNumExternalNodes( ) const ;
    
    // return connected external nodes
    const ID &getExternalNodes( ) ;
    Node **getNodePtrs( );

    // return number of dofs
    int getNumDOF( ) ;

    // commit state
    int commitState( ) ;
    
    // revert to last commit 
    int revertToLastCommit( ) ;
    
    // revert to start 
    int revertToStart( ) ;

    //print out element data
    void Print( OPS_Stream &s, int flag ) ;
	
    //return stiffness matrix 
    const Matrix &getTangentStiff( ) ;
    const Matrix &getInitialStiff( );
    const Matrix &getMass( );

    // methods for applying loads
    void zeroLoad( void );	
    int addLoad( ElementalLoad *theLoad, double loadFactor );
    int addInertiaLoadToUnbalance( const Vector &accel );

    //get residual
    const Vector &getResistingForce( ) ;
    
    //get residual with inertia terms
    const Vector &getResistingForceIncInertia( ) ;

    // public methods for element output
    int sendSelf ( int commitTag, Channel &theChannel );
    int recvSelf ( int commitTag, Channel &theChannel, FEM_ObjectBroker 
		           &theBroker );

    Response* setResponse( const char **argv, int argc, OPS_Stream &output );
    int getResponse( int responseID, Information &eleInfo );

    int setParameter(const char **argv, int argc, Parameter &param);


  private : 
    static const int ndf = 6;     // two membrane plus three bending plus one drill
    static const int nstress = 8; // three membrane, three moment, two shear
//  static const int NIP = 4;     // suppllied by quadrature
    static const int NEN = 4;

    //static data
    static Matrix stiff ;
    static Vector resid ;
    static Matrix mass ;
    static Matrix damping ;

    //node information
    ID connectedExternalNodes ;            // node tags
    std::array<Node *, NEN> theNodes;      // pointers to nodes

    // drilling stiffness
    double Ktt ;

    // material information
    std::array<SectionForceDeformation *, nip> materialPointers; // pointers to materials
					  
    //local nodal coordinates, two coordinates for each of four nodes
    //static double xl[][4] ; 
    double xl[2][4] ; 

    // shell basis vectors
    double g1[3] ;
    double g2[3] ;
    double g3[3] ;

    int applyLoad;
    double appliedB[3];     // Body forces applied with load


    //compute local coordinates and basis
    void computeBasis( ) ;
    //start Yuli Huang (yulihuang@gmail.com) & Xinzheng Lu (luxz@tsinghua.edu.cn)
    bool doUpdateBasis;
    void updateBasis( ) ;
    //end Yuli Huang (yulihuang@gmail.com) & Xinzheng Lu (luxz@tsinghua.edu.cn)
        
    // inertia terms
    void formInertiaTerms( int tangFlag ) ;

    // form residual and tangent					  
    void formResidAndTangent( int tang_flag ) ;

    // compute Bdrill matrix
    double* computeBdrill( int node, const double shp[3][4] , double Bdrill[6]);

    // assemble a B matrix 
    void
    assembleB( const Matrix &Bmembrane,
               const Matrix &Bbend, 
               const Matrix &Bshear,
               OpenSees::MatrixND<nstress, ndf, double> &B
             ) ;
  
    // compute Bmembrane matrix
    const Matrix& computeBmembrane( int node, const double shp[3][4] ) ;
  
    // compute Bbend matrix
    const Matrix& computeBbend( int node, const double shp[3][4] ) ;
  
    // shape function routine for four node quads
    void shape2d( double ss, double tt, 
		  const double x[2][4], 
		  double shp[3][4], 
		  double &xsj ) ;

    // vector for applying loads
    Vector *load;
    Matrix *Ki;
    double init_disp[4][6];
}; 

#endif
