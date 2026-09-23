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
                                                                        
// ============================================================================
// 2018 By Jose Abell @ Universidad de los Andes, Chile
// www.joseabell.com | https://github.com/jaabell | jaabell@miuandes.cl
// ============================================================================
// Implements a standard 4-node tetrahedron element.
//
// This element has one Gauss point of integration that represents the stress or
// strain field in the whole element. This is a constant stress/strain element
// because the displacement interpolation is linear. 
//
// This element is very succeptible to locking, so use with fine discretizations
// requires care. 
//
// ============================================================================

#pragma once
#include <array>
#include <stdlib.h> 
#include <math.h>

#include <ID.h> 
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <NDMaterial.h>

class FourNodeTetrahedron : public Element {

public:
  
    // full constructor
    FourNodeTetrahedron(int tag, 
      const std::array<int, 4> &nodes,
      NDMaterial &theMaterial,
      double b1, double b2, double b3);

    // destructor 
    virtual ~FourNodeTetrahedron();

    const char *getClassType() const {return "FourNodeTetrahedron";}

    void setDomain( Domain *theDomain ) ;

    // get the number of external nodes
    int getNumExternalNodes( ) const ;

    // return connected external nodes
    const ID &getExternalNodes( ) ;
    Node **getNodePtrs();

    // return number of dofs
    int getNumDOF();

    // commit state
    int commitState();
    
    // revert to last commit 
    int revertToLastCommit( ); 
    int revertToStart();
    int update();

    // print out element data
    void Print( OPS_Stream &s, int flag ) ;
	
    // return stiffness matrix 
    const Matrix &getTangentStiff();
    const Matrix &getInitialStiff();    
    const Matrix &getMass();    

    void zeroLoad( ) ;
    int  addLoad(ElementalLoad *theLoad, double loadFactor);
    int  addInertiaLoadToUnbalance(const Vector &accel);

    // get residual
    const Vector &getResistingForce( );
    
    // get residual with inertia terms
    const Vector &getResistingForceIncInertia( );

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInformation);


    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);

    void onActivate();
    void onDeactivate();

  private : 

    // Number of Gauss-points  
    static constexpr int NumGaussPoints = 1,
                         NumNodes = 4,
                         NumDOFsPerNode = 3,
                         NumStressComponents = 6,
                         NumDOFsTotal = NumNodes*NumDOFsPerNode;
    //
    // private attributes
    //
    ID connectedExternalNodes ;  // node tags
    Node *nodePointers[4] ;      // pointers to nodes

    //material information
    NDMaterial *materialPointers[1]; //pointers to eight materials

    double b[3];		// Body forces
    double appliedB[3];		// Body forces applied with load
    int applyLoad;

    Vector *load;
    Matrix *Ki;

    //
    // static attributes
    //
    static Matrix stiff ;
    static Vector resid ;
    static Matrix mass ;
    static Matrix damping ;
    static Matrix B;
  
    // quadrature data
    static const double root3 ;
    static const double one_over_root3 ;    
    static const double sg[1] ;
    static const double wg[1] ;
  
    // local nodal coordinates, three coordinates for each of four nodes
    static double xl[3][NumNodes] ; 

    Vector initDisp[NumNodes];

    int do_update;

    //
    // private methods
    //

    // inertia terms
    void formInertiaTerms( int tangFlag ) ;

    // form residual and tangent					  
    void formResidAndTangent( int tang_flag ) ;

    // compute coordinate system
    void computeBasis();

    // compute B matrix
    const Matrix& computeB( int node, const double shp[4][NumNodes] ) ;

    void shp3d( const double ss[4], double &xsj, double shp[4][4], const double xl[3][4]   );

};
