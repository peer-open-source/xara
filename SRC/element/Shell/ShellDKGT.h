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
// Written: Shuhao Zhang & Xinzheng Lu
//
// Three node flat shell element with membrane and drill DOF
//
// Ref: Plate Bending Part - DKT, thin plate element
//      Membrane Part - GT9, a membrane element with drilling DOF

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <ID.h> 
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <SectionForceDeformation.h>

class ShellDKGT : public Element {

 public:
  
  //null constructor
  ShellDKGT( ); 

  ShellDKGT( int tag, 
         int node1,
         int node2,
         int node3,
         SectionForceDeformation &theMaterial, double b1, double b2, double b3) ;
  
  // destructor 
  virtual ~ShellDKGT( ) ;

  //set domain because frank is a dumb ass 
  void setDomain( Domain *theDomain ) ;
  
  //get the number of external nodes
  int getNumExternalNodes( ) const ;
    
  //return connected external nodes
  const ID &getExternalNodes( ) ;
  Node **getNodePtrs( );

  int getNumDOF( );
  int commitState( );
  int revertToLastCommit( );
  int revertToStart( );

  // print out element data
  void Print( OPS_Stream &s, int flag ) ;
  
  //return stiffness matrix 
  const Matrix &getTangentStiff( ) ;
  const Matrix &getInitialStiff( );
  const Matrix &getMass( );

  // methods for applying loads
  void zeroLoad();    
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

  private:
    constexpr static int nip = 4;
    static constexpr int ndf         = 6;
    static constexpr int numberNodes = 3;
    static constexpr int numberGauss = 4;
    static constexpr int nShape      = 3;
    static constexpr int massIndex   = nShape - 1;

    //static data
    static Matrix stiff ;
    static Vector resid ;
    static Matrix mass ;
    static Matrix damping ;

    //quadrature data
    static const double three ;
    static const double one;
    static const double one_over_three ;
    static const double five;
    static const double one_over_five;
    static const double three_over_five;
    static const double one_over_four;
    static const double wg1;
    static const double wg2;
    static double sg[nip] ;
    static double tg[nip] ;
    static double qg[nip] ;
    static double wg[nip] ;

    //node information
    ID connectedExternalNodes ;  // node numbers
    Node *nodePointers[3] ;      // pointers to nodes

    //material information
    SectionForceDeformation *materialPointers[nip] ; // pointers to materials
                      
    //local nodal coordinates, two coordinates for each of three nodes
    double xl[2][3] ; 

    //shell basis vectors
    double g1[3] ;
    double g2[3] ;
    double g3[3] ;

    //compute local coordinates and basis
    void computeBasis( ) ;
//start Yuli Huang (yulihuang@gmail.com) & Xinzheng Lu (luxz@tsinghua.edu.cn)
    void updateBasis( ) ;
//end Yuli Huang (yulihuang@gmail.com) & Xinzheng Lu (luxz@tsinghua.edu.cn)
        
    //inertia terms
    void formInertiaTerms( int tangFlag ) ;

    //form residual and tangent                      
    void formResidAndTangent( int tang_flag ) ;

    //compute Bdrill matrix
    //double* computeBdrill( int node, const double shp[3][4] ) ;

    //assemble a B matrix 
    const Matrix& assembleB( const Matrix &Bmembrane,
                             const Matrix &Bbend, 
                             const Matrix &Bshear ) ;
  
    //compute Bmembrane matrix
    const Matrix& computeBmembrane( int node, const double shp[3][3],
                                 const double shpDrill[4][3]) ;
  
    //compute Bbend matrix
    const Matrix& computeBbend( int node, const double shpBend[6][9] ) ;
  
    //Matrix transpose
    //Matrix transpose( int dim1, int dim2, const Matrix &M ) ;

    //shape function routine for four node quads
    void shape2d( double ss, double tt,double qq, 
                  const double x[2][3], 
                  double shp[3][3], 
                  double &xsj ,double sx[2][2]) ;

    //shape function routine for membrane
    void shapeDrill(double ss, double tt,double qq,  
                    const double x[2][3],
                    double sx[2][2], double shpDrill[4][3]);
    //shape function routine for bending
    void shapeBend(double ss, double tt,double qq,  const double x[2][3],
                   double sx[2][2], double shpBend[6][9]);

    // vector for applying loads
    Vector *load;
    Matrix *Ki;

    double b[3];        // Body forces
    double appliedB[3]; // Body forces applied with load
    int applyLoad;

} ; 
