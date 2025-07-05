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
// Three node flat shell element with membrane and drill DOF
//
// Ref: Plate Bending Part - DKT, thin plate element
//      Membrane Part - GT9, a membrane element with drilling DOF
//
// Written: Shuhao Zhang & Xinzheng Lu
// 
// Modified: Yuli Huang (yulihuang@gmail.com) & Xinzheng Lu (luxz@tsinghua.edu.cn)
//
#include <ID.h> 
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <SectionForceDeformation.h>

class ShellNLDKGT : public Element {

 public:
  ShellNLDKGT( );

  ShellNLDKGT(int tag, 
              int node1,
              int node2,
              int node3,
              SectionForceDeformation &theMaterial ) ;

  virtual ~ShellNLDKGT( ) ;

  //set domain because frank is a dumb ass 
  void setDomain( Domain *theDomain ) ;
  
  //get the number of external nodes
  int getNumExternalNodes( ) const final;
  const ID &getExternalNodes( ) final;
  Node **getNodePtrs( ) final;
  int getNumDOF( ) final;

  int commitState() final;
  int revertToLastCommit() final;
  int revertToStart() final;

  // TaggedObject
  void Print( OPS_Stream &s, int flag ) final;
  
  //return stiffness matrix 
  const Vector &getResistingForce( ) final;
  const Vector &getResistingForceIncInertia( ) final;
  const Matrix &getTangentStiff( ) final;
  const Matrix &getInitialStiff( ) final;
  const Matrix &getMass( ) final;

  // methods for applying loads
  void zeroLoad( ) final;    
  int addLoad( ElementalLoad *theLoad, double loadFactor ) final;
  int addInertiaLoadToUnbalance( const Vector &accel ) final;


  // MovableObject
  int sendSelf ( int commitTag, Channel & ) final;
  int recvSelf ( int commitTag, Channel &, FEM_ObjectBroker &) final;


  Response* setResponse( const char **argv, int argc, OPS_Stream & ) final;
  int getResponse( int responseID, Information & ) final;
  int setParameter(const char **argv, int argc, Parameter &param) final;

private : 
  constexpr static int nip = 4;

  // static data
  static Matrix stiff ;
  static Vector resid ;
  static Matrix mass ;
  static Matrix damping ;

  //last resid
  // Vector CstrainGauss,TstrainGauss;
  Vector CstrainGauss,TstrainGauss;  //modify for geometric nonlinearity

  // quadrature data
  static const double three ;
  static const double one_over_three ;
  static const double five;
  static const double one_over_five;
  static const double three_over_five;
  static const double one_over_four;
  static const double wg1;
  static const double wg2;
  static double sg[4] ;
  static double tg[4] ;
  static double qg[4] ;
  static double wg[4] ;

  //node information
  ID connectedExternalNodes ;  // node numbers
  Node *nodePointers[3] ;      // pointers to nodes

  //material information
  SectionForceDeformation *materialPointers[4] ; //pointers to four materials
                    
  //local nodal coordinates, two coordinates for each node
  //static double xl[][4] ; 
  double xl[2][3] ; 

  //shell basis vectors
  double g1[3] ;
  double g2[3] ;
  double g3[3] ;

  // compute local coordinates and basis
  void computeBasis( ) ;

  void updateBasis( ) ;
//end Yuli Huang (yulihuang@gmail.com) & Xinzheng Lu (luxz@tsinghua.edu.cn)

  void formInertiaTerms( int tangFlag ) ;

  //form residual and tangent                      
  void formResidAndTangent( int tang_flag ) ;

  const Matrix& assembleB( const Matrix &Bmembrane,
                           const Matrix &Bbend, 
                           const Matrix &Bshear ) ;

  //compute Bmembrane matrix
  const Matrix& computeBmembrane( int node, const double shp[3][3],
                               const double shpDrill[4][3]) ;

  // compute Bbend matrix
  const Matrix& computeBbend( int node, const double shpBend[6][9] ) ;
  // add for geometric nonlinearity
  const Matrix& computeBG(int node, const double shpBend[6][9]);
  const Vector& computeNLdstrain(const Matrix &BG,const Vector &dispIncLocalBend);

  //shape function routine for four node quads
  void shape2d(double ss, double tt,double qq, 
               const double x[2][3], 
               double shp[3][3], 
               double &xsj ,double sx[2][2]) ;

  // shape function routine for membrane
  void shapeDrill(double ss, double tt,double qq,  
                  const double x[2][3],
                  double sx[2][2], double shpDrill[4][3]);
  // shape function routine for bending
  void shapeBend(double ss, double tt,double qq,  const double x[2][3],
                 double sx[2][2], double shpBend[6][9]);

  // vector for applying loads
  Vector *load;
  Matrix *Ki;
};

