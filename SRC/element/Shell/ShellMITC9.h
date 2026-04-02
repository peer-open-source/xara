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
// 9-node lagrandian shell element with membrane and drill
//
// Written: Leopoldo Tesser, Diego Talledo
//
#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <ID.h> 
#include <Vector.h>
#include <Matrix.h>
#include <MatrixND.h>
#include <Element.h>
#include <Node.h>
#include <SectionForceDeformation.h>
#include <quadrature/GaussLegendre2D.hpp>

class ShellMITC9 : public Element,
                   protected GaussLegendre<2,9>
{
  public:

    // null constructor
    ShellMITC9();

    // full constructor
    ShellMITC9(int tag, 
             int node1,
             int node2,
             int node3,
             int node4,
             int node5,
             int node6,
             int node7,
             int node8,
             int node9,
             SectionForceDeformation &theMaterial ) ;

    virtual ~ShellMITC9();

    const char* getClassType() const final {
      return "ShellMITC9";
    }

    //  get the number of external nodes
    int getNumExternalNodes() const ;

    const ID &getExternalNodes();
    Node **getNodePtrs();
    int getNumDOF();

    void setDomain(Domain *theDomain);
    int commitState();
    int revertToLastCommit();
    int revertToStart();

    const Matrix &getTangentStiff( ) ;
    const Matrix &getInitialStiff();
    const Matrix &getMass();

    // print out element data
    void Print( OPS_Stream &s, int flag ) ;

    // methods for applying loads
    void zeroLoad();        
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    //get residual
    const Vector &getResistingForce( ) ;
    
    //get residual with inertia terms
    const Vector &getResistingForceIncInertia( ) ;

    // public methods for element output
    int sendSelf (int commitTag, Channel &theChannel);
    int recvSelf (int commitTag, Channel &theChannel, FEM_ObjectBroker 
                  &theBroker);


    Response* setResponse(const char **argv, int argc, OPS_Stream &output);
    int getResponse(int responseID, Information &eleInfo);

    int setParameter(const char **argv, int argc, Parameter &param);

  private : 
    static constexpr int NDF = 6;      // two membrane plus three bending plus one drill
    static constexpr int nstress  = 8; // three membrane, three moment, two shear
//  static constexpr int nip   = 9;
    static constexpr int NEN = 9;

    //static data
    static Matrix stiff ;
    static Vector resid ;
    static Matrix mass ;
    static Matrix damping ;

    //quadrature data
    static const double root3 ;
    static const double root3_over_root5 ;
    static double sg[9] ;
    static double tg[9] ;
    static constexpr double wg[9] = {
        25.0 / 81.0,
        40.0 / 81.0,
        25.0 / 81.0,
        40.0 / 81.0,
        25.0 / 81.0,
        40.0 / 81.0,
        25.0 / 81.0,
        40.0 / 81.0,
        64.0 / 81.0
    };

    //node information
    ID connectedExternalNodes ;  // node tags
    Node *theNodes[NEN] ;

    // drilling stiffness
    double Ktt ;

    //material information: pointers to four materials
    SectionForceDeformation *materialPointers[9];
                                          
    // local nodal coordinates, two coordinates for each of nine nodes
    //static double xl[][9] ;
    double xl[2][9] ;

    //shell basis vectors
    double g1[3] ;
    double g2[3] ;
    double g3[3] ;

    // Compute local coordinates and basis
    void computeBasis();
    // State determination
    void formInertiaTerms( int tangFlag ) ;                                    
    void formResidAndTangent(int tang_flag) ;

    //compute Jacobian matrix and inverse at point {L1,L2}
    //void  computeJacobian( double L1, double L2,const double x[2][9], 
    //                       Matrix &JJ,Matrix &JJinv ) ;

    //compute Bdrill matrix
    double* computeBdrill( int node, const double shp[3][9] ) ;

    // assemble B matrix 
    void  assembleB( const Matrix &Bmembrane,
                     const Matrix &Bbend, 
                     const Matrix &Bshear,
                     OpenSees::MatrixND<8,6> &B ) ;
    const Matrix& computeBmembrane( int node, const double shp[3][9]) ;
    
    // compute Bbend matrix
    const Matrix& computeBbend( int node, const double shp[3][9] ) ;
    
    // compute standard Bshear matrix
    const Matrix&  computeBshear( int node, const double shp[3][9] ) ;
    
    //shape function routine for four node quads
    void shape2d( double ss, double tt, 
                  const double x[2][9], 
                  double shp[3][9], 
                  double &xsj ) ;

    // vector for applying loads
    Vector *load;
    Matrix *Ki;
}; 

