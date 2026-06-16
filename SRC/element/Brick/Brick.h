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
                                                                        
// $Revision: 1.15 $
// $Date: 2010-04-23 22:56:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/brick/Brick.h,v $

// Ed "C++" Love
//
// Eight node Brick element 
//
#pragma once

#include <array>
#include <ID.h> 
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <NDMaterial.h>


class Brick : public Element {

  public :
    Brick();
    Brick(int tag, 
          const std::array<int, 8>& node_tags,
          NDMaterial &theMaterial,
          double b1 = 0.0, double b2 = 0.0, double b3 = 0.0);
    
    // destructor 
    virtual ~Brick();

    const char *getClassType() const final {return "Brick";}

    void setDomain( Domain *) final;
    int getNumExternalNodes( ) const final;
    const ID &getExternalNodes() final;
    Node **getNodePtrs() final;
    int getNumDOF() final;


    int commitState();
    int revertToLastCommit();
    int revertToStart();
    int update();

    // return stiffness matrix 
    const Matrix &getTangentStiff();
    const Matrix &getInitialStiff();    
    const Matrix &getMass();    

    void zeroLoad();
    int addLoad(ElementalLoad *, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    const Vector &getResistingForce();
    const Vector &getResistingForceIncInertia();

    // public methods for element output

    Response *setResponse(const char **argv, int argc, OPS_Stream &);
    int getResponse(int responseID, Information &);


    int setParameter(const char **argv, int argc, Parameter &);
    int updateParameter(int parameterID, Information &);

    int sendSelf (int commitTag, Channel &);
    int recvSelf (int commitTag, Channel &, FEM_ObjectBroker  &);

    void Print( OPS_Stream &s, int flag);

  private :

    //
    // private methods
    //

    void formInertiaTerms( int tangFlag ) ;
    void formResidAndTangent( int tang_flag ) ;
    void computeBasis();

    const MatrixND<6,3>&
    computeB(int node, 
             const double shp[4][8],
             MatrixND<6,3> &B
    ) const noexcept;
  

    //
    // private attributes
    //
    constexpr static unsigned int 
                         NEN = 8,  // number of element nodes
                         NDM = 3,  // Spatial dimensions
                         NDF = 3,  // number of element dof
                         NIP = 8,  // number of integration points
                         NST = 6;  // number of stress components

    ID connectedExternalNodes ;  // node tags
    std::array<Node *, 8> theNodes;      // pointers to nodes

    // material information
    NDMaterial *materialPointers[NIP]; // pointers to materials

    double b[3];		// Body forces
    double appliedB[3];		// Body forces applied with load
    int applyLoad;

    Vector *load;
    Matrix *Ki;

    //
    // static attributes
    //

    static Matrix stiff;
    static Vector resid;
    static Matrix mass ;

    // quadrature data
    static constexpr double root3 = 1.73205080757;
    static constexpr double sg[2] = {
      -1.0/root3, 1.0/root3
    };
    static constexpr double wg[NIP] = { 
                              1.0, 1.0, 1.0, 1.0, 
                              1.0, 1.0, 1.0, 1.0  } ;
  
    // local nodal coordinates, three coordinates for each node
    double xl[NDM][NEN]; 

}; 
