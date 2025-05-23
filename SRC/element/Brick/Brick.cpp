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
// Ed "C++" Love
//
// Eight node Brick element
//
#include <stdio.h> 
#include <stdlib.h> 
#include <cmath> 

#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <Domain.h>
#include <ErrorHandler.h>
#include <Brick.h>
#include <shp3d.h>
#include <ElementResponse.h>
#include <Parameter.h>
#include <ElementalLoad.h>
#include <isoparametric.tpp>
#include <Channel.h>
#include <FEM_ObjectBroker.h>


using namespace OpenSees;
//static data
double  Brick::xl[3][8] ;

Matrix  Brick::stiff(24,24) ;
Vector  Brick::resid(24) ;
Matrix  Brick::mass(24,24) ;

    
//quadrature data
const double  Brick::root3 = sqrt(3.0) ;
const double  Brick::one_over_root3 = 1.0 / root3 ;

const double  Brick::sg[] = { -one_over_root3, one_over_root3  } ;

const double  Brick::wg[] = { 1.0, 1.0, 1.0, 1.0, 
                              1.0, 1.0, 1.0, 1.0  } ;

static Matrix B(6,3) ;


Brick::Brick() 
:Element( 0, ELE_TAG_Brick),
 connectedExternalNodes(8), applyLoad(0), load(0), Ki(0)
{
  B.Zero();

  for (int i=0; i<8; i++ ) {
    materialPointers[i] = nullptr;
    theNodes[i] = nullptr;
  }

  b[0] = 0.0;
  b[1] = 0.0;
  b[2] = 0.0;
}


Brick::Brick(int tag, 
            int node1,
            int node2,
            int node3,
            int node4,
            int node5,
            int node6,
            int node7,
            int node8,
            NDMaterial &theMaterial,
            double b1, double b2, double b3)
  : Element(tag, ELE_TAG_Brick),
   connectedExternalNodes(8), applyLoad(0), load(0), Ki(0)
{
  B.Zero();

  connectedExternalNodes(0) = node1 ;
  connectedExternalNodes(1) = node2 ;
  connectedExternalNodes(2) = node3 ;
  connectedExternalNodes(3) = node4 ;

  connectedExternalNodes(4) = node5 ;
  connectedExternalNodes(5) = node6 ;
  connectedExternalNodes(6) = node7 ;
  connectedExternalNodes(7) = node8 ;

  for (int i=0; i<8; i++ ) {
      materialPointers[i] = theMaterial.getCopy("ThreeDimensional");
      theNodes[i] = nullptr;
  }

  // Body forces
  b[0] = b1;
  b[1] = b2;
  b[2] = b3;
}
//******************************************************************


Brick::~Brick()
{

  for (int i=0 ; i<8; i++ ) {
    delete materialPointers[i] ;
  }

  if (load != 0)
    delete load;

  if (Ki != 0)
    delete Ki;
  
}


void
Brick::setDomain( Domain *theDomain ) 
{
  for (int i=0; i<8; i++ ) 
    theNodes[i] = theDomain->getNode( connectedExternalNodes(i) ) ;

  this->DomainComponent::setDomain(theDomain);

}


int
Brick::getNumExternalNodes() const
{
  return 8 ;
} 
 

const ID&
Brick::getExternalNodes() 
{
  return connectedExternalNodes ;
} 


Node **  
Brick::getNodePtrs() 
{
  return &theNodes[0];
} 


int
Brick::getNumDOF() 
{
  return 24 ;
}


int
Brick::commitState( )
{
  int success = 0 ;


  if ((success = this->Element::commitState()) != 0) {
    opserr << "Brick::commitState () - failed in base class";
  }    

  for (int i=0; i<8; i++ ) 
    success += materialPointers[i]->commitState( ) ;
  
  return success ;
}
 


int
Brick::revertToLastCommit( ) 
{
  int success = 0 ;

  for (int i=0; i<8; i++ ) 
    success += materialPointers[i]->revertToLastCommit();
  
  return success ;
}
    


int
Brick::revertToStart( ) 
{
  int success = 0 ;

  for (int i=0; i<8; i++ ) 
    success += materialPointers[i]->revertToStart();
  
  return success ;
}


void
Brick::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_ELEM_INDENT << "{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"Brick\", ";
    s << "\"nodes\": ["
      << connectedExternalNodes(0) << ", ";
    for (int i = 1; i < 7; i++)
        s << connectedExternalNodes(i) << ", ";
    s << connectedExternalNodes(7) << "], ";
    s << "\"bodyForces\": [" << b[0] << ", " << b[1] << ", " << b[2] << "], ";
    s << "\"material\": [" << materialPointers[0]->getTag() << "]}";

    return;
  }

  if (flag == 2) {
    
    s << "#Brick\n";
    
    int i;
    const int numNodes = 8;
    const int nstress = 6;

    for (i = 0; i < numNodes; i++) {
        const Vector &nodeCrd = theNodes[i]->getCrds();
        const Vector &nodeDisp = theNodes[i]->getDisp();
        s << "#NODE " << nodeCrd(0) << " " << nodeCrd(1) << " " << nodeCrd(2)
            << " " << nodeDisp(0) << " " << nodeDisp(1) << " " << nodeDisp(2) << endln;
    }
    
    // spit out the section location & invoke print on the scetion
    const int numMaterials = 8;
    
    static Vector avgStress(nstress);
    static Vector avgStrain(nstress);
    avgStress.Zero();
    avgStrain.Zero();
    for (i = 0; i < numMaterials; i++) {
        avgStress += materialPointers[i]->getStress();
        avgStrain += materialPointers[i]->getStrain();
    }
    avgStress /= numMaterials;
    avgStrain /= numMaterials;
    
    s << "#AVERAGE_STRESS ";
    for (i = 0; i < nstress; i++)
        s << avgStress(i) << " ";
    s << endln;
    
    s << "#AVERAGE_STRAIN ";
    for (i = 0; i < nstress; i++)
        s << avgStrain(i) << " ";
    s << endln;
    
    /*
    for (i=0; i<numMaterials; i++) {
      s << "#MATERIAL\n";
      //      materialPointers[i]->Print(s, flag);
      s << materialPointers[i]->getStress();
    }
    */
  }

  if (flag == OPS_PRINT_CURRENTSTATE) {
      
      s << "Standard Eight Node Brick \n";
      s << "Element Number: " << this->getTag() << endln;
      s << "Nodes: " << connectedExternalNodes;
      
      s << "Material Information : \n ";
      materialPointers[0]->Print(s, flag);
      
      s << endln;
      s << this->getTag() << " " << connectedExternalNodes(0)
          << " " << connectedExternalNodes(1)
          << " " << connectedExternalNodes(2)
          << " " << connectedExternalNodes(3)
          << " " << connectedExternalNodes(4)
          << " " << connectedExternalNodes(5)
          << " " << connectedExternalNodes(6)
          << " " << connectedExternalNodes(7)
          << endln;
      
      s << "Body Forces: " << b[0] << " " << b[1] << " " << b[2] << endln;
      s << "Resisting Force (no inertia): " << this->getResistingForce();
  }
}
 
 
const Matrix&
Brick::getTangentStiff( ) 
{
  int tang_flag = 1 ; //get the tangent 

  formResidAndTangent( tang_flag );

  return stiff;
}

const Matrix&
Brick::getInitialStiff( ) 
{
  if (Ki != 0)
    return *Ki;

  // strains ordered : eps11, eps22, eps33, 2*eps12, 2*eps23, 2*eps31 
  static const int ndm = 3 ;
  static const int ndf = 3 ; 
  static const int nstress = 6 ;
  static const int numberNodes = 8 ;
  static const int numberGauss = 8 ;
  static const int nShape = 4 ;

  int i, j, k, p, q ;
  int jj, kk ;

  static double volume ;
  static double xsj ;  // determinant jacaobian matrix 
  static double dvol[numberGauss] ; //volume element
  static double gaussPoint[ndm] ;
  static Vector strain(nstress) ;  //strain
  static double shp[nShape][numberNodes] ;  //shape functions at a gauss point
  static double Shape[nShape][numberNodes][numberGauss] ; //all the shape functions
  static Matrix stiffJK(ndf,ndf) ; //nodeJK stiffness 
  static Matrix dd(nstress,nstress) ;  //material tangent


  //---------B-matrices------------------------------------

  static Matrix BJ(nstress,ndf) ;      // B matrix node J
  static Matrix BJtran(ndf,nstress) ;
  static Matrix BK(nstress,ndf) ;      // B matrix node k
  static Matrix BJtranD(ndf,nstress) ;

  //-------------------------------------------------------

  
  //zero stiffness and residual 
  stiff.Zero( ) ;

  //compute basis vectors and local nodal coordinates
  computeBasis( ) ;

  // gauss loop to compute and save shape functions 

  int count = 0 ;
  volume = 0.0 ;

  for (int i = 0; i < 2; i++ ) {
    for (int j = 0; j < 2; j++ ) {
      for (int k = 0; k < 2; k++ ) {

        gaussPoint[0] = sg[i] ;        
        gaussPoint[1] = sg[j] ;        
        gaussPoint[2] = sg[k] ;

        //get shape functions    
        shp3d( gaussPoint, xsj, shp, xl ) ;

        //save shape functions
        for (int p = 0; p < nShape; p++ ) {
          for (int q = 0; q < numberNodes; q++ )
            Shape[p][q][count] = shp[p][q] ;
        }

        //volume element to also be saved
        dvol[count] = wg[count] * xsj ;  

        //volume += dvol[count] ;

        count++ ;
      }
    }
  }
  

  // Gauss loop 
  for (int i = 0; i < numberGauss; i++ ) {

    //extract shape functions from saved array
    for (int p = 0; p < nShape; p++ ) {
       for (int q = 0; q < numberNodes; q++ )
        shp[p][q]  = Shape[p][q][i] ;
    }


    dd = materialPointers[i]->getInitialTangent( ) ;
    dd *= dvol[i] ;
    
    jj = 0;
    for (int j = 0; j < numberNodes; j++ ) {

      BJ = computeB( j, shp ) ;
   
      //transpose 
      //BJtran = transpose( nstress, ndf, BJ ) ;
      for (p=0; p<ndf; p++) {
        for (q=0; q<nstress; q++) 
          BJtran(p,q) = BJ(q,p) ;
      }

      //BJtranD = BJtran * dd ;
      BJtranD.addMatrixProduct(0.0,  BJtran, dd, 1.0) ;
      
      kk = 0 ;
      for (int k = 0; k < numberNodes; k++ ) {

        BK = computeB( k, shp ) ;
        
        
        // stiffJK =  BJtranD * BK  ;
        stiffJK.addMatrixProduct(0.0,  BJtranD, BK, 1.0) ;
        
        for ( p = 0; p < ndf; p++ )  {
          for ( q = 0; q < ndf; q++ )
            stiff( jj+p, kk+q ) += stiffJK( p, q ) ;
        }

        kk += ndf;
      }

      jj += ndf;

    }
  }

  Ki = new Matrix(stiff);

  return stiff ;
}    


const Matrix&
Brick::getMass( ) 
{
  int tangFlag = 1 ;

  formInertiaTerms( tangFlag ) ;

  return mass ;
} 



void
Brick::zeroLoad( )
{
  if (load != 0)
    load->Zero();

  applyLoad = 0;

  appliedB[0] = 0.0;
  appliedB[1] = 0.0;
  appliedB[2] = 0.0;

  return ;
}


int 
Brick::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  int type;
  const Vector &data = theLoad->getData(type, loadFactor);

  if (type == LOAD_TAG_BrickSelfWeight) {
      applyLoad = 1;
      appliedB[0] += loadFactor * b[0];
      appliedB[1] += loadFactor * b[1];
      appliedB[2] += loadFactor * b[2];
    return 0;
  } else if (type == LOAD_TAG_SelfWeight) {
      // added compatibility with selfWeight class implemented for all continuum elements, C.McGann, U.W.
      applyLoad = 1;
      appliedB[0] += loadFactor*data(0)*b[0];
      appliedB[1] += loadFactor*data(1)*b[1];
      appliedB[2] += loadFactor*data(2)*b[2];
      return 0;
  } else {
    opserr << "Brick::addLoad() - ele with tag: " << this->getTag() << " does not deal with load type: " << type << "\n";
    return -1;
  }

  return -1;
}

int
Brick::addInertiaLoadToUnbalance(const Vector &accel)
{
  static const int numberNodes = 8 ;
  static const int numberGauss = 8 ;
  static const int ndf = 3 ; 

  int i;

  // check to see if have mass
  int haveRho = 0;
  for (i = 0; i < numberGauss; i++) {
    if (materialPointers[i]->getRho() != 0.0)
      haveRho = 1;
  }

  if (haveRho == 0)
    return 0;

  // Compute mass matrix
  int tangFlag = 1 ;
  formInertiaTerms( tangFlag ) ;

  // store computed RV for nodes in resid vector
  int count = 0;
  for (i=0; i<numberNodes; i++) {
    const Vector &Raccel = theNodes[i]->getRV(accel);
    for (int j=0; j<ndf; j++)
      resid(count++) = Raccel(j);
  }

  // create the load vector if one does not exist
  if (load == 0) 
    load = new Vector(numberNodes*ndf);

  // add -M * RV(accel) to the load vector
  load->addMatrixVector(1.0, mass, resid, -1.0);
  
  return 0;
}


const Vector&
Brick::getResistingForce( ) 
{
  int tang_flag = 0 ; // don't get the tangent

  formResidAndTangent( tang_flag ) ;

  if (load != 0)
    resid -= *load;

  return resid ;   
}


//get residual with inertia terms
const Vector&  Brick::getResistingForceIncInertia( )
{
  static Vector res(24);

  int tang_flag = 0 ; //don't get the tangent

  //do tangent and residual here 
  formResidAndTangent( tang_flag ) ;

  formInertiaTerms( tang_flag ) ;

  res = resid;

  // add the damping forces if rayleigh damping
  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
      res += this->getRayleighDampingForces();

  if (load != 0)
    res -= *load;

  return res;
}


//*********************************************************************

void
Brick::formInertiaTerms( int tangFlag ) 
{

  static constexpr int ndm = 3 ;
  static constexpr int ndf = 3 ; 
  static constexpr int numberNodes = 8 ;
  static constexpr int numberGauss = 8 ;
  static constexpr int nShape = 4 ;
  static constexpr int massIndex = nShape - 1 ;


  double dvol[numberGauss] ; //volume element
  static double Shape[nShape][numberNodes][numberGauss] ; //all the shape functions
  static double gaussPoint[ndm] ;

  static Vector momentum(ndf) ;

  int i, j, k, p, q ;
  int jj, kk ;

  double temp, rho, massJK ;


  mass.Zero();

  // compute basis vectors and local nodal coordinates
  computeBasis( ) ;

  //gauss loop to compute and save shape functions 

  int count = 0 ;

  for (int i = 0; i < 2; i++ ) {
    for (int j = 0; j < 2; j++ ) {
      for (int k = 0; k < 2; k++ ) {

        gaussPoint[0] = sg[i] ;        
        gaussPoint[1] = sg[j] ;        
        gaussPoint[2] = sg[k] ;

        // get shape functions
        double xsj;
        double shp[nShape][numberNodes];
        shp3d( gaussPoint, xsj, shp, xl ) ;

        // save shape functions
        for (int p = 0; p < nShape; p++ ) {
          for (int q = 0; q < numberNodes; q++ )
            Shape[p][q][count] = shp[p][q] ;
        }

        // volume element to also be saved
        dvol[count] = wg[count] * xsj ;  

        count++ ;

      }
    }
  } 
  


  // Gauss loop 
  for ( i = 0; i < numberGauss; i++ ) {

    // extract shape functions from saved array
    double shp[nShape][numberNodes];
    for ( p = 0; p < nShape; p++ ) {
       for ( q = 0; q < numberNodes; q++ )
        shp[p][q]  = Shape[p][q][i] ;
    }


    // node loop to compute acceleration
    momentum.Zero( ) ;
    for ( j = 0; j < numberNodes; j++ ) 
      //momentum += shp[massIndex][j] * ( theNodes[j]->getTrialAccel()  ) ; 
      momentum.addVector(1.0, theNodes[j]->getTrialAccel(), shp[massIndex][j]);


    // density
    rho = materialPointers[i]->getRho() ;


    // multiply acceleration by density to form momentum
    momentum *= rho;


    //residual and tangent calculations node loops
    jj = 0 ;
    for (int j = 0; j < numberNodes; j++ ) {

      temp = shp[massIndex][j] * dvol[i] ;

      for (int p = 0; p < ndf; p++ )
        resid( jj+p ) += ( temp * momentum(p) )  ;


      if ( tangFlag == 1 ) {

        temp *= rho ;

        // node-node mass
        kk = 0 ;
        for (int k = 0; k < numberNodes; k++ ) {

          massJK = temp * shp[massIndex][k] ;

          for (int p = 0; p < ndf; p++ )  
            mass( jj+p, kk+p ) += massJK ;
                
          kk += ndf ;
        }

      } // end if tang_flag 

      jj += ndf ;
    }
  }
}

//*********************************************************************
//form residual and tangent
int  
Brick::update() 
{

  //strains ordered : eps11, eps22, eps33, 2*eps12, 2*eps23, 2*eps31 

  static constexpr int ndm = 3 ;
  static constexpr int ndf = 3 ; 
  static constexpr int nstress = 6 ;
  static constexpr int numberNodes = 8 ;
  static constexpr int numberGauss = 8 ;
  static constexpr int nShape = 4 ;

  
  double dvol[numberGauss] ; //volume element
  static Vector strain(nstress);
  double Shape[nShape][numberNodes][numberGauss] ; //all the shape functions

  //---------B-matrices------------------------------------

  static Matrix BJ(nstress,ndf) ;      // B matrix node J
  static Matrix BJtran(ndf,nstress) ;
  static Matrix BK(nstress,ndf) ;      // B matrix node k
  static Matrix BJtranD(ndf,nstress) ;

  //-------------------------------------------------------

  
  //compute basis vectors and local nodal coordinates
  computeBasis();

  // gauss loop to compute and save shape functions 

  int count = 0 ;
  double volume = 0.0 ;

  for (int i = 0; i < 2; i++ ) {
    for (int j = 0; j < 2; j++ ) {
      for (int k = 0; k < 2; k++ ) {

        double gaussPoint[NDM] {sg[i], sg[j], sg[k]};

        //get shape functions
        double xsj;
        double shp[nShape][numberNodes];
        shp3d( gaussPoint, xsj, shp, xl ) ;

        //save shape functions
        for (int p = 0; p < nShape; p++ ) {
          for (int q = 0; q < numberNodes; q++ )
            Shape[p][q][count] = shp[p][q] ;
        }

        //volume element to also be saved
        dvol[count] = wg[count] * xsj ;  

        count++ ;

      }
    }
  }

  int success = -1;

  // Gauss loop 
  for (int i = 0; i < numberGauss; i++ ) {

    // extract shape functions from saved array
    double shp[nShape][numberNodes];
    for (int p = 0; p < nShape; p++ ) {
      for (int q = 0; q < numberNodes; q++ )
        shp[p][q]  = Shape[p][q][i] ;
    }


    //zero the strains
    strain.Zero( ) ;

    // j-node loop to compute strain 
    for (int j = 0; j < numberNodes; j++ )  {

      /**************** fmk - unwinding for performance
      //compute B matrix 
      BJ = computeB( j, shp ) ;

      //nodal displacements 
      const Vector &ul = theNodes[j]->getTrialDisp( ) ;

      // compute the strain
      // strain += (BJ*ul) ;
      ***************************************************/

      double b00 = shp[0][j];
      double b11 = shp[1][j];
      double b22 = shp[2][j];
      double b30 = shp[1][j];
      double b31 = shp[0][j];
      double b41 = shp[2][j];
      double b42 = shp[1][j];
      double b50 = shp[2][j];
      double b52 = shp[0][j];

      const Vector &ul = theNodes[j]->getTrialDisp();

      double ul0 = ul(0);
      double ul1 = ul(1);
      double ul2 = ul(2);

      strain(0) += b00 * ul0;
      strain(1) += b11 * ul1;
      strain(2) += b22 * ul2;
      strain(3) += b30 * ul0 + b31 * ul1;
      strain(4) += b41 * ul1 + b42 * ul2;
      strain(5) += b50 * ul0 + b52 * ul2;
    }
    
    //
    success = materialPointers[i]->setTrialStrain(strain);

  }

  return 0;
}


void
Brick::formResidAndTangent( int tang_flag ) 
{

  //strains ordered : eps11, eps22, eps33, 2*eps12, 2*eps23, 2*eps31 

  static constexpr int ndm = 3 ;
  static constexpr int ndf = 3 ; 
  static constexpr int nstress = 6 ;
  static constexpr int numberNodes = 8 ;
  static constexpr int numberGauss = 8 ;
  static constexpr int nShape = 4 ;



  static double dvol[numberGauss] ; //volume element
  static double shp[nShape][numberNodes] ;  //shape functions at a gauss point
  static double Shape[nShape][numberNodes][numberGauss] ; //all the shape functions

  static Vector residJ(ndf) ; //nodeJ residual 
  static Matrix stiffJK(ndf,ndf) ; //nodeJK stiffness 
  static Vector stress(nstress) ;  //stress
  static Matrix dd(nstress,nstress) ;  //material tangent


  //---------B-matrices------------------------------------

  static Matrix BJ(nstress,ndf) ;      // B matrix node J
  static Matrix BJtran(ndf,nstress) ;
  static Matrix BK(nstress,ndf) ;      // B matrix node k
  static Matrix BJtranD(ndf,nstress) ;

  //-------------------------------------------------------

  
  //zero stiffness and residual 
  stiff.Zero( ) ;
  resid.Zero( ) ;

  //compute basis vectors and local nodal coordinates
  computeBasis();

  // gauss loop to compute and save shape functions 

  int count = 0 ;

  int i, j, k, p, q ;

  for (int i = 0; i < 2; i++ ) {
    for (int j = 0; j < 2; j++ ) {
      for (int k = 0; k < 2; k++ ) {


        // Evaluate shape functions
        double xsj ;  // determinant jacaobian matrix

        const double xi[ndm]  = {sg[i], sg[j], sg[k]} ;
        shp3d( xi, xsj, shp, xl ) ;

        // save shape functions
        for ( p = 0; p < nShape; p++ ) {
          for ( q = 0; q < numberNodes; q++ )
            Shape[p][q][count] = shp[p][q] ;
        }

        //volume element to also be saved
        dvol[count] = wg[count] * xsj ;  

        count++ ;
      }
    }
  }
  
  //
  // gauss loop
  //
  for (int i = 0; i < numberGauss; i++ ) {

    //extract shape functions from saved array
    for (int p = 0; p < nShape; p++ ) {
      for (int q = 0; q < numberNodes; q++ )
        shp[p][q]  = Shape[p][q][i];
    }

    // compute the stress
    stress = materialPointers[i]->getStress( ) ;

    stress  *= dvol[i] ;
  
  
    if ( tang_flag == 1 ) {
      dd = materialPointers[i]->getTangent( ) ;
      dd *= dvol[i];
    }


    double stress0 = stress(0);
    double stress1 = stress(1);
    double stress2 = stress(2);
    double stress3 = stress(3);
    double stress4 = stress(4);
    double stress5 = stress(5);

    //residual and tangent calculations node loops

    int jj = 0 ;
    for (int j = 0; j < numberNodes; j++ ) {

      /* ************** fmk - unwinding for performance 
      ************************************************* */

      //               | N,1      0     0    | 
      //   B       =   |   0     N,2    0    |
      //               |   0      0     N,3  |   (6x3)
      //               | N,2     N,1     0   |
      //               |   0     N,3    N,2  |
      //               | N,3      0     N,1  |

      double b00 = shp[0][j];
      double b11 = shp[1][j];
      double b22 = shp[2][j];
      double b30 = shp[1][j];
      double b31 = shp[0][j];
      double b41 = shp[2][j];
      double b42 = shp[1][j];
      double b50 = shp[2][j];
      double b52 = shp[0][j];

      residJ(0) = b00 * stress0 + b30 * stress3 + b50 * stress5;
      residJ(1) = b11 * stress1 + b31 * stress3 + b41 * stress4;
      residJ(2) = b22 * stress2 + b42 * stress4 + b52 * stress5;
      
      BJ = computeB( j, shp ) ;
   
      //transpose
      for (int p=0; p<ndf; p++) {
        for (int q=0; q<nstress; q++) 
          BJtran(p,q) = BJ(q,p) ;
      }

      // residual 
      for (int p = 0; p < ndf; p++ ) {
        resid( jj + p ) += residJ(p)  ;
        if (applyLoad == 0)
          resid( jj + p ) -= dvol[i]*b[p]*shp[3][j];
        else
          resid( jj + p ) -= dvol[i]*appliedB[p]*shp[3][j];
      }

      if ( tang_flag == 1 ) {

        //BJtranD = BJtran * dd ;
        BJtranD.addMatrixProduct(0.0,  BJtran,dd,1.0) ;

        int kk = 0 ;
        for (int k = 0; k < numberNodes; k++ ) {

            BK = computeB( k, shp ) ;
            //stiffJK =  BJtranD * BK  ;
            stiffJK.addMatrixProduct(0.0,  BJtranD, BK,1.0) ;

            for (int p = 0; p < ndf; p++ )  {
              for (int q = 0; q < ndf; q++ )
                  stiff( jj+p, kk+q ) += stiffJK( p, q ) ;
            }

            kk += ndf ;
          } // end for k loop

      } // end if tang_flag 

      jj += ndf ;
    } // end for j loop


  } //end for i gauss loop 

  return ;
}


//************************************************************************
//compute local coordinates and basis

void
Brick::computeBasis( ) 
{

  //nodal coordinates 

  for (int i = 0; i < 8; i++ ) {
    const Vector &coorI = theNodes[i]->getCrds( ) ;

    xl[0][i] = coorI(0) ;
    xl[1][i] = coorI(1) ;
    xl[2][i] = coorI(2) ;

  }
}

//*************************************************************************
//compute B

const Matrix&   
Brick::computeB( int node, const double shp[4][8] )
{

//---B Matrix in standard {1,2,3} mechanics notation---------
//
//                -                   -
//               | N,1      0     0    | 
//   B       =   |   0     N,2    0    |
//               |   0      0     N,3  |   (6x3)
//               | N,2     N,1     0   |
//               |   0     N,3    N,2  |
//               | N,3      0     N,1  |
//                -                   -       
//
//-------------------------------------------------------------------
  B.Zero();
  B(0,0) = shp[0][node] ;
  B(1,1) = shp[1][node] ;
  B(2,2) = shp[2][node] ;

  B(3,0) = shp[1][node] ;
  B(3,1) = shp[0][node] ;

  B(4,1) = shp[2][node] ;
  B(4,2) = shp[1][node] ;

  B(5,0) = shp[2][node] ;
  B(5,2) = shp[0][node] ;

  return B ;

}

//**********************************************************************

int
Brick::sendSelf (int commitTag, Channel &theChannel)
{
  int res = 0;
  
  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();
  
  // Quad packs its data into a Vector and sends this to theChannel
  // along with its dbTag and the commitTag passed in the arguments

  // Now quad sends the ids of its materials
  int matDbTag;
  
  static ID idData(26);

  idData(24) = this->getTag();
  if (alphaM != 0 || betaK != 0 || betaK0 != 0 || betaKc != 0) 
    idData(25) = 1;
  else
    idData(25) = 0;
  
  int i;
  for (i = 0; i < 8; i++) {
    idData(i) = materialPointers[i]->getClassTag();
    matDbTag = materialPointers[i]->getDbTag();
    // NOTE: we do have to ensure that the material has a database
    // tag if we are sending to a database channel.
    if (matDbTag == 0) {
      matDbTag = theChannel.getDbTag();
      if (matDbTag != 0)
	materialPointers[i]->setDbTag(matDbTag);
    }
    idData(i+8) = matDbTag;
  }
  
  idData(16) = connectedExternalNodes(0);
  idData(17) = connectedExternalNodes(1);
  idData(18) = connectedExternalNodes(2);
  idData(19) = connectedExternalNodes(3);
  idData(20) = connectedExternalNodes(4);
  idData(21) = connectedExternalNodes(5);
  idData(22) = connectedExternalNodes(6);
  idData(23) = connectedExternalNodes(7);

  res += theChannel.sendID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING Brick::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return res;
  }

  static Vector dData(7);
  dData(0) = alphaM;
  dData(1) = betaK;
  dData(2) = betaK0;
  dData(3) = betaKc;
  dData(4) = b[0];
  dData(5) = b[1];
  dData(6) = b[2];

  if (theChannel.sendVector(dataTag, commitTag, dData) < 0) {
    opserr << "Brick::sendSelf() - failed to send double data\n";
    return -1;
  }    

  // Finally, quad asks its material objects to send themselves
  for (i = 0; i < 8; i++) {
    res += materialPointers[i]->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING Brick::sendSelf() - " << this->getTag() << " failed to send its Material\n";
      return res;
    }
  }
  
  return res;

}
    
int  Brick::recvSelf (int commitTag, 
		      Channel &theChannel, 
		      FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  int dataTag = this->getDbTag();

  static ID idData(26);
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING Brick::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return res;
  }

  this->setTag(idData(24));

  static Vector dData(7);
  if (theChannel.recvVector(dataTag, commitTag, dData) < 0) {
    opserr << "DispBeamColumn2d::sendSelf() - failed to recv double data\n";
    return -1;
  }    
  alphaM = dData(0);
  betaK = dData(1);
  betaK0 = dData(2);
  betaKc = dData(3);
  b[0] = dData(4);
  b[1] = dData(5);
  b[2] = dData(6);


  connectedExternalNodes(0) = idData(16);
  connectedExternalNodes(1) = idData(17);
  connectedExternalNodes(2) = idData(18);
  connectedExternalNodes(3) = idData(19);
  connectedExternalNodes(4) = idData(20);
  connectedExternalNodes(5) = idData(21);
  connectedExternalNodes(6) = idData(22);
  connectedExternalNodes(7) = idData(23);


  if (materialPointers[0] == 0) {
    for (int i = 0; i < 8; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+8);
      // Allocate new material with the sent class tag
      materialPointers[i] = theBroker.getNewNDMaterial(matClassTag);
      if (materialPointers[i] == 0) {
	opserr << "Brick::recvSelf() - Broker could not create NDMaterial of class type " << matClassTag << endln;
	return -1;
      }
      // Now receive materials into the newly allocated space
      materialPointers[i]->setDbTag(matDbTag);
      res += materialPointers[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
        opserr << "NLBeamColumn3d::recvSelf() - material " << i << "failed to recv itself\n";
        return res;
      }
    }
  }
  // materials exist , ensure materials of correct type and recvSelf on them
  else {
    for (int i = 0; i < 8; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+8);
      // Check that material is of the right type; if not,
      // delete it and create a new one of the right type
      if (materialPointers[i]->getClassTag() != matClassTag) {
        delete materialPointers[i];
        materialPointers[i] = theBroker.getNewNDMaterial(matClassTag);
        if (materialPointers[i] == 0) {
          opserr << "Brick::recvSelf() - Broker could not create NDMaterial of class type " <<
          matClassTag << endln;
          return -1;
        }
        materialPointers[i]->setDbTag(matDbTag);
      }
      // Receive the material

      res += materialPointers[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
        opserr << "Brick::recvSelf() - material " << i << "failed to recv itself\n";
        return res;
      }
    }
  }

  return res;
}
//**************************************************************************


Response*
Brick::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  Response *theResponse = 0;

  char outputData[32];

  output.tag("ElementOutput");
  output.attr("eleType","Brick");
  output.attr("eleTag",this->getTag());
  for (int i=1; i<=8; i++) {
    sprintf(outputData,"node%d",i);
    output.attr(outputData, theNodes[i-1]->getTag());
  }


  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0) {

    for (int i=1; i<=8; i++) {
      sprintf(outputData,"P1_%d",i);
      output.tag("ResponseType",outputData);
      sprintf(outputData,"P2_%d",i);
      output.tag("ResponseType",outputData);
      sprintf(outputData,"P3_%d",i);
      output.tag("ResponseType",outputData);
    }

    theResponse = new ElementResponse(this, 1, resid);
  
  }   else if (strcmp(argv[0],"material") == 0 || strcmp(argv[0],"integrPoint") == 0) {

    int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= 8) {

      output.tag("GaussPoint");
      output.attr("number",pointNum);

      theResponse =  materialPointers[pointNum-1]->setResponse(&argv[2], argc-2, output);

      output.endTag(); // GaussPoint
    }


  }
  
  else if (strcmp(argv[0],"stresses") ==0) {

    for (int i=0; i<8; i++) {
      output.tag("GaussPoint");
      output.attr("number",i+1);
      output.tag("NdMaterialOutput");
      output.attr("classType", materialPointers[i]->getClassTag());
      output.attr("tag",       materialPointers[i]->getTag());

      output.tag("ResponseType","sigma11");
      output.tag("ResponseType","sigma22");
      output.tag("ResponseType","sigma33");
      output.tag("ResponseType","sigma12");
      output.tag("ResponseType","sigma23");
      output.tag("ResponseType","sigma13");      

      output.endTag(); // NdMaterialOutput
      output.endTag(); // GaussPoint
    }
    theResponse =  new ElementResponse(this, 3, Vector(48));

  }
  
  else if (strcmp(argv[0],"strains") ==0) {

    for (int i=0; i<8; i++) {
      output.tag("GaussPoint");
      output.attr("number",i+1);
      output.tag("NdMaterialOutput");
      output.attr("classType", materialPointers[i]->getClassTag());
      output.attr("tag", materialPointers[i]->getTag());

      output.tag("ResponseType","eps11");
      output.tag("ResponseType","eps22");
      output.tag("ResponseType","eps33");
      output.tag("ResponseType","eps12");
      output.tag("ResponseType","eps23");
      output.tag("ResponseType","eps13");      

      output.endTag(); // NdMaterialOutput
      output.endTag(); // GaussPoint
    }
    theResponse =  new ElementResponse(this, 4, Vector(48));
  }

  else if (strcmp(argv[0], "stressAtNodes") == 0) {
    theResponse = new ElementResponse(this, 11, Vector(NST*NEN));
  }

  else if (strcmp(argv[0], "shape") == 0) {
    output.tag("Shape");
    output.attr("number", 1);
    output.attr("type", "Brick");
    output.attr("tag", this->getTag());
    theResponse = new ElementResponse(this, 500 + atoi(argv[1])*10 + atoi(argv[2]), Vector(8));
  }

  output.endTag(); // ElementOutput
  return theResponse;
}

int 
Brick::getResponse(int responseID, Information &eleInfo)
{
  static Vector stresses(48);

  if (responseID == 1)
    return eleInfo.setVector(this->getResistingForce());

  else if (responseID == 2)
    return eleInfo.setMatrix(this->getTangentStiff());
    
  else if (responseID == 3) {
    
    // Loop over the integration points
    int cnt = 0;
    for (int i = 0; i < 8; i++) {
      
      // Get material stress response
      const Vector &sigma = materialPointers[i]->getStress();
      stresses(cnt++) = sigma(0);
      stresses(cnt++) = sigma(1);
      stresses(cnt++) = sigma(2);
      stresses(cnt++) = sigma(3);
      stresses(cnt++) = sigma(4);
      stresses(cnt++) = sigma(5);
    }
    return eleInfo.setVector(stresses);

  } 
  else if (responseID == 4) {
    
    // Loop over the integration points
    int cnt = 0;
    for (int i = 0; i < 8; i++) {
      
      // Get material stress response
      const Vector &sigma = materialPointers[i]->getStrain();
      stresses(cnt++) = sigma(0);
      stresses(cnt++) = sigma(1);
      stresses(cnt++) = sigma(2);
      stresses(cnt++) = sigma(3);
      stresses(cnt++) = sigma(4);
      stresses(cnt++) = sigma(5);
    }
    return eleInfo.setVector(stresses);
  }

  else if (responseID == 11) {
    constexpr OpenSees::MatrixND<8,8> We {{
      { 2.549038105676660, -0.683012701892219,  0.183012701892219, -0.683012701892219, -0.683012701892219, 0.183012701892219, -0.049038105676658, 0.183012701892219}, 
      {-0.683012701892219,  2.549038105676660, -0.683012701892219,  0.183012701892219, 0.183012701892219, -0.683012701892219, 0.183012701892219, -0.049038105676658}, 
      { 0.183012701892219, -0.683012701892219,  2.549038105676660, -0.683012701892219, -0.049038105676658, 0.183012701892219, -0.683012701892219, 0.183012701892219}, 
      {-0.683012701892219,  0.183012701892219, -0.683012701892219,  2.549038105676660, 0.183012701892219, -0.049038105676658, 0.183012701892219, -0.683012701892219}, 
      {-0.683012701892219,  0.183012701892219, -0.049038105676658,  0.183012701892219, 2.54903810567666, -0.683012701892219, 0.183012701892219, -0.683012701892219}, 
      { 0.183012701892219, -0.683012701892219,  0.183012701892219, -0.049038105676658, -0.683012701892219, 2.54903810567666, -0.683012701892219, 0.183012701892219}, 
      {-0.049038105676658,  0.183012701892219, -0.683012701892219,  0.183012701892219, 0.183012701892219, -0.683012701892219, 2.54903810567666, -0.683012701892219}, 
      { 0.183012701892219, -0.049038105676658,  0.183012701892219, -0.683012701892219, -0.683012701892219, 0.183012701892219, -0.683012701892219, 2.54903810567666}
    }};

    static VectorND<NST*NEN> stressAtNodes;
    static Vector output(stressAtNodes);

    stressAtNodes.zero();
    OpenSees::StressExtrapolation<NEN,NIP,NST>(materialPointers, We, stressAtNodes);
    return eleInfo.setVector(output);
  }


  else
    return -1;
}

int
Brick::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  int res = -1;

  if ((strstr(argv[0],"material") != 0) && (strcmp(argv[0],"materialState") != 0)) {

    if (argc < 3)
      return -1;

    int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= 8)
      return materialPointers[pointNum-1]->setParameter(&argv[2], argc-2, param);
    else 
      return -1;
  }
  
  // otherwise it could be just a forall material parameter
  else {
    int matRes = res;
    for (int i=0; i<8; i++) {
      matRes =  materialPointers[i]->setParameter(argv, argc, param);
      if (matRes != -1)
        res = matRes;
    }
  }
  
  return res;
}
    
int
Brick::updateParameter(int parameterID, Information &info)
{
    int res = -1;
	int matRes = res;

    if (parameterID == res) {
        return -1;
    } else {
        for (int i = 0; i<8; i++) {
            matRes = materialPointers[i]->updateParameter(parameterID, info);
        }
		if (matRes != -1) {
			res = matRes;
		}
		return res;
    }
}

