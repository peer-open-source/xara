//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Q1/E4 - Enhanced Four Node Quadrilateral Element
//
#include <stdio.h> 
#include <string.h>
#include <stdlib.h> 
#include <math.h> 

#include <ID.h> 
#include <Vector.h>
#include <Matrix.h>
#include <MatrixND.h>
#include <Element.h>
#include <Node.h>
#include <NDMaterial.h>
#include <Domain.h>
#include <ErrorHandler.h>
#include <EnhancedQuad.h>
#include <ElementResponse.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

using namespace OpenSees;


// static data
double  EnhancedQuad::xl[2][4] ;
Matrix  EnhancedQuad::stiff(8,8) ;
Vector  EnhancedQuad::resid(8) ;
Matrix  EnhancedQuad::mass(8,8) ;



// compute Jacobian matrix and inverse at point {L1,L2} 
template <typename MatT>
static inline void 
computeJacobian(double L1, double L2,
                const double x[2][4], 
                MatT &JJ, 
                MatT &JJinv )
{
     
  static constexpr double s[] = { -0.5,  0.5, 0.5, -0.5 } ;
  static constexpr double t[] = { -0.5, -0.5, 0.5,  0.5 } ;

  double shp[2][4] ;

  double ss = L1;
  double tt = L2;

  for (int i = 0; i < 4; i++ ) {
      shp[0][i] = s[i] * ( 0.5 + t[i]*tt ) ;
      shp[1][i] = t[i] * ( 0.5 + s[i]*ss ) ;
  }

  // Construct jacobian and its inverse

  for (int i = 0; i < 2; i++ ) {
    for (int j = 0; j < 2; j++ ) {
      JJ(i,j) = 0.0;
      for (int k = 0; k < 4; k++ )
          JJ(i,j) +=  x[i][k] * shp[j][k] ;

    }
  }

  double xsj = JJ(0,0)*JJ(1,1) - JJ(0,1)*JJ(1,0) ;

  // inverse jacobian
  double jinv = 1.0 / xsj ;
  JJinv(0,0) =  JJ(1,1) * jinv ;
  JJinv(1,1) =  JJ(0,0) * jinv ;
  JJinv(0,1) = -JJ(0,1) * jinv ;
  JJinv(1,0) = -JJ(1,0) * jinv ;

  return;
}

void 
EnhancedQuad::computeB( int node, const double shp[3][4] , MatrixND<3,2> & B)
{
// ---B Matrix in standard {1,2,3} mechanics notation---------------
//
//                -             -
//               | +N,1      0   | 
// B    =        |   0     +N,2  |    (3x2)
//               | +N,2    +N,1  |
//                -             -  
//
// -------------------------------------------------------------------

// static Matrix B(3,2) ;
  B.zero( ) ;

  B(0,0) = shp[0][node];
  B(1,1) = shp[1][node];
  B(2,0) = shp[1][node];
  B(2,1) = shp[0][node];

  return;

}

//***********************************************************************
// compute enhanced strain B-matrices

template <typename MatT, typename BT>
void  
computeBenhanced(int node, 
                  double L1,
                  double L2,
                  double j, 
                  const MatT &Jinv,
                  BT &B)
{
  // static Matrix B(3,2) ;
  static double JinvTran[2][2] ;
  static double shape[2] ;


  // compute JinvTran
  JinvTran[0][0] = Jinv(0,0) ;
  JinvTran[1][1] = Jinv(1,1) ;
  JinvTran[0][1] = Jinv(1,0) ;
  JinvTran[1][0] = Jinv(0,1) ;      // residual 

  double parameter ;
  if ( node == 0 ) {
    // first column of JinvTran 
    shape[0] = JinvTran[0][0] ;
    shape[1] = JinvTran[1][0] ;

    parameter = L1 / j ;

  }
  else if ( node == 1 ) {
    // second column of JinvTran 
    shape[0] = JinvTran[0][1] ;
    shape[1] = JinvTran[1][1] ;

    parameter = L2 / j ;
  } 

  shape[0] *= parameter ;
  shape[1] *= parameter ; 

  B(0,0) = shape[0] ;
  B(1,0) = 0.0;
  B(2,0) = shape[1] ;
  B(0,1) = 0.0;
  B(1,1) = shape[1] ;
  B(2,1) = shape[0] ;
  
  return;
}

EnhancedQuad::EnhancedQuad() :
Element( 0, ELE_TAG_EnhancedQuad ),
connectedExternalNodes(4),
alpha(4), thickness(0.0), load(0), Ki(0)
{ 
  for ( int i = 0 ;  i < 4; i++ )
    materialPointers[i] = nullptr;

  // zero enhanced parameters
  alpha.Zero( ) ;
}

// full constructor
EnhancedQuad::EnhancedQuad(int tag, 
                           std::array<int,4>& nodes,
                           NDMaterial &theMaterial,
                           double thickness) 
 :
  Element( tag, ELE_TAG_EnhancedQuad ),
  connectedExternalNodes(4),
  alpha(4), thickness(thickness), load(0), Ki(0)
{

  for (int i=0; i<NEN; i++) {
    connectedExternalNodes(i) = nodes[i];
    theNodes[i] = nullptr;
  }


  for (int i = 0 ;  i < nip; i++ )
      materialPointers[i] = theMaterial.getCopy() ;

  // zero enhanced parameters
  alpha.Zero( ) ;

}

// destructor 
EnhancedQuad::~EnhancedQuad()
{
  for (int i = 0 ;  i < nip; i++ ) 
    if (materialPointers[i] != nullptr)
      delete materialPointers[i] ;

  if (load != nullptr)
    delete load;

  if (Ki != nullptr)
    delete Ki;
}


int
EnhancedQuad::getNumExternalNodes() const
{
  return NEN;
} 


const ID&
EnhancedQuad::getExternalNodes( ) 
{
  return connectedExternalNodes ;
} 

Node **
EnhancedQuad::getNodePtrs() 
{
  return theNodes;
} 


int
EnhancedQuad::getNumDOF() 
{
  return 8 ;
}



void
EnhancedQuad::setDomain( Domain *theDomain ) 
{
  if (theDomain == nullptr) {
    for (int i=0; i<NEN; i++)
      theNodes[i] = nullptr;
    return;
  }

  // node pointers
  for (int i = 0; i < NEN; i++ ) 
    theNodes[i] = theDomain->getNode( connectedExternalNodes(i) ) ;

  this->DomainComponent::setDomain(theDomain) ;
}


int
EnhancedQuad::commitState()
{
  int success = 0 ;

  // call element commitState to do any base class stuff
  if ((success = this->Element::commitState()) != 0) {
    opserr << "EnhancedQuad::commitState () - failed in base class";
  }    

  for (int i = 0; i < nip; i++ ) 
    success += materialPointers[i]->commitState( ) ;
  
  return success ;
}


int
EnhancedQuad::revertToLastCommit( ) 
{
  int success = 0 ;

  for (int i = 0; i < nip; i++ ) 
    success += materialPointers[i]->revertToLastCommit();
  
  return success ;
}


int
EnhancedQuad::revertToStart() 
{
  int success = 0 ;

  // zero enhanced parameters
  this->alpha.Zero() ;

  for (int i = 0; i < nip; i++ ) 
    success += materialPointers[i]->revertToStart();
  
  return success ;
}


// return stiffness matrix 
const Matrix&  EnhancedQuad::getTangentStiff() 
{
  int tang_flag = 1 ; // get the tangent 

  // do tangent and residual here
  formResidAndTangent( tang_flag ) ;  

  return stiff ;
}    

// return secant matrix 
const Matrix&  EnhancedQuad::getInitialStiff( ) 
{

  if (Ki != 0)
    return *Ki;

  int tang_flag = 1 ; // get the tangent

  return stiff;
}


const Matrix&
EnhancedQuad::getMass() 
{

  int tangFlag = 1 ;

  formInertiaTerms( tangFlag ) ;

  return mass ;

} 

void
EnhancedQuad::zeroLoad()
{
  if (load != nullptr)
    load->Zero();

  return ;
}

int 
EnhancedQuad::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  opserr << "EnhancedQuad::addLoad - load type unknown for ele with tag: " << this->getTag() << endln;
  return -1;
}

int
EnhancedQuad::addInertiaLoadToUnbalance(const Vector &accel)
{

  // check to see if have mass
  int haveRho = 0;
  for (int i = 0; i < nip; i++) {
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
  for (int i=0; i<numberNodes; i++) {
    const Vector &Raccel = theNodes[i]->getRV(accel);
    for (int j=0; j<NDF; j++)
      resid(count++) = Raccel(i);
  }

  // create the load vector if one does not exist
  if (load == nullptr) 
    load = new Vector(numberNodes*NDF);

  // add -M * RV(accel) to the load vector
  load->addMatrixVector(1.0, mass, resid, -1.0);
  
  return 0;
}


const Vector&
EnhancedQuad::getResistingForce( ) 
{
  int tang_flag = 0 ; // don't get the tangent

  formResidAndTangent( tang_flag ) ;

  // subtract external loads 
  if (load != nullptr)
    resid -= *load;

  return resid ;   
}

// get residual with inertia terms
const Vector&  EnhancedQuad::getResistingForceIncInertia( )
{
  int tang_flag = 0 ; // don't get the tangent

  static Vector res(8);

  // do tangent and residual here 
  formResidAndTangent( tang_flag ) ;

  // inertia terms
  formInertiaTerms( tang_flag ) ;

  res = resid;

  // add the damping forces if rayleigh damping
  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
    res += this->getRayleighDampingForces();

  // subtract external loads 
  if (load != 0)
    res -= *load;

  return res;
}

//*********************************************************************
// form inertia terms

void   EnhancedQuad::formInertiaTerms( int tangFlag ) 
{

  double shp[nShape][numberNodes] ;  // shape functions at a gauss point

  static Vector momentum(NDF);
  mass.Zero( ) ;

  // gauss loop 
  for (int i = 0; i < nip; i++ ) {

    // get shape functions  
    double xsj; // determinant of jacaobian matrix
    shape2d( pts[i][0], pts[i][1], xl, shp, xsj ) ;
    
    // volume element
    double dvol = wts[i] * xsj * thickness;

    // node loop to compute acceleration
    momentum.Zero( ) ;
    for (int j = 0; j < numberNodes; j++ ) 
      // momentum += shp[massIndex][j] * ( theNodes[j]->getTrialAccel()  ) ; 
      momentum.addVector( 1.0,
                          theNodes[j]->getTrialAccel(),
                          shp[massIndex][j]) ;

    // density
    double rho = materialPointers[i]->getRho() ;

    // multiply acceleration by density to form momentum
    momentum *= rho ;


    // residual and tangent calculations node loops
    int jj = 0 ;
    for (int j = 0; j < numberNodes; j++ ) {

      double temp = shp[massIndex][j] * dvol ;

      for (int p = 0; p < NDF; p++ )
        resid( jj+p ) += ( temp * momentum(p) )  ;

      
      if ( tangFlag == 1 ) {

         // multiply by density
         temp *= rho ;

         // node-node mass
         int kk = 0 ;
         for (int k = 0; k < numberNodes; k++ ) {

            double massJK = temp * shp[massIndex][k] ;

            for (int p = 0; p < NDF; p++ )  
              mass( jj+p, kk+p ) += massJK ;
            
            kk += NDF ;
          }
      } // end if tang_flag 

      jj += NDF ;
    }
  }
}

//*********************************************************************
// form residual and tangent
int
EnhancedQuad::formResidAndTangent( int tang_flag ) 
{

  static constexpr double tolerance = 1.0e-08 ;
  static constexpr int nIterations = 10 ;

  int success ;

  OPS_STATIC double xsj[nip] ;  // determinant jacaobian matrix 
  OPS_STATIC double dvol[nip] ; // volume element
  OPS_STATIC double Shape[nip][nShape][NEN]; // [nip] ; // all the shape functions

  static Vector residJ(NDF) ; // nodeJ residual
  // static Matrix dd(nstress,nstress) ;  // material tangent

  static Matrix Kee(nEnhanced,nEnhanced) ;

  static Vector residE(nEnhanced);
  static Vector Umode(NDF);
  static Vector dalpha(nEnhanced);
  static Matrix Kue(numberDOF,nEnhanced) ;
  static Matrix Keu(nEnhanced,numberDOF) ;
  static Matrix KeeInvKeu(nEnhanced,numberDOF);
  
  // zero stiffness and residual 
  stiff.Zero( );
  resid.Zero( );

  Kee.Zero( );
  residE.Zero( );

  Kue.Zero( );
  Keu.Zero( );

  // compute Jacobian and inverse at center
  double L1 = 0.0 ;
  double L2 = 0.0 ;

  static MatrixND<ndm,ndm> J0, J0inv; //Jacobian matrix at center of element

  computeJacobian( L1, L2, xl, J0, J0inv ) ; 

  // gauss loop to compute and save shape functions 
  double det ;
  for (int i = 0; i < nip; i++ ) {

    // get shape functions    
    shape2d( pts[i][0], pts[i][1], xl, Shape[i], det ) ;

    // save jacobian determinant
    xsj[i] = det ;

    // volume element to also be saved
    dvol[i] = wts[i] * det * thickness;  

  }

  // -------------------------------------------------------------------
  // Newton loop to solve for enhanced strain parameters
  
  VectorND<nstress> stress[nip];
  MatrixND<nstress, nstress> dd[nip];
  static MatrixND<nstress, NDF> B[numberNodes];
  int count = 0 ;
  do {

    residE.Zero( ) ;
    Kee.Zero( ) ;

    // Gauss loop
    for (int i = 0; i < nip; i++ ) {

      // Interpolate strain

      VectorND<nstress> strain {};

      // j-node loop to compute nodal strain contributions
      for (int j = 0; j < NEN; j++ )  {

        // compute B matrix 
        computeB( j, Shape[i], B[j]) ;
      
        // nodal displacements 
        const Vector &ul = theNodes[j]->getTrialDisp( ) ;

        // compute the strain
        // strain += (BJ*ul) ; 
        strain.addMatrixVector(1.0, B[j], ul, 1.0) ;

      }

      // j-node loop to compute enhanced strain contributions
      for (int j = 0; j < nModes; j++ )  {

        MatrixND<nstress,NDF> BJ ;      // B matrix node J
        // compute B matrix 
        computeBenhanced( j, pts[i][0], pts[i][1], xsj[i], J0inv, BJ) ; 
      
        // enhanced "displacements" 
        Umode(0) = this->alpha( 2*j     ) ;
        Umode(1) = this->alpha( 2*j + 1 ) ;

        // compute the strain
        // strain += (BJ*Umode) ; 
        strain.addMatrixVector(1.0, BJ, Umode, 1.0) ;

      }

      success = materialPointers[i]->setTrialStrain( strain ) ;

      // compute the stress
      stress[i] = materialPointers[i]->getStress( ) ;

      // multiply by volume element
      stress[i]  *= dvol[i] ;

      // tangent 
      dd[i] = materialPointers[i]->getTangent( ) ;

      // multiply by volume element
      dd[i] *= dvol[i] ;


      // enhanced residual and tangent calculations loops

      int jj = 0 ;
      for (int j = 0; j < nModes; j++ ) {

        MatrixND<nstress,NDF> BJ;
        computeBenhanced( j, pts[i][0], pts[i][1], xsj[i], J0inv, BJ ) ; 

        // residual
        // residJ = -BJ^T * stress) ;
        residJ.addMatrixTransposeVector(0.0, BJ, stress[i], -1.0) ;

        // residual 
        for (int p = 0; p < NDF; p++ )
          residE( jj+p ) += residJ(p)  ;

        // 
        MatrixND<NDF, nstress> BJtranD = BJ^dd[i];

        int kk = 0 ;
        for (int k = 0; k < nModes; k++ ) {

          MatrixND<nstress,NDF> BK;
          computeBenhanced( k, pts[i][0], pts[i][1], xsj[i], J0inv, BK) ;
  
          MatrixND<NDF, NDF> stiffJK =  BJtranD * BK ;

          for (int p = 0; p < NDF; p++ )  {
             for (int q = 0; q < NDF; q++ )
                Kee( jj+p, kk+q ) += stiffJK( p, q ) ;
          }

          kk += NDF ;
        }

        jj += NDF ;
      }
    }


    // solve for enhanced strain parameters
    Kee.Solve( residE, dalpha ) ;

    if (dalpha(0) > 1.0e10)  
        return -1;

    this->alpha += dalpha ;

    count++ ;
    if ( count > nIterations ) {
      opserr << "Exceeded " << nIterations
             << " iterations solving for enhanced strain parameters "
             << endln ;
      break ;
    }

    // do at least 2 iterations so saved data is good
  } while ( residE.Norm() > tolerance  ||  count < 2 ) ;

  // End enhanced strain parameters Newton loop
  // -------------------------------------------------------------------

  //
  // Gauss loop 
  //
  for (int i = 0; i < nip; i++ ) {

    // residual and tangent calculations node loops

    int jj = 0 ;
    for (int j = 0; j < NEN; j++ ) {

      computeB( j, Shape[i], B[j]) ;

      // residual
      // residJ = BJtran * stress ;
      residJ.addMatrixTransposeVector(0.0, B[j], stress[i], 1.0) ;

      for (int p = 0; p < NDF; p++ )
        resid( jj+p ) += residJ(p);


      if ( tang_flag == 1 ) {

        //BJtranD = BJtran * dd ;
        MatrixND<NDF, nstress> BJtranD = B[j]^dd[i] ;

         // node-node stiffness
         int kk = 0 ;
         for (int k = 0; k < NEN; k++ ) {

            computeB( k, Shape[i], B[k]) ;
  
            const MatrixND<NDF, NDF> stiffJK =  BJtranD * B[k];

            for (int p = 0; p < NDF; p++ )  {
               for (int q = 0; q < NDF; q++ )
                  stiff( jj+p, kk+q ) += stiffJK( p, q ) ;
            }

            kk += NDF ;
          }

         // node-enhanced stiffness Kue 
         kk = 0 ;
         for (int k = 0; k < nModes; k++ ) {

            MatrixND<nstress,NDF> BK;
            computeBenhanced( k, pts[i][0], pts[i][1], xsj[i], J0inv, BK ) ;
  
            // stiffJK =  BJtranD * BK  ;
            const MatrixND<NDF, NDF> stiffJK = BJtranD * BK;
           
            for (int p = 0; p < NDF; p++ ) {
               for (int q = 0; q < NDF; q++ )
                  Kue( jj+p, kk+q ) += stiffJK( p, q ) ;
            }

            kk += NDF ;
          }

         // enhanced-node stiffness Keu 
         kk = 0;
         for (int k = 0; k < nModes; k++ ) {

            MatrixND<nstress,NDF> BK;
            computeBenhanced( k, pts[i][0], pts[i][1], xsj[i], J0inv, BK );

            const MatrixND<NDF, NDF> stiffKJ = (BK^dd[i])*B[j];

            for (int p = 0; p < NDF; p++ )  {
               for (int q = 0; q < NDF; q++ )
                  Keu( kk+p, jj+q ) += stiffKJ( p, q ) ;
            }  

            kk += NDF ;
          }

      } // end if tang_flag 

      jj += NDF ;
    }
  } 
  // end Gauss loop 


  //
  // static condensation of enhanced parameters
  //
  if ( tang_flag == 1 ) {  
     Kee.Solve( Keu, KeeInvKeu ) ;

     // stiff -= ( Kue * KeeInvKeu ) ;
     stiff.addMatrixProduct(1.0,  Kue, KeeInvKeu, -1.0 ) ;
  }

  return 0;
}

int  
EnhancedQuad::update() 
{
  // compute basis vectors and local nodal coordinates
  computeBasis();
  return 0;
}

void
EnhancedQuad::computeBasis() 
{
  // nodal coordinates
  for (int i = 0; i < 4; i++ ) {
     const Vector &coorI = theNodes[i]->getCrds( ) ;

     xl[0][i] = coorI(0);
     xl[1][i] = coorI(1);
  }
}

//************************************************************************
// shape function routine for four node quads
void
EnhancedQuad::shape2d(double ss, double tt, 
                      const double x[2][4], 
                      double shp[3][4], 
                      double &xsj)
{

  double temp ;
     
  static constexpr double s[] = { -0.5,  0.5, 0.5, -0.5 } ;
  static constexpr double t[] = { -0.5, -0.5, 0.5,  0.5 } ;

  static Matrix xs(2,2) ;
  static Matrix sx(2,2) ;

  for (int i = 0; i < 4; i++ ) {
      shp[2][i] = ( 0.5 + s[i]*ss )*( 0.5 + t[i]*tt ) ;
      shp[0][i] = s[i] * ( 0.5 + t[i]*tt ) ;
      shp[1][i] = t[i] * ( 0.5 + s[i]*ss ) ;
  }

  
  // Construct jacobian and its inverse
  
  for (int i = 0; i < 2; i++ ) {
    for (int j = 0; j < 2; j++ ) {
      xs(i,j) = 0.0;
      for (int k = 0; k < 4; k++ )
          xs(i,j) +=  x[i][k] * shp[j][k] ;

    }
  }

  xsj = xs(0,0)*xs(1,1) - xs(0,1)*xs(1,0) ;

  // inverse jacobian
  // inverse jacobian
  double jinv = 1.0 / xsj ;
  sx(0,0) =  xs(1,1) * jinv ;
  sx(1,1) =  xs(0,0) * jinv ;
  sx(0,1) = -xs(0,1) * jinv ;
  sx(1,0) = -xs(1,0) * jinv ;


  // form global derivatives 

  for (int i = 0; i < 4; i++ ) {
    temp      = shp[0][i]*sx(0,0) + shp[1][i]*sx(1,0) ;
    shp[1][i] = shp[0][i]*sx(0,1) + shp[1][i]*sx(1,1) ;
    shp[0][i] = temp ;
  }

  return ;
}


Response*
EnhancedQuad::setResponse(const char **argv, int argc, 
                          OPS_Stream &output)
{
  Response *theResponse = nullptr;

  output.tag("ElementOutput");
  output.attr("eleType","EnhancedQuad");
  output.attr("eleTag",this->getTag());
  output.attr("node1",connectedExternalNodes[0]);
  output.attr("node2",connectedExternalNodes[1]);
  output.attr("node3",connectedExternalNodes[2]);
  output.attr("node4",connectedExternalNodes[3]);

  char dataOut[10];
  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0) {
    
    for (int i=1; i<=4; i++) {
      sprintf(dataOut,"P1_%d",i);
      output.tag("ResponseType",dataOut);
      sprintf(dataOut,"P2_%d",i);
      output.tag("ResponseType",dataOut);
    }
    
    theResponse =  new ElementResponse(this, 1, resid);

  }  else if (strcmp(argv[0],"material") == 0 || strcmp(argv[0],"integrPoint") == 0) {
    int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= 4) {

      output.tag("GaussPoint");
      output.attr("number",pointNum);
      output.attr("eta",   pts[pointNum-1][0]);
      output.attr("neta",  pts[pointNum-1][1]);

      theResponse =  materialPointers[pointNum-1]->setResponse(&argv[2], argc-2, output);
      
      output.endTag();

    } 
  }
  else if ((strcmp(argv[0],"stress") == 0) || (strcmp(argv[0],"stresses") == 0)) {

      for (int i=0; i<4; i++) {
        output.tag("GaussPoint");
        output.attr("number",i+1);
        output.attr("eta",pts[i][0]);
        output.attr("neta",pts[i][1]);

        output.tag("NdMaterialOutput");
        output.attr("classType", materialPointers[i]->getClassTag());
        output.attr("tag", materialPointers[i]->getTag());

        output.tag("ResponseType","sigma11");
        output.tag("ResponseType","sigma22");
        output.tag("ResponseType","sigma12");

        output.endTag(); // GaussPoint
        output.endTag(); // NdMaterialOutput
      }

      theResponse =  new ElementResponse(this, 3, Vector(12));
  }
  
  else if ((strcmp(argv[0],"strain") == 0) || (strcmp(argv[0],"strains") == 0)) {

      for (int i=0; i<4; i++) {
        output.tag("GaussPoint");
        output.attr("number",i+1);
        output.attr("eta",pts[i][0]);
        output.attr("neta",pts[i][1]);

        output.tag("NdMaterialOutput");
        output.attr("classType", materialPointers[i]->getClassTag());
        output.attr("tag", materialPointers[i]->getTag());

        output.tag("ResponseType", "eta11");
        output.tag("ResponseType", "eta22");
        output.tag("ResponseType", "eta12");

        output.endTag(); // GaussPoint
        output.endTag(); // NdMaterialOutput
      }

      theResponse =  new ElementResponse(this, 4, Vector(12));
  }
        
  output.endTag(); // ElementOutput

  return theResponse;
}

int 
EnhancedQuad::getResponse(int responseID, Information &eleInfo)
{
  if (responseID == 1) {

    return eleInfo.setVector(this->getResistingForce());

  } else if (responseID == 3) {

    // Loop over the integration points
    static Vector stresses(12);
    int cnt = 0;
    for (int i = 0; i < 4; i++) {

      // Get material stress response
      const Vector &sigma = materialPointers[i]->getStress();
      stresses(cnt)   = sigma(0);
      stresses(cnt+1) = sigma(1);
      stresses(cnt+2) = sigma(2);
      cnt += 3;
    }
    return eleInfo.setVector(stresses);
  }

  else if (responseID == 4) {

    // Loop over the integration points
    static Vector stresses(12);
    int cnt = 0;
    for (int i = 0; i < 4; i++) {

      // Get material stress response
      const Vector &sigma = materialPointers[i]->getStrain();
      stresses(cnt) = sigma(0);
      stresses(cnt+1) = sigma(1);
      stresses(cnt+2) = sigma(2);
      cnt += 3;
    }
    return eleInfo.setVector(stresses);
        
  }
  else
    return -1;
}

int
EnhancedQuad::sendSelf (int commitTag, Channel &theChannel)
{
  int res = 0;
  
  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();
  
  // Quad packs its data into a Vector and sends this to theChannel
  // along with its dbTag and the commitTag passed in the arguments
  static Vector data(6);
  data(0) = this->getTag();
  data(1) = thickness;

  data(2) = alphaM;
  data(3) = betaK;
  data(4) = betaK0;
  data(5) = betaKc;
  
  res += theChannel.sendVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING EnhancedQuad::sendSelf() - " << this->getTag() << " failed to send Vector\n";
    return res;
  }              
  

  // Now quad sends the ids of its materials
  int matDbTag;
  
  static ID idData(12);

  for (int i = 0; i < 4; i++) {
    idData(i) = materialPointers[i]->getClassTag();
    matDbTag = materialPointers[i]->getDbTag();
    // NOTE: we do have to ensure that the material has a database
    // tag if we are sending to a database channel.
    if (matDbTag == 0) {
      matDbTag = theChannel.getDbTag();
                        if (matDbTag != 0)
                materialPointers[i]->setDbTag(matDbTag);
    }
    idData(i+4) = matDbTag;
  }
  
  idData(8) = connectedExternalNodes(0);
  idData(9) = connectedExternalNodes(1);
  idData(10) = connectedExternalNodes(2);
  idData(11) = connectedExternalNodes(3);

  res += theChannel.sendID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING EnhancedQuad::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return res;
  }

  // Finally, quad asks its material objects to send themselves
  for (int i = 0; i < 4; i++) {
    res += materialPointers[i]->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING EnhancedQuad::sendSelf() - " << this->getTag() << " failed to send its Material\n";
      return res;
    }
  }
  
  return res;
}

int
EnhancedQuad::recvSelf(int commitTag, Channel &theChannel, 
                       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  int dataTag = this->getDbTag();

  // Quad creates a Vector, receives the Vector and then sets the 
  // internal data with the data in the Vector
  static Vector data(6);
  res += theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING EnhancedQuad::recvSelf() - failed to receive Vector\n";
    return res;
  }
  
  this->setTag((int)data(0));
  thickness = data(1);

  alphaM = data(2);
  betaK = data(3);
  betaK0 = data(4);
  betaKc = data(5);

  static ID idData(12);
  // Quad now receives the tags of its four external nodes
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING EnhancedQuad::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return res;
  }

  connectedExternalNodes(0) = idData(8);
  connectedExternalNodes(1) = idData(9);
  connectedExternalNodes(2) = idData(10);
  connectedExternalNodes(3) = idData(11);
  
  if (materialPointers[0] == 0) {
    for (int i = 0; i < 4; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+4);
      // Allocate new material with the sent class tag
      materialPointers[i] = theBroker.getNewNDMaterial(matClassTag);
      if (materialPointers[i] == 0) {
            opserr << "EnhancedQuad::recvSelf() - Broker could not create NDMaterial of class type " << matClassTag << endln;
            return -1;
      }
      // Now receive materials into the newly allocated space
      materialPointers[i]->setDbTag(matDbTag);
      res += materialPointers[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
        opserr << "EnhancedQuad::recvSelf() - material " << i << "failed to recv itself\n";
            return res;
      }
    }
  }

  // materials exist , ensure materials of correct type and recvSelf on them
  else {
    for (int i = 0; i < 4; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+4);
      // Check that material is of the right type; if not,
      // delete it and create a new one of the right type
      if (materialPointers[i]->getClassTag() != matClassTag) {
        delete materialPointers[i];
        materialPointers[i] = theBroker.getNewNDMaterial(matClassTag);
        if (materialPointers[i] == nullptr) {
          opserr << "EnhancedQuad::recvSelf() - material " << i << "failed to create\n";
          return -1;
        }
      }
      // Receive the material
      materialPointers[i]->setDbTag(matDbTag);
      res += materialPointers[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
        opserr << "EnhancedQuad::recvSelf() - material " << i << "failed to recv itself\n";
        return res;
      }
    }
  }
  
  return res;
}


void
EnhancedQuad::Print( OPS_Stream &s, int flag )
{
  if (flag == OPS_PRINT_CURRENTSTATE) {
      s << endln;
      s << "Enhanced Strain Four Node Quad \n";
      s << "Element Number: " << this->getTag() << endln;
      s << "Node 1 : " << connectedExternalNodes(0) << endln;
      s << "Node 2 : " << connectedExternalNodes(1) << endln;
      s << "Node 3 : " << connectedExternalNodes(2) << endln;
      s << "Node 4 : " << connectedExternalNodes(3) << endln;
      s << "thickness : " << thickness << endln;
      s << "Material Information : \n ";
      materialPointers[0]->Print(s, flag);
      s << endln;
  }
  
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    const ID& node_tags = this->getExternalNodes();
    s << OPS_PRINT_JSON_ELEM_INDENT << "{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"" << this->getClassType() << "\", ";

    s << "\"nodes\": [";
    for (int i=0; i < NEN-1; i++)
        s << node_tags(i) << ", ";
    s << node_tags(NEN-1) << "]";
    s << ", ";

    s << "\"thickness\": " << thickness << ", ";
    s << "\"material\": [";
    for (int i = 0; i < nip - 1; i++)
      s << materialPointers[i]->getTag() << ", ";
    s << materialPointers[nip - 1]->getTag() << "]";
    s << "}";
    return;
  }
}
