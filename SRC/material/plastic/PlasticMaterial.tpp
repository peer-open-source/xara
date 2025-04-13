/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
//

#include <Information.h>
#include <Parameter.h>
#include <Constants.h>

#include <string.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <matrix/identity.h> // IbunI, IIdev

namespace OpenSees {

template <int n, PlaneType type, typename index>
Mate<n>*
PlasticMaterial<n,type,index>::getCopy()
{
  PlasticMaterial<n,type,index> * new_mat = new PlasticMaterial<n,type,index>();
  *new_mat = *this;
  return new_mat;
}


//send back the strain
template <int n, PlaneType type, typename index>
MatrixSD<n,true>
PlasticMaterial<n,type,index>::getStrain() 
{
  return strain ;
} 


// send back the stress 
template <int n, PlaneType type, typename index>
const MatrixSD<n>&
PlasticMaterial<n,type,index>::getStress() 
{
  return stress;
}

template<int ndim, PlaneType type, typename index>
int
PlasticMaterial<ndim,type,index>::setTrialStrain(const MatrixND<ndim,ndim>& F) 
{
  // Convert deformation gradient F to small strain tensor e
  MatrixSD<ndim,true> e;
  e.addMatrixTransposeProduct(0.0, F, F, 0.5);
  e.addDiagonal(-0.5);
  return setTrialStrain(e);
}


template <int n, PlaneType type, typename index>
// set the strain and integrate plasticity equations
int 
PlasticMaterial<n,type,index>::setTrialStrain( const MatrixSD<n,true> &e) 
{

  strain = e;

  this->plastic_integrator();

  return 0 ;
}


// send back the tangent 
template <int n, PlaneType type, typename index>
const Matrix& 
PlasticMaterial<n,type,index>::getTangent() 
{
  // matrix to tensor mapping
  //  Matrix      Tensor
  // -------     -------
  //   0           0 0
  //   1           1 1
  //   2           2 2   
  //   3           0 1  ( or 1 0 )
  //   4           1 2  ( or 2 1 )
  //   5           2 0  ( or 0 2 ) 

  for (int ii = 0; ii < 6; ii++ ) {
    for (int jj = 0; jj < 6; jj++ ) {

      int i, j, k, l ;
      index::map( ii, i, j ) ;
      index::map( jj, k, l ) ;

      tangent_matrix(ii,jj) = tangent[i][j][k][l] ;

    }
  }
  return tangent_matrix ;
} 

//send back the tangent 
template <int n, PlaneType type, typename index>
const Matrix& 
PlasticMaterial<n,type,index>::getInitialTangent() 
{
  // matrix to tensor mapping
  //  Matrix      Tensor
  // -------     -------
  //   0           0 0
  //   1           1 1
  //   2           2 2   
  //   3           0 1  ( or 1 0 )
  //   4           1 2  ( or 2 1 )
  //   5           2 0  ( or 0 2 ) 

  this->doInitialTangent();

  for (int ii = 0; ii < 6; ii++ ) {
    for (int jj = 0; jj < 6; jj++ ) {

      int i, j, k, l ;
      index::map( ii, i, j ) ;
      index::map( jj, k, l ) ;

      tangent_matrix(ii,jj) = initialTangent[i][j][k][l] ;

    }
  }

  return tangent_matrix ;
} 



//zero internal variables
template <int n, PlaneType type, typename index>
void
PlasticMaterial<n,type,index>::zero()
{
  xi_n      = 0.0;
  xi_nplus1 = 0.0;

  epsilon_p_n.zero();
  epsilon_p_nplus1.zero();

  stress.zero();
  strain.zero();
}


//null constructor
template <int n, PlaneType type, typename index>
PlasticMaterial<n,type,index>::PlasticMaterial()
 : Mate<n>(0),
   parameterID(0),
  
   bulk       (0.0),
   shear      (0.0),
   sigma_0    (0.0),
   sigma_infty(0.0),
   delta      (0.0),
   Hard       (0.0),
   eta        (0.0),
   rho        (0.0)
{

  this->zero(); // or (*this).zero( )

  plastic_integrator();
}


//full constructor
template <int n, PlaneType type, typename index>
PlasticMaterial<n,type,index>::PlasticMaterial(int tag,
                                  int classTag,
                                  double K,
                                  double G,
                                  double yield0,
                                  double yield_infty,
                                  double d,
                                  double H,
                                  double viscosity,
                                  double r)
 : Mate<n>(tag),
   parameterID(0)
{
  bulk        = K;
  shear       = G;
  sigma_0     = yield0;
  sigma_infty = yield_infty;
  delta       = d;
  Hard        = H;
  eta         = viscosity;
  rho         = r;

  this->zero();

  plastic_integrator();
}


//elastic constructor
template <int n, PlaneType type, typename index>
PlasticMaterial<n,type,index>::PlasticMaterial(int tag, int classTag, double K, double G)
 : Mate<n>(tag),
   parameterID(0)
{
  bulk        = K;
  shear       = G;
  sigma_0     = 1.0e16 * shear;
  sigma_infty = sigma_0;
  delta       = 0.0;
  Hard        = 0.0;
  eta         = 0.0;

  this->zero();
}


//destructor
template <int n, PlaneType type, typename index>
PlasticMaterial<n,type,index>::~PlasticMaterial() 
{

}

//print out material data
template <int n, PlaneType type, typename index>
void
PlasticMaterial<n,type,index>::Print(OPS_Stream& s, int flag)
{
  s << "\n";
  s << "J2-Plasticity : ";
  s << this->getClassType() << "\n";
  s << "Bulk Modulus =   " << bulk << "\n";
  s << "Shear Modulus =  " << shear << "\n";
  s << "Sigma_0 =        " << sigma_0 << "\n";
  s << "Sigma_infty =    " << sigma_infty << "\n";
  s << "Delta =          " << delta << "\n";
  s << "H =              " << Hard << "\n";
  s << "Eta =            " << eta << "\n";
  s << "Rho =            " << rho << "\n";
  s << "\n";
}


//--------------------Plasticity-------------------------------------

//plasticity integration routine
template <int n, PlaneType type, typename index>
void
PlasticMaterial<n,type,index>::plastic_integrator()
{
  const double tolerance = (1.0e-8) * sigma_0;

  const double dt = ops_Dt; //time step



  double inv_norm_tau = 0.0;
  double gamma        = 0.0; //consistency parameter
  double resid        = 1.0;
  double tang         = 0.0;

  double theta     = 0.0;
  double theta_inv = 0.0;


  constexpr static int max_iterations = 25;

  // compute the deviatoric strains

  double trace = strain(0, 0) + strain(1, 1) + strain(2, 2);

  // deviatoric strain
  Matrix3D dev_strain = strain.full();
  for (int i = 0; i < 3; i++)
    dev_strain(i, i) -= (1. / 3. * trace);

  // compute the trial deviatoric stresses
  //   dev_stress = (2.0*shear) * ( dev_strain - epsilon_p_n )
  Matrix3D dev_stress; 
  dev_stress  = dev_strain;
  dev_stress -= epsilon_p_n;
  dev_stress *= 2.0 * shear;

  // compute norm of deviatoric stress

  double norm_tau = 0.0; // norm of deviatoric stress
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++)
      norm_tau += dev_stress(i, j) * dev_stress(i, j);
  }

  norm_tau = std::sqrt(norm_tau);

  static Matrix3D normal; // normal to yield surface
  if (norm_tau > tolerance) {
    inv_norm_tau = 1.0 / norm_tau;
    normal       = inv_norm_tau * dev_stress;
  } else {
    normal.zero();
    inv_norm_tau = 0.0;
  }

  // compute trial value of yield function

  double phi = norm_tau - root23 * q(xi_n);

  int iteration_counter;
  if (phi > 0.0) { //plastic

    //solve for gamma
    gamma             = 0.0;
    resid             = 1.0;
    iteration_counter = 0;
    while (fabs(resid) > tolerance) {

      resid = norm_tau - (2.0 * shear) * gamma - root23 * q(xi_n + root23 * gamma);
      if (eta > 0.0 && dt > 0.0)
        resid -= (eta / dt) * gamma;

      tang = -(2.0 * shear) - 2. / 3. * qprime(xi_n + root23 * gamma);
      if (eta > 0.0 && dt > 0.0)
        tang -= (eta / dt);

      gamma -= (resid / tang);

      iteration_counter++;

      if (iteration_counter > max_iterations) {
        opserr << "More than " << max_iterations;
        opserr << " iterations in constituive subroutine J2-plasticity \n";
        break;
      }

    } //end while resid


    gamma *= (1.0 - 1e-08);

    // update plastic internal variables

    epsilon_p_nplus1  = epsilon_p_n;
    epsilon_p_nplus1 += gamma * normal;

    xi_nplus1 = xi_n + root23 * gamma;

    // recompute deviatoric stresses

    dev_stress = (2.0 * shear) * (dev_strain - epsilon_p_nplus1);

    // compute the terms for plastic part of tangent

    theta = (2.0 * shear) + 2. / 3. * qprime(xi_nplus1);

    if (eta > 0.0 && dt > 0.0)
      theta += (eta / dt);

    theta_inv = 1.0 / theta;

  } else { 
    // elastic

    // update history variables -- they remain unchanged

    epsilon_p_nplus1 = epsilon_p_n;

    xi_nplus1 = xi_n;

    // no extra tangent terms to compute

    gamma     = 0.0;
    theta     = 0.0;
    theta_inv = 0.0;

  } // end if phi > 0


  // add on bulk part of stress

  stress = dev_stress;
  stress.addDiagonal(bulk*trace);


  // compute the tangent

  const double c1 = -4.0 * shear * shear;
  const double c2 = c1 * theta_inv;
  const double c3 = c1 * gamma * inv_norm_tau;

  for (int ii = 0; ii < 6; ii++) {
    for (int jj = 0; jj < 6; jj++) {

      int i, j, k, l;
      index::map(ii, i, j);
      index::map(jj, k, l);

      double NbunN;              // normal bun normal
      NbunN = normal(i, j) * normal(k, l);

      // elastic terms
      tangent[i][j][k][l] = bulk * IbunI[i][j][k][l];

      tangent[i][j][k][l] += (2.0 * shear) * IIdev[i][j][k][l];

      // plastic terms
      tangent[i][j][k][l] += c2 * NbunN;

      tangent[i][j][k][l] += c3 * (IIdev[i][j][k][l] - NbunN);

      // minor symmetries
      tangent[j][i][k][l] = tangent[i][j][k][l];
      tangent[i][j][l][k] = tangent[i][j][k][l];
      tangent[j][i][l][k] = tangent[i][j][k][l];

    }
  }

  return;
}


// set up for initial elastic
template <int n, PlaneType type, typename index>
void
PlasticMaterial<n,type,index>::doInitialTangent()
{

  //compute the deviatoric strains
  for (int ii = 0; ii < 6; ii++) {
    for (int jj = 0; jj < 6; jj++) {

      int i, j, k, l;
      index::map(ii, i, j);
      index::map(jj, k, l);

      //elastic terms
      initialTangent[i][j][k][l] = bulk * IbunI[i][j][k][l];
      initialTangent[i][j][k][l] += (2.0 * shear) * IIdev[i][j][k][l];

      //minor symmetries
      initialTangent[j][i][k][l] = initialTangent[i][j][k][l];
      initialTangent[i][j][l][k] = initialTangent[i][j][k][l];
      initialTangent[j][i][l][k] = initialTangent[i][j][k][l];

    }
  } 

  return;
}


//hardening function
template <int n, PlaneType type, typename index>
double
PlasticMaterial<n,type,index>::q(double xi)
{
  //  q(xi) = simga_infty + (sigma_0 - sigma_infty)*exp(-delta*xi) + H*xi

  return sigma_infty + (sigma_0 - sigma_infty) * exp(-delta * xi) + Hard * xi;
}


//hardening function derivative
template <int n, PlaneType type, typename index>
double
PlasticMaterial<n,type,index>::qprime(double xi)
{
  return (sigma_0 - sigma_infty) * (-delta) * exp(-delta * xi) + Hard;
}


template <int n, PlaneType type, typename index>
PlaneType
PlasticMaterial<n,type,index>::getType() const
{
  return type;
}


template <int n, PlaneType type, typename index>
int
PlasticMaterial<n,type,index>::getOrder() const
{
  return n*(n+1)/2;
}


template <int n, PlaneType type, typename index>
int
PlasticMaterial<n,type,index>::commitState()
{
  epsilon_p_n = epsilon_p_nplus1;
  xi_n        = xi_nplus1;

  return 0;
}


template <int n, PlaneType type, typename index>
int
PlasticMaterial<n,type,index>::revertToLastCommit()
{
  return 0;
}


template <int n, PlaneType type, typename index>
int
PlasticMaterial<n,type,index>::revertToStart()
{

  // added: C.McGann, U.Washington for InitialStateAnalysis
  if (ops_InitialStateAnalysis) {
    // do nothing, keep state variables from last step
  } else {
    // normal call for revertToStart (not initialStateAnalysis)
    this->zero();
  }

  return 0;
}

#if 0
// TODO: need to modify Parameter class to work with
// Parameterized class instead of MovableObject
template <int n, PlaneType type, typename index>
int
PlasticMaterial<n,type,index>::setParameter(const char** argv, int argc, Parameter& param)
{
  if (strcmp(argv[0], "K") == 0)
    return param.addObject(1, this);

  else if (strcmp(argv[0], "G") == 0 || 
           strcmp(argv[0], "mu") == 0)
    return param.addObject(2, this);

  else if (strcmp(argv[0], "rho") == 0)
    return param.addObject(3, this);

  return -1;
}
#endif

template <int n, PlaneType type, typename index>
int
PlasticMaterial<n,type,index>::updateParameter(int parameterID, Information& info)
{
  switch (parameterID) {
  case 1:  bulk  = info.theDouble; return 0;
  case 2:  shear = info.theDouble; return 0;
  case 3:  rho   = info.theDouble; return 0;
  default: return -1;
  }
}

template <int n, PlaneType type, typename index>
int
PlasticMaterial<n,type,index>::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}


#if 0
template <int n, PlaneType type, typename index>
int
PlasticMaterial<n,type,index>::sendSelf(int commitTag, Channel& theChannel)
{
  // place all the data needed to define material and it's state
  // int a vector object
  static Vector data(10 + 9);
  int cnt     = 0;
  data(cnt++) = this->getTag();
  data(cnt++) = bulk;
  data(cnt++) = shear;
  data(cnt++) = sigma_0;
  data(cnt++) = sigma_infty;
  data(cnt++) = delta;
  data(cnt++) = Hard;
  data(cnt++) = eta;
  data(cnt++) = rho;

  data(cnt++) = xi_n;

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      data(cnt++) = epsilon_p_n(i, j);


  // send the vector object to the channel
  if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "PlasticMaterial<n,type,index>::sendSelf - failed to send vector to channel\n";
    return -1;
  }

  return 0;
}


template <int n, PlaneType type, typename index>
int
PlasticMaterial<n,type,index>::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
  // recv the vector object from the channel which defines material param and state
  static Vector data(10 + 9);
  if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "PlasticMaterial<n,type,index>::recvSelf - failed to recv vector from channel\n";
    return -1;
  }

  // set the material parameters and state variables
  int cnt = 0;
  this->setTag(data(cnt++));
  bulk        = data(cnt++);
  shear       = data(cnt++);
  sigma_0     = data(cnt++);
  sigma_infty = data(cnt++);
  delta       = data(cnt++);
  Hard        = data(cnt++);
  eta         = data(cnt++);
  rho         = data(cnt++);

  xi_n = data(cnt++);

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      epsilon_p_n(i, j) = data(cnt++);

  epsilon_p_nplus1 = epsilon_p_n;
  xi_nplus1        = xi_n;

  return 0;
}
#endif


//
// ---------------------------------------------------
//
struct index {
  static void map(int matrix_index, int& i, int& j)
  {
    switch (matrix_index + 1) { //add 1 for standard tensor indices

    case 1:
      i = 1;
      j = 1;
      break;

    case 2:
      i = 2;
      j = 2;
      break;

    case 3:
      i = 3;
      j = 3;
      break;

    case 4:
      i = 1;
      j = 2;
      break;

    case 5:
      i = 2;
      j = 3;
      break;

    case 6:
      i = 3;
      j = 1;
      break;


    default:
      i = 1;
      j = 1;
      break;

    } //end switch

    i--; //subtract 1 for C-indexing
    j--;

    return;
  }
};

int
test()
{
  PlasticMaterial<3,PlaneType::None,index> material;
}

}
