/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
//
// Description: This file contains the class definition for LagrangeQuad.
// See https://portwooddigital.com/2022/09/11/unrolling-the-four-node-quad/
//
#include <LagrangeQuad.h>
#include <Node.h>
#include <Mate.h>
#include <Matrix.h>
#include <MatrixSD.h>
#include <Vector.h>
#include <VectorND.h>
#include <ID.h>
#include <Domain.h>
#include <string.h>
#include <Information.h>
#include <Parameter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <cassert>


namespace OpenSees {


template <int NEN, bool enh>
LagrangeQuad<NEN,enh>::LagrangeQuad(int tag, 
                                    const std::array<int, NEN>& nodes,
                                    Mate<2>& m,
                                    double thickness, 
                                    double p, double r, double b1, double b2)
 : Element(tag, ELE_TAG_LagrangeQuad),
   connectedExternalNodes(NEN),
   pressureLoad(NEN*NDF),
   thickness(thickness),
   applyLoad(0),
   pressure(p),
   rho(r),
   Ki(0)
{

  // Body forces
  b[0] = b1;
  b[1] = b2;

  for (int i = 0; i < NIP; i++) {

    // Get copies of the material model for each integration point
    theMaterial[i] = m.getCopy();

    // Check allocation
    if (theMaterial[i] == nullptr) {
      opserr << "LagrangeQuad::LagrangeQuad -- failed to get a copy of material model\n";
      return;
    }
  }

  // Set connected external node IDs
  for (int i = 0; i < NEN; i++) {
    theNodes[i]               = nullptr;
    connectedExternalNodes(i) = nodes[i];
  }

  Q.zero();
}

template <int NEN, bool enh>
LagrangeQuad<NEN,enh>::LagrangeQuad()
 : Element(0, ELE_TAG_LagrangeQuad),
   connectedExternalNodes(NEN),
   pressureLoad(NEN*NDF),
   thickness(0.0),
   applyLoad(0),
   pressure(0.0),
   Ki(nullptr)
{
  for (int i = 0; i < NEN; i++)
    theNodes[i] = nullptr;

  Q.zero();
}


template <int NEN, bool enh>
LagrangeQuad<NEN,enh>::~LagrangeQuad()
{
  for (int i = 0; i < nip; i++) {
    if (theMaterial[i])
      delete theMaterial[i];
  }

  if (Ki != nullptr)
    delete Ki;
}


template <int NEN, bool enh>
int
LagrangeQuad<NEN,enh>::getNumExternalNodes() const
{
  return NEN;
}

template <int NEN, bool enh>
const ID&
LagrangeQuad<NEN,enh>::getExternalNodes()
{
  return connectedExternalNodes;
}


template <int NEN, bool enh>
Node**
LagrangeQuad<NEN,enh>::getNodePtrs()
{
  return &theNodes[0];
}


template <int NEN, bool enh>
int
LagrangeQuad<NEN,enh>::getNumDOF()
{
  int sum = 0;

  for (const Node* node : theNodes)
    sum += node->getNumberDOF();

  return sum;
}

template <int NEN, bool enh>
void
LagrangeQuad<NEN,enh>::setDomain(Domain* theDomain)
{
  // Check Domain is not null. This happens when element is removed from a domain.
  // In this case just set null pointers to null and return.
  if (theDomain == nullptr) {
    for (int i = 0; i < NEN; i++)
      theNodes[i] = nullptr;
    return;
  }

  for (int i = 0; i < NEN; i++) {
    // Retrieve the node from the domain using its tag.
    // If no node is found, then return
    theNodes[i] = theDomain->getNode(connectedExternalNodes(i));
    if (theNodes[i] == nullptr)
      return;

    // If node is found, ensure node has the proper number of DOFs
    int dofs = theNodes[i]->getNumberDOF();
    if (dofs != NDF) {
      opserr << "WARNING element " << this->getTag() << " does not have " << NDF << " DOFs at node "
             << theNodes[i]->getTag() << ".\n";
      return;
    }
  }
  this->DomainComponent::setDomain(theDomain);

  // Compute consistent nodal loads due to pressure
  this->setPressureLoadAtNodes();
}


template <int NEN, bool enh>
int
LagrangeQuad<NEN,enh>::commitState()
{
  int retVal = 0;

  // call element commitState to do any base class stuff
  if ((retVal = this->Element::commitState()) != 0) {
    opserr << "LagrangeQuad::commitState () - failed in base class";
  }

  // Loop over the integration points and commit the material states
  for (int i = 0; i < NIP; i++)
    retVal += theMaterial[i]->commitState();

  return retVal;
}

template <int NEN, bool enh>
int
LagrangeQuad<NEN,enh>::revertToLastCommit()
{
  int retVal = 0;

  // Loop over the integration points and revert to last committed state
  for (int i = 0; i < NIP; i++)
    retVal += theMaterial[i]->revertToLastCommit();

  return retVal;
}


template <int NEN, bool enh>
int
LagrangeQuad<NEN,enh>::revertToStart()
{
  int retVal = 0;

  // Loop over the integration points and revert states to start
  for (int i = 0; i < nip; i++)
    retVal += theMaterial[i]->revertToStart();

  return retVal;
}


template <int NEN, bool enh>
int
LagrangeQuad<NEN,enh>::update()
{
  // Collect displacements at each node into a local array
  double u[NDM][NEN];

  for (int i = 0; i < NEN; i++) {
    const Vector& displ = theNodes[i]->getTrialDisp();
    for (int j = 0; j < NDM; j++) {
      u[j][i] = displ[j];
    }
  }

  int ret = 0;

  // Loop over the integration points
  for (int i = 0; i < nip; i++) {
    // Determine Jacobian for this integration point
    this->shapeFunction(pts[i][0], pts[i][1]);

    // Interpolate strains
    //   eps = B*u;
    MatrixSD<2, true> eps{{0.0}};

    for (int beta = 0; beta < NEN; beta++) {
      eps.vector[0] += shp[0][beta] * u[0][beta];
      eps.vector[1] += shp[1][beta] * u[1][beta];
      eps.vector[2] += shp[0][beta] * u[1][beta] + shp[1][beta] * u[0][beta];
    }

    // Set the material strain
    ret += theMaterial[i]->setTrialStrain(eps);
  }

  return ret;
}


template <int NEN, bool enh>
const Matrix&
LagrangeQuad<NEN,enh>::getTangentStiff()
{
  static MatrixND<NDF*NEN,NEN*NDF> K;
  K.zero();

  double DB[3][2];

  // Loop over the integration points
  for (int i = 0; i < nip; i++) {

    // Determine Jacobian for this integration point
    double dvol = this->shapeFunction(pts[i][0], pts[i][1]);
    dvol *= (thickness * wts[i]);

    // Get the material tangent
    const MatrixSD<3> D = theMaterial[i]->getTangent();

    // Perform numerical integration
    //
    // K = K + (B' * D * B) * intWt(i)*intWt(j) * detJ;

    for (int alpha = 0, ia = 0; alpha < 4; alpha++, ia += 2) {
      for (int beta = 0, ib = 0; beta < 4; beta++, ib += 2) {

        DB[0][0] = dvol * (D(0,0) * shp[0][beta] + D(0,2) * shp[1][beta]);
        DB[1][0] = dvol * (D(1,0) * shp[0][beta] + D(1,2) * shp[1][beta]);
        DB[2][0] = dvol * (D(2,0) * shp[0][beta] + D(2,2) * shp[1][beta]);
        DB[0][1] = dvol * (D(0,1) * shp[1][beta] + D(0,2) * shp[0][beta]);
        DB[1][1] = dvol * (D(1,1) * shp[1][beta] + D(1,2) * shp[0][beta]);
        DB[2][1] = dvol * (D(2,1) * shp[1][beta] + D(2,2) * shp[0][beta]);


        K(ia,     ib)     += shp[0][alpha] * DB[0][0] + shp[1][alpha] * DB[2][0];
        K(ia,     ib + 1) += shp[0][alpha] * DB[0][1] + shp[1][alpha] * DB[2][1];
        K(ia + 1, ib)     += shp[1][alpha] * DB[1][0] + shp[0][alpha] * DB[2][0];
        K(ia + 1, ib + 1) += shp[1][alpha] * DB[1][1] + shp[0][alpha] * DB[2][1];
      }
    }
  }

  static Matrix Wrapper(K);
  return Wrapper;
}


template <int NEN, bool enh>
const Matrix&
LagrangeQuad<NEN,enh>::getInitialStiff()
{
  if (Ki != 0)
    return *Ki;

  Ki = new Matrix(NEN*NDF, NEN*NDF);

  Ki->Zero();

  double DB[3][2];

  // Loop over the integration points
  for (int i = 0; i < nip; i++) {

    // Determine Jacobian for this integration point
    double dvol;
    dvol = this->shapeFunction(pts[i][0], pts[i][1]);
    dvol *= (thickness * wts[i]);

    // Get the material tangent
    MatrixSD<3> D = theMaterial[i]->getInitialTangent();

    // Perform numerical integration
    //
    // K = K + (B^ D * B) * intWt(i)*intWt(j) * detJ;
    //
    for (int beta = 0, ib = 0, colIb = 0, colIbP1 = 8; beta < 4;
         beta++, ib += 2, colIb += 16, colIbP1 += 16) {

      for (int alpha = 0, ia = 0; alpha < 4; alpha++, ia += 2) {

        DB[0][0] = dvol * (D(0,0) * shp[0][beta] + D(0,2) * shp[1][beta]);
        DB[1][0] = dvol * (D(1,0) * shp[0][beta] + D(1,2) * shp[1][beta]);
        DB[2][0] = dvol * (D(2,0) * shp[0][beta] + D(2,2) * shp[1][beta]);
        DB[0][1] = dvol * (D(0,1) * shp[1][beta] + D(0,2) * shp[0][beta]);
        DB[1][1] = dvol * (D(1,1) * shp[1][beta] + D(1,2) * shp[0][beta]);
        DB[2][1] = dvol * (D(2,1) * shp[1][beta] + D(2,2) * shp[0][beta]);

        (*Ki)(ia, ib) += shp[0][alpha] * DB[0][0] + shp[1][alpha] * DB[2][0];
        (*Ki)(ia, ib + 1) += shp[0][alpha] * DB[0][1] + shp[1][alpha] * DB[2][1];
        (*Ki)(ia + 1, ib) += shp[1][alpha] * DB[1][0] + shp[0][alpha] * DB[2][0];
        (*Ki)(ia + 1, ib + 1) += shp[1][alpha] * DB[1][1] + shp[0][alpha] * DB[2][1];
      }
    }
  }

  return *Ki;
}


template <int NEN, bool enh>
const Matrix&
LagrangeQuad<NEN,enh>::getMass()
{
  static MatrixND<NEN*NDF, NEN*NDF> M;
  static Matrix Wrapper(M);
  M.zero();

  static double rhoi[4];
  double sum = 0.0;
  for (int i = 0; i < nip; i++) {
    if (rho == 0)
      rhoi[i] = theMaterial[i]->getDensity();
    else
      rhoi[i] = rho;
    sum += rhoi[i];
  }

  if (sum == 0.0)
    return Wrapper;

  // Compute a lumped mass matrix
  for (int i = 0; i < nip; i++) {
    // Determine Jacobian for this integration point
    double rhodvol = this->shapeFunction(pts[i][0], pts[i][1]);

    // Element plus material density ... MAY WANT TO REMOVE ELEMENT DENSITY
    rhodvol *= (rhoi[i] * thickness * wts[i]);

    for (int alpha = 0, ia = 0; alpha < 4; alpha++, ia++) {
      double Nrho = shp[2][alpha] * rhodvol;
      M(ia, ia) += Nrho;
      ia++;
      M(ia, ia) += Nrho;
    }
  }

  return Wrapper;
}


template <int NEN, bool enh>
void
LagrangeQuad<NEN,enh>::zeroLoad()
{
  Q.zero();

  applyLoad = 0;

  appliedB[0] = 0.0;
  appliedB[1] = 0.0;

  return;
}

template <int NEN, bool enh>
int
LagrangeQuad<NEN,enh>::addLoad(ElementalLoad* theLoad, double loadFactor)
{
  // Added option for applying body forces in load pattern: C.McGann, U.Washington
  int type;
  const Vector& data = theLoad->getData(type, loadFactor);

  if (type == LOAD_TAG_SelfWeight) {
    applyLoad = 1;
    appliedB[0] += loadFactor * data(0) * b[0];
    appliedB[1] += loadFactor * data(1) * b[1];
    return 0;
  } else {
    opserr << "LagrangeQuad::addLoad - load type unknown for ele with tag: " << this->getTag()
           << endln;
    return -1;
  }

  return -1;
}

template <int NEN, bool enh>
int
LagrangeQuad<NEN,enh>::addInertiaLoadToUnbalance(const Vector& accel)
{
  double rhoi[nip];
  double sum = 0.0;
  for (int i = 0; i < nip; i++) {
    rhoi[i] = theMaterial[i]->getDensity();
    sum += rhoi[i];
  }

  if (sum == 0.0)
    return 0;

  static double ra[NEN*NDF];

  // Get R * accel from the nodes
  for (int i=0; i<NEN; i++) {
    const Vector& Raccel = theNodes[0]->getRV(accel);
    assert(Raccel.Size() == NDF);
    for (int j=0; j<NDF; j++)
      ra[NDF*i+j] = Raccel[j];
  }

  // Compute mass matrix
  const Matrix& M = this->getMass();

  // Add ( - fact * M R * accel ) to unbalance
  // Take advantage of lumped mass matrix
  for (int i = 0; i < NEN*NDF; i++)
      Q[i] += -M(i, i) * ra[i];

  return 0;
}


template <int NEN, bool enh>
const Vector&
LagrangeQuad<NEN,enh>::getResistingForce()
{
  thread_local VectorND<NEN*NDF> P;
  P.zero();


  // Loop over the integration points
  for (int i = 0; i < nip; i++) {
    // Determine Jacobian for this integration point
    double dvol;
    dvol = this->shapeFunction(pts[i][0], pts[i][1]);
    dvol *= (thickness * wts[i]);

    // Get material stress response
    const VectorND<3> sigma = theMaterial[i]->getStress().vector;

    // Perform numerical integration on internal force
    //P = P + (B^ sigma) * intWt(i)*intWt(j) * detJ;
    //P.addMatrixTransposeVector(1.0, B, sigma, intWt(i)*intWt(j)*detJ);
    for (int alpha = 0, ia = 0; alpha < NEN; alpha++, ia += 2) {
      P(ia)     += dvol * (shp[0][alpha] * sigma[0] + shp[1][alpha] * sigma[2]);
      P(ia + 1) += dvol * (shp[1][alpha] * sigma[1] + shp[0][alpha] * sigma[2]);

      // Subtract equiv. body forces from the nodes
      //P = P - (N^ b) * intWt(i)*intWt(j) * detJ;
      //P.addMatrixTransposeVector(1.0, N, b, -intWt(i)*intWt(j)*detJ);
      if (applyLoad == 0) {
        P(ia)     -= dvol * (shp[2][alpha] * b[0]);
        P(ia + 1) -= dvol * (shp[2][alpha] * b[1]);
      } else {
        P(ia)     -= dvol * (shp[2][alpha] * appliedB[0]);
        P(ia + 1) -= dvol * (shp[2][alpha] * appliedB[1]);
      }
    }
  }

  // Subtract pressure loading from resisting force
  if (pressure != 0.0) {
    //P = P - pressureLoad;
    P.addVector(1.0, pressureLoad, -1.0);
  }

  // Subtract other external nodal loads ... P_res = P_int - P_ext
  // P = P - Q;
  P.addVector(1.0, Q, -1.0);

  thread_local Vector wrapper(P);
  return wrapper;
}

template <int NEN, bool enh>
const Vector&
LagrangeQuad<NEN,enh>::getResistingForceIncInertia()
{
  static Vector P(NEN*NDF);

  double rhoi[nip];
  double sum = 0.0;
  for (int i = 0; i < nip; i++) {
    rhoi[i] = theMaterial[i]->getDensity();
    sum += rhoi[i];
  }

  // Compute the current resisting force
  P = this->getResistingForce();

  // if no mass terms .. just add damping terms
  if (sum == 0.0) {
    // add the damping forces if rayleigh damping
    if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
      P += this->getRayleighDampingForces();

    return P;
  }

  double a[NDF*NEN];
  for (int i=0; i<NEN; i++) {
      const Vector &accel = theNodes[i]->getTrialAccel();
      for (int j=0; j<NDF; j++)
        a[i*NDF+j] = accel[j];
  }

  // Compute the mass matrix
  const Matrix& M = this->getMass();

  // Take advantage of lumped mass matrix
  for (int i = 0; i < NEN*NDF; i++)
    P[i] += M(i, i)*a[i];

  // add the damping forces if rayleigh damping
  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
    P += this->getRayleighDampingForces();

  return P;
}


template <int NEN, bool enh>
int
LagrangeQuad<NEN,enh>::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  
  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();
  
  // Quad packs its data into a Vector and sends this to theChannel
  // along with its dbTag and the commitTag passed in the arguments
  static Vector data(9);
  data(0) = this->getTag();
  data(1) = thickness;
  data(2) = b[0];
  data(3) = b[1];
  data(4) = pressure;

  data(5) = alphaM;
  data(6) = betaK;
  data(7) = betaK0;
  data(8) = betaKc;
  
  res += theChannel.sendVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING LagrangeQuad::sendSelf() - " << this->getTag() << " failed to send Vector\n";
    return res;
  }              
  

  // Now send the ids of our materials
  
  static ID idData(NEN + 2*nip);
  
  for (int i = 0; i < nip; i++) {
#if 0 // TODO
    idData(i) = theMaterial[i]->getClassTag();
    int matDbTag = theMaterial[i]->getDbTag();
    // NOTE: we do have to ensure that the material has a database
    // tag if we are sending to a database channel.
    if (matDbTag == 0) {
      matDbTag = theChannel.getDbTag();
      if (matDbTag != 0)
        theMaterial[i]->setDbTag(matDbTag);
    }
    idData(i+4) = matDbTag;
#endif
  }
  
  for (int i=0; i<NEN; i++)
      idData(8+i) = connectedExternalNodes(i);

  res += theChannel.sendID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING LagrangeQuad::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return res;
  }

  // Finally, asks our material objects to send themselves
  for (int i = 0; i < nip; i++) {
#if 0
    res += theMaterial[i]->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING LagrangeQuad::sendSelf() - " << this->getTag() << " failed to send its Material\n";
      return res;
    }
#endif
  }
  
  return res;
}


template <int NEN, bool enh>
int
LagrangeQuad<NEN,enh>::recvSelf(int commitTag, Channel &theChannel,
                       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  int dataTag = this->getDbTag();

  // Quad creates a Vector, receives the Vector and then sets the 
  // internal data with the data in the Vector
  static Vector data(9);
  res += theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING LagrangeQuad::recvSelf() - failed to receive Vector\n";
    return res;
  }
  
  this->setTag((int)data(0));
  thickness = data(1);
  b[0] = data(2);
  b[1] = data(3);
  pressure = data(4);

  alphaM = data(5);
  betaK  = data(6);
  betaK0 = data(7);
  betaKc = data(8);

  static ID idData(12);
  // Quad now receives the tags of its four external nodes
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING LagrangeQuad::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return res;
  }

  for (int i=0; i < NEN; i++)
    connectedExternalNodes(i) = idData(8+i);
  
  if (theMaterial[0] == nullptr) {
    // Allocate new materials

    for (int i = 0; i < nip; i++) {
#if 0 // TODO
      int matClassTag = idData(i);
      int matDbTag = idData(i+4);
      // Allocate new material with the sent class tag
      theMaterial[i] = theBroker.getNewMate<2>(matClassTag);
      if (theMaterial[i] == nullptr) {
          opserr << "LagrangeQuad::recvSelf() - Broker could not create Mate<2> of class type " << matClassTag << endln;
          return -1;
      }
      // Now receive materials into the newly allocated space
      theMaterial[i]->setDbTag(matDbTag);
      res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);

      if (res < 0) {
        opserr << "LagrangeQuad::recvSelf() - material " << i << "failed to recv itself\n";
        return res;
      }
#endif
    }
  }

  // materials exist , ensure materials of correct type and recvSelf on them
  else {
    for (int i = 0; i < nip; i++) {
#if 0
      int matClassTag = idData(i);
      int matDbTag = idData(i+4);
      // Check that material is of the right type; if not,
      // delete it and create a new one of the right type
      if (theMaterial[i]->getClassTag() != matClassTag) {
          delete theMaterial[i];
//        theMaterial[i] = theBroker.getNewMate<2>(matClassTag);
          if (theMaterial[i] == nullptr) {
            opserr << "LagrangeQuad::recvSelf() - material " << i << "failed to create\n";
            return -1;
          }
      }
      // Receive the material
      theMaterial[i]->setDbTag(matDbTag);
      res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
        opserr << "LagrangeQuad::recvSelf() - material " << i << "failed to recv itself\n";
            return res;
      }
#endif
    }
  }
  
  return res;
}


template <int NEN, bool enh>
void
LagrangeQuad<NEN,enh>::Print(OPS_Stream& s, int flag)
{
  const ID& node_tags = this->getExternalNodes();

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_ELEM_INDENT << "{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"" << this->getClassType() << "\", ";

    s << "\"nodes\": [";
    for (int i = 0; i < NEN - 1; i++)
      s << node_tags(i) << ", ";
    s << node_tags(NEN - 1) << "]";
    s << ", ";

    s << "\"thickness\": " << thickness << ", ";
    s << "\"surfacePressure\": " << pressure << ", ";
    s << "\"masspervolume\": " << rho << ", ";
    s << "\"bodyForces\": [" << b[0] << ", " << b[1] << "], ";
    
    s << "\"material\": [";
    for (int i = 0; i < nip - 1; i++)
      s << theMaterial[i]->getTag() << ", ";
    s << theMaterial[nip - 1]->getTag() << "] ";
    s << "}";

    return;
  }

  if (flag == OPS_PRINT_CURRENTSTATE) {
    s << "\nLagrangeQuad, element id:  " << this->getTag() << endln;
    s << "\tConnected external nodes:  " << connectedExternalNodes;
    s << "\tthickness:  " << thickness << endln;
    s << "\tsurface pressure:  " << pressure << endln;
    s << "\tmass density:  " << rho << endln;
    s << "\tbody forces:  " << b[0] << " " << b[1] << endln;
    theMaterial[0]->Print(s, flag);
    s << "\tStress (xx yy xy)" << endln;
    for (int i = 0; i < 4; i++) {
      auto sigma = theMaterial[i]->getStress();
      s << "\t\tGauss point " << i + 1 << ": " << Vector(sigma.vector);
    }
  }
}

template <int NEN, bool enh>
Response*
LagrangeQuad<NEN,enh>::setResponse(const char** argv, int argc, OPS_Stream& output)
{
  Response* theResponse = nullptr;

  output.tag("ElementOutput");
  output.attr("eleType", "LagrangeQuad");
  output.attr("eleTag", this->getTag());
  output.attr("node1", connectedExternalNodes[0]);
  output.attr("node2", connectedExternalNodes[1]);
  output.attr("node3", connectedExternalNodes[2]);
  output.attr("node4", connectedExternalNodes[3]);

  char dataOut[10];
  if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "forces") == 0) {

    for (int i = 1; i <= 4; i++) {
      sprintf(dataOut, "P1_%d", i);
      output.tag("ResponseType", dataOut);
      sprintf(dataOut, "P2_%d", i);
      output.tag("ResponseType", dataOut);
    }

    theResponse = new ElementResponse(this, 1, Vector(NDF*NEN));
  }

  else if (strcmp(argv[0], "material") == 0 || strcmp(argv[0], "integrPoint") == 0) {

    int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= 4) {

      output.tag("GaussPoint");
      output.attr("number", pointNum);
      output.attr("eta",  pts[pointNum - 1][0]);
      output.attr("neta", pts[pointNum - 1][1]);

      theResponse = theMaterial[pointNum - 1]->setResponse(&argv[2], argc - 2, output);

      output.endTag();
    }
  } 
  
  else if ((strcmp(argv[0], "stresses") == 0) || (strcmp(argv[0], "stress") == 0)) {
    for (int i = 0; i < 4; i++) {
      output.tag("GaussPoint");
      output.attr("number", i + 1);
      output.attr("eta", pts[i][0]);
      output.attr("neta", pts[i][1]);

      output.tag("NdMaterialOutput");
      output.attr("classType", theMaterial[i]->getClassTag());
      output.attr("tag", theMaterial[i]->getTag());

      output.tag("ResponseType", "sigma11");
      output.tag("ResponseType", "sigma22");
      output.tag("ResponseType", "sigma12");

      output.endTag(); // GaussPoint
      output.endTag(); // NdMaterialOutput
    }
    theResponse = new ElementResponse(this, 3, Vector(12));
  }

  else if ((strcmp(argv[0], "stressesAtNodes") == 0) || (strcmp(argv[0], "stressAtNodes") == 0)) {
    for (int i = 0; i < 4; i++) { // nnodes
      output.tag("NodalPoint");
      output.attr("number", i + 1);
      // output.attr("eta",pts[i][0]);
      // output.attr("neta",pts[i][1]);

      // output.tag("NdMaterialOutput");
      // output.attr("classType", theMaterial[i]->getClassTag());
      // output.attr("tag", theMaterial[i]->getTag());

      output.tag("ResponseType", "sigma11");
      output.tag("ResponseType", "sigma22");
      output.tag("ResponseType", "sigma12");

      output.endTag(); // GaussPoint
      // output.endTag(); // NdMaterialOutput
    }
    theResponse = new ElementResponse(this, 11, Vector(12)); // 3 * nnodes
  }

  else if ((strcmp(argv[0], "strain") == 0) || (strcmp(argv[0], "strains") == 0)) {
    for (int i = 0; i < 4; i++) {
      output.tag("GaussPoint");
      output.attr("number", i + 1);
      output.attr("eta", pts[i][0]);
      output.attr("neta", pts[i][1]);

      output.tag("NdMaterialOutput");
      output.attr("classType", theMaterial[i]->getClassTag());
      output.attr("tag", theMaterial[i]->getTag());

      output.tag("ResponseType", "eta11");
      output.tag("ResponseType", "eta22");
      output.tag("ResponseType", "eta12");

      output.endTag(); // GaussPoint
      output.endTag(); // NdMaterialOutput
    }
    theResponse = new ElementResponse(this, 4, Vector(12));
  }

  output.endTag(); // ElementOutput

  return theResponse;
}

template <int NEN, bool enh>
int
LagrangeQuad<NEN,enh>::getResponse(int responseID, Information& eleInfo)
{
  if (responseID == 1) {

    return eleInfo.setVector(this->getResistingForce());

  } else if (responseID == 3) {

    // Loop over the integration points
    static Vector stresses(12);
    int cnt = 0;
    for (int i = 0; i < 4; i++) {
      // Get material stress response
      auto sigma        = theMaterial[i]->getStress().vector;
      stresses(cnt)     = sigma[0];
      stresses(cnt + 1) = sigma[1];
      stresses(cnt + 2) = sigma[2];
      cnt += 3;
    }

    return eleInfo.setVector(stresses);

  } else if (responseID == 11) {
    // extrapolate stress from Gauss points to element nodes
    constexpr static int ns = 3; // number of stress components

    static Vector stressGP(3*NIP);
    static Vector stressAtNodes(3*NEN);
    stressAtNodes.Zero();
    int cnt = 0;
    // first get stress components (xx, yy, xy) at Gauss points
    for (int i = 0; i < 4; i++) { // nip
      // Get material stress response
      const VectorND<3> sigma = theMaterial[i]->getStress().vector;
      stressGP(cnt)           = sigma[0];
      stressGP(cnt + 1)       = sigma[1];
      stressGP(cnt + 2)       = sigma[2];
      cnt += 3;
    }

    // [nnodes][nip]
    constexpr double We[NEN][NIP] = {{1.8660254037844386, -0.5, 0.1339745962155614, -0.5},
                                     {-0.5, 1.8660254037844386, -0.5, 0.1339745962155614},
                                     {0.1339745962155614, -0.5, 1.8660254037844386, -0.5},
                                     {-0.5, 0.1339745962155614, -0.5, 1.8660254037844386}};

    for (int i = 0; i < NEN; i++) {   // nnodes
      for (int k = 0; k < ns; k++) { // number of stress components
        int p = 3 * i + k;
        for (int j = 0; j < NIP; j++) { // nip
          // l is 
          int l = 3 * j + k;
          stressAtNodes(p) += We[i][j] * stressGP(l);
        }
      }
    }

    return eleInfo.setVector(stressAtNodes);
  }

  else if (responseID == 4) {

#if 0
    // Loop over the integration points
    static Vector stresses(12);
    int cnt = 0;
    for (int i = 0; i < nip; i++) {

      // Get material stress response
      const Vector &sigma = theMaterial[i]->getStrain();
      
      stresses(cnt) = sigma[0];
      stresses(cnt+1) = sigma[1];
      stresses(cnt+2) = sigma[2];
      cnt += 3;
    }

    return eleInfo.setVector(stresses);
#endif

  }

  return -1;
}


#if 0

template <int NEN, bool enh>
int
LagrangeQuad<NEN,enh>::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  int res = -1;

  // quad pressure loading
  if (strcmp(argv[0],"pressure") == 0) {
    return param.addObject(2, this);
  }
  // a material parameter
  else if ((strstr(argv[0],"material") != 0) && (strcmp(argv[0],"materialState") != 0)) {

    if (argc < 3)
      return -1;

    int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= nip)
      return theMaterial[pointNum-1]->setParameter(&argv[2], argc-2, param);
    else 
      return -1;
  }

  // otherwise it could be just a forall material parameter
  else {

    int matRes = res;
    for (int i=0; i<nip; i++) {

      matRes =  theMaterial[i]->setParameter(argv, argc, param);

      if (matRes != -1)
        res = matRes;
    }
  }
  
  return res;
}

template <int NEN, bool enh>
int
LagrangeQuad<NEN,enh>::updateParameter(int parameterID, Information &info)
{
  int res = -1;
  int matRes = res;
  switch (parameterID) {
    case -1:
      return -1;

    case 1:
      for (int i = 0; i < nip; i++)
        matRes = theMaterial[i]->updateParameter(parameterID, info);

      if (matRes != -1) {
        res = matRes;
      }
      return res;
  
    case 2:
      pressure = info.theDouble;
      this->setPressureLoadAtNodes();        // update consistent nodal loads
      return 0;

    default: 
      /*          
      if (parameterID >= 100) { // material parameter
        int pointNum = parameterID/100;
        if (pointNum > 0 && pointNum <= 4)
          return theMaterial[pointNum-1]->updateParameter(parameterID-100*pointNum, info);
        else
          return -1;
      } else // unknown
      */
      return -1;
  }
}


template <int NEN, bool enh>
int
LagrangeQuad<NEN,enh>::activateParameter(int param)
{
  parameterID = param;
  return 0;
}


template <int NEN, bool enh>
int
LagrangeQuad<NEN,enh>::commitSensitivity(int gradNumber, int numGrads)
{

  static double u[NDM][NEN];

  for (int i = 0; i < NEN; i++) {
    u[0][i] = theNodes[i]->getDispSensitivity(1, gradNumber);
    u[1][i] = theNodes[i]->getDispSensitivity(2, gradNumber);
  }

  static Vector eps(3);

  int ret = 0;

  // Loop over the integration points
  for (int i = 0; i < NIP; i++) {

    // Determine Jacobian for this integration point
    this->shapeFunction(pts[i][0], pts[i][1]);

    // Interpolate strains
    //eps = B*u;
    //eps.addMatrixVector(0.0, B, u, 1.0);
    eps.Zero();
    for (int beta = 0; beta < 4; beta++) {
      eps(0) += shp[0][beta] * u[0][beta];
      eps(1) += shp[1][beta] * u[1][beta];
      eps(2) += shp[0][beta] * u[1][beta] + shp[1][beta] * u[0][beta];
    }

    // Set the material strain
    //ret += theMaterial[i]->setTrialStrain(eps);
    theMaterial[i]->commitSensitivity(eps, gradNumber, numGrads);
  }

  return 0;
}


template <int NEN, bool enh>
const Vector&
LagrangeQuad<NEN,enh>::getResistingForceSensitivity(int gradNumber)
{
  P.Zero();

  // Loop over the integration points
  for (int i = 0; i < NIP; i++) {

    // Determine Jacobian for this integration point
    double dvol = this->shapeFunction(pts[i][0], pts[i][1]);
    dvol *= (thickness * wts[i]);

    // Get material stress response
    const Vector& sigma = theMaterial[i]->getStressSensitivity(gradNumber, true);

    // Perform numerical integration on internal force
    //P = P + (B^ sigma) * intWt(i)*intWt(j) * detJ;
    //P.addMatrixTransposeVector(1.0, B, sigma, intWt(i)*intWt(j)*detJ);
    for (int alpha = 0, ia = 0; alpha < 4; alpha++, ia += 2) {
      P(ia)     += dvol * (shp[0][alpha] * sigma[0] + shp[1][alpha] * sigma[2]);
      P(ia + 1) += dvol * (shp[1][alpha] * sigma[1] + shp[0][alpha] * sigma[2]);
    }
  }
  return P;
}

#endif

template <int NEN, bool enh>
double
LagrangeQuad<NEN,enh>::shapeFunction(double xi, double eta)
{
  const Vector& nd1Crds = theNodes[0]->getCrds();
  const Vector& nd2Crds = theNodes[1]->getCrds();
  const Vector& nd3Crds = theNodes[2]->getCrds();
  const Vector& nd4Crds = theNodes[3]->getCrds();

  double oneMinuseta = 1.0 - eta;
  double onePluseta  = 1.0 + eta;
  double oneMinusxi  = 1.0 - xi;
  double onePlusxi   = 1.0 + xi;

  shp[2][0] = 0.25 * oneMinusxi * oneMinuseta; // N_1
  shp[2][1] = 0.25 * onePlusxi * oneMinuseta;  // N_2
  shp[2][2] = 0.25 * onePlusxi * onePluseta;   // N_3
  shp[2][3] = 0.25 * oneMinusxi * onePluseta;  // N_4

  double J[2][2];

  J[0][0] = 0.25 * (-nd1Crds(0) * oneMinuseta + nd2Crds(0) * oneMinuseta +
                    nd3Crds(0) * (onePluseta)-nd4Crds(0) * (onePluseta));

  J[0][1] = 0.25 * (-nd1Crds(0) * oneMinusxi - nd2Crds(0) * onePlusxi + nd3Crds(0) * onePlusxi +
                    nd4Crds(0) * oneMinusxi);

  J[1][0] = 0.25 * (-nd1Crds(1) * oneMinuseta + nd2Crds(1) * oneMinuseta + nd3Crds(1) * onePluseta -
                    nd4Crds(1) * onePluseta);

  J[1][1] = 0.25 * (-nd1Crds(1) * oneMinusxi - nd2Crds(1) * onePlusxi + nd3Crds(1) * onePlusxi +
                    nd4Crds(1) * oneMinusxi);

  double detJ        = J[0][0] * J[1][1] - J[0][1] * J[1][0];
  double oneOverdetJ = 1.0 / detJ;
  double L[2][2];

  // L = inv(J)
  L[0][0] =  J[1][1] * oneOverdetJ;
  L[1][0] = -J[0][1] * oneOverdetJ;
  L[0][1] = -J[1][0] * oneOverdetJ;
  L[1][1] =  J[0][0] * oneOverdetJ;

  double L00 = 0.25 * L[0][0];
  double L10 = 0.25 * L[1][0];
  double L01 = 0.25 * L[0][1];
  double L11 = 0.25 * L[1][1];

  double L00oneMinuseta = L00 * oneMinuseta;
  double L00onePluseta  = L00 * onePluseta;
  double L01oneMinusxi  = L01 * oneMinusxi;
  double L01onePlusxi   = L01 * onePlusxi;

  double L10oneMinuseta = L10 * oneMinuseta;
  double L10onePluseta  = L10 * onePluseta;
  double L11oneMinusxi  = L11 * oneMinusxi;
  double L11onePlusxi   = L11 * onePlusxi;

  // See Cook, Malkus, Plesha p. 169 for the derivation of these terms
  shp[0][0] = -L00oneMinuseta - L01oneMinusxi; // N_1,1
  shp[0][1] = L00oneMinuseta - L01onePlusxi;   // N_2,1
  shp[0][2] = L00onePluseta + L01onePlusxi;    // N_3,1
  shp[0][3] = -L00onePluseta + L01oneMinusxi;  // N_4,1

  shp[1][0] = -L10oneMinuseta - L11oneMinusxi; // N_1,2
  shp[1][1] = L10oneMinuseta - L11onePlusxi;   // N_2,2
  shp[1][2] = L10onePluseta + L11onePlusxi;    // N_3,2
  shp[1][3] = -L10onePluseta + L11oneMinusxi;  // N_4,2

  return detJ;
}


template <int NEN, bool enh>
void
LagrangeQuad<NEN,enh>::setPressureLoadAtNodes()
{
  pressureLoad.Zero();

  if (pressure == 0.0)
    return;

  const Vector& node1 = theNodes[0]->getCrds();
  const Vector& node2 = theNodes[1]->getCrds();
  const Vector& node3 = theNodes[2]->getCrds();
  const Vector& node4 = theNodes[3]->getCrds();

  double x1 = node1(0);
  double y1 = node1(1);
  double x2 = node2(0);
  double y2 = node2(1);
  double x3 = node3(0);
  double y3 = node3(1);
  double x4 = node4(0);
  double y4 = node4(1);

  double dx12 = x2 - x1;
  double dy12 = y2 - y1;
  double dx23 = x3 - x2;
  double dy23 = y3 - y2;
  double dx34 = x4 - x3;
  double dy34 = y4 - y3;
  double dx41 = x1 - x4;
  double dy41 = y1 - y4;

  double pressureOver2 = pressure / 2.0;

  // Contribution from side 12
  pressureLoad(0) += pressureOver2 * dy12;
  pressureLoad(2) += pressureOver2 * dy12;
  pressureLoad(1) += pressureOver2 * -dx12;
  pressureLoad(3) += pressureOver2 * -dx12;

  // Contribution from side 23
  pressureLoad(2) += pressureOver2 * dy23;
  pressureLoad(4) += pressureOver2 * dy23;
  pressureLoad(3) += pressureOver2 * -dx23;
  pressureLoad(5) += pressureOver2 * -dx23;

  // Contribution from side 34
  pressureLoad(4) += pressureOver2 * dy34;
  pressureLoad(6) += pressureOver2 * dy34;
  pressureLoad(5) += pressureOver2 * -dx34;
  pressureLoad(7) += pressureOver2 * -dx34;

  // Contribution from side 41
  pressureLoad(6) += pressureOver2 * dy41;
  pressureLoad(0) += pressureOver2 * dy41;
  pressureLoad(7) += pressureOver2 * -dx41;
  pressureLoad(1) += pressureOver2 * -dx41;

  //pressureLoad = pressureLoad*thickness;
}


} // namespace OpenSees
