/* *********************************************************************
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
// by Jinchi Lu and Zhaohui Yang (May 2004)
//
// 20NodeBrick element

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <Domain.h>
#include <Logging.h>
#include <Twenty_Node_Brick.h>
#include <shp3d.h>
#include <shp3dv.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>

double Twenty_Node_Brick::xl[3][20];

Matrix Twenty_Node_Brick::stiff(60, 60);
Vector Twenty_Node_Brick::resid(60);
Matrix Twenty_Node_Brick::mass(60, 60);
Matrix Twenty_Node_Brick::damp(60, 60);

const int Twenty_Node_Brick::nintu = 27;
double Twenty_Node_Brick::shgu[4][20][27];
double Twenty_Node_Brick::shlu[4][20][27];
double Twenty_Node_Brick::wu[27];
double Twenty_Node_Brick::dvolu[27];


Twenty_Node_Brick::Twenty_Node_Brick()
 : Element(0, ELE_TAG_Twenty_Node_Brick),
   connectedExternalNodes(20),
   applyLoad(0),
   load(0),
   Ki(0) //, kc(0), rho(0)
{
  for (int i = 0; i < NEN; i++) {
    nodePointers[i] = nullptr;
  }
  b[0] = b[1] = b[2] = 0.;

  // calculate local shape functions and derivatives
  compuLocalShapeFunction();
}


Twenty_Node_Brick::Twenty_Node_Brick(int tag, int node1, int node2, int node3, int node4, int node5,
                                     int node6, int node7, int node8, int node9, int node10,
                                     int node11, int node12, int node13, int node14, int node15,
                                     int node16, int node17, int node18, int node19, int node20,
                                     NDMaterial& theMaterial, double b1, double b2, double b3)
 : Element(tag, ELE_TAG_Twenty_Node_Brick),
   connectedExternalNodes(20),
   applyLoad(0),
   load(0),
   Ki(0) //, kc(bulk), rho(rhof)
{
  connectedExternalNodes(0)  = node1;
  connectedExternalNodes(1)  = node2;
  connectedExternalNodes(2)  = node3;
  connectedExternalNodes(3)  = node4;
  connectedExternalNodes(4)  = node5;
  connectedExternalNodes(5)  = node6;
  connectedExternalNodes(6)  = node7;
  connectedExternalNodes(7)  = node8;
  connectedExternalNodes(8)  = node9;
  connectedExternalNodes(9)  = node10;
  connectedExternalNodes(10) = node11;
  connectedExternalNodes(11) = node12;
  connectedExternalNodes(12) = node13;
  connectedExternalNodes(13) = node14;
  connectedExternalNodes(14) = node15;
  connectedExternalNodes(15) = node16;
  connectedExternalNodes(16) = node17;
  connectedExternalNodes(17) = node18;
  connectedExternalNodes(18) = node19;
  connectedExternalNodes(19) = node20;

  // Allocate arrays of pointers to NDMaterials
  materialPointers = new NDMaterial*[nintu];

  for (int i = 0; i < nintu; i++) {
    materialPointers[i] = theMaterial.getCopy("ThreeDimensional");
  }

  // Body forces
  b[0] = b1;
  b[1] = b2;
  b[2] = b3;
  //	printf("b %15.6e %15.6e %15.6e \n", b1, b2,b3);
  // calculate local shape functions and derivatives
  compuLocalShapeFunction();
}
//******************************************************************


//destructor
Twenty_Node_Brick::~Twenty_Node_Brick()
{
  for (int i = 0; i < nintu; i++) {
    if (materialPointers[i])
      delete materialPointers[i];
  }

  // Delete the array of pointers to NDMaterial pointer arrays
  if (materialPointers)
    delete[] materialPointers;

  for (int i = 0; i < NEN; i++) {
    nodePointers[i] = 0;
  }

  if (load != nullptr)
    delete load;

  if (Ki != nullptr)
    delete Ki;
}


//set domain
void
Twenty_Node_Brick::setDomain(Domain* theDomain)
{
  int i, dof;
  // Check Domain is not null - invoked when object removed from a domain
  if (theDomain == nullptr) {
    for (int i = 0; i < NEN; i++)
      nodePointers[i] = nullptr;
    return;
  }
  //node pointers
  for (int i = 0; i < NEN; i++) {
    nodePointers[i] = theDomain->getNode(connectedExternalNodes(i));
    if (nodePointers[i] == 0) {
      opserr << "FATAL ERROR Twenty_Node_Brick (" << this->getTag() << "): node not found in domain"
             << endln;
      return;
    }

    dof = nodePointers[i]->getNumberDOF();
    if (dof != 3) {
      opserr << "FATAL ERROR Twenty_Node_Brick (" << this->getTag()
             << "): has wrong number of DOFs at its nodes" << endln;
      return;
    }
  }
  this->DomainComponent::setDomain(theDomain);
}



int
Twenty_Node_Brick::getNumExternalNodes() const
{
  return NEN;
}



const ID&
Twenty_Node_Brick::getExternalNodes()
{
  return connectedExternalNodes;
}


Node**
Twenty_Node_Brick::getNodePtrs()
{
  return nodePointers;
}


int
Twenty_Node_Brick::getNumDOF()
{
  return NEN*3;
}


int
Twenty_Node_Brick::commitState()
{
  int success = 0;

  // call element commitState to do any base class stuff
  if ((success = this->Element::commitState()) != 0) {
    opserr << "Twenty_Node_Brick::commitState () - failed in base class";
  }

  for (int i = 0; i < nintu; i++)
    success += materialPointers[i]->commitState();

  return success;
}


//revert to last commit
int
Twenty_Node_Brick::revertToLastCommit()
{
  int success = 0;

  for (int i = 0; i < nintu; i++)
    success += materialPointers[i]->revertToLastCommit();

  return success;
}


int
Twenty_Node_Brick::revertToStart()
{
  int success = 0;

  for (int i = 0; i < nintu; i++)
    success += materialPointers[i]->revertToStart();

  return success;
}


void
Twenty_Node_Brick::Print(OPS_Stream& s, int flag)
{

  if (flag == 2) {
    s << "#20NodeBrick\n";

    int i;
    const int numNodes = 20;
    const int nstress  = 6;

    for (i = 0; i < numNodes; i++) {
      const Vector& nodeCrd  = nodePointers[i]->getCrds();
      const Vector& nodeDisp = nodePointers[i]->getDisp();
      s << "#NODE " << nodeCrd(0) << " " << nodeCrd(1) << " " << nodeCrd(2) << " " << nodeDisp(0)
        << " " << nodeDisp(1) << " " << nodeDisp(2) << endln;
    }

    // spit out the section location & invoke print on the scetion
    const int numMaterials = nintu;

    static Vector avgStress(7);
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
    for (i = 0; i < 7; i++)
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
    s << endln;
    s << "20NodeBrick Twenty_Node_Brick \n";
    s << "Element Number: " << this->getTag() << endln;
    s << "Node 1 : " << connectedExternalNodes(0) << endln;
    s << "Node 2 : " << connectedExternalNodes(1) << endln;
    s << "Node 3 : " << connectedExternalNodes(2) << endln;
    s << "Node 4 : " << connectedExternalNodes(3) << endln;
    s << "Node 5 : " << connectedExternalNodes(4) << endln;
    s << "Node 6 : " << connectedExternalNodes(5) << endln;
    s << "Node 7 : " << connectedExternalNodes(6) << endln;
    s << "Node 8 : " << connectedExternalNodes(7) << endln;
    s << "Node 9 : " << connectedExternalNodes(8) << endln;
    s << "Node 10 : " << connectedExternalNodes(9) << endln;
    s << "Node 11 : " << connectedExternalNodes(10) << endln;
    s << "Node 12 : " << connectedExternalNodes(11) << endln;
    s << "Node 13 : " << connectedExternalNodes(12) << endln;
    s << "Node 14 : " << connectedExternalNodes(13) << endln;
    s << "Node 15 : " << connectedExternalNodes(14) << endln;
    s << "Node 16 : " << connectedExternalNodes(15) << endln;
    s << "Node 17 : " << connectedExternalNodes(16) << endln;
    s << "Node 18 : " << connectedExternalNodes(17) << endln;
    s << "Node 19 : " << connectedExternalNodes(18) << endln;
    s << "Node 20 : " << connectedExternalNodes(19) << endln;

    s << "Material Information : \n ";
    materialPointers[0]->Print(s, flag);

    s << endln;
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"20NodeBrick\", ";
    s << "\"nodes\": [" << connectedExternalNodes(0) << ", ";
    for (int i = 1; i < 19; i++)
      s << connectedExternalNodes(i) << ", ";
    s << connectedExternalNodes(19) << "], ";
    s << "\"bodyForces\": [" << b[0] << ", " << b[1] << ", " << b[2] << "], ";
    s << "\"material\": \"" << materialPointers[0]->getTag() << "\"}";
  }
}

int
Twenty_Node_Brick::update()
{
  int i, j, k, k1;
  static double u[3][20];
  static double xsj;
  static Matrix B(6, 3);
  double volume = 0.;

  for (i = 0; i < NEN; i++) {
    const Vector& disp = nodePointers[i]->getTrialDisp();
    u[0][i]            = disp(0);
    u[1][i]            = disp(1);
    u[2][i]            = disp(2);
  }

  static Vector eps(6);

  int ret = 0;

  //compute basis vectors and local nodal coordinates
  computeBasis();

  for (i = 0; i < nintu; i++) {
    // compute Jacobian and global shape functions
    Jacobian3d(i, xsj, 0);
    //volume element to also be saved
    dvolu[i] = wu[i] * xsj;
    volume += dvolu[i];
  }

  // Loop over the integration points
  for (i = 0; i < nintu; i++) {

    // Interpolate strains
    //eps = B*u;
    //eps.addMatrixVector(0.0, B, u, 1.0);
    eps.Zero();
    for (j = 0; j < NEN; j++) {

      B(0, 0) = shgu[0][j][i];
      B(0, 1) = 0.;
      B(0, 2) = 0.;
      B(1, 0) = 0.;
      B(1, 1) = shgu[1][j][i];
      B(1, 2) = 0.;
      B(2, 0) = 0.;
      B(2, 1) = 0.;
      B(2, 2) = shgu[2][j][i];
      B(3, 0) = shgu[1][j][i];
      B(3, 1) = shgu[0][j][i];
      B(3, 2) = 0.;
      B(4, 0) = 0.;
      B(4, 1) = shgu[2][j][i];
      B(4, 2) = shgu[1][j][i];
      B(5, 0) = shgu[2][j][i];
      B(5, 1) = 0.;
      B(5, 2) = shgu[0][j][i];

      //BJ = computeB( j, shp ) ;

      //nodal displacements
      const Vector& ul = nodePointers[j]->getTrialDisp();
      Vector ul3(3);
      ul3(0) = ul(0);
      ul3(1) = ul(1);
      ul3(2) = ul(2);
      //compute the strain
      //strain += (BJ*ul) ;
      eps.addMatrixVector(1.0, B, ul3, 1.0);

      /* for( k = 0; k < 6; k++) {
			for( k1 = 0; k1 < 3; k1++) {
			  eps[k] += BJ(k, k1)*u[k1][j];
			}
		}	*/
    }

    // Set the material strain
    ret += materialPointers[i]->setTrialStrain(eps);
  }

  return ret;
}


//return tangent stiffness matrix

const Matrix&
Twenty_Node_Brick::getTangentStiff()
{
  return getStiff(1);
}


const Matrix&
Twenty_Node_Brick::getInitialStiff()
{
  return getStiff(0);
}


// compute stiffness matrix

const Matrix&
Twenty_Node_Brick::getStiff(int flag)
{
  if (flag != 0 && flag != 1) {
    opserr << "FATAL Twenty_Node_Brick::getStiff() - illegal use\n";
  }


  if (flag == 0 && Ki != nullptr)
    return *Ki;


  int i, j;
  double xsj; // determinant jacaobian matrix

  double volume = 0.;

  //-------------------------------------------------------

  int j3, j3m1, j3m2, ik, ib, jk, jb;

  static Matrix B(6, NEN * 3);
  static Matrix BTDB(NEN * 3, NEN * 3);
  static Matrix D(6, 6);

  B.Zero();
  BTDB.Zero();
  stiff.Zero();

  // compute basis vectors and local nodal coordinates

  computeBasis();

  for (int i = 0; i < nintu; i++) {

    // compute Jacobian and global shape functions
    Jacobian3d(i, xsj, 0);

    // volume element to also be saved
    dvolu[i] = wu[i] * xsj;
    volume += dvolu[i];
  }

  // Loop over the integration points

  for (int i = 0; i < nintu; i++) {

    // Get the material tangent

    if (flag == 0)
      D = materialPointers[i]->getInitialTangent();
    else
      D = materialPointers[i]->getTangent();

    //const Matrix &D = materialPointers[i]->getTangent();


    for (int j = 0; j < NEN; j++) {

      j3 = 3 * j + 2;
      j3m1 = j3 - 1;
      j3m2 = j3 - 2;

      B(0, j3m2) = shgu[0][j][i];
      B(0, j3m1) = 0.;
      B(0, j3)   = 0.;

      B(1, j3m2) = 0.;
      B(1, j3m1) = shgu[1][j][i];
      B(1, j3)   = 0.;

      B(2, j3m2) = 0.;
      B(2, j3m1) = 0.;
      B(2, j3) = shgu[2][j][i];

      B(3, j3m2) = shgu[1][j][i];
      B(3, j3m1) = shgu[0][j][i];
      B(3, j3) = 0.;

      B(4, j3m2) = 0.;
      B(4, j3m1) = shgu[2][j][i];
      B(4, j3) = shgu[1][j][i];

      B(5, j3m2) = shgu[2][j][i];
      B(5, j3m1) = 0.;
      B(5, j3) = shgu[0][j][i];
    }

    // Perform numerical integration

    //K = K + (B^ D * B) * intWt(i) * detJ;

    BTDB.addMatrixTripleProduct(1.0, B, D, dvolu[i]);
  }

  for (int i = 0; i < 60; i++)
    for (int j = 0; j < 60; j++)
      stiff(i, j) = BTDB(i, j);


  if (flag == 1) {
    return stiff;
  }

  Ki = new Matrix(stiff);

  return *Ki;
}


const Matrix&
Twenty_Node_Brick::getMass()
{
  int tangFlag = 1;

  formInertiaTerms(tangFlag);

  return mass;
}


//return mass matrix

const Matrix&
Twenty_Node_Brick::getDamp()
{
  int tangFlag = 1;
  formDampingTerms(tangFlag);
  return damp;
}


void
Twenty_Node_Brick::formDampingTerms(int tangFlag)
{

  static double xsj; // determinant jacaobian matrix

  int i, j, k, m, ik, jk;
  double volume = 0.;

  // zero damp

  damp.Zero();

  if (betaK != 0.0)
    damp.addMatrix(1.0, this->getTangentStiff(), betaK);

  if (betaK0 != 0.0)
    damp.addMatrix(1.0, this->getInitialStiff(), betaK0);

  if (betaKc != 0.0)
    damp.addMatrix(1.0, *Kc, betaKc);


  if (alphaM != 0.0) {
    this->getMass();

    for (i = 0; i < 60; i++)
      for (j = 0; j < 60; j++)
        damp(i, j) += mass(i, j) * alphaM;
  }

  return;
}


void
Twenty_Node_Brick::zeroLoad()
{
  if (load != 0) {
    load->Zero();
  }

  applyLoad   = 0;

  appliedB[0] = 0.0;
  appliedB[1] = 0.0;
  appliedB[2] = 0.0;

  return;
}


int

Twenty_Node_Brick::addLoad(ElementalLoad* theLoad, double loadFactor)

{
  // Added option for applying body forces in load pattern: C.McGann, U.Washington
  int type;
  const Vector& data = theLoad->getData(type, loadFactor);

  if (type == LOAD_TAG_BrickSelfWeight) {
    applyLoad = 1;
    appliedB[0] += loadFactor * b[0];
    appliedB[1] += loadFactor * b[1];
    appliedB[2] += loadFactor * b[2];
    return 0;
  }
  
  else if (type == LOAD_TAG_SelfWeight) {
    // added compatibility with selfWeight class implemented for all continuum elements, C.McGann, U.W.
    applyLoad = 1;
    appliedB[0] += loadFactor * data(0) * b[0];
    appliedB[1] += loadFactor * data(1) * b[1];
    appliedB[2] += loadFactor * data(2) * b[2];
    return 0;

  }
  
  else {
    opserr << "Twenty_Node_Brick::addLoad - load type unknown for truss with tag: "
           << this->getTag() << endln;
    return -1;
  }

  return -1;
}


int
Twenty_Node_Brick::addInertiaLoadToUnbalance(const Vector& accel)
{

  static Vector ra(60);

  int i, j, ik;

  ra.Zero();


  for (int i = 0; i < NEN; i++) {

    const Vector& Raccel = nodePointers[i]->getRV(accel);

    if (3 != Raccel.Size()) {

      opserr << "Twenty_Node_Brick::addInertiaLoadToUnbalance matrix and vector sizes are "
                "incompatible\n";

      return -1;
    }

    ra[i * 3]     = Raccel(0);
    ra[i * 3 + 1] = Raccel(1);
    ra[i * 3 + 2] = Raccel(2);
  }


  // Compute mass matrix

  int tangFlag = 1;

  formInertiaTerms(tangFlag);

  // create the load vector if one does not exist

  if (load == 0)

    load = new Vector(60);


  // add -M * RV(accel) to the load vector

  load->addMatrixVector(1.0, mass, ra, -1.0);

  return 0;
}


//get residual

const Vector&
Twenty_Node_Brick::getResistingForce()

{

  int i, j, jk, k, k1;

  double xsj;
  static Matrix B(6, 3);
  double volume = 0.;
  resid.Zero();


  // compute basis vectors and local nodal coordinates

  computeBasis();

  // gauss loop to compute and save shape functions

  for (int i = 0; i < nintu; i++) {

    // compute Jacobian and global shape functions

    Jacobian3d(i, xsj, 0);

    // volume element to also be saved

    dvolu[i] = wu[i] * xsj;

    volume += dvolu[i];

  }


  // Loop over the integration points

  for (i = 0; i < nintu; i++) {

    // Get material stress response

    const Vector& sigma = materialPointers[i]->getStress();


    // Perform numerical integration on internal force

    //P = P + (B^ sigma) * intWt(i)*intWt(j) * detJ;

    //P.addMatrixTransposeVector(1.0, B, sigma, intWt(i)*intWt(j)*detJ);

    for (int j = 0; j < NEN; j++) {


      B(0, 0) = shgu[0][j][i];

      B(0, 1) = 0.;
      B(0, 2) = 0.;
      B(1, 0) = 0.;

      B(1, 1) = shgu[1][j][i];
      B(1, 2) = 0.;
      B(2, 0) = 0.;
      B(2, 1) = 0.;

      B(2, 2) = shgu[2][j][i];
      B(3, 0) = shgu[1][j][i];
      B(3, 1) = shgu[0][j][i];

      B(3, 2) = 0.;
      B(4, 0) = 0.;
      B(4, 1) = shgu[2][j][i];
      B(4, 2) = shgu[1][j][i];
      B(5, 0) = shgu[2][j][i];

      B(5, 1) = 0.;
      B(5, 2) = shgu[0][j][i];


      for (int k = 0; k < 3; k++) {
        for (k1 = 0; k1 < 6; k1++)
          resid(j * 3 + k) += dvolu[i] * (B(k1, k) * sigma(k1));
      }

      // Subtract equiv. body forces from the nodes

      //P = P - (N^ b) * intWt(i)*intWt(j) * detJ;

      //P.addMatrixTransposeVector(1.0, N, b, -intWt(i)*intWt(j)*detJ);

      double r = mixtureRho(i);

      if (applyLoad == 0) {
        resid(j * 3)     -= dvolu[i] * (shgu[3][j][i] * r * b[0]);
        resid(j * 3 + 1) -= dvolu[i] * (shgu[3][j][i] * r * b[1]);
        resid(j * 3 + 2) -= dvolu[i] * (shgu[3][j][i] * r * b[2]);

      } else {
        resid(j * 3)     -= dvolu[i] * (shgu[3][j][i] * r * appliedB[0]);
        resid(j * 3 + 1) -= dvolu[i] * (shgu[3][j][i] * r * appliedB[1]);
        resid(j * 3 + 2) -= dvolu[i] * (shgu[3][j][i] * r * appliedB[2]);
      }
    }
  }


  // Subtract other external nodal loads ... P_res = P_int - P_ext

  if (load != 0)
    resid -= *load;

  return resid;
}


//get residual with inertia terms

const Vector&
Twenty_Node_Brick::getResistingForceIncInertia()
{
  static Vector res(60);

  int i, j, ik;

  static double a[60];

  for (i = 0; i < NEN; i++) {

    const Vector& accel = nodePointers[i]->getTrialAccel();

    if (3 != accel.Size()) {

      opserr << "Twenty_Node_Brick::getResistingForceIncInertia matrix and vector sizes are "
                "incompatible\n";

      exit(-1);
    }

    a[i * 3]     = accel(0);
    a[i * 3 + 1] = accel(1);
    a[i * 3 + 2] = accel(2);
  }

  // Compute the current resisting force

  this->getResistingForce();

  // Compute the mass matrix

  this->getMass();


  for (i = 0; i < 60; i++) {
    for (j = 0; j < 60; j++) {
      resid(i) += mass(i, j) * a[j];
    }
  }


  for (i = 0; i < NEN; i++) {

    const Vector& vel = nodePointers[i]->getTrialVel();

    if (3 != vel.Size()) {

      opserr << "Twenty_Node_Brick::getResistingForceIncInertia matrix and vector sizes are "
                "incompatible\n";

      exit(-1);
    }

    a[i * 3] = vel(0);
    a[i * 3 + 1] = vel(1);
    a[i * 3 + 2] = vel(2);
  }


  this->getDamp();


  for (i = 0; i < 60; i++) {

    for (j = 0; j < 60; j++) {

      resid(i) += damp(i, j) * a[j];
    }
  }


  res = resid;

  return res;
}


void
Twenty_Node_Brick::formInertiaTerms(int tangFlag)

{

  static double xsj; // determinant jacaobian matrix

  int i, j, k, ik, m, jk;

  double Nrho;


  // zero mass

  mass.Zero();


  // compute basis vectors and local nodal coordinates

  computeBasis();

  //gauss loop to compute and save shape functions

  for (i = 0; i < nintu; i++) {

    // compute Jacobian and global shape functions

    Jacobian3d(i, xsj, 0);

    // volume element to also be saved

    dvolu[i] = wu[i] * xsj;
  }


  // Compute consistent mass matrix

  for (i = 0; i < NEN; i++) {

    for (j = 0; j < NEN; j++) {

      for (m = 0; m < nintu; m++) {

        Nrho = dvolu[m] * mixtureRho(m) * shgu[3][i][m] * shgu[3][j][m];

        for (k = 0; k < 3; k++) {

          mass(i * 3 + k, j * 3 + k) += Nrho;
        }
      }
    }
  }
}


double
Twenty_Node_Brick::mixtureRho(int i)
{
  double rhoi = materialPointers[i]->getRho();

  //e = 0.7;  //theMaterial[i]->getVoidRatio();
  //n = e / (1.0 + e);
  //return n * rho + (1.0-n) * rhoi;

  return rhoi;
}


void
Twenty_Node_Brick::computeBasis()
{
  for (int i = 0; i < NEN; i++) {

    const Vector& coorI = nodePointers[i]->getCrds();

    xl[0][i] = coorI(0);

    xl[1][i] = coorI(1);

    xl[2][i] = coorI(2);
  }
}


//**********************************************************************


int
Twenty_Node_Brick::sendSelf(int commitTag, Channel& theChannel)

{

  int res = 0;

  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data

  int dataTag = this->getDbTag();


  // Quad packs its data into a Vector and sends this to theChannel
  // along with its dbTag and the commitTag passed in the arguments

  int matDbTag;


  static ID idData(75);


  idData(74) = this->getTag();


  int i;

  for (i = 0; i < nintu; i++) {

    idData(i) = materialPointers[i]->getClassTag();

    matDbTag = materialPointers[i]->getDbTag();

    // NOTE: we do have to ensure that the material has a database

    // tag if we are sending to a database channel.

    if (matDbTag == 0) {
      matDbTag = theChannel.getDbTag();
      if (matDbTag != 0)
        materialPointers[i]->setDbTag(matDbTag);
    }

    idData(i + nintu) = matDbTag;
  }

  for (i = 0; i < 20; i++)
    idData(54 + i) = connectedExternalNodes(i);

  res += theChannel.sendID(dataTag, commitTag, idData);

  if (res < 0) {
    opserr << "WARNING Twenty_Node_Brick::sendSelf() - " << this->getTag()
           << " failed to send ID\n";
    return res;
  }


  // Finally, this element asks its material objects to send themselves

  for (i = 0; i < nintu; i++) {

    res += materialPointers[i]->sendSelf(commitTag, theChannel);

    if (res < 0) {

      opserr << "WARNING Twenty_Node_Brick::sendSelf() - " << this->getTag()
             << " failed to send its Material\n";

      return res;
    }
  }


  return res;
}


int
Twenty_Node_Brick::recvSelf(int commitTag,
                            Channel& theChannel,
                            FEM_ObjectBroker& theBroker)

{

  int res = 0;

  int dataTag = this->getDbTag();

  static ID idData(75);

  //  now receives the tags of its 20 external nodes

  res += theChannel.recvID(dataTag, commitTag, idData);

  if (res < 0) {

    opserr << "WARNING Twenty_Node_Brick::recvSelf() - " << this->getTag()
           << " failed to receive ID\n";

    return res;
  }


  this->setTag(idData(74));


  int i;

  for (int i = 0; i < 20; i++)
    connectedExternalNodes(i) = idData(54 + i);


  if (materialPointers[0] == 0) {
    for (i = 0; i < nintu; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i + nintu);

      // Allocate new material with the sent class tag

      materialPointers[i] = theBroker.getNewNDMaterial(matClassTag);

      if (materialPointers[i] == 0) {

        opserr
            << "Twenty_Node_Brick::recvSelf() - Broker could not create NDMaterial of class type "
            << matClassTag << endln;

        return -1;
      }

      // Now receive materials into the newly allocated space

      materialPointers[i]->setDbTag(matDbTag);

      res += materialPointers[i]->recvSelf(commitTag, theChannel, theBroker);

      if (res < 0) {

        opserr << "Twenty_Node_Brick::recvSelf() - material " << i << "failed to recv itself\n";

        return res;
      }
    }

  }

  // materials exist , ensure materials of correct type and recvSelf on them

  else {

    for (i = 0; i < nintu; i++) {

      int matClassTag = idData(i);

      int matDbTag = idData(i + nintu);

      // Check that material is of the right type; if not,

      // delete it and create a new one of the right type

      if (materialPointers[i]->getClassTag() != matClassTag) {

        delete materialPointers[i];

        materialPointers[i] = theBroker.getNewNDMaterial(matClassTag);

        if (materialPointers[i] == 0) {

          opserr
              << "Twenty_Node_Brick::recvSelf() - Broker could not create NDMaterial of class type "
              <<

              matClassTag << endln;

          exit(-1);
        }

        materialPointers[i]->setDbTag(matDbTag);
      }

      // Receive the material


      res += materialPointers[i]->recvSelf(commitTag, theChannel, theBroker);

      if (res < 0) {

        opserr << "Twenty_Node_Brick::recvSelf() - material " << i << "failed to recv itself\n";

        return res;
      }
    }
  }

  return res;
}


Response*
Twenty_Node_Brick::setResponse(const char** argv, int argc, OPS_Stream& output)
{

  Response* theResponse = 0;

  char outputData[32];

  output.tag("ElementOutput");
  output.attr("eleType", "Twenty_Node_Brick");
  output.attr("eleTag", this->getTag());
  for (int i = 1; i <= 20; i++) {
    sprintf(outputData, "node%d", i);
    output.attr(outputData, connectedExternalNodes[i - 1]);
  }

  if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "forces") == 0) {

    for (int i = 1; i <= 20; i++)
      for (int j = 1; j <= 3; j++) {
        sprintf(outputData, "P%d_%d", j, i);
        output.tag("ResponseType", outputData);
      }

    theResponse = new ElementResponse(this, 1, resid);

  } else if (strcmp(argv[0], "material") == 0 || strcmp(argv[0], "integrPoint") == 0) {
    int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= nintu) {

      output.tag("GaussPoint");
      output.attr("number", pointNum);

      theResponse = materialPointers[pointNum - 1]->setResponse(&argv[2], argc - 2, output);

      output.endTag(); // GaussPoint
    }
  } else if (strcmp(argv[0], "stresses") == 0) {

    for (int i = 0; i < 27; i++) {
      output.tag("GaussPoint");
      output.attr("number", i + 1);
      output.tag("NdMaterialOutput");
      output.attr("classType", materialPointers[i]->getClassTag());
      output.attr("tag", materialPointers[i]->getTag());

      output.tag("ResponseType", "sigma11");
      output.tag("ResponseType", "sigma22");
      output.tag("ResponseType", "sigma33");
      output.tag("ResponseType", "sigma12");
      output.tag("ResponseType", "sigma23");
      output.tag("ResponseType", "sigma13");

      output.endTag(); // NdMaterialOutput
      output.endTag(); // GaussPoint
    }
    theResponse = new ElementResponse(this, 5, Vector(162));
  }

  output.endTag(); // ElementOutput
  return theResponse;
}


int

Twenty_Node_Brick::getResponse(int responseID, Information& eleInfo)

{

  static Vector stresses(162);


  if (responseID == 1)

    return eleInfo.setVector(this->getResistingForce());


  else if (responseID == 2)

    return eleInfo.setMatrix(this->getTangentStiff());


  else if (responseID == 3)

    return eleInfo.setMatrix(this->getMass());


  else if (responseID == 4)

    return eleInfo.setMatrix(this->getDamp());


  else if (responseID == 5) {


    // Loop over the integration points

    int cnt = 0;

    for (int i = 0; i < nintu; i++) {


      // Get material stress response

      const Vector& sigma = materialPointers[i]->getStress();

      stresses(cnt++) = sigma(0);

      stresses(cnt++) = sigma(1);

      stresses(cnt++) = sigma(2);

      stresses(cnt++) = sigma(3);

      stresses(cnt++) = sigma(4);

      stresses(cnt++) = sigma(5);
    }

    return eleInfo.setVector(stresses);


  }

  else


    return -1;
}


// calculate local shape functions

void

Twenty_Node_Brick::compuLocalShapeFunction()
{


  int i, k, j;

  static double shl[4][20][27], w[27];

  // solid phase

  brcshl(shl, w, nintu, NEN);

  for (k = 0; k < nintu; k++) {

    wu[k] = w[k];

    for (j = 0; j < NEN; j++)

      for (i = 0; i < 4; i++)

        shlu[i][j][k] = shl[i][j][k];
  }
}


void

Twenty_Node_Brick::Jacobian3d(int gaussPoint, double& xsj, int mode)

{

  int i, j, k, nint, nen;
  double rxsj, c1, c2, c3;
  static double xs[3][3];
  static double ad[3][3];
  static double shp[4][NEN];


  if (mode == 0) { // solid
    nint = nintu;
    nen  = NEN;
  }
  else {
    opserr << "Twenty_Node_Brick::Jacobian3d - illegal mode: " << mode << "\n";
  }


  for (j = 0; j < NEN; j++) {

    for (i = 0; i < 4; i++) {

      if (mode == 0)

        shp[i][j] = shlu[i][j][gaussPoint];

      else {

        opserr << "Twenty_Node_Brick::Jacobian3d - illegal mode: " << mode << "\n";

      } //end if
    }
  }


  //Compute jacobian transformation


  for (j = 0; j < 3; j++) {

    for (k = 0; k < 3; k++) {

      xs[j][k] = 0;

      for (i = 0; i < nen; i++) {

        xs[j][k] += xl[j][i] * shp[k][i];
      }
    }
  }


  // Compute adjoint to jacobian

  ad[0][0] = xs[1][1] * xs[2][2] - xs[1][2] * xs[2][1];
  ad[0][1] = xs[2][1] * xs[0][2] - xs[2][2] * xs[0][1];
  ad[0][2] = xs[0][1] * xs[1][2] - xs[0][2] * xs[1][1];

  ad[1][0] = xs[1][2] * xs[2][0] - xs[1][0] * xs[2][2];
  ad[1][1] = xs[2][2] * xs[0][0] - xs[2][0] * xs[0][2];
  ad[1][2] = xs[0][2] * xs[1][0] - xs[0][0] * xs[1][2];

  ad[2][0] = xs[1][0] * xs[2][1] - xs[1][1] * xs[2][0];
  ad[2][1] = xs[2][0] * xs[0][1] - xs[2][1] * xs[0][0];
  ad[2][2] = xs[0][0] * xs[1][1] - xs[0][1] * xs[1][0];


  // Compute determinant of jacobian


  xsj = xs[0][0] * ad[0][0] + xs[0][1] * ad[1][0] + xs[0][2] * ad[2][0];

  if (xsj <= 0) {

    opserr << "Twenty_Node_Brick::Jacobian3d - Non-positive Jacobian: " << xsj << "\n";

    for (i = 0; i < nen; i++) {

      printf("%5d %15.6e %15.6e %15.6e %15.6e\n", i,

             shp[0][i], shp[1][i], shp[2][i], shp[3][i]);
    }


    exit(-1);
  }


  rxsj = 1.0 / xsj;


  //Compute jacobian inverse


  for (j = 0; j < 3; j++) {
    for (i = 0; i < 3; i++)

      xs[i][j] = ad[i][j] * rxsj;


  }


  //Compute derivatives with repect to global coords.


  for (k = 0; k < nen; k++) {

    c1 = shp[0][k] * xs[0][0] + shp[1][k] * xs[1][0] + shp[2][k] * xs[2][0];

    c2 = shp[0][k] * xs[0][1] + shp[1][k] * xs[1][1] + shp[2][k] * xs[2][1];

    c3 = shp[0][k] * xs[0][2] + shp[1][k] * xs[1][2] + shp[2][k] * xs[2][2];


    shp[0][k] = c1;

    shp[1][k] = c2;

    shp[2][k] = c3;


  } //end for k


  for (j = 0; j < nen; j++) {

    for (i = 0; i < 4; i++) {

      if (mode == 0)

        shgu[i][j][gaussPoint] = shp[i][j];

      else {

        opserr << "Twenty_Node_Brick::Jacobian3d - illegal mode: " << mode << "\n";

        exit(-1);

      } //end if
    }
  }
}
