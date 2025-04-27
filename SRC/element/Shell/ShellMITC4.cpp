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
// Bathe MITC 4 four node shell element with membrane and drill
//
// Original implementation: Ed "C++" Love
// Reimplementation: Leopoldo Tesser, Diego A. Talledo, Veronique Le Corvec
//
// [1] https://opensees.berkeley.edu/OpenSees/doc/plate.pdf
//
// [2]  Dvorkin,Bathe, A continuum mechanics based four node shell
//      element for general nonlinear analysis,
//      Eng.Comput.,1,77-88,1984
//
#include <assert.h>
#include <math.h>

#include <ID.h>
#include <Vector.h>
#include <VectorND.h>
#include <Vector3D.h>
#include <Matrix.h>
#include <MatrixND.h>
#include <Matrix3D.h>

#include <Element.h>
#include <Node.h>
#include <SectionForceDeformation.h>
#include <Domain.h>
#include <ShellMITC4.h>
#include <Renderer.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <Parameter.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>

#define min(a, b) ((a) < (b) ? (a) : (b))

using namespace OpenSees;

// static data
Matrix ShellMITC4::stiff(24, 24);
Vector ShellMITC4::resid(24);
Matrix ShellMITC4::mass(24, 24);


ShellMITC4::ShellMITC4()
    : Element(0, ELE_TAG_ShellMITC4), 
      connectedExternalNodes(4), load(0), Ki(0),
      doUpdateBasis(false)
{

  for (int i = 0; i < 4; i++)
    materialPointers[i] = nullptr;

  applyLoad = 0;

  appliedB[0] = 0.0;
  appliedB[1] = 0.0;
  appliedB[2] = 0.0;
}


ShellMITC4::ShellMITC4(int tag, int node1, int node2, int node3, int node4,
                       SectionForceDeformation &theMaterial, bool UpdateBasis)
    : Element(tag, ELE_TAG_ShellMITC4), 
      connectedExternalNodes(NEN), 
      load(0),
      Ki(0), 
      doUpdateBasis(UpdateBasis)
{

  connectedExternalNodes(0) = node1;
  connectedExternalNodes(1) = node2;
  connectedExternalNodes(2) = node3;
  connectedExternalNodes(3) = node4;

  for (int i = 0; i < nip; i++) {
    materialPointers[i] = theMaterial.getCopy();
    if (materialPointers[i] == nullptr) {
      opserr << "ShellMITC4::constructor - failed to get a material of type: "
                "ShellSection\n";
    }
  }

  applyLoad = 0;

  appliedB[0] = 0.0;
  appliedB[1] = 0.0;
  appliedB[2] = 0.0;
}


ShellMITC4::~ShellMITC4()
{
  for (int i = 0; i < 4; i++) {
    delete materialPointers[i];
    materialPointers[i] = nullptr;
    theNodes[i]     = nullptr;
  }

  if (load != nullptr)
    delete load;

  if (Ki != nullptr)
    delete Ki;
}


void ShellMITC4::setDomain(Domain *theDomain)
{
  Vector3D eig;
  Matrix3D ddMembrane;

  // node pointers
  for (int i = 0; i < 4; i++) {
    theNodes[i] = theDomain->getNode(connectedExternalNodes(i));
    if (theNodes[i] == 0) {
      opserr << "ShellMITC4::setDomain - no node " << connectedExternalNodes(i);
      opserr << " exists in the model\n";
    }
    const Vector &nodeDisp = theNodes[i]->getTrialDisp();
    assert(nodeDisp.Size() == 6);

    init_disp[i][0] = nodeDisp(0);
    init_disp[i][1] = nodeDisp(1);
    init_disp[i][2] = nodeDisp(2);
    init_disp[i][3] = nodeDisp(3);
    init_disp[i][4] = nodeDisp(4);
    init_disp[i][5] = nodeDisp(5);
  }

  // compute drilling stiffness penalty parameter
  const Matrix &dd = materialPointers[0]->getInitialTangent();

  // assemble ddMembrane ;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++)
      ddMembrane(i, j) = dd(i, j);
  }

  // eigenvalues of ddMembrane
  ddMembrane.symeig(eig);

  // set ktt
  Ktt = min(eig(2), min(eig(0), eig(1)));

  // basis vectors and local coordinates
  computeBasis();

  this->DomainComponent::setDomain(theDomain);
}


int ShellMITC4::getNumExternalNodes() const 
{
  return NEN; 
}


const ID &ShellMITC4::getExternalNodes() 
{
  return connectedExternalNodes; 
}


Node **
ShellMITC4::getNodePtrs() 
{ 
  return &theNodes[0]; 
}

// return number of dofs
int ShellMITC4::getNumDOF() 
{
  int sum = 0;

  for (const Node* node: theNodes)
    sum += node->getNumberDOF();

  return sum;
}


int ShellMITC4::commitState()
{
  int success = 0;

  // call element commitState to do any base class stuff
  if ((success = this->Element::commitState()) != 0) {
    opserr << "ShellMITC4::commitState - failed in base class";
  }

  for (auto material : materialPointers)
    success += material->commitState();

  return success;
}


int
ShellMITC4::revertToLastCommit()
{
  int success = 0;

  for (auto material : materialPointers)
    success += material->revertToLastCommit();

  return success;
}

// revert to start
int ShellMITC4::revertToStart()
{
  int success = 0;

  for (auto material : materialPointers)
    success += material->revertToStart();

  return success;
}

// print out element data
void ShellMITC4::Print(OPS_Stream &s, int flag)
{
  const ID& node_tags = this->getExternalNodes();

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_ELEM_INDENT << "{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"" << this->getClassType() << "\", ";
    s << "\"nodes\": [" 
      << connectedExternalNodes(0) << ", "
      << connectedExternalNodes(1) << ", ";
    s << connectedExternalNodes(2) << ", " 
      << connectedExternalNodes(3)
      << "], ";

    s << "\"sections\": [" << materialPointers[0]->getTag() << "]";
    s << "}";
    return;
  }

  if (flag == -1) {
    int eleTag = this->getTag();
    s << "EL_ShellMITC4\t" << eleTag << "\t";
    s << eleTag << "\t" << 1;
    s << "\t" << connectedExternalNodes(0) << "\t" << connectedExternalNodes(1);
    s << "\t" << connectedExternalNodes(2) << "\t" << connectedExternalNodes(3)
      << "\t0.00";
    s << endln;
    s << "PROP_3D\t" << eleTag << "\t";
    s << eleTag << "\t" << 1;
    s << "\t" << -1 << "\tSHELL\t1.0\0.0";
    s << endln;
  }

  if (flag < -1) {

    int counter = (flag + 1) * -1;
    int eleTag  = this->getTag();
    int i, j;
    for (i = 0; i < 4; i++) {
      const Vector &stress = materialPointers[i]->getStressResultant();

      s << "STRESS\t" << eleTag << "\t" << counter << "\t" << i << "\tTOP";
      for (j = 0; j < 6; j++)
        s << "\t" << stress(j);
      s << endln;
    }
  }

  if (flag == OPS_PRINT_CURRENTSTATE) {
    s << endln;
    s << "MITC4 Non-Locking Four Node Shell \n";
    s << "Element Number: " << this->getTag() << endln;
    s << "Node 1 : " << connectedExternalNodes(0) << endln;
    s << "Node 2 : " << connectedExternalNodes(1) << endln;
    s << "Node 3 : " << connectedExternalNodes(2) << endln;
    s << "Node 4 : " << connectedExternalNodes(3) << endln;

    s << "Material Information : \n ";
    materialPointers[0]->Print(s, flag);

    s << endln;
  }

}

Response *ShellMITC4::setResponse(const char **argv, int argc,
                                  OPS_Stream &output)
{
  Response *theResponse = nullptr;

  output.tag("ElementOutput");
  output.attr("eleType", "ShellMITC4");
  output.attr("eleTag", this->getTag());
  int numNodes    = this->getNumExternalNodes();
  const ID &nodes = this->getExternalNodes();
  static char nodeData[32];

  for (int i = 0; i < numNodes; i++) {
    sprintf(nodeData, "node%d", i + 1);
    output.attr(nodeData, nodes(i));
  }

  if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "forces") == 0 ||
      strcmp(argv[0], "globalForce") == 0 ||
      strcmp(argv[0], "globalForces") == 0) {
    const Vector &force = this->getResistingForce();
    int size            = force.Size();
    for (int i = 0; i < size; i++) {
      sprintf(nodeData, "P%d", i + 1);
      output.tag("ResponseType", nodeData);
    }
    theResponse = new ElementResponse(this, 1, this->getResistingForce());
  }

  else if (strcmp(argv[0], "material") == 0 ||
           strcmp(argv[0], "Material") == 0) {
    if (argc < 2) {
      opserr << "ShellMITC4::setResponse() - need to specify more data\n";
      return 0;
    }
    int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= 4) {

      output.tag("GaussPoint");
      output.attr("number", pointNum);
      output.attr("eta",    pts[pointNum - 1][0]);
      output.attr("neta",   pts[pointNum - 1][1]);

      theResponse = materialPointers[pointNum - 1]->setResponse(
          &argv[2], argc - 2, output);

      output.endTag();
    }

  } else if (strcmp(argv[0], "stresses") == 0) {

    for (int i = 0; i < 4; i++) {
      output.tag("GaussPoint");
      output.attr("number", i + 1);
      output.attr("eta", pts[i][0]);
      output.attr("neta", pts[i][1]);

      output.tag("SectionForceDeformation");
      output.attr("classType", materialPointers[i]->getClassTag());
      output.attr("tag", materialPointers[i]->getTag());

      output.tag("ResponseType", "p11");
      output.tag("ResponseType", "p22");
      output.tag("ResponseType", "p1212");
      output.tag("ResponseType", "m11");
      output.tag("ResponseType", "m22");
      output.tag("ResponseType", "m12");
      output.tag("ResponseType", "q1");
      output.tag("ResponseType", "q2");

      output.endTag(); // GaussPoint
      output.endTag(); // NdMaterialOutput
    }

    theResponse = new ElementResponse(this, 2, Vector(32));
  }

  else if (strcmp(argv[0], "strains") == 0) {

    for (int i = 0; i < 4; i++) {
      output.tag("GaussPoint");
      output.attr("number", i + 1);
      output.attr("eta", pts[i][0]);
      output.attr("neta", pts[i][1]);

      output.tag("SectionForceDeformation");
      output.attr("classType", materialPointers[i]->getClassTag());
      output.attr("tag", materialPointers[i]->getTag());

      output.tag("ResponseType", "eps11");
      output.tag("ResponseType", "eps22");
      output.tag("ResponseType", "gamma12");
      output.tag("ResponseType", "theta11");
      output.tag("ResponseType", "theta22");
      output.tag("ResponseType", "theta33");
      output.tag("ResponseType", "gamma13");
      output.tag("ResponseType", "gamma23");

      output.endTag(); // GaussPoint
      output.endTag(); // NdMaterialOutput
    }

    theResponse = new ElementResponse(this, 3, Vector(32));
  }

  output.endTag();
  return theResponse;
}

int ShellMITC4::getResponse(int responseID, Information &eleInfo)
{
  int cnt = 0;
  static Vector stresses(32);
  static Vector strains(32);

  switch (responseID) {
  case 1: // global forces
    return eleInfo.setVector(this->getResistingForce());
    break;

  case 2: // stresses
    for (int i = 0; i < 4; i++) {

      // Get material stress response
      const Vector &sigma = materialPointers[i]->getStressResultant();
      stresses(cnt)       = sigma(0);
      stresses(cnt + 1)   = sigma(1);
      stresses(cnt + 2)   = sigma(2);
      stresses(cnt + 3)   = sigma(3);
      stresses(cnt + 4)   = sigma(4);
      stresses(cnt + 5)   = sigma(5);
      stresses(cnt + 6)   = sigma(6);
      stresses(cnt + 7)   = sigma(7);
      cnt += 8;
    }
    return eleInfo.setVector(stresses);
    break;
  case 3: // strain
    for (int i = 0; i < 4; i++) {

      // Get section deformation
      const Vector &deformation = materialPointers[i]->getSectionDeformation();
      strains(cnt)              = deformation(0);
      strains(cnt + 1)          = deformation(1);
      strains(cnt + 2)          = deformation(2);
      strains(cnt + 3)          = deformation(3);
      strains(cnt + 4)          = deformation(4);
      strains(cnt + 5)          = deformation(5);
      strains(cnt + 6)          = deformation(6);
      strains(cnt + 7)          = deformation(7);
      cnt += 8;
    }
    return eleInfo.setVector(strains);
    break;
  default:
    return -1;
  }
  cnt = 0;
}

const Matrix &
ShellMITC4::getTangentStiff()
{
  int tang_flag = 1; // get the tangent

  // do tangent and residual here
  formResidAndTangent(tang_flag);

  return stiff;
}

// return secant matrix
const Matrix &ShellMITC4::getInitialStiff()
{
  if (Ki != 0)
    return *Ki;

  int i, j, k, p, q;
  int jj, kk;

  double volume = 0.0;

  double xsj; // determinant jacaobian matrix

  double dvol[nip]; // volume element

  double shp[3][NEN]; // shape functions at a gauss point

  //  static double Shape[3][NEN][nip] ; // all the shape functions

  static Matrix stiffJK(ndf, ndf);    // nodeJK stiffness
  static Matrix dd(nstress, nstress); // material tangent

  //---------B-matrices------------------------------------
  static Matrix BJ(nstress, ndf); // B matrix node J
  static Matrix BK(nstress, ndf); // B matrix node k
  static Matrix BJtran(ndf, nstress);
  static Matrix BJtranD(ndf, nstress);

  static Matrix Bbend(3, 3);     // bending B matrix
  static Matrix Bshear(2, 3);    // shear B matrix
  static Matrix Bmembrane(3, 2); // membrane B matrix
  OPS_STATIC double BdrillJ[ndf]; // drill B matrix
  OPS_STATIC double BdrillK[ndf];

  OPS_STATIC double saveB[nstress][ndf][NEN];
  OPS_STATIC MatrixND<nstress, ndf> B[NEN];

  //-------------------------------------------------------

  stiff.Zero();

  double dx34 = xl[0][2] - xl[0][3];
  double dy34 = xl[1][2] - xl[1][3];

  double dx21 = xl[0][1] - xl[0][0];
  double dy21 = xl[1][1] - xl[1][0];

  double dx32 = xl[0][2] - xl[0][1];
  double dy32 = xl[1][2] - xl[1][1];

  double dx41 = xl[0][3] - xl[0][0];
  double dy41 = xl[1][3] - xl[1][0];

  // TODO
  Matrix G(4, 12);
  G.Zero();
  G(0, 0)              = -0.5;
  G(0, 1)              = -dy41 * 0.25;
  G(0, 2)              =  dx41 * 0.25;
  G(0, 9)              = 0.5;
  G(0, 10)             = -dy41 * 0.25;
  G(0, 11)             =  dx41 * 0.25;
  G(1, 0)              = -0.5;
  G(1, 1)              = -dy21 * 0.25;
  G(1, 2)              =  dx21 * 0.25;
  G(1, 3)              = 0.5;
  G(1, 4)              = -dy21 * 0.25;
  G(1, 5)              =  dx21 * 0.25;
  G(2, 3)              = -0.5;
  G(2, 4)              = -dy32 * 0.25;
  G(2, 5)              =  dx32 * 0.25;
  G(2, 6)              = 0.5;
  G(2, 7)              = -dy32 * 0.25;
  G(2, 8)              =  dx32 * 0.25;
  G(3, 6)              = 0.5;
  G(3, 7)              = -dy34 * 0.25;
  G(3, 8)              =  dx34 * 0.25;
  G(3, 9)              = -0.5;
  G(3, 10)             = -dy34 * 0.25;
  G(3, 11)             =  dx34 * 0.25;

  // TODO:
  Matrix Ms(2, 4);
  Ms.Zero();
  Matrix Bs (2, 12);
  Matrix Bsv(2, 12);
  Bsv.Zero();

  double Ax = -xl[0][0] + xl[0][1] + xl[0][2] - xl[0][3];
  double Bx =  xl[0][0] - xl[0][1] + xl[0][2] - xl[0][3];
  double Cx = -xl[0][0] - xl[0][1] + xl[0][2] + xl[0][3];

  double Ay = -xl[1][0] + xl[1][1] + xl[1][2] - xl[1][3];
  double By =  xl[1][0] - xl[1][1] + xl[1][2] - xl[1][3];
  double Cy = -xl[1][0] - xl[1][1] + xl[1][2] + xl[1][3];

  double alph = atan(Ay / Ax);
  double beta = 3.141592653589793 / 2 - atan(Cx / Cy);
  Matrix Rot(2, 2);
  Rot.Zero();
  Rot(0, 0) =  sin(beta);
  Rot(0, 1) = -sin(alph);
  Rot(1, 0) = -cos(beta);
  Rot(1, 1) =  cos(alph);

  double r1 = 0;
  double r2 = 0;
  double r3 = 0;

  // gauss loop
  for (int i = 0; i < nip; i++) {

    r1 = Cx + pts[i][0] * Bx;
    r3 = Cy + pts[i][0] * By;
    r1 = r1 * r1 + r3 * r3;
    r1 = sqrt(r1);
    r2 = Ax + pts[i][1] * Bx;
    r3 = Ay + pts[i][1] * By;
    r2 = r2 * r2 + r3 * r3;
    r2 = sqrt(r2);

    // get shape functions
    shape2d(pts[i][0], pts[i][1], xl, shp, xsj);

    // volume element to also be saved
    dvol[i] = wts[i] * xsj;
    volume += dvol[i];

    Ms(1, 0) = 1 - pts[i][0];
    Ms(0, 1) = 1 - pts[i][1];
    Ms(1, 2) = 1 + pts[i][0];
    Ms(0, 3) = 1 + pts[i][1];
    Bsv      = Ms * G;

    for (int j = 0; j < 12; j++) {
      Bsv(0, j) = Bsv(0, j) * r1 / (8 * xsj);
      Bsv(1, j) = Bsv(1, j) * r2 / (8 * xsj);
    }
    Bs = Rot * Bsv;

    // j-node loop to compute strain
    for (int j = 0; j < NEN; j++) {

      // compute B matrix
      Bmembrane = computeBmembrane(j, shp);
      Bbend     = computeBbend(j, shp);

      for (int p = 0; p < 3; p++) {
        Bshear(0, p) = Bs(0, j * 3 + p);
        Bshear(1, p) = Bs(1, j * 3 + p);
      }

      assembleB(Bmembrane, Bbend, Bshear, B[j]);

      // drilling B matrix
      computeBdrill(j, shp, BdrillJ);
    }

    dd  = materialPointers[i]->getInitialTangent();
    dd *= dvol[i];

    // residual and tangent calculations node loops
    jj = 0;
    for (int j = 0; j < NEN; j++) {

      // extract BJ
      for (int p = 0; p < nstress; p++) {
        for (int q = 0; q < ndf; q++)
          BJ(p, q) = B[j](p,q); // [p][q][j];
      } 

      // multiply bending terms by (-1.0) for correct statement
      // of equilibrium
      for (int p = 3; p < 6; p++) {
        for (int q = 3; q < 6; q++)
          BJ(p, q) *= (-1.0);
      }

      // transpose
      //BJtran = transpose( 8, ndf, BJ ) ;
      for (int p = 0; p < ndf; p++) {
        for (int q = 0; q < nstress; q++)
          BJtran(p, q) = BJ(q, p);
      } 

      // drilling B matrix
      computeBdrill(j, shp, BdrillJ);

      //BJtranD = BJtran * dd ;
      BJtranD.addMatrixProduct(0.0, BJtran, dd, 1.0);

      for (int p = 0; p < ndf; p++)
        BdrillJ[p] *= (Ktt * dvol[i]);

      kk = 0;
      for (int k = 0; k < NEN; k++) {
        // extract BK
        for (int p = 0; p < nstress; p++) {
          for (int q = 0; q < ndf; q++)
            BK(p, q) = saveB[p][q][k];
        } 

        // drilling B matrix
        computeBdrill(k, shp, BdrillK);

        // stiffJK = BJtranD * BK  ;
        // +  transpose( 1,ndf,BdrillJ ) * BdrillK ;
        stiffJK.addMatrixProduct(0.0, BJtranD, BK, 1.0);

        for (int p = 0; p < ndf; p++) {
          for (int q = 0; q < ndf; q++) {
            stiff(jj + p, kk + q) += stiffJK(p, q) + (BdrillJ[p] * BdrillK[q]);
          }
        }

        kk += ndf;
      } // end for k loop

      jj += ndf;
    } // end for j loop

  } // end for i gauss loop

  Ki = new Matrix(stiff);

  return stiff;
}

// return mass matrix
const Matrix &ShellMITC4::getMass()
{
  int tangFlag = 1;

  formInertiaTerms(tangFlag);

  return mass;
}

void ShellMITC4::zeroLoad()
{
  if (load != 0)
    load->Zero();
  applyLoad = 0;

  appliedB[0] = 0.0;
  appliedB[1] = 0.0;
  appliedB[2] = 0.0;

  return;
}

int ShellMITC4::addLoad(ElementalLoad *theLoad, double loadFactor)
{

  int type;
  const Vector &data = theLoad->getData(type, loadFactor);

  if (type == LOAD_TAG_SelfWeight) {
    // added compatability with selfWeight class implemented for all continuum elements, C.McGann, U.W.
    applyLoad = 1;
    appliedB[0] += loadFactor * data(0);
    appliedB[1] += loadFactor * data(1);
    appliedB[2] += loadFactor * data(2);
    return 0;
  } else {
    opserr << "ShellMITC4::addLoad() - ele with tag: " << this->getTag()
           << " does not deal with load type: " << type << "\n";
    return -1;
  }
}

int ShellMITC4::addInertiaLoadToUnbalance(const Vector &accel)
{
  int tangFlag = 1;
  static Vector r(24);

  int allRhoZero = 0;
  for (int i = 0; i < 4; i++) {
    if (materialPointers[i]->getRho() != 0.0)
      allRhoZero = 1;
  }

  if (allRhoZero == 0)
    return 0;

  formInertiaTerms(tangFlag);

  int count = 0;
  for (int i = 0; i < 4; i++) {
    const Vector &Raccel = theNodes[i]->getRV(accel);
    for (int j = 0; j < 6; j++)
      r(count++) = Raccel(j);
  }

  if (load == 0)
    load = new Vector(24);

  load->addMatrixVector(1.0, mass, r, -1.0);

  return 0;
}

// get residual
const Vector &ShellMITC4::getResistingForce()
{
  int tang_flag = 0; // don't get the tangent

  formResidAndTangent(tang_flag);

  // subtract external loads
  if (load != 0)
    resid -= *load;

  return resid;
}

// get residual with inertia terms
const Vector &
ShellMITC4::getResistingForceIncInertia()
{
  static Vector res(24);
  int tang_flag = 0; // don't get the tangent

  // do tangent and residual here
  formResidAndTangent(tang_flag);

  formInertiaTerms(tang_flag);

  res = resid;
  // add the damping forces if rayleigh damping
  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
    res += this->getRayleighDampingForces();

  // subtract external loads
  if (load != 0)
    res -= *load;

  return res;
}

// form inertia terms
void ShellMITC4::formInertiaTerms(int tangFlag)
{
  // translational mass only
  // rotational inertia terms are neglected

  static constexpr int ndf = 6;
  static constexpr int numberNodes = 4;
  static constexpr int nShape = 3;
  static constexpr int massIndex = nShape - 1;

  double xsj; // determinant jacaobian matrix
  double dvol; // volume element
  OPS_STATIC double shp[nShape][numberNodes]; // shape functions at a gauss point

  static Vector momentum(ndf);

  double temp, rhoH, massJK;

  // zero mass
  mass.Zero();

  // Gauss loop
  for (int i = 0; i < nip; i++) {

    // get shape functions
    shape2d(pts[i][0], pts[i][1], xl, shp, xsj);

    // volume element to also be saved
    dvol = wts[i] * xsj;

    // node loop to compute accelerations
    momentum.Zero();
    for (int j = 0; j < numberNodes; j++)
      // momentum += ( shp[massIndex][j] * theNodes[j]->getTrialAccel() ) ;
      momentum.addVector(1.0, theNodes[j]->getTrialAccel(),
                         shp[massIndex][j]);

    // density
    rhoH = materialPointers[i]->getRho();

    // multiply acceleration by density to form momentum
    momentum *= rhoH;

    // residual and tangent calculations node loops
    for (int j = 0, jj = 0; j < numberNodes; j++, jj += ndf) {

      temp = shp[massIndex][j] * dvol;

      for (int p = 0; p < 3; p++)
        resid(jj + p) += (temp * momentum(p));

      if (tangFlag == 1 && rhoH != 0.0) {

        // multiply by density
        temp *= rhoH;

        // node-node translational mass
        for (int k = 0, kk = 0; k < numberNodes; k++, kk += ndf) {

          massJK = temp * shp[massIndex][k];

          for (int p = 0; p < 3; p++)
            mass(jj + p, kk + p) += massJK;

        } // end for k loop

      } // end if tang_flag

    } // end for j loop

  } // end for i gauss loop
}

//*********************************************************************

// form residual and tangent
//
//  six(6) nodal dof's ordered :
//
//    -        -
//   |    u1    |   <---plate membrane
//   |    u2    |
//   | -------- |
//   |  w = u3  |   <---plate bending
//   |  theta1  |
//   |  theta2  |
//   | -------- |
//   |  theta3  |   <---drill
//    -        -
//
// membrane strains ordered :
//
//            strain(0) =   eps00     i.e.   (11)-strain
//            strain(1) =   eps11     i.e.   (22)-strain
//            strain(2) =   gamma01   i.e.   (12)-shear
//
// curvatures and shear strains ordered  :
//
//            strain(3) =     kappa00  i.e.   (11)-curvature
//            strain(4) =     kappa11  i.e.   (22)-curvature
//            strain(5) =   2*kappa01  i.e. 2*(12)-curvature
//
//            strain(6) =     gamma02  i.e.   (13)-shear
//            strain(7) =     gamma12  i.e.   (23)-shear
//
//  same ordering for moments/shears but no 2
//
//  Then,
//              epsilon00 = -z * kappa00      +    eps00_membrane
//              epsilon11 = -z * kappa11      +    eps11_membrane
//  gamma01 = 2*epsilon01 = -z * (2*kappa01)  +  gamma01_membrane
//
//  Shear strains gamma02, gamma12 constant through cross section
//
void
ShellMITC4::formResidAndTangent(int tang_flag)
{

  int success;
  double volume = 0.0;
  double xsj;                   // determinant jacaobian matrix

  OPS_STATIC double shp[3][NEN];      // shape functions at a gauss point

  //  static double Shape[3][NEN][nip] ; // all the shape functions
  static Vector stress(nstress);      // stress resultants
  static Vector strain(nstress);      // strain
                                      //
  VectorND<ndf> residJ;
  MatrixND<nstress,nstress> dd; // material tangent


  double epsDrill = 0.0; // drilling "strain"
  double tauDrill = 0.0; // drilling "stress"

  //---------B-matrices------------------------------------
  OPS_STATIC MatrixND<nstress, ndf> B[NEN];
  static Matrix BJtranD(ndf, nstress);
  static Matrix Bbend(3, 3);           // bending B matrix
  static Matrix Bshear(2, 3);          // shear B matrix

  static Matrix Bmembrane(3, 2);       // membrane B matrix
  OPS_STATIC double BdrillJ[ndf];      // drill B matrix
  OPS_STATIC double BdrillK[ndf];
  //-------------------------------------------------------

  // zero stiffness and residual
  stiff.Zero();
  resid.Zero();

  // start Yuli Huang (yulihuang@gmail.com) & Xinzheng Lu (luxz@tsinghua.edu.cn)
  if (doUpdateBasis == true)
    updateBasis();
  // end Yuli Huang (yulihuang@gmail.com) & Xinzheng Lu (luxz@tsinghua.edu.cn)

  double dx34 = xl[0][2] - xl[0][3];
  double dy34 = xl[1][2] - xl[1][3];

  double dx21 = xl[0][1] - xl[0][0];
  double dy21 = xl[1][1] - xl[1][0];

  double dx32 = xl[0][2] - xl[0][1];
  double dy32 = xl[1][2] - xl[1][1];

  double dx41 = xl[0][3] - xl[0][0];
  double dy41 = xl[1][3] - xl[1][0];


  double Ax = -xl[0][0] + xl[0][1] + xl[0][2] - xl[0][3];
  double Bx =  xl[0][0] - xl[0][1] + xl[0][2] - xl[0][3];
  double Cx = -xl[0][0] - xl[0][1] + xl[0][2] + xl[0][3];

  double Ay = -xl[1][0] + xl[1][1] + xl[1][2] - xl[1][3];
  double By =  xl[1][0] - xl[1][1] + xl[1][2] - xl[1][3];
  double Cy = -xl[1][0] - xl[1][1] + xl[1][2] + xl[1][3];

  double alph = atan(Ay / Ax);
  double beta = 3.141592653589793 / 2 - atan(Cx / Cy);

  MatrixND<2, 2> Rot;
  Rot.zero();
  Rot(0, 0) =  sin(beta);
  Rot(0, 1) = -sin(alph);
  Rot(1, 0) = -cos(beta);
  Rot(1, 1) =  cos(alph);

  const MatrixND<4, 12> G = {{ // NOTE: initialization is transposed
   {     -0.50,        -0.50,          0.00,          0.00},
   {-dy41*0.25, -dy21 * 0.25,          0.00,          0.00},
   { dx41*0.25,  dx21 * 0.25,          0.00,          0.00},
   {      0.00,         0.50,         -0.5 ,          0.00},
   {      0.00, -dy21 * 0.25,  -dy32 * 0.25,          0.00},
   {      0.00,  dx21 * 0.25,   dx32 * 0.25,          0.00},
   {      0.00,         0.00,          0.5 ,          0.5 },
   {      0.00,         0.00,  -dy32 * 0.25,  -dy34 * 0.25},
   {      0.00,         0.00,   dx32 * 0.25,   dx34 * 0.25},
   {      0.50,         0.00,          0.00,         -0.5 },
   {-dy41*0.25,         0.00,          0.00,  -dy34 * 0.25},
   { dx41*0.25,         0.00,          0.00,   dx34 * 0.25}}};

  MatrixND<2, 4> Ms;
  Ms.zero();

  double r1 = 0;
  double r2 = 0;
  double r3 = 0;

  double dvol[nip];             // volume element
  // Gauss loop
  for (int i = 0; i < nip; i++) {

    r1 = Cx + pts[i][0] * Bx;
    r3 = Cy + pts[i][0] * By;
    r1 = r1 * r1 + r3 * r3;
    r1 = sqrt(r1);
    r2 = Ax + pts[i][1] * Bx;
    r3 = Ay + pts[i][1] * By;
    r2 = r2 * r2 + r3 * r3;
    r2 = sqrt(r2);

    // get shape functions
    shape2d(pts[i][0], pts[i][1], xl, shp, xsj);
    // volume element to also be saved
    dvol[i] = wts[i] * xsj;
    volume += dvol[i];

    Ms(1, 0) = 1 - pts[i][0];
    Ms(0, 1) = 1 - pts[i][1];
    Ms(1, 2) = 1 + pts[i][0];
    Ms(0, 3) = 1 + pts[i][1];
    auto Bsv      = Ms * G;

    for (int j = 0; j < 12; j++) {
      Bsv(0, j) = Bsv(0, j) * r1 / (8 * xsj);
      Bsv(1, j) = Bsv(1, j) * r2 / (8 * xsj);
    }
    auto Bs = Rot * Bsv;

    // zero the strains
    strain.Zero();
    epsDrill = 0.0;

    // j-node loop to compute strain
    for (int j = 0; j < NEN; j++) {

      // compute B matrix
      Bmembrane = computeBmembrane(j, shp);

      Bbend = computeBbend(j, shp);

      for (int p = 0; p < 3; p++) {
        Bshear(0, p) = Bs(0, j * 3 + p);
        Bshear(1, p) = Bs(1, j * 3 + p);
      }

      assembleB(Bmembrane, Bbend, Bshear, B[j]);

      // nodal "displacements"
      const Vector &ul_tmp = theNodes[j]->getTrialDisp();

      OPS_STATIC VectorND<6> ul;

      ul[0] = ul_tmp(0) - init_disp[j][0];
      ul[1] = ul_tmp(1) - init_disp[j][1];
      ul[2] = ul_tmp(2) - init_disp[j][2];
      ul[3] = ul_tmp(3) - init_disp[j][3];
      ul[4] = ul_tmp(4) - init_disp[j][4];
      ul[5] = ul_tmp(5) - init_disp[j][5];

      // compute the strain
      // strain += (BJ*ul) ;
      strain.addMatrixVector(1.0, B[j], ul, 1.0);

      // drilling B matrix
      computeBdrill(j, shp, BdrillJ);

      // drilling "strain"
      for (int p = 0; p < ndf; p++)
        epsDrill += BdrillJ[p] * ul[p];
    }

    // send the strain to the material
    success = materialPointers[i]->setTrialSectionDeformation(strain);

    // compute the stress
    /* VectorND<nstress> */ stress = materialPointers[i]->getStressResultant();

    // drilling "stress"
    tauDrill = Ktt * epsDrill;

    // multiply by volume element
    stress *= dvol[i];
    tauDrill *= dvol[i];

    if (tang_flag == 1) {
      dd = materialPointers[i]->getSectionTangent();
      dd *= dvol[i];
    }

    // residual and tangent calculations node loops

    int jj = 0;
    for (int j = 0; j < NEN; j++) {
      MatrixND<ndf, nstress> BJtran = B[j].transpose();

      // multiply bending terms by (-1.0) for correct statement
      // of equilibrium
      for (int p = 3; p < 6; p++) {
        for (int q = 3; q < 6; q++)
          BJtran(p, q) *= -1.0;
      }

//    residJ.addMatrixTransposeVector(0.0, B[j], stress, 1.0);
      residJ = BJtran * stress;

      // drilling B matrix
      computeBdrill(j, shp, BdrillJ);

      // residual including drill
      for (int p = 0; p < ndf; p++)
        resid(jj + p) += (residJ[p] + BdrillJ[p] * tauDrill);

      if (tang_flag == 1) {

//      BJtranD.addMatrixTransposeProduct(0.0, B[j], dd, 1.0);
        MatrixND<ndf, nstress> BJtranD = BJtran*dd;

        for (int p = 0; p < ndf; p++)
          BdrillJ[p] *= (Ktt * dvol[i]);

        int kk = 0;
        for (int k = 0; k < NEN; k++) {

          // drilling B matrix
          computeBdrill(k, shp, BdrillK);

          // stiffJK = BJtranD * BK  ;
          // +  transpose( 1,ndf,BdrillJ ) * BdrillK ;
          MatrixND<ndf,ndf> stiffJK = BJtranD * B[k]; // .addMatrixProduct(0.0, BJtranD, BK, 1.0);

          for (int p = 0; p < ndf; p++) {
            for (int q = 0; q < ndf; q++) {
              stiff(jj + p, kk + q) +=
                  stiffJK(p, q) + (BdrillJ[p] * BdrillK[q]);
            }
          }

          kk += ndf;
        } // end for k loop

      } // end if tang_flag

      jj += ndf;
    } // end for j loop

  } // end gauss loop

  if (applyLoad == 1) {
    const int nShape      = 3;
    const int numberNodes = 4;
    const int massIndex   = nShape - 1;
    double temp, rhoH;
    // If defined, apply self-weight
    static Vector momentum(ndf);
    double ddvol = 0;

    for (int i = 0; i < nip; i++) {

      // get shape functions
      shape2d(pts[i][0], pts[i][1], xl, shp, xsj);

      // volume element to also be saved
      ddvol = wts[i] * xsj;

      // node loop to compute accelerations
      momentum.Zero();
      momentum(0) = appliedB[0];
      momentum(1) = appliedB[1];
      momentum(2) = appliedB[2];

      // density
      rhoH = materialPointers[i]->getRho();

      // multiply acceleration by density to form momentum
      momentum *= rhoH;

      // residual and tangent calculations node loops
      for (int j = 0, jj = 0; j < numberNodes; j++, jj += ndf) {

        temp = shp[massIndex][j] * ddvol;

        for (int p = 0; p < 3; p++)
          resid(jj + p) += (temp * momentum(p));
      }
    }
  }
  return;
}

// ************************************************************************
// compute local coordinates and basis

void ShellMITC4::computeBasis()
{
  // could compute derivatives \frac{ \partial {\bf x} }{ \partial L_1 }
  //                     and  \frac{ \partial {\bf x} }{ \partial L_2 }
  // and use those as basis vectors but this is easier
  // and the shell is flat anyway.

  Vector3D temp;
  Vector3D v1;
  Vector3D v2;
  Vector3D v3;

  // get two vectors (v1, v2) in plane of shell by
  // nodal coordinate differences

  const Vector& coor0 = theNodes[0]->getCrds();
  const Vector& coor1 = theNodes[1]->getCrds();
  const Vector& coor2 = theNodes[2]->getCrds();
  const Vector& coor3 = theNodes[3]->getCrds();

  v1.zero();
  // v1 = 0.5 * ( coor2 + coor1 - coor3 - coor0 ) ;
  v1  = coor2;
  v1 += coor1;
  v1 -= coor3;
  v1 -= coor0;
  v1 *= 0.50;

  v2.zero();
  // v2 = 0.5 * ( coor3 + coor2 - coor1 - coor0 ) ;
  v2  = coor3;
  v2 += coor2;
  v2 -= coor1;
  v2 -= coor0;
  v2 *= 0.50;

  // normalize v1
  double length = v1.norm();
  v1 /= length;

  //Gram-Schmidt process for v2

  double alpha = v2.dot(v1);

  // v2 -= alpha*v1 ;
  temp = v1;
  temp *= alpha;
  v2 -= temp;

  // normalize v2
  length = v2.norm();
  v2 /= length;

  // cross product for v3
  v3 = v1.cross(v2);

  // local nodal coordinates in plane of shell
  for (int i = 0; i < 4; i++) {

    const Vector &coorI = theNodes[i]->getCrds();
    xl[0][i]            = v1.dot(coorI);
    xl[1][i]            = v2.dot(coorI);

  } // end for i

  // basis vectors stored as array of doubles
  for (int i = 0; i < 3; i++) {
    g1[i] = v1(i);
    g2[i] = v2(i);
    g3[i] = v3(i);
  }
}

// start Yuli Huang (yulihuang@gmail.com) & Xinzheng Lu (luxz@tsinghua.edu.cn)
//************************************************************************
// compute local coordinates and basis

void ShellMITC4::updateBasis()
{
  // could compute derivatives \frac{ \partial {\bf x} }{ \partial L_1 }
  //                     and  \frac{ \partial {\bf x} }{ \partial L_2 }
  // and use those as basis vectors but this is easier
  // and the shell is flat anyway.

  
  // get two vectors (v1, v2) in plane of shell by
  // nodal coordinate differences
  Vector3D coor[4];
  for (int i=0; i<4; i++) {
    coor[i] = theNodes[i]->getCrds();
    const Vector& displ = theNodes[i]->getTrialDisp();
    for (int j=0; j<3; j++)
      coor[i][j] += displ[j] - init_disp[i][j];
  }

  Vector3D v1;
  Vector3D v2;
  Vector3D temp;

  // v1 = 0.5 * ( coor2 + coor1 - coor3 - coor0 ) ;
  v1  = coor[2];
  v1 += coor[1];
  v1 -= coor[3];
  v1 -= coor[0];
  v1 *= 0.50;

  // v2 = 0.5 * ( coor3 + coor2 - coor1 - coor0 ) ;
  v2  = coor[3];
  v2 += coor[2];
  v2 -= coor[1];
  v2 -= coor[0];
  v2 *= 0.50;

  // normalize v1
  double length = v1.norm();
  v1 /= length;

  // Gram-Schmidt process for v2
  double alpha = v2.dot(v1);

  // v2 -= alpha*v1 ;
  temp = v1;
  temp *= alpha;
  v2 -= temp;

  // normalize v2
  length = v2.norm();
  v2 /= length;

  // cross product for v3
  Vector3D v3 = v1.cross(v2);

  // local nodal coordinates in plane of shell
  for (int i = 0; i < 4; i++) {

    const Vector &coorI = theNodes[i]->getCrds();
    xl[0][i]            = v1.dot(coorI);
    xl[1][i]            = v2.dot(coorI);

  }

  // basis vectors stored as array of doubles
  for (int i = 0; i < 3; i++) {
    g1[i] = v1(i);
    g2[i] = v2(i);
    g3[i] = v3(i);
  }
}

//*************************************************************************
// compute Bdrill
inline
double *ShellMITC4::computeBdrill(int node, const double shp[3][4], double Bdrill[6])
{

  //---Bdrill Matrix in standard {1,2,3} mechanics notation---------
  //
  //             -                                       -
  //   Bdrill = | -0.5*N,2   +0.5*N,1    0    0    0   -N |   (1x6)
  //             -                                       -
  //
  //----------------------------------------------------------------

  double B1 = -0.5 * shp[1][node];
  double B2 = +0.5 * shp[0][node];
  double B6 =       -shp[2][node];

  Bdrill[0] = B1 * g1[0] + B2 * g2[0];
  Bdrill[1] = B1 * g1[1] + B2 * g2[1];
  Bdrill[2] = B1 * g1[2] + B2 * g2[2];

  Bdrill[3] = B6 * g3[0];
  Bdrill[4] = B6 * g3[1];
  Bdrill[5] = B6 * g3[2];

  return Bdrill;
}


// assemble a B matrix
void
ShellMITC4::assembleB(const Matrix &Bmembrane,
                      const Matrix &Bbend, const Matrix &Bshear,
                      MatrixND<nstress, ndf> &B)
{
  //
  //---B Matrices in standard {1,2,3} mechanics notation---------
  //
  //            -                     _
  //           | Bmembrane  |     0    |
  //           | --------------------- |
  //    B =    |     0      |  Bbend   |   (8x6)
  //           | --------------------- |
  //           |         Bshear        |
  //            -           -         -
  //
  //-------------------------------------------------------------
  //
  static Matrix BmembraneShell(3, 3);
  static Matrix BbendShell(3, 3);
  static Matrix BshearShell(2, 6);
  static Matrix Gmem(2, 3);
  static Matrix Gshear(3, 6);

  // shell modified membrane terms

  Gmem(0, 0) = g1[0];
  Gmem(0, 1) = g1[1];
  Gmem(0, 2) = g1[2];

  Gmem(1, 0) = g2[0];
  Gmem(1, 1) = g2[1];
  Gmem(1, 2) = g2[2];

  //BmembraneShell = Bmembrane * Gmem ;
  BmembraneShell.addMatrixProduct(0.0, Bmembrane, Gmem, 1.0);

  // shell modified bending terms
  Matrix &Gbend = Gmem;

  //BbendShell = Bbend * Gbend ;
  BbendShell.addMatrixProduct(0.0, Bbend, Gbend, 1.0);

  // shell modified shear terms
  Gshear.Zero();

  Gshear(0, 0) = g3[0];
  Gshear(0, 1) = g3[1];
  Gshear(0, 2) = g3[2];

  Gshear(1, 3) = g1[0];
  Gshear(1, 4) = g1[1];
  Gshear(1, 5) = g1[2];

  Gshear(2, 3) = g2[0];
  Gshear(2, 4) = g2[1];
  Gshear(2, 5) = g2[2];

  // BshearShell = Bshear * Gshear ;
  BshearShell.addMatrixProduct(0.0, Bshear, Gshear, 1.0);

  B.zero();

  // assemble B from sub-matrices

  // membrane terms
  for (int p = 0; p < 3; p++) {
    for (int q = 0; q < 3; q++)
      B(p, q) = BmembraneShell(p, q);
  }

  // bending terms
  for (int p = 3; p < 6; p++) {
    int pp = p - 3;
    for (int q = 3; q < 6; q++)
      B(p, q) = BbendShell(pp, q - 3);
  }

  // shear terms
  for (int p = 0; p < 2; p++) {
    int pp = p + 6;

    for (int q = 0; q < 6; q++)
      B(pp, q) = BshearShell(p, q);

  }
}

//
// compute Bmembrane matrix
//
const Matrix &
ShellMITC4::computeBmembrane(int node, const double shp[3][4])
{
  //---Bmembrane Matrix in standard {1,2,3} mechanics notation---------
  //
  //                -             -
  //               | +N,1      0   |
  // Bmembrane =   |   0     +N,2  |    (3x2)
  //               | +N,2    +N,1  |
  //                -             -
  //
  //  three(3) strains and two(2) displacements (for plate)
  //-------------------------------------------------------------------

  static Matrix Bmembrane(3, 2);

  Bmembrane.Zero();

  Bmembrane(0, 0) = shp[0][node];
  Bmembrane(1, 1) = shp[1][node];
  Bmembrane(2, 0) = shp[1][node];
  Bmembrane(2, 1) = shp[0][node];

  return Bmembrane;
}

//
// compute Bbend matrix
//
const Matrix &
ShellMITC4::computeBbend(int node, const double shp[3][4])
{
//
//---Bbend Matrix in standard {1,2,3} mechanics notation---------
//
//            -             -
//   Bbend = |    0    -N,1  |
//           |  +N,2     0   |    (3x2)
//           |  +N,1   -N,2  |
//            -             -
//
//  three(3) curvatures and two(2) rotations (for plate)
//----------------------------------------------------------------
  static Matrix Bbend(3, 2);
  Bbend.Zero();

  Bbend(0, 1) = -shp[0][node];
  Bbend(1, 0) =  shp[1][node];
  Bbend(2, 0) =  shp[0][node];
  Bbend(2, 1) = -shp[1][node];

  return Bbend;
}

//************************************************************************
// shape function routine for four node quads
void ShellMITC4::shape2d(double ss, double tt, const double x[2][4],
                         double shp[3][4], double &xsj)
{

//int i, j, k;
  constexpr static const double s[] = {-0.5, 0.5, 0.5, -0.5};
  constexpr static const double t[] = {-0.5, -0.5, 0.5, 0.5};

  double xs[2][2];
  double sx[2][2];

  for (int i = 0; i < 4; i++) {
    shp[2][i] = (0.5 + s[i] * ss) * (0.5 + t[i] * tt);
    shp[0][i] = s[i] * (0.5 + t[i] * tt);
    shp[1][i] = t[i] * (0.5 + s[i] * ss);
  }

  // Construct jacobian and its inverse

  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {

      xs[i][j] = 0.0;

      for (int k = 0; k < 4; k++)
        xs[i][j] += x[i][k] * shp[j][k];

    }
  }

  xsj = xs[0][0] * xs[1][1] - xs[0][1] * xs[1][0];

  // inverse jacobian
  double jinv = 1.0 / xsj;
  sx[0][0]    =  xs[1][1] * jinv;
  sx[1][1]    =  xs[0][0] * jinv;
  sx[0][1]    = -xs[0][1] * jinv;
  sx[1][0]    = -xs[1][0] * jinv;

  // form global derivatives

  for (int i = 0; i < 4; i++) {
    double temp = shp[0][i] * sx[0][0] + shp[1][i] * sx[1][0];
    shp[1][i]   = shp[0][i] * sx[0][1] + shp[1][i] * sx[1][1];
    shp[0][i]   = temp;
  }

  return;
}

int
ShellMITC4::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;

  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();

  // Now send the ids of materials
  int matDbTag;

  static ID idData(14);

  for (int i = 0; i < 4; i++) {
    idData(i) = materialPointers[i]->getClassTag();
    matDbTag  = materialPointers[i]->getDbTag();
    // NOTE: we do have to ensure that the material has a database
    // tag if we are sending to a database channel.
    if (matDbTag == 0) {
      matDbTag = theChannel.getDbTag();
      if (matDbTag != 0)
        materialPointers[i]->setDbTag(matDbTag);
    }
    idData(i + 4) = matDbTag;
  }

  idData(8)  = this->getTag();
  idData(9)  = connectedExternalNodes(0);
  idData(10) = connectedExternalNodes(1);
  idData(11) = connectedExternalNodes(2);
  idData(12) = connectedExternalNodes(3);
  if (doUpdateBasis == true)
    idData(13) = 0;
  else
    idData(13) = 1;

  res += theChannel.sendID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING ShellMITC4::sendSelf() - " << this->getTag()
           << " failed to send ID\n";
    return res;
  }

  static Vector vectData(5 + 6 * 4);
  vectData(0) = Ktt;
  vectData(1) = alphaM;
  vectData(2) = betaK;
  vectData(3) = betaK0;
  vectData(4) = betaKc;

  int pos = 0;
  for (int node = 0; node < 4; ++node) {
    for (int dof = 0; dof < 6; ++dof) {
      vectData(5 + pos) = init_disp[node][dof];
      pos++;
    }
  }

  res += theChannel.sendVector(dataTag, commitTag, vectData);
  if (res < 0) {
    opserr << "WARNING ShellMITC4::sendSelf() - " << this->getTag()
           << " failed to send ID\n";
    return res;
  }

  // Finally, ask material objects to send themselves
  for (int i = 0; i < 4; i++) {
    res += materialPointers[i]->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING ShellMITC4::sendSelf() - " << this->getTag()
             << " failed to send its Material\n";
      return res;
    }
  }

  return res;
}

int ShellMITC4::recvSelf(int commitTag, Channel &theChannel,
                         FEM_ObjectBroker &theBroker)
{
  int res = 0;

  int dataTag = this->getDbTag();

  static ID idData(14);
  // Quad now receives the tags of its four external nodes
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING ShellMITC4::recvSelf() - " << this->getTag()
           << " failed to receive ID\n";
    return res;
  }

  this->setTag(idData(8));
  connectedExternalNodes(0) = idData(9);
  connectedExternalNodes(1) = idData(10);
  connectedExternalNodes(2) = idData(11);
  connectedExternalNodes(3) = idData(12);
  if (idData(13) == 0)
    doUpdateBasis = true;
  else
    doUpdateBasis = false;

  static Vector vectData(5 + 6 * 4);
  res += theChannel.recvVector(dataTag, commitTag, vectData);
  if (res < 0) {
    opserr << "WARNING ShellMITC4::sendSelf() - " << this->getTag()
           << " failed to send ID\n";
    return res;
  }

  Ktt    = vectData(0);
  alphaM = vectData(1);
  betaK  = vectData(2);
  betaK0 = vectData(3);
  betaKc = vectData(4);

  int pos = 0;
  for (int node = 0; node < 4; ++node) {
    for (int dof = 0; dof < 6; ++dof) {
      init_disp[node][dof] = vectData(5 + pos);
      pos++;
    }
  }

  int i;

  if (materialPointers[0] == 0) {
    for (i = 0; i < 4; i++) {
      int matClassTag = idData(i);
      int matDbTag    = idData(i + 4);
      // Allocate new material with the sent class tag
      materialPointers[i] = theBroker.getNewSection(matClassTag);
      if (materialPointers[i] == 0) {
        opserr << "ShellMITC4::recvSelf() - Broker could not create NDMaterial "
                  "of class type"
               << matClassTag << endln;
        ;
        return -1;
      }
      // Now receive materials into the newly allocated space
      materialPointers[i]->setDbTag(matDbTag);
      res += materialPointers[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
        opserr << "ShellMITC4::recvSelf() - material " << i
               << "failed to recv itself\n";
        return res;
      }
    }
  }
  // Number of materials is the same, receive materials into current space
  else {
    for (i = 0; i < 4; i++) {
      int matClassTag = idData(i);
      int matDbTag    = idData(i + 4);
      // Check that material is of the right type; if not,
      // delete it and create a new one of the right type
      if (materialPointers[i]->getClassTag() != matClassTag) {
        delete materialPointers[i];
        materialPointers[i] = theBroker.getNewSection(matClassTag);
        if (materialPointers[i] == 0) {
          opserr << "ShellMITC4::recvSelf() - Broker could not create "
                    "NDMaterial of class type"
                 << matClassTag << endln;
          exit(-1);
        }
      }
      // Receive the material
      materialPointers[i]->setDbTag(matDbTag);
      res += materialPointers[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
        opserr << "ShellMITC4::recvSelf() - material " << i
               << "failed to recv itself\n";
        return res;
      }
    }
  }

  return res;
}

int
ShellMITC4::setParameter(const char **argv, int argc, Parameter &param)
{
  int res = -1;
  // Send to all sections
  for (int i = 0; i < nip; i++) {
    int secRes = materialPointers[i]->setParameter(argv, argc, param);
    if (secRes != -1) {
      res = secRes;
    }
  }
  return res;
}
