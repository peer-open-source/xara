//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2026, OpenSees/Xara Developers
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
//
#include <ExactTruss.h>
#include <Information.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <FrameSection.h>
#include <VectorND.h>
#include <Parameter.h>
using namespace OpenSees;
#include <cmath>
#include <stdlib.h>
#include <string.h>

#include <ElementResponse.h>

Matrix ExactTruss::M2(2, 2);
Matrix ExactTruss::M4(4, 4);
Matrix ExactTruss::M6(6, 6);
Matrix ExactTruss::M12(12, 12);

Vector ExactTruss::V2(2);
Vector ExactTruss::V4(4);
Vector ExactTruss::V6(6);
Vector ExactTruss::V12(12);


ExactTruss::ExactTruss(int tag, int dim, 
                        int Nd1, int Nd2,
                        FrameSection& theSec, 
                        double r, 
                        int damp, 
                        int cm,
                        int strain)
 : Element(tag, -1),
   theSection(0),
   connectedExternalNodes(2),
   numDOF(0),
   numDIM(dim),
   Lo(0.0),
   rho(r),
   doRayleighDamping(damp),
   cMass(cm),
   R{},
   parameterID(0),
   theLoad(0),
   theMatrix(0),
   theVector(0)
{
  theSection = theSec.getFrameCopy(section_layout);

  connectedExternalNodes(0) = Nd1;
  connectedExternalNodes(1) = Nd2;

  switch (strain) {
    case 0:
      strain_type = StrainType::Linear;
      break;
    case 1:
      strain_type = StrainType::Green;
      break;
    case 2:
      strain_type = StrainType::Corotational;
      break;
    case 3:
      strain_type = StrainType::Logarithmic;
      break;
    default:
      // assert(false && "ExactTruss::ExactTruss() - invalid strain type");
      strain_type = StrainType::Corotational;
  }

  for (int i = 0; i < 2; i++)
    theNodes[i] = nullptr;
}


ExactTruss::~ExactTruss()
{
  if (theSection != nullptr)
    delete theSection;
}

int
ExactTruss::getNumExternalNodes() const
{
  return 2;
}

const ID&
ExactTruss::getExternalNodes()
{
  return connectedExternalNodes;
}


Node**
ExactTruss::getNodePtrs()
{
  return theNodes;
}

int
ExactTruss::getNumDOF()
{
  return numDOF;
}


void
ExactTruss::setDomain(Domain* theDomain)
{
  if (theDomain == nullptr) {
    theNodes[0] = 0;
    theNodes[1] = 0;
    Lo          = 0.0;
    return;
  }

  // Set node pointers
  theNodes[0] = theDomain->getNode(connectedExternalNodes(0));
  theNodes[1] = theDomain->getNode(connectedExternalNodes(1));

  if ((theNodes[0] == nullptr) || (theNodes[1] == nullptr)) {
    opserr << "ExactTruss::setDomain() - ExactTruss " << this->getTag()
           << " node does not exist in the model\n";
    numDOF = 6;
    return;
  }

  // now determine the number of dof and the dimension
  int dofNd1 = theNodes[0]->getNumberDOF();

  // if differing dof at the ends print a warning message
  if (dofNd1 != theNodes[1]->getNumberDOF()) {
    opserr << "WARNING ExactTruss::setDomain(): nodes have differing dof at ends for "
              "ExactTruss"
           << this->getTag() << "\n";
    numDOF = 6;
    return;
  }

  if (numDIM == 1 && dofNd1 == 1) {
    numDOF    = 2;
    theMatrix = &M2;
    theVector = &V2;
  } else if (numDIM == 2 && dofNd1 == 2) {
    numDOF    = 4;
    theMatrix = &M4;
    theVector = &V4;
  } else if (numDIM == 2 && dofNd1 == 3) {
    numDOF    = 6;
    theMatrix = &M6;
    theVector = &V6;
  } else if (numDIM == 3 && dofNd1 == 3) {
    numDOF    = 6;
    theMatrix = &M6;
    theVector = &V6;
  } else if (numDIM == 3 && dofNd1 == 6) {
    numDOF    = 12;
    theMatrix = &M12;
    theVector = &V12;
  } else {
    opserr << "nodal DOF not compatible with element "
           << this->getTag() << "\n";

    numDOF = 6;
    return;
  }

  // create the load vector
  if (theLoad == 0)
    theLoad = new Vector(numDOF);
  else if (theLoad->Size() != numDOF) {
    delete theLoad;
    theLoad = new Vector(numDOF);
  }

  // call the base class method
  if (theDomain != nullptr)
    this->Element::link(*theDomain);

  Vector3D dX, dx, du;
  this->locate(dX, dx, du);
  Lo = dX.norm();

  R.addDiagonal(1.0);
}

int
ExactTruss::commitState()
{
  int retVal = 0;
  if ((retVal = this->Element::commitState()) != 0) {
    opserr << "ExactTruss::commitState () - failed in base class";
  }
  retVal = theSection->commitState();
  return retVal;
}

int
ExactTruss::revertToLastCommit()
{
  return theSection->revertToLastCommit();
}

int
ExactTruss::revertToStart()
{
  return theSection->revertToStart();
}

int
ExactTruss::update()
{
  if (Lo == 0.0) { 
    // problem in setDomain(); no further warnings
    return -1;
  }

  // determine the current strain given trial displacements at nodes
  double strain = this->computeCurrentStrain();

  VectorND<1> e{strain};
  return theSection->setTrialState<1,section_layout>(e);
}


const Matrix&
ExactTruss::getTangentStiff()
{
  // Material stiffness
  //
  // Get material tangent
  double EA = theSection->getTangent<1,section_layout>(State::Pres)(0,0);
  double q  = theSection->getResultant<1,section_layout>()[0];

  Vector3D b;
  double Ln;
  double J = this->interpolate(b, Ln);

  // EA /= (Ln * Ln * Lo);
  EA *= Lo;

  Matrix3D kl{};
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      kl(i, j) = EA * b[i] * b[j];

  // Geometric stiffness
  //
  // Get material stress
  double SA = q / (Ln * Ln * Ln);
  double SL = q / Ln;

  for (int i = 0; i < 3; i++) {
    kl(i, i) += SL;
    for (int j = 0; j < 3; j++)
      kl(i, j) -= SA * b[i] * b[j];
  }

  // Compute R'*kl*R
  Matrix3D kg{};
  kg.addMatrixTripleProduct(0.0, R, kl, 1.0);

  Matrix& K = *theMatrix;
  K.Zero();

  // Copy stiffness into appropriate blocks in element stiffness
  int numDOF2 = numDOF / 2;
  for (int i = 0; i < numDIM; i++) {
    for (int j = 0; j < numDIM; j++) {
      K(i, j)                     =  kg(i, j);
      K(i, j + numDOF2)           = -kg(i, j);
      K(i + numDOF2, j)           = -kg(i, j);
      K(i + numDOF2, j + numDOF2) =  kg(i, j);
    }
  }

  return *theMatrix;
}


const Matrix&
ExactTruss::getInitialStiff()
{

  // Material stiffness
  //
  // Get material tangent
  double EA = theSection->getTangent<1,section_layout>(State::Init)(0,0);

  Matrix3D kl{};
  kl(0, 0) = EA / Lo;

  // Compute R'*kl*R
  Matrix3D kg{};
  kg.addMatrixTripleProduct(0.0, R, kl, 1.0);

  Matrix& K = *theMatrix;
  K.Zero();

  // Copy stiffness into appropriate blocks in element stiffness
  int numDOF2 = numDOF / 2;
  for (int i = 0; i < numDIM; i++) {
    for (int j = 0; j < numDIM; j++) {
      K(i, j)                     =  kg(i, j);
      K(i, j + numDOF2)           = -kg(i, j);
      K(i + numDOF2, j)           = -kg(i, j);
      K(i + numDOF2, j + numDOF2) =  kg(i, j);
    }
  }

  return *theMatrix;
}


const Matrix&
ExactTruss::getDamp()
{
  if (doRayleighDamping == 1)
    return this->Element::getDamp();

  theMatrix->Zero();
  return *theMatrix;
}


const Matrix&
ExactTruss::getMass()
{
  // zero the matrix
  Matrix& mass = *theMatrix;
  mass.Zero();

  // check for quick return
  if (Lo == 0.0 || rho == 0.0)
    return mass;

  if (cMass == 0) {
    // lumped mass matrix
    double m    = 0.5 * rho * Lo;
    int numDOF2 = numDOF / 2;
    for (int i = 0; i < numDIM; i++) {
      mass(i, i)                     = m;
      mass(i + numDOF2, i + numDOF2) = m;
    }
  } else {
    // consistent mass matrix
    double m    = rho * Lo / 6.0;
    int numDOF2 = numDOF / 2;
    for (int i = 0; i < numDIM; i++) {
      mass(i, i)                     = 2.0 * m;
      mass(i, i + numDOF2)           = m;
      mass(i + numDOF2, i)           = m;
      mass(i + numDOF2, i + numDOF2) = 2.0 * m;
    }
  }

  return *theMatrix;
}

void
ExactTruss::zeroLoad()
{
  theLoad->Zero();
}


int
ExactTruss::addLoad(ElementalLoad* theLoad, double loadFactor)
{
  opserr << "Load type unknown for truss with tag: " << this->getTag()
         << "\n";
  return -1;
}


const Vector&
ExactTruss::getResistingForce()
{
  double Ln;
  Vector3D ql{};
  double J  = this->interpolate(ql,Ln);
  double N  = (1.0/J)*theSection->getResultant<1,section_layout>()[0];

  ql *= N*Lo;

  Vector3D qg{};
  qg.addMatrixTransposeVector(0.0, R, ql, 1.0);

  Vector& P = *theVector;
  P.Zero();

  // Copy forces into appropriate places
  int numDOF2 = numDOF / 2;
  for (int i = 0; i < numDIM; i++) {
    P(i)           = -qg(i);
    P(i + numDOF2) =  qg(i);
  }

  return *theVector;
}

const Vector&
ExactTruss::getResistingForceIncInertia()
{
  Vector& P = *theVector;
  P         = this->getResistingForce();

  // subtract external load
  P -= *theLoad;

  // now include the mass portion
  if (Lo != 0.0 && rho != 0.0) {

    // add inertia forces from element mass
    const Vector& accel1 = theNodes[0]->getTrialAccel();
    const Vector& accel2 = theNodes[1]->getTrialAccel();

    int numDOF2 = numDOF / 2;

    if (cMass == 0) {
      // lumped mass matrix
      double m = 0.5 * rho * Lo;
      for (int i = 0; i < numDIM; i++) {
        P(i) += m * accel1(i);
        P(i + numDOF2) += m * accel2(i);
      }
    } else {
      // consistent mass matrix
      double m = rho * Lo / 6.0;
      for (int i = 0; i < numDIM; i++) {
        (*theVector)(i) += 2.0 * m * accel1(i) + m * accel2(i);
        (*theVector)(i + numDOF2) += m * accel1(i) + 2.0 * m * accel2(i);
      }
    }

    // add the damping forces if rayleigh damping
    if (doRayleighDamping == 1 && (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0))
      theVector->addVector(1.0, this->getRayleighDampingForces(), 1.0);
  } else {
    // add the damping forces if rayleigh damping
    if (doRayleighDamping == 1 && (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0))
      theVector->addVector(1.0, this->getRayleighDampingForces(), 1.0);
  }

  return *theVector;
}


int
ExactTruss::setParameter(const char** argv, int argc, Parameter& param)
{
  if (argc < 1)
    return -1;

  // Mass density of the truss
  if (strcmp(argv[0], "rho") == 0)
    return param.addObject(2, this);

  // Explicit specification of a material parameter
  if (strstr(argv[0], "material") != 0 || strstr(argv[0], "section") != 0) {
    if (argc < 2)
      return -1;
    else
      return theSection->setParameter(&argv[1], argc - 1, param);
  }

  // Otherwise, send it to the material
  else
    return theSection->setParameter(argv, argc, param);
}


int
ExactTruss::updateParameter(int parameterID, Information& info)
{
  switch (parameterID) {
    case 2:  
      rho = info.theDouble; 
      return 0;
    default: 
      return -1;
  }
}

int
ExactTruss::activateParameter(int passedParameterID)
{
  parameterID = passedParameterID;
  return 0;
}



void
ExactTruss::Print(OPS_Stream& s, int flag)
{
  if (flag == OPS_PRINT_CURRENTSTATE) {
    s << "\nExactTruss, tag: " << this->getTag() << "\n";
    s << "\tConnected Nodes: " << connectedExternalNodes;
    s << "\tUndeformed Length: " << Lo << "\n";
    s << "\tMass Density/Length: " << rho << "\n";
    s << "\tConsistent Mass: " << cMass << "\n";

    if (theSection) {
      s << "\tSection, tag: " << theSection->getTag() << "\n";
      theSection->Print(s, flag);
    }
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_ELEM_INDENT << "{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"ExactTruss\", ";
    s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << "], ";
    s << "\"massperlength\": " << rho << ", ";
    s << "\"section\": \"" << theSection->getTag() << "\"}";
  }
}

double
ExactTruss::computeCurrentStrain()
{
  // NOTE method will not be called if Lo == 0

  Vector3D dX, dx, du;
  double Ln = this->locate(dX, dx, du);

  // return (Ln - Lo) / Lo;
  // double egr = du[0]/Lo;
  // egr += 0.5*(du[1]*du[1] + du[2]*du[2])/(Lo*Lo);
  double C   = Ln/Lo;
  double egr = 0.5*(C*C- 1.0);
  switch (strain_type) {
    case StrainType::Linear:
      return dX.dot(du)/(Lo*Lo);
    case StrainType::Green:
      return egr;
    case StrainType::Logarithmic:
      // return 0.5*std::log(C*C);
      return std::log(Ln/Lo);
    case StrainType::Corotational:
    default:
      return 2.0*Lo*egr/(Ln + Lo);
  }
}


Response*
ExactTruss::setResponse(const char** argv, int argc, OPS_Stream& output)
{
  Response* theResponse = 0;

  output.tag("ElementOutput");
  output.attr("eleType", "Truss");
  output.attr("eleTag", this->getTag());
  output.attr("node1", connectedExternalNodes[0]);
  output.attr("node2", connectedExternalNodes[1]);

  //
  // compare argv[0] for known response types
  //

  if ((strcmp(argv[0], "force") == 0) || 
      (strcmp(argv[0], "forces") == 0) ||
      (strcmp(argv[0], "globalForce") == 0) || 
      (strcmp(argv[0], "globalForces") == 0)) {
    char outputData[32];
    int numDOFperNode = numDOF / 2;
    for (int i = 0; i < numDOFperNode; i++) {
      sprintf(outputData, "P1_%d", i + 1);
      output.tag("ResponseType", outputData);
    }
    for (int j = 0; j < numDOFperNode; j++) {
      sprintf(outputData, "P2_%d", j + 1);
      output.tag("ResponseType", outputData);
    }
    theResponse = new ElementResponse(this, 1, Vector(numDOF));

  } else if ((strcmp(argv[0], "axialForce") == 0) || (strcmp(argv[0], "basicForce") == 0) ||
             (strcmp(argv[0], "basicForces") == 0)) {
    output.tag("ResponseType", "N");
    theResponse = new ElementResponse(this, 2, 0.0);

  } else if (strcmp(argv[0], "defo") == 0 || strcmp(argv[0], "deformation") == 0 ||
             strcmp(argv[0], "deformations") == 0 || strcmp(argv[0], "basicDefo") == 0 ||
             strcmp(argv[0], "basicDeformation") == 0 ||
             strcmp(argv[0], "basicDeformations") == 0) {

    output.tag("ResponseType", "U");
    theResponse = new ElementResponse(this, 3, 0.0);

  }
  // a section quantity
  else if (strcmp(argv[0], "section") == 0) {
    if (argc > 1) {
      // we need at least one more argument otherwise
      // there is no need to forward this call to the material
      // by default assume the old call style for backward compatibility "material result"
      int offset    = 1;
      bool is_valid = true;
      // in case the user specifies the gauss point id... "section 1 result"
      if (argc > 2) {
        int sectionNum = atoi(argv[1]);
        if (sectionNum == 1) {
          // this is the only supported gauss id
          offset = 2;
        } else if (sectionNum > 1) {
          // this is a number, but not within the valid range
          is_valid = false;
        }
        // if it is 0, then it is not a number, forward it as usual...
      }
      if (is_valid) {
        output.tag("GaussPointOutput");
        output.attr("number", 1);
        output.attr("eta", 0.0);
        theResponse = theSection->setResponse(&argv[offset], argc - offset, output);
        output.endTag();
      }
    }
  }

  output.endTag();
  return theResponse;
}


int
ExactTruss::getResponse(int responseID, Information& eleInfo)
{
  double strain, force;

  switch (responseID) {
  case 1: return eleInfo.setVector(this->getResistingForce());

  case 2:
    if (Lo == 0.0) {
      strain = 0.0;
      force  = 0.0;
    } else {

      int order      = theSection->getOrder();
      const ID& code = theSection->getType();

      const Vector& s = theSection->getStressResultant();
      force           = 0.0;
      for (int i = 0; i < order; i++) {
        if (code(i) == FrameStress::N)
          force += s(i);
      }
    }
    return eleInfo.setDouble(force);

  case 3:
    if (Lo == 0.0) {
      strain = 0.0;
    } else {
      strain = this->computeCurrentStrain();
    }
    return eleInfo.setDouble(Lo * strain);

  default:
    return -1;
  }
}
