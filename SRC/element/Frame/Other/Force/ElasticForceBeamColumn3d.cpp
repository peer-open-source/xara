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
/*
 * References
 *

  State Determination Algorithm
  ---
  Neuenhofer, A. and F. C. Filippou (1997). "Evaluation of Nonlinear Frame Finite
  Element Models." Journal of Structural Engineering, 123(7):958-966.

  Spacone, E., V. Ciampi, and F. C. Filippou (1996). "Mixed Formulation of
  Nonlinear Beam Finite Element." Computers and Structures, 58(1):71-83.


  Plastic Hinge Integration
  ---
  Scott, M. H. and G. L. Fenves (2006). "Plastic Hinge Integration Methods for
  Force-Based Beam-Column Elements." Journal of Structural Engineering,
  132(2):244-252.


  Analytical Response Sensitivity (DDM)
  ---
  Scott, M. H., P. Franchin, G. L. Fenves, and F. C. Filippou (2004).
  "Response Sensitivity for Nonlinear Beam-Column Elements."
  Journal of Structural Engineering, 130(9):1281-1288.


  Software Design
  ---
  Scott, M. H., G. L. Fenves, F. T. McKenna, and F. C. Filippou (2007).
  "Software Patterns for Nonlinear Beam-Column Models."
  Journal of Structural Engineering, Approved for publication, February 2007.

 *
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include <Information.h>
#include <Parameter.h>
#include <ElasticForceBeamColumn3d.h>
#include <interpolate/cbdi.h>
#include <Domain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <math.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>

Matrix ElasticForceBeamColumn3d::theMatrix(12,12);
Vector ElasticForceBeamColumn3d::theVector(12);
double ElasticForceBeamColumn3d::workArea[200];


// invoked by a FEM_ObjectBroker, recvSelf() needs to be invoked on this object.
ElasticForceBeamColumn3d::ElasticForceBeamColumn3d(): 
  Element(0,ELE_TAG_ElasticForceBeamColumn3d), connectedExternalNodes(2), 
  beamIntegr(0), numSections(0), crdTransf(0),
  rho(0.0), initialFlag(0), Se(NEBD),
  numEleLoads(0), sizeEleLoads(0), eleLoads(0), eleLoadFactors(0),
  parameterID(0)
{
  theNodes[0] = 0;  
  theNodes[1] = 0;

  for (int i = 0; i < maxNumSections; i++)
    sections[i] = 0;
}

// constructor which takes the unique element tag, sections,
// and the node ID's of it's nodal end points. 
// allocates the necessary space needed by each object
ElasticForceBeamColumn3d::ElasticForceBeamColumn3d (int tag, 
                                                    int nodeI, int nodeJ,
                                                    int numSec, 
                                                    SectionForceDeformation **sec,
                                                    BeamIntegration &bi,
                                                    CrdTransf &coordTransf,
                                                    double massDensPerUnitLength):
  Element(tag,ELE_TAG_ElasticForceBeamColumn3d), connectedExternalNodes(2),
  beamIntegr(0), numSections(numSec), crdTransf(0),
  rho(massDensPerUnitLength),
  initialFlag(0), Se(NEBD), 
  numEleLoads(0), sizeEleLoads(0), eleLoads(0), eleLoadFactors(0),
  parameterID(0)
{
  theNodes[0] = 0;
  theNodes[1] = 0;

  connectedExternalNodes(0) = nodeI;
  connectedExternalNodes(1) = nodeJ;    
  
  beamIntegr = bi.getCopy();
  if (beamIntegr == 0) {
    opserr << "Error: ElasticForceBeamColumn3d::ElasticForceBeamColumn3d: could not create copy of beam integration object" << endln;
  }
  
  // get copy of the transformation object   
  crdTransf = coordTransf.getCopy3d(); 
  if (crdTransf == 0) {
    opserr << "Error: ElasticForceBeamColumn3d::ElasticForceBeamColumn3d: could not create copy of coordinate transformation object" << endln;
  }

  if (numSections > maxNumSections) {
    opserr << "Error: ElasticForceBeamColumn3d::ElasticForceBeamColumn3d: numSections " << numSections << " exceeds max allowed, " << maxNumSections << endln;
    numSections = maxNumSections;
  }  

  int i;
  for (i = 0; i < numSections; i++) {
    sections[i] = sec[i]->getCopy();
    if (sections[i] == 0) {
      opserr << "Error: ElasticForceBeamColumn3d::ElasticForceBeamColumn3d: could not create copy of section object " << i << endln;
    }
  }
  for ( ; i < maxNumSections; i++)
    sections[i] = 0;
}

// ~ElasticForceBeamColumn3d():
//         destructor
//      delete must be invoked on any objects created by the object
ElasticForceBeamColumn3d::~ElasticForceBeamColumn3d()
{
  for (int i=0; i < numSections; i++)
    if (sections[i] != 0)
      delete sections[i];
  
  if (sizeEleLoads != 0) {
      if (eleLoads != 0)
          delete[] eleLoads;

      if (eleLoadFactors != 0)
          delete[] eleLoadFactors;
  }
  
  if (crdTransf != 0)
    delete crdTransf;

  if (beamIntegr != 0)
    delete beamIntegr;
}

int
ElasticForceBeamColumn3d::getNumExternalNodes(void) const
{
  return 2;
}

const ID &
ElasticForceBeamColumn3d::getExternalNodes(void) 
{
  return connectedExternalNodes;
}

Node **
ElasticForceBeamColumn3d::getNodePtrs()
{
  return theNodes;
}

int
ElasticForceBeamColumn3d::getNumDOF(void) 
{
  return NEGD;
}

void
ElasticForceBeamColumn3d::setDomain(Domain *theDomain)
{
  // check Domain is not null - invoked when object removed from a domain
  if (theDomain == 0) {
    theNodes[0] = 0;
    theNodes[1] = 0;
    
    opserr << "ElasticForceBeamColumn3d::setDomain:  theDomain = 0 ";
  }

  // get pointers to the nodes
  
  int Nd1 = connectedExternalNodes(0);  
  int Nd2 = connectedExternalNodes(1);
  
  theNodes[0] = theDomain->getNode(Nd1);
  theNodes[1] = theDomain->getNode(Nd2);  
  
  if (theNodes[0] == 0) {
    opserr << "ElasticForceBeamColumn3d::setDomain: Nd1: ";
    opserr << Nd1 << "does not exist in model\n";
  }
  
  if (theNodes[1] == 0) {
    opserr << "ElasticForceBeamColumn3d::setDomain: Nd2: ";
    opserr << Nd2 << "does not exist in model\n";
  }
  
  // call the DomainComponent class method 
  this->DomainComponent::setDomain(theDomain);
  
  // ensure connected nodes have correct number of dof's
  int dofNode1 = theNodes[0]->getNumberDOF();
  int dofNode2 = theNodes[1]->getNumberDOF();
  
  if ((dofNode1 != NND) || (dofNode2 != NND)) {
    opserr << "ElasticForceBeamColumn3d::setDomain(): Nd2 or Nd1 incorrect dof ";
  }
   
  // initialize the transformation
  if (crdTransf->initialize(theNodes[0], theNodes[1])) {
    opserr << "ElasticForceBeamColumn3d::setDomain(): Error initializing coordinate transformation";  
  }
    
  // get element length
  double L = crdTransf->getInitialLength();
  if (L == 0.0) {
    opserr << "ElasticForceBeamColumn3d::setDomain(): Zero element length:" << this->getTag();  
  }
}

int
ElasticForceBeamColumn3d::commitState()
{
  return crdTransf->commitState();
}

int ElasticForceBeamColumn3d::revertToLastCommit()
{
  return crdTransf->revertToLastCommit();
}

int ElasticForceBeamColumn3d::revertToStart()
{
  return crdTransf->revertToStart();
}


const Matrix &
ElasticForceBeamColumn3d::getInitialStiff(void)
{
  static Matrix f(NEBD, NEBD);   // element flexibility matrix  
  this->getInitialFlexibility(f);

  /*
  static Matrix I(NEBD,NEBD);   // an identity matrix for matrix inverse  
  I.Zero();
  for (int i=0; i<NEBD; i++)
    I(i,i) = 1.0;
  
  // calculate element stiffness matrix
  // invert3by3Matrix(f, kv);

  static Matrix kvInit(NEBD, NEBD);
  if (f.Solve(I, kvInit) < 0)
    opserr << "ElasticForceBeamColumn3d::getInitialStiff() -- could not invert flexibility\n";
  */

  static Matrix kvInit(NEBD, NEBD);
  f.Invert(kvInit);

  static Vector SeInit(NEBD);
  SeInit.Zero();
  return crdTransf->getGlobalStiffMatrix(kvInit, SeInit);
}

const Matrix &
ElasticForceBeamColumn3d::getTangentStiff(void)
{
  crdTransf->update();        // Will remove once we clean up the corotational 2d transformation -- MHS
  return this->getInitialStiff();
}
    
void
ElasticForceBeamColumn3d::computeReactions(double *p0)
{
    int type;
    double L = crdTransf->getInitialLength();

    for (int i = 0; i < numEleLoads; i++) {

        double loadFactor = eleLoadFactors[i];
        const Vector &data = eleLoads[i]->getData(type, loadFactor);

        if (type == LOAD_TAG_Beam3dUniformLoad) {
            double wy = data(0)*loadFactor;  // Transverse
            double wz = data(1)*loadFactor;  // Transverse
            double wa = data(2)*loadFactor;  // Axial

            p0[0] -= wa*L;
            double V = 0.5*wy*L;
            p0[1] -= V;
            p0[2] -= V;
            V = 0.5*wz*L;
            p0[3] -= V;
            p0[4] -= V;
        }
        else if (type == LOAD_TAG_Beam3dPointLoad) {
            double Py = data(0)*loadFactor;
            double Pz = data(1)*loadFactor;
            double N = data(2)*loadFactor;
            double aOverL = data(3);

            if (aOverL < 0.0 || aOverL > 1.0)
                continue;

            double V1 = Py*(1.0 - aOverL);
            double V2 = Py*aOverL;
            p0[0] -= N;
            p0[1] -= V1;
            p0[2] -= V2;
            V1 = Pz*(1.0 - aOverL);
            V2 = Pz*aOverL;
            p0[3] -= V1;
            p0[4] -= V2;
        }
    }
}

const Vector &
ElasticForceBeamColumn3d::getResistingForce(void)
{
  // Will remove once we clean up the corotational 2d transformation -- MHS
//crdTransf->update();
  
  double p0[NEBD];
  Vector p0Vec(p0, NEBD);
  p0Vec.Zero();

  if (numEleLoads > 0)
    this->computeReactions(p0);

  static Matrix f(NEBD, NEBD);   // element flexibility matrix  
  this->getInitialFlexibility(f);

  static Vector Se(NEBD);
  this->computeBasicForces(Se);

  return crdTransf->getGlobalResistingForce(Se, p0Vec);
}

void
ElasticForceBeamColumn3d::computeBasicForces(Vector &q)
{
  if (q.Size() != NEBD) {
    opserr << "ElasticFBC2d::computeBasicForces -- q size not 3" << endln;
    return;
  }

  static Matrix f(NEBD, NEBD);   // element flexibility matrix  
  this->getInitialFlexibility(f);

  const Vector &v = crdTransf->getBasicTrialDisp();
  f.Solve(v, q);
}


int
ElasticForceBeamColumn3d::update()
{
  int ok = crdTransf->update();
  
  static Vector q(NEBD);
  q.Zero();
  this->computeBasicForces(q);

  double L = crdTransf->getInitialLength();
  double oneOverL = 1.0/L;

  double xi[maxNumSections];
  beamIntegr->getSectionLocations(numSections, L, xi);

  for (int i = 0; i < numSections; i++) {

      int order      = sections[i]->getOrder();
      const ID &code = sections[i]->getType();

      double xL  = xi[i];
      double xL1 = xL-1.0;
      
      static Vector s;
      s.setData(workArea, order);
      static Vector e;
      e.setData(&workArea[order], order);

      int ii;
      for (ii = 0; ii < order; ii++) {
        switch(code(ii)) {
        case SECTION_RESPONSE_P:
          s(ii) = q(0);
          break;
        case SECTION_RESPONSE_MZ:
          s(ii) =  xL1*q(1) + xL*q(2);
          break;
        case SECTION_RESPONSE_VY:
          s(ii) = oneOverL*(q(1)+q(2));
          break;
        case SECTION_RESPONSE_MY:
          s(ii) =  xL1*q(3) + xL*q(4);
          break;
        case SECTION_RESPONSE_VZ:
          s(ii) = oneOverL*(q(3)+q(4));
          break;              
            case SECTION_RESPONSE_T:
                    s(ii) = q(5);
                    break;
            default:
          s(ii) = 0.0;
          break;
        }
      }
      
      // Add the effects of element loads, if present
      // s = b*q + sp
      if (numEleLoads > 0)
        this->computeSectionForces(s, i);

      const Matrix &fs = sections[i]->getInitialFlexibility();
      e.addMatrixVector(0.0, fs, s, 1.0);

      ok += sections[i]->setTrialSectionDeformation(e);
  }

  return ok;
}

const Matrix &
ElasticForceBeamColumn3d::getMass(void)
{ 
  theMatrix.Zero();
  
  double L = crdTransf->getInitialLength();
  if (rho != 0.0)
    theMatrix(0,0) = theMatrix(1,1) = theMatrix(2,2) =
      theMatrix(6,6) = theMatrix(7,7) = theMatrix(8,8) = 0.5*L*rho;
  
  return theMatrix;
}

void 
ElasticForceBeamColumn3d::zeroLoad(void)
{
  // This is a semi-hack -- MHS
  numEleLoads = 0;

  return;
}

int
ElasticForceBeamColumn3d::addLoad(ElementalLoad *theLoad, double loadFactor)
{
    if (numEleLoads == sizeEleLoads) {

        //
        // create larger arrays, copy old, delete old & set as new
        //

        ElementalLoad ** theNextEleLoads = new ElementalLoad *[sizeEleLoads + 1];
        double *theNextEleLoadFactors = new double[sizeEleLoads + 1];
        for (int i = 0; i<numEleLoads; i++) {
            theNextEleLoads[i] = eleLoads[i];
            theNextEleLoadFactors[i] = eleLoadFactors[i];
        }

        delete[] eleLoads;
        delete[] eleLoadFactors;
        eleLoads = theNextEleLoads;
        eleLoadFactors = theNextEleLoadFactors;

        // increment array size
        sizeEleLoads += 1;
    }

    eleLoadFactors[numEleLoads] = loadFactor;
    eleLoads[numEleLoads] = theLoad;
    numEleLoads++;

    return 0;
}

void
ElasticForceBeamColumn3d::computeSectionForces(Vector &sp, int isec)
{
    int type;

    double L = crdTransf->getInitialLength();

    double xi[maxNumSections];
    beamIntegr->getSectionLocations(numSections, L, xi);
    double x = xi[isec] * L;

    int order = sections[isec]->getOrder();
    const ID &code = sections[isec]->getType();

    for (int i = 0; i < numEleLoads; i++) {

        double loadFactor = eleLoadFactors[i];
        const Vector &data = eleLoads[i]->getData(type, loadFactor);

        if (type == LOAD_TAG_Beam3dUniformLoad) {
            double wy = data(0)*loadFactor;  // Transverse
            double wz = data(1)*loadFactor;  // Transverse
            double wa = data(2)*loadFactor;  // Axial

            for (int ii = 0; ii < order; ii++) {

                switch (code(ii)) {
                case SECTION_RESPONSE_P:
                    sp(ii) += wa*(L - x);
                    break;
                case SECTION_RESPONSE_MZ:
                    sp(ii) += wy*0.5*x*(x - L);
                    break;
                case SECTION_RESPONSE_VY:
                    sp(ii) += wy*(x - 0.5*L);
                    break;
                case SECTION_RESPONSE_MY:
                    sp(ii) += wz*0.5*x*(L - x);
                    break;
                case SECTION_RESPONSE_VZ:
                    sp(ii) += wz*(x - 0.5*L);
                    break;
                default:
                    break;
                }
            }
        }
        else if (type == LOAD_TAG_Beam3dPointLoad) {
            double Py = data(0)*loadFactor;
            double Pz = data(1)*loadFactor;
            double N = data(2)*loadFactor;
            double aOverL = data(3);

            if (aOverL < 0.0 || aOverL > 1.0)
                continue;

            double a = aOverL*L;

            double Vy1 = Py*(1.0 - aOverL);
            double Vy2 = Py*aOverL;

            double Vz1 = Pz*(1.0 - aOverL);
            double Vz2 = Pz*aOverL;

            for (int ii = 0; ii < order; ii++) {

                if (x <= a) {
                    switch (code(ii)) {
                    case SECTION_RESPONSE_P:
                        sp(ii) += N;
                        break;
                    case SECTION_RESPONSE_MZ:
                        sp(ii) -= x*Vy1;
                        break;
                    case SECTION_RESPONSE_VY:
                        sp(ii) -= Vy1;
                        break;
                    case SECTION_RESPONSE_MY:
                        sp(ii) += x*Vz1;
                        break;
                    case SECTION_RESPONSE_VZ:
                        sp(ii) -= Vz1;
                        break;
                    default:
                        break;
                    }
                }
                else {
                    switch (code(ii)) {
                    case SECTION_RESPONSE_MZ:
                        sp(ii) -= (L - x)*Vy2;
                        break;
                    case SECTION_RESPONSE_VY:
                        sp(ii) += Vy2;
                        break;
                    case SECTION_RESPONSE_MY:
                        sp(ii) += (L - x)*Vz2;
                        break;
                    case SECTION_RESPONSE_VZ:
                        sp(ii) += Vz2;
                        break;
                    default:
                        break;
                    }
                }
            }
        }
        else {
            opserr << "ForceBeamColumn3d::addLoad -- load type unknown for element with tag: " <<
                this->getTag() << endln;
        }
    }

    // Don't think we need to do this anymore -- MHS
    //this->update(); // quick fix -- MHS
}

int 
ElasticForceBeamColumn3d::addInertiaLoadToUnbalance(const Vector &accel)
{
  // Check for a quick return
  if (rho == 0.0)
    return 0;

  // get R * accel from the nodes
  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);    

  double L = crdTransf->getInitialLength();
  double m = 0.5*rho*L;

  // Should be done through p0[0]
  /*
  load(0) -= m*Raccel1(0);
  load(1) -= m*Raccel1(1);
  load(3) -= m*Raccel2(0);
  load(4) -= m*Raccel2(1);
  */

  return 0;
}

const Vector &
ElasticForceBeamColumn3d::getResistingForceIncInertia()
{        
  // Compute the current resisting force
  theVector = this->getResistingForce();

  // Check for a quick return
  if (rho != 0.0) {
    const Vector &accel1 = theNodes[0]->getTrialAccel();
    const Vector &accel2 = theNodes[1]->getTrialAccel();
    
    double L = crdTransf->getInitialLength();
    double m = 0.5*rho*L;
    
    theVector(0) += m*accel1(0);
    theVector(1) += m*accel1(1);
        theVector(2) += m*accel1(2);
    theVector(6) += m*accel2(0);
    theVector(7) += m*accel2(1);
    theVector(8) += m*accel2(2);

    // add the damping forces if rayleigh damping
    if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
      theVector += this->getRayleighDampingForces();

  } else {
    // add the damping forces if rayleigh damping
    if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
      theVector += this->getRayleighDampingForces();
  }

  return theVector;
}

int
ElasticForceBeamColumn3d::sendSelf(int commitTag, Channel &theChannel)
{  
  return -1;
}    

int
ElasticForceBeamColumn3d::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return -1;
}

int
ElasticForceBeamColumn3d::getInitialFlexibility(Matrix &fe)
{
  fe.Zero();
  
  double L = crdTransf->getInitialLength();
  double oneOverL  = 1.0/L;  
  
  // Flexibility from elastic interior
  beamIntegr->addElasticFlexibility(L, fe);
  
  double xi[maxNumSections];
  beamIntegr->getSectionLocations(numSections, L, xi);
  
  double wt[maxNumSections];
  beamIntegr->getSectionWeights(numSections, L, wt);
  
  for (int i = 0; i < numSections; i++) {
    
    int order      = sections[i]->getOrder();
    const ID &code = sections[i]->getType();
    
    Matrix fb(workArea, order, NEBD);
    
    double xL  = xi[i];
    double xL1 = xL-1.0;
    double wtL = wt[i]*L;

    const Matrix &fSec = sections[i]->getInitialFlexibility();
    fb.Zero();
    double tmp;
    int ii, jj;
    for (ii = 0; ii < order; ii++) {
      switch(code(ii)) {
      case SECTION_RESPONSE_P:
        for (jj = 0; jj < order; jj++)
          fb(jj,0) += fSec(jj,ii)*wtL;
        break;
      case SECTION_RESPONSE_MZ:
        for (jj = 0; jj < order; jj++) {
          tmp = fSec(jj,ii)*wtL;
          fb(jj,1) += xL1*tmp;
          fb(jj,2) += xL*tmp;
        }
        break;
      case SECTION_RESPONSE_VY:
        for (jj = 0; jj < order; jj++) {
          tmp = oneOverL*fSec(jj,ii)*wtL;
          fb(jj,1) += tmp;
          fb(jj,2) += tmp;
        }
        break;
      case SECTION_RESPONSE_MY:
        for (jj = 0; jj < order; jj++) {
          tmp = fSec(jj,ii)*wtL;
          fb(jj,3) += xL1*tmp;
          fb(jj,4) += xL*tmp;
        }
        break;
      case SECTION_RESPONSE_VZ:
        for (jj = 0; jj < order; jj++) {
          tmp = oneOverL*fSec(jj,ii)*wtL;
          fb(jj,3) += tmp;
          fb(jj,4) += tmp;
        }
        break;
          case SECTION_RESPONSE_T:
        for (jj = 0; jj < order; jj++)
          fb(jj,5) += fSec(jj,ii)*wtL;
        break;
      default:
        break;
      }
    }
    for (ii = 0; ii < order; ii++) {
      switch (code(ii)) {
      case SECTION_RESPONSE_P:
        for (jj = 0; jj < NEBD; jj++)
          fe(0,jj) += fb(ii,jj);
        break;
      case SECTION_RESPONSE_MZ:
        for (jj = 0; jj < NEBD; jj++) {
          tmp = fb(ii,jj);
          fe(1,jj) += xL1*tmp;
          fe(2,jj) += xL*tmp;
        }
        break;
      case SECTION_RESPONSE_VY:
        for (jj = 0; jj < NEBD; jj++) {
          tmp = oneOverL*fb(ii,jj);
          fe(1,jj) += tmp;
          fe(2,jj) += tmp;
        }
        break;
      case SECTION_RESPONSE_MY:
        for (jj = 0; jj < NEBD; jj++) {
          tmp = fb(ii,jj);
          fe(3,jj) += xL1*tmp;
          fe(4,jj) += xL*tmp;
        }
        break;
      case SECTION_RESPONSE_VZ:
        for (jj = 0; jj < NEBD; jj++) {
          tmp = oneOverL*fb(ii,jj);
          fe(3,jj) += tmp;
          fe(4,jj) += tmp;
        }
        break;
      case SECTION_RESPONSE_T:
        for (jj = 0; jj < NEBD; jj++)
          fe(5,jj) += fb(ii,jj);
        break;
      default:
        break;
      }
    }
  }
  
  return 0;
}


void
ElasticForceBeamColumn3d::Print(OPS_Stream &s, int flag)
  {
  static Vector Se(NEBD);
  static Vector vp(NEBD);
  static Matrix fe(NEBD,NEBD);
  double p0[NEBD];
  Vector p0Vec(p0,NEBD);
  p0Vec.Zero();

    // flags with negative values are used by GSA
    if (flag == -1) { 
      int eleTag = this->getTag();
      s << "EL_BEAM\t" << eleTag << "\t";
      s << sections[0]->getTag() << "\t" << sections[numSections-1]->getTag(); 
      s  << "\t" << connectedExternalNodes(0) << "\t" << connectedExternalNodes(1);
      s << "\t0\t0.0000000\n";
    }  

    // flags with negative values are used by GSA  
    else if (flag < -1) {
      int eleTag = this->getTag();
      int counter = (flag +1) * -1;
          this->computeBasicForces(Se);
      double P  = Se(0);
      double MZ1 = Se(1);
      double MZ2 = Se(2);
      double MY1 = Se(3);
      double MY2 = Se(4);
      double L = crdTransf->getInitialLength();
      double VY = (MZ1+MZ2)/L;
      theVector(1) =  VY;
      theVector(4) = -VY;
      double VZ = (MY1+MY2)/L;
      double T  = Se(5);

      s << "FORCE\t" << eleTag << "\t" << counter << "\t0";
      s << "\t" << -P+p0[0] << "\t"  <<  VY+p0[1] << "\t"  << -VZ+p0[3]  << endln;
      s << "FORCE\t" << eleTag << "\t" << counter << "\t1";
      s << "\t"  << P  << ' '  << -VY+p0[2] << ' ' << VZ+p0[4] << endln;
      s << "MOMENT\t" << eleTag << "\t" << counter << "\t0";
      s << "\t" << -T << "\t"  << MY1 << "\t" << MZ1 << endln;
      s << "MOMENT\t" << eleTag << "\t" << counter << "\t1";
      s << "\t" << T << ' ' << MY2 << ' '  <<  MZ2 << endln;
    }

    // flag set to 2 used to print everything .. used for viewing data for UCSD renderer  
    else if (flag == 2) {
       static Vector xAxis(3);
       static Vector yAxis(3);
       static Vector zAxis(3);

       crdTransf->getLocalAxes(xAxis, yAxis, zAxis);

       s << "#ForceBeamColumn3D\n";
       s << "#LocalAxis " << xAxis(0) << " " << xAxis(1) << " " << xAxis(2);
       s << " " << yAxis(0) << " " << yAxis(1) << " " << yAxis(2) << " ";
       s << zAxis(0) << " " << zAxis(1) << " " << zAxis(2) << endln;

       const Vector &node1Crd = theNodes[0]->getCrds();
       const Vector &node2Crd = theNodes[1]->getCrds();        
       const Vector &node1Disp = theNodes[0]->getDisp();
       const Vector &node2Disp = theNodes[1]->getDisp();    

       s << "#NODE " << node1Crd(0) << " " << node1Crd(1) << " " << node1Crd(2)
         << " " << node1Disp(0) << " " << node1Disp(1) << " " << node1Disp(2)
         << " " << node1Disp(3) << " " << node1Disp(4) << " " << node1Disp(5) << endln;

       s << "#NODE " << node2Crd(0) << " " << node2Crd(1) << " " << node2Crd(2)
         << " " << node2Disp(0) << " " << node2Disp(1) << " " << node2Disp(2)
         << " " << node2Disp(3) << " " << node2Disp(4) << " " << node2Disp(5) << endln;

           this->computeBasicForces(Se);
       double P  = Se(0);
       double MZ1 = Se(1);
       double MZ2 = Se(2);
       double MY1 = Se(3);
       double MY2 = Se(4);
       double L = crdTransf->getInitialLength();
       double VY = (MZ1+MZ2)/L;
       theVector(1) =  VY;
       theVector(4) = -VY;
       double VZ = (MY1+MY2)/L;
       double T  = Se(5);
       s << "#END_FORCES " << -P+p0[0] << ' '  <<  VY+p0[1] << ' '  << -VZ+p0[3] << ' ' 
         << -T << ' '  << MY1 << ' ' << MZ1 << endln;
       s << "#END_FORCES "  << P  << ' '  << -VY+p0[2] << ' ' << VZ+p0[4] << ' '  
         << T << ' ' << MY2 << ' '  <<  MZ2 << endln;

       // plastic hinge rotation
       this->getInitialFlexibility(fe);
       vp = crdTransf->getBasicTrialDisp();
       vp.addMatrixVector(1.0, fe, Se, -1.0);
       s << "#PLASTIC_HINGE_ROTATION " << vp[1] << " " << vp[2] << " " << vp[3] << " " << vp[4] 
         << " " << 0.1*L << " " << 0.1*L << endln;
/*
       // allocate array of vectors to store section coordinates and displacements
       static int maxNumSections = 0;
       static Vector *coords = 0;
       static Vector *displs = 0;
       if (maxNumSections < numSections) {
         if (coords != 0) 
           delete [] coords;
         if (displs != 0)
           delete [] displs;

         coords = new Vector [numSections];
         displs = new Vector [numSections];

         if (!coords) {
           opserr << "ForceBeamColumn3d::Print() -- failed to allocate coords array";   
           exit(-1);
         }

         int i;
         for (i = 0; i < numSections; i++)
           coords[i] = Vector(NDM);

         if (!displs) {
           opserr << "ForceBeamColumn3d::Print() -- failed to allocate coords array";   
           exit(-1);
         }

         for (i = 0; i < numSections; i++)
           displs[i] = Vector(NDM);

         maxNumSections = numSections;
       }

       // compute section location & displacements
       this->compSectionDisplacements(coords, displs);

       // spit out the section location & invoke print on the scetion
       for (int i=0; i<numSections; i++) {
         s << "#SECTION " << (coords[i])(0) << " " << (coords[i])(1) << " " << (coords[i])(2);       
         s << " " << (displs[i])(0) << " " << (displs[i])(1) << " " << (displs[i])(2) << endln;
         sections[i]->Print(s, flag); 
       }
           */
    }

        if (flag == OPS_PRINT_CURRENTSTATE) {
       s << "\nElement: " << this->getTag() << " Type: ElasticForceBeamColumn3d ";
       s << "\tConnected Nodes: " << connectedExternalNodes ;
       s << "\tNumber of Sections: " << numSections;
       s << "\tMass density: " << rho << endln;
       beamIntegr->Print(s, flag);
           this->computeBasicForces(Se);
       double P  = Se(0);
       double MZ1 = Se(1);
       double MZ2 = Se(2);
       double MY1 = Se(3);
       double MY2 = Se(4);
       double L = crdTransf->getInitialLength();
       double VY = (MZ1+MZ2)/L;
       theVector(1) =  VY;
       theVector(4) = -VY;
       double VZ = (MY1+MY2)/L;
       double T  = Se(5);
       s << "\tEnd 1 Forces (P MZ VY MY VZ T): "
         << -P+p0[0] << " " << MZ1 << " " <<  VY+p0[1] << " " 
         << MY1 << " " << -VZ+p0[3] << " " << T << endln;
       s << "\tEnd 2 Forces (P MZ VY MY VZ T): "
         << P        << " " << MZ2 << " " << -VY+p0[2] << " " 
         << MY2 << " " <<  VZ+p0[4] << " " << -T << endln;

       if (flag == 1) { 
         for (int i = 0; i < numSections; i++)
           s << "\numSections "<<i<<" :" << *sections[i];
       }
    }
        
        if (flag == OPS_PRINT_PRINTMODEL_JSON) {
                s << "\t\t\t{";
                s << "\"name\": " << this->getTag() << ", ";
                s << "\"type\": \"ElasticForceBeamColumn3d\", ";
                s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << "], ";
                s << "\"sections\": [";
                for (int i = 0; i < numSections - 1; i++)
                        s << "\"" << sections[i]->getTag() << "\", ";
                s << "\"" << sections[numSections - 1]->getTag() << "\"], ";
                s << "\"integration\": ";
                beamIntegr->Print(s, flag);
                s << ", \"massperlength\": " << rho << ", ";
                s << "\"crdTransformation\": \"" << crdTransf->getTag() << "\"}";
        }
  }

  OPS_Stream &operator<<(OPS_Stream &s, ElasticForceBeamColumn3d &E)
  {
    E.Print(s);
    return s;
  }

  Response*
  ElasticForceBeamColumn3d::setResponse(const char **argv, int argc, OPS_Stream &output)
  {

    Response *theResponse = 0;

    output.tag("ElementOutput");
    output.attr("eleType","ElasticForceBeamColumn3d");
    output.attr("eleTag",this->getTag());
    output.attr("node1",connectedExternalNodes[0]);
    output.attr("node2",connectedExternalNodes[1]);

    //
    // we compare argv[0] for known response types 
    //

    // global force - 
    if (strcmp(argv[0],"forces") == 0 || strcmp(argv[0],"force") == 0
        || strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0) {

      output.tag("ResponseType","Px_1");
      output.tag("ResponseType","Py_1");
      output.tag("ResponseType","Pz_1");
      output.tag("ResponseType","Mx_1");
      output.tag("ResponseType","My_1");
      output.tag("ResponseType","Mz_1");
      output.tag("ResponseType","Px_2");
      output.tag("ResponseType","Py_2");
      output.tag("ResponseType","Pz_2");
      output.tag("ResponseType","Mx_2");
      output.tag("ResponseType","My_2");
      output.tag("ResponseType","Mz_2");


      theResponse = new ElementResponse(this, 1, theVector);

    // local force -
    }  else if (strcmp(argv[0],"localForce") == 0 || strcmp(argv[0],"localForces") == 0) {

      output.tag("ResponseType","N_1");
      output.tag("ResponseType","Vy_1");
      output.tag("ResponseType","Vz_1");
      output.tag("ResponseType","T_1");
      output.tag("ResponseType","My_1");
      output.tag("ResponseType","Mz_1");
      output.tag("ResponseType","N_2");
      output.tag("ResponseType","Vy_2");
      output.tag("ResponseType","Vz_2");
      output.tag("ResponseType","T_2");
      output.tag("ResponseType","My_2");
      output.tag("ResponseType","Mz_2");

      theResponse = new ElementResponse(this, 2, theVector);

    // chord rotation -
    }  else if (strcmp(argv[0],"chordRotation") == 0 || strcmp(argv[0],"chordDeformation") == 0 
              || strcmp(argv[0],"basicDeformation") == 0) {

      output.tag("ResponseType","eps");
      output.tag("ResponseType","thetaZ_1");
      output.tag("ResponseType","thetaZ_2");
      output.tag("ResponseType","thetaY_1");
      output.tag("ResponseType","thetaY_2");
      output.tag("ResponseType","thetaX");

      theResponse = new ElementResponse(this, 3, Vector(6));

    // plastic rotation -
    } else if (strcmp(argv[0],"plasticRotation") == 0 || strcmp(argv[0],"plasticDeformation") == 0) {

    output.tag("ResponseType","epsP");
    output.tag("ResponseType","thetaZP_1");
    output.tag("ResponseType","thetaZP_2");
    output.tag("ResponseType","thetaYP_1");
    output.tag("ResponseType","thetaYP_2");
    output.tag("ResponseType","thetaXP");

    theResponse = new ElementResponse(this, 4, Vector(6));
  
  // point of inflection
  } else if (strcmp(argv[0],"inflectionPoint") == 0) {
    theResponse = new ElementResponse(this, 5, Vector(2));
  
  // tangent drift
  } else if (strcmp(argv[0],"tangentDrift") == 0) {
    theResponse = new ElementResponse(this, 6, Vector(4));
  }
  // section response -
  else if (strstr(argv[0],"sectionX") != 0) {
    if (argc > 2) {
      float sectionLoc = atof(argv[1]);

      double xi[maxNumSections];
      double L = crdTransf->getInitialLength();
      beamIntegr->getSectionLocations(numSections, L, xi);
      
      sectionLoc /= L;

      float minDistance = fabs(xi[0]-sectionLoc);
      int sectionNum = 0;
      for (int i = 1; i < numSections; i++) {
        if (fabs(xi[i]-sectionLoc) < minDistance) {
          minDistance = fabs(xi[i]-sectionLoc);
          sectionNum = i;
        }
      }

      output.tag("GaussPointOutput");
      output.attr("number",sectionNum+1);
      output.attr("eta",xi[sectionNum]*L);
      
      if (strcmp(argv[2],"dsdh") != 0) {
        theResponse = sections[sectionNum]->setResponse(&argv[2], argc-2, output);
      } else {
        int order = sections[sectionNum]->getOrder();
        theResponse = new ElementResponse(this, 76, Vector(order));
        Information &info = theResponse->getInformation();
        info.theInt = sectionNum;
          }

    }
  }

  // section response -
  else if (strstr(argv[0],"section") != 0) { 
    if (argc > 2) {
    
      int sectionNum = atoi(argv[1]);
      if (sectionNum > 0 && sectionNum <= numSections) {
        double xi[maxNumSections];
        double L = crdTransf->getInitialLength();
        beamIntegr->getSectionLocations(numSections, L, xi);
        
        output.tag("GaussPointOutput");
        output.attr("number",sectionNum);
        output.attr("eta",2.0*xi[sectionNum-1]-1.0);
        theResponse =  sections[sectionNum-1]->setResponse(&argv[2], argc-2, output);
        
        if (strcmp(argv[2],"dsdh") != 0) {
          theResponse = sections[sectionNum-1]->setResponse(&argv[2], argc-2, output);
        } else {
          int order = sections[sectionNum-1]->getOrder();
          theResponse = new ElementResponse(this, 76, Vector(order));
          Information &info = theResponse->getInformation();
          info.theInt = sectionNum;
        }

      }
    }
  }
  
  output.endTag();
  return theResponse;
}

int 
ElasticForceBeamColumn3d::getResponse(int responseID, Information &eleInfo)
{
  static Vector vp(NEBD);
  static Matrix fe(NEBD,NEBD);
  static Vector Se(NEBD);
  this->computeBasicForces(Se);

  double p0[NEBD];
  Vector p0Vec(p0, NEBD);
  p0Vec.Zero();

  if (responseID == 1)
    return eleInfo.setVector(this->getResistingForce());
  
  else if (responseID == 2) {
    // Axial
    double N = Se(0);
    theVector(6) =  N;
    theVector(0) = -N+p0[0];
    
    // Torsion
    double T = Se(5);
    theVector(9) =  T;
    theVector(3) = -T;
    
    // Moments about z and shears along y
    double M1 = Se(1);
    double M2 = Se(2);
    theVector(5)  = M1;
    theVector(11) = M2;
    double L = crdTransf->getInitialLength();
    double V = (M1+M2)/L;
    theVector(1) =  V+p0[1];
    theVector(7) = -V+p0[2];
    
    // Moments about y and shears along z
    M1 = Se(3);
    M2 = Se(4);
    theVector(4)  = M1;
    theVector(10) = M2;
    V = (M1+M2)/L;
    theVector(2) = -V+p0[3];
    theVector(8) =  V+p0[4];
      
    return eleInfo.setVector(theVector);

  }
      
  // Chord rotation
  else if (responseID == 3) {
    vp = crdTransf->getBasicTrialDisp();
    return eleInfo.setVector(vp);
  }

  // Plastic rotation
  else if (responseID == 4) {
    this->getInitialFlexibility(fe);
    vp = crdTransf->getBasicTrialDisp();
    vp.addMatrixVector(1.0, fe, Se, -1.0);
    return eleInfo.setVector(vp);
  }

  // Point of inflection
  else if (responseID == 5) {
    static Vector LI(2);
    LI(0) = 0.0;
    LI(1) = 0.0;

    double L = crdTransf->getInitialLength();

    if (fabs(Se(1)+Se(2)) > DBL_EPSILON)
      LI(0) = Se(1)/(Se(1)+Se(2))*L;

    if (fabs(Se(3)+Se(4)) > DBL_EPSILON)
      LI(1) = Se(3)/(Se(3)+Se(4))*L;

    return eleInfo.setVector(LI);
  }
/*
  // Tangent drift
  else if (responseID == 6) {
    double d2z = 0.0;
    double d2y = 0.0;
    double d3z = 0.0;
    double d3y = 0.0;

    double L = crdTransf->getInitialLength();

    double wts[maxNumSections];
    beamIntegr->getSectionWeights(numSections, L, wts);

    double pts[maxNumSections];
    beamIntegr->getSectionLocations(numSections, L, pts);

    // Location of inflection point from node I
    double LIz = 0.0;
    if (fabs(Se(1)+Se(2)) > DBL_EPSILON)
      LIz = Se(1)/(Se(1)+Se(2))*L;

    double LIy = 0.0;
    if (fabs(Se(3)+Se(4)) > DBL_EPSILON)
      LIy = Se(3)/(Se(3)+Se(4))*L;

    int i;
    for (i = 0; i < numSections; i++) {
      double x = pts[i]*L;
      const ID &type = sections[i]->getType();
      int order = sections[i]->getOrder();
      double kappa = 0.0;
      if (x < LIz) {
        for (int j = 0; j < order; j++)
          if (type(j) == SECTION_RESPONSE_MZ)
            kappa += vs[i](j);
        double b = -LIz+x;
        d2z += (wts[i]*L)*kappa*b;
      }
      kappa = 0.0;
      if (x < LIy) {
        for (int j = 0; j < order; j++)
          if (type(j) == SECTION_RESPONSE_MY)
            kappa += vs[i](j);
        double b = -LIy+x;
        d2y += (wts[i]*L)*kappa*b;
      }
    }

    d2z += beamIntegr->getTangentDriftI(L, LIz, Se(1), Se(2));
    d2y += beamIntegr->getTangentDriftI(L, LIy, Se(3), Se(4), true);

    for (i = numSections-1; i >= 0; i--) {
      double x = pts[i]*L;
      const ID &type = sections[i]->getType();
      int order = sections[i]->getOrder();
      double kappa = 0.0;
      if (x > LIz) {
        for (int j = 0; j < order; j++)
          if (type(j) == SECTION_RESPONSE_MZ)
            kappa += vs[i](j);
        double b = x-LIz;
        d3z += (wts[i]*L)*kappa*b;
      }
      kappa = 0.0;
      if (x > LIy) {
        for (int j = 0; j < order; j++)
          if (type(j) == SECTION_RESPONSE_MY)
            kappa += vs[i](j);
        double b = x-LIy;
        d3y += (wts[i]*L)*kappa*b;
      }
    }

    d3z += beamIntegr->getTangentDriftJ(L, LIz, Se(1), Se(2));
    d3y += beamIntegr->getTangentDriftJ(L, LIy, Se(3), Se(4), true);

    static Vector d(4);
    d(0) = d2z;
    d(1) = d3z;
    d(2) = d2y;
    d(3) = d3y;

    return eleInfo.setVector(d);
  }
*/
  else
    return -1;
}


int
ElasticForceBeamColumn3d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return 0;

  if (strcmp(argv[0],"rho") == 0)
    return param.addObject(1, this);
  
  // section response -
  else if (strstr(argv[0],"sectionX") != 0) {
    if (argc > 2) {
      float sectionLoc = atof(argv[1]);

      double xi[maxNumSections];
      double L = crdTransf->getInitialLength();
      beamIntegr->getSectionLocations(numSections, L, xi);
      
      sectionLoc /= L;

      float minDistance = fabs(xi[0]-sectionLoc);
      int sectionNum = 0;
      for (int i = 1; i < numSections; i++) {
        if (fabs(xi[i]-sectionLoc) < minDistance) {
          minDistance = fabs(xi[i]-sectionLoc);
          sectionNum = i;
        }
      }

      return sections[sectionNum]->setParameter(&argv[2], argc-2, param);
    }
  }

  // If the parameter belongs to a section or lower
  else if (strstr(argv[0],"section") != 0) {
    
    if (argc < 3)
      return 0;

    // Get section number: 1...Np
    int sectionNum = atoi(argv[1]);
    
    if (sectionNum > 0 && sectionNum <= numSections)
      return sections[sectionNum-1]->setParameter(&argv[2], argc-2, param);

    else
      return 0;

    /*
    // Find the right section and call its setParameter method
    int ok = 0;
    for (int i = 0; i < numSections; i++)
      if (paramSectionTag == sections[i]->getTag())
        ok += sections[i]->setParameter(&argv[2], argc-2, param);

    return ok;
    */
  }
  
  else if (strstr(argv[0],"integration") != 0) {
    
    if (argc < 2)
      return 0;

    return beamIntegr->setParameter(&argv[1], argc-1, param);
  }

  // Default, send to everything
  else {
    int ok = 0;
    for (int i = 0; i < numSections; i++)
      ok += sections[i]->setParameter(argv, argc, param);
    ok += beamIntegr->setParameter(argv, argc, param);
    return ok;
  }
  return 0;
}

int
ElasticForceBeamColumn3d::updateParameter (int parameterID, Information &info)
{
  if (parameterID == 1) {
    rho = info.theDouble;
    return 0;
  }
  else
    return -1;
}

int
ElasticForceBeamColumn3d::activateParameter(int passedParameterID)
{
  parameterID = passedParameterID;

  return 0;  
}
