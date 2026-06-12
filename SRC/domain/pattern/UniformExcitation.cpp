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
// Purpose: This file contains the class definition for UniformExcitation.
// UniformExcitation is an abstract class.
//
// Written: fmk 11/98
// Revised:
//
#include <cassert>
#include <UniformExcitation.h>
#include <GroundMotion.h>
#include <Domain.h>
#include <NodeIter.h>
#include <Node.h>
#include <ElementIter.h>
#include <Element.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <SP_ConstraintIter.h>
#include <SP_Constraint.h>

#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <DOF_Group.h>

#include <VectorND.h>
#include <MatrixND.h>

using namespace OpenSees;


UniformExcitation::UniformExcitation()
 : LoadPattern(0, PATTERN_TAG_UniformExcitation, 1.0)
 , theMotion(0), theDof(0), vel0(0.0), fact(0.0)
 , parameterID(0)
 , currentTime(0.0)
{

}


UniformExcitation::UniformExcitation(GroundMotion *motion, 
                                     int ndm, 
                                     int dof, 
                                     int tag,
                                     double velZero, double theFactor)
 : LoadPattern(tag, PATTERN_TAG_UniformExcitation, 1.0)
 , theMotion(motion)
 , theDof(dof), vel0(velZero), fact(theFactor)
 , parameterID(0)
 , currentTime(0.0)
{

}


UniformExcitation::~UniformExcitation()
{
  delete theMotion;
}


const GroundMotion *
UniformExcitation::getGroundMotion()
{
  return theMotion;
}

int
UniformExcitation::setParameter(const char **argv, int argc, Parameter &param)
{
  return theMotion->setParameter(argv, argc, param);
}


int
UniformExcitation::updateParameter(int parameterID, Information &info)
{
  return theMotion->updateParameter(parameterID, info);
}

int
UniformExcitation::activateParameter(int pparameterID)
{
  parameterID = pparameterID;
  return theMotion->activateParameter(pparameterID);
}


void
UniformExcitation::setDomain(Domain *theDomain) 
{
  this->LoadPattern::setDomain(theDomain);

  // now we go through and set all the node velocities to be vel0 
  // for those nodes not fixed in the dirn!
  if (vel0 != 0.0) {

    SP_ConstraintIter &theSPs = theDomain->getSPs();
    SP_Constraint *theSP;
    ID constrainedNodes(0);
    int count = 0;
    while ((theSP=theSPs()) != nullptr) {
      if (theSP->getDOF_Number() == theDof) {
        constrainedNodes[count] = theSP->getNodeTag();
        count++;
      }
    }


    NodeIter &theNodes = theDomain->getNodes();
    Node *theNode;
    Vector newVel(1);
    int currentSize = 1;
    while ((theNode = theNodes()) != 0) {
      int tag = theNode->getTag();
      if (constrainedNodes.getLocation(tag) < 0) {
        int numDOF = theNode->getNumberDOF();
        if (numDOF != currentSize) 
          newVel.resize(numDOF);
        
        newVel = theNode->getVel();
        newVel(theDof) = vel0;
        theNode->setTrialVel(newVel);
        theNode->commitState();
      }
    }
  }
}


static inline double 
LeviCevita(int i, int j, int k)
{
  if (i == j || j == k || i == k)
    return 0;
  else if ((i == 0 && j == 1 && k == 2) ||
           (i == 1 && j == 2 && k == 0) ||
           (i == 2 && j == 0 && k == 1))
    return 1;
  else
    return -1;
}


static inline double
NodeAcceleration(Node& node, int NDM, int theDof, int dof, double uaccel)
{
  int irot = -1;
  if (NDM < 2)
    ;
  else if (NDM == 2 && theDof == 2)
    irot = 2;
  else if (NDM == 3 && theDof >= 3)
    irot = theDof - NDM;

  const Vector &xyz = node.getCrds();
  const int ndm = xyz.Size();
  const int numDOF = node.getNumberDOF();

  if (dof >= numDOF)
    return 0.0;
  
  else if (irot == -1) {
    // Translational excitation
    return (dof == theDof)? uaccel : 0.0;
  }
  else if (theDof == dof) {
    // Rotation is at dof
    return uaccel;
  }
  else if (dof < NDM) {
    // Rotational excitation at a translational dof
    double accel = 0.0;
    for (int j=0; j<ndm; j++)
      accel += LeviCevita(irot, j, dof)*xyz(j);

    return accel*uaccel;
  }
  return 0.0;
}


void
UniformExcitation::applyLoad(double time)
{

#if 1
  currentTime = time;
#else
  Domain *theDomain = this->getDomain();
  if (theDomain == nullptr)
      return;


  double accel = 0.0;
  if (theMotion != nullptr)
    accel = theMotion->getAccel(time)*fact;


  NodeIter &theNodes = theDomain->getNodes();
  Node *theNode;
  while ((theNode = theNodes()) != nullptr) {
      int ndf = theNode->getNumberDOF();
#if 1
      for (int i=0; i<ndf; i++) {
        double a = 0.0;
        if ((a = NodeAcceleration(*theNode, NDM, theDof, i, accel)) != 0)
          theNode->addResidualInertia(i, a);
      }
#else
      theNode->setNumColR(1);


      if (ndm == 1) {
        if (theDof < 1)
          theNode->setR(theDof, 0, fact);
      }
      else if (ndm == 2) {
        if (theDof < 2) {
          theNode->setR(theDof, 0, fact);
        }
        else if (theDof == 2) {
          double xCrd = crds(0);
          double yCrd = crds(1);
          theNode->setR(1, 0,  fact*xCrd);
          theNode->setR(0, 0, -fact*yCrd);
          theNode->setR(2, 0,  fact);
        }
      }
      else if (ndm == 3) {
        // Translational DOF
        if (theDof < 3) {
          theNode->setR(theDof, 0, fact);
        }

        // Rotational DOFs
        else if (theDof == 3) {
          double yCrd = crds(1);
          double zCrd = crds(2);
          theNode->setR(1, 0, -fact*zCrd);
          theNode->setR(2, 0,  fact*yCrd);
          theNode->setR(3, 0,  fact);
        }
        else if (theDof == 4) {
          double xCrd = crds(0);
          double zCrd = crds(2);
          theNode->setR(0, 0,  fact*zCrd);
          theNode->setR(2, 0, -fact*xCrd);
          theNode->setR(4, 0,  fact);
        }
        else if (theDof == 5) {
          double xCrd = crds(0);
          double yCrd = crds(1);
          theNode->setR(0, 0, -fact*yCrd);
          theNode->setR(1, 0,  fact*xCrd);
          theNode->setR(5, 0,  fact);
        }
      }
#endif
  }

#if 1
  static VectorND<6> an{}, fn{};
  static MatrixND<6,6> mn{}; // element mass at node


  ElementIter &theElements = theDomain->getElements();
  Element *theElement;
  while ((theElement = theElements()) != nullptr) {
    Node **theNodes = theElement->getNodePtrs();
    int numNodes = theElement->getNumExternalNodes();
    // assert(numNodes <= MaxNumNodesPerElement);

    switch (theElement->getMassType()) {
      case Element::MassType::Translation: {
        break;
      }
      case Element::MassType::General: {
        const Matrix &mass = theElement->getMass();

        for (int i=0; i<numNodes; i++) {
          // 1) Form nodal acceleration vector
          // const Vector &crds = theNodes[i]->getCrds();
          // int ndm = crds.Size();
          int ndfi = std::min(theNodes[i]->getNumberDOF(), 6);
          for (int j=0; j<ndfi; j++) {
            double a = 0.0;
            if ((a = NodeAcceleration(*theNodes[i], NDM, theDof, j, accel)) != 0)
              an[j] = a;
            else
              an[j] = 0.0;
          }

          // 2) Compute nodal forces
          for (int j=0; j<numNodes; j++) {
            int ndfj = std::min(theNodes[j]->getNumberDOF(), 6);
            fn.zero();
            mn.zero();
            for (int k=0; k<ndfj; k++) {
              for (int l=0; l<ndfi; l++) {
                mn(k,l) = mass(k + ndfj*j, l + ndfi*i);
              }
            }
            fn.addMatrixVector(0.0, mn, an, 1.0);
            theNodes[j]->addResidual(Vector(&fn[0],ndfj), -1.0);
          }
        }
      }
    }
  }
#else 
  // this->EarthquakePattern::applyLoad(time);
  {
    // see if quick return, i.e. no Ground Motions or domain set
    if (numMotions == 0)
      return;

    // check if setLoadConstant() has been called
    if (isConstant != true)
      currentTime = time;

    Domain *theDomain = this->getDomain();
    if (theDomain == nullptr)
      return;

  #if 0
    // set the vel and accel vector
    for (int i=0; i<numMotions; i++)
      (*uDotDotG)(i) = theMotions[i]->getAccel(currentTime);

    NodeIter &theNodes = theDomain->getNodes();
    Node *theNode;
    while ((theNode = theNodes()) != nullptr)
      theNode->addInertiaLoadToUnbalance(*uDotDotG, 1.0);


    ElementIter &theElements = theDomain->getElements();
    Element *theElement;
    while ((theElement = theElements()) != nullptr) 
      theElement->addInertiaLoadToUnbalance(*uDotDotG);
  #endif
  }
#endif

#endif
  return;
}

int
UniformExcitation::applyResidual(AnalysisModel &theAnalysisModel, LinearSOE &theSOE, double c)
{
#if 0
  return 0;
#else
  // Vector A(theSOE.getNumEqn()), B(theSOE.getNumEqn());
  static Vector A(236);
  A.resize(theAnalysisModel.getNumEqn());

  Domain *theDomain = this->getDomain();
  if (theDomain == nullptr)
    return -1;

  const double time = theDomain->getCurrentTime();

  double accel = 0.0;
  if (theMotion != nullptr)
    accel = theMotion->getAccel(time)*fact;


  NodeIter &theNodes = theDomain->getNodes();
  Node *theNode;
  while ((theNode = theNodes()) != nullptr) {
    const ID& idn = theNode->getDOF_GroupPtr()->getID();
    // int ndf = theNode->getNumberDOF();
    const int ndf = idn.Size();
    for (int i=0; i<ndf; i++) {
      double a = 0.0;
      if (idn(i) >=0 && (a = NodeAcceleration(*theNode, NDM, theDof, i, accel)) != 0)
        A(idn(i)) = a;
      else if (idn(i) >= 0)
        A(idn(i)) = 0.0;
    }
  }

  //
  //
  //
  theAnalysisModel.applyInertia(A, theSOE, -c);

  return 0;
#endif
}



void
UniformExcitation::applyLoadSensitivity(double time)
{
  Domain *theDomain = this->getDomain();
  if (theDomain == nullptr)
    return;

  // TODO(cmp)
  // NodeIter &theNodes = theDomain->getNodes();
  // Node *theNode;
  // while ((theNode = theNodes()) != nullptr) {
  //   theNode->setNumColR(1);
  //   theNode->setR(theDof, 0, 1.0);
  // }

  // this->EarthquakePattern::applyLoadSensitivity(time);
  {


    Domain *theDomain = this->getDomain();
    if (theDomain == nullptr)
      return;


    // set the vel and accel vector
    double uDotG = theMotion->getVel(time);
    double uDotDotG;
    if (parameterID != 0) { // Something is random in the motions
      uDotDotG = theMotion->getAccelSensitivity(time);
    }
    else {
      uDotDotG = theMotion->getAccel(time);
    }

    bool somethingRandomInMotions = false;
    if (parameterID != 0) {
      somethingRandomInMotions = true;
    }

#if 0
    NodeIter &theNodes = theDomain->getNodes();
    Node *theNode;
    while ((theNode = theNodes()) != 0) 
      theNode->addInertiaLoadSensitivityToUnbalance(*uDotDotG, 1.0,  somethingRandomInMotions);


    ElementIter &theElements = theDomain->getElements();
    Element *theElement;
    while ((theElement = theElements()) != 0) 
      theElement->addInertiaLoadSensitivityToUnbalance(*uDotDotG,  somethingRandomInMotions);
#endif
  }

  return;
}



int
UniformExcitation::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}


int 
UniformExcitation::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return -1;
}


void 
UniformExcitation::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_ELEM_INDENT << "{";
    s << "\"name\": \"" << this->getTag() << "\", ";
    s << "\"type\": \"" << this->getClassType() <<"\", ";
    s << "\"dof\": " << theDof+1 << ", ";
    s << "\"vel0\": " << vel0 << ", ";
    s << "\"fact\": " << fact;
    // if (theMotion != 0) {
    //   s << ", \"groundMotion\": ";
    //   theMotion->Print(s, flag);
    // }
    s << "}";

  }
  else {
    s << "UniformExcitation  " << this->getTag() << "\n";
  }
}
