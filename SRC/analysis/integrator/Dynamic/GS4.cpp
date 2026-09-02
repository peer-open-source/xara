//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, Claudio M. Perez
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
//===----------------------------------------------------------------------===//
//
// Description: This file contains the implementation of the GS4 class.
//
// Written : cmp
// Created : 06/2024
// Revision: A
//
#include <stdexcept>
#include <GS4.h>
#include <FE_Element.h>
#include <LinearSOE.h>
#include <AnalysisModel.h>
#include <Vector.h>
#include <DOF_Group.h>
#include <AnalysisModel.h>
// #include <string.h>
// for sensitivity
#include <Domain.h>
#include <NodeIter.h>
#include <Node.h>
#include <LoadPattern.h>
#include <Parameter.h>
#include <ParameterIter.h>


enum {
  eu=0,
  ev=1,
  ea=2
};

int 
GS4::setConstants(int unknown, 
                  double deltaT, double beta, double gamma,
                  std::array<double,3> const &alpha,
                  GammaScheme& scheme)
{
  switch (unknown) {
  case GS4::Displacement:
    if (beta == 0.0)  {
      opserr << "Invalid beta, requires beta != 0.0\n";
      return -1;
    }
    scheme.g[eu] = 1.0;
    scheme.g[ev] = gamma/(beta*deltaT);
    scheme.g[ea] = 1.0/(beta*deltaT*deltaT);


    scheme.G[eu][eu] =  1.0 - alpha[eu];
    scheme.G[eu][ev] =  0.0;
    scheme.G[eu][ea] =  0.0;

    scheme.G[ev][eu] = -gamma*alpha[ev]/(beta*deltaT);
    scheme.G[ev][ev] = 1.0 -  alpha[ev]*gamma/beta; 
    scheme.G[ev][ea] = deltaT*(1.0 - 0.5*gamma/beta);

    scheme.G[ea][eu] = -1.0/(beta*deltaT*deltaT);
    scheme.G[ea][ev] = -1.0/(beta*deltaT);
    scheme.G[ea][ea] =  1.0 - 0.5/beta;
    return 0;

  case GS4::Velocity:
    if (gamma == 0.0)  {
      opserr << "Invalid gamma, requires gamma != 0.0\n";
      return -1;
    }
    scheme.g[eu] = deltaT*beta/gamma;
    scheme.g[ev] = 1.0;
    scheme.g[ea] = 1.0/(gamma*deltaT);


    scheme.G[eu][eu] = 1.0;
    scheme.G[eu][ev] = -deltaT*beta/gamma*(1.0 - gamma/beta);
    scheme.G[eu][ea] =  deltaT*deltaT*beta/gamma*(gamma*0.5/beta - 1.0);

    scheme.G[ev][eu] = 0.0;
    scheme.G[ev][ev] = 0.0;
    scheme.G[ev][ea] = 0.0;

    scheme.G[ea][eu] = 0.0;
    scheme.G[ea][ev] = -1/(gamma*deltaT);
    scheme.G[ea][ea] =  1 - 1/gamma;
    return 0;

  case GS4::Acceleration:
    scheme.g[eu] = beta*deltaT*deltaT;
    scheme.g[ev] = gamma*deltaT;
    scheme.g[ea] = 1.0;


    scheme.G[eu][eu] = 1.0;
    scheme.G[eu][ev] = deltaT;
    scheme.G[eu][ea] = deltaT*deltaT*(0.5 - beta);

    scheme.G[ev][eu] = 0.0;
    scheme.G[ev][ev] = 1.0;
    scheme.G[ev][ea] = deltaT*(1.0 - gamma);

    scheme.G[ea][eu] = 0.0;
    scheme.G[ea][ev] = 0.0;
    scheme.G[ea][ea] = 0.0;
    return 0;

  default:
    opserr << "GS4::SetConstants -- unknown type " << unknown << endln;
    return -1;
  }
  return 0;
}

GS4::GS4(double gamma,  double beta, 
          double alphaF, double alphaM,
          int uFlag, int iFlag, bool aFlag)
  : TransientIntegrator(0),
    gamma(gamma), beta(beta), 
    alphaF(1.0), alphaM(1.0),
    alpha{1, 1, 1},
    unknown(uFlag), unknown_initialize(iFlag),
    step(0),
    dt(0.0),
    cu(0.0), cv(0.0), ca(0.0),
    G1(), G2(),
    Uo(nullptr), Vo(nullptr), Ao(nullptr),
    Ua(nullptr), Va(nullptr), Aa(nullptr),
    Un(nullptr), Vn(nullptr), An(nullptr),
    determiningMass(false),
    isSensitivityResidual(0), gradNumber(0), 
    dAa(0),
    dVa(0), 
    assemblyFlag(aFlag),
    dUn(), dVn(), dAn()
{

}


GS4::~GS4()
{
  if (Uo != nullptr)
    delete Uo;
  if (Vo != nullptr)
    delete Vo;
  if (Ao != nullptr)
    delete Ao;
  if (Ua != nullptr)
    delete Ua;
  if (Va != nullptr)
    delete Va;
  if (Aa != nullptr)
    delete Aa;
  if (Un != nullptr)
    delete Un;
  if (Vn != nullptr)
    delete Vn;
  if (An != nullptr)
    delete An;

  // clean up sensitivity
  if (dAa != nullptr)
    delete dAa;
  
  if (dVa != nullptr)
    delete dVa;
}


int
GS4::newStep(double deltaT)
{
  if (deltaT <= 0.0)  {
    opserr << "GS4::newStep() - error in variable\n";
    opserr << "dT = " << deltaT << endln;
    return -2;  
  }
  
  if (Un == nullptr)
    return -3;

  // mark step as bootstrap or not
  if (deltaT != dt)
    step = 0;
  else
    step++;

  dt = deltaT;

  //
  // Set response at t to be that at t+deltaT of previous step
  //
  (*Uo) = *Un;        
  (*Vo) = *Vn;  
  (*Ao) = *An;

  // set scheme constants
  GS4::setConstants(unknown, dt, beta, gamma, alpha, G1);
  cu = G1.g[eu];
  cv = G1.g[ev];
  ca = G1.g[ea];
  double buu = G1.G[eu][eu],
         buv = G1.G[eu][ev],
         bua = G1.G[eu][ea],

         bvu = G1.G[ev][eu],
         bvv = G1.G[ev][ev],
         bva = G1.G[ev][ea],

         bau = G1.G[ea][eu],
         bav = G1.G[ea][ev],
         baa = G1.G[ea][ea];

  //
  // Set initial guesses for {u,v,a}_{t + dt}
  //
  int init = unknown_initialize; // (step < 2 && unknown==Displacement) ? unknown : unknown_initialize;

  if (unknown == Displacement) {

    switch (init) {
      case Displacement:
    //  Vn->addVector(bvv, *Ao, bva);
        Vn->addVector(0.0, *Uo,  bvu  + cv/cu*(1.0 - buu)); //  = 0
        Vn->addVector(1.0, *Vo,  bvv  + cv/cu*(    - buv)); // += bvv*Vo
        Vn->addVector(1.0, *Ao,  bva  + cv/cu*(    - bua)); // += bva*Ao

        An->addVector(baa, *Vo, bav);
        break;

      case Acceleration:
        Un->addVector(0.0, *Uo,  buu  + cu/ca*(    - bau));
        Un->addVector(1.0, *Vo,  buv  + cu/ca*(    - bav));
        Un->addVector(1.0, *Ao,  bua  + cu/ca*(1.0 - baa));

        Vn->addVector(0.0, *Uo,  bvu  + cv/ca*(    - bau));
        Vn->addVector(1.0, *Vo,  bvv  + cv/ca*(    - bav));
        Vn->addVector(1.0, *Ao,  bva  + cv/ca*(1.0 - baa));
        break;

      case Velocity:
        Un->addVector(0.0, *Uo,  buu  + cu*(    - bvu)/cv);
        Un->addVector(1.0, *Vo,  buv  + cu*(1.0 - bvv)/cv);
        Un->addVector(1.0, *Ao,  bua  + cu*(    - bva)/cv);

        An->addVector(0.0, *Uo,  bau  + ca*(    - bvu)/cv);
        An->addVector(1.0, *Vo,  bav  + ca*(1.0 - bvv)/cv);
        An->addVector(1.0, *Ao,  baa  + ca*(    - bva)/cv);
        break;
    }
  }

  else if (unknown == Velocity) {
    switch (init) {
      case Displacement:
        Vn->addVector(0.0, *Uo,  bvu + cv*(1.0 - buu)/cu);
        Vn->addVector(1.0, *Vo,  bvv + cv*(    - buv)/cu);
        Vn->addVector(1.0, *Ao,  bva + cv*(    - bua)/cu);

        // a += c3*a_{n+1}
        An->addVector(0.0, *Uo,  bau + ca*(1.0 - buu)/cu);
        An->addVector(1.0, *Vo,  bav + ca*(    - buv)/cu);
        An->addVector(1.0, *Ao,  baa + ca*(    - bua)/cu);
        break;

      case Velocity:
        // TODO: Check
        Un->addVector(0.0, *Uo,  buu + cu*(    - bvu)/cv);
        Un->addVector(1.0, *Vo,  buv + cu*(1.0 - bvv)/cv);
        Un->addVector(1.0, *Ao,  bua + cu*(    - bva)/cv);

        An->addVector(0.0, *Uo,  bau + ca*(    - bvu)/cv);
        An->addVector(1.0, *Vo,  bav + ca*(1.0 - bvv)/cv);
        An->addVector(1.0, *Ao,  baa + ca*(    - bva)/cv);
        break;

      case Acceleration:
        Un->addVector(0.0, *Uo,  buu + cu*(    - bau)/ca);
        Un->addVector(1.0, *Vo,  buv + cu*(    - bav)/ca);
        Un->addVector(1.0, *Ao,  bua + cu*(1.0 - baa)/ca);

        Vn->addVector(0.0, *Uo,  bvu + cv*(    - bau)/ca);
        Vn->addVector(1.0, *Vo,  bvv + cv*(    - bav)/ca);
        Vn->addVector(1.0, *Ao,  bva + cv*(1.0 - baa)/ca);
        break;
    }
  }
  else {
    switch (init) {
      case Displacement:
        // Initialize: U == Ut
        // implying   Da = -vc/(beta dt) - ac/(2 beta)

        Vn->addVector(0.0, *Uo,  bvu + cv*(1.0 - buu)/cu); // 0
        Vn->addVector(1.0, *Vo,  bvv + cv*(    - buv)/cu); // (beta*deltaT));
        Vn->addVector(1.0, *Ao,  bva + cv*(    - bua)/cu);

        An->addVector(0.0, *Uo,  bau + ca*(1.0 - buu)/cu); // 0
        An->addVector(1.0, *Vo,  bav + ca*(    - buv)/cu);
        An->addVector(1.0, *Ao,  baa + ca*(    - bua)/cu);
        break; 

      case Velocity: // TODO: Check
        // Dv = 0
        Un->addVector(0.0, *Uo,  buu + cu*(    - bvu)/cv);
        Un->addVector(1.0, *Vo,  buv + cu*(1.0 - bvv)/cv);
        Un->addVector(1.0, *Ao,  bua + cu*(    - bva)/cv);

        // a += c3*a_{n+1}
        An->addVector(0.0, *Uo,  bau + ca*(    - bvu)/cv);
        An->addVector(1.0, *Vo,  bav + ca*(1.0 - bvv)/cv);
        An->addVector(1.0, *Ao,  baa + ca*(    - bva)/cv);
        break;

      case Acceleration:
        // Da = 0
        Un->addVector(buu, *Vo,  buv);
        Un->addVector(1.0, *Ao,  bua + cu);

        Vn->addVector(0.0, *Uo,  bvu);
        Vn->addVector(1.0, *Vo,  bvv);
        Vn->addVector(1.0, *Ao,  bva + cv); // deltaT
        break;
    }
  }

  //
  // set the trial response quantities
  //
  // determine the displacements at t+alphaF*deltaT
  (*Ua) = *Uo;
  Ua->addVector((1.0-alphaF), *Un, alphaF);

  // determine the velocities at t+alphaF*deltaT
  (*Va) = *Vo;
  Va->addVector((1.0-alphaF), *Vn, alphaF);

  // determine the velocities at t+alphaM*deltaT
  (*Aa) = *Ao;
  Aa->addVector((1.0-alphaM), *An, alphaM);

  AnalysisModel *theModel = this->getAnalysisModel();
  
  theModel->setResponse(*Ua, *Va, *Aa);

  //
  // increment the time and apply the load
  //
  double time = theModel->getCurrentDomainTime();
  time += alphaF*deltaT;
  if (theModel->updateDomain(time, deltaT) < 0)  {
    opserr << "GS4::newStep - failed to update\n";
    return -4;
  }

  return 0;
}


int
GS4::update(const Vector &deltaX)
{
  AnalysisModel *theModel = this->getAnalysisModel();
  assert(theModel != nullptr);
  
  // Check domainChanged() has been called, i.e. Ut is not null
  assert(Uo != nullptr);

  // check deltaX is of correct size
  if (deltaX.Size() != Un->Size())  {
    opserr << "WARNING GS4::update() - Vectors of incompatible size ";
    opserr << " expecting " << Un->Size() << " obtained " << deltaX.Size() << endln;
    return -3;
  }
  
  //  determine the response at t+deltaT
  Un->addVector(1.0, deltaX, G1.g[eu]);
  Vn->addVector(1.0, deltaX, G1.g[ev]);
  An->addVector(1.0, deltaX, G1.g[ea]);

  // determine state at t + alpha*deltaT
  Ua->addVector(1.0, deltaX, alphaF*G1.g[eu]);
  Va->addVector(1.0, deltaX, alphaF*G1.g[ev]);
  Aa->addVector(1.0, deltaX, alphaM*G1.g[ea]);

//  (*Ua) = *Uo;
//  Ua->addVector((1.0-alphaF), *Un, alphaF);

//  (*Va) = *Vo;
//  Va->addVector((1.0-alphaF), *Vn, alphaF);

//  (*Aa) = *Ao;
//  Aa->addVector((1.0-alphaM), *An, alphaM);

  // update the response at the DOFs
  theModel->setResponse(*Ua,*Va,*Aa);
  if (theModel->updateDomain() < 0)  {
    opserr << "GS4::update - failed to update the domain\n";
    return -4;
  }
  
  return 0;
}    


const Vector &
GS4::getVel()
{
  return *Va;
}

int
GS4::revertToLastStep()
{
  // set response at t+deltaT to be that at t .. for next newStep
  if (Un != nullptr)  {
    (*Un) = *Uo;        
    (*Vn) = *Vo;  
    (*An) = *Ao;  
  }
  return 0;
}


int
GS4::formEleTangent(FE_Element *theEle)
{
  if (determiningMass == true)
    return 0;

  theEle->zeroTangent();
  
  switch (statusFlag) {
  case CURRENT_TANGENT:
    theEle->addKtToTang(alphaF*G1.g[eu]);
    theEle->addCtoTang( alphaF*G1.g[ev]);
    theEle->addMtoTang( alphaM*G1.g[ea]);
    break;
  case INITIAL_TANGENT:
    theEle->addKiToTang(alphaF*G1.g[eu]);
    theEle->addCtoTang( alphaF*G1.g[ev]);
    theEle->addMtoTang( alphaM*G1.g[ea]);
    break;
  case HALL_TANGENT:
    theEle->addKtToTang(cu*cFactor);
    theEle->addKiToTang(cu*iFactor);
    theEle->addCtoTang(cv);
    theEle->addMtoTang(ca);
    break;
  }
  return 0;
}

int
GS4::formNodTangent(DOF_Group *theDof)
{
  if (determiningMass == true)
    return 0;

  theDof->zeroTangent();
  theDof->addCtoTang(alphaF*G1.g[ev]);
  theDof->addMtoTang(alphaM*G1.g[ea]);

  return 0;
}


int
GS4::domainChanged()
{
  AnalysisModel *myModel = this->getAnalysisModel();
  LinearSOE *theLinSOE = this->getLinearSOE();
  const Vector &x = theLinSOE->getX();
  int size = x.Size();

  // create the new Vector objects
  if (Uo == nullptr || Uo->Size() != size)  {
    // delete the old
    if (Uo != nullptr)
      delete Uo;
    if (Vo != nullptr)
      delete Vo;
    if (Ao != nullptr)
      delete Ao;
    if (Ua != nullptr)
      delete Ua;
    if (Va != nullptr)
      delete Va;
    if (Aa != nullptr)
      delete Ao;
    if (Un != nullptr)
      delete Un;
    if (Vn != nullptr)
      delete Vn;
    if (An != nullptr)
      delete An;

    // perform the allocations
    Uo = new Vector(size);
    Vo = new Vector(size);
    Ao = new Vector(size);
    Ua = new Vector(size);
    Va = new Vector(size);
    Aa = new Vector(size);
    Un = new Vector(size);
    Vn = new Vector(size);
    An = new Vector(size);
    dUn.resize(size); 
    dUn.Zero();
    dVn.resize(size); 
    dVn.Zero();
    dAn.resize(size); 
    dAn.Zero(); 
  }        
  
  // now go through and populate U, Udot and Udotdot by iterating through
  // the DOF_Groups and getting the last committed velocity and accel
  DOF_GrpIter &theDOFs = myModel->getDOFs();
  DOF_Group *group;
  while ((group = theDOFs()) != nullptr)  {
    const ID &id = group->getID();
    int idSize = id.Size();

    const Vector &disp = group->getCommittedDisp();  
    for (int i=0; i < idSize; i++)  {
      int loc = id(i);
      if (loc >= 0)  {
        (*Un)(loc) = disp(i);    
      }
    }
    
    const Vector &vel = group->getCommittedVel();
    for (int i=0; i < idSize; i++)  {
      int loc = id(i);
      if (loc >= 0) {
        (*Vn)(loc) = vel(i);
      }
    }

    const Vector &accel = group->getCommittedAccel();  
    for (int i=0; i < idSize; i++)  {
      int loc = id(i);
      if (loc >= 0) {
        (*An)(loc) = accel(i);
      }
    }

    // The remaining get**Sensitivity methods cause seg faults with Lagrange constraint
    // handler in dynamic (transient) analysis even when there is no sensitivity algorithm.
    // However, I don't think these methods need to be called in domainChanged -- MHS
    continue;
    
    const Vector &dispSens = group->getDispSensitivity(gradNumber);  
    for (int i=0; i < idSize; i++) {
        int loc = id(i);
        if (loc >= 0) {
          dUn(loc) = dispSens(i);    
        }
    }

    const Vector &velSens = group->getVelSensitivity(gradNumber);
    for (int i=0; i < idSize; i++) {
        int loc = id(i);
        if (loc >= 0) {
          dVn(loc) = velSens(i);
        }
    }

    const Vector &accelSens = group->getAccSensitivity(gradNumber);  
    for (int i=0; i < idSize; i++) {
        int loc = id(i);
        if (loc >= 0) {
          dAn(loc) = accelSens(i);
        }
    }
  }    
  
  return 0;
}


int
GS4::revertToStart()
{
  if (Uo != nullptr) 
    Uo->Zero();
  if (Vo != nullptr) 
    Vo->Zero();
  if (Ao != nullptr) 
    Ao->Zero();
  if (Un != nullptr) 
    Un->Zero();
  if (Vn != nullptr) 
    Vn->Zero();
  if (An != nullptr) 
    An->Zero();

  return 0;
}

int
GS4::formEleResidual(FE_Element* elem)
{
  if (isSensitivityResidual == false) {
    this->TransientIntegrator::formEleResidual(elem);
  }
  else {

    elem->zeroResidual();


    // So, the constants can be computed as follows:
    if (unknown != Displacement) {
      opserr << "ERROR: GS4::formEleResidual() -- the implemented"
        << " scheme only works if the displ variable is set to true." << endln;
    }

    double dt = gamma/(beta*cv);

    double a2 = -ca;
    double a3 = -cv/gamma;
    double a4 = 1.0 - 1.0/(2.0*beta);
    double bvu = -cv;
    double bvv = 1.0 - gamma/beta;
    double bva = dt*(1.0 - gamma/(2.0*beta));

    // Pre-compute the vectors involving a2, a3, etc.
    //Vector tmp1 = V*a2 + Vdot*a3 + Vdotdot*a4;
    const int n = Un->Size();
    Vector dUn(n);
    Vector dVn(n);
    Vector dAn(n);

    AnalysisModel *myModel = this->getAnalysisModel();
    myModel->getStateGradient(dUn, dVn, dAn, gradNumber);

    // Pre-compute the vectors involving a2, a3, etc.
    // Vector tmp1 = V*a2 + Vdot*a3 + Vdotdot*a4;
    Vector tmp1(n);
    tmp1.addVector(0.0, dUn, a2);
    tmp1.addVector(1.0, dVn, a3);
    tmp1.addVector(1.0, dAn, a4);
    //Vector tmp2 = V*bvu + Vdot*bvv + Vdotdot*bva;
    Vector tmp2(n);
    tmp2.addVector(0.0, dUn, bvu);
    tmp2.addVector(1.0, dVn, bvv);
    tmp2.addVector(1.0, dAn, bva);

    if (dAa == nullptr)
      dAa = new Vector(tmp1.Size());

    if (dVa == nullptr)
      dVa = new Vector(tmp2.Size());

    (*dAa) = tmp1;
    (*dVa) = tmp2;


    // Now we're ready to make calls to the FE Element:

    // The term -dPint/dh|u fixed
    elem->addResistingForceSensitivity(gradNumber); 

    // The term -dM/dh*acc
    elem->addM_ForceSensitivity(gradNumber, *An, -1.0);

    // The term -M*(a2*v + a3*vdot + a4*vdotdot)
    elem->addM_Force(*dAa, -1.0);

    // The term -C*(bvu*v + bvv*vdot + bva*vdotdot)
    elem->addD_Force(*dVa, -1.0);

    // The term -dC/dh*vel
    elem->addD_ForceSensitivity(gradNumber, *Vn, -1.0);
  }

  return 0;
}

int
GS4::formNodUnbalance(DOF_Group *theDof)
{

  if (isSensitivityResidual == 0) {
    this->TransientIntegrator::formNodUnbalance(theDof);
  }
  else {
    // Assemble sensitivity residual

    theDof->zeroUnbalance();

    // The term -M*(a2*v + a3*vdot + a4*vdotdot)
    theDof->addM_Force(*dAa,-1.0);

    // The term -dM/dh*acc
    theDof->addM_ForceSensitivity(*An, -1.0);

    // The term -C*(bvu*v + bvv*vdot + bva*vdotdot)
    theDof->addD_Force(*dVa,-1.0);

    // The term -dC/dh*vel
    theDof->addD_ForceSensitivity(*Vn,-1.0);

    // In case of random loads (have already been formed by 'applyLoadSensitivity')
    theDof->addPtoUnbalance();

  }
  return 0;
}

int 
GS4::formSensitivityRHS(int grad)
{
  // Set to sensitivity mode. This will change the behaviour of formEleResidual
  // and formNodUnbalance so that they add the sensitivity terms to the RHS
  // rather than the standard terms.
  isSensitivityResidual = true;
  gradNumber = grad;

  LinearSOE *theSOE = this->getLinearSOE();
  AnalysisModel *theModel = this->getAnalysisModel();

  // 3)

  // element/material contributions
  FE_Element *elem;
  FE_EleIter &theEles = theModel->getFEs();    
  while ((elem = theEles()) != nullptr)
    theSOE->addB(elem->getResidual(this), elem->getID());


  // Loop through DOF groups (IT IS IMPORTANT THAT THIS IS DONE LAST!)
  DOF_Group *group;
  DOF_GrpIter &theDOFs = theModel->getDOFs();
  while ((group = theDOFs()) != nullptr)
    theSOE->addB(group->getUnbalance(this),  group->getID());


  // Reset the sensitivity flag
  isSensitivityResidual = false;

  return 0;
}


int 
GS4::updateGradient(const Vector & vNew,int gradNum,int numGrads)
{
  // Compute GS4 parameters
  double dt = gamma/(beta*cv);
  double a2 = -ca;
  double a3 = -cv/gamma;
  double a4 = 1.0 - 1.0/(2.0*beta);
  double a5 = cv;
  double bvu = -cv;
  double bvv = 1.0 - gamma/beta;
  double bva = dt*(1.0 - gamma/(2.0*beta));


  // Obtain sensitivity vectors from previous step modified by lei July 2018
  int vectorSize = Un->Size();
  Vector dUn(vectorSize);
  Vector dVn(vectorSize);
  Vector dAn(vectorSize);

  AnalysisModel *myModel = this->getAnalysisModel();
  myModel->getStateGradient(dUn, dVn, dAn, gradNum);

  // Compute new acceleration and velocity vectors:
  Vector vdotNew(vectorSize);
  Vector vdotdotNew(vectorSize);
  //(*vdotdotNewPtr) = vNew*c3 + V*a2 + Vdot*a3 + Vdotdot*a4;
  vdotdotNew.addVector(0.0, vNew, ca);
  vdotdotNew.addVector(1.0, dUn,  a2);
  vdotdotNew.addVector(1.0, dVn,  a3);
  vdotdotNew.addVector(1.0, dAn,  a4);
  
  //(*vdotNewPtr) = vNew*a5 + V*bvu + Vdot*bvv + Vdotdot*bva;
  vdotNew.addVector(0.0, vNew, a5);
  vdotNew.addVector(1.0, dUn, bvu);
  vdotNew.addVector(1.0, dVn, bvv);
  vdotNew.addVector(1.0, dAn, bva);

  // update
  dUn = vNew;
  dVn = vdotNew;
  dAn = vdotdotNew;

  myModel->setStateGradient(vNew, vdotNew, vdotdotNew, gradNum, numGrads);

  return 0;
}



double
GS4::getCFactor() 
{
  return alphaF*cv;
}


int 
GS4::computeSensitivities()
{
  LinearSOE *theSOE = this->getLinearSOE();

  // Zero out the old right-hand side of the SOE
  theSOE->zeroB();

  AnalysisModel *theModel = this->getAnalysisModel();  //Abbas 
  Domain *theDomain=theModel->getDomainPtr();//Abbas

  // De-activate all parameters
  Parameter *theParam;
  ParameterIter &paramIter = theDomain->getParameters();
  while ((theParam = paramIter()) != nullptr)
    theParam->activate(false);
  
  // Compute sensitivity wrt each parameter
  const int numGrads = theDomain->getNumParameters();
  paramIter = theDomain->getParameters();
  while ((theParam = paramIter()) != nullptr) {
    
    // Activate this parameter
    theParam->activate(true);
    
    // Zero the RHS vector
    theSOE->zeroB();
    
    // Get the grad index for this parameter
    int gradIndex = theParam->getGradIndex();

    // Form the residual
    theModel->applyLoadGradient();

    this->formSensitivityRHS(gradIndex);
    
    // Solve for displacement sensitivity 
    theSOE->solve();

    // Save sensitivity to nodes
    this->updateGradient( theSOE->getX(), gradIndex, numGrads );
    

    // Commit unconditional history variables (also for elastic problems; strain sens may be needed anyway)
    this->commitGradient(gradIndex, numGrads);
    theModel->commitGradient(gradIndex, numGrads);
    
    // De-activate this parameter for next parameter sensitivity
    theParam->activate(false);
  }
  
  return 0;
}



void
GS4::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON)
    return;

  AnalysisModel *theModel = this->getAnalysisModel();
  if (theModel != nullptr) {
      double currentTime = theModel->getCurrentDomainTime();
      s << "\t GS4 - currentTime: " << currentTime;
  }

  s << "\t gamma: " << gamma << "  beta: " << beta << "\n";
  s << "\t alphaF: " << alphaF << "  alphaM: " << alphaM << "\n";
  s << "\t unknown: " << unknown << "  initialization: " << unknown_initialize << "\n";
}

