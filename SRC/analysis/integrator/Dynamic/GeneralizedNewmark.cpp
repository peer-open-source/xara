//===----------------------------------------------------------------------===//
//
//                                   xara
//
//===----------------------------------------------------------------------===//
//                              https://xara.so
//===----------------------------------------------------------------------===//
//
// Description: This file contains the implementation of the GeneralizedNewmark class.
//
// Written : cmp
// Created : 06/2024
// Revision: A
//
#include <stdexcept>
#include <GeneralizedNewmark.h>
#include <FE_Element.h>
#include <FE_EleIter.h>
#include <LinearSOE.h>
#include <AnalysisModel.h>
#include <Vector.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <AnalysisModel.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <string.h>
#include <NodeIter.h>
#include <Domain.h>
// for sensitivity
#include <Node.h>
#include <LoadPattern.h>
#include <LoadPatternIter.h>
#include <Parameter.h>
#include <ParameterIter.h>


GeneralizedNewmark::GeneralizedNewmark(double gamma,  double beta, 
                                       double alphaF, double alphaM,
                                       int uFlag, int iFlag, bool aFlag)
  : TransientIntegrator(0),
    gamma(gamma), beta(beta), 
    alphaF(1.0), alphaM(1.0), 
    unknown(uFlag), unknown_initialize(iFlag),
    step(0),
    dt(0.0),
    cu(0.0), cv(0.0), ca(0.0), 
    Uo(nullptr), Vo(nullptr), Ao(nullptr),
    Ua(nullptr), Va(nullptr), Aa(nullptr),
    Un(nullptr), Vn(nullptr), An(nullptr),
    determiningMass(false),
    sensitivityFlag(0), gradNumber(0), 
    dAa(0),
    dVa(0), 
    assemblyFlag(aFlag), 
    independentRHS(),
    dUn(), dVn(), dAn()
{

}


GeneralizedNewmark::~GeneralizedNewmark()
{
  // clean up the memory created
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
GeneralizedNewmark::newStep(double deltaT)
{
  if (deltaT <= 0.0)  {
      opserr << "GeneralizedNewmark::newStep() - error in variable\n";
      opserr << "dT = " << deltaT << endln;
      return -2;  
  }
  
  if (Un == nullptr)
    throw std::invalid_argument("domainChange failed or not called");
    // return -3;

  // mark step as bootstrap or not
  if (deltaT != dt)
      step = 0;
  else
      step++;

  dt = deltaT;

  // Set response at t to be that at t+deltaT of previous step
  (*Uo) = *Un;        
  (*Vo) = *Vn;  
  (*Ao) = *An;

  // set the constants
  switch (unknown) {
  case Displacement:
    if (beta == 0.0)  {
        opserr << "Invalid beta, requires beta != 0.0\n";
        return -1;
    }
    cu = 1.0;
    cv = gamma/(beta*deltaT);
    ca = 1.0/(beta*deltaT*deltaT);
    break;

  case Velocity:
    if (gamma == 0.0)  {
        opserr << "Invalid gamma, requires gamma != 0.0\n";
        return -1;
    }
    cu = deltaT*beta/gamma;
    cv = 1.0;
    ca = 1.0/(gamma*deltaT);
    break;

  case Acceleration:
    cu = beta*deltaT*deltaT;
    cv = gamma*deltaT;
    ca = 1.0;
    break;
  }

  //
  // Set initial guesses for {u,v,a}_{t + dt}
  //
  int init = (step < 2 && unknown==Displacement) ? unknown : unknown_initialize;

  if (unknown == Displacement) {
    double buu =  0.0;
    double buv =  0.0;
    double bua =  0.0;

    double bvu = -gamma/(beta*deltaT);
    double bvv = 1.0 - gamma/beta; 
    double bva = deltaT*(1.0 - 0.5*gamma/beta);

    double bau = -1.0/(beta*deltaT*deltaT);
    double bav = -1.0/(beta*deltaT);
    double baa =  1.0 - 0.5/beta;

    switch (init) {
      case Displacement:
    //  Vn->addVector(bvv, *Ao, bva);
        Vn->addVector(0.0, *Uo,  bvu  + cv/cu*(1.0 - buu)); //  = 0
        Vn->addVector(1.0, *Vo,  bvv  + cv/cu*(    - buv)); // += bvv*Vo
        Vn->addVector(1.0, *Ao,  bva  + cv/cu*(    - bua)); // += bva*Ao

        // An->addVector(baa-1/beta, *Vo, bav);
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
    double buu = 1.0;
    double buv = -deltaT*beta/gamma*(1.0 - gamma/beta);
    double bua =  deltaT*deltaT*beta/gamma*(gamma*0.5/beta - 1.0);

    double bvu = 0.0;
    double bvv = 0.0;
    double bva = 0.0;

    double bau = 0.0;
    double bav = -1/(gamma*deltaT);
    double baa =  1 - 1/gamma;

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
    double buu = 1.0;
    double buv = deltaT;
    double bua = deltaT*deltaT*(0.5 - beta);

    double bvu = 0.0;
    double bvv = 1.0;
    double bva = deltaT*(1.0 - gamma);

    double bau = 0.0;
    double bav = 0.0;
    double baa = 0.0;

    // Choose how to initialize state
    // u += c1*Da
//    U->addVector(*Utdot,                   -c1*(      buv/c1));
//    U->addVector(*Utdotdot,                -c1*(1.0 + bua/c1));

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
    opserr << "GeneralizedNewmark::newStep() - failed to update the domain\n";
    return -4;
  }

  return 0;
}


int
GeneralizedNewmark::update(const Vector &deltaX)
{
  AnalysisModel *theModel = this->getAnalysisModel();
  if (theModel == nullptr)  {
      opserr << "WARNING GeneralizedNewmark::update() - no AnalysisModel set\n";
      return -1;
  }
  
  // check domainChanged() has been called, i.e. Ut will not be null
  if (Uo == nullptr)  {
      opserr << "WARNING GeneralizedNewmark::update() - domainChange() failed or not called\n";
      return -2;
  }  

  // check deltaX is of correct size
  if (deltaX.Size() != Un->Size())  {
      opserr << "WARNING GeneralizedNewmark::update() - Vectors of incompatible size ";
      opserr << " expecting " << Un->Size() << " obtained " << deltaX.Size() << endln;
      return -3;
  }
  
  //  determine the response at t+deltaT
  switch (unknown) {
  case Displacement:
    (*Un) += deltaX;
    Vn->addVector(1.0, deltaX, cv);
    An->addVector(1.0, deltaX, ca);
    break;

  case Velocity:
    Un->addVector(1.0, deltaX, cu);
    (*Vn) += deltaX;
    An->addVector(1.0, deltaX, ca);
    break;

  case Acceleration:
    Un->addVector(1.0, deltaX, cu);
    Vn->addVector(1.0, deltaX, cv);        
    (*An) += deltaX;
    break;
  }

  // determine displacement and velocity at t + alphaF*deltaT
  (*Ua) = *Uo;
  Ua->addVector((1.0-alphaF), *Un, alphaF);

  (*Va) = *Vo;
  Va->addVector((1.0-alphaF), *Vn, alphaF);

  // determine the velocities at t+alphaM*deltaT
  (*Aa) = *Ao;
  Aa->addVector((1.0-alphaM), *An, alphaM);

  // update the response at the DOFs
  theModel->setResponse(*Ua,*Va,*Aa);
  if (theModel->updateDomain() < 0)  {
    return -4;
  }
  
  return 0;
}    


const Vector &
GeneralizedNewmark::getVel()
{
  return *Vn;
}

int
GeneralizedNewmark::revertToLastStep()
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
GeneralizedNewmark::formEleTangent(FE_Element *theEle)
{
    if (determiningMass == true)
        return 0;

    theEle->zeroTangent();
    
    switch (statusFlag) {
    case CURRENT_TANGENT:
        theEle->addKtToTang(alphaF*cu);
        theEle->addCtoTang(alphaF*cv);
        theEle->addMtoTang(alphaM*ca);
        break;
    case INITIAL_TANGENT:
        theEle->addKiToTang(alphaF*cu);
        theEle->addCtoTang(alphaF*cv);
        theEle->addMtoTang(alphaM*ca);
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
GeneralizedNewmark::formNodTangent(DOF_Group *theDof)
{
  if (determiningMass == true)
      return 0;
  
  theDof->zeroTangent();
  theDof->addCtoTang(alphaF*cv);
  theDof->addMtoTang(alphaM*ca);
  
  return 0;
}    


int
GeneralizedNewmark::domainChanged()
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
    DOF_Group *dofPtr;
    while ((dofPtr = theDOFs()) != nullptr)  {
      const ID &id = dofPtr->getID();
      int idSize = id.Size();

      const Vector &disp = dofPtr->getCommittedDisp();  
      for (int i=0; i < idSize; i++)  {
          int loc = id(i);
          if (loc >= 0)  {
              (*Un)(loc) = disp(i);    
          }
      }
      
      const Vector &vel = dofPtr->getCommittedVel();
      for (int i=0; i < idSize; i++)  {
          int loc = id(i);
          if (loc >= 0) {
              (*Vn)(loc) = vel(i);
          }
      }
      
      const Vector &accel = dofPtr->getCommittedAccel();  
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
      
      const Vector &dispSens = dofPtr->getDispSensitivity(gradNumber);  
      for (int i=0; i < idSize; i++) {
          int loc = id(i);
          if (loc >= 0) {
            dUn(loc) = dispSens(i);    
          }
      }

      const Vector &velSens = dofPtr->getVelSensitivity(gradNumber);
      for (int i=0; i < idSize; i++) {
          int loc = id(i);
          if (loc >= 0) {
            dVn(loc) = velSens(i);
          }
      }

      const Vector &accelSens = dofPtr->getAccSensitivity(gradNumber);  
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
GeneralizedNewmark::sendSelf(int cTag, Channel &theChannel)
{
    Vector data(3);
    data(0) = gamma;
    data(1) = beta;
    data(2) = unknown;

    
    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING GeneralizedNewmark::sendSelf() - could not send data\n";
        return -1;
    }

    return 0;
}


int
GeneralizedNewmark::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  Vector data(3);
  if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0)  {
    opserr << "WARNING GeneralizedNewmark::recvSelf() - could not receive data\n";
    gamma = 0.5;
    beta = 0.25; 
    return -1;
  }
  
  gamma  = data(0);
  beta   = data(1);
  unknown  = data(2);

  return 0;
}


void
GeneralizedNewmark::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON)
    return;

  AnalysisModel *theModel = this->getAnalysisModel();
  if (theModel != nullptr) {
      double currentTime = theModel->getCurrentDomainTime();
      s << "\t GeneralizedNewmark - currentTime: " << currentTime;
  }
  s << "\t gamma: " << gamma << "  beta: " << beta << "\n";
  s << "\t alphaF: " << alphaF << "  alphaM: " << alphaM << "\n";
  s << "\t unknown: " << unknown << "  initialization: " << unknown_initialize << "\n";
}



int
GeneralizedNewmark::revertToStart()
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
GeneralizedNewmark::formEleResidual(FE_Element* theEle)
{
  if (sensitivityFlag == 0) {  // no sensitivity
      this->TransientIntegrator::formEleResidual(theEle);

  } else {
  
      theEle->zeroResidual();

      // Compute the time-stepping parameters of the form
      // udotdot = c3*ui+1 + a2*ui + a3*udoti + a4*udotdoti
      // udot    = a5*ui+1 + bvu*ui + bvv*udoti + bva*udotdoti
      // (see p. 166 of Chopra)

      // The constants are:
      // c3 = 1.0/(beta*dt*dt)
      // a2 = -1.0/(beta*dt*dt)
      // a3 = -1.0/beta*dt
      // a4 = 1.0 - 1.0/(2.0*beta)
      // a5 = gamma/(beta*dt)
      // bvu = -gamma/(beta*dt)
      // bvv = 1.0 - gamma/beta
      // bva = 1.0 - gamma/(2.0*beta)

      // So, the constants can be computed as follows:
      if (unknown != Displacement) {
          opserr << "ERROR: GeneralizedNewmark::formEleResidual() -- the implemented"
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
      int vectorSize = Un->Size();
      Vector dUn(vectorSize);
      Vector dVn(vectorSize);
      Vector dAn(vectorSize);

      AnalysisModel *myModel = this->getAnalysisModel();
      DOF_GrpIter &theDOFs = myModel->getDOFs();
      DOF_Group *dofPtr;
      while ((dofPtr = theDOFs()) != nullptr) {

        const ID &id = dofPtr->getID();
        const int idSize = id.Size();
        const Vector &dispSens = dofPtr->getDispSensitivity(gradNumber);
        for (int i = 0; i < idSize; i++)
          if (int loc = id(i); loc >= 0)
            dUn(loc) = dispSens(i);


        const Vector &velSens = dofPtr->getVelSensitivity(gradNumber);
        for (int i = 0; i < idSize; i++)
          if (int loc = id(i); loc >= 0)
            dVn(loc) = velSens(i);

        const Vector &accelSens = dofPtr->getAccSensitivity(gradNumber);
        for (int i = 0; i < idSize; i++)
          if (int loc = id(i); loc >= 0)
            dAn(loc) = accelSens(i);
      }

      // Pre-compute the vectors involving a2, a3, etc.
      // Vector tmp1 = V*a2 + Vdot*a3 + Vdotdot*a4;
      Vector tmp1(vectorSize);
      tmp1.addVector(0.0, dUn, a2);
      tmp1.addVector(1.0, dVn, a3);
      tmp1.addVector(1.0, dAn, a4);
      //Vector tmp2 = V*bvu + Vdot*bvv + Vdotdot*bva;
      Vector tmp2(vectorSize);
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
      theEle->addResistingForceSensitivity(gradNumber); 

      // The term -dM/dh*acc
      theEle->addM_ForceSensitivity(gradNumber, *An, -1.0);

      // The term -M*(a2*v + a3*vdot + a4*vdotdot)
      theEle->addM_Force(*dAa, -1.0);

      // The term -C*(bvu*v + bvv*vdot + bva*vdotdot)
      theEle->addD_Force(*dVa, -1.0);

      // The term -dC/dh*vel
      theEle->addD_ForceSensitivity(gradNumber, *Vn, -1.0);
    }

    return 0;
}

int
GeneralizedNewmark::formNodUnbalance(DOF_Group *theDof)
{

    if (sensitivityFlag == 0) {
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
GeneralizedNewmark::formSensitivityRHS(int passedGradNumber)
{
  // Set a couple of data members
  sensitivityFlag = 1;
  gradNumber = passedGradNumber;


  LinearSOE *theSOE = this->getLinearSOE();

  // Possibly set the independent part of the RHS
  if (assemblyFlag != 0)
    theSOE->setB(independentRHS);

  // Get the analysis model
  AnalysisModel *theModel = this->getAnalysisModel();

  //
  // Randomness in external load (including randomness in time series)
  //

  Domain *theDomain = theModel->getDomainPtr();

  // Loop through nodes to zero the unbalaced load
  Node *nodePtr;
  NodeIter &theNodeIter = theDomain->getNodes();
  while ((nodePtr = theNodeIter()) != nullptr)
    nodePtr->zeroUnbalancedLoad();

  // Loop through load patterns to add external load sensitivity
  LoadPattern *loadPatternPtr;
  LoadPatternIter &thePatterns = theDomain->getLoadPatterns();

  while ((loadPatternPtr = thePatterns()) != nullptr)
    loadPatternPtr->applyLoadSensitivity(theDomain->getCurrentTime());


  // Randomness in element/material contributions
  // Loop through FE elements
  FE_Element *elePtr;
  FE_EleIter &theEles = theModel->getFEs();    
  while ((elePtr = theEles()) != nullptr) {
    theSOE->addB(  elePtr->getResidual(this),  elePtr->getID()  );
  }

  // Loop through DOF groups (IT IS IMPORTANT THAT THIS IS DONE LAST!)
  DOF_Group *dofPtr;
  DOF_GrpIter &theDOFs = theModel->getDOFs();
  while ((dofPtr = theDOFs()) != nullptr) {
    theSOE->addB(  dofPtr->getUnbalance(this),  dofPtr->getID()  );
  }

  // Reset the sensitivity flag
  sensitivityFlag = 0;

  return 0;
}

int 
GeneralizedNewmark::formIndependentSensitivityRHS()
{
    // For now; don't use this
/*
  sensitivityFlag = 2; // Tell subsequent methods what to be assembled

  // Get pointer to the SOE
  LinearSOE *theSOE = this->getLinearSOEPtr();

  
  // Loop through FE elements
  FE_Element *elePtr;
  FE_EleIter &theEles = theModel->getFEs();    
  while((elePtr = theEles()) != 0) {
  theSOE->addB(  elePtr->getResidual(this),  elePtr->getID()  );
  }


  // Loop through DOF groups (IT IS IMPORTANT THAT THIS IS DONE LAST!)
  DOF_Group *dofPtr;
  DOF_GrpIter &theDOFs = theModel->getDOFs();
  while((dofPtr = theDOFs()) != 0) {
  theSOE->addB(  dofPtr->getUnbalance(this),  dofPtr->getID()  );
  }


  // Set the data member of this class
  independentRHS = theSOE->getB();


  // Reset the sensitivity flag
  sensitivityFlag = 0;
*/

    return 0;
}

int 
GeneralizedNewmark::saveSensitivity(const Vector & vNew,int gradNum,int numGrads)
{
    // Compute GeneralizedNewmark parameters
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
  DOF_GrpIter &theDOFs = myModel->getDOFs();
  DOF_Group *dofPtr;
  while ((dofPtr = theDOFs()) != nullptr) {

    const ID &id = dofPtr->getID();
    int idSize = id.Size();
    const Vector &dispSens = dofPtr->getDispSensitivity(gradNumber);
    for (int i = 0; i < idSize; i++)
      if (int loc = id(i); loc >= 0)
        dUn(loc) = dispSens(i);

    const Vector &velSens = dofPtr->getVelSensitivity(gradNumber);
    for (int i = 0; i < idSize; i++)
      if (int loc = id(i); loc >= 0)
        dVn(loc) = velSens(i);

    const Vector &accelSens = dofPtr->getAccSensitivity(gradNumber);
    for (int i = 0; i < idSize; i++)
      if (int loc = id(i); loc >= 0)
        dAn(loc) = accelSens(i);
  }

    // Compute new acceleration and velocity vectors:
    Vector vdotNew(vectorSize);
    Vector vdotdotNew(vectorSize);
    //(*vdotdotNewPtr) = vNew*c3 + V*a2 + Vdot*a3 + Vdotdot*a4;
    vdotdotNew.addVector(0.0, vNew, ca);
    vdotdotNew.addVector(1.0, dUn, a2);
    vdotdotNew.addVector(1.0, dVn, a3);
    vdotdotNew.addVector(1.0, dAn, a4);
    
    //(*vdotNewPtr) = vNew*a5 + V*bvu + Vdot*bvv + Vdotdot*bva;
    vdotNew.addVector(0.0, vNew, a5);
    vdotNew.addVector(1.0, dUn, bvu);
    vdotNew.addVector(1.0, dVn, bvv);
    vdotNew.addVector(1.0, dAn, bva);

    // update
    dUn = vNew;
    dVn = vdotNew;
    dAn = vdotdotNew;

    // Now we can save vNew, vdotNew and vdotdotNew
    DOF_GrpIter &theDOFGrps = myModel->getDOFs();
    DOF_Group   *dofPtr1;
    while ((dofPtr1 = theDOFGrps()) != nullptr)
      dofPtr1->saveSensitivity(vNew,vdotNew,vdotdotNew,gradNum,numGrads);

    return 0;
}

int 
GeneralizedNewmark::commitSensitivity(int gradNum, int numGrads)
{
    // Loop through the FE_Elements and set unconditional sensitivities
    AnalysisModel *theAnalysisModel = this->getAnalysisModel();
    FE_Element *elePtr;
    FE_EleIter &theEles = theAnalysisModel->getFEs();    
    while ((elePtr = theEles()) != nullptr)
      elePtr->commitSensitivity(gradNum, numGrads);

    return 0;
}


double
GeneralizedNewmark::getCFactor() 
{
  return cv;
}


int 
GeneralizedNewmark::computeSensitivities()
{

  LinearSOE *theSOE = this->getLinearSOE();

  // Zero out the old right-hand side of the SOE
  theSOE->zeroB();
  
  // Form the part of the RHS which are indepent of parameter
  this->formIndependentSensitivityRHS();
  AnalysisModel *theModel = this->getAnalysisModel();  //Abbas 
  Domain *theDomain=theModel->getDomainPtr();//Abbas
  ParameterIter &paramIter = theDomain->getParameters();

  Parameter *theParam;
  // De-activate all parameters
  while ((theParam = paramIter()) != nullptr)
    theParam->activate(false);
  
  // Now, compute sensitivity wrt each parameter
  int numGrads = theDomain->getNumParameters();

  paramIter = theDomain->getParameters();
  
  while ((theParam = paramIter()) != nullptr) {
    
    // Activate this parameter
    theParam->activate(true);
    
    // Zero the RHS vector
    theSOE->zeroB();
    
    // Get the grad index for this parameter
    int gradIndex = theParam->getGradIndex();

    // Form the RHS
    this->formSensitivityRHS(gradIndex);
    
    // Solve for displacement sensitivity 
    theSOE->solve();

    // Save sensitivity to nodes
    this->saveSensitivity( theSOE->getX(), gradIndex, numGrads );
    

    // Commit unconditional history variables (also for elastic problems; strain sens may be needed anyway)
    this->commitSensitivity(gradIndex, numGrads);
    
    // De-activate this parameter for next parameter sensitivity
    theParam->activate(false);
  }
  
  return 0;
}

