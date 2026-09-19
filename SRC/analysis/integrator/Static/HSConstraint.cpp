//===============================================================================
//# COPYRIGHT (C): Woody's license (by BJ):
//                 ``This    source  code is Copyrighted in
//                 U.S.,  for  an  indefinite  period,  and anybody
//                 caught  using it without our permission, will be
//                 mighty good friends of ourn, cause we don't give
//                 a  darn.  Hack it. Compile it. Debug it. Run it.
//                 Yodel  it.  Enjoy it. We wrote it, that's all we
//                 wanted to do.''
//
//# PROJECT:           Object Oriented Finite Element Program
//# PURPOSE:           Hyper-spherical Constraint
//# CLASS:             HSConstraint
//#
//# VERSION:           0.61803398874989 (golden section)
//# LANGUAGE:          C++
//# TARGET OS:         all...
//# DESIGN:            Ritu Jain, Boris Jeremic
//# PROGRAMMER(S):     Ritu, Boris Jeremic
//#
//#
//# DATE:              14Mar2003
//# UPDATE HISTORY:
//#
//#
//===============================================================================
//
#include <Logging.h>
#include <HSConstraint.h>
#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <Vector.h>
#include <math.h>


HSConstraint::HSConstraint(double arcLength, double psi_u, double psi_f, double u_ref)
: StaticIntegrator(INTEGRATOR_TAGS_HSConstraint),
  arcLength2(arcLength*arcLength), /*alpha2(alpha*alpha),*/
  psi_u2(psi_u*psi_u), 
  psi_f2(psi_f*psi_f),
  u_ref2(u_ref*u_ref),
  deltaUhat(0), 
  deltaUbar(0), 
  deltaU(0), 
  deltaUstep(0),
  phat(0), 
  deltaLambdaStep(0.0), 
  currentLambda(0.0),
  signLastDeltaLambdaStep(1)
{

}

HSConstraint::~HSConstraint()
{
    // delete any vector object created
    if (deltaUhat != 0)
        delete deltaUhat;
    if (deltaU != 0)
        delete deltaU;
    if (deltaUstep != 0)
        delete deltaUstep;
    if (deltaUbar != 0)
        delete deltaUbar;
    if (phat != 0)
        delete phat;
}

int
HSConstraint::newStep()
{
    // get pointers to AnalysisModel and LinearSOE
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();
    if (theModel == nullptr || theLinSOE == nullptr) {
        opserr << "WARNING HSConstraint::newStep() ";
        opserr << "No AnalysisModel or LinearSOE has been set\n";
        return -1;
    }

    // get the current load factor
    currentLambda = theModel->getCurrentDomainTime();

    if (deltaLambdaStep < 0)
        signLastDeltaLambdaStep = -1;
    else
        signLastDeltaLambdaStep = +1;



    // determine dUhat
    this->formTangent();
    theLinSOE->setB(*phat);
    theLinSOE->solve();

    (*deltaUhat) = theLinSOE->getX();
    Vector &dUhat = *deltaUhat;

    Vector f_ext = *phat;


    // determine delta lambda(1) == dlambda
//    double dLambda = sqrt(arcLength2/((dUhat^dUhat)+alpha2));
// out temp BJ 
//    double dLambda = sqrt(arcLength2/((psi_u2/u_ref2*fabs(dUhat^dUhat))+psi_f2));
// old version with fext
    double dLambda = std::sqrt(
                      arcLength2/( (psi_u2/u_ref2*std::fabs(dUhat^dUhat) ) + psi_f2*(f_ext^f_ext) ));
    dLambda *= signLastDeltaLambdaStep; // base sign of load change
                                        // on what was happening last step
    deltaLambdaStep = dLambda;
    currentLambda += dLambda;

    // determine delta U(1) == dU
    (*deltaU) = dUhat;
    (*deltaU) *= dLambda;
    (*deltaUstep) = (*deltaU);

    // update model with delta lambda and delta U
    theModel->incrDisp(*deltaU);
    theModel->applyLoadDomain(currentLambda);
    theModel->updateDomain();

    return 0;
}

int
HSConstraint::update(const Vector &dU)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();
    if (theModel == 0 || theLinSOE == 0) {
        opserr << "WARNING ArcLength::update() ";
        opserr << "No AnalysisModel or LinearSOE has been set\n";
        return -1;
    }

    (*deltaUbar) = dU; // have to do this as the SOE is gonna change

    // determine dUhat
    theLinSOE->solve(*phat, *deltaUhat);

    Vector f_ext = *phat;

    // determine the coeeficients of our quadratic equation
    double a1 = 
           psi_u2/u_ref2*((*deltaUhat)^(*deltaUhat)) 
           + 
           psi_f2 * (f_ext^f_ext);
    
    double a2 = 2.0 *(
           psi_u2/u_ref2*(((*deltaUhat)^(*deltaUbar))+((*deltaUhat)^(*deltaUstep)))
           +
           psi_f2*deltaLambdaStep  * (f_ext^f_ext));
    
    double a3 = 
           psi_u2/u_ref2 * (((*deltaUstep)+(*deltaUbar))^((*deltaUstep)+(*deltaUbar)))
           - 
           arcLength2 
           + 
           (deltaLambdaStep*deltaLambdaStep)*psi_f2 * (f_ext^f_ext) ;

    // check for a solution to quadratic
    double b24ac = a2*a2 - a1*a3;
    if (b24ac < 0) {
      opserr << "HSConstraint::update() - imaginary roots due to multiple instability";
      opserr << " directions - initial load increment was too large\n";
      opserr << "a1: " << a1 << " a2: " << a2 << " a3: " << a3 << " b24ac: " << b24ac << endln;
      return -1;
    }
    double dLambda;
    if (a1 == 0.0) {
     // opserr << "HSConstraint::update() - zero denominator";
     // opserr << "\n";
     // return -2;
        dLambda = -a3/(2.0*a2);
    }
    else
    {
        // determine the roots of the quadratic
        double sqrtb24ac = sqrt(b24ac);
        double dlambda1 = (-a2 + sqrtb24ac)/a1;
        double dlambda2 = (-a2 - sqrtb24ac)/a1;

        //Vector deltaU1 = (*deltaUbar);
        //deltaU1->addVector(1.0, *deltaUhat,dlambda1);
        //double costheta1 = (*deltaUstep)^((*deltaUstep)+(*deltaU1));

        //Vector deltaU2 = (*deltaUbar);
        //deltaU2->addVector(1.0, *deltaUhat,dlambda2);
        //double costheta2 = (*deltaUstep)^((*deltaUstep)+(*deltaU2));

        double val = (*deltaUhat)^(*deltaUstep);
            double costheta1 = ((*deltaUstep)^(*deltaUstep)) + ((*deltaUbar)^(*deltaUstep));
            double costheta2 = costheta1 + dlambda2*val;

            costheta1 += dlambda1*val;

            // choose dLambda based on angle between incremental displacement before
            // and after this step -- want positive
            /*double dLambda;*/
            if (costheta1 > costheta2)
                dLambda = dlambda1;
            else
                dLambda = dlambda2;
    }


    // determine delta U(i)
    (*deltaU) = (*deltaUbar);
    deltaU->addVector(1.0, *deltaUhat,dLambda);

    // update dU and dlambda
    (*deltaUstep) += *deltaU;
    deltaLambdaStep += dLambda;
    currentLambda += dLambda;

    // update the model
    theModel->incrDisp(*deltaU);
    theModel->applyLoadDomain(currentLambda);
    theModel->updateDomain();

    // set the X soln in linearSOE to be deltaU for convergence Test
    theLinSOE->setX(*deltaU);

    return 0;
}



int
HSConstraint::domainChanged()
{
    // we first create the Vectors needed
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();
    if (theModel == 0 || theLinSOE == 0) {
        opserr << "WARNING HSConstraint::domainChanged() ";
        opserr << "No AnalysisModel or LinearSOE has been set\n";
        return -1;
    }
    int size = theModel->getNumEqn(); // ask model in case N+1 space

    if (deltaUhat == 0 || deltaUhat->Size() != size) { // create new Vector
        if (deltaUhat != 0)
            delete deltaUhat;   // delete the old
        deltaUhat = new Vector(size);
        if (deltaUhat == 0 || deltaUhat->Size() != size) { // check got it
            opserr << "FATAL HSConstraint::domainChanged() - ran out of memory for";
            opserr << " deltaUhat Vector of size " << size << endln;
            return -1;
        }
    }

    if (deltaUbar == 0 || deltaUbar->Size() != size) { // create new Vector
        if (deltaUbar != 0)
            delete deltaUbar;   // delete the old
        deltaUbar = new Vector(size);
    }


    if (deltaU == 0 || deltaU->Size() != size) { // create new Vector
        if (deltaU != 0)
            delete deltaU;   // delete the old
        deltaU = new Vector(size);
    }

    if (deltaUstep == 0 || deltaUstep->Size() != size) {
        if (deltaUstep != 0)
            delete deltaUstep;
        deltaUstep = new Vector(size);
    }

    if (phat == 0 || phat->Size() != size) {
        if (phat != 0)
            delete phat;
        phat = new Vector(size);
    }

    // now we have to determine phat
    // do this by incrementing lambda by 1, applying load
    // and getting phat from unbalance.
    currentLambda = theModel->getCurrentDomainTime();
    currentLambda += 1.0;
    theModel->applyLoadDomain(currentLambda);
    this->formUnbalance(*phat); // NOTE: this assumes unbalance at last was 0
    currentLambda -= 1.0;
    theModel->setCurrentDomainTime(currentLambda);


    // check there is a reference load
    int haveLoad = 0;
    for (int i=0; i<size; i++)
      if ( (*phat)(i) != 0.0 ) {
        haveLoad = 1;
        i = size;
      }

    if (haveLoad == 0) {
      opserr << "WARNING HSConstraint::domainChanged() - zero reference load";
      return -1;
    }

    return 0;
}


void
HSConstraint::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel != 0) {
        double cLambda = theModel->getCurrentDomainTime();
        s << "\t HSConstraint - currentLambda: " << cLambda;
        s << "  HSConstraint: " << sqrt(arcLength2) /*<<  "  alpha: ";
        s << sqrt(alpha2) */ << endln;
    } else
        s << "\t HSConstraint - no associated AnalysisModel\n";
}
