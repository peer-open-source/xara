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
// Written: Ed Love
// Created: 06/01
//
#include <BFGS.h>
#include <IncrementalIntegrator.h>
#include <LinearSOE.h>
#include <ConvergenceTest.h>
#include <ID.h>
#include <LineSearch.h>


BFGS::BFGS(IncrementalIntegrator::TangentFlagType theTangentToUse, int n , LineSearch* theSearch)
:EquiSolnAlgo(EquiALGORITHM_TAGS_BFGS),
 tangent(theTangentToUse), numberLoops(n),
 du(0),
 b(0),
 Go(0),
 Gn(0),
 search(theSearch),
 action(n, this)
{
  theTest = nullptr;

  s  = new Vector*[numberLoops+3]{};
  z  = new Vector*[numberLoops+3]{};
}


BFGS::~BFGS()
{
  for (int i =0; i < numberLoops+3; i++ ) {
    if (s[i] != nullptr )
      delete s[i];
    if (z[i] != nullptr )
      delete z[i];
    s[i] = nullptr;
    z[i] = nullptr;
  }

  if ( s != nullptr ) 
    delete[] s; 
  if ( z != nullptr ) 
    delete[] z;

  s = nullptr;
  z = nullptr;
}



int 
BFGS::solveCurrentStep()
{
  // set up some pointers and check they are valid
  // NOTE this could be taken away if we set Ptrs as protecetd in superclass

  IncrementalIntegrator *theIntegrator = this->getIncrementalIntegratorPtr();

  LinearSOE  *theSOE = this->getLinearSOEptr();

  if ((theIntegrator == 0) || (theSOE == 0) || (theTest == 0)) {
    return SolutionAlgorithm::BadAlgorithm;
  }        


  if (theTest->start(*theSOE) < 0)
    return SolutionAlgorithm::BadTestStart;


  ConvergenceTest* localTest = theTest->getCopy(this->numberLoops);


  const int systemSize = theSOE->getNumEqn();
  Go.resize(systemSize);
  Gn.resize(systemSize);
  du.resize(systemSize);
  b.resize(systemSize);

  action.link(theIntegrator, theSOE);

  int result = ConvergenceTest::Continue;
  int count = 0;
  do {

    // Form the initial tangent
    if (theIntegrator->formTangent(tangent) < 0)
      return SolutionAlgorithm::BadFormTangent;

    // form the initial residual 
    if (theIntegrator->formUnbalance(b) < 0) {
      return SolutionAlgorithm::BadFormResidual;  
    }

    // solve
    if (theSOE->solve(b, du) < 0)
      return SolutionAlgorithm::BadLinearSolve;

    // update
    if (theIntegrator->update(du) < 0)
      return SolutionAlgorithm::BadStepUpdate;

    // initial displacement increment
    if (s[1] == nullptr)
      s[1] = new Vector(systemSize);
    else
      s[1]->resize(systemSize);

    *s[1] = du;

    Go.addVector(0.0, b, -1.0);


    // form the residual again
    if (theIntegrator->formUnbalance(b) < 0)
      return SolutionAlgorithm::BadFormResidual;

    theSOE->setB(b);
    localTest->start(*theSOE);

    int nBFGS = 1;
    do {

      // save residual
      Gn.addVector(0.0, b, -1.0);

      // solve
      if (theSOE->solve(b, du) < 0)
        return SolutionAlgorithm::BadLinearSolve;

      // BFGS modifications to du
      BFGSUpdate(theIntegrator, theSOE, du, b, nBFGS);

      // if (search != nullptr) {
      //   const double s0 = - (du ^ theSOE->getB());
      //   if (search->search(du, *theSOE, *theIntegrator) < 0) {
      //     return -1;// SolutionAlgorithm::BadLineSearch;
      //   }
      // }
  
      if ( theIntegrator->update(du) < 0 )
        return SolutionAlgorithm::BadStepUpdate;
      
      // increment broyden counter
      nBFGS += 1;

      // save displacement increment
      if (s[nBFGS] == nullptr )
        s[nBFGS] = new Vector(systemSize);
      else
        s[nBFGS]->resize(systemSize);

      *s[nBFGS] = du;

      // swap residuals
      Go = Gn;

      // form the residual again
      if (theIntegrator->formUnbalance(b) < 0)
        return SolutionAlgorithm::BadFormResidual;

      result = localTest->test(b, du); 
      
    } while (result == ConvergenceTest::Continue && nBFGS <= numberLoops);

    result = theTest->test(b, du);

    this->record(count++);

  } while (result == ConvergenceTest::Continue);


  if (result == ConvergenceTest::Failure)
    return SolutionAlgorithm::TestFailed;

  // if positive result we are returning what the convergence test returned
  // which should be the number of iterations
  return result;
}



int 
BFGS::BFGSUpdate(IncrementalIntegrator *theIntegrator, 
                 LinearSOE *theSOE, 
                 Vector &du, 
                 const Vector &b,
                 int nBFGS) 
{

  static constexpr double eps = 1.0e-16;

  int systemSize = theSOE->getNumEqn();
  std::vector<double>& rdotz = action.rdotz;
  std::vector<double>& sdotr = action.sdotr;
  Vector& temp = action.temp;

  // compute z
  //  theSOE->setB( (*r[nBFGS]) - (*r[nBFGS-1]) );
  //    theSOE->setB( (*residNew) - (*residOld) );
  temp.addVector(0.0, Gn, 1.0);
  temp.addVector(1.0, Go, -1.0);
  theSOE->setB(temp);


  if ( z[nBFGS] == nullptr ) 
    z[nBFGS] = new Vector(systemSize);


  if (theSOE->solve(temp, *z[nBFGS]) < 0)
    return SolutionAlgorithm::BadLinearSolve;

  //  *z[nBFGS] *= (-1.0);


  for (int i=1; i<=(nBFGS-1); i++ ) {

    if ( sdotr[i] < eps )
      break; 

    double fact1 = 1.0 + ( rdotz[i] / sdotr[i] );

    fact1 /= sdotr[i];

    double pdotb = (*s[i]) ^ ( theSOE->getB() );

    fact1 *= pdotb;

    //    *z[nBFGS] +=  fact1 * ( *s[i] );
    z[nBFGS]->addVector(1.0, *s[i], fact1);

    double bdotz = (*z[i])^(theSOE->getB());

    //    *z[nBFGS] -= (1.0/sdotr[i]) * 
    //             ( bdotz * (*s[i])   +  pdotb * (*z[i]) );   
    temp = *s[i];
    temp *= bdotz;
    temp /= sdotr[i];
    *z[nBFGS] -= temp;

    temp = *z[i];
    temp *= pdotb;
    temp /= sdotr[i];
    *z[nBFGS] -= temp;
 
  } // end for i


  //sdotr[nBFGS] = *s[nBFGS] ^ ( *residNew - *residOld );

  //rdotz[nBFGS] = *z[nBFGS] ^ ( *residNew - *residOld );   
  temp  = Gn;
  temp -= Go;

  sdotr[nBFGS] = *s[nBFGS] ^ (temp);

  rdotz[nBFGS] = *z[nBFGS] ^ (temp);


  // BFGS modifications to du
  for (int i=1; i<=nBFGS; i++ ) {

    if ( sdotr[i] < eps )
      break;

    double fact1 = 1.0 + ( rdotz[i] / sdotr[i] );

    fact1 /= sdotr[i];

    double sdotb = (*s[i]) ^ b;

    fact1 *= sdotb;

    //du +=  fact1 * ( *s[i] );
    du.addVector(1.0, *s[i], fact1);


    double bdotz = (*z[i]) ^ b;  

    //du -= (1.0/sdotr[i]) * 
    //             ( bdotz * (*s[i])   +  sdotb * (*z[i]) );
    du.addVector(1.0, *s[i], -bdotz/sdotr[i]);

    du.addVector(1.0, *z[i], -sdotb/sdotr[i]);

  } // end for i

  return 0;
}


void
BFGS::Print(OPS_Stream &s, int flag) const
{
  if (flag == 0) {
    s << "BFGS" << "\n";
    s << "  Number of Iterations = " << numberLoops << "\n";
  }
}
