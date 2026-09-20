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
// Description: This file contains the class definition implementation of
// Broyden.  
//
// Written: Ed C++ Love
// Created: 04/01
//
#include <Broyden.h>
#include <IncrementalIntegrator.h>
#include <LinearSOE.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ConvergenceTest.h>
#include <ID.h>
#include <cmath>


Broyden::Broyden(int theTangentToUse, int n )
:EquiSolnAlgo(EquiALGORITHM_TAGS_Broyden),
 tangent(theTangentToUse), numberLoops(n),
 du(0),
 temp(0)
{
  s  = new Vector*[numberLoops+3] ;
  z  = new Vector*[numberLoops+3] ;

  residOld = nullptr ;
  residNew = nullptr ;

  for ( int i =0; i < numberLoops+3; i++ ) {
    s[i] = 0 ;
    z[i] = 0 ;
    //r[i] = 0 ;
  }
}


// Destructor
Broyden::~Broyden()
{
  if ( residOld != 0 )
    delete residOld ;  
  residOld = 0 ;

  if ( residNew != 0 ) 
    delete residNew ;
  residNew = 0 ;


  for (int i =0; i < numberLoops+3; i++ ) {
    if ( s[i] != 0 )
      delete s[i] ;
    if ( z[i] != 0 )
      delete z[i] ;
    s[i] = 0 ;
    z[i] = 0 ;
  } // end for i

  if ( s != 0 ) delete[] s ; 

  if ( z != 0 ) delete[] z ;
}


int 
Broyden::solveCurrentStep()
{
  // set up some pointers and check they are valid
  // NOTE this could be taken away if we set Ptrs as protecetd in superclass

  IncrementalIntegrator *theIntegrator = this->getIncrementalIntegratorPtr();

  LinearSOE  *theSOE = this->getLinearSOEptr();


  if ((theIntegrator == nullptr) 
    || (theSOE == nullptr)
    || (theTest == nullptr)) {
    return SolutionAlgorithm::BadAlgorithm;
  }        

  if (theTest->start(*theSOE) < 0) {
    return SolutionAlgorithm::BadTestStart;
  }

  ConvergenceTest *localTest = theTest->getCopy(this->numberLoops);

  const int systemSize = theSOE->getNumEqn();
  temp.resize(systemSize);
  du.resize(systemSize);

  // initial displacement increment
  if ( s[1] == nullptr ) 
    s[1] = new Vector(systemSize);

  int result = -1 ;
  int count = 0 ;
  do {

    // form the initial tangent
    if (theIntegrator->formTangent(tangent) < 0)
      return SolutionAlgorithm::BadFormTangent;


    // form the initial residual 
    if (theIntegrator->formUnbalance(temp) < 0) {
      opserr << "WARNING Broyden::solveCurrentStep() -";
      opserr << "the Integrator failed in formUnbalance\n";
    }            

    // solve
    if (theSOE->solve(temp, *s[1]) < 0)
      return SolutionAlgorithm::BadLinearSolve;

    // update
    if ( theIntegrator->update(*s[1]) < 0)
      return SolutionAlgorithm::BadStepUpdate;


    // initial residual
    if ( residOld == nullptr )
      residOld = new Vector(systemSize) ;

    *residOld = theSOE->getB();
    *residOld *= (-1.0 ) ;

    //form the residual again
    if (theIntegrator->formUnbalance(temp) < 0) {
      opserr << "WARNING Broyden::solveCurrentStep() -";
      opserr << "the Integrator failed in formUnbalance\n";        
    }            

    if ( residNew == nullptr ) 
      residNew = new Vector(systemSize) ;


    localTest->start(*theSOE) ;

    int nBroyden = 1 ;
    do {

      // save residual
      /*    if ( r[nBroyden] == 0 ) r[nBroyden] = new Vector(systemSize) ;
      *r[nBroyden] =  theSOE->getB( ) ; 
      *r[nBroyden] *= (-1.0 ) ; 
      */

      *residNew =  -1.0*temp;

      // solve
      if (theSOE->solve(temp, du) < 0)
        return SolutionAlgorithm::BadLinearSolve;

      // broyden modifications to du
      BroydenUpdate(*theIntegrator, *theSOE, du, nBroyden );

      if ( theIntegrator->update(du) < 0 )
        return SolutionAlgorithm::BadStepUpdate;
      
      //increment broyden counter
      nBroyden += 1 ;

      // save displacement increment
      if ( s[nBroyden] == 0 ) 
        s[nBroyden] = new Vector(systemSize) ;

      *s[nBroyden] = du ;

      // swap residuals
      *residOld = *residNew ;

      //form the residual again
      if (theIntegrator->formUnbalance(temp) < 0) {
        ;
      }            
      
      result = localTest->test(*theSOE); ;
      
    } while ( result == -1 && nBroyden <= numberLoops );


    result = theTest->test(*theSOE);
    this->record(count++);

  } while (result == ConvergenceTest::Continue);


  if (result == ConvergenceTest::Failure)
    return SolutionAlgorithm::TestFailed;


  // note - if positive result we are returning what the convergence test returned
  // which should be the number of iterations
  return result;
}



int
Broyden::BroydenUpdate( IncrementalIntegrator &theIntegrator, 
                        LinearSOE &theSOE, 
                        Vector &du, 
                        int nBroyden ) 
{

  static constexpr double eps = 1.0e-16 ;


  //compute z
  //  theSOE->setB( (*r[nBroyden]) - (*r[nBroyden-1]) ) ;
  //    theSOE->setB( (*residNew) - (*residOld) ) ;
  temp  = (*residNew);
  temp -= (*residOld);
  theSOE.setB(temp);

  if ( z[nBroyden] == 0 ) 
    z[nBroyden] = new Vector(du.Size());

  if (theSOE.solve(temp, *z[nBroyden]) < 0)
    return SolutionAlgorithm::BadLinearSolve;


  *z[nBroyden] *= (-1.0) ;


  for (int i=1; i<=(nBroyden-1); i++ ) {

    double p = - ( (*s[i]) ^ (*z[i]) ) ;

    if ( std::fabs(p) < eps ) break ;

    double sdotz = (*s[i]) ^ (*z[nBroyden]) ;

    //*z[nBroyden] += (1.0/p) * sdotz * ( *s[i] + *z[i] ) ;
    z[nBroyden]->addVector(1.0, *s[i], (1.0/p) * sdotz);
    z[nBroyden]->addVector(1.0, *z[i], (1.0/p) * sdotz);
  }


  // broyden modifications to du
  for (int i=1; i<=nBroyden; i++ ) {

    double p = - ( (*s[i]) ^ (*z[i]) ) ;

    if ( std::fabs(p) < eps )
      break ;

    double sdotdu = (*s[i]) ^ du ;

    //du += (1.0/p) * sdotdu * ( *s[i] + *z[i] ) ;
    du.addVector(1.0, *s[i], (1.0/p) * sdotdu);
    du.addVector(1.0, *z[i], (1.0/p) * sdotdu);

  }

  return 0;
}


void
Broyden::Print(OPS_Stream &s, int flag) const
{
  if (flag == 0) {
    s << "Broyden" << endln ;
    s << "  Number of Iterations = " << numberLoops << endln ;
  }
}


//const Vector &pi = *p[i] ;

//    //residual at this iteration before next solve 
//    Resid0 = theSOE->getB() ;

//    //line search direction 
//    dx0 = theSOE->getX() ;





/*

        //first solve step

         if (theIntegrator->formTangent(tangent) < 0){
            opserr << "WARNING Broyden::solveCurrentStep() -";
            opserr << "the Integrator failed in formTangent()\n";
            return -1;
        }                    
        
        if (theSOE->solve() < 0) {
            opserr << "WARNING Broyden::solveCurrentStep() -";
            opserr << "the LinearSysOfEqn failed in solve()\n";        
            return -3;
        }            


        if (theIntegrator->update(theSOE->getX()) < 0) {
            opserr << "WARNING Broyden::solveCurrentStep() -";
            opserr << "the Integrator failed in update()\n";        
            return -4;
        }                


        if (theIntegrator->formUnbalance() < 0) {
            opserr << "WARNING Broyden::solveCurrentStep() -";
            opserr << "the Integrator failed in formUnbalance\n";        
            return -2;
        }        
        
        result = theTest->test();
        this->record(nBroyden++);

      const Vector &du = BroydengetX( theIntegrator, theSOE, nBroyden )  ;

*/
