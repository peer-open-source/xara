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
// Written: fmk 
// Created: 11/01
//
#include <RegulaFalsiLineSearch.h>
#include <IncrementalResidual.h>
#include <Vector.h>
#include <cmath>

RegulaFalsiLineSearch::RegulaFalsiLineSearch(double tol, int mIter, double mnEta, double mxEta, int pFlag)
: LineSearch(LINESEARCH_TAGS_RegulaFalsiLineSearch),
  tolerance(tol), maxIter(mIter), minEta(mnEta), maxEta(mxEta), printFlag(pFlag)
{   

}

RegulaFalsiLineSearch::~RegulaFalsiLineSearch()
{

}


int 
RegulaFalsiLineSearch::newStep(const Vector &Go)
{
  return 0;
}


int 
RegulaFalsiLineSearch::search(double s0, 
                              double s1, 
                              const Vector& dU,
                              Vector& G,
                              Vector& Xs,
                              IncrementalResidual &theIntegrator)
{
  // Initialize residual ratio
  double r0 = 0.0;
  if ( s0 != 0.0 ) 
    r0 = std::fabs( s1 / s0 );

  if (r0 <= tolerance )
    // Line Search Not Required Residual Decrease Less Than Tolerance
    return 0;

  if (s1 == s0)
    return 0;  // RegulaFalsi will have a divide-by-zero error if continue

  //
  // 1) Initialize search
  //
  double r      = r0;
  double s      = s1;
  double eta    = 1.0;
  double etaU   = 1.0;
  double etaL   = 0.0;
  double sU     = s1;
  double sL     = s0;
  double etaJ   = 1.0;
  double compoundFactor = 0.0;

  Xs = dU;

  if (printFlag == 0) {
    opserr << "        Line Search - initial: "
           << "      eta(0) : " 
           << eta << " , Ratio |s/s0| = " << r0 << "\n";
  }

  //
  // 2) Search for a bracket
  //

  // we first search for a bracket to a solution, i.e. we want sU * sL < 0.0
  int count = 0;
  while ((sU * sL > 0.0) && (etaU < maxEta)) {

    count++;

    /*
    if (count == 1)
      etaU = 0.5;
    else
    */
    etaU = etaJ * 4.0;
    if (etaU > maxEta)
      etaU = maxEta;

    // update the incremental difference in response and determine new unbalance
    double factor = etaU - etaJ;
    Xs.addVector(0, dU, factor);
    compoundFactor += factor;

    etaJ = etaU;

    // new value of sU
    if (theIntegrator.update(Xs) < 0)
      return -1;
    G.Zero();
    if (theIntegrator.formUnbalance(G) < 0)
      return -2;

    sU = dU ^ G;

    // check if we have a solution we are happy with
    r = std::fabs( sU / s0 ); 
    if (r < tolerance) {
      Xs.addVector(0, dU, etaJ);
      return 0;
    }

    if (printFlag == 0) {
      opserr << "Bisection Line Search - bracketing: " << count 
             << " , eta(j) : " << etaU << " , Ratio |sj/s0| = " << r << "\n";
    }
  }

  // return if no bracket for a solution found, reset to initial values
  if (sU * sL > 0.0) {
    Xs = dU;
    // theSOE.setX(Xs);
    Xs *= -compoundFactor;
    if (theIntegrator.update(Xs) < 0)
      return -1;
    if (theIntegrator.formUnbalance(G) < 0)
      return -2;
    return 0; 
  }

  //
  // 3) Perform the secant iterations:
  //
  //                eta(j+1) = eta(u) -  s(u) * (eta(l) -eta(u))
  //                                     ------------------------
  //                                           s(l) - s(u)

  count = 0; //initial value of iteration counter 
  while ( r > tolerance  &&  count < maxIter ) {
    
    count++;

    eta = etaU - sU * (etaL-etaU) / (sL - sU);


    // Put limits on eta(i)
    if (eta > maxEta)  eta = maxEta;
    if (  r >  r0   )  eta =  1.0;
    if (eta < minEta)  eta = minEta;

    if (eta == etaJ) // break if going to have a zero *x
      break;
    
    //update the incremental difference in response and determine new unbalance
    Xs = dU;
    Xs *= eta-etaJ;
            
    if (theIntegrator.update(Xs) < 0) {  
      return -1;
    }
    G.Zero();
    if (theIntegrator.formUnbalance(G) < 0) { 
      return -2;
    }

    // new value of s
    s = dU ^ G;

    // new value of r 
    r = std::fabs( s / s0 ); 


    if (printFlag == 0) {
      opserr << "RegulaFalsi Line Search - iteration: " << count 
           << " , eta(j) : " << eta << " , Ratio |sj/s0| = " << r << "\n";
    }

    if (etaJ == eta)
      count = maxIter;

    // set variables for next iteration
    etaJ = eta;
    
    if (s*sU < 0.0) {
      etaL = eta;
      sL   = s;
    } else if (s*sU == 0.0)
      count = maxIter;
    else {
      etaU = eta;
      sU   = s;
    } 

    if (sL == sU)
      count = maxIter;

  } // end while

  // set X in the SOE for the revised dU, needed for convergence tests

  if (eta == 0.0)
    eta = 1.0;

  Xs.addVector(0, dU, etaJ);
  // theSOE.setX(Xs);

  return 0;
}



void
RegulaFalsiLineSearch::Print(OPS_Stream &s, int flag)
{
  if (flag == 0) {
    s << "RegulaFalsiLineSearch :: Line Search Tolerance = " << tolerance << endln; 
    s << "                         max num Iterations = " << maxIter << endln;
    s << "                         max value on eta = " << maxEta << endln;
  }
}



