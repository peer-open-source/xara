//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, OpenSees/Xara Developers
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
//
//
// Written: cmp
// Adapted from BisectionLineSearch.cpp by fmk, dated 11/01
//
#include <BisectionLineSearch.h>
#include <IncrementalResidual.h>
#include <Vector.h>
#include <cmath>

BisectionLineSearch::BisectionLineSearch(double tol, int mIter, double mnEta, double mxEta, int pFlag)
 : LineSearch(LINESEARCH_TAGS_BisectionLineSearch)
 , tolerance(tol)
 , maxIter(mIter)
 , minEta(mnEta)
 , maxEta(mxEta)
 , printFlag(pFlag)
{   

}

BisectionLineSearch::~BisectionLineSearch()
{

}


int 
BisectionLineSearch::newStep(const Vector &Go)
{
  return 0;
}


int
BisectionLineSearch::search(double s0, 
                            double s1, 
                            const Vector& dU,
                            Vector& G,
                            Vector& Xs,
                            IncrementalResidual &theIntegrator)
{
  double r0 = 0.0;
  if ( s0 != 0.0 ) 
    r0 = std::fabs( s1 / s0 );

  Xs = dU;

  if  (r0 <= tolerance )
    return 0; // Line Search Not Required Residual Decrease Less Than Tolerance

  if (s1 == s0)
    return 0;  // Bisection will have a divide-by-zero error if continue

  // set some variables
  double eta    = 1.0;
  double s      = s1;
  double etaU   = 1.0;
  double etaL   = 0.0;
  double sU     = s1;
  double sL     = s0;
  double r      = r0;
  double etaJ   = 1.0;
  double compoundFactor = 0.0;


  if (printFlag == 0) {
    opserr << "           Line Search: " << 0
         << "    eta : " << eta 
         << " , Ratio |s/s0| = " << r0 
         << "\n";
  }

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

    // update the incremental difference in response and determine new unbalance
    double factor = etaU - etaJ;
    Xs.addVector(0, dU, factor);
    compoundFactor += factor;

    etaJ = etaU;

    if (theIntegrator.update(Xs) < 0)
      return -1;
    if (theIntegrator.formUnbalance(G) < 0)
      return -2;

    sU = dU ^ G;

    // check if we have a solution we are happy with
    r = std::fabs( sU / s0 ); 
    if (r < tolerance)
      return 0;

    if (printFlag == 0) {
      opserr << "           Line Search: " << count 
            << " ,  eta: " << eta
            << " , |s/s0| = " << r 
            << "\n";
    }
  }

  // return if no bracket for a solution found, resetting to initial values
  if (sU * sL > 0.0) {
    Xs = dU;
    Xs *= -compoundFactor;
    theIntegrator.update(Xs);
    theIntegrator.formUnbalance(G);
    return 0; 
  }

  // perform the secant iterations:
  //
  //                eta(j+1) = eta(l) + eta(u)
  //                           ---------------
  //                                2.0

  count = 0; //initial value of iteration counter 
  while ( r > tolerance  &&  count < maxIter ) {
    
    count++;

    eta = (etaU + etaL) * 0.5;

    //-- want to put limits on eta(i)
    //    if (r   > r0    )  eta =  1.0;
    
    // update the incremental difference in response and determine new unbalance
    Xs = dU;
    double fact = eta-etaJ;

    if (fact == 0)
      break;

    Xs *= fact;
            
    if (theIntegrator.update(Xs) < 0) { 
      return -1;
    }
    if (theIntegrator.formUnbalance(G) < 0) {
      return -2;
    }
    //new value of s
    s = Xs ^ G;
    
    //new value of r 
    r = std::fabs( s / s0 ); 

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

    if (printFlag == 0) {
      opserr << "Bisection Line Search - iteration: " << count 
           << " , eta(j) : " << eta << " , Ratio |sj/s0| = " << r << "\n";
    }
    
  } // end while

  // set X in the SOE for the revised dU, needed for convergence tests
  Xs = dU;
  if (eta != 0.0) 
    Xs *= eta;

  return 0;
}



void
BisectionLineSearch::Print(OPS_Stream &s, int flag)
{
  if (flag == 0) {
    s << "BisectionLineSearch :: Line Search Tolerance = " << tolerance << "\n";
    s << "                         max num Iterations = " << maxIter << "\n";
    s << "                         max value on eta = " << maxEta << "\n";
  }
}

