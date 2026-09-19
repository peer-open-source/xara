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
#include <SecantLineSearch.h>
#include <IncrementalResidual.h>
#include <Vector.h>
#include <cmath>

SecantLineSearch::SecantLineSearch(double tol, int mIter, double mnEta, double mxEta, int pFlag)
:LineSearch(LINESEARCH_TAGS_SecantLineSearch),
 tolerance(tol), maxIter(mIter), minEta(mnEta), maxEta(mxEta), printFlag(pFlag)
{   

}

SecantLineSearch::~SecantLineSearch()
{

}


int 
SecantLineSearch::newStep(const Vector &Go)
{
  return 0;
}

int 
SecantLineSearch::search(double s0, 
                          double s1, 
                          const Vector& dU,
                          Vector& G,
                          Vector& Xs,
                          IncrementalResidual &theIntegrator)
{
  double r0 = 0.0;
  if ( s0 != 0.0 ) 
    r0 = std::fabs( s1 / s0 );
        
  if  (r0 <= tolerance )
    return 0; // Line Search Not Required Residual Decrease Less Than Tolerance

  if (s1 == s0)
    return 0;  // Secant will have a divide-by-zero if continue

  // set some variables
  double eta    = 1.0;
  double s      = s1;
  double etaJ   = 1.0;
  double etaJm1 = 0.0;
  double sJ     = s1;
  double sJm1   = s0;
  double r = r0;

  // Xo = theSOE.getX();
  // const Vector &dU = Xo;

  if (printFlag == 0) {
    opserr << "Secant Line Search - initial: "
         << "      eta(0) : " << eta << " , Ratio |s/s0| = " << r0 << "\n";
  }

  // perform the secant iterations:
  //
  //                eta(j+1) = eta(j) -  s(j) * (eta(j-1)-eta(j))
  //                                     ------------------------
  //                                           s(j-1) - s(j)
  Xs = dU;
  int count = 0; //initial value of iteration counter 
  while ( r > tolerance  &&  count < maxIter ) {
    
    count++;

    eta = etaJ - sJ * (etaJm1-etaJ) / (sJm1 - sJ);

    //-- want to put limits on eta and stop solution blowing up
    if (eta > maxEta)  eta = maxEta;
    if (r   > r0    )  eta =  1.0;
    if (eta < minEta)  eta = minEta;
    
    // update the incremental difference in response and determine new unbalance
    if (eta == etaJ) 
      break; // no change in response

    // set Xs = dU * (eta-etaJ)
    Xs.addVector(0, dU, eta-etaJ);
            
    if (theIntegrator.update(Xs) < 0)
      return -1;

    G.Zero();
    if (theIntegrator.formUnbalance(G) < 0)
      return -2;

    // new value of s
    s = Xs ^ G;

    // new value of r 
    r = std::fabs( s / s0 ); 

    if (printFlag == 0) {
      opserr << "Secant Line Search - iteration: " << count 
           << " , eta(j) : " << eta << " , Ratio |sj/s0| = " << r << "\n";
    }

    if (etaJ == eta)
      count = maxIter;

    // set variables for next iteration
    etaJm1 = etaJ;
    etaJ = eta;
    sJm1 = sJ;
    sJ = s;

    if (sJm1 == sJ)
      count = maxIter;
  } // end while

  // set X in the SOE for the revised dU, needed for convergence tests

  if (eta == 0.0)
    eta = 1.0;
    
  Xs.addVector(0, dU, eta);
  // theSOE.setX(Xs);

  return 0;
}


void
SecantLineSearch::Print(OPS_Stream &s, int flag)
{
  if (flag == 0) {
    s << "SecantLineSearch :: Line Search Tolerance = " << tolerance << "\n"; 
    s << "                       max num Iterations = " << maxIter << "\n";
    s << "                         max value on eta = " << maxEta << "\n";
  }
}









