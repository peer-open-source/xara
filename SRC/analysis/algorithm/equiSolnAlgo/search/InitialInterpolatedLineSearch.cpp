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
// [1] M.A. Crissfield, Nonlinear Finite Element Analysis of Solid and Structures,
//     Wiley. 97
//
#include <InitialInterpolatedLineSearch.h>
#include <SolutionAlgorithm.h>
#include <IncrementalResidual.h>
#include <Vector.h>
#include <cmath>


InitialInterpolatedLineSearch::InitialInterpolatedLineSearch(double tol, int mIter, double mnEta,
                                                             double mxEta, int pFlag)
 : LineSearch(LINESEARCH_TAGS_InitialInterpolatedLineSearch)
 , tolerance(tol)
 , maxIter(mIter)
 , minEta(mnEta)
 , maxEta(mxEta)
 , printFlag(pFlag)
{

}

InitialInterpolatedLineSearch::~InitialInterpolatedLineSearch()
{

}


int 
InitialInterpolatedLineSearch::newStep(const Vector &Go)
{
  return 0;
}

int 
InitialInterpolatedLineSearch::search(double s0, 
                                      double s1,
                                      const Vector& dU,
                                      Vector& G,
                                      Vector& Xs,
                                      IncrementalResidual &theIntegrator)
{

  // initialize r = ratio of residuals 
  double r0 = 0.0;
  if ( s0 != 0.0 ) 
    r0 = std::fabs( s1 / s0 );

  if (r0 <= tolerance ) {
    // Line Search Not Required Residual Decrease Less Than Tolerance
    return 0;
  }

  //
  // 1) Initialize search
  //
  double r = r0;
  double s = s1;
  double eta = 1.0;     // initial value of line search parameter
  double etaPrev = 1.0;

  // Xo = theSOE.getX();
  // const Vector &dU = Xo;


  if (printFlag == 0) {
    opserr << "           Line Search: " << 0
         << "    eta : " << eta 
         << " , Ratio |s/s0| = " << r0 
         << "\n";
  }


  // Solution procedure follows the one in Crissfields book [1].
  // NOTE: it is not quite linear interpolation/false-position/regula-falsi as eta(0) = 0.0
  // does not change. uses eta(i) = eta(i-1)*s0
  //                                -----------
  //                                s0 - s(i-1)  to compute eta(i)


  int count = 0;
  while ( r > tolerance  &&  count < maxIter ) {

    count++;

    eta = eta * s0 / (s0 - s);
    // eta = eta * s0 / (s - s0);

    // Put limits on eta(i)
    if (eta > maxEta)  eta = maxEta;
    if (r   > r0    )  eta =  1.0;
    if (eta < minEta)  eta = minEta;

    if (eta == etaPrev)
      break; // no change in response break

    // dx = ( eta * dx0 );
    Xs.addVector(0, dU, eta-etaPrev);

    if (theIntegrator.update(Xs) < 0)
      return SolutionAlgorithm::BadStepUpdate;

    if (theIntegrator.formUnbalance(G) < 0)
      return SolutionAlgorithm::BadFormResidual;

    // new value of s
    s = dU ^ G;

    // new value of r 
    r = std::fabs( s / s0 ); 

    if (printFlag == 0) {
      opserr << "           Line Search: " << count 
            << " ,  eta: " << eta
            << " , |s/s0| = " << r 
            << "\n";
    }

    // reset the variables, also check not just hitting bounds over and over
    if (eta == etaPrev)
      count = maxIter;
    else
      etaPrev = eta;

  } // end while

  //
  // set final X in the SOE for the revised dU, needed for convergence tests
  //
  if (eta == 0.0)
    eta = 1.0;

  Xs.addVector(0, dU, eta);

  return 0;
}



void
InitialInterpolatedLineSearch::Print(OPS_Stream &s, int flag)
{
  if (flag == 0)
    s << "InitialInterpolatedLineSearch :: Line Search Tolerance = " << tolerance << "\n"; 
}

