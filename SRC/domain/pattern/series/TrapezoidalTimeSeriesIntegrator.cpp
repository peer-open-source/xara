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
// Description: This file contains the class definition for 
// a TrapezoidalTimeSeriesIntegrator, which integrates a
// ground motion TimeSeries using the trapezoidal rule.
//
// Written: MHS
// Created: 10/99
// Revision: A
//
#include <TrapezoidalTimeSeriesIntegrator.h>
#include <Vector.h>
#include <Channel.h>
#include <PathSeries.h>

TrapezoidalTimeSeriesIntegrator::TrapezoidalTimeSeriesIntegrator() 
  :TimeSeriesIntegrator(TIMESERIES_INTEGRATOR_TAG_Trapezoidal)
{

}

TrapezoidalTimeSeriesIntegrator::~TrapezoidalTimeSeriesIntegrator()
{

}


TimeSeries*
TrapezoidalTimeSeriesIntegrator::integrate(TimeSeries *theSeries, double delta)
{	
  // Check for zero time step, before dividing to get number of steps
  if (delta <= 0.0) {
    opserr << "TrapezoidalTimeSeriesIntegrator::integrate() Attempting to integrate using time step" <<
      delta << "<= 0\n";
    return 0;
   }

  // check a TimeSeries object was passed
  if (theSeries == 0) {
    opserr << "TrapezoidalTimeSeriesIntegrator::integrate() - - no TimeSeries passed\n";
    return 0;
  }

  // Add one to get ceiling out of type cast
  long long numSteps = (long long)(theSeries->getDuration() / delta + 1.0);

  Vector *theInt = new Vector (numSteps);


  double dummyTime = theSeries->getStartTime();

  double F = 0.0; // integral value
  double fi = 0.0;

  for (long long i = 0; i < numSteps; i++, dummyTime += delta) {
    double fj = theSeries->getFactor(dummyTime);

    // Apply the trapezoidal rule to update the integral
    F = F + 0.5 * delta * (fi + fj);

    (*theInt)[i] = F;

    fi = fj;
  }

  // Set the method return value
  PathSeries *returnSeries = new PathSeries (0, *theInt, delta, 1.0, true, false, theSeries->getStartTime());
  delete theInt;

  return returnSeries;
}

TimeSeries*
TrapezoidalTimeSeriesIntegrator::differentiate(TimeSeries *theSeries, double delta)
{	
  // Check for zero time step, before dividing to get number of steps
  if (delta <= 0.0) {
    opserr << "TrapezoidalTimeSeriesIntegrator::differentiate() Attempting to differentiate using time step" <<
      delta << "<= 0\n";
    return 0;
   }

  // check a TimeSeries object was passed
  if (theSeries == 0) {
    opserr << "TrapezoidalTimeSeriesIntegrator::differentiate() - - no TimeSeries passed\n";
    return 0;
  }

  // Add one to get ceiling out of type cast
  long long numSteps = (long long)(theSeries->getDuration() / delta + 1.0);

  Vector *theDif = new Vector (numSteps);

  double f = 0.0;       // derivative value
  
  // Dummy variable for integrating
  double dummyTime = theSeries->getStartTime();


  double Fi = 0.0;

  //opserr<<"differentiate()\n";
  for (long long i = 0; i < numSteps; i++, dummyTime += delta) {
    double Fj = theSeries->getFactor(dummyTime);

    // Apply the trapezoidal rule to update the derivative
    f = 2.0 * (Fj - Fi) / delta - f;

    (*theDif)[i] = f;

    Fi = Fj;
  }

  // Set the method return value
  PathSeries *returnSeries = new PathSeries(0, *theDif, delta, 1.0, true, false, theSeries->getStartTime());
  delete theDif;

  return returnSeries;
}


int
TrapezoidalTimeSeriesIntegrator::sendSelf(int commitTag, Channel &theChannel)
{
  return 0;
}

int
TrapezoidalTimeSeriesIntegrator::recvSelf(int commitTag, Channel &theChannel, 
		           FEM_ObjectBroker &theBroker)
{
  return 0;
}

void
TrapezoidalTimeSeriesIntegrator::Print(OPS_Stream &s, int flag)
{
   // Need to implement, return for now
   return;
}
