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
                                                                        
// $Revision: 1.3 $
// $Date: 2010-02-04 00:34:11 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/RectangularSeries.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 07/99
// Revision: A
//
// Purpose: This file contains the class implementation of RectangularSeries.
//
// What: "@(#) RectangularSeries.C, revA"

#include <string.h>
#include <RectangularSeries.h>
#include <Vector.h>
#include <Channel.h>


RectangularSeries::RectangularSeries(int tag, 
				     double startTime, 
				     double finishTime,
				     double theFactor)
  :TimeSeries(tag, TSERIES_TAG_RectangularSeries),
   tStart(startTime),tFinish(finishTime),cFactor(theFactor)
{
  // does nothing
}


RectangularSeries::RectangularSeries()
  :TimeSeries(TSERIES_TAG_RectangularSeries),
   tStart(0),tFinish(0),cFactor(0)
{
  // does nothing
}


RectangularSeries::~RectangularSeries()
{
  // does nothing
}

TimeSeries *
RectangularSeries::getCopy()
{
  return new RectangularSeries(this->getTag(), tStart, tFinish, cFactor);
}

double
RectangularSeries::getFactor(double pseudoTime)
{	
  if (pseudoTime >= tStart && pseudoTime <= tFinish)
    return cFactor;
  else
    return 0;
}

int
RectangularSeries::sendSelf(int commitTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();
  Vector data(3);
  data(0) = cFactor;
  data(1) = tStart;	
  data(2) = tFinish;
  int result = theChannel.sendVector(dbTag,commitTag, data);
  if (result < 0) {
    opserr << "RectangularSeries::sendSelf() - channel failed to send data\n";
    return result;
  }
  return 0;
}


int 
RectangularSeries::recvSelf(int commitTag, Channel &theChannel, 
		       FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();
  Vector data(3);
  int result = theChannel.recvVector(dbTag,commitTag, data);
  if (result < 0) {
    opserr << "RectangularSeries::sendSelf() - channel failed to receive data\n";
    cFactor = 1.0;
    tStart= 0.0;
    tFinish = 0.0;
    return result;
  }
  cFactor = data(0);
  tStart = data(1);
  tFinish = data(2);

  return 0;    
}

void
RectangularSeries::Print(OPS_Stream &s, int flag)
{
    s << "Linear Series: constant factor: " << cFactor;
    s << "  tStart: " << tStart << "  tFinish: " << tFinish << endln;

}
