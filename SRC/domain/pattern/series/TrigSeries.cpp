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
// Purpose: This file contains the class implementation of TrigSeries.
//
// Written: fmk 
// Created: 07/99
// Revision: A
//
#include <TrigSeries.h>
#include <Vector.h>
#include <Channel.h>
#include <classTags.h>

#include <string.h>
#include <math.h>


TrigSeries::TrigSeries(int tag,
    double startTime, 
    double finishTime,
    double T, 
    double phaseshift, 
    double theFactor,
    double zeroshift)
    : TimeSeries(tag, TSERIES_TAG_TrigSeries),
    tStart(startTime), tFinish(finishTime),
    period(T), phaseShift(phaseshift),
    cFactor(theFactor), zeroShift(zeroshift)
{

}


TrigSeries::TrigSeries()
    : TimeSeries(TSERIES_TAG_TrigSeries),
    tStart(0.0), tFinish(0.0),
    period(1.0), phaseShift(0.0),
    cFactor(1.0), zeroShift(0.0)
{
    // does nothing
}


TrigSeries::~TrigSeries()
{
    // does nothing
}

TimeSeries *TrigSeries::getCopy()
{
    return new TrigSeries(this->getTag(), tStart, tFinish, period,
        phaseShift, cFactor, zeroShift);
}


double TrigSeries::getFactor(double pseudoTime)
{
    static double twopi = 4*asin(1.0);

    if (pseudoTime >= tStart && pseudoTime <= tFinish)  {
        double phi = phaseShift - period/twopi*asin(zeroShift/cFactor);
        return cFactor*sin(twopi*(pseudoTime-tStart)/period + phi) + zeroShift;
    }
    else
        return 0.0;
}


int TrigSeries::sendSelf(int commitTag, Channel &theChannel)
{
    int dbTag = this->getDbTag();
    Vector data(6);
    data(0) = cFactor;
    data(1) = tStart;
    data(2) = tFinish;
    data(3) = period;
    data(4) = phaseShift;
    data(5) = zeroShift;

    int result = theChannel.sendVector(dbTag,commitTag, data);
    if (result < 0) {
        opserr << "TrigSeries::sendSelf() - channel failed to send data\n";
        return result;
    }

    return 0;
}


int TrigSeries::recvSelf(int commitTag, Channel &theChannel, 
    FEM_ObjectBroker &theBroker)
{
    int dbTag = this->getDbTag();
    Vector data(6);
    int result = theChannel.recvVector(dbTag,commitTag, data);
    if (result < 0) {
        opserr << "TrigSeries::recvSelf() - channel failed to receive data\n";
        cFactor    = 1.0;
        tStart     = 0.0;
        tFinish    = 0.0;
        period     = 1.0;
        phaseShift = 0.0;
        zeroShift  = 0.0;
        return result;
    }
    cFactor    = data(0);
    tStart     = data(1);
    tFinish    = data(2);
    period     = data(3);
    phaseShift = data(4);
    zeroShift  = data(5);

    return 0;
}


void TrigSeries::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_MATE_INDENT << "{";
    s << "\"name\": \"" << this->getTag() << "\", ";
    s << "\"type\": \"TrigSeries\", ";
    s << "\"factor\": " << cFactor << ", ";
    s << "\"tStart\": " << tStart << ", ";
    s << "\"tFinish\": " << tFinish << ", ";
    s << "\"period\": " << period << ", ";
    s << "\"phaseShift\": " << phaseShift << ", ";
    s << "\"zeroShift\": " << zeroShift;
    s << "}";

  }
  else {
    s << "Trig Series" << "\n";
    s << "\tFactor: " << cFactor << "\n";
    s << "\ttStart: " << tStart << "\n";
    s << "\ttFinish: " << tFinish << "\n";
    s << "\tPeriod: " << period << "\n";
    s << "\tPhase Shift: " << phaseShift << "\n";
    s << "\tZero Shift: " << zeroShift << "\n";
  }
}

