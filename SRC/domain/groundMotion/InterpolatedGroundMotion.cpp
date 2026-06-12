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
// InterpolatedGroundMotion. 
//
// Written: fmk 
// Created: 11/00
// Revision: A
// 
#include <InterpolatedGroundMotion.h>
#include <stdlib.h>
#include <math.h>
#include <classTags.h>
#include <Vector.h>
#include <Channel.h>

InterpolatedGroundMotion::InterpolatedGroundMotion()
:GroundMotion(GROUND_MOTION_TAG_InterpolatedGroundMotion),
 theMotions(0), factors(0),
 destroyMotions(0), data(3), deltaPeak(0.0)
{
    
}

InterpolatedGroundMotion::InterpolatedGroundMotion(GroundMotion **groundMotions, 
						   const Vector &fact,
						   bool destroyMotions,
						   double dT)
:GroundMotion(GROUND_MOTION_TAG_InterpolatedGroundMotion),
 theMotions(0), factors(0),
 destroyMotions(0), data(3), deltaPeak(dT)
{
  factors = new Vector(fact);
  theMotions = new GroundMotion *[fact.Size()];

  for (int i=0; i<fact.Size(); i++)
    theMotions[i] = groundMotions[i];
  
  if (destroyMotions == true)
    this->destroyMotions = 1;
}


InterpolatedGroundMotion::~InterpolatedGroundMotion()
{
  if (destroyMotions == 1) {
      for (int i=0; i<factors->Size(); i++)
	  delete theMotions[i];
  }  
  
  delete [] theMotions;
  delete factors;
}

double 
InterpolatedGroundMotion::getDuration()
{
  double value = 0.0;
  int numMotions = factors->Size();
  for (int i=0; i<numMotions; i++) {
      double motionValue = theMotions[i]->getDuration();
      if (motionValue > value)
	  value = motionValue;
  }
  return value;
}

double 
InterpolatedGroundMotion::getPeakAccel()
{
  double value = 0.0;
  double duration = this->getDuration();
  double time = 0.0;
  while (time < duration) {
    double accel = this->getAccel(time);
    if (accel > value)
      value = accel;
    time += deltaPeak;
  }
  return value;
}

double 
InterpolatedGroundMotion::getPeakVel()
{
  double value = 0.0;
  double duration = this->getDuration();
  double time = 0.0;
  while (time < duration) {
    double accel = this->getVel(time);
    if (accel > value)
      value = accel;
    time += deltaPeak;
  }
  return value;
}

double 
InterpolatedGroundMotion::getPeakDisp()
{
  double value = 0.0;
  double duration = this->getDuration();
  double time = 0.0;
  while (time < duration) {
    double accel = this->getDisp(time);
    if (accel > value)
      value = accel;
    time += deltaPeak;
  }
  return value;
}

double 
InterpolatedGroundMotion::getAccel(double time)
{
  if (time < 0.0)
    return 0.0;

  double value = 0.0;
  int numMotions = factors->Size();
  for (int i=0; i<numMotions; i++) {
      value += (*factors)(i) * theMotions[i]->getAccel(time);
  }

  return value;

}     

double 
InterpolatedGroundMotion::getVel(double time)
{
  if (time < 0.0)
    return 0.0;

  double value = 0.0;
  int numMotions = factors->Size();
  for (int i=0; i<numMotions; i++) {
      value += (*factors)(i) * theMotions[i]->getVel(time);
  }


  return value;
}

double 
InterpolatedGroundMotion::getDisp(double time)
{
  if (time < 0.0)
    return 0.0;

  double value = 0.0;
  int numMotions = factors->Size();
  for (int i=0; i<numMotions; i++) {
      value += (*factors)(i) * theMotions[i]->getDisp(time);
  }

  return value;
}

const Vector &
InterpolatedGroundMotion::getDispVelAccel(double time)
{
  if (time < 0.0) {
    data(0) = 0.0;
    data(1) = 0.0;
    data(2) = 0.0;
    return data;
  }

  data.Zero();
  static Vector motionData(3);

  int numMotions = factors->Size();
  for (int i=0; i<numMotions; i++) {
      motionData = theMotions[i]->getDispVelAccel(time);
      motionData *= (*factors)(i);
      data += motionData;
  }
  
  return data;
}


int 
InterpolatedGroundMotion::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}


int 
InterpolatedGroundMotion::recvSelf(int commitTag, Channel &theChannel,
				   FEM_ObjectBroker &theBroker)
{
  return -1;
}



