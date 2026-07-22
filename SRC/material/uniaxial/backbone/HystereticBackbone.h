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
// Written: MHS
// Created: Aug 2000
//
// Description: This file contains the interface for HystereticBackbone,
// which represents a backbone curve for hysteretic models.

#ifndef HystereticBackbone_h
#define HystereticBackbone_h

#include <TaggedObject.h>
#include <MovableObject.h>

class Information;
class Parameter;

class HystereticBackbone : public TaggedObject//, public MovableObject
{
 public:
  HystereticBackbone(int tag, int classTag);
  virtual ~HystereticBackbone();
  
  virtual double getStress(double strain) = 0;
  virtual double getTangent(double strain) = 0;
  virtual double getEnergy(double strain) = 0;
  
  virtual double getYieldStrain(void) = 0;
  
  virtual HystereticBackbone *getCopy() = 0;
  
  virtual int setVariable(char *argv) {return -1;}
  virtual int getVariable(int varID, double &theValue) {return -1;}

  virtual int setParameter(char **argv, int argc, Parameter &eleInformation);
  virtual int updateParameter(int responseID, Information &eleInformation);	
  
 // CMP: removed MovableObject inheritance to avoid shadowing of get/setVariable
 int sendSelf(int commitTag, Channel &) {return 0;} // cmp
 int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &) {return 0;} // cmp
 int getDbTag() {return 0;} // cmp
 int setDbTag(int dbTag) {return 0;} // cmp  
 int getClassTag() {return 0;} // cmp
};

extern bool OPS_addHystereticBackbone(HystereticBackbone *newComponent);
extern HystereticBackbone *OPS_getHystereticBackbone(int tag);
extern void OPS_clearAllHystereticBackbone(void);

#endif
