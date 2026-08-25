/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
#ifndef Parameter_h
#define Parameter_h

#include <Information.h>
#include <TaggedObject.h>

class MovableObject;
class Channel;
class FEM_ObjectBroker;
class Domain;

class Parameter : public TaggedObject, public MovableObject
{
 public:
  Parameter(int tag, 
	    MovableObject *theObject,
	    const char **argv, 
	    int argc);

  Parameter(const Parameter &);
  Parameter(int tag, int classTag = PARAMETER_TAG_Parameter);
  Parameter();
  virtual ~Parameter();
  
  virtual void Print(OPS_Stream &s, int flag);
  
  virtual int update(int newValue); 
  virtual int update(double newValue);
  virtual int activate(bool active);
  virtual double getValue() {return theInfo.theDouble;}
  virtual void setValue(double newValue) {theInfo.theDouble = newValue;}

  virtual int addComponent(MovableObject *, const char **argv, int argc);  
  virtual int addComponent(int, const char **argv, int argc);  
  virtual int addObject(int parameterID, MovableObject *object);

  virtual int clean();

  void setGradIndex(int gradInd) {gradIndex = gradInd;}
  int getGradIndex() {return gradIndex;}

  virtual bool isImplicit() {return true;}
  virtual double getSensitivity(int index);
  virtual double getPerturbation() {return 0.001;}
  virtual const char *getType() {return "FEModel";}
  virtual int getPointerTag() {return -1;}

  virtual void setDomain(Domain *theDomain);

  virtual int sendSelf(int commitTag, Channel &);  
  virtual int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &);

 protected:
  int *parameterID;

  MovableObject **theObjects;
  int numObjects;
  int maxNumObjects;
  
 private:
  Information theInfo;
  double currentValue;

  static constexpr int initialSize = 64;
  static constexpr int expandSize = 128;

  MovableObject **theComponents;
  int numComponents;
  int maxNumComponents;

  int gradIndex; // 0,...,nparam-1
};

#endif
