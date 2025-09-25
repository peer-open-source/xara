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
#ifndef BeamIntegration_h
#define BeamIntegration_h

#include <MovableObject.h>
#include <TaggedObject.h>
#include <ID.h>

class Matrix;
class ElementalLoad;
class Information;

class BeamIntegration : public MovableObject
{
 public:
  BeamIntegration(int classTag);
  virtual ~BeamIntegration();

  virtual void getSectionLocations(int nIP, double L, double *xi) const = 0;
  virtual void getSectionWeights(int nIP, double L, double *wt) const = 0;

  virtual BeamIntegration *getCopy(void) = 0;

//
//
// (cmp: changed from virtual to inline; nothing was implementing?)
  inline  void addElasticDeformations(ElementalLoad *theLoad,
				      double loadFactor,
				      double L, double *v0) const {return;}
  // Return 0 if there is no elastic interior, -1 otherwise
  inline  int addElasticFlexibility(double L, Matrix &fe) const {return 0;}

  // Return 0 if there is no elastic interior, -1 otherwise
  inline  int addElasticFlexDeriv(double L, Matrix &dfedh,
				  double dLdh = 0.0) const {return 0;}
//
//
//
  virtual double getTangentDriftI(double L, double LI, double q2,
				  double q3, bool yAxis = false) {return 0.0;}
  virtual double getTangentDriftJ(double L, double LI, double q2,
				  double q3, bool yAxis = false) {return 0.0;}


  virtual void getLocationsDeriv(int nIP, double L, double dLdh, double *dptsdh);
  virtual void getWeightsDeriv(int nIP, double L, double dLdh, double *dwtsdh);


  virtual void Print(OPS_Stream &s, int flag = 0) = 0;

#if 0
  static void
  getData(int n, const double **const x, const double **const w) {
    *x = nullptr;
    *w = nullptr;
  }

  template<template<int, int> class Gauss> static void
  getData(int n, const double **const x, const double **const w)
  {
    switch (n) {
      case  1: *x = Gauss<1, 1>::pts; *w = Gauss<1, 1>::wts; break;
      case  2: *x = Gauss<1, 2>::pts; *w = Gauss<1, 2>::wts; break;
      case  3: *x = Gauss<1, 3>::pts; *w = Gauss<1, 3>::wts; break;
      case  4: *x = Gauss<1, 4>::pts; *w = Gauss<1, 4>::wts; break;
      case  5: *x = Gauss<1, 5>::pts; *w = Gauss<1, 5>::wts; break;
      case  6: *x = Gauss<1, 6>::pts; *w = Gauss<1, 6>::wts; break;
      case  7: *x = Gauss<1, 7>::pts; *w = Gauss<1, 7>::wts; break;
      case  8: *x = Gauss<1, 8>::pts; *w = Gauss<1, 8>::wts; break;
      case  9: *x = Gauss<1, 9>::pts; *w = Gauss<1, 9>::wts; break;
      case 10: *x = Gauss<1,10>::pts; *w = Gauss<1,10>::wts; break;
    }
  }
#endif
};

// a BeamIntegrationRule store BeamIntegration and section tags
class BeamIntegrationRule : public TaggedObject
{
public:
    BeamIntegrationRule(int tag, BeamIntegration* bi, const ID& stags)
	:TaggedObject(tag),theInt(bi),secTags(stags){}

    ~BeamIntegrationRule(){
	if (theInt != 0) 
          delete theInt;
    }

    BeamIntegration* getBeamIntegration(){return theInt;}
    const ID& getSectionTags() const {return secTags;}

    void Print(OPS_Stream &s, int flag) {
      theInt->Print(s);
    }
private:
    BeamIntegration* theInt;
    ID secTags;
};

#endif
