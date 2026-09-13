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
                                                                        
/*
Source: /DEVELOPER/element/cpp/Macroelement3d/Macroelement3d/NoTensionSection3d.h
Written by Francesco Vanin (francesco.vanin@epfl.ch)
Ecole Polytechnique Federale de Lausanne, Switzerland,
Earthquake Engineering and Structural Dynamics laboratory, 2019

Reference: Vanin F., Penna A., Beyer K.;"A three dimensional macro-element
for modelling of the in-plane and out-of-plane response of masonry walls",
submitted to Earthquake Engineering and Structural Dynamics (2019)

Last edit: 26 Feb 2019
*/

#ifndef NoTensionSection3d_h
#define NoTensionSection3d_h

#include <SectionForceDeformation.h>
#include <Matrix.h>
#include <Vector.h>

class Channel;
class FEM_ObjectBroker;
class Information;
class Parameter;

class NoTensionSection3d : public SectionForceDeformation
{
 public:
  NoTensionSection3d(int tag, double _k, double _kg, double _t, double _L, double _J, double _fc, int _nSections, 
	                  bool stronger=false, bool elastic=false, bool crushing=true, bool spandrel=false, double r=1.0, double muMax=1e6, bool triangular=false);
  NoTensionSection3d(void);
  ~NoTensionSection3d(void);
  
  const char *getClassType(void) const {return "NoTensionSection3d";}
  
  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);
  
  int setTrialSectionDeformation(const Vector&);
  const Vector &getSectionDeformation(void);
  
  const Vector &getStressResultant();
  const Matrix &getSectionTangent();
  const Matrix &getInitialTangent();
  
  SectionForceDeformation *getCopy();
  const ID &getType();
  int getOrder() const;
  
  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel,
	       FEM_ObjectBroker &theBroker);
  
  void Print(OPS_Stream &s, int flag = 0);
  
  int setParameter(const char **argv, int argc, Parameter &param);
  int updateParameter(int parameterID, Information &info);
  int activateParameter(int parameterID);

  Response *setResponse(const char **argv, int argc, OPS_Stream &s);
  int getResponse(int responseID, Information &info);
  Vector getCrushingVariables(int sliceNum);

 protected:
  
 private:
  
  double k, kg, L, t, J, fc;
  
  int nSections;
  double  IPfactor;        // factor defining the ratio between in-plane and out-of-plane deformation
  Matrix muZ, muY;         // ductility demands for all section discretisations (converged)
  Matrix zetaZ, zetaY;     // damaged portions for all section discretisations (converged)
  Matrix muZt, muYt;       // ductility demands for all section discretisations (trial)
  Matrix zetaZt, zetaYt;   // damaged portions for all section discretisations  (trial)
  Matrix zetaZ2, zetaZ2t;  // crushed portions for all section discretisations  (converged, trial)
  double r, muMax;         // residual strength attained at max ductility
  bool triangular;         // flag for triangular softening response

  Vector pos;              // position of the sections
  Vector weight;           // weight of the sections
  
  Vector eCommitted, e, sCommitted, s, s_zero;     // strain and stress vectors (committed and trial) 
  Matrix D, Dcommitted, Dtrial;            // stiffness matrix (elastic, committed, trial)

  static ID code;

  bool stronger;           // option for adding a small linear elastic contribution for helping convergence
  double factorStronger;   // factor defining the amount of the added elastic contribution, when used
  bool elastic;            // option for a linear elastic formulation
  bool crushing;           // option for accounting or not for crushing in compression
  bool spandrel;           // option for impleting a different behaviour for spandrels (not used now) 
 
  int parameterID;
  int sliceOutput;
  double torsionalStiffnessFactor;  // torsional stiffness modifier, if applied
};

#endif
