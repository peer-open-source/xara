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
                                                                        
// Written: Chris McGann
//          February 2011

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <NDMaterial.h>
#include <Vector.h>
#include <Matrix.h>

#include "BoundingCamClay.h"

class BoundingCamClayPlaneStrain : public BoundingCamClay {

  public : 

  //null constructor
  BoundingCamClayPlaneStrain();

  //full constructor
  BoundingCamClayPlaneStrain(int tag, double mDen, double c, double bulk, double OCR, double mu_o,
							          double alpha, double lambda, double h, double m);


  //destructor
  ~BoundingCamClayPlaneStrain();

  NDMaterial* getCopy();
  const char* getType() const;
  int getOrder() const;

  int setTrialStrain(const Vector &strain_from_element);

  // Unused trialStrain functions
  int setTrialStrain(const Vector &v, const Vector &r);
    
  const Vector& getStrain();

  const Vector& getStress();

  const Matrix& getTangent();
  const Matrix& getInitialTangent();

  private :

  // static vectors and matrices
  static Vector strain;
  static Vector stress;
  static Matrix tangent;


};


