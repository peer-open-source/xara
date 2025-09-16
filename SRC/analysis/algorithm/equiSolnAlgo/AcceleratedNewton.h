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
// AcceleratedNewton.  AcceleratedNewton is a class which uses a Krylov
// subspace accelerator on the modified Newton method.
// The accelerator is described by Carlson and Miller in
//
// - "Design and Application of a 1D GWMFE Code"
//   from SIAM Journal of Scientific Computing (Vol. 19, No. 3,
//   pp. 728-765, May 1998)
//
// Written: MHS
// Created: Oct 2001
//
#ifndef AcceleratedNewton_h
#define AcceleratedNewton_h

#include <EquiSolnAlgo.h>
#include <Vector.h>

class Accelerator;

class AcceleratedNewton: public EquiSolnAlgo
{
 public:
  AcceleratedNewton(Accelerator *theAccel, int tangent);
  ~AcceleratedNewton();
  
  int solveCurrentStep(); 
  
  int getNumFactorizations() const override {return numFactorizations;}
  int getNumIterations()  const override {return numIterations;}
  //double getTotalTimeCPU(void)   {return totalTimeCPU;}
  //double getTotalTimeReal(void)  {return totalTimeReal;}
  //double getSolveTimeCPU(void)   {return solveTimeCPU;}
  //double getSolveTimeReal(void)  {return solveTimeReal;}
  //double getAccelTimeCPU(void)   {return accelTimeCPU;}
  //double getAccelTimeReal(void)  {return accelTimeReal;}
  
  virtual int sendSelf(int commitTag, Channel &);
  virtual int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &);
  void Print(OPS_Stream &, int flag) const final;    
  
 protected:
  
 private:
  int tangent;
  
  Accelerator *theAccelerator;
  
  // Storate for accelerated mod-Newton prediction
  Vector *vAccel;
  
  int numFactorizations;
  int numIterations;
};

#endif
