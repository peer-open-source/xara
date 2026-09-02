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
                                                                        
// File: ~/OOP/analysis/algorithm/Linear.h 
// 
// Written: fmk 
// Created: Sun Sept 15 15:06:47: 1996 
// Revision: A 
//

// Description: This file contains the class definition for 
// Linear. Linear is a class which performs a linear solution algorithm
// to solve the equations.
// 
// What: "@(#)Linear.h, revA"
#pragma once
#include <EquiSolnAlgo.h>

class Linear: public EquiSolnAlgo
{
public:
  Linear(int theTangent = CURRENT_TANGENT, int factorOnce = 0);
  ~Linear();

  int solveCurrentStep() final;
  void Print(OPS_Stream &, int flag) const final;    

private:
  int incrTangent;
  int factorOnce;
};
