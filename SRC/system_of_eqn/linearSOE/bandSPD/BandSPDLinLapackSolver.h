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
// $Revision: 1.2 $
// $Date: 2009-05-11 20:55:31 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/bandSPD/BandSPDLinLapackSolver.h,v $
//
// Written: fmk 
// Created: 11/96
//
// Description: This file contains the class definition for 
// BandSPDLinLapackSolver. It solves the BandSPDLinSOE object by calling
// Lapack routines.
//
// What: "@(#) BandSPDLinLapackSolver.h, revA"
//
#pragma once
#include <BandSPDLinSolver.h>

class BandSPDLinLapackSolver : public BandSPDLinSolver
{
public:
  BandSPDLinLapackSolver();    
  ~BandSPDLinLapackSolver();

  int solve() override;
  int solve(const Vector& B, Vector& X) override;
  int setSize() override;
  bool requireDeterminant() override { return false; }
};

