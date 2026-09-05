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
// File: ~/analysis/analysis/Analysis.h
// 
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the interface for the Analysis class.
// Analysis is an abstract class, i.e. no objects of it's type can be created. 
//
// What: "@(#) Analysis.h, revA"
#pragma once

class Domain;

class Analysis
{
  public:
    Analysis(Domain &);
    virtual ~Analysis();

    // pure virtual functions
    virtual int domainChanged() = 0;
    
  protected:
    Domain *getDomainPtr();
    
  private:
    Domain *theDomain;
};

