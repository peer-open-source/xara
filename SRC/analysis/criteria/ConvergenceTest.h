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
// Purpose: This file contains the class definition for ConvergenceTest,
// which is an abstract class. Objects of concrete subclasses can be used
// to test the convergence of an algorithm.
//
// Written: fmk
// Date: 09/98
//
#ifndef ConvergenceTest_h
#define ConvergenceTest_h

#define OPS_MAXTOL 1.7e307

#include <MovableObject.h>
#include <StandardStream.h>
#include <Vector.h>
#include <stdbool.h>
#include <string>

class LinearSOE;


class ConvergenceTest: public MovableObject
{
  public:
    enum Status {
      Continue =-1,
      Failure  =-2
    };
    enum Protocol {         // is | was
      Silent        = 0<<0, //  1    0 print nothing
      PrintTest     = 1<<1, //  2    1 print information on norms on test()
      PrintSuccess  = 1<<2, //  4    2 print information on norms and number of iterations at end of successful test
      PrintFailure  = 1<<3, //  .    . 
      PrintTest02   = 1<<4, //  .    4 More verbose test() output
      AlwaysSucceed = 1<<5, //       5 if it fails to converge at end of $numIter it will
                            //         print an error message BUT RETURN A SUCEESSFULL test
      // TODO: add output option 7:
      //       print current iterations dx and du vectors (see commit 9cd8104)
    };

    ConvergenceTest(int classTag);
    virtual ~ConvergenceTest();

    virtual ConvergenceTest *getCopy( int iterations ) = 0 ;

    virtual int start(LinearSOE&) =0;
    virtual int test(LinearSOE&) =0;

    virtual int getNumTests() =0;
    virtual int getMaxNumTests() =0;
    virtual double getRatioNumToMax() =0;
    virtual const Vector &getNorms() =0;


  protected:
    std::string pad(double x);
    std::string pad(int i);
    StandardStream pstream;
};

#endif

