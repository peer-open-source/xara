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
// Description: This file contains the class definition for BFGS.
//
// See also:
// - FEAP: feap/program/iterat.f
//
// [1] 
//
// Written: Ed Love
// Created: 06/01
//
#pragma once
#include <vector>
#include <EquiSolnAlgo.h>
#include <Vector.h>
#include <IncrementalIntegrator.h>
#include <LinearSOE.h> 
#include <LinearAction.h>
class LineSearch;

class BFGS: public EquiSolnAlgo
{
  public:

    BFGS(IncrementalIntegrator::TangentFlagType tangent, int n, LineSearch* search=nullptr);

    ~BFGS();

    int solveCurrentStep() final;

    void Print(OPS_Stream &, int flag) const final;    
    

    
  private:
    int BFGSUpdate(IncrementalIntegrator *,
                    LinearSOE *theSOE,
                    Vector &du, 
                    const Vector &b, 
                    int count);


    struct ApplyBFGS //: public LinearAction
    {
      public:
        ApplyBFGS(int n, BFGS *theAlgo) : 
          rdotz(n+3), sdotr(n+3), temp(0)
        {
        }

        int link(IncrementalIntegrator* integrator, LinearSOE* soe) {
          theIntegrator = integrator;
          theSOE = soe;
          temp.resize(soe->getNumEqn());
          return 0;
        }

        // void reset(int n) {
        //   nBFGS = 0;
        //   this->systemSize = n;
        // }

        // int update(const Vector& Go, const Vector& Gn, const Vector& b) {

        //   // compute z
        //   //  theSOE->setB( (*r[nBFGS]) - (*r[nBFGS-1]) );
        //   //    theSOE->setB( (*residNew) - (*residOld) );
        //   temp.addVector(0.0, Gn, 1.0);
        //   temp.addVector(1.0, Go, -1.0);
        //   theSOE->setB(temp);


        //   if ( z[nBFGS] == nullptr ) 
        //     z[nBFGS] = new Vector(systemSize);


        //   if (theSOE->solve(temp, *z[nBFGS]) < 0)
        //     return -1;

        //   //  *z[nBFGS] *= (-1.0);


        //   for (int i=1; i<=(nBFGS-1); i++ ) {

        //     if ( sdotr[i] < eps )
        //       break; 

        //     double fact1 = 1.0 + ( rdotz[i] / sdotr[i] );

        //     fact1 /= sdotr[i];

        //     double pdotb = (*s[i]) ^ b;

        //     fact1 *= pdotb;

        //     //    *z[nBFGS] +=  fact1 * ( *s[i] );
        //     z[nBFGS]->addVector(1.0, *s[i], fact1);

        //     double bdotz = (*z[i])^b;

        //     //    *z[nBFGS] -= (1.0/sdotr[i]) * 
        //     //             ( bdotz * (*s[i])   +  pdotb * (*z[i]) );   
        //     temp = *s[i];
        //     temp *= bdotz;
        //     temp /= sdotr[i];
        //     *z[nBFGS] -= temp;

        //     temp = *z[i];
        //     temp *= pdotb;
        //     temp /= sdotr[i];
        //     *z[nBFGS] -= temp;
        
        //   } // end for i


        //   //sdotr[nBFGS] = *s[nBFGS] ^ ( *residNew - *residOld );

        //   //rdotz[nBFGS] = *z[nBFGS] ^ ( *residNew - *residOld );   
        //   temp  = Gn;
        //   temp -= Go;

        //   sdotr[nBFGS] = *s[nBFGS] ^ (temp);

        //   rdotz[nBFGS] = *z[nBFGS] ^ (temp);


        //   nBFGS++;
        //   return 0;
        // }

        // int apply(const Vector& x, Vector&b) final {
        //   return -1;
        // }

        // int solve(const Vector &b, Vector &du) final {
        //   // BFGS modifications to du
        //   for (int i=1; i<=nBFGS; i++ ) {

        //     if ( sdotr[i] < eps )
        //       break;

        //     double fact1 = 1.0 + ( rdotz[i] / sdotr[i] );

        //     fact1 /= sdotr[i];

        //     double sdotb = (*s[i]) ^ b;

        //     fact1 *= sdotb;

        //     //du +=  fact1 * ( *s[i] );
        //     du.addVector(1.0, *s[i], fact1);


        //     double bdotz = (*z[i]) ^ b;  

        //     //du -= (1.0/sdotr[i]) * 
        //     //             ( bdotz * (*s[i])   +  sdotb * (*z[i]) );
        //     du.addVector(1.0, *s[i], -bdotz/sdotr[i]);

        //     du.addVector(1.0, *z[i], -sdotb/sdotr[i]);
        //   }
        //   return 0;
        // }
      public:
        IncrementalIntegrator* theIntegrator;
        LinearSOE* theSOE;
        std::vector<double> rdotz;
        std::vector<double> sdotr;
        Vector temp; // temporary vector 

        Vector **s;  // displacement increments
        Vector **z;
    } action;


    IncrementalIntegrator::TangentFlagType tangent;
    int numberLoops;
    LineSearch* search;

    Vector **s;  // displacement increments
    Vector **z;
    Vector Go;  // residuals
    Vector Gn;

    Vector du; // displacement increment
    Vector b;  // current right-hand side
};
