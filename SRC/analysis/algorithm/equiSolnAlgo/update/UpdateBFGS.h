//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, OpenSees/Xara Developers
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
//
#pragma once
#include <LinearAction.h>
#include <Vector.h>
#include <vector>


struct UpdateBFGS : public LinearAction
{
public:
  UpdateBFGS(int n) : 
      rdotz(n+3), sdotr(n+3), bz(0),
      nBFGS(0)
  {
  }


  void reset(int n) {
    nBFGS = 0;
    this->systemSize = n;
  }

  int update(const Vector& Go, const Vector& Gn, LinearSOE& theSOE) {

    // compute z
    //    theSOE->setB( (*residNew) - (*residOld) );
    bz.addVector(0.0, Gn,  1.0);
    bz.addVector(1.0, Go, -1.0);


    if ( z[nBFGS] == nullptr ) 
      z[nBFGS] = new Vector(systemSize);


    if (theSOE.solve(bz, *z[nBFGS]) < 0)
      return -1;

    //  *z[nBFGS] *= (-1.0);


    for (int i=1; i<=(nBFGS-1); i++ ) {

      if ( sdotr[i] < eps )
          break; 

      double fact1 = 1.0 + ( rdotz[i] / sdotr[i] );

      fact1 /= sdotr[i];

      double pdotb = (*s[i]) ^ bz;

      fact1 *= pdotb;

      //    *z[nBFGS] +=  fact1 * ( *s[i] );
      z[nBFGS]->addVector(1.0, *s[i], fact1);

      double bdotz = (*z[i])^bz;

      //    *z[nBFGS] -= (1.0/sdotr[i]) * 
      //             ( bdotz * (*s[i])   +  pdotb * (*z[i]) ); 
      z[nBFGS]->addVector(1.0, *s[i], -bdotz/sdotr[i]);
      z[nBFGS]->addVector(1.0, *z[i], -pdotb/sdotr[i]);

    } // end for i


    // sdotr[nBFGS] = *s[nBFGS] ^ ( *residNew - *residOld );
    // rdotz[nBFGS] = *z[nBFGS] ^ ( *residNew - *residOld );

    sdotr[nBFGS] = (*s[nBFGS] ^ (Go)) - (*s[nBFGS] ^ (Gn));
    rdotz[nBFGS] = (*z[nBFGS] ^ (Go)) - (*z[nBFGS] ^ (Gn));

    nBFGS++;
    return 0;
  }

  int solve(const Vector& x, Vector&b) final {
    return -1;
  }

  int apply(const Vector &b, Vector &du) final {
      // BFGS modifications to du
      for (int i=1; i<=nBFGS; i++ ) {

          if ( sdotr[i] < eps )
              break;

          double fact1 = 1.0 + ( rdotz[i] / sdotr[i] );

          fact1 /= sdotr[i];

          double sdotb = (*s[i]) ^ b;

          fact1 *= sdotb;

          //du +=  fact1 * ( *s[i] );
          du.addVector(1.0, *s[i], fact1);


          double bdotz = (*z[i]) ^ b;  

          //du -= (1.0/sdotr[i]) * 
          //             ( bdotz * (*s[i])   +  sdotb * (*z[i]) );
          du.addVector(1.0, *s[i], -bdotz/sdotr[i]);

          du.addVector(1.0, *z[i], -sdotb/sdotr[i]);
      }
      return 0;
  }
private:
  int nBFGS;
  int systemSize;
  std::vector<double> rdotz;
  std::vector<double> sdotr;
  Vector bz; // temporary vector 

  Vector **s;  // displacement increments
  Vector **z;

  static constexpr double eps = 1.0e-16;
};
