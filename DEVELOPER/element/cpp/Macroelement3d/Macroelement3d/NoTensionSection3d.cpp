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
Source: /DEVELOPER/element/cpp/Macroelement3d/Macroelement3d/NoTensionSection3d.cpp

Written by Francesco Vanin (francesco.vanin@epfl.ch)
Ecole Polytechnique Federale de Lausanne, Switzerland, 
Earthquake Engineering and Structural Dynamics laboratory, 2019

Reference: Vanin F., Penna A., Beyer K.; "A three dimensional macro-element 
for modelling of the in-plane and out-of-plane response of masonry walls", 
submitted to Earthquake Engineering and Structural Dynamics (2019)

Last edit: 26 Feb 2019
*/

/* -------------------------------------------------------------------------- */
#include "NoTensionSection3d.h"
/* -------------------------------------------------------------------------- */
#include <MaterialResponse.h>
#include <classTags.h>
#include <elementAPI.h>
#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <MatrixUtil.h>
#include <Parameter.h>
/* -------------------------------------------------------------------------- */
#include <limits>
#include <cmath>
/* -------------------------------------------------------------------------- */
#ifndef DBL_EPSILON
#define DBL_EPSILON (std::numeric_limits<double>::epsilon())
#endif

#include <stdlib.h>
#include <string.h>

ID NoTensionSection3d::code(4);

NoTensionSection3d::NoTensionSection3d(void) 
	:SectionForceDeformation(0, 0), k(0.0), kg(0.0), t(0.0), L(0.0), J(0.0), fc(0.0), e(4), r(0), muMax(0), triangular(false),
	eCommitted(4), sCommitted(4), s(4), D(4,4), Dtrial(4,4), Dcommitted(4,4), nSections(0), IPfactor(1.0), stronger(false), elastic(false), crushing(true)
{
    if (code(0) != SECTION_RESPONSE_P) {
      code(0) = SECTION_RESPONSE_P;	    // axial force
      code(1) = SECTION_RESPONSE_MZ;	// in-plane moment  
      code(2) = SECTION_RESPONSE_MY;	// out-of-plane moment 
      code(3) = SECTION_RESPONSE_T;	    // Torsion 
    }
    s_zero.Zero();
}

NoTensionSection3d::NoTensionSection3d (int tag, double _k, double _kg, double _t, double _L, double _J, double _fc, int _nSections, bool stronger, bool elastic, bool crushing, bool spandrel, double r, double muMax, bool triangular)
   :SectionForceDeformation(tag, 0),
   k(_k), kg(_kg), t(_t), L(_L), J(_J), fc(_fc), e(4), eCommitted(4), s(4), s_zero(4), sCommitted(4), D(4,4), Dtrial(4,4), Dcommitted(4,4), nSections(_nSections), IPfactor(0.0),
   muZ(_nSections,2),  muY(_nSections,2),  zetaZ(_nSections,2),  zetaY(_nSections,2), pos(_nSections), weight(_nSections), sliceOutput(0),
   muZt(_nSections,2), muYt(_nSections,2), zetaZt(_nSections,2), zetaYt(_nSections,2), zetaZ2(_nSections, 2), zetaZ2t(_nSections, 2),
   stronger(stronger), elastic(elastic), crushing(crushing), spandrel(spandrel), factorStronger(0.001), torsionalStiffnessFactor(1.0), r(r), muMax(muMax), triangular(triangular)
{  
  s_zero.Zero();
	// calculate torsional stiffness if a negative stiffness is provided (assume relatively thin rectangular section)
	if (J<0.0) {
		if (t>L) {
			//J = 1/3.*t*L*L*L;                                            // good for thin sections
			J = t*L*L*L * (1./3. - 0.21*L/t*(1.-std::pow(L,4)/(12*std::pow(t,4))));   // more refined formulation for squatter sections
		}
		else {
			//J = 1/3.*t*t*t*L;
			J = L*t*t*t * (1. / 3. - 0.21*t / L*(1. - std::pow(t, 4) / (12 * std::pow(L, 4))));
		}
	}
	
	// initialise elastic stiffness matrix 
	D(0,0) = k*L*t;
    D(1,1) = k*t*L*L*L /12.0;
    D(2,2) = k*t*t*t*L /12.0;

	D(3,3) = kg*J;

	// change torsional stiffness is a modifier is introduced
	D(3,3)*= torsionalStiffnessFactor;

	// sum a (small) elastic contribution if the section is defined as "stronger"
	// do nothing is the section is already defined as "elastic"
	if (stronger && !elastic) {
		D(0,0) += factorStronger*k*L*t;
        D(1,1) += factorStronger*k*t*L*L*L /12.0;
        D(2,2) += factorStronger*k*t*t*t*L /12.0;
		D(3,3) += factorStronger*kg*J;
	} 

	Dtrial = D;       // trial stiffness matrix (updated at each displacement update)
	Dcommitted= D;    // committed stiffness matrix (updated at each converged step)

	s.Zero();
	sCommitted.Zero();

	// initialise response identifiers (standard ordering)
    if (code(0) != SECTION_RESPONSE_P) {
	    code(0) = SECTION_RESPONSE_P;	// axial force
	    code(1) = SECTION_RESPONSE_MZ;	// in-plane moment  
	    code(2) = SECTION_RESPONSE_MY;	// out-of-plane moment 
	    code(3) = SECTION_RESPONSE_T;	// Torsion 
    }

	// initialise slice weights for nonlinearity in compression
    if (nSections==1) {
	    pos(0) = 0.0;
	    weight(0) = 1;
    } else {
	   double slice = 1./(nSections);
	   for (int i=0; i<nSections; i++) {
		    pos(i) = -0.5 + slice/2. + slice * i;
		    weight(i) = slice;
	   }
    }

	// inititalise ductility demand parameters to 1 
    for (int i=0; i<nSections; i++) {
	  muZ(i,0) = 1.;
	  muZ(i,1) = 1.;
	  muY(i,0) = 1.;
	  muY(i,1) = 1.;
    }

}

NoTensionSection3d::~NoTensionSection3d(void)
{
    return;
}

int 
NoTensionSection3d::commitState(void)
{
	eCommitted = e;
	sCommitted = s;

	Dcommitted = Dtrial;

    muZ = muZt;
	muY = muYt;
	zetaZ = zetaZt;
	zetaZ2 = zetaZ2t;
	zetaY = zetaYt;

  return 0;
}

int 
NoTensionSection3d::revertToLastCommit(void)
{
	e = eCommitted;
	s = sCommitted;
	Dtrial = Dcommitted;

	muZt = muZ;
	muYt = muY;
	zetaZt = zetaZ;
	zetaZ2t = zetaZ2;
	zetaYt = zetaY;


  return 0;
}

int 
NoTensionSection3d::revertToStart(void)
{
	muZ.Zero();
	muY.Zero();
	zetaZ.Zero();
	zetaZ2.Zero();
	zetaY.Zero();

	e.Zero();
	s.Zero();
	eCommitted.Zero();
	sCommitted.Zero();
	Dtrial = D;
	Dcommitted = D;

  return 0;
}

int
NoTensionSection3d::setTrialSectionDeformation (const Vector &def)
{
	Vector increment(4);
	increment = e - def;

	if (std::abs(increment.Norm()) < DBL_EPSILON) {
      return 0;
	}

    e = def;
	increment = e - eCommitted;

	if (elastic) {  // return a linear elastic solution
		s = D*e;
		Dtrial = D;
		return 0;
	}
	
	// initialise trial damage variables to last converged value
	muZt = muZ;
	muYt = muY;
	zetaZt = zetaZ;
	zetaYt = zetaY;


	//---------------------------------------------------------------------------------------------------------//
	// set stress
	//---------------------------------------------------------------------------------------------------------//

	s.Zero();
	Dtrial.Zero();

	if (std::abs(e(2)) < DBL_EPSILON) {
		// zero in plane moment
		if (std::abs(e(1)) + 2.0*e(0)/L > DBL_EPSILON) {
			if  (-std::abs(e(1)) + 2.0*e(0)/L > DBL_EPSILON) {
				// all tension, all zeroes, do nothing
			} else {				
					// case 5
					double eps0 = e(0);    // axial deformation at the origin
					double chiZ = e(1);    // curvature around axis z
					double chiY = 0.0;     // curvature around axis y

					if (chiZ<0.0) {
						// case 4a,c: both out from different sides, chiZ>=0
							 s(0) = (std::pow(chiY,2)*std::pow(t,3) + 3*t*std::pow(chiZ*L + 2*eps0,2))/(24.*chiZ);
							 s(1) = -(std::pow(chiY,2)*std::pow(t,3)*eps0 + t*(-(chiZ*L) + eps0)*std::pow(chiZ*L + 2*eps0,2))/(24.*std::pow(chiZ,2));
							 s(2) = (chiY*std::pow(t,3)*(chiZ*L + 2*eps0))/(24.*chiZ);

					} else {
						// case 4b,d: both out from different sides, chiZ<0
							 s(0) = -(std::pow(chiY,2)*std::pow(t,3) + 3*t*std::pow(chiZ*L - 2*eps0,2))/(24.*chiZ);
							 s(1) = (std::pow(chiY,2)*std::pow(t,3)*eps0 + t*std::pow(chiZ*L - 2*eps0,2)*(chiZ*L + eps0))/(24.*std::pow(chiZ,2));
							 s(2) = (chiY*std::pow(t,3)*(chiZ*L - 2*eps0))/(24.*chiZ);
					}
			}

		} else {
			// all compression, elastic solution
			s(0) = L*t *e(0);
            s(1) = t*L*L*L /12.0 * e(1);
            s(2) = t*t*t*L /12.0 * e(2);
		}


	} else {
		double eps0 = e(0);    // axial deformation at the origin
		double chiZ = e(1);    // curvature around axis z
		double chiY = e(2);    // curvature around axis y

		double t1 = (2.0*eps0 + L*chiZ)/(2.0*chiY);   // x coordinate of the neutral axis for y = +L/2
		double t2 = (2.0*eps0 - L*chiZ)/(2.0*chiY);   // x coordinate of the neutral axis for y = -L/2

		int inside1 = 0;
		if (t1<-t/2.0)   inside1=-1;
		if (t1> t/2.0)   inside1=+1;

	    int inside2 = 0;
		if (t2<-t/2.0)   inside2=-1;
		if (t2> t/2.0)   inside2=+1;
		
		if (inside1==0 && inside2==0) {
			if (chiY>0.0) {
			   // case 1a,b: both neutral points inside the long side, fi_y>0
				 s(0) = -(L*(std::pow(chiZ,2)*std::pow(L,2) + 3*std::pow(chiY*t - 2*eps0,2)))/(24.*chiY);
				 s(1) = (chiZ*std::pow(L,3)*(chiY*t - 2*eps0))/(24.*chiY);
				 s(2) = (L*(std::pow(chiY,3)*std::pow(t,3) - 3*std::pow(chiY,2)*std::pow(t,2)*eps0 + std::pow(chiZ,2)*std::pow(L,2)*eps0 + 4*std::pow(eps0,3)))/(24.*std::pow(chiY,2));

			} else {
               // case 1c,d: both neutral points inside the long side, fi_y<0
				 s(0) = (L*(std::pow(chiZ,2)*std::pow(L,2) + 3*std::pow(chiY*t + 2*eps0,2)))/(24.*chiY);
				 s(1) = (chiZ*std::pow(L,3)*(chiY*t + 2*eps0))/(24.*chiY);
				 s(2) = (L*(std::pow(chiY,3)*std::pow(t,3) + 3*std::pow(chiY,2)*std::pow(t,2)*eps0 - std::pow(chiZ,2)*std::pow(L,2)*eps0 - 4*std::pow(eps0,3)))/(24.*std::pow(chiY,2));
			}

		} else {
			if (inside1*inside2>0) {
				if (eps0 - t/2.0*chiY + L/2.0*chiZ < DBL_EPSILON) {
				    // case 0a: both neutral points ouside the same side, all compressed
					s(0) = L*t *e(0);
                    s(1) = t*L*L*L /12.0 * e(1);
                    s(2) = t*t*t*L /12.0 * e(2);
				} else {
					// case 0b: both neutral points ouside the same side, all tension
				    s(0) = 0.0;
                    s(1) = 0.0;
                    s(2) = 0.0;
				}
			} else {
				if (inside1*inside2==0) {
					// one in and one out. Cases 2-3
					if (inside2>0) {
						if (chiY>0.0) {
							// case 2a
							 s(0) = -std::pow(-(chiY*t) + chiZ*L + 2*eps0,3)/(48.*chiZ*chiY);
							 s(1) = -((chiY*t + 3*chiZ*L - 2*eps0)*std::pow(-(chiY*t) + chiZ*L + 2*eps0,3))/(384.*std::pow(chiZ,2)*chiY);
							 s(2) = (std::pow(-(chiY*t) + chiZ*L + 2*eps0,3)*(3*chiY*t + chiZ*L + 2*eps0))/(384.*chiZ*std::pow(chiY,2));

						} else {
							// case 3a
						     s(0) = t*L*eps0 + std::pow(-(chiY*t) + chiZ*L + 2*eps0,3)/(48.*chiZ*chiY);
							 s(1) = (chiZ*t*std::pow(L,3))/12. + ((chiY*t + 3*chiZ*L - 2*eps0)*std::pow(-(chiY*t) + chiZ*L + 2*eps0,3))/(384.*std::pow(chiZ,2)*chiY);
							 s(2) = (chiY*std::pow(t,3)*L)/12. - (std::pow(-(chiY*t) + chiZ*L + 2*eps0,3)*(3*chiY*t + chiZ*L + 2*eps0))/(384.*chiZ*std::pow(chiY,2));
						}

					} else {
						if (inside2<0) {
							if (chiY<0.0) {
								// case 2c
								 s(0) = std::pow(chiY*t + chiZ*L + 2*eps0,3)/(48.*chiZ*chiY);
								 s(1) = -((chiY*t - 3*chiZ*L + 2*eps0)*std::pow(chiY*t + chiZ*L + 2*eps0,3))/(384.*std::pow(chiZ,2)*chiY);
								 s(2) = -((-3*chiY*t + chiZ*L + 2*eps0)*std::pow(chiY*t + chiZ*L + 2*eps0,3))/(384.*chiZ*std::pow(chiY,2));;

							} else {
								// case 3c
									 s(0) = t*L*eps0 - std::pow(chiY*t + chiZ*L + 2*eps0,3)/(48.*chiZ*chiY);
									 s(1) = (chiZ*t*std::pow(L,3))/12. + ((chiY*t - 3*chiZ*L + 2*eps0)*std::pow(chiY*t + chiZ*L + 2*eps0,3))/(384.*std::pow(chiZ,2)*chiY);
									 s(2) = (chiY*std::pow(t,3)*L)/12. + ((-3*chiY*t + chiZ*L + 2*eps0)*std::pow(chiY*t + chiZ*L + 2*eps0,3))/(384.*chiZ*std::pow(chiY,2));

							}

						} else {
							if (inside1>0) {
								if (chiY>0.0) {
									// case 2b
								    s(0) = -std::pow(chiY*t + chiZ*L - 2*eps0,3)/(48.*chiZ*chiY);
									s(1) = (std::pow(chiY*t + chiZ*L - 2*eps0,3)*(-(chiY*t) + 3*chiZ*L + 2*eps0))/(384.*std::pow(chiZ,2)*chiY);
									s(2) = (std::pow(chiY*t + chiZ*L - 2*eps0,3)*(3*chiY*t - chiZ*L + 2*eps0))/(384.*chiZ*std::pow(chiY,2));

								} else {
									// case 3b
								     s(0) = std::pow(chiY*t + chiZ*L - 2*eps0,3)/(48.*chiZ*chiY) + t*L*eps0;
									 s(1) = (chiZ*t*std::pow(L,3))/12. + ((chiY*t - 3*chiZ*L - 2*eps0)*std::pow(chiY*t + chiZ*L - 2*eps0,3))/(384.*std::pow(chiZ,2)*chiY);
									 s(2) = (chiY*std::pow(t,3)*L)/12. + ((-3*chiY*t + chiZ*L - 2*eps0)*std::pow(chiY*t + chiZ*L - 2*eps0,3))/(384.*chiZ*std::pow(chiY,2));
								}

							} else {
								if (chiY<0.0) {
									// case 2d
									 s(0) = std::pow(-(chiY*t) + chiZ*L - 2*eps0,3)/(48.*chiZ*chiY);
									 s(1) = (std::pow(chiY*t - chiZ*L + 2*eps0,3)*(chiY*t + 3*chiZ*L + 2*eps0))/(384.*std::pow(chiZ,2)*chiY);
									 s(2) = (std::pow(-(chiY*t) + chiZ*L - 2*eps0,3)*(3*chiY*t + chiZ*L - 2*eps0))/(384.*chiZ*std::pow(chiY,2));
								} else {
									// case 3d
							         s(0) = t*L*eps0 + std::pow(chiY*t - chiZ*L + 2*eps0,3)/(48.*chiZ*chiY);
							         s(1) = (chiZ*t*std::pow(L,3))/12. + (std::pow(-(chiY*t) + chiZ*L - 2*eps0,3)*(chiY*t + 3*chiZ*L + 2*eps0))/(384.*std::pow(chiZ,2)*chiY);
							         s(2) = (chiY*std::pow(t,3)*L)/12. + ((3*chiY*t + chiZ*L - 2*eps0)*std::pow(chiY*t - chiZ*L + 2*eps0,3))/(384.*chiZ*std::pow(chiY,2));
								}
							}
						}
					}

				} else {
					if (chiZ<0.0) {
						// case 4a,c: both out from different sides, chiZ>=0
							 s(0) = (std::pow(chiY,2)*std::pow(t,3) + 3*t*std::pow(chiZ*L + 2*eps0,2))/(24.*chiZ);
							 s(1) = -(std::pow(chiY,2)*std::pow(t,3)*eps0 + t*(-(chiZ*L) + eps0)*std::pow(chiZ*L + 2*eps0,2))/(24.*std::pow(chiZ,2));
							 s(2) = (chiY*std::pow(t,3)*(chiZ*L + 2*eps0))/(24.*chiZ);

					} else {
						// case 4b,d: both out from different sides, chiZ<0							 
							 s(0) = -(std::pow(chiY,2)*std::pow(t,3) + 3*t*std::pow(chiZ*L - 2*eps0,2))/(24.*chiZ);
							 s(1) = (std::pow(chiY,2)*std::pow(t,3)*eps0 + t*std::pow(chiZ*L - 2*eps0,2)*(chiZ*L + eps0))/(24.*std::pow(chiZ,2));
							 s(2) = (chiY*std::pow(t,3)*(chiZ*L - 2*eps0))/(24.*chiZ);
					}
				}
			}
		}

		}

	    s *= k;

		// add elastic torsional response
	    s(3) = kg*J*e(3);

		// apply stiffness modifier (default: 1)
		s(3)*= torsionalStiffnessFactor;

		// ----------------------------------------------------------------------------------------------------------------------
		// CRUSHING CORRECTION
		// ----------------------------------------------------------------------------------------------------------------------

		// slices along direction z (higher accuracy in the in-plane direction; 
		// for out-of-plane loading works like a fibre discretisation) 
		//---------------------------------------------------------------------------------------------------------//
		double muNew;
		double zetaNew;
		double zeta2New;

		double muNew_p;
		double zetaNew_p;
		double zeta2New_p;
		double muNew_n;
		double zetaNew_n;
		double zeta2New_n;

		double f1_p, f2_p, fy_p;
		double f1_n, f2_n, fy_n;
		double d_ny, d_n2, d_n;
		double d_py, d_p2, d_p;
		double corr;
		double OOPfactor;

		Vector bi(4);
		Vector c(4);

		Vector dmu0_de(4);
		Vector dmu1_de(4);
		Vector dzeta0_de(4);
		Vector dzeta1_de(4);
		Vector dCorr_de(4);
		


		// enter only if there is some stress (anf therefore the section is not open in tension)
		if ((s.Norm() > DBL_EPSILON) && crushing) {
			/*
		    if ((std::sqrt(std::pow(e(1),2) + std::pow(e(2),2)) ) > DBL_EPSILON )  // if there is some rotation define the proportion of IP/OOP rotation
     		  if (std::abs(e(1))*L/2. > 1./100. *(fc/k) || std::abs(e(2))*t/2. > 1./100. *(fc/k) ) 
			  	IPfactor = std::abs(e(1)) / std::sqrt(std::pow(e(1),2.) + std::pow(e(2),2.));     
		    */ 

			IPfactor = 1.0;  // commented previous lines: applies only slice discretisation parallel to the dimension L 

            OOPfactor = 1.0-IPfactor;  // not used at the moment

			// initialise damage variables
			muZt = muZ;
			zetaZt = zetaZ;
			zetaZ2t = zetaZ2;

			// initialise constant values of vectors b_i and c
			bi(0) = 1.0;
			bi(3) = 0.0;
			c.Zero();

			for (int i=0; i<nSections; i++) {
				
				// zero derivatives
				dmu0_de.Zero();
				dmu1_de.Zero();
				dzeta0_de.Zero();
				dzeta1_de.Zero();
				dCorr_de.Zero();
				double x_p = 0.0;
				double x_n = 0.0;
					
				
					
					// check negative side
					// --------------------------------------------------------------------------------------------------
					// define vector b_i for the i-th slice (1, -y, z)
					bi(1) = L / 2.0;
					bi(2) = t*pos(i);

					// define vector c
					c(1) = -L / 3.;

					// trial values for damage variables
					muNew_n = -k / fc *(bi^e);
					if (e(1) < 0) {
						zetaNew_n = (muNew_n - 1.0)*fc / k / (3.0*c^e);
						zeta2New_n = (muNew_n - muMax)*fc / k / (3.0*c^e);
					}
					else {
						zetaNew_n = 0;
						zeta2New_n = 0;
					}

					// updates
					if (muNew_n > muZt(i, 0)) {
						muZt(i, 0) = muNew_n;
					}

					if (zetaNew_n > zetaZt(i, 0)) {
						if (zetaNew_n <= 1.0) {
							zetaZt(i, 0) = zetaNew_n;
							//dzeta0_de = -1.0*((bi*(3.0*c^e) - 3.0*c*(bi^e + fc / k)) / pow(3.0*c^e, 2.0));
						}
						else {
							zetaZt(i, 0) = 1.0;
							zetaZt(i, 1) = 1.0;
						}
					}

					if (zeta2New_n > zetaZ2t(i, 0)) {
						if (zeta2New_n <= 1.0) {
							zetaZ2t(i, 0) = zeta2New_n;
						}

						else {
							zetaZ2t(i, 0) = 1.0;
							zetaZ2t(i, 1) = 1.0;
						}
					}


					d_ny = e(0) + L*(1. / 2 - zetaZt(i, 0))  * e(1) + t*pos(i)*e(2);
					d_n2 = e(0) + L*(1. / 2 - zetaZ2t(i, 0)) * e(1) + t*pos(i)*e(2);
					d_n = e(0) + L / 2.*e(1) + t*pos(i)*e(2);

					if (zetaZ2t(i, 0) > 0 && d_n2 < 0) {
						fy_n = -k*(d_ny);
						f1_n = r*fc* d_n2 / (-fc / k) / muMax;
						if (zetaZ2t(i, 0) > 0.0) {
							f2_n = fy_n + (f1_n - fy_n) / (zetaZt(i, 0) - zetaZ2t(i, 0))*zetaZt(i, 0);
						}
						if (d_n < 0.0) {
							f1_n = r*fc* d_n / (-fc / k) / muZt(i, 0);
							x_n = 0.0;
						}
						else {
							x_n = d_n / e(1);
							f1_n = -f1_n / (zetaZ2t(i, 0)*L - x_n)*x_n;
						}

					}
					else {
						x_n = f1_n = f2_n = fy_n = 0;
					}



					// check positive side
					// --------------------------------------------------------------------------------------------------
					// define vector b_i for the i-th slice (1, -y, z)
					bi(1) = -L / 2.0;
					c(1) = L / 3.;

					// trial values for damage variables
					muNew_p = -k / fc *(bi^e);
					if (e(1)>0) {
						zetaNew_p = (muNew_p - 1.0)*fc / k / (3.0*c^e);
						zeta2New_p = (muNew_p - muMax)*fc / k / (3.0*c^e);
					}
					else {
						zetaNew_p = 0;
						zeta2New_p = 0;
					}

					// updates
					if (muNew_p > muZt(i, 1)) {
						muZt(i, 1) = muNew_p;
					}

					if (zetaNew_p > zetaZt(i, 1)) {
						if (zetaNew_p <= 1.0) {
							zetaZt(i, 1) = zetaNew_p;
							//dzeta0_de = -1.0*((bi*(3.0*c^e) - 3.0*c*(bi^e + fc / k)) / pow(3.0*c^e, 2.0));
						}
						else {
							zetaZt(i, 0) = 1.0;
							zetaZt(i, 1) = 1.0;
						}
					}

					if (zeta2New_p > zetaZ2t(i, 1)) {
						if (zeta2New_p <= 1.0) {
							zetaZ2t(i, 1) = zeta2New_p;
						}

						else {
							zetaZ2t(i, 0) = 1.0;
							zetaZ2t(i, 1) = 1.0;
						}
					}


					d_py = e(0) - L*(1. / 2 - zetaZt(i, 1))  * e(1) + t*pos(i)*e(2);
					d_p2 = e(0) - L*(1. / 2 - zetaZ2t(i, 1)) * e(1) + t*pos(i)*e(2);
					d_p  = e(0) - L / 2.*e(1) + t*pos(i)*e(2);

					if (zetaZ2t(i, 1) > 0 && d_p2 < 0) {
						fy_p = -k*(d_py);
						f1_p = r*fc* d_p2 / (-fc / k) / muMax;
						if (zetaZ2t(i, 1) > 0.0) {
							f2_p = fy_p + (f1_p - fy_p) / (zetaZt(i, 1) - zetaZ2t(i, 1))*zetaZt(i, 1);
						}
						if (d_p < 0.0) {
							f1_p = r*fc* d_p / (-fc / k) / muZt(i, 1);
							x_p = 0.0;
						}
						else {
							x_p = -d_p / e(1);
							f1_p = -f1_p / (zetaZ2t(i, 1)*L - x_p)*x_p;
						}

					}
					else {
						x_p = f1_p = f2_p = fy_p = 0;
					}

				

					// apply correction 
					// -------------------------------------------------------------------------

					if (muZt(i, 0)>1.0) {     // negative side correction term

						bi(1) = L / 2.0;
						bi(2) = t*pos(i);
						c(1) = -L / 3.;
						double beta = (muMax - r) / (muMax - 1);

						if (d_n < 0.0 && muZt(i, 0) < muMax) {  // apply correction as long as the the edge is compressed
							corr = -k*beta*(muZt(i, 0) - 1.0) / (2.0*muZt(i, 0))  *t*weight(i) * (zetaZt(i, 0)*L) * (bi^e);
							s += (bi + zetaZt(i, 0)*c)*corr;

							dCorr_de.Zero();
							dCorr_de += -k*(muZt(i, 0) - 1.0) / (2.0*muZt(i, 0))  *t*weight(i) * (zetaZt(i, 0)*L) * (bi);
							dCorr_de += corr / (muZt(i, 0) - 1.0) / (2.0*muZt(i, 0))* pow(muZt(i, 0), -2.0) * dmu0_de;
							dCorr_de += (corr / zetaZt(i, 0))  * dzeta0_de;

							//	Dtrial += (bi + zetaZt(i, 0)*c) % (dCorr_de);
							//	Dtrial += corr* (c % dzeta0_de);
						}

						if (d_n2 < 0.0 && muZt(i, 0) > muMax) {
							corr = 0.5*(-k* (bi^e) - f2_n) *t*weight(i) * (zetaZt(i, 0)*L);
							s += (bi + zetaZt(i, 0)*c)*corr;

							c(1) = -L / 3.;
							corr = 0.5*(f1_n - f2_n)  *t*weight(i) * zetaZ2t(i, 0)*L;
							s -= (bi + zetaZ2t(i, 0)*c)*corr;

							if (d_n>0) {
								c(1) = -x_n / 3.;
								corr = 0.5*(f1_n - (-k* (bi^e)))  *t*weight(i) * x_n;
								s += (bi + c)*corr;
							}
						}

					}

					if (muZt(i, 1) > 1.0) {     // positive side correction term

						bi(1) = -L / 2.0;
						bi(2) = t*pos(i);
						c(1) = L / 3.;
						double beta = (muMax - r) / (muMax - 1);

						if (d_p < 0.0 && muZt(i, 1) <= muMax) {  // apply correction as long as the the edge is compressed

							//corr = -k*(muZt(i, 1) - 1.0) / (2.0*muZt(i, 1))  *t*weight(i) * (zetaZt(i, 1)*L) * (bi^e);
							//s += (bi + zetaZt(i, 1)*c)*corr;

							dCorr_de.Zero();
							dCorr_de += -k*(muZt(i, 1) - 1.0) / (2.0*muZt(i, 1))  *t*weight(i) * (zetaZt(i, 1)*L) * (bi);
							dCorr_de += corr / (muZt(i, 1) - 1.0) / (2.0*muZt(i, 1))* pow(muZt(i, 1), -2.0) * dmu1_de;
							dCorr_de += corr / zetaZt(i, 1)   * dzeta1_de;

							//	Dtrial += (bi + zetaZt(i, 1)*c) % (dCorr_de);
							//	Dtrial += corr* (c % dzeta1_de);
							corr = -k*beta*(muZt(i, 1) - 1.0) / (2.0*muZt(i, 1))  *t*weight(i) * (zetaZt(i, 1)*L) * (bi^e);
							s += (bi + zetaZt(i, 1)*c)*corr;
						}

						if (d_p2 < 0.0 && muZt(i, 1) > muMax) {
							bi(1) = -L / 2.0;
							corr = 0.5*(-k* (bi^e) - f2_p) *t*weight(i) * (zetaZt(i, 1)*L);
							s += (bi + zetaZt(i, 1)*c)*corr;

							c(1) = L / 3.;
							corr = 0.5*(f1_p - f2_p)  *t*weight(i) * zetaZ2t(i, 1)*L;
							s -= (bi + zetaZ2t(i, 1)*c)*corr;

							if (d_p > 0) {
								c(1) = x_p / 3.;
								corr = 0.5*(f1_p - (-k* (bi^e)))  *t*weight(i) * x_p;
								s += (bi + c)*corr;
							}
						}
					}

				
				
				
				
				
				
				
				/*
				if (e(1) < -DBL_EPSILON) {  // negative side (more) compressed. 

					// define vector b_i for the i-th slice (1, -y, z)
					bi(1) = L / 2.0;
					bi(2) = t*pos(i);

					// define vector c
					c(1) = -L / 3.;


					// trial values for damage variables
					muNew = -k/fc *(bi^e);
					zetaNew = (muNew-1.0)*fc/k / (3.0*c^e);
					zeta2New = (muNew - muMax)*fc / k / (3.0*c^e);
					

					if (muNew > muZt(i, 0)) {
						muZt(i, 0) = muNew;
						dmu0_de = -k / fc *bi;
					}

					if (zetaNew > zetaZt(i, 0)) {
						if (zetaNew <= 1.0) {
							zetaZt(i, 0) = zetaNew;
							dzeta0_de = -1.0*((bi*(3.0*c^e) - 3.0*c*(bi^e + fc/k)) / std::pow(3.0*c^e, 2.0));
						}

						else {
							zetaZt(i, 0) = 1.0;
							zetaZt(i, 1) = 1.0;

							// update ductility demand at the other side
							double muOtherSide = muNew *(zetaNew - 1.0) / zetaNew;
							if (muOtherSide > muZt(i, 1)) {
								muZt(i, 1) = muOtherSide;
								dmu1_de = ((-1. + zetaNew)) / zetaNew  * dmu0_de;
							}
						}
					}

					if (zeta2New > zetaZ2t(i, 0)) {
						if (zeta2New <= 1.0) {
							zetaZ2t(i, 0) = zeta2New;
						}

						else {
							zetaZ2t(i, 0) = 1.0;
							zetaZ2t(i, 1) = 1.0;
						}
					}

					fy_n = -k*((bi^e) - (e(1))*zetaZt(i, 0) *L);
					f1_n = r*fc * ((bi^e) - (e(1))*zetaZ2t(i, 0) *L)/(-fc/k) / muMax;
					f2_n = fy_n + (f1_n - fy_n) / (zetaZt(i, 0) - zetaZ2t(i, 0))*zetaZt(i, 0);
					f1_n = r*fc *muNew / muZt(i, 0);


					double d_py = e(0) - L*(1. / 2 - zetaZt(i, 1))  * e(1) + t*pos(i)*e(2);
					double d_p2 = e(0) - L*(1. / 2 - zetaZ2t(i, 1)) * e(1) + t*pos(i)*e(2);
					double d_p  = e(0) - L / 2.*e(1) + t*pos(i)*e(2);

					if (d_p < 0.0 && zetaZ2t(i, 1)>0) {
						fy_p = -k*(d_py);
						f1_p = r*fc* d_p2 / (-fc / k) / muMax;
						if (zetaZ2t(i, 1) > 0.0) {
							f2_p = fy_p + (f1_p - fy_p) / (zetaZt(i, 1) - zetaZ2t(i, 1))*zetaZt(i, 1);
						}
						f1_p = r*fc* d_p / (-fc / k) / muZt(i, 1);

					}
					else {
						if (d_p2 < 0 && zetaZ2t(i, 1)>0) {
							fy_p = -k*(d_py);
							f1_p = r*fc* d_p2 / (-fc / k) / muMax;
							if (zetaZ2t(i, 1) > 0.0) {
								f2_p = fy_p + (f1_p - fy_p) / (zetaZt(i, 1) - zetaZ2t(i, 1))*zetaZt(i, 1);
							}
							x_p = -d_p / e(1);
							f1_p = -f1_p / (zetaZ2t(i, 1)*L - x_p)*x_p;

						}
						else {
							f2_p = 0.0;
							f1_p = 0.0;
						}
					}
	
				} else  {               

					// check  positive side  displ[z_] := wLoc - fiLoc*z; 

					// define vector b_i for the i-th slice
					bi(1) = -L / 2.0;
					bi(2) = t*pos(i);

					// define vector c
					c(1) = L / 3.;

					// trial values for damage variables
					muNew = -k/fc *bi^e;
					zetaNew = (muNew - 1.0)*fc / k / (3.0*c^e);
					zeta2New = (muNew - muMax)*fc / k / (3.0*c^e);

					if (muNew > muZt(i,1))  {
						muZt(i,1) = muNew;
						dmu1_de = -k / fc *bi;
					}

					if (zetaNew > zetaZt(i,1)) {
						if (zetaNew<=1.0) {
							zetaZt(i,1) = zetaNew;
							dzeta1_de = -1.0*((bi*(3.0*c^e) - 3.0*c*(bi^e + fc/k))/ std::pow(3.0*c^e, 2.0));

						} else {						
							zetaZt(i,0) = 1.0;
							zetaZt(i,1) = 1.0;
							
							// update ductility demand at the other side
							double muOtherSide = muNew *(zetaNew-1.0)/zetaNew;

							if (muOtherSide > muZt(i,0)) {
								muZt(i,0) = muOtherSide;
								dmu0_de = ((-1. + zetaNew)) / zetaNew  * dmu1_de;
							}				
						}
					}


					if (zeta2New > zetaZ2t(i, 1)) {
						if (zeta2New <= 1.0) {
							zetaZ2t(i, 1) = zeta2New;
						}
						else {
							zetaZ2t(i, 0) = 1.0;
							zetaZ2t(i, 1) = 1.0;
						}
					}
					
					fy_p = -k*(e(0) - L*(1. / 2 - zetaZt(i, 1))  * e(1) + t*pos(i)*e(2));
					f1_p = r*fc * (e(0) - L*(1. / 2 - zetaZ2t(i, 1))  * e(1) + t*pos(i)*e(2)) / (-fc / k) / muMax;
					f2_p = fy_p + (f1_p - fy_p) / (zetaZt(i, 1) - zetaZ2t(i, 1))*zetaZt(i, 1);
					f1_p = r*fc *muNew / muZt(i, 1);
					

					double d_ny = e(0) + L*(1. / 2 - zetaZt(i, 0))  * e(1) + t*pos(i)*e(2);
					double d_n2 = e(0) + L*(1. / 2 - zetaZ2t(i, 0)) * e(1) + t*pos(i)*e(2);
					double d_n  = e(0) + L / 2.*e(1) + t*pos(i)*e(2);
					
					if (d_n < 0.0 && zetaZ2t(i, 0)>0) {
						fy_n = -k*(d_ny);
						f1_n = r*fc* d_n2 / (-fc / k) / muMax;
						if (zetaZ2t(i, 0) > 0.0) {
							f2_n =  fy_n + (f1_n - fy_n) / (zetaZt(i, 0) - zetaZ2t(i, 0))*zetaZt(i, 0);
						}
						f1_n = r*fc* d_n / (-fc / k) / muZt(i, 0);

					} else {
						if (d_n2 < 0 && zetaZ2t(i, 0)>0) {
							fy_n = -k*(d_ny);
							f1_n = r*fc* d_n2 / (-fc / k) / muMax;
							if (zetaZ2t(i, 0) > 0.0) {
								f2_n = fy_n + (f1_n - fy_n) / (zetaZt(i, 0) - zetaZ2t(i, 0))*zetaZt(i, 0);
							}
							x_n = d_n/e(1);
							f1_n = -f1_n/(zetaZ2t(i, 0)*L-x_n)*x_n;
							
						} else {
							f2_n = 0.0;
							f1_n = 0.0;
						}
					}
					
					


				} else {              // the section curvature along axis z is basically null
					bi(1) = 0;        // any value would be fine
					bi(2) = t*pos(i);
					
					muNew = -k / fc *bi^e;

					// zeta is (in theory) infinite

					if (muNew > muZt(i,0))  {
						muZt(i,0) = muNew;
						dmu0_de = -k / fc *bi;
					}

					if (muNew > muZt(i,1))  {
						muZt(i,1) = muNew;
						dmu1_de = dmu0_de;
					}

					if (zetaZt(i,0) != 1.0 && muNew>1.) {  // was not already all crushed and there is some crushing now
						zetaZt(i,0) = 1.0;
					}
			        if (zetaZt(i,1) != 1.0 && muNew>1.) {  // was not already all crushed and there is some crushing now
						zetaZt(i,1) = 1.0;
					}

					if (zetaZ2t(i, 0) != 1.0 && muNew>muMax) {  // was not already all crushed and there is some crushing now
						zetaZ2t(i, 0) = 1.0;
					}
					if (zetaZ2t(i, 1) != 1.0 && muNew>muMax) {  // was not already all crushed and there is some crushing now
						zetaZ2t(i, 1) = 1.0;
					}

					double d_ny = e(0) + L*(1. / 2 - zetaZt(i, 0))  * e(1) + t*pos(i)*e(2);
					double d_n2 = e(0) + L*(1. / 2 - zetaZ2t(i, 0)) * e(1) + t*pos(i)*e(2);
					double d_n = e(0) + L / 2.*e(1) + t*pos(i)*e(2);

					if (d_n < 0.0 && zetaZ2t(i, 0)>0) {
						fy_n = -k*(d_ny);
						f1_n = r*fc* d_n2 / (-fc / k) / muMax;
						if (zetaZ2t(i, 0) > 0.0) {
							f2_n = fy_n + (f1_n - fy_n) / (zetaZt(i, 0) - zetaZ2t(i, 0))*zetaZt(i, 0);
						}
						f1_n = r*fc* d_n / (-fc / k) / muZt(i, 0);
					}

					double d_py = e(0) - L*(1. / 2 - zetaZt(i, 1))  * e(1) + t*pos(i)*e(2);
					double d_p2 = e(0) - L*(1. / 2 - zetaZ2t(i, 1)) * e(1) + t*pos(i)*e(2);
					double d_p = e(0) - L / 2.*e(1) + t*pos(i)*e(2);

					if (d_p < 0.0 && zetaZ2t(i, 1)>0) {
						fy_p = -k*(d_py);
						f1_p = r*fc* d_p2 / (-fc / k) / muMax;
						if (zetaZ2t(i, 1) > 0.0) {
							f2_p = fy_p + (f1_p - fy_p) / (zetaZt(i, 1) - zetaZ2t(i, 1))*zetaZt(i, 1);
						}
						f1_p = r*fc* d_p / (-fc / k) / muZt(i, 1);

					}

				}
				*/


				// apply correction term to all sectional forces
/*

				if (muZt(i,0)>1.0) {     // negative side correction term

					bi(1) = L / 2.0;
					bi(2) = t*pos(i);
					c(1) = -L / 3.;
					double beta = (muMax - r) / (muMax - 1);

					if ((bi^e) < 0.0 && muZt(i, 0) < muMax) {  // apply correction as long as the the edge is compressed
							corr = -k*beta*(muZt(i, 0) - 1.0) / (2.0*muZt(i, 0))  *t*weight(i) * (zetaZt(i, 0)*L) * (bi^e);
							s += (bi + zetaZt(i, 0)*c)*corr;

							dCorr_de.Zero();
							dCorr_de += -k*(muZt(i, 0) - 1.0) / (2.0*muZt(i, 0))  *t*weight(i) * (zetaZt(i, 0)*L) * (bi);
							dCorr_de += corr / (muZt(i, 0) - 1.0) / (2.0*muZt(i, 0))* pow(muZt(i, 0), -2.0) * dmu0_de;
							dCorr_de += (corr / zetaZt(i, 0))  * dzeta0_de;

							//	Dtrial += (bi + zetaZt(i, 0)*c) % (dCorr_de);
							//	Dtrial += corr* (c % dzeta0_de);
					}

					bi(1) = L / 2.0 - zetaZ2t(i, 0)*L ;

					if ((bi^e) < 0.0 && muZt(i, 0) > muMax) {
						bi(1) = L / 2.0;
						corr = 0.5*(-k* (bi^e)-f2_n) *t*weight(i) * (zetaZt(i, 0)*L) ;
						s += (bi + zetaZt(i, 0)*c)*corr;

						c(1) = -L / 3.;
						corr = 0.5*(f1_n - f2_n)  *t*weight(i) * zetaZ2t(i, 0)*L;
						s -= (bi + zetaZ2t(i, 0)*c)*corr;

						if (x_n > 0.0) {
							c(1) = -x_n / 3.; 
							corr = 0.5*(f1_n - (-k* (bi^e)))  *t*weight(i) * x_n;
							s += (bi + c)*corr;
						}
					}
											
				}
				
				if (muZt(i,1)>1.0) {     // positive side correction term

					bi(1) = -L / 2.0;
					bi(2) = t*pos(i);
					c(1) =  L / 3.;
					double beta = (muMax - r) / (muMax - 1);

					if ((bi^e) < 0.0 && muZt(i, 1) <= muMax) {  // apply correction as long as the the edge is compressed

						corr = -k*(muZt(i, 1) - 1.0) / (2.0*muZt(i, 1))  *t*weight(i) * (zetaZt(i, 1)*L) * (bi^e);
						//s += (bi + zetaZt(i, 1)*c)*corr;

						dCorr_de.Zero();
						dCorr_de += -k*(muZt(i, 1) - 1.0) / (2.0*muZt(i, 1))  *t*weight(i) * (zetaZt(i, 1)*L) * (bi);
						dCorr_de += corr / (muZt(i, 1) - 1.0) / (2.0*muZt(i, 1))* std::pow(muZt(i, 1), -2.0) * dmu1_de;
						dCorr_de += corr / zetaZt(i, 1)   * dzeta1_de;

						//	Dtrial += (bi + zetaZt(i, 1)*c) % (dCorr_de);
						//	Dtrial += corr* (c % dzeta1_de);
						corr = -k*beta*(muZt(i, 1) - 1.0) / (2.0*muZt(i, 1))  *t*weight(i) * (zetaZt(i, 1)*L) * (bi^e);
						s += (bi + zetaZt(i, 1)*c)*corr;
					}

					bi(1) = -L / 2.0 + zetaZ2t(i, 1)*L;
					if ((bi^e) < 0.0 && muZt(i, 1) > muMax) {
						bi(1) = -L / 2.0;
						corr = 0.5*(-k* (bi^e) - f2_p) *t*weight(i) * (zetaZt(i, 1)*L);
						s += (bi + zetaZt(i, 1)*c)*corr;

						c(1) = L / 3.;
						corr = 0.5*(f1_p - f2_p)  *t*weight(i) * zetaZ2t(i, 1)*L;
						s -= (bi + zetaZ2t(i, 1)*c)*corr;

						if (x_p > 0.0) {
							c(1) = x_p / 3.;
							corr = 0.5*(f1_p - (-k* (bi^e)))  *t*weight(i) * x_p;
							s += (bi + c)*corr;
							//opserr << e;
							//opserr << x_p << ", " << zetaZ2t(i, 1) << ", " << zetaZt(i, 1) << endln;
						}

						}					
				}
				*/


			}  // end for sections
	} // end if crushing and no all tension
	// end crushing correction


	// add an elastic contribution if defined "stronger"
	if (stronger) {	
		s(0) += factorStronger* k*L*t           * e(0);
        s(1) += factorStronger* k*t*L*L*L /12.0 * e(1);
        s(2) += factorStronger* k*t*t*t*L /12.0 * e(2);
		s(3) += factorStronger* kg*J            * e(3);	
	}


	//---------------------------------------------------------------------------------------------------------//
	// set tangent
	//---------------------------------------------------------------------------------------------------------//
	// terms corresponding to the crushing correction are already added in the previous section

	if (std::abs(e(2)) < DBL_EPSILON) {
		// zero in plane moment
		if (std::abs(e(1)) + 2.0*e(0)/L > DBL_EPSILON) {
			if  (-std::abs(e(1)) + 2.0*e(0)/L > DBL_EPSILON) {
				// all tension, all zeroes, do nothing
				Dtrial.Zero();
			} else {
					// case 5
					double eps0 = e(0);    // axial deformation at the origin
					double chiZ = e(1);    // curvature around axis z
					double chiY = 0.0;     // curvature around axis y

				  if (chiZ<0.0) {
						// case 4a,c: both out from different sides, chiZ>=0
						 Dtrial(0,0) += k*t*(L/2. + eps0/chiZ); 
						 Dtrial(1,0) += -k*(t*(std::pow(chiY,2)*std::pow(t,2) - 3*std::pow(chiZ,2)*std::pow(L,2) + 12*std::pow(eps0,2)))/(24.*std::pow(chiZ,2)); 
						 Dtrial(2,0) += k*(chiY*std::pow(t,3))/(12.*chiZ); 
						 Dtrial(0,1) += -k*(t*(std::pow(chiY,2)*std::pow(t,2) - 3*std::pow(chiZ,2)*std::pow(L,2) + 12*std::pow(eps0,2)))/(24.*std::pow(chiZ,2)); 
						 Dtrial(1,1) += k*(t*(std::pow(chiZ,3)*std::pow(L,3) + 2*std::pow(chiY,2)*std::pow(t,2)*eps0 + 8*std::pow(eps0,3)))/(24.*std::pow(chiZ,3)); 
						 Dtrial(2,1) += -k*(chiY*std::pow(t,3)*eps0)/(12.*std::pow(chiZ,2)); 
						 Dtrial(0,2) += k*(chiY*std::pow(t,3))/(12.*chiZ); 
						 Dtrial(1,2) += -k*(chiY*std::pow(t,3)*eps0)/(12.*std::pow(chiZ,2)); 
						 Dtrial(2,2) += k*(std::pow(t,3)*(chiZ*L + 2*eps0))/(24.*chiZ); 

					} else {
						// case 4b,d: both out from different sides, chiZ<0
						 Dtrial(0,0) += k*(t*(L - (2*eps0)/chiZ))/2.; 
						 Dtrial(1,0) += k*(t*(std::pow(chiY,2)*std::pow(t,2) - 3*std::pow(chiZ,2)*std::pow(L,2) + 12*std::pow(eps0,2)))/(24.*std::pow(chiZ,2)); 
						 Dtrial(2,0) += -k*(chiY*std::pow(t,3))/(12.*chiZ); 
						 Dtrial(0,1) += k*(t*(std::pow(chiY,2)*std::pow(t,2) - 3*std::pow(chiZ,2)*std::pow(L,2) + 12*std::pow(eps0,2)))/(24.*std::pow(chiZ,2)); 
						 Dtrial(1,1) += k*(t*(std::pow(chiZ,3)*std::pow(L,3) - 2*std::pow(chiY,2)*std::pow(t,2)*eps0 - 8*std::pow(eps0,3)))/(24.*std::pow(chiZ,3)); 
						 Dtrial(2,1) += k*(chiY*std::pow(t,3)*eps0)/(12.*std::pow(chiZ,2)); 
						 Dtrial(0,2) += -k*(chiY*std::pow(t,3))/(12.*chiZ); 
						 Dtrial(1,2) += k*(chiY*std::pow(t,3)*eps0)/(12.*std::pow(chiZ,2)); 
						 Dtrial(2,2) += k*(std::pow(t,3)*(chiZ*L - 2*eps0))/(24.*chiZ); 
				    }
			}

		} else {
			Dtrial(0,0) += k*L*t;
            Dtrial(1,1) += k*t*L*L*L /12.0;
            Dtrial(2,2) += k*t*t*t*L /12.0;
		}
	} else {
		double eps0 = e(0);    // axial deformation at the origin
		double chiZ = e(1);    // curvature around axis z
		double chiY = e(2);    // curvature around axis y

		double t1 = (2.0*eps0 + L*chiZ)/(2.0*chiY);   // x coordinate of the neutral axis for y = +L/2
		double t2 = (2.0*eps0 - L*chiZ)/(2.0*chiY);   // x coordinate of the neutral axis for y = -L/2

		int inside1 = 0;
		if (t1<-t/2.0)   inside1=-1;
		if (t1> t/2.0)   inside1=+1;

	    int inside2 = 0;
		if (t2<-t/2.0)   inside2=-1;
		if (t2> t/2.0)   inside2=+1;
		
		if (inside1==0 && inside2==0) {
			if (chiY>0.0) {
			   // case 1a,b: both neutral points inside the long side, fi_y>0
			     Dtrial(0,0) += k*(L*(t - (2*eps0)/chiY))/2.; 
				 Dtrial(1,0) += -k*(chiZ*std::pow(L,3))/(12.*chiY); 
				 Dtrial(2,0) += k*(L*(-3*std::pow(chiY,2)*std::pow(t,2) + std::pow(chiZ,2)*std::pow(L,2) + 12*std::pow(eps0,2)))/(24.*std::pow(chiY,2)); 
				 Dtrial(0,1) += -k*(chiZ*std::pow(L,3))/(12.*chiY); 
				 Dtrial(1,1) += k*(std::pow(L,3)*(chiY*t - 2*eps0))/(24.*chiY); 
				 Dtrial(2,1) += k*(chiZ*std::pow(L,3)*eps0)/(12.*std::pow(chiY,2)); 
				 Dtrial(0,2) += k*(L*(-3*std::pow(chiY,2)*std::pow(t,2) + std::pow(chiZ,2)*std::pow(L,2) + 12*std::pow(eps0,2)))/(24.*std::pow(chiY,2)); 
				 Dtrial(1,2) += k*(chiZ*std::pow(L,3)*eps0)/(12.*std::pow(chiY,2)); 
				 Dtrial(2,2) += k*(L*(std::pow(chiY,3)*std::pow(t,3) - 2*std::pow(chiZ,2)*std::pow(L,2)*eps0 - 8*std::pow(eps0,3)))/(24.*std::pow(chiY,3)); 

			} else {
               // case 1c,d: both neutral points inside the long side, fi_y<0
				 Dtrial(0,0) += k*(t*L)/2. + k*(L*eps0)/chiY; 
				 Dtrial(1,0) += k*(chiZ*std::pow(L,3))/(12.*chiY); 
				 Dtrial(2,0) += -k*(L*(-3*std::pow(chiY,2)*std::pow(t,2) + std::pow(chiZ,2)*std::pow(L,2) + 12*std::pow(eps0,2)))/(24.*std::pow(chiY,2)); 
				 Dtrial(0,1) += k*(chiZ*std::pow(L,3))/(12.*chiY); 
				 Dtrial(1,1) += k*(std::pow(L,3)*(chiY*t + 2*eps0))/(24.*chiY); 
				 Dtrial(2,1) += -k*(chiZ*std::pow(L,3)*eps0)/(12.*std::pow(chiY,2)); 
				 Dtrial(0,2) += -k*(L*(-3*std::pow(chiY,2)*std::pow(t,2) + std::pow(chiZ,2)*std::pow(L,2) + 12*std::pow(eps0,2)))/(24.*std::pow(chiY,2)); 
				 Dtrial(1,2) += -k*(chiZ*std::pow(L,3)*eps0)/(12.*std::pow(chiY,2)); 
				 Dtrial(2,2) += k*(L*(std::pow(chiY,3)*std::pow(t,3) + 2*std::pow(chiZ,2)*std::pow(L,2)*eps0 + 8*std::pow(eps0,3)))/(24.*std::pow(chiY,3)); 
			}

		} else {
			if (inside1*inside2>0) {
				if (eps0 - t/2.0*chiY + L/2.0*chiZ < DBL_EPSILON) {
				    // case 0a: both neutral points ouside the same side, all compressed
					Dtrial(0,0) += k*L*t;
                    Dtrial(1,1) += k*t*L*L*L /12.0;
                    Dtrial(2,2) += k*t*t*t*L /12.0;
				} else {
					// case 0b: both neutral points ouside the same side, all tension
				    // all zeroes, do nothing
				}
			} else {
				if (inside1*inside2==0) {
					// one in and one out. Cases 2-3
					if (inside2>0) {
						if (chiY>0.0) {
							// case 2a
							Dtrial(0,0) += -k*std::pow(-(chiY*t) + chiZ*L + 2*eps0,2)/(8.*chiZ*chiY); 
							Dtrial(1,0) += -k*((chiY*t + 2*chiZ*L - 2*eps0)*std::pow(-(chiY*t) + chiZ*L + 2*eps0,2))/(48.*std::pow(chiZ,2)*chiY); 
							Dtrial(2,0) += k*(std::pow(-(chiY*t) + chiZ*L + 2*eps0,2)*(2*chiY*t + chiZ*L + 2*eps0))/(48.*chiZ*std::pow(chiY,2)); 
							Dtrial(0,1) += -k*((chiY*t + 2*chiZ*L - 2*eps0)*std::pow(-(chiY*t) + chiZ*L + 2*eps0,2))/(48.*std::pow(chiZ,2)*chiY); 
							Dtrial(1,1) += -k*(std::pow(-(chiY*t) + chiZ*L + 2*eps0,2)*(std::pow(chiY,2)*std::pow(t,2) + 2*chiZ*chiY*t*L + 3*std::pow(chiZ,2)*std::pow(L,2) - 4*(chiY*t + chiZ*L)*eps0 + 4*std::pow(eps0,2)))/(192.*std::pow(chiZ,3)*chiY); 
							Dtrial(2,1) += k*(std::pow(-(chiY*t) + chiZ*L + 2*eps0,2)*(3*std::pow(chiY*t + chiZ*L,2) - 4*chiY*t*eps0 + 4*chiZ*L*eps0 - 4*std::pow(eps0,2)))/(384.*std::pow(chiZ,2)*std::pow(chiY,2)); 
							Dtrial(0,2) += k*(std::pow(-(chiY*t) + chiZ*L + 2*eps0,2)*(2*chiY*t + chiZ*L + 2*eps0))/(48.*chiZ*std::pow(chiY,2)); 
							Dtrial(1,2) += k*(std::pow(-(chiY*t) + chiZ*L + 2*eps0,2)*(3*std::pow(chiY*t + chiZ*L,2) - 4*chiY*t*eps0 + 4*chiZ*L*eps0 - 4*std::pow(eps0,2)))/(384.*std::pow(chiZ,2)*std::pow(chiY,2)); 
							Dtrial(2,2) += -k*(std::pow(-(chiY*t) + chiZ*L + 2*eps0,2)*(3*std::pow(chiY,2)*std::pow(t,2) + 2*chiY*t*(chiZ*L + 2*eps0) + std::pow(chiZ*L + 2*eps0,2)))/(192.*chiZ*std::pow(chiY,3)); 

						} else {
							// case 3a
							 Dtrial(0,0) += k*t*L + std::pow(-(chiY*t) + chiZ*L + 2*eps0,2)/(8.*chiZ*chiY); 
							 Dtrial(1,0) += k*((chiY*t + 2*chiZ*L - 2*eps0)*std::pow(-(chiY*t) + chiZ*L + 2*eps0,2))/(48.*std::pow(chiZ,2)*chiY); 
							 Dtrial(2,0) += -k*(std::pow(-(chiY*t) + chiZ*L + 2*eps0,2)*(2*chiY*t + chiZ*L + 2*eps0))/(48.*chiZ*std::pow(chiY,2)); 
							 Dtrial(0,1) += k*((chiY*t + 2*chiZ*L - 2*eps0)*std::pow(-(chiY*t) + chiZ*L + 2*eps0,2))/(48.*std::pow(chiZ,2)*chiY); 
							 Dtrial(1,1) += k*(std::pow(chiY,4)*std::pow(t,4) + 3*std::pow(chiZ,4)*std::pow(L,4) - 8*std::pow(chiY,3)*std::pow(t,3)*eps0 + 8*std::pow(chiZ,3)*std::pow(L,3)*eps0 + 24*std::pow(chiY,2)*std::pow(t,2)*std::pow(eps0,2) + 16*std::pow(eps0,4) + 4*chiY*t*(3*std::pow(chiZ,3)*std::pow(L,3) - 8*std::pow(eps0,3)))/(192.*std::pow(chiZ,3)*chiY); 
							 Dtrial(2,1) += -k*(std::pow(-(chiY*t) + chiZ*L + 2*eps0,2)*(3*std::pow(chiY*t + chiZ*L,2) - 4*chiY*t*eps0 + 4*chiZ*L*eps0 - 4*std::pow(eps0,2)))/(384.*std::pow(chiZ,2)*std::pow(chiY,2)); 
							 Dtrial(0,2) += -k*(std::pow(-(chiY*t) + chiZ*L + 2*eps0,2)*(2*chiY*t + chiZ*L + 2*eps0))/(48.*chiZ*std::pow(chiY,2)); 
							 Dtrial(1,2) += -k*(std::pow(-(chiY*t) + chiZ*L + 2*eps0,2)*(3*std::pow(chiY*t + chiZ*L,2) - 4*chiY*t*eps0 + 4*chiZ*L*eps0 - 4*std::pow(eps0,2)))/(384.*std::pow(chiZ,2)*std::pow(chiY,2)); 
							 Dtrial(2,2) += k*(3*std::pow(chiY,4)*std::pow(t,4) + 4*std::pow(chiY,3)*std::pow(t,3)*(3*chiZ*L - 2*eps0) + std::pow(chiZ*L + 2*eps0,4))/(192.*chiZ*std::pow(chiY,3)); 
						}

					} else {
						if (inside2<0) {
							if (chiY<0.0) {
								// case 2c
								 Dtrial(0,0) += k*std::pow(chiY*t + chiZ*L + 2*eps0,2)/(8.*chiZ*chiY); 
								 Dtrial(1,0) += -k*((chiY*t - 2*chiZ*L + 2*eps0)*std::pow(chiY*t + chiZ*L + 2*eps0,2))/(48.*std::pow(chiZ,2)*chiY); 
								 Dtrial(2,0) += -k*((-2*chiY*t + chiZ*L + 2*eps0)*std::pow(chiY*t + chiZ*L + 2*eps0,2))/(48.*chiZ*std::pow(chiY,2)); 
								 Dtrial(0,1) += -k*((chiY*t - 2*chiZ*L + 2*eps0)*std::pow(chiY*t + chiZ*L + 2*eps0,2))/(48.*std::pow(chiZ,2)*chiY); 
								 Dtrial(1,1) += k*(std::pow(chiY*t + chiZ*L + 2*eps0,2)*(std::pow(chiY,2)*std::pow(t,2) - 2*chiZ*chiY*t*L + 3*std::pow(chiZ,2)*std::pow(L,2) + 4*chiY*t*eps0 - 4*chiZ*L*eps0 + 4*std::pow(eps0,2)))/(192.*std::pow(chiZ,3)*chiY); 
								 Dtrial(2,1) += -k*(std::pow(chiY*t + chiZ*L + 2*eps0,2)*(3*std::pow(chiY*t - chiZ*L,2) + 4*(chiY*t + chiZ*L)*eps0 - 4*std::pow(eps0,2)))/(384.*std::pow(chiZ,2)*std::pow(chiY,2)); 
								 Dtrial(0,2) += -k*((-2*chiY*t + chiZ*L + 2*eps0)*std::pow(chiY*t + chiZ*L + 2*eps0,2))/(48.*chiZ*std::pow(chiY,2)); 
								 Dtrial(1,2) += -k*(std::pow(chiY*t + chiZ*L + 2*eps0,2)*(3*std::pow(chiY*t - chiZ*L,2) + 4*(chiY*t + chiZ*L)*eps0 - 4*std::pow(eps0,2)))/(384.*std::pow(chiZ,2)*std::pow(chiY,2)); 
								 Dtrial(2,2) += k*(std::pow(chiY*t + chiZ*L + 2*eps0,2)*(3*std::pow(chiY,2)*std::pow(t,2) - 2*chiY*t*(chiZ*L + 2*eps0) + std::pow(chiZ*L + 2*eps0,2)))/(192.*chiZ*std::pow(chiY,3)); 

							} else {
								// case 3c
									 Dtrial(0,0) += k*t*L - k*std::pow(chiY*t + chiZ*L + 2*eps0,2)/(8.*chiZ*chiY); 
									 Dtrial(1,0) += k*((chiY*t - 2*chiZ*L + 2*eps0)*std::pow(chiY*t + chiZ*L + 2*eps0,2))/(48.*std::pow(chiZ,2)*chiY); 
									 Dtrial(2,0) += k*((-2*chiY*t + chiZ*L + 2*eps0)*std::pow(chiY*t + chiZ*L + 2*eps0,2))/(48.*chiZ*std::pow(chiY,2)); 
									 Dtrial(0,1) += k*((chiY*t - 2*chiZ*L + 2*eps0)*std::pow(chiY*t + chiZ*L + 2*eps0,2))/(48.*std::pow(chiZ,2)*chiY); 
									 Dtrial(1,1) += -k*(std::pow(chiY,4)*std::pow(t,4) + 3*std::pow(chiZ,4)*std::pow(L,4) + 8*std::pow(chiY,3)*std::pow(t,3)*eps0 + 8*std::pow(chiZ,3)*std::pow(L,3)*eps0 + 24*std::pow(chiY,2)*std::pow(t,2)*std::pow(eps0,2) + 16*std::pow(eps0,4) + 4*chiY*t*(-3*std::pow(chiZ,3)*std::pow(L,3) + 8*std::pow(eps0,3)))/(192.*std::pow(chiZ,3)*chiY); 
									 Dtrial(2,1) += k*(std::pow(chiY*t + chiZ*L + 2*eps0,2)*(3*std::pow(chiY*t - chiZ*L,2) + 4*(chiY*t + chiZ*L)*eps0 - 4*std::pow(eps0,2)))/(384.*std::pow(chiZ,2)*std::pow(chiY,2)); 
									 Dtrial(0,2) += k*((-2*chiY*t + chiZ*L + 2*eps0)*std::pow(chiY*t + chiZ*L + 2*eps0,2))/(48.*chiZ*std::pow(chiY,2)); 
									 Dtrial(1,2) += k*(std::pow(chiY*t + chiZ*L + 2*eps0,2)*(3*std::pow(chiY*t - chiZ*L,2) + 4*(chiY*t + chiZ*L)*eps0 - 4*std::pow(eps0,2)))/(384.*std::pow(chiZ,2)*std::pow(chiY,2)); 
									 Dtrial(2,2) += -k*(3*std::pow(chiY,4)*std::pow(t,4) + 4*std::pow(chiY,3)*std::pow(t,3)*(-3*chiZ*L + 2*eps0) + std::pow(chiZ*L + 2*eps0,4))/(192.*chiZ*std::pow(chiY,3)); 

							}

						} else {
							if (inside1>0) {
								if (chiY>0.0) {
									// case 2b
									Dtrial(0,0) += k*std::pow(chiY*t + chiZ*L - 2*eps0,2)/(8.*chiZ*chiY); 
									Dtrial(1,0) += k*(std::pow(chiY*t + chiZ*L - 2*eps0,2)*(chiY*t - 2*(chiZ*L + eps0)))/(48.*std::pow(chiZ,2)*chiY); 
									Dtrial(2,0) += k*((-2*chiY*t + chiZ*L - 2*eps0)*std::pow(chiY*t + chiZ*L - 2*eps0,2))/(48.*chiZ*std::pow(chiY,2)); 
									Dtrial(0,1) += k*(std::pow(chiY*t + chiZ*L - 2*eps0,2)*(chiY*t - 2*(chiZ*L + eps0)))/(48.*std::pow(chiZ,2)*chiY); 
									Dtrial(1,1) += k*(std::pow(chiY*t + chiZ*L - 2*eps0,2)*(std::pow(chiY,2)*std::pow(t,2) + 3*std::pow(chiZ,2)*std::pow(L,2) + 4*chiZ*L*eps0 + 4*std::pow(eps0,2) - 2*chiY*t*(chiZ*L + 2*eps0)))/(192.*std::pow(chiZ,3)*chiY); 
									Dtrial(2,1) += -k*(std::pow(chiY*t + chiZ*L - 2*eps0,2)*(3*std::pow(chiY*t - chiZ*L,2) - 4*(chiY*t + chiZ*L)*eps0 - 4*std::pow(eps0,2)))/(384.*std::pow(chiZ,2)*std::pow(chiY,2)); 
									Dtrial(0,2) += k*((-2*chiY*t + chiZ*L - 2*eps0)*std::pow(chiY*t + chiZ*L - 2*eps0,2))/(48.*chiZ*std::pow(chiY,2)); 
									Dtrial(1,2) += -k*(std::pow(chiY*t + chiZ*L - 2*eps0,2)*(3*std::pow(chiY*t - chiZ*L,2) - 4*(chiY*t + chiZ*L)*eps0 - 4*std::pow(eps0,2)))/(384.*std::pow(chiZ,2)*std::pow(chiY,2)); 
									Dtrial(2,2) += k*(std::pow(chiY*t + chiZ*L - 2*eps0,2)*(3*std::pow(chiY,2)*std::pow(t,2) + std::pow(chiZ*L - 2*eps0,2) + chiY*(-2*chiZ*t*L + 4*t*eps0)))/(192.*chiZ*std::pow(chiY,3)); 


								} else {
									// case 3b
									 Dtrial(0,0) += k*t*L - std::pow(chiY*t + chiZ*L - 2*eps0,2)/(8.*chiZ*chiY); 
									 Dtrial(1,0) += k*(std::pow(chiY*t + chiZ*L - 2*eps0,2)*(-(chiY*t) + 2*(chiZ*L + eps0)))/(48.*std::pow(chiZ,2)*chiY); 
									 Dtrial(2,0) += k*(std::pow(chiY*t + chiZ*L - 2*eps0,2)*(2*chiY*t - chiZ*L + 2*eps0))/(48.*chiZ*std::pow(chiY,2)); 
									 Dtrial(0,1) += k*(std::pow(chiY*t + chiZ*L - 2*eps0,2)*(-(chiY*t) + 2*(chiZ*L + eps0)))/(48.*std::pow(chiZ,2)*chiY); 
									 Dtrial(1,1) += -k*(std::pow(chiY,4)*std::pow(t,4) + 3*std::pow(chiZ,4)*std::pow(L,4) - 8*std::pow(chiY,3)*std::pow(t,3)*eps0 - 8*std::pow(chiZ,3)*std::pow(L,3)*eps0 + 24*std::pow(chiY,2)*std::pow(t,2)*std::pow(eps0,2) + 16*std::pow(eps0,4) - 4*chiY*t*(3*std::pow(chiZ,3)*std::pow(L,3) + 8*std::pow(eps0,3)))/(192.*std::pow(chiZ,3)*chiY); 
									 Dtrial(2,1) += k*(std::pow(chiY*t + chiZ*L - 2*eps0,2)*(3*std::pow(chiY*t - chiZ*L,2) - 4*(chiY*t + chiZ*L)*eps0 - 4*std::pow(eps0,2)))/(384.*std::pow(chiZ,2)*std::pow(chiY,2)); 
									 Dtrial(0,2) += k*(std::pow(chiY*t + chiZ*L - 2*eps0,2)*(2*chiY*t - chiZ*L + 2*eps0))/(48.*chiZ*std::pow(chiY,2)); 
									 Dtrial(1,2) += k*(std::pow(chiY*t + chiZ*L - 2*eps0,2)*(3*std::pow(chiY*t - chiZ*L,2) - 4*(chiY*t + chiZ*L)*eps0 - 4*std::pow(eps0,2)))/(384.*std::pow(chiZ,2)*std::pow(chiY,2)); 
									 Dtrial(2,2) += -k*(3*std::pow(chiY,4)*std::pow(t,4) + std::pow(chiZ*L - 2*eps0,4) - 4*std::pow(chiY,3)*std::pow(t,3)*(3*chiZ*L + 2*eps0))/(192.*chiZ*std::pow(chiY,3)); 
								}

							} else {
								if (chiY<0.0) {
									// case 2d
									 Dtrial(0,0) += -k*std::pow(chiY*t - chiZ*L + 2*eps0,2)/(8.*chiZ*chiY); 
									 Dtrial(1,0) += k*(std::pow(chiY*t - chiZ*L + 2*eps0,2)*(chiY*t + 2*(chiZ*L + eps0)))/(48.*std::pow(chiZ,2)*chiY); 
									 Dtrial(2,0) += -k*((2*chiY*t + chiZ*L - 2*eps0)*std::pow(chiY*t - chiZ*L + 2*eps0,2))/(48.*chiZ*std::pow(chiY,2)); 
									 Dtrial(0,1) += k*(std::pow(chiY*t - chiZ*L + 2*eps0,2)*(chiY*t + 2*(chiZ*L + eps0)))/(48.*std::pow(chiZ,2)*chiY); 
									 Dtrial(1,1) += -k*(std::pow(chiY*t - chiZ*L + 2*eps0,2)*(std::pow(chiY,2)*std::pow(t,2) + 3*std::pow(chiZ,2)*std::pow(L,2) + 4*chiZ*L*eps0 + 4*std::pow(eps0,2) + 2*chiY*t*(chiZ*L + 2*eps0)))/(192.*std::pow(chiZ,3)*chiY); 
									 Dtrial(2,1) += k*(std::pow(chiY*t - chiZ*L + 2*eps0,2)*(3*std::pow(chiY*t + chiZ*L,2) + 4*(chiY*t - chiZ*L)*eps0 - 4*std::pow(eps0,2)))/(384.*std::pow(chiZ,2)*std::pow(chiY,2)); 
									 Dtrial(0,2) += -k*((2*chiY*t + chiZ*L - 2*eps0)*std::pow(chiY*t - chiZ*L + 2*eps0,2))/(48.*chiZ*std::pow(chiY,2)); 
									 Dtrial(1,2) += k*(std::pow(chiY*t - chiZ*L + 2*eps0,2)*(3*std::pow(chiY*t + chiZ*L,2) + 4*(chiY*t - chiZ*L)*eps0 - 4*std::pow(eps0,2)))/(384.*std::pow(chiZ,2)*std::pow(chiY,2)); 
									 Dtrial(2,2) += -k*((3*std::pow(chiY,2)*std::pow(t,2) + 2*chiY*t*(chiZ*L - 2*eps0) + std::pow(chiZ*L - 2*eps0,2))*std::pow(chiY*t - chiZ*L + 2*eps0,2))/(192.*chiZ*std::pow(chiY,3)); 
								} else {
									// case 3d
									 Dtrial(0,0) += k*t*L + k*std::pow(chiY*t - chiZ*L + 2*eps0,2)/(8.*chiZ*chiY); 
									 Dtrial(1,0) += -k*(std::pow(chiY*t - chiZ*L + 2*eps0,2)*(chiY*t + 2*(chiZ*L + eps0)))/(48.*std::pow(chiZ,2)*chiY); 
									 Dtrial(2,0) += k*((2*chiY*t + chiZ*L - 2*eps0)*std::pow(chiY*t - chiZ*L + 2*eps0,2))/(48.*chiZ*std::pow(chiY,2)); 
									 Dtrial(0,1) += -k*(std::pow(chiY*t - chiZ*L + 2*eps0,2)*(chiY*t + 2*(chiZ*L + eps0)))/(48.*std::pow(chiZ,2)*chiY); 
									 Dtrial(1,1) += k*(std::pow(chiY,4)*std::pow(t,4) + 3*std::pow(chiZ,4)*std::pow(L,4) + 8*std::pow(chiY,3)*std::pow(t,3)*eps0 - 8*std::pow(chiZ,3)*std::pow(L,3)*eps0 + 24*std::pow(chiY,2)*std::pow(t,2)*std::pow(eps0,2) + 16*std::pow(eps0,4) + 4*chiY*t*(3*std::pow(chiZ,3)*std::pow(L,3) + 8*std::pow(eps0,3)))/(192.*std::pow(chiZ,3)*chiY); 
									 Dtrial(2,1) += -k*(std::pow(chiY*t - chiZ*L + 2*eps0,2)*(3*std::pow(chiY*t + chiZ*L,2) + 4*(chiY*t - chiZ*L)*eps0 - 4*std::pow(eps0,2)))/(384.*std::pow(chiZ,2)*std::pow(chiY,2)); 
									 Dtrial(0,2) += k*((2*chiY*t + chiZ*L - 2*eps0)*std::pow(chiY*t - chiZ*L + 2*eps0,2))/(48.*chiZ*std::pow(chiY,2)); 
									 Dtrial(1,2) += -k*(std::pow(chiY*t - chiZ*L + 2*eps0,2)*(3*std::pow(chiY*t + chiZ*L,2) + 4*(chiY*t - chiZ*L)*eps0 - 4*std::pow(eps0,2)))/(384.*std::pow(chiZ,2)*std::pow(chiY,2)); 
									 Dtrial(2,2) += k*(3*std::pow(chiY,4)*std::pow(t,4) + std::pow(chiZ*L - 2*eps0,4) + 4*std::pow(chiY,3)*std::pow(t,3)*(3*chiZ*L + 2*eps0))/(192.*chiZ*std::pow(chiY,3)); 
								}
							}
						}
					}

				} else {
					if (chiZ<0.0) {
						// case 4a,c: both out from different sides, chiZ>=0
						 Dtrial(0,0) += k*t*(L/2. + eps0/chiZ); 
						 Dtrial(1,0) += -k*(t*(std::pow(chiY,2)*std::pow(t,2) - 3*std::pow(chiZ,2)*std::pow(L,2) + 12*std::pow(eps0,2)))/(24.*std::pow(chiZ,2)); 
						 Dtrial(2,0) += k*(chiY*std::pow(t,3))/(12.*chiZ); 
						 Dtrial(0,1) += -k*(t*(std::pow(chiY,2)*std::pow(t,2) - 3*std::pow(chiZ,2)*std::pow(L,2) + 12*std::pow(eps0,2)))/(24.*std::pow(chiZ,2)); 
						 Dtrial(1,1) += k*(t*(std::pow(chiZ,3)*std::pow(L,3) + 2*std::pow(chiY,2)*std::pow(t,2)*eps0 + 8*std::pow(eps0,3)))/(24.*std::pow(chiZ,3)); 
						 Dtrial(2,1) += -k*(chiY*std::pow(t,3)*eps0)/(12.*std::pow(chiZ,2)); 
						 Dtrial(0,2) += k*(chiY*std::pow(t,3))/(12.*chiZ); 
						 Dtrial(1,2) += -k*(chiY*std::pow(t,3)*eps0)/(12.*std::pow(chiZ,2)); 
						 Dtrial(2,2) += k*(std::pow(t,3)*(chiZ*L + 2*eps0))/(24.*chiZ); 

					} else {
						// case 4b,d: both out from different sides, chiZ<0
						 Dtrial(0,0) += k*(t*(L - (2*eps0)/chiZ))/2.; 
						 Dtrial(1,0) += k*(t*(std::pow(chiY,2)*std::pow(t,2) - 3*std::pow(chiZ,2)*std::pow(L,2) + 12*std::pow(eps0,2)))/(24.*std::pow(chiZ,2)); 
						 Dtrial(2,0) += -k*(chiY*std::pow(t,3))/(12.*chiZ); 
						 Dtrial(0,1) += k*(t*(std::pow(chiY,2)*std::pow(t,2) - 3*std::pow(chiZ,2)*std::pow(L,2) + 12*std::pow(eps0,2)))/(24.*std::pow(chiZ,2)); 
						 Dtrial(1,1) += k*(t*(std::pow(chiZ,3)*std::pow(L,3) - 2*std::pow(chiY,2)*std::pow(t,2)*eps0 - 8*std::pow(eps0,3)))/(24.*std::pow(chiZ,3)); 
						 Dtrial(2,1) += k*(chiY*std::pow(t,3)*eps0)/(12.*std::pow(chiZ,2)); 
						 Dtrial(0,2) += -k*(chiY*std::pow(t,3))/(12.*chiZ); 
						 Dtrial(1,2) += k*(chiY*std::pow(t,3)*eps0)/(12.*std::pow(chiZ,2)); 
						 Dtrial(2,2) += k*(std::pow(t,3)*(chiZ*L - 2*eps0))/(24.*chiZ); 

					}
				}
			}
		}
	}

	Dtrial(3,3) = kg*J;
	Dtrial(3,3)*= torsionalStiffnessFactor;

    if (stronger) {		
		Dtrial(0,0) += factorStronger*k*L*t;
        Dtrial(1,1) += factorStronger*k*t*L*L*L /12.0;
        Dtrial(2,2) += factorStronger*k*t*t*t*L /12.0;
		Dtrial(3,3) += factorStronger*kg*J;

	} else {
		// fake stiffness instead of zero
		
		double factor = 0.000001; 
		if  (std::sqrt(std::pow(s(0),2) + std::pow(s(1),2) + std::pow(s(1),2))<DBL_EPSILON) {

			Dtrial(0,0) += factor *k* L*t;
			Dtrial(1,1) += factor *k* t*L*L*L /12.0;
			Dtrial(2,2) += factor *k* t*t*t*L /12.0;			
	     }
	}

    // maybe implement different conditions for spandrel elements?
	if (spandrel) {
	}

	return 0;
}

const Vector &
NoTensionSection3d::getSectionDeformation (void)
{
    return e;
}

const Vector &
NoTensionSection3d::getStressResultant (void)
{
  	if (std::abs(s(0))>-DBL_EPSILON && 
		std::abs(s(1))>-DBL_EPSILON &&
	    std::abs(s(2))>-DBL_EPSILON ) {
		return s;
	}else {
	  //return s*0.0; // This return a reference to a temporary (undefined
	  //behavior)
	  return s_zero;
	}

  return s;
}

const Matrix &
NoTensionSection3d::getSectionTangent(void)
{
	if (std::abs(Dtrial(0,0))>-DBL_EPSILON && 
		std::abs(Dtrial(0,1))>-DBL_EPSILON &&
		std::abs(Dtrial(0,2))>-DBL_EPSILON &&
				
		std::abs(Dtrial(1,0))>-DBL_EPSILON &&
		std::abs(Dtrial(1,1))>-DBL_EPSILON &&
		std::abs(Dtrial(1,2))>-DBL_EPSILON &&
		
		std::abs(Dtrial(2,0))>-DBL_EPSILON &&
		std::abs(Dtrial(2,1))>-DBL_EPSILON &&
		std::abs(Dtrial(2,2))>-DBL_EPSILON )		
	
     	return Dtrial;
	else
		return D;

}

const Matrix &
NoTensionSection3d::getInitialTangent(void)
{  
  return D;
}


SectionForceDeformation*
NoTensionSection3d::getCopy ()
{
    // Make a copy of the hinge
    NoTensionSection3d *theCopy = new NoTensionSection3d (this->getTag(), k, kg, t, L, J, fc, nSections, stronger, elastic, crushing, spandrel, r, muMax,triangular);
    theCopy->parameterID = parameterID;

    return theCopy;
}

const ID&
NoTensionSection3d::getType ()
{
    return code;
}

int
NoTensionSection3d::getOrder () const
{
    return 4;
}

int
NoTensionSection3d::sendSelf(int commitTag, Channel &theChannel)
{
    int res = 0;

    static Vector data(6);

    int dataTag = this->getDbTag();
    
	data(0) = this->getTag();
    data(1) = k;
    data(2) = kg;    
    data(3) = t;
    data(4) = L;
    data(5) = J;
    
    res += theChannel.sendVector(dataTag, commitTag, data);
    if(res < 0) {
      opserr << "ElasticSection3d::sendSelf -- failed to send data\n";
      return res;
    }
    
    return res;
}

int
NoTensionSection3d::recvSelf(int commitTag, Channel &theChannel,
					 FEM_ObjectBroker &theBroker)
{
    int res = 0;
    
	static Vector data(6);

    int dataTag = this->getDbTag();

    res += theChannel.recvVector(dataTag, commitTag, data);
    if(res < 0) {
      opserr << "ElasticSection3d::recvSelf -- failed to receive data\n";
      return res;
    }

	this->setTag((int)data(0));
    k = data(1);
    kg = data(2);    
    t = data(3);
    L = data(4);
    J = data(5);

    return res;
}
 
void
NoTensionSection3d::Print(OPS_Stream &s, int flag)
{
  if (flag == 2) {

  } else {
    s << "NoTensionSection3d, tag: " << this->getTag() << endln;
    s << " k: " << k << endln;
    s << " kg: " << kg << endln;
    s << " t: " << t << endln;
    s << " L: " << L << endln;
    s << " J: " << J << endln;
  }
}

int
NoTensionSection3d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  if (strcmp(argv[0],"k") == 0) {
    param.setValue(k);
    return param.addObject(1, this);
  }
  if (strcmp(argv[0],"kg") == 0) {
    param.setValue(kg);
    return param.addObject(2, this);
  }
  if (strcmp(argv[0],"t") == 0) {
    param.setValue(t);
    return param.addObject(3, this);
  }
  if (strcmp(argv[0],"L") == 0) {
    param.setValue(L);
    return param.addObject(4, this);
  }
  if (strcmp(argv[0],"J") == 0) {
    param.setValue(J);
    return param.addObject(5, this);
  }
  return -1;
}

int
NoTensionSection3d::updateParameter(int paramID, Information &info)
{
  if (paramID == 1)
    k = info.theDouble;
  if (paramID == 2)
    kg = info.theDouble;
  if (paramID == 3)
    t = info.theDouble;
  if (paramID == 4)
    L = info.theDouble;
  if (paramID == 5)
    J = info.theDouble;

  return 0;
}

int
NoTensionSection3d::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}

// recorder methods
Vector
NoTensionSection3d::getCrushingVariables(int sliceNum) {
	Vector crushingVarSet(4);

	crushingVarSet(0) = muZ(sliceNum-1, 0);
	crushingVarSet(1) = zetaZ(sliceNum-1, 0);
	crushingVarSet(2) = muZ(sliceNum-1, 1);
	crushingVarSet(3) = zetaZ(sliceNum-1, 1);

	return crushingVarSet;
}

Response*
NoTensionSection3d::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  const ID &type = this->getType();
  int typeSize = this->getOrder();

  Response *theResponse = 0;

  output.tag("SectionOutput");
  output.attr("secType", this->getClassType());
  output.attr("secTag", this->getTag());

  // deformations
  if (strcmp(argv[0],"deformations") == 0 || strcmp(argv[0],"deformation") == 0) {

         theResponse =  new MaterialResponse(this, 1, this->getSectionDeformation());
  
  // forces
  } else if (strcmp(argv[0],"forces") == 0 || strcmp(argv[0],"force") == 0) {
    
         theResponse =  new MaterialResponse(this, 2, this->getStressResultant());
  
  // force and deformation
  } else if (strcmp(argv[0],"forceAndDeformation") == 0) { 

         theResponse =  new MaterialResponse(this, 4, Vector(2*this->getOrder()));
  // state parameters
  } else if (strcmp(argv[0],"state") == 0 || strcmp(argv[0],"muZeta")==0 || strcmp(argv[0],"crushing") == 0) {
	  
	int sliceNum = nSections/2 +1;
	
	if (argc > 1) 
		sliceNum = atoi(argv[1]);
      
	if (sliceNum > 0 && sliceNum <= nSections) {
			   
		output.tag("SliceOutput");
		output.attr("number",sliceNum);

		sliceOutput = sliceNum;
		theResponse = new MaterialResponse(this, 5, this->getCrushingVariables(sliceNum) );
		output.endTag();
      
    } else if (sliceNum == 0) { // argv[1] was not an int, we want the central slice only
	
		sliceNum = nSections/2 +1;

		output.tag("SliceOutput");
		output.attr("number",sliceNum);

		sliceOutput = sliceNum;
		theResponse = new MaterialResponse(this, 5, this->getCrushingVariables(sliceNum) );
		output.endTag();

	}
   

  }

  output.endTag();
  return theResponse;
}

int 
NoTensionSection3d::getResponse(int responseID, Information &sectInfo)
{
  if (responseID == 5) { 
    Vector data(4);
		return sectInfo.setVector(this->getCrushingVariables(sliceOutput));	
  } else {  
		return SectionForceDeformation::getResponse(responseID, sectInfo);
  }
}
