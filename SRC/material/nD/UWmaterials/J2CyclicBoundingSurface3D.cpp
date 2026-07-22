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
                                                                       
// Written: UW Computational Geomechanics Group
//          Pedro Arduino (*), ALborz Ghofrani(*)
//            (*)  University of Washington
//          March 2020
//
// Description: This file contains the implementation of the J2CyclicBoundingSurface3D class.
#include "J2CyclicBoundingSurface3D.h"



#include <elementAPI.h>

void * OPS_ADD_RUNTIME_VPV(OPS_J2CyclicBoundingSurfaceMaterial)
{
    int numdata = OPS_GetNumRemainingInputArgs();

    if (numdata < 10) {
        opserr << "WARNING: Insufficient arguements\n";
        opserr << "Want: nDMaterial J2CyclicBoundingSurface tag? G? K? su? rho? h? m? k_in?  chi? beta?\n";
        return 0;
    }

    int tag;

    numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) {
        opserr << "WARNING invalid J2CyclicBoundingSurface tag\n";
        return 0;
    }

    double data[9] = { 0,0,0,0,0,0,0,0,0 };
    numdata = OPS_GetNumRemainingInputArgs();
    if (numdata != 9) {
        opserr << "WARNING error in  J2CyclicBoundingSurface number of arg incorrect\n";
        return 0;
    }
    if (OPS_GetDoubleInput(&numdata, data)) {
        opserr << "WARNING invalid J2CyclicBoundingSurface double inputs\n";
        return 0;
    }

    return new J2CyclicBoundingSurface3D(tag, data[0], data[1], data[2], data[3], data[4], data[5], data[6], data[7], data[8]);

}


Matrix J2CyclicBoundingSurface3D::tangent(3, 3);


// full constructor
J2CyclicBoundingSurface3D::J2CyclicBoundingSurface3D( int tag, double G, double K, double su, double rho, double h, double m, double h0, double chi, double beta)
    :J2CyclicBoundingSurface (tag, ND_TAG_J2CyclicBoundingSurface3D, G, K, su, rho, h, m, h0, chi, beta)
{
}

// null constructor
J2CyclicBoundingSurface3D::J2CyclicBoundingSurface3D()
  : J2CyclicBoundingSurface()
{  
}

// destructor
J2CyclicBoundingSurface3D::~J2CyclicBoundingSurface3D()
{ 
} 

// make a clone of this material
NDMaterial* 
J2CyclicBoundingSurface3D::getCopy()
{
    J2CyclicBoundingSurface3D *clone = new J2CyclicBoundingSurface3D(this->getTag(), 
                        m_shear, m_bulk, m_su, m_density, m_h_par, m_m_par, m_h0_par, m_chi, m_beta);
    // *clone = *this;
    return clone;
}

// send back type of material
const char* 
J2CyclicBoundingSurface3D::getType() const
{
    return "ThreeDimensional";
}

// send back order of strain
int 
J2CyclicBoundingSurface3D::getOrder() const
{ 
    return 6; 
} 

// get the strain and integrate plasticity equations
int 
J2CyclicBoundingSurface3D::setTrialStrain(const Vector &strain_from_element)
{
    for (int i = 0; i < 6; ++i)
        m_strain_np1(i) = strain_from_element(i);
    // m_strain_np1 = strain_from_element;
    this->integrate();

    return 0;
}

// unused trial strain functions
int 
J2CyclicBoundingSurface3D::setTrialStrain(const Vector &v, const Vector &r)
{
    m_strainRate_n1 = r;
    m_strain_np1 = v;
    this->integrate();

    return 0;
}

// send back the strain
const Vector& 
J2CyclicBoundingSurface3D::getStrain()
{
    return_vector.setData(m_strain_np1);
    return return_vector;
} 

// send back the stress 
const Vector& 
J2CyclicBoundingSurface3D::getStress()
{
    //return m_stress_np1;
    // return m_stress_t_n1;
    return_vector.setData(m_stress_t_n1);
    return return_vector;
}

// send back the tangent 
const Matrix& 
J2CyclicBoundingSurface3D::getTangent()
{
    return calcTangent();
} 

// send back the tangent 
const Matrix& 
J2CyclicBoundingSurface3D::getInitialTangent()
{
    return m_Ce;
}