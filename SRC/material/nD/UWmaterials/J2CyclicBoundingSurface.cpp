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
** ****************************************************************** */

// Written: Diego Turello(*), Alborz Ghofrani and Pedro Arduino
//            Sep 2017, University of Washington
//          (*) Universidad Nacional de C�rdoba, FCEFyN. Depto Estructuras.
//              Universidad Tecnol�gica Nacional, GIMNI.
//              CONICET
// 
// Description: This file contains the implementation for the Borja material class.
// MultiaxialCyclicPlasticity for clays 
//
// Borja R.I, Amies, A.P.Multiaxial Cyclic Plasticity Model for Clays,
//   ASCE J.Geotech.Eng.Vol 120, No 6, 1051 - 1070
//            

#include <J2CyclicBoundingSurface.h>
#include <J2CyclicBoundingSurface3D.h>
#include <J2CyclicBoundingSurfacePlaneStrain.h>
#include <Information.h>
#include <Parameter.h>
#include <string.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <plastic/Voigt.hpp>

using namespace OpenSees;

//parameters

char  unsigned      J2CyclicBoundingSurface::m_ElastFlag = 1;  // --> default visco-elasto -plastic   
//char  unsigned      J2CyclicBoundingSurface::m_ElastFlag = 2; // default visco-elastic




//null constructor
J2CyclicBoundingSurface::J2CyclicBoundingSurface() :
    NDMaterial()
{

}


// full constructor
J2CyclicBoundingSurface::J2CyclicBoundingSurface(int  tag, int classTag,
    double G,
    double K,
    double su,
    double rho,
    double h,
    double m,
    double h0,
    double chi,
    double beta)
    :
    NDMaterial(tag, classTag),
    // m_sigma0_n(6), m_sigma0_np1(6), m_stress_n(6), m_stress_np1(6), 
    m_sigma0_n{}, m_sigma0_np1{}, m_stress_n{}, m_stress_np1{}, 
    m_strain_np1{}, m_strain_n{}, m_strainRate_n{}, m_strainRate_n1{}, 
    m_Cep(6, 6), m_Ce(6, 6), 
    m_D(6, 6), 
    // m_stress_vis_n(6), m_stress_vis_n1(6), m_stress_t_n1(6)
    m_stress_vis_n{}, m_stress_vis_n1{}, m_stress_t_n1{}
{
    double m_poiss = (3.*K - 2.*G) / 2. / (3. * K + G);

    if (m_poiss > 0.5) {
        opserr << "\n ERROR! J2CyclicBoundingSurface Poiss can not be grater than 0.50!" << "\n";
        exit(-1);
        return;
    }

    m_su = su;
    m_bulk = K; // E / 3 / (1-2*nu)
    m_shear = G;  // E / 2 / (1 + nu);
    m_R = sqrt(8. / 3.)*su;
    m_density = rho;
    m_kappa_inf = 1.0e10;
    m_h_par = h;
    m_m_par = m;
    m_h0_par = h0;
    m_beta = beta;
    m_kappa_n = m_kappa_inf;
    m_kappa_np1 = m_kappa_inf;
    m_psi_n = 2.*m_shear;
    m_chi = chi;

    m_isElast2Plast = false;
    debugFlag = false;
    small = 1.0e-10;

    calcInitialTangent();
}


// destructor
J2CyclicBoundingSurface::~J2CyclicBoundingSurface()
{

}


//--------------------Plasticity-------------------------------------

void 
J2CyclicBoundingSurface::zero()
{
}

void
J2CyclicBoundingSurface::integrate()
{
    if (m_ElastFlag == 0) // Force elastic response
        elastic_integrator();
    else if (m_ElastFlag == 1)  // ElastoPlastic response
        plastic_integrator();
    else if (m_ElastFlag == 2)
        viscoElastic_integrator();
}

void
J2CyclicBoundingSurface::elastic_integrator()
{
    VectorND<6> dStrain = m_strain_np1 - m_strain_n; //delta strain for the step
    m_stress_np1 = m_stress_n + m_Ce*dStrain;
    m_stress_t_n1 = m_stress_np1;
}

void
J2CyclicBoundingSurface::viscoElastic_integrator()
{
    // Update using strain rate from element
    //Vector dStrain = m_strain_np1 - m_strain_n; //delta strain for the step
    //m_stress_vis_n1 = m_D * (m_beta * m_strainRate_n1 + (1.0 - m_beta) * m_strainRate_n);
    //m_stress_np1 = m_stress_n + m_Ce * dStrain + m_D * ((1 - m_beta) * m_strainRate_n + m_beta * m_strainRate_n1);

    //Update using strain rate approx as dStrain/dt
    VectorND<6> dStrain = m_strain_np1 - m_strain_n; //delta strain for the step
    //m_stress_vis_n1 = m_D * (m_beta * m_strainRate_n1 + (1.0 - m_beta) * m_strainRate_n);
    double dt = 0.0;
    if (this->getTimeStep(dt) == 0) { 
        m_stress_vis_n1 = m_D * (dStrain) / dt;
    }
    else {
        m_stress_vis_n1 = m_stress_vis_n;
    }
    //m_stress_np1 = m_stress_n + m_Ce * dStrain + (m_stress_vis_n1 - m_stress_vis_n);
    //m_stress_t_n1 = m_stress_np1;
    m_stress_np1 = m_stress_n + m_Ce * dStrain;
    m_stress_t_n1 = m_stress_np1 + m_stress_vis_n1;
}


// plastic integration routine
void
J2CyclicBoundingSurface::plastic_integrator()
{  
    const double tol_rel = 1.0e-10;
    // static Vector eye(Voigt::ivol);
    VectorND<6> eye{};
    eye(0) = 1.0;
    eye(1) = 1.0;
    eye(2) = 1.0;

    // incremental strain for the step
    VectorND<6> dStrain = m_strain_np1 - m_strain_n;
    VectorND<6> dStrain_dev = Voigt::Dev(dStrain);    // incremental deviatoric strain
    double dStrain_vol = Voigt::Trace(dStrain);       // incremental volumetric strain



    m_sigma0_np1 = m_sigma0_n;

    // VectorND<6> dev_stress_n{};   // deviatoric stress
    VectorND<6> dev_stress_np1{};    // deviatoric stress
    // VectorND<6> dev_sigma0_np1{};    // deviatoric stress
    VectorND<6> dev_stress_n = Voigt::Dev(m_stress_n);
    VectorND<6> dev_sigma0_np1 = Voigt::Dev(m_sigma0_np1);

    double loadingCond = -1;
    double norm_dev_stress_n   = sqrt(Voigt::Dot(dev_stress_n,   dev_stress_n, 1));
    double norm_dev_sigma0_np1 = sqrt(Voigt::Dot(dev_sigma0_np1, dev_sigma0_np1, 1));

    m_kappa_np1 = m_kappa_n;
    m_psi_np1 = m_psi_n;

    m_stress_np1.zero();

    // check loading/unloading
#if 0
    // this is how it's written in the paper
    double temp_numerator = inner_product(-(1 + m_kappa_n) * dev_stress_n - m_kappa_n * (1 + m_kappa_n)*(dev_stress_n - dev_sigma0_np1), dStrain_dev, 3);
    double temp_denominator = inner_product(  (1 + m_kappa_n) * dev_stress_n - m_kappa_n * dev_sigma0_np1, (dev_stress_n - dev_sigma0_np1), 1 );

    if (abs(temp_denominator) < small)
        loadingCond = 0.0;
    else
        loadingCond = temp_numerator / temp_denominator;
#else
    // this is how we do it
    loadingCond = Voigt::Dot(m_kappa_n/(1 + m_kappa_n) * dev_sigma0_np1 - dev_stress_n, dStrain_dev, 3);
#endif

    if (loadingCond > 0.0)
    {
        m_sigma0_np1 = m_stress_n;
        dev_sigma0_np1 = Voigt::Dev(m_sigma0_np1);
        loadingCond = 0.0;
    }

    double H_n;
    double H_np1;
    if (vector_norm(dev_stress_n - dev_sigma0_np1, 1) == 0.0)
    {
        //
        // unloading (or beginning of loading) just happened
        //
        if (debugFlag)
            opserr << "Initial loading." << "\n";

        VectorND<6> devStrainDir{};
        double devStrainNorm = vector_norm(dStrain_dev, 2);
        if (fabs(devStrainNorm) < small) {
            m_stress_np1 = m_stress_n + m_bulk * dStrain_vol * eye + m_psi_np1 * convert_to_stressLike(dStrain_dev);
            return;
        }
        else
            devStrainDir = dStrain_dev / devStrainNorm;

        m_psi_np1 = 2 * m_shear;
        m_kappa_np1 = 1.0e10;

        H_np1 = H(m_kappa_np1);

        VectorND<2> res{}; 
        res(0) = m_psi_np1 * (1.0 + 3.0 * m_shear *  m_beta / H_np1) / (2.0 * m_shear) - 1.0;
        res(1) = vector_norm(dev_stress_n + (1.0 + m_kappa_np1) * m_psi_np1 * convert_to_stressLike(dStrain_dev), 1) / m_R - 1.0;

        double res_norm = res.norm(); //vector_norm(res, 3);

        // Initialize variables for the Newton
        int       iteration_counter = 0;
        constexpr static int max_iterations = 50;
        MatrixND<2,2>    Ktan{};

        VectorND<2> incVar{};
        double tol_material = tol_rel * res_norm;

        for (iteration_counter = 0; iteration_counter < max_iterations; iteration_counter++)
        {
            if (debugFlag) opserr << "iteration " << iteration_counter << " , norm = " << res_norm << "\n";

            if (res_norm < tol_material + small)
            {
                m_stress_np1 = m_stress_n + m_bulk*dStrain_vol * eye + m_psi_np1 * convert_to_stressLike(dStrain_dev);
                break;
            }

            VectorND<6> temp = (dev_stress_n + (1.0 + m_kappa_np1) * m_psi_np1 * convert_to_stressLike(dStrain_dev));
            temp /= vector_norm(temp, 1);
            Ktan(0, 0) = (1.0 + 3.0 * m_shear *  m_beta / H_np1) / (2.0 * m_shear);
            Ktan(0, 1) = (-3.0 * m_shear * m_psi_np1 * m_beta * m_m_par / m_h_par / pow(m_kappa_np1, m_m_par + 1.0)) / (2.0 * m_shear);
            Ktan(1, 0) = ((1 + m_kappa_np1) * inner_product(temp, dStrain_dev, 3)) / m_R;
            Ktan(1, 1) = (inner_product(temp, m_psi_np1 * convert_to_stressLike(dStrain_dev), 1)) / m_R;

            // Solve the system
            Ktan.rsolve(res, incVar);

            m_psi_np1   -= incVar(0);
            m_kappa_np1 -= incVar(1);

            H_np1 = H(m_kappa_np1);

            // calculate new residual
            res(0) = m_psi_np1 * (1.0 + 3.0 * m_shear * (+m_beta / H_np1))/(2.0 * m_shear) - 1.0;
            res(1) = vector_norm(dev_stress_n + (1.0 + m_kappa_np1) * m_psi_np1 * convert_to_stressLike(dStrain_dev), 1) / m_R - 1.0;

            res_norm = res.norm(); //vector_norm(res, 3);
        }
    }
    else
    {
        if (debugFlag) opserr << "Loading continues..." << "\n";

        VectorND<6> devStrainDir{};
        double devStrainNorm = vector_norm(dStrain_dev, 2);
        if (fabs(devStrainNorm) < small) {
            m_stress_np1 = m_stress_n + m_bulk * dStrain_vol * eye + m_psi_np1 * convert_to_stressLike(dStrain_dev);
            return;
        }

        // calculate the initial residual
        H_n = H(m_kappa_n);
        H_np1 = H(m_kappa_np1);

        VectorND<2> res{};
        res(0) = m_psi_np1 * (1.0 + 3.0 * m_shear * ((1 - m_beta) / H_n + m_beta / H_np1)) / (2.0 * m_shear) - 1.0;
        res(1) = vector_norm(dev_stress_n + (1.0 + m_kappa_np1) * m_psi_np1 * convert_to_stressLike(dStrain_dev) + m_kappa_np1 * (dev_stress_n - dev_sigma0_np1), 1) / m_R - 1.0;;

        double res_norm = res.norm(); //vector_norm(res, 3);

        // Initialize variables for the Newton
        int iteration_counter = 0;
        const int max_iterations = 50;
        MatrixND<2,2> Ktan{};

        VectorND<2> incVar{};
        double tol_material = tol_rel * res_norm;

        for (iteration_counter = 0; iteration_counter < max_iterations; iteration_counter++)
        {
            if (debugFlag) opserr << "iteration " << iteration_counter << " , norm = " << res_norm << "\n";

            if (res_norm < tol_material + small)
            {
                m_stress_np1 = m_stress_n + m_bulk * dStrain_vol * eye + m_psi_np1 * convert_to_stressLike(dStrain_dev);
                break;
            }

            VectorND<6> temp = (dev_stress_n + (1.0 + m_kappa_np1) * m_psi_np1 * convert_to_stressLike(dStrain_dev) + m_kappa_np1 * (dev_stress_n - dev_sigma0_np1));
            temp = temp / vector_norm(temp, 1);
            Ktan(0, 0) = (1.0 + 3.0 * m_shear * ((1 - m_beta) / H_n + m_beta / H_np1)) / (2.0 * m_shear);
            Ktan(0, 1) = (-3.0 * m_shear * m_psi_np1 * m_beta * m_m_par / m_h_par / pow(m_kappa_np1, m_m_par + 1.0)) / (2.0 * m_shear);
            Ktan(1, 0) = ((1 + m_kappa_np1) * inner_product(temp, dStrain_dev, 3)) / m_R;
            Ktan(1, 1) = (inner_product(temp, dev_stress_n + m_psi_np1 * convert_to_stressLike(dStrain_dev) - dev_sigma0_np1, 1)) / m_R;

            // Solve the system
            Ktan.rsolve(res, incVar);

            m_psi_np1   = m_psi_np1 - incVar(0);
            m_kappa_np1 = m_kappa_np1 - incVar(1);

            H_np1 = H(m_kappa_np1);

            // calculate new residual
            res(0) = m_psi_np1 * (1.0 + 3.0 * m_shear * ((1 - m_beta) / H_n + m_beta / H_np1)) / (2.0 * m_shear) - 1.0;
            res(1) = vector_norm(dev_stress_n + (1.0 + m_kappa_np1) * m_psi_np1 * convert_to_stressLike(dStrain_dev) + m_kappa_np1 * (dev_stress_n - dev_sigma0_np1), 1) / m_R - 1.0;

            res_norm = res.norm(); //vector_norm(res, 3);
        }
    }

    //m_stress_vis_n1 = m_D * (m_beta * m_strainRate_n1 + (1.0 - m_beta) * m_strainRate_n);
    double dt = 0.0;
    if (this->getTimeStep(dt) == 0) {
        m_stress_vis_n1 = m_D * (dStrain) / dt; 
    }
    else {
        m_stress_vis_n1 = m_stress_vis_n; 
    }

    // m_stress_np1 += (m_stress_vis_n1 - m_stress_vis_n);  //stress update using incremental form

    m_stress_t_n1 = m_stress_np1 + m_stress_vis_n1;

    return;
}


// Trace Operator
double 
J2CyclicBoundingSurface::trace(const VectorND<6>& V)
{
    return V(0) + V(1) + V(2);
}


// Inner product
double 
J2CyclicBoundingSurface::inner_product(const VectorND<6> &x, const VectorND<6> &y, int type)
{
    double modifier = 1.0;
    double inner = 0.0;
    switch (type)
    {
    case 1: // stress : stress
        modifier = 2.0;
        break;
    case 2: // strain : strain
        modifier = 0.5;
        break;
    case 3: // stress : strain
        modifier = 1.0;
        break;
    default:
        break;
    }

    for (int i = 0; i < x.size(); i++)
        inner += x(i) * y(i) + (i > 2) * (modifier - 1.0) * x(i) * y(i);

    return inner;
}

// Norm of a vector
double
J2CyclicBoundingSurface::vector_norm(const VectorND<6>& x, int type)
{
    return sqrt(inner_product(x, x, type));
}

VectorND<6>
J2CyclicBoundingSurface::convert_to_stressLike(const VectorND<6>& v)
{
    VectorND<6> res = v;
    for (int ii = 3; ii < 6; ii++)
        res(ii) *= 0.5;

    return res;
}

const Matrix &
J2CyclicBoundingSurface::getDampTangent()
{
    return m_D;
}


// set up for initial elastic
void
J2CyclicBoundingSurface::calcInitialTangent()
{
    Matrix I2xI2(6, 6), I4dev(6, 6), eye(6, 6);
    for (unsigned int i = 0; i < 3; i++)
        for (unsigned int j = 0; j < 3; j++)
        {
            I2xI2(i, j) = 1;
        }

    eye(0, 0) = 1.0;
    eye(1, 1) = 1.0;
    eye(2, 2) = 1.0;
    eye(3, 3) = 1.0;
    eye(4, 4) = 1.0;
    eye(5, 5) = 1.0;

    I4dev = eye - 1 / 3 * I2xI2;

    //m_Ce = m_bulk * I2xI2 + 2 * m_shear*I4dev;
    m_Ce = m_bulk * I2xI2 + m_shear * I4dev;
    m_D  = m_chi * m_Ce;

    return;
}


// hardening function
double 
J2CyclicBoundingSurface::H(double kappa)
{
    if (kappa < 0) 
      return 1.0e-10;

    return m_h_par * pow(kappa, m_m_par);
}




NDMaterial*
J2CyclicBoundingSurface::getCopy(const char* type)
{
    if (strcmp(type, "ThreeDimensional") == 0 || strcmp(type, "3D") == 0) {
        J2CyclicBoundingSurface3D* clone;
        clone = new J2CyclicBoundingSurface3D(this->getTag(), m_shear, m_bulk, m_su, m_density, m_h_par, m_m_par, m_h0_par, m_chi, m_beta);
        return clone;
    }
    else if (strcmp(type, "PlaneStrain2D") == 0 || strcmp(type, "PlaneStrain") == 0) {
        J2CyclicBoundingSurfacePlaneStrain* clone;
        clone = new J2CyclicBoundingSurfacePlaneStrain(this->getTag(), m_shear, m_bulk, m_su, m_density, m_h_par, m_m_par, m_h0_par, m_chi, m_beta);
        return clone;
    }
    else {
        return this->NDMaterial::getCopy(type);
    }
}



int
J2CyclicBoundingSurface::commitState()
{
    m_sigma0_n = m_sigma0_np1;
    m_stress_n = m_stress_np1;
    m_kappa_n  = m_kappa_np1;
    m_psi_n    = m_psi_np1;
    m_strain_n = m_strain_np1;
    m_strainRate_n = m_strainRate_n1;
    m_stress_vis_n = m_stress_vis_n1;

    return 0;
}

int
J2CyclicBoundingSurface::revertToLastCommit()
{
    return 0;
}


int
J2CyclicBoundingSurface::revertToStart() {

    // added: C.McGann, U.Washington for InitialStateAnalysis
    if (ops_InitialStateAnalysis) {
        // do nothing, keep state variables from last step
    }
    else {
        // normal call for revertToStart (not initialStateAnalysis)
        this->zero();
    }

    return 0;
}

int
J2CyclicBoundingSurface::setParameter(const char **argv, int argc,
    Parameter &param)
{
    int theMaterialTag;
    theMaterialTag = atoi(argv[1]);

    if (theMaterialTag == this->getTag()) {

        if (strcmp(argv[0], "updateMaterialStage") == 0) {     // enforce elastic/elastoplastic response
            return param.addObject(1, this);
        }
        else if (strcmp(argv[0], "materialState") == 0) {     // enforce elastic/elastoplastic response
            return param.addObject(2, this);
        }
    }
    return -1;
}

int
J2CyclicBoundingSurface::updateParameter(int responseID, Information &info)
{
    // called updateMaterialStage
    if (responseID == 1) {
        m_ElastFlag = info.theInt;
        m_isElast2Plast = true;
        return 0;
    }
    // called materialState in tcl file
    else if (responseID == 2) {
        m_ElastFlag = (int)info.theDouble;
        m_isElast2Plast = true;
        return 0;
    }
    return -1;
}

int
J2CyclicBoundingSurface::activateParameter(int paramID)
{
    // TODO : implement this
    return -1;
}

int
J2CyclicBoundingSurface::sendSelf(int commitTag, Channel &theChannel)
{
    // TODO : implement this
    return 0;
}

int
J2CyclicBoundingSurface::recvSelf(int commitTag, Channel &,
    FEM_ObjectBroker &)
{
    // TODO : implement this
    return -1;
}


// send back the tangent 
const Matrix&
J2CyclicBoundingSurface::calcTangent()
{
    // Force elastic response
    if (m_ElastFlag == 0) {
        return    m_Ce;
    }
    // ElastoPlastic response
    else if (m_ElastFlag == 1) {
        // only the symetric part
        Matrix I2xI2(6, 6), I4dev(6, 6), eye(6, 6);
        for (unsigned int i = 0; i < 3; i++)
            for (unsigned int j = 0; j < 3; j++)
            {
                I2xI2(i, j) = 1;
            }

        for (unsigned int i = 0; i < 6; i++)
            for (unsigned int j = 0; j < 6; j++)
            {
                if (i == j) { eye(i, j) = 1; }
            }
        I4dev = eye - 1 / 3 * I2xI2;

        m_Cep = m_bulk * I2xI2 + 0.5 * m_psi_np1 * I4dev;
        
        double dt = 0.0;
        if (this->getTimeStep(dt) == 0) {
            m_Cep += (1.0 / dt) * m_D; // Adding viscous damping contribution to Cep
        }

        return m_Cep;
        //return m_Ce;
    }
    // ViscoElastic response
    else if (m_ElastFlag == 2) {

        m_Cep = m_Ce;
        double dt = 0.0;   
        if (this->getTimeStep(dt) == 0) {
            m_Cep += (1.0 / dt) * m_D; // Adding viscous damping contribution to Cep
        }
        return m_Cep;
    }
    else {
        opserr << "\n ERROR! J2CyclicBoundingSurface m_ElastFlag not valid - returning Ce" << "\n";
        return    m_Ce;
    }

}

// send back the tangent 
//const Matrix&
//J2CyclicBoundingSurface::calcInitialTangent()
//{
//    return m_Ce;
//}


void
J2CyclicBoundingSurface::Print(OPS_Stream & s, int flag)
{
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << OPS_PRINT_JSON_MATE_INDENT << "{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"" << this->getClassType() << "\", ";
        s << "\"G\": " << m_shear << ", ";
        s << "\"K\": " << m_bulk << ", ";
        s << "\"Su\": " << m_su << ", ";
        s << "\"density\": " << m_density << ", ";
        s << "\"H\": " << m_h_par << ", ";
        s << "\"m\": " << m_m_par << ", ";
        s << "\"Ho\": " << m_h0_par << ", ";
        s << "\"eta\": " << m_chi << ", ";
        s << "\"beta\": " << m_beta;
        s << "}";
    }
}

