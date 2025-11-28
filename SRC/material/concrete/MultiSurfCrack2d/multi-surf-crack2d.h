#pragma once
#include <string>
#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <algorithm>
#include <cstdlib>

#include "Vector.h"
#include "Matrix.h"
#include <ID.h> 
#include <NDMaterial.h>

class MultiSurfCrack2D: public NDMaterial
{
public:

    // FUNCTIONS REQUIRED BY OTHER OPENSEES ROUTINES
    // ---------------------------------------------

	// empty constructor
    MultiSurfCrack2D();

	// full constructor
    MultiSurfCrack2D(int tag,
          Matrix De,                // elastic loading stiffness matrix (MPa/mm)
          Matrix DeU,               // elastic unloading stiffness matrix (MPa/mm)
          double user_fc,           // concrete strength (MPa)
          double user_ag,           // max aggregate diameter (mm)
          double user_fcl,          // crack closure stress magnitude (MPa)
          double user_Acr,          // nominal crack area (thickness * element discretization length)
          double user_rho_lok,      // interlock dilation parameter
          double user_chi_lok,      // cohesion ratio for interlock surface
          double user_rho_act,      // rubble dilation parameter
          double user_mu,           // unloading friction coefficient
          double user_chi,          // unloading cohesion ratio (c/vcimax)  -- activation surface
          double user_zeta,         // slip break point (pinching)
          double user_kappa,        // roughness transition parameter
          double user_theta,        // average contact angle (radians)
          double user_w,            // initial crack width (mm)
          int cPath);              // 1 if load path is width-constant, else 0

    //Crack(const Crack&);
    
    // destructor
	~MultiSurfCrack2D();

	// clone the material
	NDMaterial* getCopy();
    NDMaterial* getCopy(const char* type);

    int getOrder() const { return 2; }
    const char* getType() const { return "MultiSurfCrack2D"; }

    int commitState();
    int revertToLastCommit();
    int revertToStart();

    //get the strain 
    int setTrialStrain(const Vector& strainFromElement);

    //send back the strain
    const Vector& getStrain();

    //send back the stress 
    const Vector& getStress();

    //send back the tangent 
    const Matrix& getTangent();
    const Matrix& getInitialTangent();

    void Print(OPS_Stream& s, int flag);
    int sendSelf(int commitTag, Channel& theChannel);
    int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);
	

    void PrintState();

private:

	// INTERNAL "HELPER" FUNCTIONS
    // --------------------------------------------
    void setLoadingDirection();
    void setDamagedStiffness();
    void trialTangentStiffness();
    void takeElasticTrialStep();
    void checkActiveYieldSurface();
    void checkActiveYieldSurface(std::string& controlling_surface);
    int overlapYieldSurface();
    double getF_interlock();
    double getF_activation();
    void plasticReturn();
    double incrementPlasticMultiplier();
    void trialStressReturn();
    void updateFlowDir();
    void updateInterlockFlowDir();
    void updateActivationFlowDir();
    void initializeSurfacesShape();
    void updateInterlockShape();
    void updateDerivedShapeParameters();
    double get_w_tilde();
    void setvcimax();                                       // set vcimax, and drag rest of yield surfaces along with MCFT surface
    double d_vcimax_d_w_tilde();
    double d_vcimax_dnu();
    double d_F_dvcimax();
    double d_F_dfcc();
    double d_F_dw();
    double d_F_dnu();                                       // partial of yield surface w/rsp crush strain
    double d_F_dchi();
    double d_F_dk();
    Vector strain_gradient_w_tilde();
    Vector strain_gradient_nu();
    Vector strain_gradient_chi();
    Vector strain_gradient_k();
    double d_r_ds();
    double get_postPeakCompressionSoftening();


    // CLASS MEMBERS
    // ---------------------------------------------
    // input elastic properties (initial) ...............
    const Matrix tangent_elastic_loading_initial;
    const Matrix tangent_elastic_unloading_initial;

//public:

    // committed state ...................................
    Vector strain_committed;
    Vector strain_plastic_committed;
    Vector stress_committed;
    Matrix tangent_committed;
    double sMAX;                            // maximum (positive) plastic slip
    double sMIN;                            // maximum (negative) plastic slip
    double sDAM;                            // maximum (pos or neg) total slip -- used for sliding damage calculations
    double wMAX;

    // current state .....................................
    Matrix tangent_elastic_damaged;
    Matrix tangent_current;                 // Stress/length, for internal calculations
    Matrix tangentF;                        // Force/length, for handing off to element
    Vector strain_increment;
    Vector strain_current;                  // strain from element (not counting initial residual strain)
    Vector strain_plastic_current;
    Vector stress_current;                  // crack stresses used for internal calculations
    Vector force_current;                   // crack forces, passed to zeroLengthND element
    double vcimax;                          // peak attainable shear stress (positive valued)
    double crushing_strain;                 // accumulated crushing strain (crack closure, therefore negative) while F_c^lok is active
    double dLambda;

    // crack properties ..................................
    double fc;                              // concrete strength (MPa)
    double ag;                              // max aggregate diameter (mm)
    double fcl;                             // closure stress magnitude (MPa)
    double Acr;                             // area of this crack section (thickness x discretization length)
    double w_res;                           // residual crack width at start of analysis (mm)

    // MCFT aggregate interlock (lok) parabola ...........
    double c_lok;                           // cohesion
    double chi_lok;                         // cohesion ratio
    double chi_lok_r;                       // bounding value of chi_lok for rough surface
    double chi_lok_s;                       // bounding value of chi_lok for smooth surface
    double d;                               // parameters defining yield surface parabola
    double h;                               //   --
    double k;                               //   --
    double k_rough;                         // bounding value of k for rough surface
    double k_smooth;                        // bounding value of k for smooth surface
    double rho_lok;                         // dilation parameter for aggregate interlock
    Vector dF_lok;                          // flow direction for aggregate interlock
    Vector dG_lok;                          // non-associated flow direction for interlock
    double kappa;                           // roughness factor for interlock
    double theta_tilde;                     // average contact angle (radians) for sliding damage calculation
    int criticalPath;                       // 1 if loading path is constant-width, else 0
    double shift;                           // amount of tensile strength available for constant-width surface (i.e., horizontal shift from origin)

    // partial Elliptical tension yield surface.............
    double ft;                              // peak traction
    double chi_t;                           // ratio of ft / vcimax
    double xc;                              // horizontal position of tensile ellipse center
    double chi_c;                           // ratio of xc / vcimax
    double a_t;                             // half horizontal axis length for tensile ellipse
    double b_t;                             // half vertical axis length for tensile ellipse
    double m;                               // parameter enforcing size of tensile ellipse (for real ellipse, need 0 < m < 1)
    double l;                               // survival algebra parameter
    double n;                               // survival algebra parameter

    // semi-Elliptical compression yield surface............
    double fcc;                             // fc reduced to account for crushing damage (positive-valued)
    double a_c;                             // half horizontal axis length for compression ellipse
    double b_c;                             // half vertical axis length for compression ellipse

    // Hyperbolic unloading activation (act) surface .....
    double c_act;                           // cohesion for activtion branch (may differ from c_lok)
    double chi_act;                         // cohesion ratio for activation branch (may differ from chi_lok)
    double mu;                              // friction ratio
    double rho_act;                         // dilation parameter for frictional unloading
    Vector dF_act;                          // flow direction for frictional unloading
    Vector dG_act;                          // non-associated flow direction for unloading
    double zeta;                            // break point (aka pinching point)

    // crack stage .......................................
    bool isElastic;
    std::string dir;
    double yield_criterion;
    std::string active_surface;
    double last_vci_lok;                    // shear value of most recent stress point that converged to the combined interlock surface
    double* rho_engaged;                    // pointers to the currently engaged yield surface
    Vector* dF_engaged;                     //  --
    Vector* dG_engaged;                     //  --

    // opensees printout
    bool verbose = false;
};