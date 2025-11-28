#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <algorithm>
#include <cstdlib>

#include <MultiSurfCrack2D.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <MaterialResponse.h>
#include <elementAPI.h>


// Parse user inputs and create material
void* OPS_MultiSurfCrack2D()
{
	int argc = OPS_GetNumRemainingInputArgs();
	if (argc < 22) {
		opserr << "WARNING wrong number of arguments\n";
		opserr << "Want: nDMaterial MultiSurfCrack2D tag? E? H? M? K? Eunl? Hunl? Munl? Kunl? fc? ag? fcl? Acr? rho_lok? chi_lok? rho_act? mu? chi_act? zeta? kappa? theta? w? <cPath?>" << endln;
		return 0;
	}

	int tag;
	int numdata = 1;
	if (OPS_GetIntInput(&numdata, &tag) < 0) {
		opserr << "WARNING invalid nDMaterial MultiSurfCrack2D tag" << endln;
		return 0;
	}

	double data[21];
	numdata = 21;
	if (OPS_GetDoubleInput(&numdata, data) < 0) {
		opserr << "WARNING invalid double inputs" << endln;
		opserr << "MultiSurfCrack2D: " << tag << endln;
		return 0;
	}

	int path = 0;
	if (argc == 23)
	{
		numdata = 1;
		if (OPS_GetIntInput(&numdata, &path) < 0) {
			opserr << "WARNING invalid int input" << endln;
			opserr << "MultiSurfCrack2D: " << tag << endln;
			return 0;
		}
	}

	Matrix De_user(2, 2);
	De_user(0, 0) = data[0];
	De_user(0, 1) = data[1];
	De_user(1, 0) = data[2];
	De_user(1, 1) = data[3];
	Matrix DeU_user(2, 2);
	DeU_user(0, 0) = data[4];
	DeU_user(0, 1) = data[5];
	DeU_user(1, 0) = data[6];
	DeU_user(1, 1) = data[7];

	return new MultiSurfCrack2D(tag, De_user, DeU_user, data[8], data[9], data[10], data[11], data[12], data[13], data[14], data[15], data[16], data[17], data[18], data[19], data[20], path);
}


MultiSurfCrack2D::MultiSurfCrack2D(int tag,
			 Matrix De,
			 Matrix DeU,
			 double user_fc,
			 double user_ag,
			 double user_fcl,
			 double user_Acr,
			 double user_rho_lok,
			 double user_chi_lok,
			 double user_rho_act,
			 double user_mu,
			 double user_chi,
			 double user_zeta,
			 double user_kappa,
			 double user_theta,
			 double user_w,
			 int cPath = 0)
  : NDMaterial(tag, ND_TAG_MultiSurfCrack2D),
	tangent_elastic_loading_initial(De), tangent_elastic_unloading_initial(DeU),
	strain_committed(2),strain_plastic_committed(2),stress_committed(2),tangent_committed(De),
	sMAX(0.001),sMIN(-0.001),sDAM(0.),wMAX(0.1),
	tangent_elastic_damaged(De),tangent_current(De),tangentF(2,2),
	strain_increment(2),strain_current(2),strain_plastic_current(2),stress_current(2),force_current(2),
	vcimax(15.),crushing_strain(0.),dLambda(0.),fc(user_fc),ag(user_ag),fcl(user_fcl),Acr(user_Acr),w_res(user_w),
	c_lok(0.),chi_lok(user_chi_lok),chi_lok_r(user_chi_lok),chi_lok_s(user_chi),h(-1.),k(1.),k_rough(1.),k_smooth(0.6),
	rho_lok(user_rho_lok),dF_lok(2),dG_lok(2),kappa(user_kappa),theta_tilde(user_theta),criticalPath(cPath),shift(0.005),m(0.5),fcc(user_fc),
	c_act(0.),chi_act(user_chi),mu(user_mu),rho_act(user_rho_act),
	dF_act(2),dG_act(2),zeta(user_zeta),isElastic(true),dir("loading"),
	yield_criterion(-1.),active_surface(""),last_vci_lok(0.0),
	rho_engaged(&rho_lok),dF_engaged(&dF_lok),dG_engaged(&dG_lok)
{

	// some error handling to ensure that combined yield surface is smooth
	initializeSurfacesShape();

	// set residual crack width (Maybe set residual stresses instead??)
	// 
	//    ! residual crack can be treated like a prestrain !
	//strain_current[0] = user_w;			// should be strain_committed so that getCopy() sends this to the element?
	//strain_committed[0] = user_w;
	wMAX = std::max(wMAX, user_w);

	// set initial size of MCFT parabola, cohesion and compression/tension ellipses
	setvcimax();

	// set residual crack stresses
	// ToDo...

	// initialize flow directions
	updateFlowDir();
}


MultiSurfCrack2D::MultiSurfCrack2D(): 
	NDMaterial(0, ND_TAG_MultiSurfCrack2D),
	tangent_elastic_loading_initial(2,2), tangent_elastic_unloading_initial(2,2),
	strain_committed(2), strain_plastic_committed(2), stress_committed(2),tangent_committed(2,2),
	sMAX(0.001),sMIN(-0.001),sDAM(0.),wMAX(0.1),
	tangent_elastic_damaged(2,2),tangent_current(2,2),tangentF(2,2),
	strain_increment(2),strain_current(2),strain_plastic_current(2),stress_current(2),force_current(2),
	vcimax(15.),crushing_strain(0.),dLambda(0.),fc(30.),ag(16.),
	c_lok(0.),chi_lok(0.18),chi_lok_r(0.18),chi_lok_s(0.0),h(-1.),k(1.),k_rough(1.),k_smooth(0.6),
	rho_lok(1.0),dF_lok(2),dG_lok(2),kappa(0.0),theta_tilde(0.0),criticalPath(0),shift(0.005),m(0.5),fcc(30.),
	c_act(0.),chi_act(0.18),mu(0.3),
	rho_act(1.0),dF_act(2),dG_act(2),zeta(0.5),isElastic(true),dir("loading"),
	yield_criterion(-1.),active_surface(""),last_vci_lok(0.0),
	rho_engaged(&rho_lok),dF_engaged(&dF_lok),dG_engaged(&dG_lok)
{
	// do nothing
}


// ToDo: inherit from NDMaterial
//Crack::Crack(const Crack& mat)
//{
//
//}



MultiSurfCrack2D::~MultiSurfCrack2D()
{

}


void MultiSurfCrack2D::initializeSurfacesShape()
{
	// activation cohesion point not to exceed interlock cohesion point
	chi_act = std::min(chi_lok, std::max(chi_act, 0.0001));

	// initialize "survival algebra variables"
	l = m * m / (8 * (1 - m));
	n = m / 4;

	// calculate derived shape parameters
	updateDerivedShapeParameters();
}


// solve a trial displacement step
int MultiSurfCrack2D::setTrialStrain(const Vector& strainFromElement)
{
	//opserr << "MultiSurfCrack2D::setTrialStrain" << endln;
	//opserr << "  strainFromElement: [" << strainFromElement[0] << "," << strainFromElement[1] << "]" << endln;

	strain_current = strainFromElement;
	strain_increment = strain_current - strain_committed;

	// updating size of yield surface due to trial width increment
	setvcimax();							

	// decide if crack is loading or unloading
	setLoadingDirection();

	// update elastic stiffness terms (loading or unloading, damage)
	setDamagedStiffness();


	// GET CRACK TRIAL STRESSES FOR GIVEN STRAIN
	// ============================================================================
	// =	ELASTIC TRIAL --> F < 0? ---> YES ---------------------------->  END  =
	// =                       |                                       ^          =
	// =                       v                                       |          =
	// =                        --------> NO --> PLASTIC CORRECTION ---           =
	// ============================================================================

	// get elastic trial stresses and check yield function
	takeElasticTrialStep();

	// plastic return process if elastic stresses violate yield surface
	if (yield_criterion > 0)
		plasticReturn();

	// check for surface overlap after return to activation surface (if relevant)
	int overlap = overlapYieldSurface();
	if (overlap)
	{
		// end the activation process and reset loading stiffness
		dir = "loading";
		setDamagedStiffness();

		// return to previous converged state and redo elastic-trial, plastic correction
		takeElasticTrialStep();
		if (yield_criterion > 0)
			plasticReturn();
	}

	// tangent stiffness at trial stress that has been returned to yield surface
	trialTangentStiffness();

	return 0;
}


void MultiSurfCrack2D::setLoadingDirection()
{	
	// check direction of slip increment
	double s = strain_current[1];
	double vci = stress_current[1];
	double s_prev = strain_committed[1];
	double ds = strain_increment[1];
	bool increasingSlipMagnitude = ((std::abs(s) - std::abs(s_prev)) > 0);
	

	if (increasingSlipMagnitude)
	{	
		// check if pinching point is exceeded
		if ((s > zeta * sMAX) || (s < zeta * sMIN))
			dir = "loading";

		else
		{
			dir = "unloading";
			//if (vci * last_vci_lok > 0)
			//{
			//	// single-sided cycles, in which case no activation reversal
			//	dir = "loading";
			//}
			//else
			//{
			//	// pinching point not exceeded
			//	dir = "unloading";
			//}
		}	
	}
	else
	{
		if (ds == 0)
		{
			// case of pure mode-I loading
			dir = "loading";
		}
		else
		{
			// slding reversal
			dir = "unloading";
		}
	}
}


void MultiSurfCrack2D::setDamagedStiffness()
{
	// crack stiffness in MPa/mm

	double& sigTR = stress_committed[0];									// current normal stress (note, still haven't taken trial step)
	double wCURR = strain_current[0] + w_res;								// current crack width
	double& wInc = strain_increment[0];										// crack opening direction
	double E_penalty = 1000.;												// penalty stiffness against crack overlap
	double E_ruggiero;


	// Update Crack Normal Stiffness -- Ruggiero model (simplified)
	if (wCURR >= 0.0)		// separated crack nodes
	{
		E_ruggiero = std::min(E_penalty, fcl / wMAX);						// ruggiero model with q = infty
	}
	else					// overlapping crack nodes
	{
		if (wInc > 0)
		{
			// reload to closure point
			double cohesive_stress = 0.0;									// ToDo: update cohesion as fxn of wMAX
			E_ruggiero = (-fcl + cohesive_stress - sigTR) / (-wCURR);
		}
		else
		{
			// penalty against further overlap
			E_ruggiero = E_penalty;
		}
	}

	tangent_elastic_damaged(0, 0) = E_ruggiero;
	//tangent_elastic_damaged(0, 0) = std::min(E_ruggiero, 1.0);


	// Update Crack Shear Stiffness 

	//		-- User Input w/ sliding damage
	double shearDamage = std::max(0.1, (1 - 2 * sMAX / ag));


	if (dir == "loading")
	{
		tangent_elastic_damaged(1, 1) = tangent_elastic_loading_initial(1, 1) * shearDamage;
		tangent_elastic_damaged(0, 1) = tangent_elastic_loading_initial(0, 1);
		tangent_elastic_damaged(1, 0) = tangent_elastic_loading_initial(1, 0);
	}
	else
	{						// undamaged stiffness if unloading
		tangent_elastic_damaged(1, 1) = tangent_elastic_unloading_initial(1, 1);
		tangent_elastic_damaged(0, 1) = tangent_elastic_unloading_initial(0, 1);
		tangent_elastic_damaged(1, 0) = tangent_elastic_unloading_initial(1, 0);
	}
}



void MultiSurfCrack2D::takeElasticTrialStep()
{
	// clear any plastic residual from previous non-solution-path steps
	isElastic = true;
	strain_plastic_current = strain_plastic_committed;
	stress_current = stress_committed;
	dLambda = 0.0;
	active_surface = "";

	double& dwe = strain_increment[0];
	double& dse = strain_increment[1];	
	double& s = strain_current[1];
	double sign_s = (s > 0) ? 1 : ((s < 0) ? -1 : 0);

	Matrix& De = tangent_elastic_damaged;
	double dsig = De(0, 0) * dwe + De(0, 1) * dse * sign_s;
	double dvci = De(1, 0) * dwe * sign_s + De(1, 1) * dse;
	stress_current[0] = stress_committed[0] + dsig;
	stress_current[1] = stress_committed[1] + dvci;

	checkActiveYieldSurface();
}


void MultiSurfCrack2D::checkActiveYieldSurface()
{
	double& ds = strain_increment[1];
	double& vci = stress_current[1];

	if (dir == "loading")
	{
		// check if interlock surface breached
		yield_criterion = getF_interlock();
		if (yield_criterion > 0)
		{
			isElastic = false;
			active_surface = "interlock";
			rho_engaged = &rho_lok;

			updateInterlockFlowDir();
			dF_engaged = &dF_lok;
			dG_engaged = &dG_lok;
		}
	}
	// else if (dir == "unloading" && (ds * vci > 0 && vci * last_vci_lok < 0))
	else if ((dir == "unloading") && (ds * vci > 0))
	{
		// check if activation surface breached
		yield_criterion = getF_activation();
		if (yield_criterion > 0)
		{
			isElastic = false;
			active_surface = "activation";
			rho_engaged = &rho_act;

			updateActivationFlowDir();
			dF_engaged = &dF_act;
			dG_engaged = &dG_act;
		}

		// check possible case that activation path crosses tensile branch of interlock surface
		// ....after activation surface has been returned
	}
	else
	{
		yield_criterion = -1.;
	}
}


void MultiSurfCrack2D::checkActiveYieldSurface(std::string& controlling_surface)
{
	if (controlling_surface == "interlock")
		yield_criterion = getF_interlock();
	else
		yield_criterion = getF_activation();
}


int MultiSurfCrack2D::overlapYieldSurface()
{
	int overlap = 0;
	double& sigTR = stress_current[0];

	// check if activation surface overlaps with tensile de-cohesion surface
	if ((active_surface == "activation") && (sigTR > 0))
	{
		// check if interlock surface breached for current stress state, which has already been returned to the activation surface
		yield_criterion = getF_interlock();

		// signal if breached
		if (yield_criterion > 0)
		{
			overlap = 1;
		}
	}
	else 
	{
		overlap = 0;
	}

	return overlap;
}


double MultiSurfCrack2D::getF_interlock()
{
	double& sigTR = stress_current[0];
	double& vciTR = stress_current[1];

	double F_interlock;
	
	if (!criticalPath)
	{
		// yield surface in tensile region
		if (sigTR > 0)
		{
			//F_interlock = std::pow(sigTR-xc, 2) + 2*chi_c*d/(h*chi_lok) * std::pow(vciTR, 2) - std::pow(k*(chi_t-chi_c)*vcimax, 2);

			double hkxv = h * k * chi_lok * vcimax;
			F_interlock = l / 2 * std::pow(h * vciTR, 2) + std::pow((k - chi_lok) * sigTR - l * hkxv, 2) - std::pow((n + l) * hkxv, 2);
		}

		else if (sigTR < h * vcimax)
		{
			// yield surface in compressive region
			F_interlock = std::pow((sigTR - h * vcimax) / a_c, 2) + std::pow(vciTR / b_c, 2) - 1;
		}

		else
		{
			// Parabolic (MCFT) yield surface
			//F_interlock = 4 * d * (std::abs(vciTR) / vcimax - k) + std::pow((sigTR / vcimax - h), 2);
			F_interlock = h * h / (k - chi_lok) * (std::abs(vciTR) / vcimax - k) + std::pow((sigTR / vcimax - h), 2);
		}
	}

	else
	{
		// Linear (constant crack width) yield surface
		F_interlock = std::abs(vciTR) + (sigTR - shift * vcimax);

		// ToDo: add compressive override here as well??
	}

	return F_interlock;
}


double MultiSurfCrack2D::getF_activation()
{
	double& sigTR = stress_current[0];
	double& vciTR = stress_current[1];

	// Hyperbolic yield surface for cohesion-frictional sliding reversal
	double F_activation = vciTR * vciTR - std::pow(k * chi_act * vcimax, 2) - std::pow(mu * sigTR, 2);

	return F_activation;
}


// stress-return algorithm
void MultiSurfCrack2D::plasticReturn()
{
	// if entering this method, FTR > 0, hence do-while loop

	// newton-raphson method with updating flow direction (general cutting plane algorithm)
	const double F_TOL = 1e-8;
	int itercnt = 0;
	const int CNT_MAX = 1000;

	do
	{
		// update flow direction at each iteration
		if (active_surface == "interlock")
			updateFlowDir();
		
		// estimate plastic multiplier to return stress state to yield surface
		double dlam = incrementPlasticMultiplier();

		// add the incremental flow to the accumulated flow during the current iterative plastic return step
		dLambda += dlam;

		// updated stresses (still trial)
		trialStressReturn();

		// check active yield function at new approximation
		checkActiveYieldSurface(active_surface);
		//checkActiveYieldSurface();

		// emergency exit (non-convergence of newton-raphson solver)
		itercnt += 1;
		if (itercnt > CNT_MAX)
		{
			// PRINT opserr WARNING
			break;
		}
	} while (std::abs(yield_criterion) > F_TOL);		// (yield_criterion > F_TOL);

	std::cout << "crack return-map iterations: " << itercnt << "; f: " << yield_criterion << " dir: " << dir << std::endl;

}


double MultiSurfCrack2D::incrementPlasticMultiplier()
{
	Matrix& De = tangent_elastic_damaged;
	Vector& dF = *dF_engaged;
	Vector& dG = *dG_engaged;

	// generalized cutting plane algorithm
	double dlam = yield_criterion / (dF ^ (De * dG));

	return dlam;
}


void MultiSurfCrack2D::trialStressReturn()
{
	// updated plastic strains
	Vector depsP = dLambda * (*dG_engaged);
	strain_plastic_current = strain_plastic_committed + depsP;

	// toDo: update crushing_strain here?? Or wait 'til after each commit?

	// some shortcut parameters
	double& s = strain_current[1];
	double& ds = strain_increment[1];
	double sgn_s = (s > 0) ? 1 : ((s < 0 ? -1 : 0));
	double& dw = strain_increment[0];
	double& dsP = depsP[1];
	double& dwP = depsP[0];
	Matrix& De = tangent_elastic_damaged;
	
	// updated stress
	double dsig = De(0, 0) * (dw - dwP) + De(0, 1) * (ds - dsP) * sgn_s;
	double dvci = De(1, 0) * (dw - dwP) * sgn_s + De(1, 1) * (ds - dsP);
	stress_current[0] = stress_committed[0] + dsig;
	stress_current[1] = stress_committed[1] + dvci;
}


// get flow vectors
void MultiSurfCrack2D::updateFlowDir()
{
	updateInterlockFlowDir();
	updateActivationFlowDir();
}


void MultiSurfCrack2D::updateInterlockFlowDir()
{
	double sigTR = stress_current[0];
	double vciTR = stress_current[1];
	
	// get sign of vci
	double sign_vci = (vciTR > 0) ? 1 : ((vciTR < 0) ? -1 : 1);

	if (!criticalPath)
	{
		// override for tensile region
		if (sigTR > 0)
		{
			dF_lok[0] = 2 * std::pow(k - chi_lok, 2) * (sigTR - l * h * k * chi_lok * vcimax / (k - chi_lok));
			dF_lok[1] = l * h * h * vciTR;
		}

		// override for compressive region
		else if (sigTR < h * vcimax)
		{
			dF_lok[0] = 2 * (sigTR - h * vcimax) / (a_c * a_c);
			dF_lok[1] = 2 * vciTR / (b_c * b_c);
		}

		else
		{
			// parabolic MCFT yield parabola
			dF_lok[0] = 2.0 / vcimax * (sigTR / vcimax - h);
			dF_lok[1] = h * h / vcimax / (k - chi_lok) * sign_vci;
		}
	}

	else
	{
		// linear "critical path" parabola
		dF_lok[0] = 1.0;
		dF_lok[1] = sign_vci;

	}

	// also, non-associated flow
	dG_lok[0] = rho_lok * dF_lok[0];
	dG_lok[1] =           dF_lok[1];
}


void MultiSurfCrack2D::updateActivationFlowDir()
{
	double sigTR = stress_current[0];
	double vciTR = stress_current[1];

	// use current stresses to calculate flow direction
	dF_act[0] = -2 * mu * mu * sigTR;
	dF_act[1] = 2 * vciTR;

	// also non-associated flow
	dG_act[0] = rho_act * dF_act[0];
	dG_act[1] =           dF_act[1];
}


void MultiSurfCrack2D::setvcimax()
{
	// update crushing damage to f'c
	double Ecp = get_postPeakCompressionSoftening();					// don't need to keep repeating this function call
	fcc = std::max(0.2*fc, fc + Ecp * crushing_strain);					// recall fc is positive values

	// get effective crack width, including sliding damage
	double w_tilde = get_w_tilde();

	// update vcimax (both crushing and sliding damage)
	vcimax = std::sqrt(fcc) / (0.31 + 24 * w_tilde / (ag + 16));

	// update compression ellipse
	a_c = -fcc - h * vcimax;
	b_c = k * vcimax;

	// update interlock cohesion
	c_lok = k * chi_lok * vcimax;

	// update activation cohesion
	c_act = k * chi_act * vcimax;

	// update traction ellipse
	ft = k * chi_t * vcimax;
	xc = k * chi_c * vcimax;
	a_t = ft - xc;
	b_t = a_t * std::sqrt(h * c_lok / (2 * xc * d));
}


double MultiSurfCrack2D::get_w_tilde()
{
	double w = strain_current[0] + w_res;													// total crack width
	double pi = 2 * std::acos(0.0);
	double w_tilde = std::max(0.0, w) + 8 / (3 * pi) * std::sin(theta_tilde) * sDAM;		// effective crack width, including sliding damage

	return w_tilde;
}


double MultiSurfCrack2D::get_postPeakCompressionSoftening()
{
	// Softening Modulus (MPa/mm) -- positive valued
	double Ecp = 100.0;				

	// ToDo: calculate Ecp with user-input G_f^c and characteristic length


	return Ecp;
}


void MultiSurfCrack2D::trialTangentStiffness()
{
	if (isElastic)
	{
		tangent_current = tangent_elastic_damaged;
	}
	else
	{
		Matrix& De = tangent_elastic_damaged;

		// stress vector contributions
		Vector dF = *dF_engaged;							// stress gradient for active yield surface
		Vector dG = *dG_engaged;							// flow vector
 
		// strain-hardening contributions
		double f_w = d_F_dw();								// yield surface sensitivity to crack width
		double f_nu = d_F_dnu();							// yield surface sensitivity to crushing -- contributes to elastoplastic tangent if...
		Vector w_eps = strain_gradient_w_tilde();			// strain gradient of effective width, omega_tilde
		Vector nu_eps = strain_gradient_nu();				// strain gradient of crushing strain, nu

		// shape change contributions
		double f_chi = d_F_dchi();							// variation of yield surface with shape parameter chi_lok
		double f_k = d_F_dk();								// variation of yield surface with shape parameter k
		Vector chi_eps = strain_gradient_chi();				// strain gradient of cohesion shape parameter chi_lok
		Vector k_eps = strain_gradient_k();					// strain gradient of size shape parameter k
 
		Vector Dg = De * dG;
		Vector fD = De ^ dF;
		Vector v_strain = w_eps * f_w + nu_eps * f_nu;
		Vector v_shape = f_chi * chi_eps + f_k * k_eps;
		tangent_current = De - (Dg % (fD + v_strain + v_shape)) / (dF ^ (De * dG));
	}
}


// partial derivative of yield surface w/rsp crack width
double MultiSurfCrack2D::d_F_dw()
{
	double f_w = d_F_dvcimax() * d_vcimax_d_w_tilde();
	return f_w;
}


// partial derivative of yield surface w/rsp crushing strain (nu in crack formulation paper)
double MultiSurfCrack2D::d_F_dnu()
{
	double Ecp = get_postPeakCompressionSoftening();

	double f_nu = d_F_dvcimax() * d_vcimax_dnu() + d_F_dfcc() * Ecp;

	return f_nu;
}


// partial derivative of yield surface w/rsp size parameter d
double MultiSurfCrack2D::d_F_dchi()
{
	double f_chi = 0.0;

	double& sigTR = stress_current[0];
	double& vciTR = stress_current[1];

	if (active_surface == "interlock")
	{
		// compressive surface
		if (sigTR < h * vcimax)
		{
			f_chi = 0.0;
		}

		// tensile surface
		else if (sigTR > 0)
		{
			//f_chi = 2 / h * chi_c / chi_lok * vciTR * vciTR;

			double hkv = h * k * vcimax;
			f_chi = -2 * (sigTR + l * hkv) * ((k - chi_lok) * sigTR - l * chi_lok * hkv) - 2 * chi_lok * std::pow((n + l) * hkv, 2);
		}

		// MCFT surface
		else
		{
			//f_chi = 4.0 * (std::abs(vciTR) / vcimax - k);
			f_chi = std::pow(h / (k - chi_lok), 2) * (std::abs(vciTR) / vcimax - k);
		}
	}

	else if (active_surface == "activation")
	{
		f_chi = 0.0;
	}

	return f_chi;
}


// partial derivative of yield surface w/rsp size parameter k
double MultiSurfCrack2D::d_F_dk()
{
	double f_k = 0.0;

	double& sigTR = stress_current[0];
	double& vciTR = stress_current[1];

	if (active_surface == "interlock")
	{
		// compressive surface
		if (sigTR < h * vcimax)
		{
			f_k = -2 / k * std::pow(vciTR / (k * vcimax), 2);
		}

		// tensile surface
		else if (sigTR > 0)
		{
			//f_k = -2 * chi_c * vcimax * (sigTR - k * chi_c * vcimax) - 2 * k * std::pow((chi_t - chi_c)*vcimax, 2);

			double hxv = h * chi_lok * vcimax;
			f_k = 2 * (sigTR - l * hxv) * ((k - chi_lok) * sigTR - l * k * hxv) - 2 * k * std::pow((n + l) * hxv, 2);
		}

		// MCFT surface
		else
		{
			//f_k = -4 * d;
			f_k = -h * h / (k - chi_lok);
		}
	}

	else if (active_surface == "activation")
	{
		f_k = -2 * k * std::pow(chi_act * vcimax, 2);
	}

	return f_k;
}


// derivative of vcimax w/rsp effective crack width
double MultiSurfCrack2D::d_vcimax_d_w_tilde()
{
	double w_tilde = get_w_tilde();

	double denom = 0.31 + 24 * w_tilde / (ag + 16.0);
	double deriv = -24.0 / (ag + 16.0) * std::sqrt(fcc) / (std::pow(denom, 2));

	return deriv;
}


// derivative of vcimax w/rsp crushing strain
double MultiSurfCrack2D::d_vcimax_dnu()
{
	double w_tilde = get_w_tilde();								// effective crack width

	double denom = 0.31 + 24 * w_tilde / (ag + 16.0);
	double Ecp = get_postPeakCompressionSoftening();
	double deriv = 1. / 2. / std::sqrt(fcc) * Ecp / denom;

	return deriv;
}


// derivative of yield surface w/rsp vcimax
double MultiSurfCrack2D::d_F_dvcimax()
{
	double deriv;
	double& sigTR = stress_current[0];
	double& vciTR = stress_current[1];

	if (active_surface == "interlock")
	{
		if (!criticalPath)
		{
			// override for compressive surface
			if (sigTR < h * vcimax)
			{
				//deriv = 2 * a_c * (sigTR + vcimax) * (a_c - (sigTR + vcimax)) / std::pow(a_c, 4) - 2 * vciTR * vciTR / std::pow(b_c, 3);
				deriv = 2 * h * (sigTR - h * vcimax) * ((sigTR - h * vcimax) - a_c) / std::pow(a_c, 3) - 2 / vcimax * std::pow(vciTR / b_c, 2);
			}

			// override for tensile surface
			else if (sigTR > 0)
			{
				//deriv = -2 * k * chi_c * (sigTR - xc) - 2 * vcimax * std::pow(k*(chi_t - chi_c), 2);

				double hkx = h * k * chi_lok;
				double k_minus_x = k - chi_lok;
				deriv = -2 * hkx * k_minus_x * (l*sigTR + (n*n + 2*n*l) * hkx/k_minus_x * vcimax);
			}

			else 
			{
				// MCFT surface
					//deriv = -2 / std::pow(vcimax, 2) * (2 * d * std::abs(vciTR) + sigTR * (sigTR / vcimax - h));
				deriv = -1 / std::pow(vcimax, 2) * (h * h / (k - chi_lok) * std::abs(vciTR) + 2 * sigTR * (sigTR / vcimax - h));
			}
		}

		else
		{
			deriv = -shift;
		}
	}
	else if (active_surface == "activation")
	{
		double radical = std::sqrt(1 + std::pow(mu * sigTR / (chi_act * vcimax), 2));
		//deriv = -radical + std::pow(mu * sigTR / (chi * vcimax), 2) / radical;
		//deriv *= chi;

		//deriv = std::abs(vciTR) / vcimax * std::pow(mu * sigTR / c_act, 2) / std::pow(radical, 3);
		//deriv -= chi_act;

		deriv = -2 * vcimax * std::pow(k * chi_act, 2);
	}

	return deriv;
}


// derivative of yield surface w/rsp reduced compression strength fcc
double MultiSurfCrack2D::d_F_dfcc()
{
	double deriv = 0.0;

	double& sigTR = stress_current[0];

	// this derivative only non-zero when compression surface is active	
	if (active_surface == "interlock" && sigTR < h * vcimax)
	{	
		deriv = 2 * std::pow(sigTR - h * vcimax, 2) / std::pow(a_c, 3);
	}

	return deriv;
}


// strain gradient of effective crack width omega_tilde
Vector MultiSurfCrack2D::strain_gradient_w_tilde()
{
	Vector w_tilde_eps(2);
	
	double w = strain_current[0] + w_res;					// total crack width
	w += 0.0000001;											// avoid div.by.zero

	double t_w = std::max(0., w) / w;						// toggle width part off if w < 0
	w_tilde_eps[0] = t_w;

	double& s = strain_current[1];
	double t_s;
	t_s = std::floor(std::abs(s) / (sDAM + 0.00001));		// toggle sliding part off if |s| < sDAM    (avoid div.by.zero)
	t_s = std::min(1.0, std::max(0.0, t_s));
	
	double pi = 2 * std::acos(0.0);
	w_tilde_eps[1] = 8 / (3 * pi) * std::sin(theta_tilde) * t_s;

	return w_tilde_eps;
}


// strain gradient of crushing strain nu
Vector MultiSurfCrack2D::strain_gradient_nu()
{
	Vector nu_eps(2);

	// toggle function (1 if compression surface engaged, else 0)
	double t_c;

	Vector dF = *dF_engaged;								// stress gradient for active yield surface

	dF[0] += 0.00000001;									// avoid div.by.zero
	t_c = std::min(0., dF[0]) / dF[0];						// ...outward normal points in negative crack direction (t_c = 1)
	t_c = (active_surface == "interlock") ? t_c : 0.0;		// ...and interlock surface is active

	nu_eps[0] = t_c;

	return nu_eps;
}


// strain gradient of shape parameter d
Vector MultiSurfCrack2D::strain_gradient_chi()
{
	Vector chi_eps(2);

	double d_chi_r = chi_lok_r - chi_lok_s;
	double d_chi_s = d_chi_r * d_r_ds();

	chi_eps[1] = d_chi_s;

	return chi_eps;
}


// strain gradient of shape parameter k
Vector MultiSurfCrack2D::strain_gradient_k()
{
	Vector k_eps(2);

	double k_r = k_rough - k_smooth;
	double k_sl = k_r * d_r_ds();

	k_eps[1] = k_sl;

	return k_eps;
}


int MultiSurfCrack2D::commitState()
{
	double& sP = strain_plastic_current[1];
	double& sCurr = strain_current[1];
	double wCurr = strain_current[0] + w_res;

	// update pinching break points (aggregate re-engagment)
	sMAX += (sP - sMAX) * (sP > sMAX);
	sMIN += (sP - sMIN) * (sP < sMIN);
	//sMAX += (sCurr - sMAX) * (sCurr > sMAX);
	//sMIN += (sCurr - sMIN) * (sCurr < sMIN);


	// update max crack width (for ruggiero stiffness update)
	wMAX += (wCurr - wMAX) * (wCurr > wMAX);

	// update max crack slip (for sliding damage calcs)
	sDAM += (std::abs(sCurr) - sDAM) * (std::abs(sCurr) > sDAM);

	// update crushing strain accumulation
	if ((active_surface == "interlock") && (stress_current[0] < h * vcimax))
	{
		// use total strains so post-peak response isn't affected by user's choice of rho^lok
		double& wPrev = strain_committed[0];
		double crush_increment = wCurr - wPrev;
		if (crush_increment < 0)
		{
			crushing_strain += crush_increment;
		}
	}									

	// update committed state
	strain_committed = strain_current;
	strain_plastic_committed = strain_plastic_current;
	stress_committed = stress_current;
	tangent_committed = tangent_current;
	dLambda = 0.0;

	// update shape of interlocking yield surface (only after committing state)
	if (sP && (!criticalPath))
	{
		updateInterlockShape();
	}
	
	//updateInterlockShape();					// updated even during elastic steps

	// update yield surface (with updated sliding damage, updated crushing damage and updated interlock shape)
	setvcimax();

	// store last shear point converged to interlock surface
	if (active_surface == "interlock")
	{
		last_vci_lok = stress_committed[1];
	}


	return 0;
}


int MultiSurfCrack2D::revertToLastCommit()
{
	strain_current = strain_committed;
	stress_current = stress_committed;
	strain_plastic_current = strain_plastic_committed;
	tangent_current = tangent_committed;
	dLambda = 0.;
	setvcimax();

	return 0;
}


void MultiSurfCrack2D::updateInterlockShape()
{
	double& vci = stress_committed[1];
	double& s = strain_committed[1];
	double& sP = strain_plastic_committed[1];
	
	double s_rough;

	if (vci > 0)
	//if (s > 0)
	{
		s_rough = 1.0 * sMAX;
	}
	else
	{
		s_rough = 1.0 * sMIN;
	}

	// intermediate slip variable
	double x = std::abs(s / s_rough);
	x = std::min(1.0, x);
	x = std::max(0.0, x);

	// update Roughness parameter based on current slip
	double r = (x - kappa * x) / (kappa - 2 * kappa * x + 1);

	// update shape of interlocking parabola, F^lok
	k = k_rough * r + k_smooth * (1 - r);
	chi_lok = chi_lok_r * r + chi_lok_s * (1 - r);
	h = -1;

	// update dependent shape parameters (d, chi_t, chi_c)
	updateDerivedShapeParameters();

	// don't let chi_lok underpass chi_act (avoid surface overlap)  --> set chi_lok_smooth = chi_act
}



void MultiSurfCrack2D::updateDerivedShapeParameters()
{
	// calculate derived shape parameters from updated k and chi_lok
	d = h * h / 4 / (k - chi_lok);
	chi_t = -m * h / 4 * chi_lok / (k - chi_lok);
	chi_c = m * m * h / (8 * (1 - m)) * chi_lok / (k - chi_lok);
}



double MultiSurfCrack2D::d_r_ds()
{
	double& vci = stress_committed[1];
	double& s = strain_committed[1];
	double& sP = strain_plastic_committed[1];
	
	double s_rough;
	if (vci >= 0)
	{
		s_rough = 1.0 * sMAX;
	}
	else
	{
		s_rough = 1.0 * sMIN;
	}

	// intermediate slip variable
	double x = std::abs(s / s_rough);
		//double x = s / s_rough;
	x = std::min(1.0, x);
	x = std::max(0.0, x);

	// toggle off if s = s_rough
	double t_d = 1 - std::floor(x);

	// avoid div.by.zero if kappa = -1 or kappa = 1
	x = std::min(0.999, x);
	x = std::max(0.001, x);

	double denom = kappa - 2 * kappa * x + 1;
	double r_s = 1 / s_rough * (1 - std::pow(kappa, 2)) / std::pow(denom,2);

	return r_s * t_d;
}


void MultiSurfCrack2D::PrintState()
{
	double& w = strain_current[0];
	double& s = strain_current[1];
	double& sig = stress_current[0];
	double& vci = stress_current[1];
	
	std::cout << "(w, s): (" << w << ", " << s << ")" << std::endl;
	std::cout << "\tdir: " << dir << std::endl;
	std::cout << "\tPHASE: " << active_surface << std::endl;
	std::cout << "\tsMAX_plus: " << sMAX << std::endl;
	std::cout << "\tsMAX_minus: " << sMIN << std::endl;
	std::cout << "\t(sig, vci): " << sig << ", " << vci << ")" << std::endl;
}


const Vector& MultiSurfCrack2D::getStrain()
{
	return strain_current;
}


const Vector& MultiSurfCrack2D::getStress()
{
	if (verbose)
	{
		opserr << "MultiSurfCrack2D::getStress()" << endln;
		opserr << "\tcrack forces: " << stress_current[0] * Acr << ", " << stress_current[1] * Acr << endln;
	}
	
	// send force back to zLND element
	force_current = stress_current * Acr;

	return force_current;
}


const Matrix& MultiSurfCrack2D::getTangent()
{
	if (verbose)
		opserr << "MultiSurfCrack2D::getTangent()" << endln;

	// want units of Force/Disp for zLND element tangent, so multiply by this crack element's area
	tangentF = tangent_current * Acr;

	//opserr << " matTangent: [[" << tangentF(0, 0) << ", " << tangentF(0, 1) << "],[" << tangentF(1, 0) << ", " << tangentF(1, 1) << "]]" << endln;
	
	return tangentF;
}


const Matrix& MultiSurfCrack2D::getInitialTangent()
{
	// multiply by this crack element's area to get units of Force/Disp
	tangentF = tangent_elastic_loading_initial * Acr;
	
	return tangentF;
}


int MultiSurfCrack2D::revertToStart()
{
	// ToDo: update to handle built-in stresses
	stress_current.Zero();
	stress_committed.Zero();
	
	strain_current.Zero();
	strain_committed.Zero();
	strain_plastic_current.Zero();
	strain_plastic_committed.Zero();
	setvcimax();
	sMAX = 0.0;
	sMIN = 0.0;
	wMAX = std::min(0.1, w_res);
	sDAM = 0.0;

	tangent_elastic_damaged = tangent_elastic_loading_initial;
	tangent_current = tangent_elastic_loading_initial;
	tangent_committed = tangent_elastic_loading_initial;

	return 0;
}


NDMaterial* MultiSurfCrack2D::getCopy() {

	MultiSurfCrack2D* theCrackMaterial = new MultiSurfCrack2D(this->getTag(),
										tangent_elastic_loading_initial,
										tangent_elastic_unloading_initial,
										fc,ag,fcl,Acr,rho_lok,chi_lok,
										rho_act,mu,chi_act,zeta,kappa,theta_tilde,
										strain_committed[0]+w_res,criticalPath);
	return theCrackMaterial;
};


NDMaterial* MultiSurfCrack2D::getCopy(const char* type) {

	// This crack model is limited to 2D and meant for use with zLND element -- warn users if type != "PlaneStrain2D", which comes from zeroLengthND element parsing

	if (strcmp(type, "PlaneStrain2D") != 0)
	{
		opserr << "WARNING: element must be of type PlaneStrain2D to handle crack material " << endln;
	}

	MultiSurfCrack2D* theCrackMaterial = new MultiSurfCrack2D(this->getTag(),
										tangent_elastic_loading_initial,
										tangent_elastic_unloading_initial,
										fc,ag,fcl,Acr,rho_lok,chi_lok,
										rho_act,mu,chi_act,zeta,kappa,theta_tilde,
										strain_committed[0]+w_res,criticalPath);
	return theCrackMaterial;
};


int MultiSurfCrack2D::sendSelf(int commitTag, Channel& theChannel)
{
	int res = 0;

	static Vector data(16+2+2+2+4+1);

	data(0) = this->getTag();

	// package up elastic stiffness components in both loading and unloading
	//data(1) = tangent_elastic_loading_initial(0, 0);
	//data(2) = tangent_elastic_loading_initial(0, 1);
	//data(3) = tangent_elastic_loading_initial(1, 0);
	//data(4) = tangent_elastic_loading_initial(1, 1);
	//data(5) = tangent_elastic_unloading_initial(0, 0);
	//data(6) = tangent_elastic_unloading_initial(0, 1);
	//data(7) = tangent_elastic_unloading_initial(1, 0);
	//data(8) = tangent_elastic_unloading_initial(1, 1);
	// 
	// package up concrete properties to be able to define size of yield surface
	data(1) = fc;
	data(2) = ag;
	data(3) = fcl;
	// package up crack area to convert stress --> force for use with zLND ele
	data(4) = Acr;
	// package up model parameters
	data(5) = rho_lok;
	data(6) = chi_lok;
	data(7) = rho_act;
	data(8) = mu;
	data(9) = chi_act;
	data(10) = zeta;
	data(11) = kappa;
	data(12) = theta_tilde;
	// package up slip history variables
	data(13) = sMAX;
	data(14) = sMIN;
	data(15) = sDAM;

	// need to also recreate current stress and strain state
	for (int i = 0; i < 2; i++) {
		data(16 +     i) = strain_committed[i];
		data(16 + 2 + i) = strain_plastic_committed[i];
		data(16 + 4 + i) = stress_committed[i];
	}
	
	// package up current tangent
	data(16+6) = tangent_committed(0, 0);
	data(16+6+1) = tangent_committed(0, 1);
	data(16+6+2) = tangent_committed(1, 0);
	data(16+6+3) = tangent_committed(1, 1);

	data(16+6+4) = criticalPath;

	res = theChannel.sendVector(this->getDbTag(), commitTag, data);
	if (res < 0)
		opserr << "Crack::sendSelf() - failed to send data" << endln;

	return res;
}


int MultiSurfCrack2D::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
	static Vector data(16+2+2+2+4+1);
	if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
		opserr << "Crack::recvSelf - failed to recv vector from channel\n";
		return -1;
	}

	this->setTag((int)data(0));
	 
	//tangent_elastic_loading_initial(0, 0) = data(1);
	//tangent_elastic_loading_initial(0, 1) = data(2);
	//tangent_elastic_loading_initial(1, 0) = data(3);
	//tangent_elastic_loading_initial(1, 1) = data(4);
	//tangent_elastic_unloading_initial(0, 0) = data(5);
	//tangent_elastic_unloading_initial(0, 1) = data(6);
	//tangent_elastic_unloading_initial(1, 0) = data(7);
	//tangent_elastic_unloading_initial(1, 1) = data(8);

	fc = data(1);
	ag = data(2);
	fcl = data(3);
	Acr = data(4);
	rho_lok = data(5);
	chi_lok = data(6);
	rho_act = data(7);
	mu = data(8);
	chi_act = data(9);
	zeta = data(10);
	kappa = data(11);
	theta_tilde = data(12);
	sMAX = data(13);
	sMIN = data(14);
	sDAM = data(15);

	for (int i = 0; i < 2; i++) {
		strain_committed[i] = data(16 + i);
		strain_plastic_committed[i] = data(16 + 2 + i);
		stress_committed[i] = data(16 + 4 + i);
	}

	tangent_committed(0, 0) = data(16 + 6);
	tangent_committed(0, 1) = data(16 + 6 + 1);
	tangent_committed(1, 0) = data(16 + 6 + 2);
	tangent_committed(1, 1) = data(16 + 6 + 3);

	criticalPath = (int)data(16+6+4);

	// set trial values to committed values
	strain_current = strain_committed;
	strain_plastic_current = strain_plastic_committed;
	stress_current = stress_committed;
	tangent_current = tangent_committed;

	// also set current value of cohesion
	this->setvcimax();

	// re-initialize flow direction
	this->updateFlowDir();

	return 0;
}


void MultiSurfCrack2D::Print(OPS_Stream& s, int flag)
{
	s << "MultiSurfCrack2D Material tag: " << this->getTag() << endln;
	s << " ...printout not yet implemented ";
}