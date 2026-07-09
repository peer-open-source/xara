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
Source: /DEVELOPER/element/cpp/Macroelement3d/Macroelement3d/Macroelement3d.cpp
Written by Francesco Vanin (francesco.vanin@epfl.ch)
Ecole Polytechnique Federale de Lausanne, Switzerland,
Earthquake Engineering and Structural Dynamics laboratory, 2019

Reference: Vanin F., Penna A., Beyer K.;"A three dimensional macro-element
for modelling of the in-plane and out-of-plane response of masonry walls",
submitted to Earthquake Engineering and Structural Dynamics (2019)

Last edit: 27 Feb 2019
*/

/* -------------------------------------------------------------------------- */
#include "Macroelement3d.h"
#include "DamageShearInterface.h"
#include "GambarottaLagomarsinoModel.h"
#include "NoTensionSection3d.h"
#include "GenericDamagePlasticityShear.h"
#include "WrappedMaterial.h"
/* -------------------------------------------------------------------------- */
#include <Channel.h>
#include <CompositeResponse.h>
#include <Domain.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <FEM_ObjectBroker.h>
#include <ID.h>
#include <Information.h>
#include <Matrix.h>
#include <Node.h>
#include <Parameter.h>
#include <Renderer.h>
#include <SectionForceDeformation.h>
#include <TransientIntegrator.h>
#include <Vector.h>
#include <elementAPI.h>
/* -------------------------------------------------------------------------- */
#include <algorithm>
#include <limits>
#include <string>

#include <cmath>
/* -------------------------------------------------------------------------- */
#ifndef DBL_EPSILON
#define DBL_EPSILON (std::numeric_limits<double>::epsilon())
#endif
#include <TransientIntegrator.h>
#include <CompositeResponse.h>
#include <typeinfo>
#include <string>



Matrix Macroelement3d::K(18,18);
Vector Macroelement3d::P(18);
double Macroelement3d::workArea[200];

#ifdef _USRDLL
#include <windows.h>
#define OPS_Export extern "C" _declspec(dllexport)
#elif _MACOSX
#define OPS_Export extern "C" __attribute__((visibility("default")))
#else
#define OPS_Export extern "C"
#endif

static int numMyMacroelement = 0;

OPS_Export 	void * 
OPS_Macroelement3d()
{
	Vector intLength(3);
	intLength(0) = 1.0/6;
	intLength(1) = 2.0/3;
	intLength(2) = 1.0/6;

	Vector intLengthMasses(3);
	/*
	intLengthMasses(0) = 0.1666;
	intLengthMasses(1) = 0.6666;
	intLengthMasses(2) = 0.1666;
	*/

	intLengthMasses(0) = 0.25;
	intLengthMasses(1) = 0.50;
	intLengthMasses(2) = 0.25;

	Vector driftModelF(0);
	Vector driftModelF_ALR(0);
	Vector driftModelS(0);
	Vector driftModelS_ALR(0);

	double Ltfc = 0.;
	double alphaAC_HC = -4.0/3.0;
	double failureFactorF = 0.001;
	double failureFactorS = 0.001;
	double betaShearSpanF = 0.0;
	double betaShearSpanS = 0.0;
	
	
    if (numMyMacroelement  == 0) {
        opserr << "Macroelement3d loaded from external library - Written by Francesco Vanin, EPFL, 2019" << endln;
        numMyMacroelement++;
    }

    int remaining=OPS_GetNumRemainingInputArgs();

    if (remaining < 6) {
		opserr << "WARNING: insufficient arguments for element geometry." << endln;
		opserr << "Required: eleTag, iNode, jNode, eNode, axisX, axisY, axisZ, oopX, oopY, oopZ <-flags>" << endln;
	return 0;
    }

    // inputs: 
    int iData[4];
    int numData = 4;
    if(OPS_GetIntInput(&numData,&iData[0]) < 0) {
	opserr<<"WARNING: invalid integer inputs, element geometry definition." << endln;
	return 0;
    }

    // inputs: 
    double dData[6];
    numData = 6;
    if(OPS_GetDoubleInput(&numData,&dData[0]) < 0) {
	opserr<<"WARNING: invalid double inputs, element geometry definition (tag: "<< iData[0] << ")" << endln;
	return 0;
    }

    Vector axis(3); 
    axis(0) = dData[0];
    axis(1) = dData[1];
    axis(2) = dData[2];

	Vector oop(3);
    oop(0) = dData[3]; 
    oop(1) = dData[4]; 
    oop(2) = dData[5]; 

	SectionForceDeformation* theSectionI{nullptr};
    SectionForceDeformation* theSectionE{nullptr};
	SectionForceDeformation* theSectionJ{nullptr};


    NDMaterial* theShearModel{nullptr};
	NDMaterial* theShearModelOOP{nullptr};
    double E_, h; 
	
	bool gable = false;
	double b = 0.;
	double t = 0.;

    // options for input
	const char* inputStructure = OPS_GetString();
	if (strcmp(inputStructure,"-tremuri") == 0)  {
		// standard input, model equivalent to Tremuri (in-plane)
		// arguments: h,b,t,E,G,fc,mu,tau0,Gc,beta
		double dData2[10];
	    numData = 10;
		if (OPS_GetDoubleInput(&numData,&dData2[0]) < 0) {
			opserr << "WARNING Macroelement3d (tag: " << iData[0] << "):'-tremuri' input structure incorrect, invalid double input(s)." << endln;
			opserr << "Required structure : h, b, t, E, G, fc, mu, tau0, Gc, beta <-flags>" << endln;
			return 0;
		}

		intLength(0) = 0.495;
		intLength(1) = 0.010;
		intLength(2) = 0.495;

		h = dData2[0];
		b = dData2[1];
		t = dData2[2];
		E_ = dData2[3];
		double G = dData2[4];
		double fc = dData2[5];
		double  mu = dData2[6];
		double  tau0 = dData2[7];
		double Gc = dData2[8];
		double beta = dData2[9];

		Ltfc = std::abs(fc*b*t);
		
		// No tension 3d section model
		// true for stronger, true for elastic, true for crushing
		theSectionI = new NoTensionSection3d(0, E_, G, t, b, -1.0, dData2[5], 5, false, false,  true);   	
		theSectionE = new NoTensionSection3d(0, E_, G, t, b, -1.0, dData2[5], 5, true,  false,  true);
		theSectionJ = new NoTensionSection3d(0, E_, G, t, b, -1.0, dData2[5], 5, false, false,  true);

		// Gambarotta Lagomarsino model for shear
		// true-false for elastic solution or not. The 5/6 factor is dropped to make it equivalent to Tremuri
		theShearModelOOP = new GambarottaLagomarsinoModel(0, E_, G, tau0, mu, Gc/G, beta, b, t, h, true);  
		theShearModel    = new GambarottaLagomarsinoModel(0, E_, G, tau0, mu, Gc/G, beta, b, t, h, false);
			
     } 




		else if (strcmp(inputStructure,"-tremuriIP") == 0) {
		// standard input, model equivalent to Tremuri (thinner element with only in-plane response)
		// arguments: h,b,t, E,G, fc, mu,tau0, Gc,beta
		double dData2[10];
	    numData = 10;
		if (OPS_GetDoubleInput(&numData,&dData2[0]) < 0) {
			opserr << "WARNING Macroelement3d (tag: " << iData[0] << "):'-tremuriIP' input structure incorrect, invalid double input(s)." << endln;
			opserr << "Required structure : h, b, t, E, G, fc, mu, tau0, Gc, beta <-flags>" << endln;
			return 0;
		}

		intLength(0) = 0.495;
		intLength(1) = 0.010;
		intLength(2) = 0.495;

		intLengthMasses = intLength;

		h = dData2[0];
		b = dData2[1];
		t = dData2[2];
		E_ = dData2[3];
		double G = dData2[4];
		double fc = dData2[5];
		double  mu = dData2[6];
		double  tau0 = dData2[7];
		double Gc = dData2[8];
		double beta = dData2[9];

		Ltfc = std::abs(fc*b*t);

		// true for stronger, true for elastic, true for crushing
		double thinner = 100.;

		theSectionI = new NoTensionSection3d(0, thinner*E_, thinner*G, t / thinner, b, -1.0, thinner*dData2[5], 5, false, false, true);
		theSectionE = new NoTensionSection3d(0, thinner*E_, thinner*G, t / thinner, b, -1.0, thinner*dData2[5], 5, false, true, true);
		theSectionJ = new NoTensionSection3d(0, thinner*E_, thinner*G, t / thinner, b, -1.0, thinner*dData2[5], 5, false, false, true);

		// Gambarotta Lagomarsino model for shear
		// true-false for elastic solution or not. The 5/6 factor is dropped to make it equivalent to Tremuri
		theShearModelOOP = new GambarottaLagomarsinoModel(0, E_, G, tau0, mu, Gc / G, beta, b, t, h, true);
		theShearModel    = new GambarottaLagomarsinoModel(0, E_, G, tau0, mu, Gc / G, beta, b, t, h, false);

     } 




		else if (strcmp(inputStructure,"-fiberSection") == 0) {
		// input for the use of one or more fiber section models for the rocking behaviour, coupled to a standard shear model
	    // two-line input for standard rectangular prismatic elements-> line 1: section definition, line 2: macroelement definition 
	    // arguments: sectionI, sectionJ, sectionK, h,b,t, E,G, fc, mu,tau0, Gc,beta

		int iData2[3];
		int numData = 3;
		if (OPS_GetIntInput(&numData,&iData2[0]) < 0) {
			opserr << "WARNING Macroelement3d (tag: " << iData[0] << "): '-fiberSection' input structure incorrect, invalid integer input(s)." << endln;
			opserr << "Required structure : sectionI, sectionE, sectionJ, h, b, t, E, G, fc, mu, tau0, Gc, beta <-flags>" << endln;
			return 0;
		}

		
		// create copies of the sectional models
		int secTag = iData2[0];
	    theSectionI = OPS_GetSectionForceDeformation(secTag);
		if (theSectionI == 0) {
			opserr<<"SectionI: section model tagged " << secTag << " not found." << endln;
			return 0;
		}

		secTag = iData2[1];
		theSectionE = OPS_GetSectionForceDeformation(secTag);
		if (theSectionE == 0) {
			opserr<<"SectionE: section model tagged " << secTag << " not found." << endln;
			return 0;
		}

		secTag = iData2[2];
	    theSectionJ = OPS_GetSectionForceDeformation(secTag);
		if (theSectionJ == 0) {
			opserr<<"SectionJ: section model tagged " << secTag << " not found." << endln;
			return 0;
		}


		double dData2[10];
	    numData = 10;
		if (OPS_GetDoubleInput(&numData,&dData2[0]) < 0) {
			opserr << "WARNING Macroelement3d (tag: " << iData[0] << "): '-fiberSection' input structure incorrect, invalid integer input(s)." << endln;
			opserr << "Required structure : sectionI, sectionE, sectionJ, h, b, t, E, G, fc, mu, tau0, Gc, beta <-flags>" << endln;
			return 0;
		}


		h = dData2[0];
		b = dData2[1];
		t = dData2[2];
		E_ = dData2[3];
		double G = dData2[4];
		double fc = dData2[5];
		double  mu = dData2[6];
		double  tau0 = dData2[7];
		double Gc = dData2[8];
		double beta = dData2[9];

		Ltfc = std::abs(fc*b*t);

		// Gambarotta Lagomarsino model for shear
		// true-false for elastic solution or not. The 5/6 factor is dropped to make it equivalent to Tremuri
		theShearModelOOP = new GambarottaLagomarsinoModel(0, E_, G, tau0, mu, Gc/G, beta, b, t, h, true);  		
		theShearModel    = new GambarottaLagomarsinoModel(0, E_, G, tau0, mu, Gc/G, beta, b, t, h, false);
     } 



	else if (strcmp(inputStructure,"-fiberSectionShearModel1d") == 0) {
		// input for the use of one or more fiber section models for the rocking behaviour, plus a uniaxial shear model, uncoupled from the axial load
	    // arguments: sectionI, sectionJ, sectionK, ShearModel1dIP, ShearModel1dOOP

		int iData2[5];
		int numData = 5;
		if (OPS_GetIntInput(&numData,&iData2[0]) < 0) {
			opserr << "WARNING Macroelement3d (tag: " << iData[0] << "): '-fiberSectionShearModel1d' input structure incorrect, invalid integer input(s)." << endln;
			opserr << "Required structure : sectionI, sectionE, sectionJ, shearModel1dIP, shearModel1dOOP, alpha, h <-flags>" << endln;
			return 0;
		}

		// create copies of the sectional models
		int secTag = iData2[0];
	    theSectionI = OPS_GetSectionForceDeformation(secTag);
		if (theSectionI == 0) {
			opserr<<"SectionI: section model tagged " << secTag << " not found." << endln;
			return 0;
		}

		secTag = iData2[1];
		theSectionE = OPS_GetSectionForceDeformation(secTag);
		if (theSectionE == 0) {
			opserr<<"SectionE: section model tagged " << secTag << " not found." << endln;
			return 0;
		}

		secTag = iData2[2];
	    theSectionJ = OPS_GetSectionForceDeformation(secTag);
		if (theSectionJ == 0) {
			opserr<<"SectionJ: section model tagged " << secTag << " not found." << endln;
			return 0;
		}

		E_ = (theSectionE->getInitialTangent())(0,0);
		

		double dData2[1];
	    numData = 1;
		if (OPS_GetDoubleInput(&numData,&dData2[0]) < 0) {
			opserr << "WARNING Macroelement3d (tag: " << iData[0] << "):'-fiberSectionShearModel1d' input structure incorrect, invalid double input(s)." << endln;
			opserr << "Required structure : sectionI, sectionE, sectionJ, shearModel1dIP, shearModel1dOOP, h <-flags>" << endln;
			return 0;
		}

		intLength(0) = 0.495;
		intLength(1) = 0.010;
		intLength(2) = 0.495;

		t = 0.;
		h = dData2[0];
			
		// create pointers to uniaxial shear model
		UniaxialMaterial* pointerShear1dIP;
		UniaxialMaterial* pointerShear1dOOP;

		int matTag = iData2[3];
	    pointerShear1dIP = OPS_GetUniaxialMaterial(matTag);
		if (pointerShear1dIP== 0) {
			opserr<<"Shear1dIP: material model tagged " << secTag << " not found." << endln;
			return 0;
		}

		matTag = iData2[4];
	    pointerShear1dOOP = OPS_GetUniaxialMaterial(matTag);
		if (pointerShear1dOOP== 0) {
			opserr<<"Shear1dOOP: material model tagged " << secTag << " not found." << endln;
			return 0;
		}

		theShearModel     = new WrappedMaterial(0, E_, pointerShear1dIP);
		theShearModelOOP  = new WrappedMaterial(0, E_, pointerShear1dOOP);
     } 



	
	else if (strcmp(inputStructure,"-pier") == 0)  {
		// standard input for piers, with my shear model and 3 inelastic sections
		// arguments: h, b, t, E, G, fc, mu, c, Gc, dropDrift, muR 

		double dData2[11];
	    numData = 11;
		if (OPS_GetDoubleInput(&numData,&dData2[0]) < 0) {
			opserr << "WARNING Macroelement3d (tag: " << iData[0] << "):'-pier' input structure incorrect, invalid double input(s)." << endln;
			opserr << "Required structure : h, b, t, E, G, fc, mu, c, Gc, dropDrift, muR <-flags>" << endln;
			return 0;
		}

		intLength(0) = 1./6.;
		intLength(1) = 2./3.;
		intLength(2) = 1./6.;

		h = dData2[0];	

		b = dData2[1];
		t = dData2[2];
		E_ = dData2[3];
		double G = dData2[4];
		double fc = dData2[5];
		double  mu = dData2[6];
		double  tau0 = dData2[7];
		double Gc = 1.0 + dData2[8];
		double dropDrift = dData2[9];
		double muR = dData2[10];

		Ltfc = std::abs(fc*b*t);

		// true for stronger, true for elastic, true for crushing
		theSectionI = new NoTensionSection3d(0, E_, G, t, b, -1.0, fc, 11, false, false, true);   	
		theSectionE = new NoTensionSection3d(0, E_, G, t, b, -1.0, fc, 11, false, false, true);
		theSectionJ = new NoTensionSection3d(0, E_, G, t, b, -1.0, fc, 11, false, false, true);

		theShearModelOOP = new DamageShearInterface(0, E_*b*t, G/(h)*5./6.*b*t, tau0*b*t, mu, muR, Gc, dropDrift*h, true);
		theShearModel    = new DamageShearInterface(0, E_*b*t, G/(h)*5./6.*b*t, tau0*b*t, mu, muR, Gc, dropDrift*h, false);										
     }
	


	else if (strcmp(inputStructure, "-pierGeneric") == 0) {
		// input for piers, with generic shear model to include in the end. Needs to substitute in the long term the -pier definition
		// left for back compatibility
		// arguments: h, b, t, E, G, fc "shearModel" <variable shear model parameters>

		double dData2[6];
		numData = 6;
		if (OPS_GetDoubleInput(&numData, &dData2[0]) < 0) {
			opserr << "WARNING Macroelement3d (tag: " << iData[0] << "):'-pierGeneric' input structure incorrect, invalid double input(s)." << endln;
			opserr << "Required structure : h, b, t, E, G, fc shearModel <variable shear model parameters> <-flags>" << endln;
			return 0;
		}

		intLength(0) = 1. / 6.;
		intLength(1) = 2. / 3.;
		intLength(2) = 1. / 6.;

		h = dData2[0];

		b = dData2[1];
		t = dData2[2];
		E_ = dData2[3];
		double G = dData2[4];
		double fc = dData2[5];

		Ltfc = abs(fc*b*t);

		// true for stronger, true for elastic, true for crushing
		theSectionI = new NoTensionSection3d(0, E_, G, t, b, -1.0, fc, 11, false, false, true);
		theSectionE = new NoTensionSection3d(0, E_, G, t, b, -1.0, fc, 11, false, false, true);
		theSectionJ = new NoTensionSection3d(0, E_, G, t, b, -1.0, fc, 11, false, false, true);

		const char* type = NULL;
        type = OPS_GetString();

		if (strcmp(type, "-MohrCoulomb") == 0) {
			double dataShear[5];
			numData = 5;

			if (OPS_GetDoubleInput(&numData, &dataShear[0]) < 0) {
				opserr << "WARNING Macroelement3d (tag: " << iData[0] << "):unable to read shear model input(s)." << endln;
				opserr << "Required structure: -MohrCoulomb mu c Gc dropDrift alpha <-flags>" << endln;
				return 0;
			}

			Vector props(4);
			props(0) = b;
			props(1) = t;
			props(2) = dataShear[1];
			props(3) = dataShear[0];

			theShearModel = new GenericDamagePlasticityShear(0, G/(h)*5./6.*b*t, 1, props, 1.0 + dataShear[2], dataShear[3] *h, dataShear[4], false, 0.0);
			theShearModelOOP = new GenericDamagePlasticityShear(0, G / (h)*5. / 6.*b*t, 1, props, 1.0 + dataShear[2], dataShear[3] * h, dataShear[4], true, 0.0);
		}
		else if (strcmp(type, "-TurnsekCacovic")==0) {
			double dataShear[5];
			numData = 5;

			if (OPS_GetDoubleInput(&numData, &dataShear[0]) < 0) {
				opserr << "WARNING Macroelement3d (tag: " << iData[0] << "):unable to read shear model input(s)." << endln;
				opserr << "Required structure: -MohrCoulomb ft b Gc dropDrift alpha <-flags>" << endln;
				return 0;
			}

			Vector props(4);
			props(0) = b;
			props(1) = t;
			props(2) = dataShear[1];
			props(3) = dataShear[0];

			theShearModel = new GenericDamagePlasticityShear(0, G / (h)*5. / 6.*b*t, 2, props, 1.0+dataShear[2], dataShear[3] * h, dataShear[4], false, 0.0);
			theShearModelOOP = new GenericDamagePlasticityShear(0, G / (h)*5. / 6.*b*t, 2, props, 1.0+dataShear[2], dataShear[3] * h, dataShear[4], true, 0.0);

		}
		else if ((strcmp(type, "-MohrCoulombCompressed") == 0) || (strcmp(type, "-Eurocode") == 0)) {
			// [L, t, c, mu, 0.065fm]
			double dataShear[6];
			numData = 6;

			if (OPS_GetDoubleInput(&numData, &dataShear[0]) < 0) {
				opserr << "WARNING Macroelement3d (tag: " << iData[0] << "):unable to read shear model input(s)." << endln;
				opserr << "Required structure: -MohrCoulomb mu c 0.065fm Gc dropDrift alpha <-flags>" << endln;
				return 0;
			}

			Vector props(5);
			props(0) = b;
			props(1) = t;
			props(2) = dataShear[1];
			props(3) = dataShear[0];
			props(4) = dataShear[2];

			theShearModel = new GenericDamagePlasticityShear(0, G / (h)*5. / 6.*b*t, 3, props, 1.0 + dataShear[3], dataShear[4] * h, dataShear[5], false, 0.0);
			theShearModelOOP = new GenericDamagePlasticityShear(0, G / (h)*5. / 6.*b*t, 3, props, 1.0 + dataShear[3], dataShear[4] * h, dataShear[5], true, 0.0);

		}

	}



	
	else if (strcmp(inputStructure,"-spandrel") == 0) {
        // standard input for spandrels, with my shear model and the middle section linear elastic (end sections slightly stronger)
		// arguments: h, b, t, E, G, fc, mu, c, Gc, dropDrift, muR  

		double dData2[11];
	    numData = 11;
		if (OPS_GetDoubleInput(&numData,&dData2[0]) < 0) {
			opserr << "WARNING Macroelement3d (tag: " << iData[0] << "):'-spandrel' input structure incorrect, invalid double input(s)." << endln;
			opserr << "Required structure : h, b, t, E, G, fc, mu, c, Gc, dropDrift, muR <-flags>" << endln;
			return 0;
		}

		intLength(0) = 1./6.;
		intLength(1) = 2./3.;
		intLength(2) = 1./6.;

		h = dData2[0];

		b = dData2[1];
		t = dData2[2];
		E_ = dData2[3];
		double G = dData2[4];
		double fc = dData2[5];
		double  mu = dData2[6];
		double  tau0 = dData2[7];
		double Gc = 1.0 + dData2[8];
		double dropDrift = dData2[9];
		double muR = dData2[10];

		Ltfc = std::abs(fc*b*t);

		// true for stronger, true for elastic, true for crushing
		theSectionI = new NoTensionSection3d(0, E_, G, t, b, -1.0, fc, 5, true, false, false);  	
		theSectionE = new NoTensionSection3d(0, E_, G, t, b, -1.0, fc, 5, false, true, false);
		theSectionJ = new NoTensionSection3d(0, E_, G, t, b, -1.0, fc, 5, true, false, false);

		theShearModelOOP = new DamageShearInterface(0, E_*b*t, G/(h)*5./6.*b*t, tau0*b*t, mu, muR, Gc, dropDrift*h, true);
		theShearModel    = new DamageShearInterface(0, E_*b*t, G/(h)*5./6.*b*t, tau0*b*t, mu, muR, Gc, dropDrift*h, false);									
     } 



	else if (strcmp(inputStructure, "-spandrelGeneric") == 0) {
		// input for piers, with generic shear model to include in the end. Needs to substitute in the long term the -pier definition
		// left for back compatibility
		// arguments: h, b, t, E, G, fc "shearModel" <variable shear model parameters>

		double dData2[6];
		numData = 6;
		if (OPS_GetDoubleInput(&numData, &dData2[0]) < 0) {
			opserr << "WARNING Macroelement3d (tag: " << iData[0] << "):'-pierGeneric' input structure incorrect, invalid double input(s)." << endln;
			opserr << "Required structure : h, b, t, E, G, fc shearModel <variable shear model parameters> <-flags>" << endln;
			return 0;
		}

		intLength(0) = 1. / 6.;
		intLength(1) = 2. / 3.;
		intLength(2) = 1. / 6.;

		h = dData2[0];

		b = dData2[1];
		t = dData2[2];
		E_ = dData2[3];
		double G = dData2[4];
		double fc = dData2[5];

		Ltfc = abs(fc*b*t);

		// true for stronger, true for elastic, true for crushing
		theSectionI = new NoTensionSection3d(0, E_, G, t, b, -1.0, fc, 5, true, false, false);
		theSectionE = new NoTensionSection3d(0, E_, G, t, b, -1.0, fc, 5, false, true, false);
		theSectionJ = new NoTensionSection3d(0, E_, G, t, b, -1.0, fc, 5, true, false, false);

		const char* type = NULL;
		type = OPS_GetString();

		if (strcmp(type, "-MohrCoulomb") == 0) {
			double dataShear[5];
			numData = 5;

			if (OPS_GetDoubleInput(&numData, &dataShear[0]) < 0) {
				opserr << "WARNING Macroelement3d (tag: " << iData[0] << "):unable to read shear model input(s)." << endln;
				opserr << "Required structure: -MohrCoulomb mu c Gc dropDrift alpha <-flags>" << endln;
				return 0;
			}

			Vector props(4);
			props(0) = b;
			props(1) = t;
			props(2) = dataShear[1];
			props(3) = dataShear[0];

			theShearModel = new GenericDamagePlasticityShear(0, G / (h)*5. / 6.*b*t, 1, props, 1.0 + dataShear[2], dataShear[3] * h, dataShear[4], false, 0.0);
			theShearModelOOP = new GenericDamagePlasticityShear(0, G / (h)*5. / 6.*b*t, 1, props, 1.0 + dataShear[2], dataShear[3] * h, dataShear[4], true, 0.0);
		}
		else if (strcmp(type, "-TurnsekCacovic") == 0) {
			double dataShear[5];
			numData = 5;

			if (OPS_GetDoubleInput(&numData, &dataShear[0]) < 0) {
				opserr << "WARNING Macroelement3d (tag: " << iData[0] << "):unable to read shear model input(s)." << endln;
				opserr << "Required structure: -MohrCoulomb ft b Gc dropDrift alpha <-flags>" << endln;
				return 0;
			}

			Vector props(4);
			props(0) = b;
			props(1) = t;
			props(2) = dataShear[1];
			props(3) = dataShear[0];

			theShearModel = new GenericDamagePlasticityShear(0, G / (h)*5. / 6.*b*t, 2, props, 1.0 + dataShear[2], dataShear[3] * h, dataShear[4], false, 0.0);
			theShearModelOOP = new GenericDamagePlasticityShear(0, G / (h)*5. / 6.*b*t, 2, props, 1.0 + dataShear[2], dataShear[3] * h, dataShear[4], true, 0.0);

		}
		else if ((strcmp(type, "-MohrCoulombCompressed") == 0) || (strcmp(type, "-Eurocode") == 0)) {
			// [L, t, c, mu, 0.065fm]
			double dataShear[6];
			numData = 6;

			if (OPS_GetDoubleInput(&numData, &dataShear[0]) < 0) {
				opserr << "WARNING Macroelement3d (tag: " << iData[0] << "):unable to read shear model input(s)." << endln;
				opserr << "Required structure: -MohrCoulomb mu c 0.065fm Gc dropDrift alpha <-flags>" << endln;
				return 0;
			}

			Vector props(5);
			props(0) = b;
			props(1) = t;
			props(2) = dataShear[1];
			props(3) = dataShear[0];
			props(4) = dataShear[2];

			theShearModel = new GenericDamagePlasticityShear(0, G / (h)*5. / 6.*b*t, 3, props, 1.0 + dataShear[3], dataShear[4] * h, dataShear[5], false, 0.0);
			theShearModelOOP = new GenericDamagePlasticityShear(0, G / (h)*5. / 6.*b*t, 3, props, 1.0 + dataShear[3], dataShear[4] * h, dataShear[5], true, 0.0);

		}

	}




	
	else if (strcmp(inputStructure,"-gable") == 0)  {

        // standard input for gables, with my shear model
		// arguments: h, b, t, E, G, fc, mu, c, Gc, dropDrift, muR 
		gable = true;

		double dData2[11];
	    numData = 11;
		if (OPS_GetDoubleInput(&numData,&dData2[0]) < 0) {
			opserr << "WARNING Macroelement3d (tag: " << iData[0] << "):'-gable' input structure incorrect, invalid double input(s)." << endln;
			opserr << "Required structure : h, b, t, E, G, fc, mu, c, Gc, dropDrift, muR <-flags>" << endln;
			return 0;
		}

		intLength(0) = 1./6.;
		intLength(1) = 2./3.;
		intLength(2) = 1./6.;

		h = dData2[0];
		
		b = dData2[1];
		t = dData2[2];
		E_ = dData2[3];
		double G = dData2[4];
		double fc = dData2[5];
		double  mu = dData2[6];
		double  tau0 = dData2[7];
		double Gc = 1.0 + dData2[8];
		double dropDrift = dData2[9];
		double muR = dData2[10];

		Ltfc = std::abs(fc*0.5*b*t);

		// true for stronger, true for elastic, true for crushing
		theSectionI = new NoTensionSection3d(0, E_, 0.5*G, t,     b, -1.0, fc, 11, false, false, true);   
		theSectionE = new NoTensionSection3d(0, E_,     G, t, 0.5*b, -1.0, fc, 11, false, false, true);
		theSectionJ = new NoTensionSection3d(0, E_,     G, t, 0.1*b, -1.0, fc, 11, false, false, true);

		theShearModelOOP = new DamageShearInterface(0, E_*0.5*b*t, G/(h)*5./6.*0.5*b*t, tau0*0.5*b*t, mu, muR, Gc, dropDrift*h, true);
		theShearModel    = new DamageShearInterface(0, E_*0.5*b*t, G/(h)*5./6.*0.5*b*t, tau0*0.5*b*t, mu, muR, Gc, dropDrift*h, true);	
     }

	else if (strcmp(inputStructure, "-gableGeneric") == 0) {

		// standard input for gables, with my shear model
		// arguments: h, b, t, E, G, fc, mu, c, Gc, dropDrift, muR 
		gable = true;

		double dData2[6];
		numData = 6;
		if (OPS_GetDoubleInput(&numData, &dData2[0]) < 0) {
			opserr << "WARNING Macroelement3d (tag: " << iData[0] << "):'-gable' input structure incorrect, invalid double input(s)." << endln;
			opserr << "Required structure : h, b, t, E, G, fc, mu, c, Gc, dropDrift, muR <-flags>" << endln;
			return 0;
		}

		intLength(0) = 1. / 6.;
		intLength(1) = 2. / 3.;
		intLength(2) = 1. / 6.;

		h = dData2[0];

		b = dData2[1];
		t = dData2[2];
		E_ = dData2[3];
		double G = dData2[4];
		double fc = dData2[5];

		Ltfc = abs(fc*b*t);

		// true for stronger, true for elastic, true for crushing
		theSectionI = new NoTensionSection3d(0, E_, G, t, b,     -1.0, fc, 11, false, false, true);
		theSectionE = new NoTensionSection3d(0, E_, G, t, 0.5*b, -1.0, fc, 11, false, false, true);
		theSectionJ = new NoTensionSection3d(0, E_, G, t, 0.1*b, -1.0, fc, 11, false, false, true);

		const char* type = NULL;
		type = OPS_GetString();

		if (strcmp(type, "-MohrCoulomb") == 0) {
			double dataShear[5];
			numData = 5;

			if (OPS_GetDoubleInput(&numData, &dataShear[0]) < 0) {
				opserr << "WARNING Macroelement3d (tag: " << iData[0] << "):unable to read shear model input(s)." << endln;
				opserr << "Required structure: -MohrCoulomb mu c Gc dropDrift alpha <-flags>" << endln;
				return 0;
			}

			Vector props(4);
			props(0) = b;
			props(1) = t;
			props(2) = dataShear[1];
			props(3) = dataShear[0];

			theShearModel = new GenericDamagePlasticityShear(0,    G / (h)*5. / 6.*b*t*0.5, 1, props, 1.0 + dataShear[2], dataShear[3] * h, dataShear[4], false, 0.0);
			theShearModelOOP = new GenericDamagePlasticityShear(0, G / (h)*5. / 6.*b*t*0.5, 1, props, 1.0 + dataShear[2], dataShear[3] * h, dataShear[4], true, 0.0);
		}
		else if (strcmp(type, "-TurnsekCacovic") == 0) {
			double dataShear[5];
			numData = 5;

			if (OPS_GetDoubleInput(&numData, &dataShear[0]) < 0) {
				opserr << "WARNING Macroelement3d (tag: " << iData[0] << "):unable to read shear model input(s)." << endln;
				opserr << "Required structure: -MohrCoulomb ft b Gc dropDrift alpha <-flags>" << endln;
				return 0;
			}

			Vector props(4);
			props(0) = b;
			props(1) = t;
			props(2) = dataShear[1];
			props(3) = dataShear[0];

			theShearModel = new GenericDamagePlasticityShear(0,    G / (h)*5. / 6.*b*t*0.5, 2, props, 1.0 + dataShear[2], dataShear[3] * h, dataShear[4], false, 0.0);
			theShearModelOOP = new GenericDamagePlasticityShear(0, G / (h)*5. / 6.*b*t*0.5, 2, props, 1.0 + dataShear[2], dataShear[3] * h, dataShear[4], true, 0.0);

		}
		else if ((strcmp(type, "-MohrCoulombCompressed") == 0) || (strcmp(type, "-Eurocode") == 0)) {
			// [L, t, c, mu, 0.065fm]
			double dataShear[6];
			numData = 6;

			if (OPS_GetDoubleInput(&numData, &dataShear[0]) < 0) {
				opserr << "WARNING Macroelement3d (tag: " << iData[0] << "):unable to read shear model input(s)." << endln;
				opserr << "Required structure: -MohrCoulomb mu c 0.065fm Gc dropDrift alpha <-flags>" << endln;
				return 0;
			}

			Vector props(5);
			props(0) = b;
			props(1) = t;
			props(2) = dataShear[1];
			props(3) = dataShear[0];
			props(4) = dataShear[2];

			theShearModel = new GenericDamagePlasticityShear(0, G / (h)*5. / 6.*b*t*0.5, 3, props, 1.0 + dataShear[3], dataShear[4] * h, dataShear[5], false, 0.0);
			theShearModelOOP = new GenericDamagePlasticityShear(0, G / (h)*5. / 6.*b*t*0.5, 3, props, 1.0 + dataShear[3], dataShear[4] * h, dataShear[5], true, 0.0);

		}
	}
	

	
	else {
		// standard input, opensees like. Custom sectional models have to be created before, and their tags are passed to the macrolement.
		// arguments: sectionTag, sectionTagE, shearModelTag, h, E 
		// inputs: 

		intLength(0) = 1. / 6.;
		intLength(1) = 2. / 3.;
		intLength(2) = 1. / 6.;

		int iData2[5];
		int numData = 5;
		if (OPS_GetIntInput(&numData,&iData2[0]) < 0) {
			opserr<<"WARNING: invalid integer inputs, sectional models. \nRequired structure: sectionItag, sectionEtag, sectionJtag, shearModelIP, shearModelOOP, h, E <-flags>\n";
			return 0;
		}

		numData = 1;
		if (OPS_GetDoubleInput(&numData,&h) < 0) {
			opserr<<"WARNING: invalid double inputs, sectional models\n";
			return 0;
		}

		numData = 1;
		if (OPS_GetDoubleInput(&numData,&E_) < 0) {
			opserr<<"WARNING: invalid double inputs, sectional models\n";
			return 0;
		}

		// create copies of the sectional models
		int secTag = iData2[0];
	    theSectionI = OPS_GetSectionForceDeformation(secTag);
		if (theSectionI == 0) {
			opserr<<"SectionI: section model tagged " << secTag << " not found." << endln;
			return 0;
		}

		secTag = iData2[1];
		theSectionE = OPS_GetSectionForceDeformation(secTag);
		if (theSectionE == 0) {
			opserr<<"SectionE: section model tagged " << secTag << " not found." << endln;
			return 0;
		}

		secTag = iData2[2];
	    theSectionJ = OPS_GetSectionForceDeformation(secTag);
		if (theSectionJ == 0) {
			opserr<<"SectionJ: section model tagged " << secTag << " not found." << endln;
			return 0;
		}

		int matTag = iData2[3];
		theShearModel = OPS_GetNDMaterial(matTag);
		if (theShearModel == 0) {
			opserr<<"Shear model IP: NDMaterial " << matTag << " not found." << endln;
			return 0;
		}

		int matTagOOP = iData2[4];
		theShearModelOOP = OPS_GetNDMaterial(matTagOOP);
		if (theShearModelOOP == 0) {
			opserr<<"Shear model OOP: NDMaterial " << matTag << " not found." << endln;
			return 0;
		}

		t=0.0;

	 }
	 
	// optional inputs

    double mass = 0.0;
	double weights[3];
	double massGlobal[3];
    int cmass = 0;
    int PDelta = 0;
	int dampingModel = 0;
	Vector massDir(3);
	// default: mass applied along all global directions
	massDir(0) = massDir(1) = massDir(2) = 1.;

	bool read_Ltfc=false;
	if (Ltfc==0) { read_Ltfc=true;}
	bool alreadReadType = false;
	const char* type = NULL;


    while(OPS_GetNumRemainingInputArgs() > 0) {
	if (!alreadReadType)	type = OPS_GetString();

	if ((strcmp(type,"-cMass") == 0) || (strcmp(type, "-cmass") == 0)) {
	    cmass = 1;
		alreadReadType = false;
	} 
	else if (strcmp(type,"-mass") == 0) {
	    if (OPS_GetNumRemainingInputArgs() > 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData,&mass) < 0) {
				opserr<<"WARNING: Macroelement " << iData[0] <<", invalid mass\n";
				return 0;
			}
			alreadReadType = false;
	    }
	} 
	else if (strcmp(type,"-AxialCollapseRatio") == 0) {
	    if (OPS_GetNumRemainingInputArgs() > 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData,&alphaAC_HC) < 0) {
				opserr<<"WARNING: Macroelement " << iData[0] <<", invalid AxialCollapseRatio\n";
			}
			if (alphaAC_HC < 1) {
				opserr<<"WARNING: Macroelement " << iData[0] <<", AxialCollapseRatio must be bigger than 1.0 (set equal to 1.0)\n";
				alphaAC_HC = 1.0;
			}
			alreadReadType = false;
	    }
	}

	else if (strcmp(type,"-failureFactor") == 0) {
	    if (OPS_GetNumRemainingInputArgs() > 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData,&failureFactorF) < 0) {
				opserr<<"WARNING: Macroelement " << iData[0] <<", invalid failureFactor\n";
			}
			if ((failureFactorF > 1.0) || (failureFactorF<0.0)) {
				opserr<<"WARNING: Macroelement " << iData[0] <<", failureFactor must be between 0.0 and 1.0 (set equal to 0.001)\n";
				failureFactorF = 0.001;
			}
			failureFactorS = failureFactorF;
			alreadReadType = false;
	    }
	}

	else if (strcmp(type,"-failureFactorFlexure") == 0) {
	    if (OPS_GetNumRemainingInputArgs() > 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData,&failureFactorF) < 0) {
				opserr<<"WARNING: Macroelement " << iData[0] <<", invalid failureFactor\n";
			}
			if ((failureFactorF > 1.0) || (failureFactorF<0.0)) {
				opserr<<"WARNING: Macroelement " << iData[0] <<", failureFactorFlexure must be between 0.0 and 1.0 (set equal to 0.001)\n";
				failureFactorF = 0.001;
			}
			alreadReadType = false;
	    }
	}

     else if (strcmp(type,"-failureFactorShear") == 0) {
	    if (OPS_GetNumRemainingInputArgs() > 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData,&failureFactorS) < 0) {
				opserr<<"WARNING: Macroelement " << iData[0] <<", invalid failureFactor\n";
			}
			if ((failureFactorS > 1.0) || (failureFactorS<0.0)) {
				opserr<<"WARNING: Macroelement " << iData[0] <<", failureFactorShear must be between 0.0 and 1.0 (set equal to 0.001)\n";
				failureFactorS = 0.001;
			}
			alreadReadType = false;
	    }
	}


	else if (strcmp(type,"-density") == 0) {
	    if (OPS_GetNumRemainingInputArgs() > 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData,&mass) < 0) {
				opserr << "WARNING: Macroelement " << iData[0] <<", invalid mass.\n";
				return 0;
			}
			mass *= b*t; 
			if (b*t<DBL_EPSILON)
				opserr << "WARNING: Macroelement " << iData[0] <<" has zero area, use -mass flag instead of -density.\n";
	    }
		alreadReadType = false;
	} 
	else if (strcmp(type,"-PDelta") == 0  || strcmp(type,"-pDelta") == 0) {
             PDelta = 1;
			 alreadReadType = false;
	       }

	else if (strcmp(type, "-secantDamping") == 0 || strcmp(type, "-SecantDamping") == 0) {
		dampingModel = 1;
		alreadReadType = false;
	}

	else if (strcmp(type, "-secantDampingCycle") == 0 || strcmp(type, "-SecantDampingCycle") == 0) {
		dampingModel = 2;
		alreadReadType = false;
	}
	       
	else  if (strcmp(type,"-intWeights") == 0) {
              if(OPS_GetNumRemainingInputArgs() > 2) {
				  numData = 3;
				  if (OPS_GetDoubleInput(&numData,&weights[0]) < 0) {
					  opserr << "WARNING: Macroelement " << iData[0] <<", invalid integration weights, standard weights are applied.\n";
				  } else {
					  intLength(0) = weights[0];
					  intLength(1) = weights[1];
					  intLength(2) = weights[2];
					  // maybe check that they sum to 1
					  double sum = 0.;
					  sum = intLength(0) + intLength(1) + intLength(2);
					  if (std::abs(sum-1.0) > 0.01)
						  opserr << "WARNING: Macroelement " << iData[0] <<", specified integration weights do not sum exactly to 1 (sum=" << sum << ").\n";
				  }
			  }
			  alreadReadType = false;
	       }
		else  if (strcmp(type,"-massDir") == 0) {
              if (OPS_GetNumRemainingInputArgs() > 2) {
				  numData = 3;
				  if (OPS_GetDoubleInput(&numData,&massGlobal[0]) < 0) {
					  opserr << "WARNING: Macroelement " << iData[0] <<", invalid mass directions, mass applied in all global directions.\n";
				  } else {
					  massDir(0) = massGlobal[0];
					  massDir(1) = massGlobal[1];
					  massDir(2) = massGlobal[2];
					  // check that they are not bigger than 1, but allow the choice anyway
					  if ( massDir(0)>1.0 || massDir(0)<0.0 || massDir(1)>1.0 || massDir(1)<0.0 || massDir(2)>1.0 || massDir(2)<0.0 )
						  opserr << "WARNING: Macroelement " << iData[0] <<", specified mass directions are negative or bigger than one.\n";
				  }
			  }
			  alreadReadType = false;
	       }
		//drift flexure
		else  if (strcmp(type,"-driftFlexure") == 0) {

			int otherInputs=OPS_GetNumRemainingInputArgs();
			std::vector<double> parsed;

			bool goOn = true;
			int i=0;
			//int numData = 1;
			double singleInput;

			const char* stringRead = NULL;
			char* pEnd;
			
			while (goOn && i<otherInputs) {
				stringRead = OPS_GetString();
				
				char *err;
				double dnotUsed = strtod(stringRead, &err);
				if (*err == 0) { 
					singleInput = strtod(stringRead, &pEnd);
					parsed.push_back(singleInput);
					i += 1;
				} else if (!isspace((unsigned char)*err)) { goOn = false; }

			}

			if (read_Ltfc) {
				Ltfc = parsed[parsed.size()-1];
				parsed.pop_back();
			}

			betaShearSpanF = parsed[parsed.size()-1];
			parsed.pop_back();

			// if only one value is given, constant drfit model
			if (parsed.size() ==1) {
				double constantDrift  = parsed[0];
				parsed[0] = 0.0;
				parsed.push_back(constantDrift);
				parsed.push_back(1.0);
				parsed.push_back(constantDrift);
			}

			// check that the law is defined between ALR=0 and ALR = 1, if not extend it (constant up to the next defined point, or to the end)
			if (parsed[0]>0.0) {
				parsed.insert(parsed.begin(), parsed[1]);
				parsed.insert(parsed.begin(), 0.0);
			}

			if (parsed[parsed.size()-2]<1.0) {
				parsed.push_back(1.0),
				parsed.push_back(parsed[parsed.size()-2]);
			}

			int nPoints = (int) (parsed.size()/2);
			driftModelF_ALR.resize(nPoints);
			driftModelF.resize(nPoints);

			for (int i=0; i<nPoints; i++) {
				driftModelF_ALR(i) = parsed[2*i];
				driftModelF(i)     = parsed[2*i+1];
			}

			alreadReadType = true;
			type = stringRead;

		}

		//drift shear
		else  if (strcmp(type,"-driftShear") == 0) {
            int otherInputs=OPS_GetNumRemainingInputArgs();
			std::vector<double> parsedS;

			bool goOn = true;
			int i=0;
			//		int numData = 1;
			double singleInput;

			const char* stringRead = NULL;
			char* pEnd;
			
			while (goOn && i<otherInputs) {
				stringRead = OPS_GetString();

				char *err;
				double dnotUsed = strtod(stringRead, &err);
				if (*err == 0) {
					singleInput = strtod(stringRead, &pEnd);
					parsedS.push_back(singleInput);
					i += 1;
				}
				else if (!isspace((unsigned char)*err)) { goOn = false; }
			}

			if (read_Ltfc) {
				Ltfc = parsedS[parsedS.size()-1];
				parsedS.pop_back();
			}

			betaShearSpanS = parsedS[parsedS.size()-1];
			parsedS.pop_back();

			// if only one value is given, constant drfit model
			if (parsedS.size() ==1) {
				double constantDrift  = parsedS[0];
				parsedS[0] = 0.0;
				parsedS.push_back(constantDrift);
				parsedS.push_back(1.0);
				parsedS.push_back(constantDrift);
			}

			// check that the law is defined between ALR=0 and ALR = 1, if not extend it (constant up to the next defined point, or to the end)
			if (parsedS[0]>0.0) {
				parsedS.insert(parsedS.begin(), parsedS[1]);
				parsedS.insert(parsedS.begin(), 0.0);
			}

			if (parsedS[parsedS.size()-2]<1.0) {
				parsedS.push_back(1.0),
				parsedS.push_back(parsedS[parsedS.size()-2]);
			}

			
			int nPoints = (int) (parsedS.size()/2);
			driftModelS_ALR.resize(nPoints);
			driftModelS.resize(nPoints);

			for (int i=0; i<nPoints; i++) {
				driftModelS_ALR(i) = parsedS[2*i];
				driftModelS(i)     = parsedS[2*i+1];
			}

			alreadReadType = true;
			type = stringRead;

			//opserr << "driftS_ALR = " << driftModelS_ALR;
			//opserr << "driftS = " << driftModelS;
	       }
	
	}  // end while other inputs


	 Element *theEle =  new Macroelement3d(iData[0],iData[1],iData[2],iData[3],theSectionI,theSectionE,theSectionJ,theShearModel,theShearModelOOP,h,E_,
		                                    driftModelF,driftModelF_ALR, driftModelS, driftModelS_ALR, Ltfc, alphaAC_HC, betaShearSpanF, betaShearSpanS, failureFactorF, failureFactorS,
											axis,oop, intLength, intLengthMasses, massDir, PDelta,mass,cmass, gable, dampingModel);

	 if ((strcmp(inputStructure,"-tremuri")==0) || 
		 (strcmp(inputStructure,"-tremuriIP")==0)  || 
		 (strcmp(inputStructure,"-pier")==0) || 
		 (strcmp(inputStructure,"-spandrel")==0) || 
		 (strcmp(inputStructure,"-gable")==0) )  {

	    if (theSectionI!=NULL)  delete theSectionI;
	    if (theSectionE!=NULL)  delete theSectionE;
	    if (theSectionJ!=NULL)  delete theSectionJ;
	    if (theShearModel!=NULL)  delete theShearModel;
	    if (theShearModelOOP!=NULL)  delete theShearModelOOP;
	 }
	 
     return theEle;
}

Macroelement3d::Macroelement3d(
    int tag, int nd1, int nd2, int ndE, SectionForceDeformation *sI,
    SectionForceDeformation *sE, SectionForceDeformation *sJ,
    NDMaterial *shearModel, NDMaterial *shearModelOOP, double h, double E_,
    Vector driftF, Vector driftF_ALR, Vector driftS, Vector driftS_ALR,
    double Ltfc, double alphaAC_HC, double betaShearSpanF,
    double betaShearSpanS, double failureFactorF, double failureFactorS,
    Vector axis, Vector oop, Vector _intLength, Vector _intLengthMasses,
    Vector massDir, int PDelta, double rho, int cm, double _isGable, int dampingModel)
    : Element(tag, 0), numSections(3), numShearModels(2), theSections(0),
      theShearModel(0), connectedExternalNodes(3), GammaC(12, 18), Tgl(18, 18),
      Tgl6(6, 6), Q(18), q(12), uBasic(12), uBasicCommitted(12), rho(rho),
      cMass(cm), massglobalDir(massDir), PDelta(PDelta), deltaW1(0.),
      deltaV1(0.), deltaW3(0.), deltaV3(0.), parameterID(0),
      nodeIInitialDisp(0), nodeJInitialDisp(0), nodeIOffset(0), nodeJOffset(0),
      xAxis(3), yAxis(3), zAxis(3), intLength(_intLength),
      intLengthMasses(_intLengthMasses), L(h / 2.0), E(E_), driftF(0.),
      driftS(0.), Ltfc(Ltfc), betaShearSpanF(betaShearSpanF),
      betaShearSpanS(betaShearSpanS), limitDriftF(driftF), limitDriftS(driftS),
      limitDriftF_ALR(driftF_ALR), limitDriftS_ALR(driftS_ALR),
      alphaAC_HC(alphaAC_HC), failedF(false), failedS(false),
      failedFcommitted(false), failedScommitted(false), collapsedF(false),
      collapsedS(false), collapsedFcommitted(false), collapsedScommitted(false),
      failureFactorF(failureFactorF), failureFactorS(failureFactorS),
      isGable(_isGable), wx(0.0), wy(0.0), wz(0.0), ductilityDemand(0.0), ductilityDemandCommitted(1.0), 
	  ductilityDemandPeak(1.0), ductilityDemandCycle(1.), triggerDuctility(true), dampingModel(dampingModel) {

  // Allocate arrays of pointers to SectionForceDeformations
  theSections = new SectionForceDeformation *[numSections];
  if (theSections == 0) {
    opserr << "Macroelement3d::Macroelement3d - failed to allocate section model pointer\n";
    exit(-1);
  }
  
  theSections[0] = sI->getCopy();       if (theSections[0] == 0) { opserr << "Macroelement3d::Macroelement3d -- failed to get a copy of section model\n";  exit(-1);  }
  theSections[1] = sE->getCopy();       if (theSections[1] == 0) { opserr << "Macroelement3d::Macroelement3d -- failed to get a copy of section model\n";  exit(-1);  }
  theSections[2] = sJ->getCopy();       if (theSections[2] == 0) { opserr << "Macroelement3d::Macroelement3d -- failed to get a copy of section model\n";  exit(-1);  }

  CommittedDeformation[0] = (theSections[0]->getSectionDeformation());
  CommittedDeformation[1] = (theSections[1]->getSectionDeformation());
  CommittedDeformation[2] = (theSections[2]->getSectionDeformation());

  theShearModel = new NDMaterial *[numShearModels];
  theShearModel[0] = shearModel->getCopy("BeamFiber");
  theShearModel[1] = shearModelOOP->getCopy("BeamFiber");
  if (theShearModel == 0) {
    opserr << "Macroelement3d::Macroelement3d - failed to allocate shear model pointer\n";
    exit(-1);
  }
 
  // Set connected external node IDs
  connectedExternalNodes(0) = nd1;
  connectedExternalNodes(1) = nd2;
  connectedExternalNodes(2) = ndE;

  theNodes[0] = 0;
  theNodes[1] = 0;
  theNodes[2] = 0;

  for (int i=0; i<12; i++) {
	  q0[i] = 0.0;
      p0[i] = 0.0;
  }


  // define orientation
  Vector vAxis(3);
  for (int i=0; i<3; i++) {
      xAxis(i) = axis(i);
      vAxis(i) = oop(i);
  }

  // make xAxis unitary
  double xnorm = xAxis.Norm();
  if (xnorm == 0) {
        opserr << "Macroelement3d ("<< tag << "). Defined between two coincident nodes, check geometry." << endln;
        exit(-1);
  }
  
  // y = v cross x
  yAxis(0) = vAxis(1)*xAxis(2) - vAxis(2)*xAxis(1);
  yAxis(1) = vAxis(2)*xAxis(0) - vAxis(0)*xAxis(2);
  yAxis(2) = vAxis(0)*xAxis(1) - vAxis(1)*xAxis(0);

  double ynorm = yAxis.Norm();
  if (ynorm == 0) {
        opserr << "Macroelement3d ("<< tag << ").The vector v that defines plane xz is parallel to x axis." << endln;
        exit(-1);
  }
    
  yAxis /= ynorm;

  // z = x cross y. Already unitary
  zAxis(0) = xAxis(1)*yAxis(2) - xAxis(2)*yAxis(1);
  zAxis(1) = xAxis(2)*yAxis(0) - xAxis(0)*yAxis(2);
  zAxis(2) = xAxis(0)*yAxis(1) - xAxis(1)*yAxis(0);

  // Fill in transformation matrix
  R[0][0] = xAxis(0);
  R[0][1] = xAxis(1);
  R[0][2] = xAxis(2);
    
  R[1][0] = yAxis(0);
  R[1][1] = yAxis(1);
  R[1][2] = yAxis(2);
    
  R[2][0] = zAxis(0);
  R[2][1] = zAxis(1);
  R[2][2] = zAxis(2);

  // setup transformation matrix from global to local
    Tgl.Zero();
    Tgl(0,0) = Tgl(3,3) = Tgl(6,6) = Tgl(9,9)   = Tgl(12,12) = Tgl(15,15)  = R[0][0];
    Tgl(0,1) = Tgl(3,4) = Tgl(6,7) = Tgl(9,10)  = Tgl(12,13) = Tgl(15,16)  = R[0][1];
    Tgl(0,2) = Tgl(3,5) = Tgl(6,8) = Tgl(9,11)  = Tgl(12,14) = Tgl(15,17)  = R[0][2];
    Tgl(1,0) = Tgl(4,3) = Tgl(7,6) = Tgl(10,9)  = Tgl(13,12) = Tgl(16,15)  = R[1][0];
    Tgl(1,1) = Tgl(4,4) = Tgl(7,7) = Tgl(10,10) = Tgl(13,13) = Tgl(16,16)  = R[1][1];
    Tgl(1,2) = Tgl(4,5) = Tgl(7,8) = Tgl(10,11) = Tgl(13,14) = Tgl(16,17)  = R[1][2];
    Tgl(2,0) = Tgl(5,3) = Tgl(8,6) = Tgl(11,9)  = Tgl(14,12) = Tgl(17,15)  = R[2][0];
    Tgl(2,1) = Tgl(5,4) = Tgl(8,7) = Tgl(11,10) = Tgl(14,13) = Tgl(17,16)  = R[2][1];
    Tgl(2,2) = Tgl(5,5) = Tgl(8,8) = Tgl(11,11) = Tgl(14,14) = Tgl(17,17)  = R[2][2];
    

    Tgl6(0,0) = Tgl6(3,3)  = R[0][0];
    Tgl6(0,1) = Tgl6(3,4)  = R[0][1];
    Tgl6(0,2) = Tgl6(3,5)  = R[0][2];
    Tgl6(1,0) = Tgl6(4,3)  = R[1][0];
    Tgl6(1,1) = Tgl6(4,4)  = R[1][1];
    Tgl6(1,2) = Tgl6(4,5)  = R[1][2];
    Tgl6(2,0) = Tgl6(5,3)  = R[2][0];
    Tgl6(2,1) = Tgl6(5,4)  = R[2][1];
    Tgl6(2,2) = Tgl6(5,5)  = R[2][2];

    
	this->intLength *= 2*L;
	this->intLengthMasses *= 2*L;

}

Macroelement3d::Macroelement3d()
    : Element(0, 0), numSections(0), numShearModels(0), theSections(0),
      theShearModel(0), connectedExternalNodes(3), Q(18), q(12), uBasic(12),
      uBasicCommitted(12), Tgl(18, 18), Tgl6(6, 6), GammaC(12, 18), rho(0.0),
      cMass(0), parameterID(0), E(0.), intLength(3), intLengthMasses(3),
      driftF(0.), driftS(0.), limitDriftF(0), limitDriftS(0),
      limitDriftF_ALR(0), limitDriftS_ALR(0), alphaAC_HC(0.0),
      betaShearSpanF(0.0), betaShearSpanS(0.0), failureFactorF(0.0),
      failureFactorS(0.0), failedF(false), failedS(false),
      failedFcommitted(false), failedScommitted(false), collapsedF(false),
      collapsedS(false), collapsedFcommitted(false), collapsedScommitted(false),
      nodeIInitialDisp(0), nodeJInitialDisp(0), nodeIOffset(0), nodeJOffset(0),
      xAxis(3), yAxis(3), zAxis(3), L(0.0), PDelta(0), deltaW1(0.), deltaV1(0.),
      deltaW3(0.), deltaV3(0.), ductilityDemand(0.0), ductilityDemandCommitted(1.0), 
	  ductilityDemandPeak(1.0), ductilityDemandCycle(1.0), triggerDuctility(true), dampingModel(0) {
  for (int i=0; i<12; i++) {
	  q0[i] = 0.0;
      p0[i] = 0.0;
  }

  for (int i=0; i<3; i++) {
      R[i][0] = 0.0;
      R[i][1] = 0.0;
      R[i][2] = 0.0;
   }

  theNodes[0] = 0;
  theNodes[1] = 0;
  theNodes[2] = 0;
}

Macroelement3d::~Macroelement3d()
{    
  for (int i = 0; i < numSections; i++) {
    if (theSections[i])
      delete theSections[i];
  }

  for (int i = 0; i < 2; i++) {
    if (theShearModel[i])
      delete theShearModel[i];
  }
  
  // Delete the array of pointers to SectionForceDeformation pointer arrays and to the Shear model
  if (theSections)
    delete [] theSections;

  if (theShearModel)
    delete [] theShearModel;
}

int
Macroelement3d::getNumExternalNodes() const {
    return 3;
}

const ID&
Macroelement3d::getExternalNodes() {
    return connectedExternalNodes;
}

Node **
Macroelement3d::getNodePtrs() {
    return theNodes;
}

int
Macroelement3d::getNumDOF() {
    return 18;
}

void
Macroelement3d::setDomain(Domain *theDomain) {
	// Check Domain is not null - invoked when object removed from a domain

	if (theDomain == 0) {
		theNodes[0] = 0;
		theNodes[1] = 0;
        theNodes[2] = 0;
	return;
    }

    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);
    int NdE = connectedExternalNodes(2);

    theNodes[0] = theDomain->getNode(Nd1);
    theNodes[1] = theDomain->getNode(Nd2);
    theNodes[2] = theDomain->getNode(NdE);

	

    if (theNodes[0] == 0 || theNodes[1] == 0) {
		return;
    }

    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();
    int dofNdE = theNodes[2]->getNumberDOF();
    
    if (dofNd1 != 6 || dofNd2 != 6 || dofNdE != 6) {
		return;
    }


    // element projection
    static Vector dx(3);
    
    const Vector &ndICoords = theNodes[0]->getCrds();
    const Vector &ndJCoords = theNodes[1]->getCrds();
    
    dx(0) = ndJCoords(0) - ndICoords(0);
    dx(1) = ndJCoords(1) - ndICoords(1);
    dx(2) = ndJCoords(2) - ndICoords(2);



    // see if there is some initial displacements at nodes
    const Vector &nodeIDisp = theNodes[0]->getDisp();
    const Vector &nodeJDisp = theNodes[0]->getDisp();
    bool addInitialDisplacement = false;

    for (int i=0; i<6; i++) {
            if (nodeIDisp(i) != 0.0 || nodeJDisp(i) != 0.0) {
			     addInitialDisplacement = true;
			}
	}

	if (addInitialDisplacement) {
          nodeIInitialDisp = new double [6];
          for (int j=0; j<6; j++) {
               nodeIInitialDisp[j] = nodeIDisp(j);
		  }
            
          nodeJInitialDisp = new double [6];
          for (int j=0; j<6; j++) {
				  nodeJInitialDisp[j] = nodeJDisp(j);
          }
	}

    if (nodeIInitialDisp != 0) {
        dx(0) -= nodeIInitialDisp[0];
        dx(1) -= nodeIInitialDisp[1];
        dx(2) -= nodeIInitialDisp[2];
    }
    
    if (nodeJInitialDisp != 0) {
        dx(0) += nodeJInitialDisp[0];
        dx(1) += nodeJInitialDisp[1];
        dx(2) += nodeJInitialDisp[2];
    }



    // define position
    Vector inode(3);
    Vector jnode(3);
    const Vector &ndECoords = theNodes[2]->getCrds();

    inode = ndECoords - L*xAxis;
    jnode = ndECoords + L*xAxis;

     Vector offset(3);
     offset =  inode - ndICoords;
	 //opserr << "Macroelement " << this->getTag() << ": offset I = " << offset;
	 if (offset.Norm()>DBL_EPSILON) {
            nodeIOffset = new double[3];
            nodeIOffset[0] = offset(0);
            nodeIOffset[1] = offset(1);
            nodeIOffset[2] = offset(2);
      }


     offset = jnode - ndJCoords;
	 //opserr << "Macroelement " << this->getTag() << ": offset J = " << offset;

	 if (offset.Norm()>DBL_EPSILON) {
            nodeJOffset = new double[3];
            nodeJOffset[0] = offset(0);
            nodeJOffset[1] = offset(1);
            nodeJOffset[2] = offset(2);
      }

    
    if (nodeJOffset != 0) {
        dx(0) += nodeJOffset[0];
        dx(1) += nodeJOffset[1];
        dx(2) += nodeJOffset[2];
		
    }
    
    if (nodeIOffset != 0) {
        dx(0) -= nodeIOffset[0];
        dx(1) -= nodeIOffset[1];
        dx(2) -= nodeIOffset[2];
		
    }

	
    
    // calculate the element length
    L = dx.Norm();

	// we call L the height of each of the two halves
	L /= 2.;
   
    this->DomainComponent::setDomain(theDomain);
	this->getIncrementalCompatibilityMatrix(false);
	this->update();
	
}

int
Macroelement3d::commitState() {

    int retVal = 0;

    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
		opserr << "Macroelement3d::commitState () - failed in base class";
    }    

    // Loop over the interfaces and commit the material states
	for (int i = 0; i < numSections; i++) {
		retVal += theSections[i]->commitState();
		CommittedDeformation[i] = (theSections[i]->getSectionDeformation());
	}

    // commit shear model
	for (int i = 0; i < numShearModels; i++)
		retVal += theShearModel[i]->commitState();


	failedFcommitted = failedF;
	failedScommitted = failedS;

	collapsedFcommitted = collapsedF;
	collapsedScommitted = collapsedS;

	// commit basic displacements
	uBasicCommitted = uBasic;

	ductilityDemandCommitted = ductilityDemand;

	// cycle trigger
	if (ductilityDemandPeak < ductilityDemandCommitted && triggerDuctility)  ductilityDemandPeak = ductilityDemandCommitted;
	if (ductilityDemandPeak > ductilityDemandCommitted && !triggerDuctility)  ductilityDemandPeak = ductilityDemandCommitted;
	if (ductilityDemandCycle < ductilityDemandPeak) ductilityDemandCycle = ductilityDemandPeak;
	if (ductilityDemandCommitted < 0.95*ductilityDemandPeak && triggerDuctility) {
		triggerDuctility = false;
		ductilityDemandPeak = ductilityDemandCommitted;
		ductilityDemandCycle = ductilityDemandPeak;
	}
	if (ductilityDemandCommitted > 1.05*ductilityDemandPeak && !triggerDuctility) {
		triggerDuctility = true;
		ductilityDemandPeak = ductilityDemandCommitted;
	}

	//opserr << "Macroelement " << this->getTag() << ": ductility demand: " << ductilityDemandCommitted << endln;

	return retVal;
}

int
Macroelement3d::revertToLastCommit() {

    int retVal = 0;

    // Loop over the integration points and revert to last committed state
    for (int i = 0; i < numSections; i++)
		retVal += theSections[i]->revertToLastCommit();

    // revert shear model
	for (int i = 0; i < numShearModels; i++)
		retVal += theShearModel[i]->revertToLastCommit();

	failedF = failedFcommitted;
	failedS = failedScommitted;

	collapsedF = collapsedFcommitted;
	collapsedS = collapsedScommitted;

    return retVal;
}

int
Macroelement3d::revertToStart()
{
    int retVal = 0;

    // Loop over the integration points and revert states to start
    for (int i = 0; i < numSections; i++)
		retVal += theSections[i]->revertToStart();
    
	// revert shear model
	for (int i = 0; i < numShearModels; i++)
		retVal += theShearModel[i]->revertToStart();

	failedS = false;
	failedF = false;
	failedScommitted = false;
	failedFcommitted = false;

	collapsedS = false;
	collapsedF = false;
	collapsedScommitted = false;
	collapsedFcommitted = false;

    return retVal;
}

int
Macroelement3d::update(void)
{
  int err = 0;

  const Vector &dispI = theNodes[0]->getTrialDisp();
  const Vector &dispJ = theNodes[1]->getTrialDisp();
  const Vector &dispE = theNodes[2]->getTrialDisp();

  uBasic = this->getBasicDisplacement(dispI, dispJ, dispE);

  double oneOverL = 1./L;
  double M1, M3;

  
    // Set section deformations
	// divide by integration length
  double trialDuctilityDemand = 1.0;
  ductilityDemand = 1.0;
  double OOPStiffness = 1.;
  double secantStiffness = 1.;
  double tangentStiffness = 1.;

	// get ordering of sectional outputs (standard: P, Mz, My, T; a different ordering though can be implemented for some section models)
	// assume that different section models with different orderings can be applied to the three sections
	int order = theSections[0]->getOrder();
	const ID &code0 = theSections[0]->getType();
	Vector e(order);  // allow for longer vectors if the section order is higher (section aggregator, sections with shear deformations)

	int ordering[4];
	int j;
	for (j = 0; j < 4; j++)
		ordering[j] = -1;

	for (j = 0; j < order; j++) {
		switch (code0(j)) {
		case SECTION_RESPONSE_P:     ordering[0] = j;	break;
		case SECTION_RESPONSE_MZ:    ordering[1] = j;   break;
		case SECTION_RESPONSE_MY:    ordering[2] = j;	break;
		case SECTION_RESPONSE_T:     ordering[3] = j;	break;
		default: break;
		}
	}

	// apply to the first section 
	for (int i = 0; i<3; i++)
		if ( ordering[i]>=0)
			e( ordering[i] ) = uBasic(i) / intLength(0);

	if (ordering[3] >= 0)
		e( ordering[3] ) = uBasic(3) / (2.*L);

	err += theSections[0]->setTrialSectionDeformation(e);


	M1 = (theSections[0]->getStressResultant())(ordering[1]);

	if (dampingModel > 0) {
		if (abs(e(ordering[2])) > DBL_EPSILON) {
			OOPStiffness = (theSections[0]->getInitialTangent())(ordering[2], ordering[2]);
			tangentStiffness = ((theSections[0]->getSectionTangent())*(CommittedDeformation[0] - e))[2] / ((CommittedDeformation[0])[2] - e[2]);
			secantStiffness = (theSections[0]->getStressResultant())(ordering[2]) / e(ordering[2]);

			if (secantStiffness > OOPStiffness)   secantStiffness = OOPStiffness;
			if (secantStiffness < tangentStiffness)  secantStiffness = tangentStiffness;
			trialDuctilityDemand = OOPStiffness / secantStiffness;
			if (ductilityDemand < trialDuctilityDemand) {
				ductilityDemand = trialDuctilityDemand;
			}
		}

		if (abs(e(ordering[1])) > DBL_EPSILON) {
			OOPStiffness = (theSections[0]->getInitialTangent())(ordering[1], ordering[1]);
			tangentStiffness = ((theSections[0]->getSectionTangent())*(CommittedDeformation[0] - e))[1] / ((CommittedDeformation[0])[1] - e[1]);
			secantStiffness = (theSections[0]->getStressResultant())(ordering[1]) / e(ordering[1]);

			if (secantStiffness > OOPStiffness)   secantStiffness = OOPStiffness;
			if (secantStiffness < tangentStiffness)  secantStiffness = tangentStiffness;
			trialDuctilityDemand = OOPStiffness / secantStiffness;
			if (ductilityDemand < trialDuctilityDemand) {
				ductilityDemand = trialDuctilityDemand;
			}
		}
	}
	
	//opserr << "; ductility: " << trialDuctilityDemand << endln;


	// second section
	order = theSections[1]->getOrder();
	const ID &code1 = theSections[1]->getType();

	for (j = 0; j < 4; j++)
		ordering[j] = -1;

	for (j = 0; j < order; j++) {
		switch (code1(j)) {
		case SECTION_RESPONSE_P:     ordering[0] = j;	break;
		case SECTION_RESPONSE_MZ:    ordering[1] = j;   break;
		case SECTION_RESPONSE_MY:    ordering[2] = j;	break;
		default: break;
		}
	}

	e.resize(order); // different section orders might be used (unlikely..)
    e.Zero();
    for (int i=0; i<3; i++)
		if (ordering[i]>=0)
			e( ordering[i] ) = uBasic(i+4)   / intLength(1);
    err += theSections[1]->setTrialSectionDeformation(e);

	OOPStiffness = (theSections[1]->getInitialTangent())(ordering[2], ordering[2]);
	
	if (dampingModel > 0) {
		if (abs(e(ordering[2])) > DBL_EPSILON) {
			tangentStiffness = ((theSections[1]->getSectionTangent())*(CommittedDeformation[1] - e))[2] / ((CommittedDeformation[1])[2] - e[2]);
			secantStiffness = (theSections[1]->getStressResultant())(ordering[2]) / e(ordering[2]);
			if (secantStiffness > OOPStiffness)   secantStiffness = OOPStiffness;
			if (secantStiffness < tangentStiffness)  secantStiffness = tangentStiffness;
			trialDuctilityDemand = OOPStiffness / secantStiffness;
			if (ductilityDemand < trialDuctilityDemand) {
				ductilityDemand = trialDuctilityDemand;
			}
		}

		if (abs(e(ordering[1])) > DBL_EPSILON) {
			tangentStiffness = ((theSections[1]->getSectionTangent())*(CommittedDeformation[1] - e))[1] / ((CommittedDeformation[1])[1] - e[1]);
			secantStiffness = (theSections[1]->getStressResultant())(ordering[1]) / e(ordering[1]);
			if (secantStiffness > OOPStiffness)   secantStiffness = OOPStiffness;
			if (secantStiffness < tangentStiffness)  secantStiffness = tangentStiffness;
			trialDuctilityDemand = OOPStiffness / secantStiffness;
			if (ductilityDemand < trialDuctilityDemand) {
				ductilityDemand = trialDuctilityDemand;
			}
		}
	}

	//third section
	order = theSections[2]->getOrder();
	const ID &code2 = theSections[2]->getType();

	for (j = 0; j < 4; j++)
		ordering[j] = -1;

	for (j = 0; j < order; j++) {
		switch (code2(j)) {
		case SECTION_RESPONSE_P:     ordering[0] = j;	break;
		case SECTION_RESPONSE_MZ:    ordering[1] = j;   break;
		case SECTION_RESPONSE_MY:    ordering[2] = j;	break;
		default: break;
		}
	}

	e.resize(order); // different section orders might be used (unlikely..)
    e.Zero();
    for (int i=0; i<3; i++)
		if (ordering[i]>=0)
			e( ordering[i] ) = uBasic(i+7)   / intLength(2);
    err += theSections[2]->setTrialSectionDeformation(e);
	M3 = (theSections[2]->getStressResultant())(ordering[1]);

	OOPStiffness = (theSections[2]->getInitialTangent())(ordering[2], ordering[2]);
	
	if (dampingModel > 0) {
		if (abs(e(ordering[2])) > DBL_EPSILON) {
			tangentStiffness = ((theSections[2]->getSectionTangent())*(CommittedDeformation[2] - e))[2] / ((CommittedDeformation[2])[2] - e[2]);
			secantStiffness = (theSections[2]->getStressResultant())(ordering[2]) / e(ordering[2]);
			if (secantStiffness > OOPStiffness)   secantStiffness = OOPStiffness;
			if (secantStiffness < tangentStiffness)  secantStiffness = tangentStiffness;
			trialDuctilityDemand = OOPStiffness / secantStiffness;
			if (ductilityDemand < trialDuctilityDemand) {
				ductilityDemand = trialDuctilityDemand;
			}
		}

		if (abs(e(ordering[1])) > DBL_EPSILON) {
			tangentStiffness = ((theSections[2]->getSectionTangent())*(CommittedDeformation[2] - e))[1] / ((CommittedDeformation[2])[1] - e[2]);
			secantStiffness = (theSections[2]->getStressResultant())(ordering[1]) / e(ordering[1]);
			if (secantStiffness > OOPStiffness)   secantStiffness = OOPStiffness;
			if (secantStiffness < tangentStiffness)  secantStiffness = tangentStiffness;
			trialDuctilityDemand = OOPStiffness / secantStiffness;
			if (ductilityDemand < trialDuctilityDemand) {
				ductilityDemand = trialDuctilityDemand;
			}
		}
	}
	Vector N(4);
	N = theSections[1]->getStressResultant();


	// check for in-plane drift failure
	driftS = std::abs(-uBasic(10)) /2.0 *oneOverL;
// #ifdef _WINDOWS_
// 	driftF = std::max(std::abs(-(uBasic(1)*(2 * L) + uBasic(5)*L) / (2.0*L)),
//                 std::abs(-(uBasic(8)*(2 * L) + uBasic(5)*L) / (2.0*L)));
// #else
	driftF = std::max(std::abs(-(uBasic(1)*(2 * L) + uBasic(5)*L) / (2.0*L)),
		std::abs(-(uBasic(8)*(2 * L) + uBasic(5)*L) / (2.0*L)));
//#endif // _WINDOWS_


	double axialLoadRatio;
	if (Ltfc>0.0)
		axialLoadRatio = -N(ordering[0])/Ltfc;
	else
	    axialLoadRatio = 0.0;

	double shearSpanRatio;
	if (std::abs(M1 - M3) < DBL_EPSILON) {
		shearSpanRatio = 10.;
	}
	else {

// #ifdef _WINDOWS_
// 		shearSpanRatio = max(M1 / (M1 - M3),
// 			                 M3 / (M3 - M1));
// #else
		shearSpanRatio = std::max(M1 / (M1 - M3),
			                      M3 / (M3 - M1));
//#endif // _WINDOWS_

		if (shearSpanRatio > 10)
			shearSpanRatio = 10.;
	}



	Matrix EshearModel = theShearModel[0]->getInitialTangent();
	e(0) = N(ordering[0]) / EshearModel(0, 0);

	Vector s(2);
	s(0) = e(0);
	s(1) = uBasic(10);

	GenericDamagePlasticityShear* genericShearType = dynamic_cast<GenericDamagePlasticityShear*>(theShearModel[0]);
	if (genericShearType != NULL) {
		Vector params(2);
		params(0) = shearSpanRatio*2.*L;
		params(1) = axialLoadRatio;
		err += genericShearType->setTrialStrainExtraParams(s, params);
	}
	else {
		err += theShearModel[0]->setTrialStrain(s);
	}

	EshearModel = theShearModel[1]->getInitialTangent();
	e(0) = N(ordering[0]) / EshearModel(0, 0);

	s(0) = e(0);
	s(1) = uBasic(11);
	
	GenericDamagePlasticityShear* genericShearTypeOOP = dynamic_cast<GenericDamagePlasticityShear*>(theShearModel[1]);
	if (genericShearTypeOOP != NULL) {
		Vector params(2);
		params(0) = 2.*L;
		params(1) = axialLoadRatio;
		err += genericShearTypeOOP->setTrialStrainExtraParams(s, params);
	}
	else {
		err += theShearModel[1]->setTrialStrain(s);
	}

	
	this->driftModel(driftF, driftS, axialLoadRatio, shearSpanRatio);


    if (err != 0) {
        opserr << "Macroelement3d::update() - failed setTrialSectionDeformations()\n";
    return err;
    }

    return 0;
}

void 
Macroelement3d::driftModel(double /*currentDriftF*/, double /*currentDriftS*/, double axialLoadRatio, double H0overL) {

	if (axialLoadRatio<DBL_EPSILON)  axialLoadRatio = 0.0;
	if (axialLoadRatio>1.0)  axialLoadRatio = 0.999;

	// flexure model
	if (limitDriftF.Size()>0) {

		int index = 0;
		bool ALR_found = false;

		while (index<limitDriftF_ALR.Size() && !ALR_found) {
			if (axialLoadRatio >= limitDriftF_ALR(index) ) {
				index++;
			} else { 
				ALR_found = true;
			}
		}

		double driftF_capacity = limitDriftF(index-1) + (limitDriftF(index)-limitDriftF(index-1)) / (limitDriftF_ALR(index)-limitDriftF_ALR(index-1)) 
		                     * (axialLoadRatio - limitDriftF_ALR(index-1));
		driftF_capacity *= std::pow(H0overL, betaShearSpanF);
	    
		//check for loss of lateral force capacity
		if (std::abs(driftF)>driftF_capacity)    failedF = true;

		if (alphaAC_HC>0) {  // check for axial collapse
			if (std::abs(driftF)>driftF_capacity*alphaAC_HC)    collapsedF = true;
		}


	}
	 
	//shear model
	if (limitDriftS.Size()>0) {

		int index = 0;
		bool ALR_found = false;

		while (index<limitDriftS_ALR.Size() && !ALR_found) {
			if (axialLoadRatio >= limitDriftS_ALR(index) ) 
				index++;
			else 
				ALR_found = true;
		}

		double driftS_capacity = limitDriftS(index-1) + (limitDriftS(index)-limitDriftS(index-1)) / (limitDriftS_ALR(index)-limitDriftS_ALR(index-1)) 
		                     * (axialLoadRatio - limitDriftS_ALR(index-1));
		driftS_capacity *= std::pow(H0overL, betaShearSpanS);
	    
		//check for loss of lateral force capacity
		if (std::abs(driftS)>driftS_capacity)    failedS = true;

		if (alphaAC_HC>0) {  // check for axial collapse
			if (std::abs(driftS)>driftS_capacity*alphaAC_HC)    collapsedS = true;
		}
	}

	return;
}

Vector 
Macroelement3d::getBasicDisplacement(Vector dispI, Vector dispJ, Vector dispE) {

  // passes from global displacements/velocities to local
  Vector uGlobal(18);
  Vector uILocal(6);
  Vector uELocal(6);
  Vector uJLocal(6);

  for (int i=0; i<6; i++) {
       uGlobal(i)    = dispI(i);
       uGlobal(i+6)  = dispJ(i);
       uGlobal(i+12) = dispE(i);

       uILocal(i) = dispI(i);
       uJLocal(i) = dispJ(i);
       uELocal(i) = dispE(i);
  }

    
  if (nodeIInitialDisp != 0) {
        for (int j=0; j<6; j++) {
            uGlobal(j) -= nodeIInitialDisp[j];
            uILocal(j) -= nodeIInitialDisp[j];
		}
    }
    
    if (nodeJInitialDisp != 0) {
        for (int j=0; j<6; j++) {
            uGlobal(j+6) -= nodeJInitialDisp[j];
            uJLocal(j)   -= nodeIInitialDisp[j];
		}
    }

	// rotate to local
    Vector uLocal(18);
	uILocal = Tgl6*uILocal; 
    uJLocal = Tgl6*uJLocal; 
    uELocal = Tgl6*uELocal; 


	// correct if offset
	static double Wu[3];
    
    if (nodeIOffset) {
        Wu[0] =  nodeIOffset[2]*uGlobal(4) - nodeIOffset[1]*uGlobal(5);
        Wu[1] = -nodeIOffset[2]*uGlobal(3) + nodeIOffset[0]*uGlobal(5);
        Wu[2] =  nodeIOffset[1]*uGlobal(3) - nodeIOffset[0]*uGlobal(4);
        
		uILocal(0) += R[0][0]*Wu[0] + R[0][1]*Wu[1] + R[0][2]*Wu[2];
        uILocal(1) += R[1][0]*Wu[0] + R[1][1]*Wu[1] + R[1][2]*Wu[2];
        uILocal(2) += R[2][0]*Wu[0] + R[2][1]*Wu[1] + R[2][2]*Wu[2];
    }
    
    if (nodeJOffset) {
        Wu[0] =  nodeJOffset[2]*uGlobal(10) - nodeJOffset[1]*uGlobal(11);
        Wu[1] = -nodeJOffset[2]*uGlobal(9)  + nodeJOffset[0]*uGlobal(11);
        Wu[2] =  nodeJOffset[1]*uGlobal(9)  - nodeJOffset[0]*uGlobal(10);
        
		uJLocal(0) += R[0][0]*Wu[0] + R[0][1]*Wu[1] + R[0][2]*Wu[2];
        uJLocal(1) += R[1][0]*Wu[0] + R[1][1]*Wu[1] + R[1][2]*Wu[2];
        uJLocal(2) += R[2][0]*Wu[0] + R[2][1]*Wu[1] + R[2][2]*Wu[2];
    }

	
    for (int i=0; i<6; i++) {
       uLocal(i)    = uILocal(i);
       uLocal(i+6)  = uJLocal(i);
       uLocal(i+12) = uELocal(i);
    }

	double dV1, dV2, dV3;
	double dW1, dW2, dW3;

    // P-Delta transformation: update displacements
	if (PDelta) {   
		dV1 = uELocal(1) - uILocal(1);
		dW1 = uELocal(2) - uILocal(2);

		dV2 = uELocal(4) - uELocal(1);
		dW2 = uELocal(5) - uELocal(2);

		dV3 = uJLocal(1) - uELocal(4);
		dW3 = uJLocal(2) - uELocal(5);
	}

    // drop subscript Local for brevity
    Vector &uI = uILocal;
    Vector &uE = uELocal;
    Vector &uJ = uJLocal;

    double oneOverL = 1./L;

	Vector updatedBasic(12);

    // Get basic deformations. Total compatibility relations (linear) can be expressed also in matrix form as uBasic = Gamma_c * uLocal. 
	if (PDelta) {
		// non-linear total compatibility relations: OpenSees does not incluide this
		// first order approximationd of rotations, second order approximations of axial strain
		updatedBasic(0)  =  uE(0) - uI(0)  + 0.5*oneOverL*(std::pow(dV1,2) + std::pow(dW1,2));
		updatedBasic(1)  = -uI(5)  + oneOverL*dV1 ;
		updatedBasic(2)  = -uI(4)  - oneOverL*dW1;
		updatedBasic(3)  =  uJ(3) - uI(3);
		updatedBasic(4)  = -uE(0) + uE(3); // + 0.5*oneOverL*(-deltaV2*(uI(1)-uJ(1)) -deltaW2*(uI(2)-uJ(2)));  term to be included if exact formulation is applied
		updatedBasic(5)  =  oneOverL*(-dV1 + dV3);
		updatedBasic(6)  =  oneOverL*( dW1 - dW3);
		updatedBasic(7)  =  uJ(0) - uE(3) + 0.5*oneOverL*(std::pow(dV3,2) + std::pow(dW3,2));
		updatedBasic(8)  =  uJ(5) - oneOverL*dV3;
		updatedBasic(9)  =  uJ(4) + oneOverL*dW3;
		updatedBasic(10) =  dV2;  // the terms for the exact formulation are not included as they need an update of deltaU also
		updatedBasic(11) =  dW2;  // same here
	} else {
		// linear total compatibility relations
		updatedBasic(0)  =  uE(0) - uI(0);
		updatedBasic(1)  = -oneOverL*uI(1) - uI(5)  +oneOverL*uE(1);
		updatedBasic(2)  =  oneOverL*uI(2) - uI(4)  -oneOverL*uE(2);
		updatedBasic(3)  =  uJ(3) - uI(3);
		updatedBasic(4)  = -uE(0) + uE(3);
		updatedBasic(5)  =  oneOverL*(uI(1) + uJ(1) -uE(1) - uE(4));
		updatedBasic(6)  = oneOverL*(-uI(2) - uJ(2) +uE(2) + uE(5));
		updatedBasic(7)  = uJ(0) - uE(3);
		updatedBasic(8)  = -oneOverL*uJ(1) + uJ(5) + oneOverL*uE(4);
		updatedBasic(9)  =  oneOverL*uJ(2) + uJ(4) - oneOverL*uE(5);
		updatedBasic(10) = -uE(1) + uE(4);
		updatedBasic(11) = -uE(2) + uE(5);
	}


	if (PDelta) {
		deltaV1 = dV1;
		deltaW1 = dW1;

		deltaV3 = dV3;
		deltaW3 = dW3;

		getIncrementalCompatibilityMatrix(true);
	} 
	
	return updatedBasic;
}

const Matrix&
 Macroelement3d::getIncrementalCompatibilityMatrix(bool flagIncremental) {

	GammaC.Zero();
    double oneOverL = 1./L;

	    GammaC(0,0) = -1;
		GammaC(0,12) = 1;
		GammaC(1,1) = -oneOverL;
		GammaC(1,5) = -1;
		GammaC(1,13) = oneOverL;
		GammaC(2,2) = oneOverL;
		GammaC(2,4) = -1;
		GammaC(2,14) = -oneOverL;
		GammaC(3,3) = -1;
		GammaC(3,9) = 1;
		GammaC(4,12) = -1;
		GammaC(4,15) = 1;
		GammaC(5,1) = oneOverL;
		GammaC(5,7) = oneOverL;
		GammaC(5,13) = -oneOverL;
		GammaC(5,16) = -oneOverL;
		GammaC(6,2) = -oneOverL;
		GammaC(6,8) = -oneOverL;
		GammaC(6,14) = oneOverL;
		GammaC(6,17) = oneOverL;
		GammaC(7,6) = 1;
		GammaC(7,15) = -1;
		GammaC(8,7) = -oneOverL;
		GammaC(8,11) = 1;
		GammaC(8,16) = oneOverL;
		GammaC(9,8) = oneOverL;
		GammaC(9,10) = 1;
		GammaC(9,17) = -oneOverL;
		GammaC(10,13) = -1;
		GammaC(10,16) = 1;
		GammaC(11,14) = -1;
		GammaC(11,17) = 1;


    if ((PDelta) && flagIncremental) {

		GammaC(0,1) = -(deltaV1*oneOverL);
		GammaC(0,2) = -(deltaW1*oneOverL);
		GammaC(0,13) = deltaV1*oneOverL;
		GammaC(0,14) = deltaW1*oneOverL;
		GammaC(7,7) = deltaV3*oneOverL;
		GammaC(7,8) = deltaW3*oneOverL;
		GammaC(7,16) = -(deltaV3*oneOverL);
		GammaC(7,17) = -(deltaW3*oneOverL);
    }

    return GammaC;
  };

const Matrix&
Macroelement3d::getTangentStiff()
{
  static Matrix kb(12,12);
  
  kb.Zero();
  q.Zero();
  
  double oneOverL = 1.0/L;

  // first interface
  int order = theSections[0]->getOrder();
  const ID &code0 = theSections[0]->getType();

  int ordering[4];
  for (int j = 0; j < 4; j++)
	  ordering[j] = -1;

  for (int j = 0; j < order; j++) {
	  switch (code0(j)) {
	  case SECTION_RESPONSE_P:     ordering[0] = j;	break;
	  case SECTION_RESPONSE_MZ:    ordering[1] = j; break;
	  case SECTION_RESPONSE_MY:    ordering[2] = j;	break;
	  case SECTION_RESPONSE_T:     ordering[3] = j;	break;
	  default: break;
	  }
  }

  const Vector &s0  = theSections[0]->getStressResultant();
  const Matrix &ks0 = theSections[0]->getSectionTangent();

  // i index: line;  j index: column
  for (int i=0; i<4; i++) {
	  if (ordering[i] >= 0) {
		  q(i) = s0(ordering[i]);
		  for (int j = 0; j < 4; j++)
			  if (i < 3)
				  kb(i, j) = ks0(ordering[i], ordering[j]) / intLength(0);
			  else
				  kb(i, j) = ks0(ordering[i], ordering[j]) / (2.*L);
	  }
  }

  // second interface
  order = theSections[1]->getOrder();
  const ID &code1 = theSections[1]->getType();

  for (int j = 0; j < 4; j++)
	  ordering[j] = -1;

  for (int j = 0; j < order; j++) {
	  switch (code1(j)) {
	  case SECTION_RESPONSE_P:     ordering[0] = j;	break;
	  case SECTION_RESPONSE_MZ:    ordering[1] = j; break;
	  case SECTION_RESPONSE_MY:    ordering[2] = j;	break;
	  default: break;
	  }
  }

  const Vector &s1  = theSections[1]->getStressResultant();
  const Matrix &ks1 = theSections[1]->getSectionTangent();

  // i index: line;  j index: column
  for (int i = 0; i<3; i++) {
	  if (ordering[i] >= 0) {
		  q(4 + i) = s1(ordering[i]);
		  for (int j = 0; j < 3; j++)
			  kb(4 + i, 4 + j) = ks1(ordering[i], ordering[j]) / intLength(1);
	  }
  }


  // third interface
  order = theSections[2]->getOrder();
  const ID &code2 = theSections[2]->getType();

  for (int j = 0; j < 4; j++)
	  ordering[j] = -1;

  for (int j = 0; j < order; j++) {
	  switch (code2(j)) {
	  case SECTION_RESPONSE_P:     ordering[0] = j;	break;
	  case SECTION_RESPONSE_MZ:    ordering[1] = j; break;
	  case SECTION_RESPONSE_MY:    ordering[2] = j;	break;
	  default: break;
	  }
  }

  const Vector &s2  = theSections[2]->getStressResultant();
  const Matrix &ks2 = theSections[2]->getSectionTangent();

  // i index: line;  j index: column
  for (int i = 0; i<3; i++) {
	  if (ordering[i] >= 0) {
		  q(7 + i) = s2(ordering[i]);
		  for (int j = 0; j < 3; j++)
			  kb(7 + i, 7 + j) = ks2(ordering[i], ordering[j]) / intLength(2);
	  }
  }
 
  // shear interface 1 and 2
  for (int sect=0; sect<2; sect++) {

	Matrix EshearModel = theShearModel[sect]->getInitialTangent();
	//double N = q(4);
	double Esm = EshearModel(0,0);

    const Matrix &ks = theShearModel[sect]->getTangent();
    const Vector &s  = theShearModel[sect]->getStress();

    kb(10+sect,10+sect) = ks(1,1);
	q(10+sect) = s(1);

	// update non-diagonal terms of the stiffness matrix (dependency on N2)
	for (int j=0; j<3; j++)
		kb(10+sect, 4+j)  += ks(1,0) /Esm * kb(4, 4+j);
   
  }
 
  // kill lateral capacity if failed
  if (collapsedS || collapsedF) {
	  double failureFactor;
	  if (collapsedS) 
		  failureFactor = failureFactorS;
	  else
		  failureFactor = failureFactorF;

	  // dofs to kill: all
	  int dofsToKill[12] = {0,1,2,3,4,5,6,7,8,9,10,11}; 
	  int dof;

	  for (int kDof=0; kDof<12; kDof++) {
		  dof = dofsToKill[kDof];
		  q(dof) *= failureFactor;
		  for (int j=0; j<kb.noCols(); j++)
			  kb(dof,j) *= failureFactor;
	  }
  } else if (failedS || failedF) {

	  double failureFactor;
	  if (failedS) 
		  failureFactor = failureFactorS;
	  else
		  failureFactor = failureFactorF;

	  // dofs to kill: 1,2, 5,6, 8,9,  10,11  (moments of every section, shear of every section)
	  int dofsToKill[8] = {1,2, 5,6, 8,9,  10,11}; 
	  int dof;

	  for (int kDof=0; kDof<8; kDof++) {
		  dof = dofsToKill[kDof];
		  q(dof) *= failureFactor;
		  for (int j=0; j<kb.noCols(); j++)
			  kb(dof,j) *= failureFactor;
	  }
  }
  
   
  // Transform to local stiffness matrix
  
  Matrix Klocal(18,18); 
  //if (PDelta==1)   getIncrementalCompatibilityMatrix(true);
  Klocal.addMatrixTripleProduct(0.0, GammaC, kb, 1.0); 

  //opserr << "local stiffness = " << Klocal;
  
  // add geometric local stiffness matrix  
  if (PDelta==1) {
        Matrix kgeom(18,18);

	    // add fixed end forces. Update
	    for (int i=0; i<12; i++) {
		   q(i) += q0[i];
	    }

	    double qq0 = q(0);
	    double qq7 = q(7);
	    //double qq4 = q(4);

	  	kgeom(1,1) = qq0*oneOverL;
		kgeom(1,13) = -(qq0*oneOverL);
		kgeom(2,2) = qq0*oneOverL;
		kgeom(2,14) = -(qq0*oneOverL);
		kgeom(7,7) = qq7*oneOverL;
		kgeom(7,16) = -(qq7*oneOverL);
		kgeom(8,8) = qq7*oneOverL;
		kgeom(8,17) = -(qq7*oneOverL);
		kgeom(13,1) = -(qq0*oneOverL);
		kgeom(13,13) = qq0*oneOverL;
		kgeom(14,2) = -(qq0*oneOverL);
		kgeom(14,14) = qq0*oneOverL;
		kgeom(16,7) = -(qq7*oneOverL);
		kgeom(16,16) = qq7*oneOverL;
		kgeom(17,8) = -(qq7*oneOverL);
		kgeom(17,17) = qq7*oneOverL;

        Klocal += kgeom;
  }
   

   // Transform to global stiffness matrix   Kglobal = Tgl' * Klocal * Tgl
   // done now otherwise the secant stiffness matrix would be rrturned
	trasformMatrixToGlobal(Klocal);

  return K;
}

const Matrix&
Macroelement3d::getInitialBasicStiff()
{

  static Matrix kb(12,12);  
  kb.Zero();
  
  //double oneOverL = 1.0/L;

  // first interface
  int order = theSections[0]->getOrder();
  const ID &code0 = theSections[0]->getType();

  int ordering[4];
  for (int j = 0; j < 4; j++)
	  ordering[j] = -1;

  for (int j = 0; j < order; j++) {
	  switch (code0(j)) {
	  case SECTION_RESPONSE_P:     ordering[0] = j;	break;
	  case SECTION_RESPONSE_MZ:    ordering[1] = j; break;
	  case SECTION_RESPONSE_MY:    ordering[2] = j;	break;
	  case SECTION_RESPONSE_T:     ordering[3] = j;	break;
	  default: break;
	  }
  }
  
  const Matrix &ks0 = theSections[0]->getInitialTangent();
  
  // i index: line;  j index: column
  for (int i = 0; i<4; i++) {
	  if (ordering[i] >= 0) {
		  for (int j = 0; j < 4; j++)
			  if (i < 3)
				  kb(i, j) = ks0(ordering[i], ordering[j]) / intLength(0);
			  else
				  kb(i, j) = ks0(ordering[i], ordering[j]) / (2.*L);
	  }
  }

  // second interface
  order = theSections[1]->getOrder();
  const ID &code1 = theSections[1]->getType();

  for (int j = 0; j < 4; j++)
	  ordering[j] = -1;

  for (int j = 0; j < order; j++) {
	  switch (code1(j)) {
	  case SECTION_RESPONSE_P:     ordering[0] = j;	break;
	  case SECTION_RESPONSE_MZ:    ordering[1] = j; break;
	  case SECTION_RESPONSE_MY:    ordering[2] = j;	break;
	  default: break;
	  }
  }

  const Matrix &ks1 = theSections[1]->getInitialTangent();

  // i index: line;  j index: column
  for (int i = 0; i<3; i++) {
	  if (ordering[i] >= 0) {
		  for (int j = 0; j < 3; j++)
			  kb(4 + i, 4 + j) = ks1(ordering[i], ordering[j]) / intLength(1);
	  }
  }


  // third interface
  order = theSections[2]->getOrder();
  const ID &code2 = theSections[2]->getType();

  for (int j = 0; j < 4; j++)
	  ordering[j] = -1;

  for (int j = 0; j < order; j++) {
	  switch (code2(j)) {
	  case SECTION_RESPONSE_P:     ordering[0] = j;	break;
	  case SECTION_RESPONSE_MZ:    ordering[1] = j; break;
	  case SECTION_RESPONSE_MY:    ordering[2] = j;	break;
	  default: break;
	  }
  }

  const Matrix &ks2 = theSections[2]->getInitialTangent();

  // i index: line;  j index: column
  for (int i = 0; i<3; i++) {
	  if (ordering[i] >= 0) {
		  for (int j = 0; j < 3; j++)
			  kb(7 + i, 7 + j) = ks2(ordering[i], ordering[j]) / intLength(2);
	  }
  }

  // shear interface 1 and 2
  for (int sect=0; sect<2; sect++) {
    const Matrix &ks = theShearModel[sect]->getInitialTangent();
    kb(10+sect,10+sect) = ks(1,1);
	// no need to update diagonal terms->they are already uncoupled
  }


  // kill lateral capacity if failed - added even if the initial stiffness is returned to help convergence qhen one element fails
  if (collapsedS || collapsedF) {
	  double failureFactor;
	  if (collapsedS)
		  failureFactor = failureFactorS;
	  else
		  failureFactor = failureFactorF;

	  // dofs to kill: all
	  int dofsToKill[12] = { 0,1,2,3,4,5,6,7,8,9,10,11 };
	  int dof;

	  for (int kDof = 0; kDof<12; kDof++) {
		  dof = dofsToKill[kDof];
		  q(dof) *= failureFactor;
		  for (int j = 0; j<kb.noCols(); j++)
			  kb(dof, j) *= failureFactor;
	  }
  }
  else if (failedS || failedF) {

	  double failureFactor;
	  if (failedS)
		  failureFactor = failureFactorS;
	  else
		  failureFactor = failureFactorF;

	  // dofs to kill: 1,2, 5,6, 8,9,  10,11  (moments of every section, shear of every section)
	  int dofsToKill[8] = { 1,2, 5,6, 8,9,  10,11 };
	  int dof;

	  for (int kDof = 0; kDof<8; kDof++) {
		  dof = dofsToKill[kDof];
		  q(dof) *= failureFactor;
		  for (int j = 0; j<kb.noCols(); j++)
			  kb(dof, j) *= failureFactor;
	  }
  }

  return kb;
}

const Matrix&
Macroelement3d::getInitialStiff()
{
  const Matrix &kb = this->getInitialBasicStiff();

  Matrix Klocal(18,18);
  Klocal.addMatrixTripleProduct(0.0, GammaC, kb, 1.0); 
  trasformMatrixToGlobal(Klocal);
  
  return K;
}

const Matrix&
Macroelement3d::getMass()
{
  K.Zero();
  if (rho == 0.0)
    return K;
  
  if (cMass == 0)  {
    // lumped mass matrix
	  if (isGable) {
		double m = 0.5*rho*2.*L;
		K(0,0) = 5./12.*m *massglobalDir(0);
		K(1,1) = 5./12.*m *massglobalDir(1);
		K(2,2) = 5./12.*m *massglobalDir(2);
		
		K(6,6) = 1./12.*m *massglobalDir(0);
		K(7,7) = 1./12.*m *massglobalDir(1);
		K(8,8) = 1./12.*m *massglobalDir(2);

		K(12,12) = 1./3.*m *massglobalDir(0);
		K(13,13) = 1./3.*m *massglobalDir(1);
		K(14,14) = 1./3.*m *massglobalDir(2);

		K(15,15) = 1./6.*m *massglobalDir(0);
		K(16,16) = 1./6.*m *massglobalDir(1);
		K(17,17) = 1./6.*m *massglobalDir(2);
		
	  } else {
		double m = rho;
		K(0,0) = m *intLengthMasses(0) *massglobalDir(0);
		K(1,1) = m *intLengthMasses(0) *massglobalDir(1);
		K(2,2) = m *intLengthMasses(0) *massglobalDir(2);
		
		K(6,6) = m *intLengthMasses(2) *massglobalDir(0);
		K(7,7) = m *intLengthMasses(2) *massglobalDir(1);
		K(8,8) = m *intLengthMasses(2) *massglobalDir(2);

		K(12,12) = m *intLengthMasses(1)/2. *massglobalDir(0);
		K(13,13) = m *intLengthMasses(1)/2. *massglobalDir(1);
		K(14,14) = m *intLengthMasses(1)/2. *massglobalDir(2);
		
		K(15,15) = m *intLengthMasses(1)/2. *massglobalDir(0);
		K(16,16) = m *intLengthMasses(1)/2. *massglobalDir(1);
		K(17,17) = m *intLengthMasses(1)/2. *massglobalDir(2);
	  }

  } else  {
    // consistent mass matrix	  

	 Matrix massI_over12(3,3);
	 double m = rho*L / 6.;  // 1/12 mTot, 1/6 mTot if gable

	 if (isGable) {  
		 m /= 4.;     // mTotGable/24

		 massI_over12(0, 0) = m *massglobalDir(0);
		 massI_over12(1, 1) = m *massglobalDir(1);
		 massI_over12(2, 2) = m *massglobalDir(2);

		 // node I
		 K(0, 0) = 7.*massI_over12(0, 0);
		 K(1, 1) = 7.*massI_over12(1, 1);
		 K(2, 2) = 7.*massI_over12(2, 2);

		 // node J
		 K(6, 6) = massI_over12(0, 0);
		 K(7, 7) = massI_over12(1, 1);
		 K(8, 8) = massI_over12(2, 2);

		 // node E
		 K(12, 12) = 5.*massI_over12(0, 0);
		 K(13, 13) = 5.*massI_over12(1, 1);
		 K(14, 14) = 5.*massI_over12(2, 2);

		 K(15, 15) = 3.*massI_over12(0, 0);
		 K(16, 16) = 3.*massI_over12(1, 1);
		 K(17, 17) = 3.*massI_over12(2, 2);

		 // coupled dofs node I
		 K(0, 12) = 3.*massI_over12(0, 0);
		 K(1, 13) = 3.*massI_over12(1, 1);
		 K(2, 14) = 3.*massI_over12(2, 2);

		 K(12, 0) = 3.*massI_over12(0, 0);
		 K(13, 1) = 3.*massI_over12(1, 1);
		 K(14, 2) = 3.*massI_over12(2, 2);

		 // coupled dofs node J
		 K(6, 15) = massI_over12(0, 0);
		 K(7, 16) = massI_over12(1, 1);
		 K(8, 17) = massI_over12(2, 2);

		 K(15, 6) = massI_over12(0, 0);
		 K(16, 7) = massI_over12(1, 1);
		 K(17, 8) = massI_over12(2, 2);

	 }
	 else {
		 massI_over12(0, 0) = m *massglobalDir(0);
		 massI_over12(1, 1) = m *massglobalDir(1);
		 massI_over12(2, 2) = m *massglobalDir(2);

		 // node I
		 K(0, 0) = 2.*massI_over12(0, 0);
		 K(1, 1) = 2.*massI_over12(1, 1);
		 K(2, 2) = 2.*massI_over12(2, 2);

		 // node J
		 K(6, 6) = 2.*massI_over12(0, 0);
		 K(7, 7) = 2.*massI_over12(1, 1);
		 K(8, 8) = 2.*massI_over12(2, 2);

		 // node E
		 K(12, 12) = 2.*massI_over12(0, 0);
		 K(13, 13) = 2.*massI_over12(1, 1);
		 K(14, 14) = 2.*massI_over12(2, 2);

		 K(15, 15) = 2.*massI_over12(0, 0);
		 K(16, 16) = 2.*massI_over12(1, 1);
		 K(17, 17) = 2.*massI_over12(2, 2);

		 // coupled dofs node I
		 K(0, 12) = massI_over12(0, 0);
		 K(1, 13) = massI_over12(1, 1);
		 K(2, 14) = massI_over12(2, 2);

		 K(12, 0) = massI_over12(0, 0);
		 K(13, 1) = massI_over12(1, 1);
		 K(14, 2) = massI_over12(2, 2);

		 // coupled dofs node J
		 K(6, 15) = massI_over12(0, 0);
		 K(7, 16) = massI_over12(1, 1);
		 K(8, 17) = massI_over12(2, 2);

		 K(15, 6) = massI_over12(0, 0);
		 K(16, 7) = massI_over12(1, 1);
		 K(17, 8) = massI_over12(2, 2);
	 }
	 

	 // apply offset
	 if (nodeIOffset) {

		 K(0, 4) = K(4, 0) =  nodeIOffset[2] * K(0, 0);
		 K(0, 5) = K(5, 0) = -nodeIOffset[1] * K(0, 0);
		 K(1, 3) = K(3, 1) = -nodeIOffset[2] * K(1, 1);
		 K(1, 5) = K(5, 1) =  nodeIOffset[0] * K(1, 1);
		 K(2, 3) = K(3, 2) =  nodeIOffset[1] * K(2, 2);
		 K(2, 4) = K(4, 2) = -nodeIOffset[0] * K(2, 2);

		 K(3, 3) = (std::pow(nodeIOffset[1], 2) + std::pow(nodeIOffset[2], 2))* K(0, 0);
		 K(4, 4) = (std::pow(nodeIOffset[0], 2) + std::pow(nodeIOffset[2], 2))* K(1, 1);
		 K(5, 5) = (std::pow(nodeIOffset[0], 2) + std::pow(nodeIOffset[1], 2))* K(2, 2);

		 K(3, 4) = K(4, 3) = -nodeIOffset[0] * nodeIOffset[1] * K(2, 2);
		 K(3, 5) = K(5, 3) = -nodeIOffset[0] * nodeIOffset[2] * K(1, 1);
		 K(4, 5) = K(5, 4) = -nodeIOffset[1] * nodeIOffset[2] * K(0, 0);


		 K(12, 4) = K(4, 12) =  nodeIOffset[2] * K(12, 0);
		 K(12, 5) = K(5, 12) = -nodeIOffset[1] * K(12, 0); 
		 K(13, 3) = K(3, 13) = -nodeIOffset[2] * K(13, 1);
		 K(13, 5) = K(5, 13) =  nodeIOffset[0] * K(13, 1);
		 K(14, 3) = K(3, 14) =  nodeIOffset[1] * K(14, 2);
		 K(14, 4) = K(4, 14) = -nodeIOffset[0] * K(14, 2);
	 }


	 if (nodeJOffset) {

		 K(6, 10) = K(10, 6) =  nodeJOffset[2] * K(6, 6);
		 K(6, 11) = K(11, 6) = -nodeJOffset[1] * K(6, 6);
		 K(7, 9)  = K(9,  6) = -nodeJOffset[2] * K(7, 7);
		 K(7, 11) = K(11, 7) =  nodeJOffset[0] * K(7, 7);
		 K(8, 9)  = K(9,  8) =  nodeJOffset[1] * K(8, 8);
		 K(8, 10) = K(10, 8) = -nodeJOffset[0] * K(8, 8);

		 K(9,  9)  = (std::pow(nodeJOffset[1], 2) + std::pow(nodeJOffset[2], 2))* K(6, 6);
		 K(10, 10) = (std::pow(nodeJOffset[0], 2) + std::pow(nodeJOffset[2], 2))* K(7, 7);
		 K(11, 11) = (std::pow(nodeJOffset[0], 2) + std::pow(nodeJOffset[1], 2))* K(8, 8);

		 K(9,  10) = K(6,  9)  = -nodeJOffset[0] * nodeJOffset[1] * K(8, 8);
		 K(9,  11) = K(11, 9)  = -nodeJOffset[0] * nodeJOffset[2] * K(7, 7);
		 K(10, 11) = K(11, 10) = -nodeJOffset[1] * nodeJOffset[2] * K(6, 6);

		 K(15, 10) = K(10, 15) =  nodeIOffset[2] * K(15, 6);
		 K(15, 11) = K(11, 15) = -nodeIOffset[1] * K(15, 6);
		 K(16, 9)  = K(9,  16) = -nodeIOffset[2] * K(16, 7);
		 K(16, 11) = K(11, 16) =  nodeIOffset[0] * K(16, 7);
		 K(17, 9)  = K(9,  17) =  nodeIOffset[1] * K(17, 8);
		 K(17, 10) = K(10, 17) = -nodeIOffset[0] * K(17, 8);

	 } 	
  }
  
  return K;
}

void
Macroelement3d::zeroLoad(void)
{
 Q.Zero();

  for (int i=0; i<12; i++) {
	 q0[i] = 0.0;
	 p0[i] = 0.0;
  }
  return;
}

int 
Macroelement3d::addLoad(ElementalLoad *theLoad, double loadFactor)
{
   int type;
   const Vector &data = theLoad->getData(type, loadFactor);

   if (type == LOAD_TAG_SelfWeight) {

	   Vector gravity(6);
	   gravity(0) = data(0) * loadFactor;   // gravity direction in global coordinates
	   gravity(1) = data(1) * loadFactor;
	   gravity(2) = data(2) * loadFactor;

	   // rotate to local
	   gravity = Tgl6*gravity; 

	    double wx = gravity(0)*rho;  // Axial (+direction from node I to J)
		double wy = gravity(1)*rho;  // Transverse y
		double wz = gravity(2)*rho;  // Transverse z

		if (isGable) {
			double Vy = wy*L;
			double Vz = wz*L;
			double Fx = wx*L;   

			// Reactions in basic system
			p0[0] -= 5./12.*Fx;
			p0[1] -= 1./12.*Fx;
			p0[2] -= 1./3. *Fx;
			p0[3] -= 1./6. *Fx;

		    p0[4] -= 5./12.*Vy;
			p0[5] -= 1./12.*Vy;
			p0[6] -= 1./3. *Vy;
			p0[7] -= 1./6. *Vy;

			p0[8]  -= 5./12.*Vz;
			p0[9]  -= 1./12.*Vz;
			p0[10] -= 1./3. *Vz;
			p0[11] -= 1./6. *Vz;
		
			// Fixed end forces in basic system
			q0[0] -= 5./12.*Fx;
			q0[7] += 1./12.*Fx;

		} else {

			double Vy = 0.5*wy*L;
			double Vz = 0.5*wz*L;
			double Fx = wx*L*2.;   

			// Reactions in basic system
			p0[0] -= 0.25*Fx;
			p0[1] -= 0.25*Fx;
			p0[2] -= 0.25*Fx;
			p0[3] -= 0.25*Fx;

			p0[4] -= Vy;
			p0[5] -= Vy;
			p0[6] -= Vy;
			p0[7] -= Vy;

			p0[8] -= Vz;
			p0[9] -= Vz;
			p0[10] -= Vz;
			p0[11] -= Vz;

			// Fixed end forces in basic system
			q0[0] -= 0.25*Fx;
			q0[7] += 0.25*Fx;
		}
  }
   else if (type == LOAD_TAG_Beam3dUniformLoad) {
	    double wy = data(0)*loadFactor;  // Transverse y
		double wz = data(1)*loadFactor;  // Transverse z
		double wx = data(2)*loadFactor;  // Axial 

			double Vy = 0.5*wy*L;
			double Vz = 0.5*wz*L;
			double Fx = wx*L*2.;   

			// Reactions in basic system
			p0[0] -= 0.25*Fx;
			p0[1] -= 0.25*Fx;
			p0[2] -= 0.25*Fx;
			p0[3] -= 0.25*Fx;
		
			p0[4] -= Vy;
			p0[5] -= Vy;
			p0[6] -= Vy;
			p0[7] -= Vy;
		
			p0[8] -= Vz;
			p0[9] -= Vz;
			p0[10] -= Vz;
			p0[11] -= Vz;

			// Fixed end forces in basic system
			q0[0] -= 0.25*Fx;
			q0[7] += 0.25*Fx;
   }
   else {
		opserr << "Macroelement3d::addLoad() -- load type unknown for element with tag: " << this->getTag() << endln;
		opserr << "Available elemental load types: -selfWeight $gX $gY $gZ (global system) and -beamUniform $Wy $Wz $Wx (local system)\n";
		return -1;
  }

   return 0;
 }

int 
Macroelement3d::addInertiaLoadToUnbalance(const Vector &accel)
{
  if (rho == 0.0) {
    return 0;
  }
  
  // Get R * accel from the nodes
  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);
  const Vector &RaccelE = theNodes[2]->getRV(accel);
  
  if (6 != Raccel1.Size() || 6 != Raccel2.Size()) {
    opserr << "Macroelement3d::addInertiaLoadToUnbalance matrix and vector sizes are incompatible\n";
    return -1;
  }

  
  // want to add ( - fact * M R * accel ) to unbalance
  if (cMass == 0)  {
    // take advantage of lumped mass matrix
	if (isGable) {
		double m = 0.5*rho*2.*L;
		Q(0)  -= 5./12.*m *Raccel1(0) *massglobalDir(0);
		Q(1)  -= 5./12.*m *Raccel1(1) *massglobalDir(1);
		Q(2)  -= 5./12.*m *Raccel1(2) *massglobalDir(2);
		
		Q(6)  -= 1./12.*m *Raccel2(0) *massglobalDir(0);
		Q(7)  -= 1./12.*m *Raccel2(1) *massglobalDir(1);
		Q(8)  -= 1./12.*m *Raccel2(2) *massglobalDir(2);

		Q(12) -= 1./3.*m *RaccelE(0) *massglobalDir(0);
		Q(13) -= 1./3.*m *RaccelE(1) *massglobalDir(1);
		Q(14) -= 1./3.*m *RaccelE(2) *massglobalDir(2);

		Q(15) -= 1./6.*m *RaccelE(3) *massglobalDir(0);
		Q(16) -= 1./6.*m *RaccelE(4) *massglobalDir(1);
		Q(17) -= 1./6.*m *RaccelE(5) *massglobalDir(2);
		
	  } else {
		double m = rho;
		Q(0)  -= m *intLengthMasses(0) *Raccel1(0) *massglobalDir(0);
		Q(1)  -= m *intLengthMasses(0) *Raccel1(1) *massglobalDir(1);
		Q(2)  -= m *intLengthMasses(0) *Raccel1(2) *massglobalDir(2);
		
		Q(6)  -= m *intLengthMasses(2) *Raccel2(0) *massglobalDir(0);
		Q(7)  -= m *intLengthMasses(2) *Raccel2(1) *massglobalDir(1);
		Q(8)  -= m *intLengthMasses(2) *Raccel2(2) *massglobalDir(2);

		Q(12) -= m *intLengthMasses(1)/2. *RaccelE(0) *massglobalDir(0);
		Q(13) -= m *intLengthMasses(1)/2. *RaccelE(1) *massglobalDir(1);
		Q(14) -= m *intLengthMasses(1)/2. *RaccelE(2) *massglobalDir(2);

		Q(15) -= m *intLengthMasses(1)/2. *RaccelE(3) *massglobalDir(0);
		Q(16) -= m *intLengthMasses(1)/2. *RaccelE(4) *massglobalDir(1);
		Q(17) -= m *intLengthMasses(1)/2. *RaccelE(5) *massglobalDir(2);
	  }


  } else  {
    // use matrix vector multip. for consistent mass matrix
    static Vector Raccel(18);
    for (int i=0; i<6; i++)  {
      Raccel(i)    = Raccel1(i);
      Raccel(i+6)  = Raccel2(i);
	  Raccel(i+12) = RaccelE(i);
    }

    Q.addMatrixVector(1.0, this->getMass(), Raccel, -1.0);
  }
  
  return 0;
}

const Vector&
Macroelement3d::getResistingForce()
{ 
  P.Zero();
  q.Zero();
  
  //double oneOverL = 1.0/L;

  // first interface
  int order, j;
  int ordering[4];
  const Vector &s0  = theSections[0]->getStressResultant();

  order = theSections[0]->getOrder();
  const ID &code0 = theSections[0]->getType();

  for (j = 0; j < 4; j++)
	  ordering[j] = -1;

  for (j = 0; j < order; j++) {
	  switch (code0(j)) {
	  case SECTION_RESPONSE_P:     ordering[0] = j;	break;
	  case SECTION_RESPONSE_MZ:    ordering[1] = j; break;
	  case SECTION_RESPONSE_MY:    ordering[2] = j;	break;
	  case SECTION_RESPONSE_T:     ordering[3] = j;	break;
	  default: break;
	  }
  }

  for (int i = 0; i<4; i++) {
	  if (ordering[i] >= 0) {
     	  q(i) = s0(ordering[i]);
	  }
  }

  // second interface
  order = theSections[1]->getOrder();
  const ID &code1 = theSections[1]->getType();

  for (j = 0; j < 4; j++)
	  ordering[j] = -1;

  for (j = 0; j < order; j++) {
	  switch (code1(j)) {
	  case SECTION_RESPONSE_P:     ordering[0] = j;	break;
	  case SECTION_RESPONSE_MZ:    ordering[1] = j; break;
	  case SECTION_RESPONSE_MY:    ordering[2] = j;	break;
	  default: break;
	  }
  }

  const Vector &s1 = theSections[1]->getStressResultant();
  for (int i = 0; i<3; i++) {
	  if (ordering[i] >= 0) {
		  q(4 + i) = s1(ordering[i]);
	  }
  }


  // third interface
  order = theSections[2]->getOrder();
  const ID &code2 = theSections[2]->getType();

  for (j = 0; j < 4; j++)
	  ordering[j] = -1;

  for (j = 0; j < order; j++) {
	  switch (code2(j)) {
	  case SECTION_RESPONSE_P:     ordering[0] = j;	break;
	  case SECTION_RESPONSE_MZ:    ordering[1] = j; break;
	  case SECTION_RESPONSE_MY:    ordering[2] = j;	break;
	  default: break;
	  }
  }

  const Vector &s2 = theSections[2]->getStressResultant();
  for (int i = 0; i<3; i++) {
	  if (ordering[i] >= 0) {
		  q(7 + i) = s2(ordering[i]);
	  }
  }

  // shear interface 1 and 2
  for (int sect=0; sect<2; sect++) {
    const Vector &s  = theShearModel[sect]->getStress();
    q(10+sect) = s(1);
  }

    if (collapsedS || collapsedF) {

		double failureFactor;
	  if (collapsedS) 
		  failureFactor = failureFactorS;
	  else
		  failureFactor = failureFactorF;

	  // dofs to kill: all
	  int dofsToKill[12] = {0,1,2,3,4,5,6,7,8,9,10,11}; 
	  int dof;

	  for (int kDof=0; kDof<12; kDof++) {
		  dof = dofsToKill[kDof];
		  q(dof) *= failureFactor;
	  }
   } else if (failedS || failedF) {

	   double failureFactor;
	  if (failedS) 
		  failureFactor = failureFactorS;
	  else
		  failureFactor = failureFactorF;

	  // dofs to kill: 1,2, 5,6, 8,9,  10,11  (moments of every section, shear of every section)
	  int dofsToKill[8] = {1,2, 5,6, 8,9,  10,11}; 
	  int dof;

	  for (int kDof=0; kDof<8; kDof++) {
		  dof = dofsToKill[kDof];
		  q(dof) *= failureFactor;
	  }
   }

  // add fixed end forces . Update
  for (int i=0; i<12; i++) {
	  q(i) += q0[i];
  }
     
  // Transform forces to local reference system
  Vector pLocal(18);

  pLocal = GammaC ^ q;  // pLocal = GammaE * q = GammaC'*q;

  // add fixed reactions
  pLocal(0)  += p0[0];
  pLocal(6)  += p0[1];
  pLocal(12) += p0[2];
  pLocal(15) += p0[3];

  pLocal(1)  += p0[4];
  pLocal(7)  += p0[5];
  pLocal(13) += p0[6];
  pLocal(16) += p0[7];

  pLocal(2)  += p0[8];
  pLocal(8)  += p0[9];
  pLocal(14) += p0[10];
  pLocal(17) += p0[11];
  

  // Transform forces to global reference system
  P = Tgl ^ pLocal;

  // add offset effect
    if (nodeIOffset) {
        P(3) += -nodeIOffset[2]*P(1) + nodeIOffset[1]*P(2);
        P(4) +=  nodeIOffset[2]*P(0) - nodeIOffset[0]*P(2);
        P(5) += -nodeIOffset[1]*P(0) + nodeIOffset[0]*P(1);
    }
    
    if (nodeJOffset) {
        P(9)  += -nodeJOffset[2]*P(7) + nodeJOffset[1]*P(8);
        P(10) +=  nodeJOffset[2]*P(6) - nodeJOffset[0]*P(8);
        P(11) += -nodeJOffset[1]*P(6) + nodeJOffset[0]*P(7);
    }

  return P;
}

const Vector&
Macroelement3d::getResistingForceIncInertia()
{
  P = this->getResistingForce();

  // Subtract other external nodal loads ... P_res = P_int - P_ext
  P.addVector(1.0, Q, -1.0);

  if (rho != 0.0) {
    const Vector &accel1 = theNodes[0]->getTrialAccel();
    const Vector &accel2 = theNodes[1]->getTrialAccel();
	const Vector &accelE = theNodes[2]->getTrialAccel();
    
  if (cMass == 0)  {
    // take advantage of lumped mass matrix
	if (isGable) {
		double m = 0.5*rho*2.*L;
		P(0)  += 5./12.*m* accel1(0) *massglobalDir(0);
		P(1)  += 5./12.*m *accel1(1) *massglobalDir(1);
		P(2)  += 5./12.*m *accel1(2) *massglobalDir(2);
		
		P(6)  += 1./12.*m *accel2(0) *massglobalDir(0);
		P(7)  += 1./12.*m *accel2(1) *massglobalDir(1);
		P(8)  += 1./12.*m *accel2(2) *massglobalDir(2);

		P(12) += 1./3.*m *accelE(0) *massglobalDir(0);
		P(13) += 1./3.*m *accelE(1) *massglobalDir(1);
		P(14) += 1./3.*m *accelE(2) *massglobalDir(2);

		P(15) += 1./6.*m *accelE(3) *massglobalDir(0);
		P(16) += 1./6.*m *accelE(4) *massglobalDir(1);
		P(17) += 1./6.*m *accelE(5) *massglobalDir(2);
		
	} else {
		double m = rho;
		P(0)  += m *intLengthMasses(0) *accel1(0) *massglobalDir(0);
		P(1)  += m *intLengthMasses(0) *accel1(1) *massglobalDir(1);
		P(2)  += m *intLengthMasses(0) *accel1(2) *massglobalDir(2);
		
		P(6)  += m *intLengthMasses(2) *accel2(0) *massglobalDir(0);
		P(7)  += m *intLengthMasses(2) *accel2(1) *massglobalDir(1);
		P(8)  += m *intLengthMasses(2) *accel2(2) *massglobalDir(2);

		P(12) += m *intLengthMasses(1)/2. *accelE(0) *massglobalDir(0);
		P(13) += m *intLengthMasses(1)/2. *accelE(1) *massglobalDir(1);
		P(14) += m *intLengthMasses(1)/2. *accelE(2) *massglobalDir(2);

		P(15) += m *intLengthMasses(1)/2. *accelE(3) *massglobalDir(0);
		P(16) += m *intLengthMasses(1)/2. *accelE(4) *massglobalDir(1);
		P(17) += m *intLengthMasses(1)/2. *accelE(5) *massglobalDir(2);
	  }


  } else  {
    // use matrix vector multip. for consistent mass matrix
    static Vector accel(18);
    for (int i=0; i<6; i++)  {
      accel(i)   = accel1(i);
      accel(i+6) = accel2(i);
	  accel(i+12) = accelE(i);
    }
    P.addMatrixVector(1.0, this->getMass(), accel, 1.0);
  }

   
   // add the damping forces if rayleigh damping
  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0) {
	     Vector dampingForces = this->getRayleighDampingForces();
         P.addVector(1.0, dampingForces, 1.0);
  }

  } else {
    // add the damping forces if rayleigh damping
    if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0 ) {
      P.addVector(1.0, this->getRayleighDampingForces(), 1.0);
	  
	}
  }

  return P;
}

int
Macroelement3d::sendSelf(int /*commitTag*/, Channel &/*theChannel*/)
{
	opserr << "Macroelement3d::sendSelf: method not implemented\n";
	return -1;
}

int
Macroelement3d::recvSelf(int /*commitTag*/, Channel &/*theChannel*/,
						 FEM_ObjectBroker &/*theBroker*/)
{
	opserr << "Macroelement3d::recvSelf: method not implemented\n";
	return -1;
}

void
Macroelement3d::Print(OPS_Stream &s, int /*flag*/)
{
  s << "\nMacroelement3d, element id:  " << this->getTag() << endln;
  s << "\tConnected external nodes:  " << connectedExternalNodes;
}

Response*
Macroelement3d::setResponse(const char **argv, int argc, OPS_Stream &output)
{
    Response *theResponse = 0;

    output.tag("ElementOutput");
    output.attr("eleType","Macroelement3d");
    output.attr("eleTag",this->getTag());
    output.attr("node1",connectedExternalNodes[0]);
    output.attr("node2",connectedExternalNodes[1]);


    // global force - 
    if (strcmp(argv[0],"forces") == 0 || strcmp(argv[0],"force") == 0
	|| strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0) {

      output.tag("ResponseType","Px_1");
      output.tag("ResponseType","Py_1");
      output.tag("ResponseType","Pz_1");
      output.tag("ResponseType","Mx_1");
      output.tag("ResponseType","My_1");
      output.tag("ResponseType","Mz_1");
      output.tag("ResponseType","Fe_11");
      output.tag("ResponseType","Fe_21");
      output.tag("ResponseType","Fe_31");
      output.tag("ResponseType","Fe_12");
      output.tag("ResponseType","Fe_22");
      output.tag("ResponseType","Fe_32");
      output.tag("ResponseType","Px_2");
      output.tag("ResponseType","Py_2");
      output.tag("ResponseType","Pz_2");
      output.tag("ResponseType","Mx_2");
      output.tag("ResponseType","My_2");
      output.tag("ResponseType","Mz_2");

      theResponse = new ElementResponse(this, 1, P);

    // basic force - check if correct (allows confusion between basic forces and local forces)
    }  else if (strcmp(argv[0],"basicForces") == 0 || strcmp(argv[0],"localForces") == 0 || strcmp(argv[0], "localForce") == 0) {

      output.tag("ResponseType","N_1");
      output.tag("ResponseType","Mz_1");
      output.tag("ResponseType","My_1");
      output.tag("ResponseType","T_1");
      output.tag("ResponseType","N_2");
      output.tag("ResponseType","Mz_2");
      output.tag("ResponseType","My_2");
	  output.tag("ResponseType","N_3");
      output.tag("ResponseType","Mz_3");
      output.tag("ResponseType","My_3");
	  output.tag("ResponseType","V_y");
      output.tag("ResponseType","V_z");

      theResponse = new ElementResponse(this, 2, q);

    // in-plane drifts - shear and flexure
    }  else if (strcmp(argv[0],"drifts") == 0 || strcmp(argv[0],"drift") == 0 
	      || strcmp(argv[0],"driftSF") == 0) {

      output.tag("ResponseType","driftShear");
      output.tag("ResponseType","driftFlexure");

      theResponse = new ElementResponse(this, 3, Vector(2));


	// generic sectional response
	}  else if (strstr(argv[0],"section") != 0) {

    if (argc > 1) {
      
      int sectionNum = atoi(argv[1]);
      
      if (sectionNum > 0 && sectionNum <= 3 && argc > 2) {

		output.tag("SectionOutput");
		output.attr("number",sectionNum);

		theResponse = theSections[sectionNum-1]->setResponse(&argv[2], argc-2, output);
		output.endTag();
      
      } else if (sectionNum == 0) { // argv[1] was not an int, we want all sections, 
	
		CompositeResponse *theCResponse = new CompositeResponse();
		int numResponse = 0;

		for (int i=0; i<numSections; i++) {
	  
			output.tag("SectionOutput");
			output.attr("number",i+1);
	  
			Response *theSectionResponse = theSections[i]->setResponse(&argv[1], argc-1, output);

			output.endTag();	  

			if (theSectionResponse != 0) {
				numResponse = theCResponse->addResponse(theSectionResponse);
			}
		}
	
		if (numResponse == 0) // no valid responses found
			delete theCResponse;
		else
			theResponse = theCResponse;
	  }
    }
  
	 
	// in-plane shear state  
    }  else if (strcmp(argv[0],"shear") == 0 || strcmp(argv[0],"shearSpring") == 0) {

		output.tag("Shear");
		theResponse = theShearModel[0]->setResponse(&argv[1], argc-1, output);
		output.endTag();



  } else if (strcmp(argv[0],"RayleighForces") == 0 || strcmp(argv[0],"rayleighForces") == 0) {

    theResponse =  new ElementResponse(this, 12, P);

  }   

  output.endTag();
  return theResponse;
}

int 
Macroelement3d::getResponse(int responseID, Information &eleInfo)
{
  double oneOverL = 1.0/L;

  if (responseID == 1)
    return eleInfo.setVector(this->getResistingForce());

  else if (responseID == 12)
    return eleInfo.setVector(this->getRayleighDampingForces());
    
  else if (responseID == 2) {

	this->getResistingForce();
	for (int i=0; i<12; i++) {
	  q(i) -= q0[i];
    }
	Vector basicForces(12);
	basicForces = q;

	double tmp = basicForces(3);
	for (int ii = 3; ii<11; ii++)
		basicForces(ii) = basicForces(ii + 1);
	basicForces(11) = tmp;

    return eleInfo.setVector(basicForces);
  }

  // Drifts
  else if (responseID == 3) {

     Vector driftsVec(2);
	 driftsVec(0) = driftS;
	 driftsVec(1) = driftF;
     
    return eleInfo.setVector(driftsVec);
  }

  // Crushing state variables
  else if (responseID == 4) {
     Vector shearStateVec(2);

	 Vector sectionResponse(4);

    return eleInfo.setVector(shearStateVec);
  }

  else
    return -1;
}

int 
Macroelement3d::trasformMatrixToGlobal(Matrix& A) {

	int m,j;
	double tmp[18][18];
	double kl[18][18];

	for (m=0; m<18; m++)
		for (j=0;j<18;j++) 
			kl[m][j] = A(m,j);


    static double RWI[3][3];
        
        if (nodeIOffset) {
            // Compute RWI
            RWI[0][0] = -R[0][1]*nodeIOffset[2] + R[0][2]*nodeIOffset[1];
            RWI[1][0] = -R[1][1]*nodeIOffset[2] + R[1][2]*nodeIOffset[1];
            RWI[2][0] = -R[2][1]*nodeIOffset[2] + R[2][2]*nodeIOffset[1];
            
            RWI[0][1] =  R[0][0]*nodeIOffset[2] - R[0][2]*nodeIOffset[0];
            RWI[1][1] =  R[1][0]*nodeIOffset[2] - R[1][2]*nodeIOffset[0];
            RWI[2][1] =  R[2][0]*nodeIOffset[2] - R[2][2]*nodeIOffset[0];
            
            RWI[0][2] = -R[0][0]*nodeIOffset[1] + R[0][1]*nodeIOffset[0];
            RWI[1][2] = -R[1][0]*nodeIOffset[1] + R[1][1]*nodeIOffset[0];
            RWI[2][2] = -R[2][0]*nodeIOffset[1] + R[2][1]*nodeIOffset[0];
        }
        
        static double RWJ[3][3];
        
        if (nodeJOffset) {
            // Compute RWJ
            RWJ[0][0] = -R[0][1]*nodeJOffset[2] + R[0][2]*nodeJOffset[1];
            RWJ[1][0] = -R[1][1]*nodeJOffset[2] + R[1][2]*nodeJOffset[1];
            RWJ[2][0] = -R[2][1]*nodeJOffset[2] + R[2][2]*nodeJOffset[1];
            
            RWJ[0][1] =  R[0][0]*nodeJOffset[2] - R[0][2]*nodeJOffset[0];
            RWJ[1][1] =  R[1][0]*nodeJOffset[2] - R[1][2]*nodeJOffset[0];
            RWJ[2][1] =  R[2][0]*nodeJOffset[2] - R[2][2]*nodeJOffset[0];
            
            RWJ[0][2] = -R[0][0]*nodeJOffset[1] + R[0][1]*nodeJOffset[0];
            RWJ[1][2] = -R[1][0]*nodeJOffset[1] + R[1][1]*nodeJOffset[0];
            RWJ[2][2] = -R[2][0]*nodeJOffset[1] + R[2][1]*nodeJOffset[0];
        }
        
        // Transform local stiffness to global system
        // First compute kl*T_{lg}

        for (m = 0; m < 18; m++) {
            tmp[m][0] = kl[m][0]*R[0][0] + kl[m][1]*R[1][0]  + kl[m][2]*R[2][0];
            tmp[m][1] = kl[m][0]*R[0][1] + kl[m][1]*R[1][1]  + kl[m][2]*R[2][1];
            tmp[m][2] = kl[m][0]*R[0][2] + kl[m][1]*R[1][2]  + kl[m][2]*R[2][2];
            
            tmp[m][3] = kl[m][3]*R[0][0] + kl[m][4]*R[1][0]  + kl[m][5]*R[2][0];
            tmp[m][4] = kl[m][3]*R[0][1] + kl[m][4]*R[1][1]  + kl[m][5]*R[2][1];
            tmp[m][5] = kl[m][3]*R[0][2] + kl[m][4]*R[1][2]  + kl[m][5]*R[2][2];
            
            if (nodeIOffset) {
                tmp[m][3]  += kl[m][0]*RWI[0][0]  + kl[m][1]*RWI[1][0]  + kl[m][2]*RWI[2][0];
                tmp[m][4]  += kl[m][0]*RWI[0][1]  + kl[m][1]*RWI[1][1]  + kl[m][2]*RWI[2][1];
                tmp[m][5]  += kl[m][0]*RWI[0][2]  + kl[m][1]*RWI[1][2]  + kl[m][2]*RWI[2][2];
            }
            
            tmp[m][6] = kl[m][6]*R[0][0] + kl[m][7]*R[1][0]  + kl[m][8]*R[2][0];
            tmp[m][7] = kl[m][6]*R[0][1] + kl[m][7]*R[1][1]  + kl[m][8]*R[2][1];
            tmp[m][8] = kl[m][6]*R[0][2] + kl[m][7]*R[1][2]  + kl[m][8]*R[2][2];
            
            tmp[m][9]  = kl[m][9]*R[0][0] + kl[m][10]*R[1][0] + kl[m][11]*R[2][0];
            tmp[m][10] = kl[m][9]*R[0][1] + kl[m][10]*R[1][1] + kl[m][11]*R[2][1];
            tmp[m][11] = kl[m][9]*R[0][2] + kl[m][10]*R[1][2] + kl[m][11]*R[2][2];
            
            if (nodeJOffset) {
                tmp[m][9]   += kl[m][6]*RWJ[0][0]  + kl[m][7]*RWJ[1][0]  + kl[m][8]*RWJ[2][0];
                tmp[m][10]  += kl[m][6]*RWJ[0][1]  + kl[m][7]*RWJ[1][1]  + kl[m][8]*RWJ[2][1];
                tmp[m][11]  += kl[m][6]*RWJ[0][2]  + kl[m][7]*RWJ[1][2]  + kl[m][8]*RWJ[2][2];
            }

			tmp[m][12] = kl[m][12]*R[0][0] + kl[m][13]*R[1][0]  + kl[m][14]*R[2][0];
            tmp[m][13] = kl[m][12]*R[0][1] + kl[m][13]*R[1][1]  + kl[m][14]*R[2][1];
            tmp[m][14] = kl[m][12]*R[0][2] + kl[m][13]*R[1][2]  + kl[m][14]*R[2][2];
            
            tmp[m][15] = kl[m][15]*R[0][0] + kl[m][16]*R[1][0] + kl[m][17]*R[2][0];
            tmp[m][16] = kl[m][15]*R[0][1] + kl[m][16]*R[1][1] + kl[m][17]*R[2][1];
            tmp[m][17] = kl[m][15]*R[0][2] + kl[m][16]*R[1][2] + kl[m][17]*R[2][2];
            
        }
        
        // Now compute T'_{lg}*(kl*T_{lg})
        for (m = 0; m < 18; m++) {
            K(0,m) = R[0][0]*tmp[0][m] + R[1][0]*tmp[1][m]  + R[2][0]*tmp[2][m];
            K(1,m) = R[0][1]*tmp[0][m] + R[1][1]*tmp[1][m]  + R[2][1]*tmp[2][m];
            K(2,m) = R[0][2]*tmp[0][m] + R[1][2]*tmp[1][m]  + R[2][2]*tmp[2][m];
            
            K(3,m) = R[0][0]*tmp[3][m] + R[1][0]*tmp[4][m]  + R[2][0]*tmp[5][m];
            K(4,m) = R[0][1]*tmp[3][m] + R[1][1]*tmp[4][m]  + R[2][1]*tmp[5][m];
            K(5,m) = R[0][2]*tmp[3][m] + R[1][2]*tmp[4][m]  + R[2][2]*tmp[5][m];
            
            if (nodeIOffset) {
                K(3,m) += RWI[0][0]*tmp[0][m]  + RWI[1][0]*tmp[1][m] + RWI[2][0]*tmp[2][m];
                K(4,m) += RWI[0][1]*tmp[0][m]  + RWI[1][1]*tmp[1][m] + RWI[2][1]*tmp[2][m];
                K(5,m) += RWI[0][2]*tmp[0][m]  + RWI[1][2]*tmp[1][m] + RWI[2][2]*tmp[2][m];
            }
            
            K(6,m) = R[0][0]*tmp[6][m] + R[1][0]*tmp[7][m]  + R[2][0]*tmp[8][m];
            K(7,m) = R[0][1]*tmp[6][m] + R[1][1]*tmp[7][m]  + R[2][1]*tmp[8][m];
            K(8,m) = R[0][2]*tmp[6][m] + R[1][2]*tmp[7][m]  + R[2][2]*tmp[8][m];
            
            K(9,m)  = R[0][0]*tmp[9][m] + R[1][0]*tmp[10][m] + R[2][0]*tmp[11][m];
            K(10,m) = R[0][1]*tmp[9][m] + R[1][1]*tmp[10][m] + R[2][1]*tmp[11][m];
            K(11,m) = R[0][2]*tmp[9][m] + R[1][2]*tmp[10][m] + R[2][2]*tmp[11][m];
            
            if (nodeJOffset) {
                K(9,m)  += RWJ[0][0]*tmp[6][m]  + RWJ[1][0]*tmp[7][m] + RWJ[2][0]*tmp[8][m];
                K(10,m) += RWJ[0][1]*tmp[6][m]  + RWJ[1][1]*tmp[7][m] + RWJ[2][1]*tmp[8][m];
                K(11,m) += RWJ[0][2]*tmp[6][m]  + RWJ[1][2]*tmp[7][m] + RWJ[2][2]*tmp[8][m];
            }

			K(12,m) = R[0][0]*tmp[12][m] + R[1][0]*tmp[13][m]  + R[2][0]*tmp[14][m];
            K(13,m) = R[0][1]*tmp[12][m] + R[1][1]*tmp[13][m]  + R[2][1]*tmp[14][m];
            K(14,m) = R[0][2]*tmp[12][m] + R[1][2]*tmp[13][m]  + R[2][2]*tmp[14][m];
            
            K(15,m) = R[0][0]*tmp[15][m] + R[1][0]*tmp[16][m] + R[2][0]*tmp[17][m];
            K(16,m) = R[0][1]*tmp[15][m] + R[1][1]*tmp[16][m] + R[2][1]*tmp[17][m];
            K(17,m) = R[0][2]*tmp[15][m] + R[1][2]*tmp[16][m] + R[2][2]*tmp[17][m];

        }

	return 0;
};









// damping methods. 
int
Macroelement3d::setRayleighDampingFactors(double alpham, double betak, double betak0, double betakc)
{
	alphaM = alpham;
	betaK = betak;
	betaK0 = betak0;
	betaKc = betakc;

	// check that memory has been allocated to store compute/return
	// damping matrix & residual force calculations
	if (index == -1) {
		int numDOF = this->getNumDOF();

		for (int i = 0; i<numMatrices; i++) {
			Matrix *aMatrix = theMatrices[i];
			if (aMatrix->noRows() == numDOF) {
				index = i;
				i = numMatrices;
			}
		}
		if (index == -1) {
			Matrix **nextMatrices = new Matrix *[numMatrices + 1];
			if (nextMatrices == 0) {
				opserr << "Element::getTheMatrix - out of memory\n";
			}
			int j;
			for (j = 0; j<numMatrices; j++)
				nextMatrices[j] = theMatrices[j];
			Matrix *theMatrix = new Matrix(numDOF, numDOF);
			if (theMatrix == 0) {
				opserr << "Element::getTheMatrix - out of memory\n";
				exit(-1);
			}
			nextMatrices[numMatrices] = theMatrix;

			Vector **nextVectors1 = new Vector *[numMatrices + 1];
			Vector **nextVectors2 = new Vector *[numMatrices + 1];
			if (nextVectors1 == 0 || nextVectors2 == 0) {
				opserr << "Element::getTheVector - out of memory\n";
				exit(-1);
			}

			for (j = 0; j<numMatrices; j++) {
				nextVectors1[j] = theVectors1[j];
				nextVectors2[j] = theVectors2[j];
			}

			Vector *theVector1 = new Vector(numDOF);
			Vector *theVector2 = new Vector(numDOF);
			if (theVector1 == 0 || theVector2 == 0) {
				opserr << "Element::getTheVector - out of memory\n";
				exit(-1);
			}

			nextVectors1[numMatrices] = theVector1;
			nextVectors2[numMatrices] = theVector2;

			if (numMatrices != 0) {
				delete[] theMatrices;
				delete[] theVectors1;
				delete[] theVectors2;
			}
			index = numMatrices;
			numMatrices++;
			theMatrices = nextMatrices;
			theVectors1 = nextVectors1;
			theVectors2 = nextVectors2;
		}
	}

	// if need storage for Kc go get it
	if (betaKc != 0.0) {
		if (Kc == 0)
			Kc = new Matrix(this->getTangentStiff());
		if (Kc == 0) {
			opserr << "WARNING - ELEMENT::setRayleighDampingFactors - out of memory\n";
			betaKc = 0.0;
		}

		// if don't need storage for Kc & have allocated some for it, free the memory
	}
	else if (Kc != 0) {
		delete Kc;
		Kc = 0;
	}

	return 0;
}



const Matrix &
Macroelement3d::getDamp(void)
{
	if (index == -1) {
		this->setRayleighDampingFactors(alphaM, betaK, betaK0, betaKc);
	}

	// now compute the damping matrix
	Matrix *theMatrix = theMatrices[index];
	theMatrix->Zero();
	if (alphaM != 0.0)
		theMatrix->addMatrix(0.0, this->getMass(), alphaM);
	if (betaK != 0.0)
		theMatrix->addMatrix(1.0, this->getTangentStiff(), betaK);
	if (betaK0 != 0.0)
		if (dampingModel == 0) {
			theMatrix->addMatrix(1.0, this->getInitialStiff(), betaK0);
		}
		else {
			if (ductilityDemandCommitted < 1) {
				theMatrix->addMatrix(1.0, this->getInitialStiff(), betaK0);
			}
			else {
				if (dampingModel == 1) {
					theMatrix->addMatrix(1.0, this->getInitialStiff(), betaK0 / ductilityDemandCommitted);
				}
				else {
					theMatrix->addMatrix(1.0, this->getInitialStiff(), betaK0 / ductilityDemandCycle);
				}
			}
		}
	if (betaKc != 0.0)
		theMatrix->addMatrix(1.0, this->getSecantStiff(), betaKc);

	// return the computed matrix
	return *theMatrix;
}


const Vector &
Macroelement3d::getRayleighDampingForces(void)
{

	if (index == -1) {
		this->setRayleighDampingFactors(alphaM, betaK, betaK0, betaKc);
	}

	Matrix *theMatrix = theMatrices[index];
	Vector *theVector = theVectors2[index];
	Vector *theVector2 = theVectors1[index];

	//
	// perform: R = (alphaM * M + betaK0 * K0 + betaK * K) * v
	//            = D * v
	//

	// determine the vel vector from ele nodes
	Node **theNodes = this->getNodePtrs();
	int numNodes = this->getNumExternalNodes();
	int loc = 0;
	for (int i = 0; i<numNodes; i++) {
		const Vector &vel = theNodes[i]->getTrialVel();
		for (int i = 0; i<vel.Size(); i++) {
			(*theVector2)(loc++) = vel[i];
		}
	}

	// now compute the damping matrix
	theMatrix->Zero();
	if (alphaM != 0.0)
		theMatrix->addMatrix(0.0, this->getMass(), alphaM);
	if (betaK != 0.0)
		theMatrix->addMatrix(1.0, this->getTangentStiff(), betaK);
	if (betaK0 != 0.0)
		if (dampingModel == 0) {
			theMatrix->addMatrix(1.0, this->getInitialStiff(), betaK0);
		}
		else {
			if (ductilityDemandCommitted < 1.0) {
				theMatrix->addMatrix(1.0, this->getInitialStiff(), betaK0);
			}
			else {
				if (dampingModel == 1) {
					theMatrix->addMatrix(1.0, this->getInitialStiff(), betaK0 / ductilityDemandCommitted);
				}
				else {
					theMatrix->addMatrix(1.0, this->getInitialStiff(), betaK0 / ductilityDemandCycle);
				}
			}
		}
	if (betaKc != 0.0)
		theMatrix->addMatrix(1.0, this->getSecantStiff(), betaKc);

	// finally the D * v
	theVector->addMatrixVector(0.0, *theMatrix, *theVector2, 1.0);

	return *theVector;
}



const Matrix&
Macroelement3d::getSecantStiff()
{
	static Matrix kb(12, 12);

	kb.Zero();
	q.Zero();

	double oneOverL = 1.0 / L;

	// first interface
	int order = theSections[0]->getOrder();
	const ID &code0 = theSections[0]->getType();

	int ordering[4];
	for (int j = 0; j < 4; j++)
		ordering[j] = -1;

	for (int j = 0; j < order; j++) {
		switch (code0(j)) {
		case SECTION_RESPONSE_P:     ordering[0] = j;	break;
		case SECTION_RESPONSE_MZ:    ordering[1] = j;   break;
		case SECTION_RESPONSE_MY:    ordering[2] = j;	break;
		case SECTION_RESPONSE_T:     ordering[3] = j;	break;
		default: break;
		}
	}

	const Vector &s0 = theSections[0]->getStressResultant();
	const Matrix &ks0 = theSections[0]->getSectionTangent();

	// i index: line;  j index: column
	for (int i = 0; i<4; i++) {
		if (ordering[i] >= 0) {
			q(i) = s0(ordering[i]);
			for (int j = 0; j < 4; j++)
				if (i < 3)
					kb(i, j) = ks0(ordering[i], ordering[j]) / intLength(0);
				else
					kb(i, j) = ks0(ordering[i], ordering[j]) / (2.*L);
		}
	}

	// second interface
	order = theSections[1]->getOrder();
	const ID &code1 = theSections[1]->getType();

	for (int j = 0; j < 4; j++)
		ordering[j] = -1;

	for (int j = 0; j < order; j++) {
		switch (code1(j)) {
		case SECTION_RESPONSE_P:     ordering[0] = j;	break;
		case SECTION_RESPONSE_MZ:    ordering[1] = j; break;
		case SECTION_RESPONSE_MY:    ordering[2] = j;	break;
		default: break;
		}
	}

	const Vector &s1 = theSections[1]->getStressResultant();
	const Matrix &ks1 = theSections[1]->getSectionTangent();

	// i index: line;  j index: column
	for (int i = 0; i<3; i++) {
		if (ordering[i] >= 0) {
			q(4 + i) = s1(ordering[i]);
			for (int j = 0; j < 3; j++)
				kb(4 + i, 4 + j) = ks1(ordering[i], ordering[j]) / intLength(1);
		}
	}


	// third interface
	order = theSections[2]->getOrder();
	const ID &code2 = theSections[2]->getType();

	for (int j = 0; j < 4; j++)
		ordering[j] = -1;

	for (int j = 0; j < order; j++) {
		switch (code2(j)) {
		case SECTION_RESPONSE_P:     ordering[0] = j;	break;
		case SECTION_RESPONSE_MZ:    ordering[1] = j; break;
		case SECTION_RESPONSE_MY:    ordering[2] = j;	break;
		default: break;
		}
	}

	const Vector &s2 = theSections[2]->getStressResultant();
	const Matrix &ks2 = theSections[2]->getSectionTangent();

	// i index: line;  j index: column
	for (int i = 0; i<3; i++) {
		if (ordering[i] >= 0) {
			q(7 + i) = s2(ordering[i]);
			for (int j = 0; j < 3; j++)
				kb(7 + i, 7 + j) = ks2(ordering[i], ordering[j]) / intLength(2);
		}
	}

	// shear interface 1 and 2
	for (int sect = 0; sect<2; sect++) {

		Matrix EshearModel = theShearModel[sect]->getInitialTangent();
		double N = q(4);
		double Esm = EshearModel(0, 0);

		const Matrix &ks = theShearModel[sect]->getTangent();
		const Vector &s = theShearModel[sect]->getStress();

		kb(10 + sect, 10 + sect) = ks(1, 1);
		q(10 + sect) = s(1);

		// update non-diagonal terms of the stiffness matrix (dependency on N2)
		for (int j = 0; j<3; j++)
			kb(10 + sect, 4 + j) += ks(1, 0) / Esm * kb(4, 4 + j);

	}

	// kill lateral capacity if failed
	if (collapsedS || collapsedF) {
		double failureFactor;
		if (collapsedS)
			failureFactor = failureFactorS;
		else
			failureFactor = failureFactorF;

		// dofs to kill: all
		int dofsToKill[12] = { 0,1,2,3,4,5,6,7,8,9,10,11 };
		int dof;

		for (int kDof = 0; kDof<12; kDof++) {
			dof = dofsToKill[kDof];
			q(dof) *= failureFactor;
			for (int j = 0; j<kb.noCols(); j++)
				kb(dof, j) *= failureFactor;
		}
	}
	else if (failedS || failedF) {

		double failureFactor;
		if (failedS)
			failureFactor = failureFactorS;
		else
			failureFactor = failureFactorF;

		// dofs to kill: 1,2, 5,6, 8,9,  10,11  (moments of every section, shear of every section)
		int dofsToKill[8] = { 1,2, 5,6, 8,9,  10,11 };
		int dof;

		for (int kDof = 0; kDof<8; kDof++) {
			dof = dofsToKill[kDof];
			q(dof) *= failureFactor;
			for (int j = 0; j<kb.noCols(); j++)
				kb(dof, j) *= failureFactor;
		}
	}


	// Transform to local stiffness matrix

	Matrix Klocal(18, 18);
	//if (PDelta==1)   getIncrementalCompatibilityMatrix(true);
	Klocal.addMatrixTripleProduct(0.0, GammaC, kb, 1.0);

	//opserr << "local stiffness = " << Klocal;

	// add geometric local stiffness matrix  
	if (PDelta == 1) {
		Matrix kgeom(18, 18);

		// add fixed end forces. Update
		for (int i = 0; i<12; i++) {
			q(i) += q0[i];
		}

		double qq0 = q(0);
		double qq7 = q(7);
		double qq4 = q(4);

		kgeom(1, 1) = qq0*oneOverL;
		kgeom(1, 13) = -(qq0*oneOverL);
		kgeom(2, 2) = qq0*oneOverL;
		kgeom(2, 14) = -(qq0*oneOverL);
		kgeom(7, 7) = qq7*oneOverL;
		kgeom(7, 16) = -(qq7*oneOverL);
		kgeom(8, 8) = qq7*oneOverL;
		kgeom(8, 17) = -(qq7*oneOverL);
		kgeom(13, 1) = -(qq0*oneOverL);
		kgeom(13, 13) = qq0*oneOverL;
		kgeom(14, 2) = -(qq0*oneOverL);
		kgeom(14, 14) = qq0*oneOverL;
		kgeom(16, 7) = -(qq7*oneOverL);
		kgeom(16, 16) = qq7*oneOverL;
		kgeom(17, 8) = -(qq7*oneOverL);
		kgeom(17, 17) = qq7*oneOverL;

		Klocal += kgeom;
	}


	// Transform to global stiffness matrix   Kglobal = Tgl' * Klocal * Tgl
	// done now otherwise the secant stiffness matrix would be rrturned
	trasformMatrixToGlobal(Klocal);

	return K;
}

