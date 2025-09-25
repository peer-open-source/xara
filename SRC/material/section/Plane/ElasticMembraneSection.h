// Code written/implemented by: Maria Jose Nunez
//
// Created: 03/2024
// 
// Description: This file contains the ElasticMembraneSection class definition
// A ElasticMembraneSection is a subclass of the sectionForceDeformation class and corresponds to the abstract representation
// for the stress-strain behavior for a elastic membrane section in the Finite Element Method or Structural Analysis. 
//
// Reference:
// 1. Rojas, F., Anderson, J. C., Massone, L. M. (2016). A nonlinear quadrilateral layered membrane element with drilling degrees of freedom for 
// the modeling of reinforced concrete walls. Engineering Structures, 124, 521-538.
//
// Source: \OpenSees\SRC\material\section\LayeredMembraneSection
//
// Rev: 1.0

#ifndef ElasticMembraneSection_h
#define ElasticMembraneSection_h

#include <stdlib.h>
#include <math.h>

#include <Vector.h>
#include <Matrix.h>
#include <ID.h>

#include <SectionForceDeformation.h>

class ElasticMembraneSection : public SectionForceDeformation {
public:

	ElasticMembraneSection(int tag,								// section tag
		double elasticModulus,									// Young's Modulus
		double poisson,										    // Poisson ratio
		double thickness,                                       // section thickness
	    double rho = 0.0);                                      // mass density	(mass per unit 3D volume)									

	ElasticMembraneSection();

	~ElasticMembraneSection();

	int setTrialSectionDeformation(const Vector& newTrialSectionStrain);

	// Public methods to obtain strain, stress, tangent and residual information
	const Vector& getSectionDeformation();
	const Vector& getStressResultant();
	const Matrix& getSectionTangent();
	const Matrix& getInitialTangent();

	// Public methods to obtain a copy of section and other information
	SectionForceDeformation* getCopy();
	const char* getClassType() const { return "ElasticMembraneSection"; }
	const ID& getType();
	int getOrder() const;
	
	// Public methods to set the state of the section
	int commitState();
	int revertToLastCommit();
	int revertToStart();

	// Public methods for output
	int sendSelf(int commitTag, Channel&);
	int recvSelf(int commitTag, Channel&, FEM_ObjectBroker&);

	void Print(OPS_Stream&, int flag);

	// Functions used for recorders
	const Vector& getCommittedStrain();
	const Vector& getCommittedStress();

	// Density per unit area
	double getRho();															

private:

	double E;																        // Store the elastic modulus
	double nu;                                                                      // Store the Poisson ratio
	double t;																		// Store the thickness of the section 
	double rho;																		// Store the mass per unit 2D area 

	// Committed State Variables
	Vector CSectionStrain;															// Store the commit strains of the membrane section
	Vector CSectionStress;															// Store the commit stress of the membrane section
	Matrix CSectionTangent;															// Store the commit tangent of the membrane section
	// Trial State Variables
	Vector TSectionStrain;															// Store the trial strains of the membrane section
	Vector TSectionStress;															// Store the trial stress of the membrane section
	Matrix TSectionTangent;															// Store the trial tangent of the membrane section
	Matrix InitialTangent;															// Store the initial tangent of the membrane section

	static ID array;
};

#endif // !ElasticMembraneSection_h
