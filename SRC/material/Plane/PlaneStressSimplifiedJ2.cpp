// -------------------
//
// Written: Quan Gu and Zhijian Qiu
// Created: 2013.7
//
// Reference: JP Conte, MK. Jagannath, Seismic relibility analysis of concrete 
// gravity dams, A Report on Research, Rice University, 1995. 
//           EA de Souza Neto, D Peri´c, DRJ Owen, Computational methods for 
// plasticity, Theory and applications (see pages 357 to 366), 2008.
// 
// 3D J2 plasticity model with linear isotropic and kinematic hardening
//
// -------------------

#include <math.h>   
#include <stdlib.h>
#include <PlaneStressSimplifiedJ2.h>
#include <Information.h>
#include <ID.h>
#include <MaterialResponse.h>
#include <Parameter.h>

#include <Hash.h>
using namespace OpenSees::Hash::literals;
// #include <classTags.h>
#define ND_TAG_PlaneStress   3452


Matrix PlaneStressSimplifiedJ2::tmpMatrix(3,3);
Vector PlaneStressSimplifiedJ2::tmpVector(3);

// --- element: eps(1,1),eps(2,2),eps(3,3),2*eps(1,2),2*eps(2,3),2*eps(1,3) ----
// --- material strain: eps(1,1),eps(2,2),eps(3,3),eps(1,2),eps(2,3),eps(1,3) , same sign ----

// be careful! Here we use  eps(1,1),eps(2,2),2*eps(1,2). i.e., the same as that of element. 

#include <SimplifiedJ2.h>


PlaneStressSimplifiedJ2::PlaneStressSimplifiedJ2 (int pTag, 
						   int nd, 
						   NDMaterial &passed3DMaterial)
  : NDMaterial(pTag, ND_TAG_PlaneStress), stress(3),
    strain(3), Cstress(3), Cstrain(3),theTangent(3,3)
{
  this->ndm = 2;    
  the3DMaterial = passed3DMaterial.getCopy();  
  
  stress.Zero();
  strain.Zero();
  
  Cstress.Zero();
  Cstrain.Zero();
  
  savedStrain33=0.0;
  CsavedStrain33 = 0.0;
}



PlaneStressSimplifiedJ2::~PlaneStressSimplifiedJ2() {
	
	return; 
};

     

	
int PlaneStressSimplifiedJ2::plastIntegrator(){

	int maxIter = 25;
	double tol = 1e-12;
	double e33 = CsavedStrain33;

	static Vector strain3D(6);
	static Vector stress3D(6);
	static Matrix tangent3D(6,6);
	
	strain3D(0) = strain(0);
	strain3D(1) = strain(1);
	strain3D(2) = e33;
	strain3D(3) = strain(2);
	strain3D(4) = 0.0;
	strain3D(5) = 0.0;

	the3DMaterial->setTrialStrain(strain3D);
	stress3D = the3DMaterial->getStress();
	tangent3D = the3DMaterial->getTangent();


	int i =0;

	double e33_old=e33+1.0;

	while (( fabs(e33-e33_old)>tol) &&( fabs(stress3D(2))>tol) &&(i<maxIter)) {
	    e33_old = e33;		
		e33 -= stress3D(2)/tangent3D(2,2);
		strain3D(2) = e33;
	    the3DMaterial->setTrialStrain(strain3D);
	    stress3D = the3DMaterial->getStress();
		tangent3D = the3DMaterial->getTangent();
		i++;
	} 

	if (( fabs(e33-e33_old)>tol) &&(fabs(stress3D(2))>tol)) {
		opserr<<"PlaneStressSimplifiedJ2::plastIntegrator() can not find e33"<<endln;
		return -1;
	}

	// --------- update the stress and tangent -----
	savedStrain33 = e33;

   stress(0) = stress3D(0);
   stress(1) = stress3D(1);
   stress(2) = stress3D(3);
   

   double D22 = tangent3D(2,2);
   static Vector D12(3);
   static Vector D21(3);
   static Matrix D11(3,3);

 D11(0,0)=tangent3D(0,0);
 D11(0,1)=tangent3D(0,1);
 D11(0,2)=tangent3D(0,3);
 D11(1,0)=tangent3D(1,0);
 D11(1,1)=tangent3D(1,1);
 D11(1,2)=tangent3D(1,3);
 D11(2,0)=tangent3D(3,0);
 D11(2,1)=tangent3D(3,1);
 D11(2,2)=tangent3D(3,3);

D12(0) = tangent3D(0,2);
D12(1) = tangent3D(1,2);
D12(2) = tangent3D(3,2);

D21(0) = tangent3D(2,0);
D21(1) = tangent3D(2,1);
D21(2) = tangent3D(2,3);

for( int i=0; i<3; i++)
  for (int j=0; j<3; j++)
	  theTangent(i,j) = D11(i,j)-1.0/D22*D12(i)*D21(j);

 
  return 0;

}
 

int
PlaneStressSimplifiedJ2::setTrialStrain (const Vector &pStrain)
{

    strain = pStrain;
 
  // ----- change to real strain instead of eng. strain

  // strain[2] /=2.0;     be careful!           

  return this->plastIntegrator();
}

int
PlaneStressSimplifiedJ2::setTrialStrain(const Vector &v, const Vector &r)
{
	return this->setTrialStrain ( v);
}

int PlaneStressSimplifiedJ2::setTrialStrainIncr(const Vector &v){
	
	// ----- change to real strain instead of eng. strain
   // ---- since all strain in material is the true strain, not eng.strain. 

		strain[0] = Cstrain[0]+v[0];
		strain[1] = Cstrain[1]+v[1];
		strain[2] = Cstrain[2]+v[2];     //  no need to divide by 2.0;
	  
	 this->plastIntegrator();

	 return 0;

};

int
PlaneStressSimplifiedJ2::setTrialStrainIncr(const Vector &v, const Vector &r)
{
	return this->setTrialStrainIncr(v);
}

// Calculates current tangent stiffness.
const Matrix & 
PlaneStressSimplifiedJ2::getTangent()
{
	return theTangent;
}

const Matrix & PlaneStressSimplifiedJ2::getInitialTangent()
{
	return this->getTangent();
}
        
     
const Vector & PlaneStressSimplifiedJ2::getStress()
{
  return stress;
}

const Vector & PlaneStressSimplifiedJ2::getStrain()
{
  return strain; 
}

const Vector & PlaneStressSimplifiedJ2::getCommittedStress()
{
    return Cstress;
}

const Vector & PlaneStressSimplifiedJ2::getCommittedStrain()
{
    return Cstrain;
}


int PlaneStressSimplifiedJ2::commitState()
{
	CsavedStrain33 = savedStrain33; 
	Cstress = stress;
	Cstrain = strain;
	the3DMaterial->commitState();
	//CcumPlastStrainDev = cumPlastStrainDev;

	return 0;

}

int PlaneStressSimplifiedJ2::revertToLastCommit()
{
// -- to be implemented.
	return 0;
}

int PlaneStressSimplifiedJ2::revertToStart()
{
	// -- to be implemented.
	return 0;
}



NDMaterial * 
PlaneStressSimplifiedJ2::getCopy()
{
    PlaneStressSimplifiedJ2 * theJ2 = new PlaneStressSimplifiedJ2(this->getTag(),this->ndm, *the3DMaterial);
    return theJ2;
}

NDMaterial * PlaneStressSimplifiedJ2::getCopy (const char *type){
  if (strcmp(type,"PlaneStress") == 0) {
    PlaneStressSimplifiedJ2 * theJ2 = new PlaneStressSimplifiedJ2(this->getTag(),this->ndm, *the3DMaterial);
    return theJ2;
  } else {
    return 0;
  }
}
 


int
PlaneStressSimplifiedJ2::sendSelf(int commitTag, Channel &theChannel){
	// -- to be implemented.


	return 0;
}

int PlaneStressSimplifiedJ2::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker){
	// -- to be implemented.
	return 0;
}
  
     
Response * PlaneStressSimplifiedJ2::setResponse (const char **argv, int argc, OPS_Stream &s){


  if (strcmp(argv[0],"stress") == 0 || strcmp(argv[0],"stresses") == 0)
		return new MaterialResponse(this, 1, stress);

  else if (strcmp(argv[0],"strain") == 0 || strcmp(argv[0],"strains") == 0)
		return new MaterialResponse(this, 2, strain);

  else if (strcmp(argv[0],"tangent") == 0 || strcmp(argv[0],"Tangent") == 0)
		return new MaterialResponse(this, 3, theTangent);

   else if (strcmp(argv[0],"strain33") == 0 || strcmp(argv[0],"Strain33") == 0)
		return new MaterialResponse(this, 4, savedStrain33 );

  else
		return 0;
	
}



int
PlaneStressSimplifiedJ2::getResponse(int responseID, Information &matInfo)
{
		
	switch (responseID) {
		case -1:
			return -1;
		case 1:
			if (matInfo.theVector != 0)
				*(matInfo.theVector) =stress;
			return 0;

		case 2:
			if (matInfo.theVector != 0)
				*(matInfo.theVector) = strain;
			return 0;

		case 3:
			if (matInfo.theMatrix != 0)
				*(matInfo.theMatrix) = theTangent;
			return 0;

	 	case 4:
		  //if (matInfo.theDouble != 0)
			    matInfo.setDouble (savedStrain33);
			return 0;



		}
		
 

	return 0;
}

void PlaneStressSimplifiedJ2::Print(OPS_Stream &s, int flag)
{
  the3DMaterial->Print(s,flag);
	// if (flag == OPS_PRINT_PRINTMODEL_JSON) {
	// 	s << TaggedObject::JsonPropertyIndent;
	// 	s << "{";
	// 	s << "\"name\": \"" << this->getTag() << "\", ";
	// 	s << "\"type\": \"PlaneStressSimplifiedJ2\", ";
	// 	s << "\"G\":"
	// 	s << "}";
	// 	return;
	// }
	// return;
}


int PlaneStressSimplifiedJ2::setParameter(const char **argv, int argc, Parameter &param){
  // -- to be implemented.
  return 0;
};

int PlaneStressSimplifiedJ2::updateParameter(int responseID, Information &eleInformation){
  return 0;
};

