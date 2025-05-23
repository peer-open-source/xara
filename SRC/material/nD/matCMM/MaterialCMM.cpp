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

#include "MaterialCMM.h"
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <MaterialResponse.h>

#ifdef _HAVE_MAT_CMM

#ifdef _WIN32

#else
#define MATCMM matcmm_
#endif

extern "C" 
{

  void MAT_CMM(double *Stress,double *Strain, double *dStrain, int *kLayer,
			double *Mat_Par, double *Stress1, double *dsdePl, double *ustatev);

}

#else
  void MAT_CMM(double *Stress,double *Strain, double *dStrain, int *kLayer,
	       double *Mat_Par, double *Stress1, double *dsdePl, double *ustatev)
   
{
  opserr << "MAT_CMM - NOT DEFINED IN THIS VERSION, SOURCE CODE RESTRICTED\n";
}
#endif

#include <elementAPI.h>

void * OPS_ADD_RUNTIME_VPV(OPS_MaterialCMM)
{
  NDMaterial *theMaterial = 0;
  
  int numArgs = OPS_GetNumRemainingInputArgs();
  
  if (numArgs < 2 + MaterialCMM_NumParameters) {
    opserr << "Want: insufficient args: nDMaterial MaterialCMMc $tag $layer 71 parameters!" << endln;
    return 0;	
  }
  
  int iData[2];
  double dData[MaterialCMM_NumParameters];
  
  int numData = 2;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag\n";
    return 0;
  }
  
  numData = MaterialCMM_NumParameters;
  
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING invalid data: " << iData[0] <<"\n";
    return 0;
  }  
  
  theMaterial = new MaterialCMM(iData[0], iData[1], dData);
	
  return theMaterial;
}


//null constructor
MaterialCMM::MaterialCMM( ) : 
NDMaterial(0, ND_TAG_MaterialCMM ), 
stress(5), strain(5), tangent(5,5)
{ 

}


//full constructor
MaterialCMM::MaterialCMM(int tag, int theLayer, double *matParameters) :
NDMaterial( tag, ND_TAG_MaterialCMM ),
stress(5), strain(5), tangent(5,5), layer(theLayer)
{

  for (int i=0; i<MaterialCMM_NumParameters; i++)
    matPar[i] = matParameters[i];

  for (int i=0; i<MaterialCMM_NumStateVar; i++) {
    stateVarC[i] = 0.;
    stateVarT[i] = 0.;
  }

  for (int i=0; i<5; i++) {
    stressC[i] =0.;
    stressT[i] =0.;
    strainC[i] =0.;
    strainT[i] =0.;
  }
  // set initial
  for (int i=0; i<9; i++) {      
    tangentC[i]=0.;
    tangentT[i]=0.;
  }
}

//destructor
MaterialCMM::~MaterialCMM( ) 
{ 

} 


NDMaterial*
MaterialCMM::getCopy( ) 
{
  MaterialCMM *clone ;   //new instance of this class

  clone = new MaterialCMM(this->getTag(), layer, matPar);

  return clone ;
}


NDMaterial* 
MaterialCMM::getCopy( const char *type ) 
{
  if (strcmp(type, this->getType()) == 0)
    return this->getCopy( ) ;
  else
    return 0;
}


int 
MaterialCMM::getOrder( ) const
{
  return 5 ;
}


const char*
MaterialCMM::getType( ) const 
{
  return "PlaneStress" ; 
}



//swap history variables
int 
MaterialCMM::commitState( ) 
{
  for (int i=0; i<5; i++) {
    stressC[i] = stressT[i];
    strainC[i] = strainT[i];
  }

  for (int i=0; i<9; i++)
    tangentC[i] = tangentT[i];

  for (int i=0; i<MaterialCMM_NumStateVar; i++)
    stateVarC[i] = stateVarT[i];

  return 0;
}


int 
MaterialCMM::revertToLastCommit( )
{
  return 0;
}


int
MaterialCMM::revertToStart( )
{
  for (int i=0; i<MaterialCMM_NumStateVar; i++) {
    stateVarC[i] = 0.;
    stateVarT[i] = 0.;
  }

  for (int i=0; i<5; i++) {
    stressC[i] =0.;
    stressT[i] =0.;
    strainC[i] =0.;
    strainT[i] =0.;
  }
  // set initial
  for (int i=0; i<9; i++) {      
    tangentC[i]=0.;
    tangentT[i]=0.;
  }
  return 0;
}

//receive the strain
int 
MaterialCMM::setTrialStrain( const Vector &strainIN )
{
  strain = strainIN;

  for (int i=0; i<5; i++) {
    strainT[i] = strainIN(i);
    stressT[i] = stressC[i];
    dStrain[i] = strainT[i] - strainC[i];
  }
  for (int i=0; i< MaterialCMM_NumStateVar; i++) {
    stateVarT[i] = stateVarC[i];
  }


  MAT_CMM(stressC, strainC, dStrain, &layer, matPar, stressT, tangentT, stateVarT);

  
  for (int i = 0; i < 5; i++) {
    stress(i) = stressT[i];
    
    for (int j = 0; j < 5; j++)
      tangent(i,j) = tangentT[j + 5 * i];
  }

  return 0;
}


const Vector& 
MaterialCMM::getStrain( )
{
  return strain ;
}


const Vector&  
MaterialCMM::getStress( )
{
  return stress ;
}


const Matrix&  
MaterialCMM::getTangent( )
{
  
  tangent.setData(tangentT,5,5);
  return tangent;
}

const Matrix&  
MaterialCMM::getInitialTangent( )
{
  tangent.setData(tangentI,5,5);
  return tangent;
}


void  
MaterialCMM::Print( OPS_Stream &s, int flag )
{
  s << "MaterialCMM Material tag: " << this->getTag() << endln ;
}


int 
MaterialCMM::sendSelf(int commitTag, Channel &theChannel) 
{
  int res = 0;

  return res;
}

int 
MaterialCMM::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;

  return res;
}
 
