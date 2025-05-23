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
//
// Written: QGU. UCSD Group for response sensitivity analysis
//  Contact: Quan Gu(guquan2005@gmail.com)  
//           Joel P. Conte (jpconte@ucsd.edu)
//           Michele Barbato (mbarbato@lsu.edu)
//
// Created: May 2009
//
// refer to paper:
// Conte, J. P., Vijalapura, P., and Meghella, M., "Consistent Finite Element Response Sensitivities Analysis," 
// Journal of Engineering Mechanics, ASCE, Vol. 129, No. 12, pp. 1380-1393, 2003.
//
#include <UniaxialJ2Plasticity.h>
#include <Vector.h>
#include <Channel.h>
#include <Matrix.h>
#include <Information.h>
#include <Parameter.h>
#include <string.h>
#include <math.h>
#include <float.h>


UniaxialJ2Plasticity::UniaxialJ2Plasticity(int pTag, double pE, double pSigmaY,
                                           double pHkin, double pHiso)
:UniaxialMaterial(pTag,MAT_TAG_UniaxialJ2Plasticity),
 E(pE), sigmaY(pSigmaY), Hiso(pHiso), Hkin(pHkin)
{
   parameterID = 0;
   SHVs = 0;

   // Initialize variables
   this->revertToStart();
}

UniaxialJ2Plasticity::~UniaxialJ2Plasticity()
{
    if (SHVs != 0) 
        delete SHVs;
    SHVs =0;
}

int 
UniaxialJ2Plasticity::setTrialStrain(double strain, double strainRate)
{
    t_strain = strain;

    // ------ Elastic trial -------

    TPlasticStrain = CPlasticStrain;
    TBackStress = CBackStress;
    TAccumulatedPlasticStrain = CAccumulatedPlasticStrain;
    t_stress = E * (t_strain - CPlasticStrain);
    double CSigmaY = sigmaY+Hiso*CAccumulatedPlasticStrain;

    // Compute trial stress relative to committed back stress
    double xsi = t_stress - TBackStress;

    // Compute yield criterion
    double f = fabs(xsi) - (sigmaY + Hiso*TAccumulatedPlasticStrain);

    // Elastic step ... no updates required
    if (f <= -DBL_EPSILON * E) {
        // Set trial tangent
        TTangent = E;
    }

    // ------- Plastic step ... perform return mapping algorithm ---
    else {
      double deltaLambda = (fabs(xsi) - CSigmaY)/(E+Hkin+Hiso);

      // Find sign of xsi
      int sign = (xsi < 0) ? -1 : 1;

      TPlasticStrain = CPlasticStrain +deltaLambda*sign;

      TBackStress = CBackStress + Hkin*deltaLambda*sign;

      TAccumulatedPlasticStrain = CAccumulatedPlasticStrain + deltaLambda;

      t_stress = E * (t_strain - TPlasticStrain);
    
      TTangent = E*(Hkin+Hiso) / (E+Hkin+Hiso);
    }


    return 0;
}

double 
UniaxialJ2Plasticity::getStress()
{
    return t_stress;
}

double 
UniaxialJ2Plasticity::getTangent()
{
    return TTangent;
}

double 
UniaxialJ2Plasticity::getStrain()
{
    return t_strain;
}

int 
UniaxialJ2Plasticity::commitState()
{
    // Commit trial history variables
    CPlasticStrain = TPlasticStrain;
    CBackStress = TBackStress;
    CAccumulatedPlasticStrain = TAccumulatedPlasticStrain;

    c_strain = t_strain;
    c_stress = t_stress;
    CTangent = TTangent;
    
    return 0;
}

int 
UniaxialJ2Plasticity::revertToLastCommit()
{
  return 0;
}

int 
UniaxialJ2Plasticity::revertToStart()
{
    // Reset committed history variables
    CPlasticStrain = 0.0;
    CBackStress = 0.0;
    CAccumulatedPlasticStrain = 0.0;

    // Reset committed history variables
    TPlasticStrain = 0.0;
    TBackStress = 0.0;
    TAccumulatedPlasticStrain = 0.0;

    // Initialize state variables
    t_strain = 0.0;
    t_stress = 0.0;
    TTangent = E;

    // Initialize committed variables
    c_strain = 0.0;
    c_stress = 0.0;
    CTangent = E;

    if (SHVs != 0) 
        SHVs->Zero();

    return 0;
}

UniaxialMaterial *
UniaxialJ2Plasticity::getCopy()
{
    UniaxialJ2Plasticity *theCopy =
        new UniaxialJ2Plasticity(this->getTag(), E, sigmaY, Hkin, Hiso);

    // Copy committed history variables
    theCopy->CPlasticStrain = CPlasticStrain;
    theCopy->CBackStress = CBackStress;
    theCopy->CAccumulatedPlasticStrain = CAccumulatedPlasticStrain;

    // Copy trial history variables
    theCopy->TPlasticStrain = TPlasticStrain;
    theCopy->TBackStress = TBackStress;
    theCopy->TAccumulatedPlasticStrain = TAccumulatedPlasticStrain;

    // Copy trial state variables
    theCopy->t_strain = t_strain;
    theCopy->t_stress = t_stress;
    theCopy->TTangent = TTangent;
    theCopy->c_strain = c_strain;
    theCopy->c_stress = c_stress;
    theCopy->CTangent = CTangent;   
    return theCopy;
}

int 
UniaxialJ2Plasticity::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(12);
  
  data(0) = this->getTag();
  data(1) = E;
  data(2) = sigmaY;
  data(3) = Hiso;
  data(4) = Hkin;
  data(5) = CPlasticStrain;
  data(6) = CBackStress;
  data(7) = CAccumulatedPlasticStrain;
  data(8) = t_strain;
  data(9) = t_stress;
  data(10) = TTangent;
  data(11) = c_strain;
  data(12) = c_stress;
  data(13) = CTangent;  
  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "UniaxialJ2Plasticity::sendSelf() - failed to send data\n";

  return res;
}

int 
UniaxialJ2Plasticity::recvSelf(int cTag, Channel &theChannel, 
                   FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static Vector data(12);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
  if (res < 0) {
      opserr << "UniaxialJ2Plasticity::recvSelf() - failed to receive data\n";
      E = 0; 
      this->setTag(0);      
  }
  else {
    this->setTag((int)data(0));
    E = data(1);
    sigmaY = data(2);
    Hiso = data(3);
    Hkin = data(4);
    CPlasticStrain = data(5);
    CBackStress = data(6);
    CAccumulatedPlasticStrain = data(7);
    t_strain = data(8);
    t_stress = data(9);
    TTangent = data(10);
    c_strain = data(11);
    c_stress = data(12);
    CTangent = data(13);  
  }
    
  return res;
}

void 
UniaxialJ2Plasticity::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
      s << OPS_PRINT_JSON_ELEM_INDENT << "{";
      s << "\"name\": \"" << this->getTag() << "\", ";
      s << "\"type\": \"" << this->getClassType() <<"\", ";
      s << "\"E\": " << E << ", ";
      s << "\"Y\": " << sigmaY << ", ";
      s << "\"Hiso\": " << Hiso << ", ";
      s << "\"Hkin\": " << Hkin;
      s << "}";
      return;
    }
    else {
      s << "UniaxialJ2Plasticity, tag: " << this->getTag() << endln;
      s << "  E: " << E << endln;
      s << "  sigmaY: " << sigmaY << endln;
      s << "  Hiso: " << Hiso << endln;
      s << "  Hkin: " << Hkin << endln;
    }
   
}


int
UniaxialJ2Plasticity::setParameter(const char **argv, int argc, Parameter &param)
{
  if (strcmp(argv[0],"sigmaY") == 0 || strcmp(argv[0],"fy") == 0)
    return param.addObject(1, this);

  if (strcmp(argv[0],"E") == 0)
    return param.addObject(2, this);

  if ((strcmp(argv[0],"H_kin") == 0)||(strcmp(argv[0],"Hkin") == 0))
    return param.addObject(3, this);

  if ((strcmp(argv[0],"H_iso") == 0)||(strcmp(argv[0],"Hiso") == 0))
    return param.addObject(4, this);

  return -1;
}


int
UniaxialJ2Plasticity::updateParameter(int parameterID, Information &info)
{
    switch (parameterID) {
    case -1:
        return -1;
    case 1:
        this->sigmaY = info.theDouble;
        break;
    case 2:
        this->E = info.theDouble;
        break;
    case 3:
        this->Hkin = info.theDouble;
        break;
    case 4:
        this->Hiso = info.theDouble;
        break;
    default:
        return -1;
    }

    return 0;
}



int
UniaxialJ2Plasticity::activateParameter(int passedParameterID)
{
    parameterID = passedParameterID;

    return 0;
}




double
UniaxialJ2Plasticity::getStrainSensitivity(int gradIndex)
{
    
    if (SHVs ==0) return 0.0;
    else{
        double sensitivity =(*SHVs)(4,gradIndex-1); // unconditional stress sensitivity
        return sensitivity;
    }
}

double
UniaxialJ2Plasticity::getStressSensitivity(int gradIndex, bool conditional)
{
    
    if (conditional == false) {  // return stress sensitivity for recorder purpose
        if (SHVs ==0) return 0.0;
        else {
            double sensitivity =(*SHVs)(3,gradIndex-1); // unconditional stress sensitivity
            return sensitivity;
        }
    }

    // First set values depending on what is random
    double SigmaYSensitivity = 0.0;
    double ESensitivity = 0.0;
    double HkinSensitivity = 0.0;
    double HisoSensitivity = 0.0;

    if (parameterID == 1) {  // sigmaY
        SigmaYSensitivity = 1.0;
    }
    else if (parameterID == 2) {  // E
        ESensitivity = 1.0;
    }
    else if (parameterID == 3) {  // Hkin
        HkinSensitivity = 1.0;
    }
    else if (parameterID == 4) {  // Hiso
        HisoSensitivity = 1.0;
    }
    else {
        // Nothing random here, but may have to return something in any case
    }

    double t_strainSensitivity = 0.0; 

    // Then pick up history variables for this gradient number
    double CPlasticStrainSensitivity = 0.0;
    double CBackStressSensitivity     = 0.0;
    double CAccumulatedPlasticStrainSensitivity     = 0.0;
    if (SHVs != 0) {
        CPlasticStrainSensitivity = (*SHVs)(0,gradIndex);
        CBackStressSensitivity     = (*SHVs)(1,gradIndex);
        CAccumulatedPlasticStrainSensitivity = (*SHVs)(2,gradIndex);
    }

    // ------ Elastic trial -------

    double deltaLambda = 0.0;
    TPlasticStrain = CPlasticStrain;
    TBackStress = CBackStress;
    TAccumulatedPlasticStrain = CAccumulatedPlasticStrain;
    t_stress = E * (t_strain-CPlasticStrain);
    double t_stressSensitivity = E*(t_strainSensitivity-CPlasticStrainSensitivity)+ESensitivity*(t_strain-CPlasticStrain);
    double CSigmaY = sigmaY+Hiso*CAccumulatedPlasticStrain;

    // Compute trial stress relative to committed back stress
    double xsi = t_stress - TBackStress;

    // Compute yield criterion
    double f = fabs(xsi) - (sigmaY + Hiso*TAccumulatedPlasticStrain);

    double sensitivity;


    // Elastic step ... no updates required
    if (f <= -DBL_EPSILON * E) {
        //if (f <= 1.0e-8) {
        // Set trial tangent
        TTangent = E;
        sensitivity = t_stressSensitivity;
    }

    // ------- Plastic corrector ... perform return mapping algorithm ---
    else {
      deltaLambda = (fabs(xsi)-CSigmaY)/(E+Hkin+Hiso);
      
      // Find sign of xsi
      int sign = (xsi < 0) ? -1 : 1;
      TPlasticStrain = CPlasticStrain +deltaLambda*sign;

      TBackStress = CBackStress + Hkin*deltaLambda*sign;

      TAccumulatedPlasticStrain = CAccumulatedPlasticStrain + deltaLambda;

      t_stress = E * (t_strain - TPlasticStrain);
    
      TTangent = E*(Hkin+Hiso) / (E+Hkin+Hiso);

      double CSigmaYSensitivity = SigmaYSensitivity + HisoSensitivity*CAccumulatedPlasticStrain
           + Hiso*CAccumulatedPlasticStrainSensitivity;

      double deltaLambdaSensitivity = ((t_stressSensitivity-CBackStressSensitivity)*sign-CSigmaYSensitivity)/(E+Hkin+Hiso)
          - (ESensitivity+HkinSensitivity+HisoSensitivity)*((E * (t_strain-CPlasticStrain)-CBackStress)*sign-CSigmaY)/pow((E+Hkin+Hiso),2.0);

      double TPlasticStrainSensitivity = CPlasticStrainSensitivity + deltaLambdaSensitivity*sign;

      sensitivity = E * (t_strainSensitivity - TPlasticStrainSensitivity)+ESensitivity * (t_strain - TPlasticStrain);

   }
    return sensitivity;
}



double
UniaxialJ2Plasticity::getInitialTangentSensitivity(int gradIndex)
{
    // For now, assume that this is only called for initial stiffness 
    if (parameterID == 2) {
        return 1.0; 
    }
    else {
        return 0.0;
    }
}


int
UniaxialJ2Plasticity::commitSensitivity(double t_strainSensitivity, int gradIndex, int numGrads)
{

        
    if (SHVs == 0) {
        SHVs = new Matrix(5,numGrads);
        SHVs->Zero();
    }



    // First set values depending on what is random
    double SigmaYSensitivity = 0.0;
    double ESensitivity = 0.0;
    double HkinSensitivity = 0.0;
    double HisoSensitivity = 0.0;

    if (parameterID == 1) {  // sigmaY
        SigmaYSensitivity = 1.0;
    }
    else if (parameterID == 2) {  // E
        ESensitivity = 1.0;
    }
    else if (parameterID == 3) {  // Hkin
        HkinSensitivity = 1.0;
    }
    else if (parameterID == 4) {  // Hiso
        HisoSensitivity = 1.0;
    }
    else {
        // Nothing random here, but may have to return something in any case
    }

    
    // Then pick up history variables for this gradient number

    double CPlasticStrainSensitivity = (*SHVs)(0,gradIndex);
    double CBackStressSensitivity     = (*SHVs)(1,gradIndex);
    double CAccumulatedPlasticStrainSensitivity = (*SHVs)(2,gradIndex);


    // ------ Elastic trial -------

    double deltaLambda = 0.0;
    TPlasticStrain = CPlasticStrain;
    TBackStress = CBackStress;
    TAccumulatedPlasticStrain = CAccumulatedPlasticStrain;
    t_stress = E * (t_strain-CPlasticStrain);
    double t_stressSensitivity = E*(t_strainSensitivity-CPlasticStrainSensitivity)+ESensitivity*(t_strain-CPlasticStrain);
    double CSigmaY = sigmaY+Hiso*CAccumulatedPlasticStrain;

    // Compute trial stress relative to committed back stress
    double xsi = t_stress - TBackStress;

    // Compute yield criterion
    double f = fabs(xsi) - (sigmaY + Hiso*TAccumulatedPlasticStrain);

    double sensitivity;


    // Elastic step ... no updates required
    if (f <= -DBL_EPSILON * E) {
        //if (f <= 1.0e-8) {
        // Set trial tangent
        TTangent = E;
        sensitivity = t_stressSensitivity;
    }

    // ------- Plastic step ... perform return mapping algorithm ---
    else {
      deltaLambda = (fabs(xsi)-CSigmaY)/(E+Hkin+Hiso);
      
      // Find sign of xsi
      int sign = (xsi < 0) ? -1 : 1;
      TPlasticStrain = CPlasticStrain +deltaLambda*sign;

      TBackStress = CBackStress + Hkin*deltaLambda*sign;

      TAccumulatedPlasticStrain = CAccumulatedPlasticStrain + deltaLambda;

      t_stress = E * (t_strain - TPlasticStrain);
    
      TTangent = E*(Hkin+Hiso) / (E+Hkin+Hiso);

      double CSigmaYSensitivity = SigmaYSensitivity + HisoSensitivity*CAccumulatedPlasticStrain
           + Hiso*CAccumulatedPlasticStrainSensitivity;

      double deltaLambdaSensitivity = ((t_stressSensitivity-CBackStressSensitivity)*sign-CSigmaYSensitivity)/(E+Hkin+Hiso)
          - (ESensitivity+HkinSensitivity+HisoSensitivity)*((E * (t_strain-CPlasticStrain)-CBackStress)*sign-CSigmaY)/pow((E+Hkin+Hiso),2.0);

      double TPlasticStrainSensitivity = CPlasticStrainSensitivity + deltaLambdaSensitivity*sign;

      sensitivity = E * (t_strainSensitivity - TPlasticStrainSensitivity)+ESensitivity * (t_strain - TPlasticStrain);

      double TAccumulatedPlasticStrainSensitivity = CAccumulatedPlasticStrainSensitivity + deltaLambdaSensitivity;

      double TBackStressSensitivity = CBackStressSensitivity + HkinSensitivity*deltaLambda*sign + Hkin*deltaLambdaSensitivity*sign;


      (*SHVs)(0,gradIndex) = TPlasticStrainSensitivity;
      (*SHVs)(1,gradIndex) = TBackStressSensitivity;
      (*SHVs)(2,gradIndex) = TAccumulatedPlasticStrainSensitivity;
      (*SHVs)(3,gradIndex) = sensitivity;      // for recorder purpose
      (*SHVs)(4,gradIndex) = t_strainSensitivity;  // for recorder purpose
    }

      

    return 0;
}

