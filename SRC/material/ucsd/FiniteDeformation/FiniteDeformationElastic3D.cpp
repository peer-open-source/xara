//===============================================================================
//# COPYRIGHT (C): Woody's license (by BJ):
//                 ``This    source  code is Copyrighted in
//                 U.S.,  for  an  indefinite  period,  and anybody
//                 caught  using it without our permission, will be
//                 mighty good friends of ourn, cause we don't give
//                 a  darn.  Hack it. Compile it. Debug it. Run it.
//                 Yodel  it.  Enjoy it. We wrote it, that's all we
//                 wanted to do.''
//
//# PROJECT:           Object Oriented Finite Element Program
//# PURPOSE:           Finite Deformation Hyper-Elastic classes
//# CLASS:
//#
//# VERSION:           0.6_(1803398874989) (golden section)
//# LANGUAGE:          C++
//# TARGET OS:         all...
//# DESIGN:            Zhao Cheng, Boris Jeremic (jeremic@ucdavis.edu)
//# PROGRAMMER(S):     Zhao Cheng, Boris Jeremic
//#
//#
//# DATE:              19AUg2003
//# UPDATE HISTORY:    Sept 2003
//#                    May28, 2004
//#
//===============================================================================

#include <FiniteDeformationElastic3D.h>

//-----------------------------------------------------------------------------------------------------------------------------------------------
FiniteDeformationElastic3D::FiniteDeformationElastic3D(int tag,
                                                       int classTag,
                                                       double rho_in= 0.0)
:BJMaterial(tag, classTag), rho(rho_in)
{

}

//------------------------------------------------------------------------------------------------------------------------------------------------
FiniteDeformationElastic3D::FiniteDeformationElastic3D()
:BJMaterial(0, 0), rho(0.0)
{

}

//------------------------------------------------------------------------------------------------------------------------------------------------
FiniteDeformationElastic3D::~FiniteDeformationElastic3D()
{

}

//-------------------------------------------------------------------------------------------------------------------------------------------------
double FiniteDeformationElastic3D::getRho()
{
   return rho;
}

//--------------------------------------------------------------------------------------------------------------------------------------------------
int
FiniteDeformationElastic3D::setTrialF(const straintensor &f)
{
   opserr << "FiniteDeformationElastic3D-- subclass responsibility\n";
   exit (-1);
   return 0;
}

//---------------------------------------------------------------------------------------------------------------------------------------------------
int FiniteDeformationElastic3D::setTrialFIncr(const straintensor &df)
{
   opserr << "FiniteDeformationElastic3D-- subclass responsibility\n";
   exit (-1);
   return 0;
}

//---------------------------------------------------------------------------------------------------------------------------------------------------
int FiniteDeformationElastic3D::setTrialC(const straintensor &c)
{
   opserr << "FiniteDeformationElastic3D-- subclass responsibility\n";
   exit (-1);
   return 0;
}

//---------------------------------------------------------------------------------------------------------------------------------------------------
int FiniteDeformationElastic3D::setTrialCIncr(const straintensor &dc)
{
   opserr << "FiniteDeformationElastic3D-- subclass responsibility\n";
   exit (-1);
   return 0;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------
const straintensor& FiniteDeformationElastic3D::getF()
{
   opserr << "FiniteDeformationElastic3D-- subclass responsibility\n";
   exit (-1);
   // Just to make it compile
   straintensor *ret = new straintensor;
   return *ret;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------
const straintensor& FiniteDeformationElastic3D::getC()
{
   opserr << "FiniteDeformationElastic3D-- subclass responsibility\n";
   exit (-1);
   // Just to make it compile
   straintensor *ret = new straintensor;
   return *ret;

}

//------------------------------------------------------------------------------------------------------------------------------------------------------------------------
const Tensor& FiniteDeformationElastic3D::getTangentTensor()
{
   exit (-1);
   // Just to make it compile
   Tensor *ret = new Tensor;
   return *ret;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
const Tensor
&FiniteDeformationElastic3D::getInitialTangentTensor()
{
   opserr << "FiniteDeformationElastic3D-- subclass responsibility\n";
   exit (-1);
   // Just to make it compile
   Tensor *ret = new Tensor;
   return *ret;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
const straintensor& 
FiniteDeformationElastic3D::getStrainTensor()
{
   opserr << "FiniteDeformationElastic3D-- subclass responsibility\n";
   exit (-1);
   // Just to make it compile
   straintensor *ret = new straintensor;
   return *ret;
}
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
const stresstensor& FiniteDeformationElastic3D::getStressTensor()
{
   opserr << "FiniteDeformationElastic3D-- subclass responsibility\n";
   exit (-1);
   stresstensor *ret = new stresstensor; 
   return *ret;
}
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
const stresstensor& FiniteDeformationElastic3D::getPK1StressTensor()
{
   opserr << "FiniteDeformationElastic3D-- subclass responsibility\n";
   exit (-1);
   stresstensor *ret = new stresstensor; 
   return *ret;
}
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
const stresstensor& FiniteDeformationElastic3D::getCauchyStressTensor(void)
{
   opserr << "FiniteDeformationElastic3D-- subclass responsibility\n";
   exit (-1);
   stresstensor *ret = new stresstensor; 
   return *ret;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int FiniteDeformationElastic3D::commitState ()
{
   opserr << "FiniteDeformationElastic3D-- subclass responsibility\n";
   exit (-1);
   return -1;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int FiniteDeformationElastic3D::revertToLastCommit ()
{
   opserr << "FiniteDeformationElastic3D-- subclass responsibility\n";
   exit (-1);
   return -1;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int FiniteDeformationElastic3D::revertToStart ()
{
   opserr << "FiniteDeformationElastic3D-- subclass responsibility\n";
   exit (-1);
   return -1;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
BJMaterial * FiniteDeformationElastic3D::getCopyBJ()
{
   opserr << "FiniteDeformationElastic3D-- subclass responsibility\n";
   exit (-1);
   return 0;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
BJMaterial * FiniteDeformationElastic3D::getCopyBJ (const char *type)
{
   opserr << "FiniteDeformationElastic3D-- subclass responsibility\n";
   exit (-1);
   return 0;
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
const char* FiniteDeformationElastic3D::getType () const
{
   opserr << "FiniteDeformationElastic3D-- subclass responsibility\n";
   exit (-1);
   return 0;
}
////--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//int FiniteDeformationElastic3D::getOrder (void) const
//{
//   opserr << "FiniteDeformationElastic3D-- subclass responsibility\n";
//   exit (-1);
//   return 0;
//}
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int FiniteDeformationElastic3D::sendSelf (int commitTag, Channel &theChannel)
{
   opserr << "FiniteDeformationElastic3D-- subclass responsibility\n";
   exit (-1);
   return 0;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int FiniteDeformationElastic3D::recvSelf (int commitTag,
                                          Channel &theChannel,
                                          FEM_ObjectBroker &theBroker)
{
   opserr << "FiniteDeformationElastic3D-- subclass responsibility\n";
   exit (-1);
   return 0;
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void FiniteDeformationElastic3D::Print (OPS_Stream &s, int flag)
{
   opserr << "FiniteDeformationElastic3D-- subclass responsibility\n";
   exit (-1);
   return;
}
