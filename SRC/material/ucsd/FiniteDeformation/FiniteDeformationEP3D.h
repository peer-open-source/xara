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
//# DATE:              July 2004
//# UPDATE HISTORY:
//#
//===============================================================================
#pragma once

#include "BJMaterial.h"

#include "FDEPState.h"
#include "FiniteDeformationElastic3D.h"

#include <fdYield.h>
#include <fdFlow.h>
#include <fdEvolution_S.h>
#include <fdEvolution_T.h>

#include <FEM_ObjectBroker.h>
#include <Channel.h>
#include <OPS_Globals.h>
#include <ConsoleErrorHandler.h>

#include <stresst.h>
#include <straint.h>
#include <BJmatrix.h>
#include <BJvector.h>

class FiniteDeformationEP3D : public BJMaterial
{
public:
  // Constructor 00
  FiniteDeformationEP3D( );
  // Constructor 01
  FiniteDeformationEP3D(int tag,
                        BJMaterial *fde3d_in,
                        fdYield *fdy_in,
                        fdFlow *fdf_in,
                        fdEvolution_S *fdEvolutionS_in,
                        fdEvolution_T *fdEvolutionT_in);
  // Constructor 02
  FiniteDeformationEP3D(int tag,
                        BJMaterial *fde3d_in,
                        fdYield *fdy_in,
                        fdFlow *fdf_in,
                        fdEvolution_S *fdEvolutionS_in);
  // Constructor 03
  FiniteDeformationEP3D(int tag,
                        BJMaterial *fde3d_in,
                        fdYield *fdy_in,
                        fdFlow *fdf_in,
                        fdEvolution_T *fdEvolutionT_in);
  // Constructor 04
  FiniteDeformationEP3D(int tag,
                        BJMaterial *fde3d_in,
                        fdYield *fdy_in,
                        fdFlow *fdf_in);

  virtual ~FiniteDeformationEP3D( );
  const char *getClassType() const {return "FiniteDeformationEP3D";}
  
  double getRho();  

  int setTrialF(const straintensor &f);
  int setTrialFIncr(const straintensor &df);

  const Tensor& getTangentTensor() ;

  const  straintensor& getStrainTensor() ;  // Default Green Strain
  const  stresstensor& getStressTensor() ;  // Default 2nd Piola Kirchhoff Stress
  const  straintensor& getF();
  const  straintensor& getFp();

  int commitState(void) ;
  int revertToLastCommit(void) ;
  int revertToStart(void) ;

  BJMaterial *getCopyBJ();
  BJMaterial *getCopyBJ(const char *type);

  const char *getType() const;

  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

  void Print(OPS_Stream &s, int flag);

  const  stresstensor& getCauchyStressTensor();
  const  stresstensor& getPK1StressTensor() ;

private:

  BJMaterial *getFDE3D() const;
  fdYield *getFDY() const;
  fdFlow *getFDF() const;
  fdEvolution_S *getFDEvolutionS() const;
  fdEvolution_T *getFDEvolutionT() const;  
  FDEPState *getFDEPState() const; 

  int ImplicitAlgorithm();
  int SemiImplicitAlgorithm();

private:
  BJMaterial *fde3d;
  fdYield *fdy;
  fdFlow *fdf;
  fdEvolution_S *fdEvolutionS;
  fdEvolution_T *fdEvolutionT;

  // material input
  straintensor F;

  // material response
  straintensor iniGreen;
  stresstensor iniPK2;
  tensor iniTangent;
  
  stresstensor B_PK2;
  straintensor Fe;
  //stresstensor cauchystress;
  
  FDEPState *fdeps;
  
  static stresstensor static_stress; //Only for reference return
  static straintensor static_strain; //Only for reference return
};
