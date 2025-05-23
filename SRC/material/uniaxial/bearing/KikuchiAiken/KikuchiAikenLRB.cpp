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
// Written: Masaru Kikuchi
// Created: August 10, 2012
//
// Kikuchi&Aiken model for lead rubber bearing
//
// Description: This file contains the function to parse the TCL input
// uniaxialMaterial KikuchiAikenLRB matTag? type? ar? hr? gr? ap? tp? alph? beta? <-T temp? > <-coKQ rk? rq?> <-coMSS rs? rf?>


#include <string.h>
#include <ID.h>
#include <Logging.h>
#include <Vector.h>
#include "KikuchiAikenLRB.h"
#include <Channel.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <float.h>


void
KikuchiAikenLRB::setType(int Tp) {
  //parameter function for each rubber
  switch (Tp) {
  case 1: // LRB
    trgStrain = 0.05;
    lmtStrain = 4.10;
    calcN   = KikuchiAikenLRB::calcNType1;
    calcP   = KikuchiAikenLRB::calcPType1;
    calcA   = KikuchiAikenLRB::calcAType1;
    calcB   = KikuchiAikenLRB::calcBType1;
    calcC   = KikuchiAikenLRB::calcCType1;
    calcCQd = KikuchiAikenLRB::calcCQdType1;
    calcCKd = KikuchiAikenLRB::calcCKdType1;
    calcCHeq= KikuchiAikenLRB::calcCHeqType1;
    break;
  }
}

KikuchiAikenLRB::KikuchiAikenLRB(int tag, int type, double ar, double hr, double gr, double ap, double tp, 
		  double alph, double beta, double temp, double rk, double rq, double rs, double rf)
  :UniaxialMaterial(tag,MAT_TAG_KikuchiAikenLRB),Type(type),Ar(ar),Hr(hr),Gr(gr),Ap(ap),Tp(tp),
   Alph(alph),Beta(beta),Temp(temp),Rk(rk),Rq(rq),Rs(rs),Rf(rf)
{

  this->setType(Type);

  //initialize
  qd100 = Tp*Ap * exp(-0.00879*(Temp-15.0)) * Rq;
  kd100 = (Gr*Ar/Hr + Alph*Ap/Hr) * exp(-0.00271*(Temp-15.0)) * Rk;
  ku100 = Beta * kd100;

  qd = (this->calcCQd)(trgStrain)*qd100;
  kd = (this->calcCKd)(trgStrain)*kd100;
  ku = (this->calcCKd)(trgStrain)*ku100;
  initialStiff = compKeq(fabs(trgStrain*Hr),qd,kd);

  //
  numIdx = 500;
  revXBgn   = new double [numIdx];
  revQ2Bgn  = new double [numIdx];
  revXEnd   = new double [numIdx];
  revQ2End  = new double [numIdx];
  revB      = new double [numIdx];
  revAlpha  = new double [numIdx];

  for (int i = 0; i < numIdx; i++) {
    revXBgn[i] = 0.0;
    revQ2Bgn[i] = 0.0;
    revXEnd[i] = 0.0;
    revQ2End[i] = 0.0;
    revB[i] = 0.0;
    revAlpha[i] = 0.0;
  }
  
  trialDeform  = 0.0;
  trialForce   = 0.0;
  trialStiff   = initialStiff;
  trialStrain  = 0.0;
  trialIfElastic = true;
  trialQ1 = 0.0;
  trialQ2 = 0.0;
  trialMaxStrain = 0.0;
  trialDDeform = 0.0;
  trialDDeformLastSign = 0;
  trialIdxRev=0;

  commitDeform  = 0.0;
  commitForce   = 0.0;
  commitStiff   = initialStiff;
  commitStrain  = 0.0;
  commitIfElastic = true;
  commitQ1 = 0.0;
  commitQ2 = 0.0;
  commitMaxStrain = 0.0;
  commitDDeform = 0.0;
  commitDDeformLastSign = 0;
  commitIdxRev=0;

  revB[0] = 0.0;


}

KikuchiAikenLRB::KikuchiAikenLRB()
  :UniaxialMaterial(0,MAT_TAG_KikuchiAikenLRB)
{
  revXBgn   = 0;
  revQ2Bgn  = 0;
  revXEnd   = 0;
  revQ2End  = 0;
  revB      = 0;
  revAlpha  = 0;
  
  //

  trialDeform  = 0.0;
  trialForce   = 0.0;
  trialStiff   = 0.0;
  commitDeform = 0.0;
  commitForce  = 0.0;
  commitStiff  = 0.0;


}

KikuchiAikenLRB::~KikuchiAikenLRB()
{

  if (revXBgn != 0)
    delete [] revXBgn;

  if (revQ2Bgn != 0)
    delete [] revQ2Bgn;

  if (revXEnd != 0)
    delete [] revXEnd;

  if (revQ2End != 0)
    delete [] revQ2End;

  if (revB != 0)
    delete [] revB;

  if (revAlpha != 0)
    delete [] revAlpha;


}

int 
KikuchiAikenLRB::setTrialStrain(double strain, double strainRate)
{

  // (input)                             (output)
  // deformation -> strain -> stress  -> force
  //                strain -> modulus -> spring constant


  //(input: deformation->strain)
  trialDeform = strain;
  trialDeform = trialDeform * (Rs/Rf); //for MSS model
  trialStrain = trialDeform/Hr;


  //incremental deformation
  trialDDeform = trialDeform - commitDeform;

  //if incremental deformation is zero
  if ( fabs(trialDDeform) < DBL_EPSILON ) { // if ( trialDDeform == 0.0 )
    trialForce   = commitForce;
    trialStiff   = commitStiff;
    return 0;
  }

  //sign of incremental deformation (used in next caluculation step)
  if (trialDDeform > 0) {
    trialDDeformLastSign = +1;
  } else if (trialDDeform < 0) {
    trialDDeformLastSign = -1;
  } else {
    trialDDeformLastSign = commitDDeformLastSign;
  }

  //copy reversal point index
  trialIdxRev = commitIdxRev;


  //application limit strain
  if ( fabs(trialStrain) > lmtStrain ) {
    opserr << "uniaxialMaterial KikuchiAikenLRB: \n";
    opserr << "   Response value exceeded limited strain.\n";
    // return -1;
    // warning, extrapolation
  }

  //elastic limit strain
  if ( fabs(trialStrain) > trgStrain ) {
      trialIfElastic = false;
  }

  //maximum strain
  if ( fabs(trialStrain) > commitMaxStrain ) {
    trialMaxStrain = fabs(trialStrain);
  }

  //parameters for Kikuchi&Aiken model

  if ( trialIfElastic || fabs(trialStrain) == trialMaxStrain ) {

    //if elactic, tmpstrain is max(abs(trialStrain),trgStrain)
    tmpStrain = (fabs(trialStrain)>trgStrain) ? fabs(trialStrain) : trgStrain ; //max(abs(trialStrain),trgStrain)
    tmpDeform = tmpStrain * Hr;

    qd = (this->calcCQd)(tmpStrain)*qd100;
    kd = (this->calcCKd)(tmpStrain)*kd100;
    ku = (this->calcCKd)(tmpStrain)*ku100;

    keq = compKeq(tmpDeform,qd,kd);
    heq = (this->calcCHeq)(tmpStrain) * compHeq(tmpDeform,qd,kd,ku);
    u = qd / (keq*tmpDeform);

    xm  = fabs(trialDeform);
    fm  = keq*xm;

    n = (this->calcN)(fabs(trialStrain));
    p = (this->calcP)(fabs(trialStrain));
    c = (this->calcC)(fabs(trialStrain));
    a = (this->calcA)(fabs(trialStrain),heq,u);

  }  // if ( trialIfElastic || abs(trialStrain) == trialMaxStrain )  ... end


  // normalized deformation
  x = (xm>0.0) ? trialDeform/xm : 0.0;

  // reversal points
  if (!trialIfElastic) {

    // unload or reverse
    if (trialDDeform*commitDDeformLastSign < 0) {


      if ( trialIdxRev == 0 ) { // unload

	trialIdxRev = 1;

	revXBgn[1]  = commitDeform/xm;
	revQ2Bgn[1] = commitQ2;
	
	//b
	double tqd = (this->calcCQd)(fabs(commitStrain))*qd100;
	double tkd = (this->calcCKd)(fabs(commitStrain))*kd100;
	double tku = (this->calcCKd)(fabs(commitStrain))*ku100;
	double tkeq = compKeq(fabs(commitDeform),tqd,tkd);
	double theq = (this->calcCHeq)(fabs(commitStrain)) * compHeq(fabs(commitDeform),tqd,tkd,tku);
	double tu = tqd/(tkeq*fabs(commitDeform));
	b = (this->calcB)(fabs(commitStrain),a,c,theq,tu);

	revB[1] = b;
	revAlpha[1] = 1.0;


      } else { //reverse

	trialIdxRev++;


	if (trialIdxRev >= numIdx ) { //memorize reversal points of hysteretic curve
	  int newIdx;
	  double *newXBgn;
	  double *newQ2Bgn;
	  double *newXEnd;
	  double *newQ2End;
	  double *newB;
	  double *newAlpha;
	  
	  newIdx = numIdx + 500;
	  newXBgn  = new double [newIdx];
	  newQ2Bgn = new double [newIdx];
	  newXEnd  = new double [newIdx];
	  newQ2End = new double [newIdx];
	  newB     = new double [newIdx];
	  newAlpha = new double [newIdx];
	  
	  for(int i = 0; i<numIdx; i++){
	    newXBgn[i]  = revXBgn[i];
	    newQ2Bgn[i] = revQ2Bgn[i];
	    newXEnd[i]  = revXEnd[i];
	    newQ2End[i] = revQ2End[i];
	    newB[i]     = revB[i];
	    newAlpha[i] = revAlpha[i];
	  }
	  
	  numIdx = newIdx;

	  delete [] revXBgn;
	  delete [] revQ2Bgn;
	  delete [] revXEnd;
	  delete [] revQ2End;
	  delete [] revB;
	  delete [] revAlpha;

	  revXBgn  = newXBgn;
	  revQ2Bgn = newQ2Bgn;
	  revXEnd  = newXEnd;
	  revQ2End = newQ2End;
	  revB     = newB;
	  revAlpha = newAlpha;
	}


	revXEnd[trialIdxRev]  = revXBgn[trialIdxRev-1];
	revQ2End[trialIdxRev] = revQ2Bgn[trialIdxRev-1];
	revXBgn[trialIdxRev]  = commitDeform/xm;
	revQ2Bgn[trialIdxRev] = commitQ2;
	
	//b
	if ( revB[trialIdxRev-1] == 0 ) { // after reversal point with b=0
	  b = 0;
	} else if (revXEnd[trialIdxRev]*revXBgn[trialIdxRev] > 0) { //consecutive reverse at same sign of "x"
	  b = 0;
	} else { //reverse at different sign of "x"
	  double tqd = (this->calcCQd)(fabs(commitStrain))*qd100;
	  double tkd = (this->calcCKd)(fabs(commitStrain))*kd100;
	  double tku = (this->calcCKd)(fabs(commitStrain))*ku100;
	  double tkeq = compKeq(fabs(commitDeform),tqd,tkd);
	  double theq = (this->calcCHeq)(fabs(commitStrain)) * compHeq(fabs(commitDeform),tqd,tkd,tku);
	  double tu = tqd/(tkeq*fabs(commitDeform)); 
	  b = (this->calcB)(fabs(commitStrain),a,c,theq,tu);
	}

	
	//alpha
	if (trialDDeform > 0) {
	  alpha = (this->compAlpha)(a,revB[trialIdxRev-1],b,c,revXEnd[trialIdxRev],revXBgn[trialIdxRev],revAlpha[trialIdxRev-1]);
	} else {
	  alpha = (this->compAlpha)(a,revB[trialIdxRev-1],b,c,-revXEnd[trialIdxRev],-revXBgn[trialIdxRev],revAlpha[trialIdxRev-1]);
	}

	revB[trialIdxRev] = b;
	revAlpha[trialIdxRev] = alpha;

      }

    }


    //forget reversal points
    if  (fabs(trialStrain) == trialMaxStrain) { //if reach skeleton curve
      trialIdxRev = 0;
    } else if (trialIdxRev >= 2 ) { //if pass through reversal point
      while (trialIdxRev >= 2 && (x-revXBgn[trialIdxRev])*(x-revXEnd[trialIdxRev]) > 0) {
	trialIdxRev--;
      }
    }

  }  //   if (!trialIfElastic)  ... end

  // calculate stress

  // Q1 component
  if (trialStrain > 0) {
    trialQ1 = (this->compQ1)(u,n,p,fm,x);
    q1Stf   = (this->compQ1Derivertive)(u,n,p,keq,x);
  } else {
    trialQ1 = (this->compQ1)(u,n,p,-fm,-x);
    q1Stf   = (this->compQ1Derivertive)(u,n,p,keq,-x);
  }

  //Q2 component
  if (trialIdxRev == 0) {// skeleton curve

    if (trialDDeform > 0) {
      trialQ2 = (this->compQ2Unload)(u,a,revB[trialIdxRev],c,-fm,-x);
      q2Stf   = (this->compQ2UnloadDerivertive)(u,a,revB[trialIdxRev],c,keq,x);
    } else {
      trialQ2 = (this->compQ2Unload)(u,a,revB[trialIdxRev],c,fm,x);
      q2Stf   = (this->compQ2UnloadDerivertive)(u,a,revB[trialIdxRev],c,keq,-x);
    }

  } else if (trialIdxRev == 1) { //unload

    if (trialDDeform > 0) {
      trialQ2 = (this->compQ2Unload)(u,a,revB[trialIdxRev],c,fm,x);
      q2Stf   = (this->compQ2UnloadDerivertive)(u,a,revB[trialIdxRev],c,keq,x);
    } else {
      trialQ2 = (this->compQ2Unload)(u,a,revB[trialIdxRev],c,-fm,-x);
      q2Stf   = (this->compQ2UnloadDerivertive)(u,a,revB[trialIdxRev],c,keq,-x);
    }

  } else { //reverse

    if (trialDDeform > 0) {
      trialQ2 = (this->compQ2Masing)(u,a,revB[trialIdxRev],c,fm,x,revXBgn[trialIdxRev],revQ2Bgn[trialIdxRev],revAlpha[trialIdxRev]);
      q2Stf   = (this->compQ2MasingDerivertive)(u,a,revB[trialIdxRev],c,keq,x,revXBgn[trialIdxRev],revAlpha[trialIdxRev]);
    } else {
      trialQ2 = (this->compQ2Masing)(u,a,revB[trialIdxRev],c,-fm,-x,-revXBgn[trialIdxRev],revQ2Bgn[trialIdxRev],revAlpha[trialIdxRev]);
      q2Stf   = (this->compQ2MasingDerivertive)(u,a,revB[trialIdxRev],c,keq,-x,-revXBgn[trialIdxRev],revAlpha[trialIdxRev]);
    }

  }


  //force
  trialForce  = trialQ1 + trialQ2;

  //stiffness
  //trialStiff = (trialForce-commitForce)/trialDDeform;
  if (trialIfElastic) {
    trialStiff = initialStiff;
  } else {
    trialStiff = q1Stf + q2Stf;
  }

  trialForce = trialForce * Rf; //for MSS model
  trialStiff = trialStiff * Rs; //for MSS model

  return 0;

}



double 
KikuchiAikenLRB::getStress(void)
{
  return trialForce;
}

double 
KikuchiAikenLRB::getTangent(void)
{
  return trialStiff;
}

double 
KikuchiAikenLRB::getInitialTangent(void)
{
  return initialStiff;
}

double 
KikuchiAikenLRB::getStrain(void)
{
  return trialDeform;
}

int 
KikuchiAikenLRB::commitState(void)
{
  commitDeform  = trialDeform;
  commitForce   = trialForce;
  commitStiff   = trialStiff;
  commitStrain  = trialStrain;
  commitIfElastic = trialIfElastic;
  commitQ1 = trialQ1;
  commitQ2 = trialQ2;
  commitMaxStrain = trialMaxStrain;
  commitDDeform   = trialDDeform;
  commitDDeformLastSign = trialDDeformLastSign;
  commitIdxRev = trialIdxRev;

  return 0;
}

int 
KikuchiAikenLRB::revertToLastCommit(void)
{

  trialDeform  = commitDeform;
  trialForce   = commitForce;
  trialStiff   = commitStiff;
  trialStrain  = commitStrain;
  trialIfElastic = commitIfElastic;
  trialQ1 = commitQ1;
  trialQ2 = commitQ2;
  trialMaxStrain = commitMaxStrain;
  trialDDeform = commitDDeform;
  trialDDeformLastSign = commitDDeformLastSign;
  trialIdxRev = commitIdxRev;

  return 0;
}

int 
KikuchiAikenLRB::revertToStart(void)
{

  trialDeform  = 0.0;
  trialForce   = 0.0;
  trialStiff   = initialStiff;
  trialStrain  = 0.0;
  trialIfElastic = true;
  trialQ1 = 0.0;
  trialQ2 = 0.0;
  trialMaxStrain = 0.0;
  trialDDeform = 0.0;
  trialDDeformLastSign = 0;
  trialIdxRev=0;

  commitDeform  = 0.0;
  commitForce   = 0.0;
  commitStiff   = initialStiff;
  commitStrain  = 0.0;
  commitIfElastic = true;
  commitQ1 = 0.0;
  commitQ2 = 0.0;
  commitMaxStrain = 0.0;
  commitDDeform = 0.0;
  commitDDeformLastSign = 0;
  commitIdxRev=0;

  revB[0] = 0.0;

  return 0;
}

UniaxialMaterial *
KikuchiAikenLRB::getCopy(void)
{

  KikuchiAikenLRB *theCopy = new KikuchiAikenLRB(this->getTag(), Type, Ar, Hr, Gr, Ap, Tp, 
			      Alph, Beta, Temp, Rk, Rq, Rs, Rf );

  // Copy temporary variables
  theCopy->tmpStrain = tmpStrain;
  theCopy->keq = keq;
  theCopy->heq = heq;
  theCopy->u = u;
  theCopy->n = n;
  theCopy->a = a;
  theCopy->b = b;
  theCopy->c = c;
  theCopy->xm = xm;
  theCopy->fm = fm;
  theCopy->x = x;
  theCopy->alpha = alpha;

  // Copy trial variables
  theCopy->trialDeform = trialDeform;
  theCopy->trialForce = trialForce;
  theCopy->trialStiff = trialStiff;
  theCopy->trialStrain = trialStrain;
  theCopy->trialIfElastic = trialIfElastic;
  theCopy->trialQ1 = trialQ1;
  theCopy->trialQ2 = trialQ2;
  theCopy->trialMaxStrain = trialMaxStrain;
  theCopy->trialDDeform = trialDDeform;
  theCopy->trialDDeformLastSign = trialDDeformLastSign;
  theCopy->trialIdxRev = trialIdxRev;

  // Copy commit variables
  theCopy->commitDeform = commitDeform;
  theCopy->commitForce = commitForce;
  theCopy->commitStiff = commitStiff;
  theCopy->commitStrain = commitStrain;
  theCopy->commitIfElastic = commitIfElastic;
  theCopy->commitQ1 = commitQ1;
  theCopy->commitQ2 = commitQ2;
  theCopy->commitMaxStrain = commitMaxStrain;
  theCopy->commitDDeform = commitDDeform;
  theCopy->commitDDeformLastSign = commitDDeformLastSign;
  theCopy->commitIdxRev = commitIdxRev;

  return theCopy;
}

int 
KikuchiAikenLRB::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;

  ID idata(3);
  idata(0) = this->getTag();
  idata(1) = Type;
  idata(2) = numIdx;

  res = theChannel.sendID(this->getDbTag(), cTag, idata);
  if (res < 0) {
    opserr << "KikuchiAikenLRB::sendSelf - failed to send ID data" << endln;
    return -1;
  }  

  Vector data(18+3+11+6*numIdx);
  data(0) = Ar;
  data(1) = Hr;
  data(2) = Gr;
  data(3) = Ap;
  data(4) = Tp;
  data(5) = Alph;
  data(6) = Beta;
  data(7) = Temp;
  data(8) = Rk;
  data(9) = Rq;
  data(10) = Rs;
  data(11) = Rf;

  data(12) = qd100;
  data(13) = kd100;
  data(14) = ku100;
  data(15) = qd;
  data(16) = kd;
  data(17) = ku;

  data(18) = trgStrain;
  data(19) = lmtStrain;
  data(20) = initialStiff;

  data(21) = commitDeform;
  data(22) = commitForce;
  data(23) = commitStiff;
  data(24) = commitStrain;
  data(25) = commitIfElastic ? 1.0 : -1.0;
  data(26) = commitQ1;
  data(27) = commitQ2;
  data(28) = commitMaxStrain;
  data(29) = commitDDeform;
  data(30) = commitDDeformLastSign;
  data(31) = commitIdxRev;          

  for (int i = 0; i < numIdx; i++) {
    data(32 +            i) = revXBgn[i];
    data(32 +   numIdx + i) = revQ2Bgn[i];
    data(32 + 2*numIdx + i) = revXEnd[i];
    data(32 + 3*numIdx + i) = revQ2End[i];
    data(32 + 4*numIdx + i) = revB[i];
    data(32 + 5*numIdx + i) = revAlpha[i];                    
  }
  
  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) {
    opserr << "KikuchiAikenLRB::sendSelf - failed to send vector data" << endln;
    return -2;
  }
  
  return res;
}

int 
KikuchiAikenLRB::recvSelf(int cTag, Channel &theChannel, 
			      FEM_ObjectBroker &theBroker)
{
  int res = 0;

  ID idata(3);
  res = theChannel.recvID(this->getDbTag(), cTag, idata);
  if (res < 0) {
    opserr << "KikuchiAikenLRB::recvSelf - failed to receive ID data" << endln;
    return -1;
  }

  this->setTag(idata(0));
  Type = idata(1);
  numIdx = idata(2);

  this->setType(Type);
  
  Vector data(18+3+11+6*numIdx);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  if (res < 0) {
    opserr << "KikuchiAikenLRB::recvSelf - failed to receive vector data" << endln;
    return -2;
  }

  Ar = data(0);
  Hr = data(1);
  Gr = data(2);
  Ap = data(3);
  Tp = data(4);
  Alph = data(5);
  Beta = data(6);
  Temp = data(7);
  Rk = data(8);
  Rq = data(9);
  Rs = data(10);
  Rf = data(11);

  qd100 = data(12);
  kd100 = data(13);
  ku100 = data(14);
  qd = data(15);
  kd = data(16);
  ku = data(17);
  
  trgStrain = data(18);
  lmtStrain = data(19);
  initialStiff = data(20);

  commitDeform = data(21);
  commitForce = data(22);
  commitStiff = data(23);
  commitStrain = data(24);
  commitIfElastic = data(25) > 0.0 ? true : false;
  commitQ1 = data(26);
  commitQ2 = data(27);
  commitMaxStrain = data(28);
  commitDDeform = data(29);
  commitDDeformLastSign = int(data(30));
  commitIdxRev = int(data(31));          

  if (revXBgn != 0) delete [] revXBgn;
  revXBgn = new double[numIdx];
  if (revQ2Bgn != 0) delete [] revQ2Bgn;
  revQ2Bgn = new double[numIdx];
  if (revXEnd != 0) delete [] revXEnd;
  revXEnd = new double[numIdx];
  if (revQ2End != 0) delete [] revQ2End;
  revQ2End = new double[numIdx];
  if (revB != 0) delete [] revB;
  revB = new double[numIdx];
  if (revAlpha != 0) delete [] revAlpha;
  revAlpha = new double[numIdx];        
  
  for (int i = 0; i < numIdx; i++) {
    revXBgn[i]  = data(32 +            i);
    revQ2Bgn[i] = data(32 +   numIdx + i);
    revXEnd[i]  = data(32 + 2*numIdx + i);
    revQ2End[i] = data(32 + 3*numIdx + i);
    revB[i]     = data(32 + 4*numIdx + i);
    revAlpha[i] = data(32 + 5*numIdx + i);
  }

  this->revertToLastCommit();
  
  return res;
}

void 
KikuchiAikenLRB::Print(OPS_Stream &s, int flag)
{
  s << "KikuchiAikenLRB : " << this->getTag() << endln;
  s << "  Type: " << Type << endln;
  s << "  Ar: " << Ar << endln;
  s << "  Hr: " << Hr << endln;
  s << "  Gr: " << Gr << endln;
  s << "  Ap: " << Ap << endln;
  s << "  Tp: " << Tp << endln;
  s << "  Alph: " << Alph << endln;
  s << "  Beta: " << Beta << endln;
  s << "  Temp: " << Temp << endln;
  s << "  Rk: " << Rk << endln;
  s << "  Rq: " << Rq << endln;
  s << "  Rs: " << Rs << endln;
  s << "  Rf: " << Rf << endln;

  return;
}



//------------------------------------------------------------
//formulae of Kikuchi&Aiken model
//------------------------------------------------------------

double
KikuchiAikenLRB::compQ1(double u, double n, double p, double fm, double x)
{
//  new F1 model
  return (1-u)*fm*( (1-p)*x + p*pow(x,n) );
}

double
KikuchiAikenLRB::compQ2Unload(double u, double a, double b, double c, double fm, double x)
{
  return u*fm*(1-2*exp(-a*(1+x))+b*(1+x)*exp(-c*(1+x)));
}

double
KikuchiAikenLRB::compQ2Masing(double u, double a, double b, double c, double fm, double x1, double x2, double q2i, double alpha)
{
  return q2i + alpha*u*fm*(2-2*exp(-a*(x1-x2))+b*(x1-x2)*exp(-c*(x1-x2)));
}


double
KikuchiAikenLRB::compAlpha(double a, double b1, double b2, double c, double x1, double x2, double alpha0)
{
  return alpha0*(2-2*exp(-a*(x1-x2))+b1*(x1-x2)*exp(-c*(x1-x2)))/(2-2*exp(-a*(x1-x2))+b2*(x1-x2)*exp(-c*(x1-x2)));
}


double
KikuchiAikenLRB::compQ1Derivertive(double u, double n, double p, double keq, double x)
{
//  new F1 model
    return (1-u)*keq*( (1-p) + p*n*pow(x,n-1) );
}

double
KikuchiAikenLRB::compQ2UnloadDerivertive(double u, double a, double b, double c, double keq, double x)
{
  return u*keq*(2*a*exp(-a*(1+x))+b*exp(-c*(1+x))-b*c*(1+x)*exp(-c*(1+x)));
}


double
KikuchiAikenLRB::compQ2MasingDerivertive(double u, double a, double b, double c, double keq, double x1, double x2,double alpha)
{
  return alpha*u*keq*(2*a*exp(-a*(x1-x2))+b*exp(-c*(x1-x2))-b*c*(x1-x2)*exp(-c*(x1-x2)));
}

//bisection method to calculate parameter "a"
double
KikuchiAikenLRB::compABisection(double heq, double u, double min, double max, double tol, double lim)
{
  double aMin, aMax, aTmp ;
  double RHS,LHS ;

  RHS = (2*u-M_PI*heq)/(2*u);

  aMin = min;
  aMax = max;

  while (true) {
    aTmp = (aMin+aMax)/2;
    LHS = (1-exp(-2*aTmp))/(aTmp);
    if (fabs((LHS-RHS)/RHS) < tol) {
      break;
    } else if (LHS < RHS) {
      aMax = aTmp;
    } else {
      aMin = aTmp;
    }
  }

  return (aTmp < lim) ? aTmp : lim ; //min(aTmp,lim)

}

//Keq
double
KikuchiAikenLRB::compKeq(double xm, double qd, double kd)
{
  return (qd + kd*xm)/xm;
}

//Heq
double
KikuchiAikenLRB::compHeq(double xm, double qd, double kd, double ku)
{
  double dy,qy,fm,heq;

    dy  = qd/(ku-kd);
    qy  = ku*dy;
    fm  = qd + kd*xm;
    heq = 2.0*qd/M_PI/fm*(1.0-qy/ku/xm);

  return heq ;

}



//------------------------------------------------------------
//parameters for each rubber
//------------------------------------------------------------

//--------------------------
// type1
//--------------------------

double KikuchiAikenLRB::calcNType1(double gm)
{return 8.0;}

double KikuchiAikenLRB::calcPType1(double gm)
{
  if      (gm<2.0)  {return 0.0;}
  else              {return -0.30226 + 0.15113*gm;}
}

double KikuchiAikenLRB::calcCType1(double gm)
{return 6.0;}

double KikuchiAikenLRB::calcCQdType1(double gm)
{
  if      (gm<0.1)  {return 2.036*pow(gm,0.410);}
  else if (gm<0.5)  {return 1.106*pow(gm,0.145);}
  else              {return 1.0;}
}

double KikuchiAikenLRB::calcCKdType1(double gm)
{
  if      (gm<0.25) {return 0.779*pow(gm,-0.43 );}
  else if (gm<1.0)  {return       pow(gm,-0.25 );}
  else if (gm<2.0)  {return       pow(gm,-0.12 );}
  else              {return 0.0482025*pow((gm-2.0),2.0)+0.92019;}
}

double KikuchiAikenLRB::calcCHeqType1(double gm)
{
  if      (gm<2.0)  {return 1.0;}
  else              {return 0.036375*pow((gm-2.0),2.0)+1.0;}
}

double KikuchiAikenLRB::calcAType1(double gm, double heq, double u)
{
  if      (gm<2.0)  {return compABisection(heq,u,0.0,50.0,1e-6,26.501472);}
  else              {return 26.501472;}
}

double KikuchiAikenLRB::calcBType1(double gm, double a, double c, double heq, double u)
{
  if      (gm<2.0)  {return 0.0;}
  else              {return c*c*(M_PI*heq/u-(2+2/a*(exp(-2*a)-1)));}
}


