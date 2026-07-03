/* ============================================================================
 * FEDEASLab - Release 5.1, July 2020
 * Matlab Finite Elements for Design, Evaluation and Analysis of Structures
 * Professor Filip C. Filippou (filippou@berkeley.edu)
 * Department of Civil and Environmental Engineering, UC Berkeley
 * Copyright(c) 1998-2021. The Regents of the University of California.
 * All Rights Reserved.
 * ============================================================================
 * Adapted from MATLAB source code by Claudio Perez                     07/2021
 */



#include "DegradingUniaxialWrapper.h"
#include "distributions.h"

/* `struct DmgEvol`, `struct DmgFrac` and `struct DmgResp` all look
 * like they'll be small enough for passing by value to
 * be efficient, but maybe this should change.
 */

struct DmgResp
dmglib_none(double x, struct DmgEvol dp, struct DmgFrac Frac)
{
  struct DmgResp resp = {0.0, 0.0};
  return resp;
}

struct DmgResp
dmglib_MBeta(double x, struct DmgEvol dp, struct DmgFrac Frac)
{
  double y, dydx;
  double a, b, psiF, psiU;

  if (dp.theParam[0] > 1) {
    a = dp.theParam[1];
    b = dp.theParam[0] * dp.theParam[1];
  } else {
    a = dp.theParam[1] / dp.theParam[0];
    b = dp.theParam[1];
  }

  // Fracture damage
  if (Frac.Activ) {
    psiF = Frac.psiF;
    psiU = Frac.psiU;
    y = 0.0;
    dydx = y;
    if (x <= psiF) {
      y = betacdf(x, a, b);
      dydx = betapdf(x, a, b);

    } else if (x <= psiU) {
      y = betacdf(psiF, a, b) +
          (x - psiF) / (psiU - psiF) * (1 - betacdf(psiF, a, b));
      dydx = (1 - betacdf(psiF, a, b)) / (psiU - psiF);

    } else {      // x > psiU
      y = 1.0;
      dydx = 0.0;
    }

  } else {
    // no fracture damage
    y = betacdf(x, a, b);
    dydx = betapdf(x, a, b);
  }

  struct DmgResp resp = {y, dydx};
  return resp;
}

struct DmgResp
dmglib_OBeta(double x, struct DmgEvol dp, struct DmgFrac Frac)
{
  double y, dydx;

  y = betacdf(x, dp.theParam[0], dp.theParam[1]);
  dydx = betapdf(x, dp.theParam[0], dp.theParam[1]);

  struct DmgResp resp = {y, dydx};
  return resp;
}

struct DmgResp
dmglib_Wbl(double x, struct DmgEvol dp, struct DmgFrac Frac)
{
  double y, dydx;

  y = wblcdf(x, dp.theParam[0], dp.theParam[1]);
  dydx = wblpdf(x, dp.theParam[0], dp.theParam[1]);

  struct DmgResp resp = {y, dydx};
  return resp;
}

struct DmgResp
dmglib_Logn(double x, struct DmgEvol dp, struct DmgFrac Frac)
{
  double y, dydx;

  y = logncdf(x, dp.theParam[0], dp.theParam[1]);
  dydx = lognpdf(x, dp.theParam[0], dp.theParam[1]);

  struct DmgResp resp = {y, dydx};
  return resp;
}


/* Bilin
 * normalized bilinear damage evolution */
struct DmgResp
dmglib_Bilin(double x, struct DmgEvol dp, struct DmgFrac Frac)
{
  double y, dydx;
  double ylim, slop;

  /* consider residual strength ylim if third damage
   * parameter is present */
  if (dp.nParam > 2) {
    ylim = dp.theParam[2];
  } else {
    ylim = 1;
  }

  if (x <= 0) {
    y = 0;
    dydx = 0;
  } else if (x <= dp.theParam[0]) {
    slop = dp.theParam[1] / dp.theParam[0];
    y = x * slop;
    dydx = slop;
  } else if (x <= 1) {
    slop = (ylim - dp.theParam[1]) / (1 - dp.theParam[0]);
    y = dp.theParam[1] + (x - dp.theParam[0]) * slop;
    dydx = slop;
  } else {
    y = ylim;
    dydx = 0.0;
  }

  struct DmgResp resp = {y, dydx};
  return resp;
}


/* function Trilin
   normalized bilinear damage evolution */
struct DmgResp
dmglib_Trilin(double x, struct DmgEvol dp, struct DmgFrac Frac)
{
  double y, dydx;
  double ylim, slop;

  /* consider residual strength ylim if fifth damage parameter is present */
  if (dp.nParam > 4) {
    ylim = dp.theParam[4];
  } else {
    ylim = 1;
  }

  if (x <= 0.0) {
    y = 0.0;
    dydx = 0.0;
  } else if (x <= dp.theParam[0]) {
    slop = dp.theParam[1] / dp.theParam[0];
    y = x * slop;
    dydx = slop;
  } else if (x <= dp.theParam[2]) {
    slop =
        (dp.theParam[3] - dp.theParam[1]) / (dp.theParam[2] - dp.theParam[0]);
    y = dp.theParam[1] + (x - dp.theParam[0]) * slop;
    dydx = slop;
  } else if (x <= 1) {
    slop = (ylim - dp.theParam[3]) / (1 - dp.theParam[2]);
    y = dp.theParam[3] + (x - dp.theParam[2]) * slop;
    dydx = slop;
  } else {
    y = ylim;
    dydx = 0;
  }

  struct DmgResp resp = {y, dydx};
  return resp;
}
