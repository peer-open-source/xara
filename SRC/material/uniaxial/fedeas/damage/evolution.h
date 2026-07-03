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

#ifdef __cplusplus
extern "C" {
#endif

/* Damage response; definitely rename */
struct DmgResp {double y, dydx;};

struct DmgResp dmglib_Wbl(double, struct DmgEvol, struct DmgFrac);
struct DmgResp dmglib_Logn(double, struct DmgEvol, struct DmgFrac);
struct DmgResp dmglib_none(double, struct DmgEvol, struct DmgFrac);
struct DmgResp dmglib_MBeta(double, struct DmgEvol, struct DmgFrac);
struct DmgResp dmglib_OBeta(double, struct DmgEvol, struct DmgFrac);
struct DmgResp dmglib_Bilin(double, struct DmgEvol, struct DmgFrac);
struct DmgResp dmglib_Trilin(double, struct DmgEvol, struct DmgFrac);

#ifdef __cplusplus
}
#endif
