/***********************************************************************
*  direction_logic.h   – single copy of the X/Y partitioning algorithm *
***********************************************************************/
#pragma once
#include <cmath>
inline int sgn(double x) { return (x > 0) - (x < 0); }

struct DirectionState
{
    /*=== I/O STATE (all passed by reference) =======================*/
    double& V;                        // current shear force (was f)
    double&    slopeSign;                // +1 up, -1 down
    double&    loadingFlag;              // 1 after first “bottom→top” pass
    double&    unloadingFlag;            // 1 after first “top→bottom” pass
    double&    signV;                    // saved sign of first cycle (+/-)

    /* displacement DOFs */
    double& u1; double& u2; double& u3; double& u4;

    /* running extreme / turning‐point data */
    double& V_tr;     double& V_tr_u;      // last stored turning forces
    double& V_ref;                         // reference force for current loop

    double& u1_tr; double& u2_tr; double& u3_tr; double& u4_tr;
    double& u1_tr_u; double& u2_tr_u; double& u3_tr_u; double& u4_tr_u;

    /* stored, pre-step displacements */
    double& u1_st; double& u2_st; double& u3_st; double& u4_st;
};

/* utility lambdas local to this TU */
namespace _dirdetail {
    inline double lim(double val, double lim) {
        return std::fabs(val) > lim ? sgn(val) * lim : val;
    }
}

/**
 * Implements the entire “regime 0 … regime 5” logic for one direction.
 * All parameters come from the caller – the function is completely stateless.
 */
inline void updateDirection(
        DirectionState s,
        /* thresholds */ double F_f1, double F_f2, double F_f4,
        double F_dr1, double F_dr4,
        /* geometry */  double Wavg, double L1, double L2, double L3,
        /* capacity  */ double Ubar1, double Ubar2, double Ubar3)
{
    using _dirdetail::lim;

    /*-------------------------------------------------------------
     * 1.  Capture the sign on the very first excursion
     *------------------------------------------------------------*/
    if (std::fabs(s.V) <= F_f2)         // before regime 1
        s.signV = s.slopeSign;          // remember direction of first push

    /*-------------------------------------------------------------
     * 2.  FIRST CYCLE  (loadingFlag==0 && unloadingFlag==0)
     *------------------------------------------------------------*/
    if (s.unloadingFlag == 0 && s.loadingFlag == 0 &&
        s.slopeSign      == s.signV)
    {
        s.signV = sgn(s.slopeSign);     // freeze sign for later loops

        /* === REGIME-selection for first cycle ===================*/
        const double Sv = s.signV * s.V;

        if (Sv <= F_f1) {                       // Regime 1
            s.u1 = s.u1_st;
            s.u4 = s.u4_st;
            s.u2 = (Sv - F_f2) * L1 / Wavg;
            s.u3 = s.u2;

            if (Sv <= F_f2) s.u2 = s.u3 = 0.0;
        }
        else if (Sv <= F_f4) {                  // Regime 2
            s.u1 = (Sv - F_f1) * L2 / Wavg;
            s.u4 = s.u4_st;
            s.u2 = s.u2_st;
            s.u3 = (Sv - F_f2) * L1 / Wavg;
        }
        else if (Sv <= F_dr1) {                 // Regime 3
            s.u1 = (Sv - F_f1) * L2 / Wavg;
            s.u4 = (Sv - F_f4) * L3 / Wavg;
            s.u2 = s.u2_st;
            s.u3 = s.u3_st;
        }
        else if (Sv <= F_dr4) {                 // Regime 4
            s.u1 = s.u1_st;
            s.u4 = (Sv - F_f4) * L3 / Wavg;
            s.u2 = ((Sv - F_f2) / Wavg - Ubar2 / L2) * L1;
            s.u3 = s.u3_st;
        }
        else {                                  // Regime 5
            s.u1 = s.u1_st;  s.u4 = s.u4_st;
            s.u2 = ((Sv - F_f2) / Wavg - Ubar2 / L2) * L1;
            s.u3 = ((Sv - F_f2) / Wavg - Ubar3 / L3) * L1;
            s.u2 = lim(s.u2, Ubar1);
            s.u3 = lim(s.u3, Ubar1);
        }

        /* store turning point for the upcoming unloading branch */
        s.V_tr  = s.V;   s.u1_tr  = s.u1;  s.u2_tr  = s.u2;
        s.u3_tr = s.u3;  s.u4_tr  = s.u4;
        return;  // first-cycle done
    }

    /*-------------------------------------------------------------
     * 3.  SUBSEQUENT LOOPS:   s.slopeSign  (+1 up or -1 down)
     *------------------------------------------------------------*/
    /* Which half of the hysteresis loop are we on? */
    const bool climbing   = (s.slopeSign == +1);   // bottom→top
    const bool descending = (s.slopeSign == -1);   // top→bottom

    /* Mark that we have completed at least one half-loop */
    if (climbing)   s.loadingFlag   = 1;
    if (descending) s.unloadingFlag = 1;

    /* ------------------------------------------------------------
     * Step A : choose the reference (turning) point for this half
     * -----------------------------------------------------------*/
    if (climbing) {
        if (s.signV == -1) {                           // first half climbed in −ve direction
            s.V_ref  = s.V_tr;
            s.u1_tr_u = s.u1_tr; s.u4_tr_u = s.u4_tr;
            s.u2_tr_u = s.u2_tr; s.u3_tr_u = s.u3_tr;
        } else {
            s.V_ref  = s.V_tr_u;
            s.u1_tr  = s.u1_tr_u; s.u4_tr  = s.u4_tr_u;
            s.u2_tr  = s.u2_tr_u; s.u3_tr  = s.u3_tr_u;
        }
    } else { // descending
        if (s.signV == -1) {
            s.V_ref  = s.V_tr_u;
        } else {
            s.V_ref  = s.V_tr;
        }
    }

    /* ------------------------------------------------------------
     * Step B :   common lambda that writes u1–u4 & history capture
     * -----------------------------------------------------------*/
    const auto SAVE_TR = [&](bool upper) {
        if (upper) {                         // store to _u  (descending loop)
            s.V_tr_u = s.V;
            s.u1_tr_u = s.u1; s.u4_tr_u = s.u4;
            s.u2_tr_u = s.u2; s.u3_tr_u = s.u3;
        } else {                             // store to    (climbing loop)
            s.V_tr = s.V;
            s.u1_tr = s.u1; s.u4_tr = s.u4;
            s.u2_tr = s.u2; s.u3_tr = s.u3;
        }
    };

    /* ------------------------------------------------------------
     * Step C :    Evaluate REGIME for current V, relative to V_ref
     *             A sign-agnostic formulation is easier:
     *    ΔV = (climbing) ? V - V_ref : V_ref - V
     * -----------------------------------------------------------*/
    const double dV = climbing ? (s.V - s.V_ref)
                               : (s.V_ref - s.V);

    /* signless thresholds for ΔV */
    const double R0 = 2*F_f2;
    const double R1 = 2*F_f1;
    const double R4 = 2*F_f4;

    /* helper to translate sign back */
    const double σ = (s.signV == 0) ? 1.0 : s.signV;

    /* ------------------------------------------------------------
     * REGIME MAP  (matches original tables 0‒5)
     * -----------------------------------------------------------*/
    if ( (climbing  && s.V <= s.V_ref + R0) ||
         (descending && s.V >= s.V_ref - R0) )
    {
        /* Regime 0 : purely unloading — stay where you are */
        s.u1 = s.u1_st; s.u2 = s.u2_st;
        s.u3 = s.u3_st; s.u4 = s.u4_st;
        SAVE_TR(descending);
        return;
    }

    /* ---- Regime 1 ------------------------------------------------*/
    if ( (climbing  && dV <=  R1) ||
         (descending && dV <=  R1) )
    {
        s.u1 = s.u1_st;  s.u4 = s.u4_st;
        s.u2 = s.u2_st + σ * dV * L1 / Wavg;
        s.u3 = s.u3_st + σ * dV * L1 / Wavg;
        s.u2 = lim(s.u2, Ubar1);
        s.u3 = lim(s.u3, Ubar1);
        SAVE_TR(descending);
        return;
    }

    /* ---- Regime 2 ------------------------------------------------*/
    if (dV <= R4) {
        s.u1 = s.u1_st + σ * (dV - R1) * L2 / Wavg;
        s.u4 = s.u4_st;
        s.u2 = s.u2_st;
        s.u3 = s.u3_st + σ * dV * L1 / Wavg;
        s.u3 = lim(s.u3, Ubar1);
        s.u1 = lim(s.u1, Ubar2);
        SAVE_TR(descending);
        return;
    }

    /* ---- Regime 3 ------------------------------------------------*/
    if (( climbing && ((s.signV*s.V) <= F_dr1) ) ||
        (descending && ((s.signV*s.V) >= -F_dr1)) )
    {
        s.u1 = s.u1_st + σ * (dV - R1) * L2 / Wavg;
        s.u4 = s.u4_st + σ * (dV - R4) * L3 / Wavg;
        s.u2 = s.u2_st;  s.u3 = s.u3_st;
        s.u1 = lim(s.u1, Ubar2);
        s.u4 = lim(s.u4, Ubar3);
        SAVE_TR(descending);
        return;
    }

    /* ---- Regime 4 ------------------------------------------------*/
    if (( s.signV*s.V ) * climbing >= F_dr1 * σ &&
        ( s.signV*s.V ) * climbing <= F_dr4 * σ)
    {
        s.u1 = s.u1_st;
        s.u4 = s.u4_st + σ * (s.V - σ*F_f4) * L3 / Wavg;
        s.u2 = s.u2_st + σ*((s.V - σ*F_f2)/Wavg - Ubar2/L2)*L1;
        s.u3 = s.u3_st;
        s.u2 = lim(s.u2, Ubar1);
        s.u4 = lim(s.u4, Ubar3);
        SAVE_TR(descending);
        return;
    }

    /* ---- Regime 5 ------------------------------------------------*/
    /* Anything above F_dr4 (climbing) or below −F_dr4 (descending) */
    s.u1 = s.u1_st;  s.u4 = s.u4_st;
    s.u2 = s.u2_st + σ * (dV - R0) * L1 / Wavg;
    s.u3 = s.u3_st + σ * (dV - R0) * L1 / Wavg;
    s.u2 = lim(s.u2, Ubar1);
    s.u3 = lim(s.u3, Ubar1);
    SAVE_TR(descending);
}
