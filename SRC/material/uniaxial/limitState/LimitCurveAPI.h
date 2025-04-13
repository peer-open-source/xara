#ifndef LIMIT_CURVE_API_H
#define LIMIT_CURVE_API_H
#include <LimitCurve.h>
#include <tcl.h>
#include <elementAPI.h>
#include <OpenSeesFFI.h>

struct limCrvObject;

typedef void (*limCrvFunct)(struct limCrvObject*, modelState*, double* strain, double* tang, double* stress, int* isw, int* error);

typedef struct limitCurveFunction {
    char* funcName;
    limCrvFunct theFunct;
    struct limitCurveFunction* next;
} LimitCurveFunction;

struct limCrvObject {
    int tag;
    int nParam;
    int nState;
    double* theParam;
    double* cState;
    double* tState;
    limCrvFunct limCrvFunctPtr;
    void* limCrvObjectPtr;
};

typedef struct limCrvObject limCrvObj;

#define OPS_AllocateLimitCurve ops_allocatelimitcurve_
#define OPS_GetLimitCurveType ops_getlimitcurvetype_

// // extern LimitCurveFunction *theLimitCurveFunctions;
// limCrv* OPS_GetLimitCurveType(char* type, int sizeType);//**MRL
// int     OPS_AllocateLimitCurve(limCrvObj*);//**MRL

#ifdef __cplusplus
extern "C" limCrvObj*  OPS_GetLimitCurveType(char* type, int sizeType);
extern "C" int         OPS_AllocateLimitCurve(limCrvObject * theLimCrv);
#endif

#endif // LIMIT_CURVE_API_H
