#include <tcl.h>
#include <LimitCurve.h>
#include <elementAPI.h>
#include <LimitCurveAPI.h>
#include <packages.h>

static LimitCurveFunction *theLimitCurveFunctions = NULL;

#if 0
LimitCurve *
Tcl_addWrapperLimitCurve(limCrvObj *theLimCrv, ClientData clientData,
                         Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  theInterp = interp;

  //  theModelBuilder = builder;
  currentArgv = argv;
  currentArg = 2;
  maxArg = argc;

  // get the current load factor
  static modelState theModelState;
  if (theDomain != 0) {
    double time = theDomain->getCurrentTime();
    double dt = theDomain->getCurrentTime() - time;
    theModelState.time = time;
    theModelState.dt = dt;
  }

  // invoke the limit curve function with isw = 0
  int isw = ISW_INIT;
  int result;
  theLimCrv->limCrvFunctPtr(theLimCrv, &theModelState, 0, 0, 0, &isw, &result);

  if (result != 0) {
    opserr << "Tcl_addWrapperLimitCurve - failed in limit curve function "
           << result << endln;
    return 0;
  }

  WrapperLimitCurve *theLimitCurve = new WrapperLimitCurve(argv[1], theLimCrv);

  return theLimitCurve;
}
#endif

/*
typedef struct limitCurveFunction {
  char *funcName;
  limCrvFunct theFunct;
  struct limitCurveFunction *next;
} LimitCurveFunction;
*/

extern "C" limCrvObj *
OPS_GetLimitCurveType(char *type, int sizeType)
{

  // try existing loaded routines
  LimitCurveFunction *limCrvFunction = theLimitCurveFunctions;
  bool found = false;
  while (limCrvFunction != NULL && found == false) {
    if (strcmp(type, limCrvFunction->funcName) == 0) {

      // create a new eleObject, set the function ptr &  return it

      limCrvObj *theLimCrvObject = new limCrvObj;
      theLimCrvObject->limCrvFunctPtr = limCrvFunction->theFunct;
      /* opserr << "limCrvObj *OPS_GetLimitCurveType() - FOUND " << endln;  */
      return theLimCrvObject;
    } else
      limCrvFunction = limCrvFunction->next;
  }

  // try to load new routine from dynamic library in load path
  limCrvFunct limCrvFunctPtr;
  void *libHandle;
  int res =
      getLibraryFunction(type, type, &libHandle, (void **)&limCrvFunctPtr);

  if (res == 0) {
    // add the routine to the list of possible elements
    char *funcName = new char[strlen(type) + 1];
    strcpy(funcName, type);
    limCrvFunction = new LimitCurveFunction;
    limCrvFunction->theFunct = limCrvFunctPtr;
    limCrvFunction->funcName = funcName;
    limCrvFunction->next = theLimitCurveFunctions;
    theLimitCurveFunctions = limCrvFunction;

    // create a new eleObject, set the function ptr &  return it
    limCrvObj *theLimCrvObject = new limCrvObj;
    theLimCrvObject->limCrvFunctPtr = limCrvFunction->theFunct;
    return theLimCrvObject;
  }

  return 0;
}

extern "C" int
OPS_AllocateLimitCurve(limCrvObject *theLimCrv)
{

  /*fprintf(stderr,"allocateLimitCurve Address %p\n",theLimCrv);*/

  if (theLimCrv->nParam > 0)
    theLimCrv->theParam = new double[theLimCrv->nParam];

  int nState = theLimCrv->nState;

  if (nState > 0) {
    theLimCrv->cState = new double[nState];
    theLimCrv->tState = new double[nState];
    for (int i = 0; i < nState; i++) {
      theLimCrv->cState[i] = 0;
      theLimCrv->tState[i] = 0;
    }
  } else {
    theLimCrv->cState = 0;
    theLimCrv->tState = 0;
  }

  return 0;
}
