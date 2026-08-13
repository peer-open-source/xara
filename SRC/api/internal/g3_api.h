/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
//
// Written: cmp
//
#if !defined(G3_API_H_) && defined(OPS_USE_RUNTIME)
#define G3_API_H_
   struct Tcl_Interp;
   typedef void* ClientData;
#  include <api/elementAPI.h>
//
#  ifndef OPS_Export
#    define OPS_Export
#  endif
//
// ANALYSIS RUNTIME
//
#  undef  OPS_GetDomain
#  define OPS_GetDomain()        G3_getDomain(rt)
#  undef  OPS_GetAnalysisModel
#  define OPS_GetAnalysisModel() G3_getAnalysisModelPtr(rt)
//
// MODELBUILDER
//
// Materials
#  undef  OPS_GetUniaxialMaterial
#  define OPS_GetUniaxialMaterial(tag) G3_getUniaxialMaterialInstance(rt, (tag))
#  define OPS_getUniaxialMaterial(tag) G3_getUniaxialMaterialInstance(rt, (tag))

// Time series
#  define OPS_getTimeSeries(tag) G3_getTimeSeries(rt, (tag))


// Coordinate Transforms
#  undef  OPS_getCrdTransf
#  define OPS_getCrdTransf(tag) G3_getCrdTransf(rt, (tag))
#  undef  OPS_GetCrdTransf
#  define OPS_GetCrdTransf(tag) G3_getCrdTransf(rt, (tag))

// Section
SectionForceDeformation *G3_getSectionForceDeformation(G3_Runtime *, int);
#  undef  OPS_getSectionForceDeformation
#  define OPS_getSectionForceDeformation(tag) G3_getSectionForceDeformation(rt, (tag))
#  undef  OPS_GetSectionForceDeformation
#  define OPS_GetSectionForceDeformation(tag) G3_getSectionForceDeformation(rt, (tag))

// NDMaterial
NDMaterial *G3_GetNDMaterial(G3_Runtime *, int);
#  undef  OPS_getNDMaterial
#  define OPS_getNDMaterial(tag) G3_GetNDMaterial(rt, (tag))
#  undef  OPS_GetNDMaterial
#  define OPS_GetNDMaterial(tag) G3_GetNDMaterial(rt, (tag))

typedef int G3_Tag;
#  define G3_Char TCL_Char
#  define G3_getDouble Tcl_GetDoubleFromObj
class G3_Runtime;
class ModelBuilder;
#  include <runtime/runtime/ModelRegistry.h>
class AnalysisModel;
class EquiSolnAlgo;
class ConstraintHandler;
class DOF_Numberer;
class LinearSOE;
class EigenSOE;
class AnalysisModel;
class StaticAnalysis;
class DirectIntegrationAnalysis;
class StaticIntegrator;
class TransientIntegrator;
class ConvergenceTest;

class TimeSeries;

class UniaxialMaterial;
class NDMaterial;
class SectionForceDeformation;
class CrdTransf;
class FrictionModel;
class LimitCurve;
class Domain;
class FE_Datastore;


#ifdef __cplusplus
extern "C" {
#endif
//
// Runtime
//
G3_Runtime *G3_getRuntime(Tcl_Interp *);
Tcl_Interp *G3_getInterpreter(G3_Runtime*);

// Domain
Domain *G3_getDomain(G3_Runtime *);

ModelRegistry *G3_getSafeBuilder(G3_Runtime *);
ModelRegistry *G3_getModelBuilder(G3_Runtime *);
int             G3_setModelBuilder(G3_Runtime *, ModelRegistry*);
int     G3_getNDM(G3_Runtime *);

// Materials
UniaxialMaterial *G3_getUniaxialMaterialInstance(G3_Runtime *, G3_Tag);
int G3_addUniaxialMaterial(G3_Runtime *, UniaxialMaterial *);
// Coordinate Transforms
CrdTransf *G3_getCrdTransf(G3_Runtime *, G3_Tag);


// Systems and Solvers
LinearSOE **G3_getLinearSoePtr(G3_Runtime* );

//
// RUNTIME
//
// Time Series
TimeSeries *G3_getTimeSeries(G3_Runtime *, G3_Tag);
// int         G3_removeTimeSeries(G3_Runtime *, G3_Tag);

// Analysis
AnalysisModel **G3_getAnalysisModelPtr(G3_Runtime *);
StaticIntegrator *G3_getStaticIntegrator(G3_Runtime *);
int G3_setStaticIntegrator(G3_Runtime *, StaticIntegrator *);

#ifdef __cplusplus
}
#endif
#endif // G3_API_H_
