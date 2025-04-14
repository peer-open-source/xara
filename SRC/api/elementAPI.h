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
// Written: fmk
// $Date: 2010-03-05 22:32:36 $
// /usr/local/cvs/OpenSees/SRC/api/elementAPI.h
//
#ifndef _eleAPI
#define _eleAPI

#define OPS_Error ops_error_
#define OPS_GetIntInput ops_getintinput_
#define OPS_GetDoubleInput ops_getdoubleinput_
#define OPS_SetDoubleListsOutput ops_setdoublelistsoutput_
#define OPS_SetDoubleDictOutput ops_setdoubledictoutput_
#define OPS_SetDoubleDictListOutput ops_setdoubledictlistoutput_
//
#define OPS_AllocateMaterial ops_allocatematerial_
#define OPS_AllocateElement ops_allocateelement_
#define OPS_GetMaterialType ops_getmaterialtype_
#define OPS_GetMaterial ops_getmaterial_
#define OPS_GetMaterialPtr ops_getmaterialptr_
#define OPS_GetCrdTransf ops_getcrdtransf_
#define OPS_GetNodeCrd ops_getnodecrd_
#define OPS_GetNodeDisp ops_getnodedisp_
#define OPS_GetNodeVel ops_getnodevel_
#define OPS_GetNodeAccel ops_getnodeaccel_
#define OPS_GetNodeIncrDisp ops_getnodeincrdisp_
#define OPS_GetNodeIncrDeltaDisp ops_getnodeincrdeltadisp_
#define OPS_GetInt ops_getintinput_
#define OPS_GetDouble ops_getdoubleinput_
#define OPS_GetString ops_getstring
#define OPS_GetStringFromAll ops_getstringfromall_
#define OPS_SetString ops_setstring
#define OPS_GetNDM ops_getndm_
#define OPS_GetNDF ops_getndf_
#define OPS_GetFEDatastore ops_getfedatastore_
#define OPS_GetInterpPWD ops_getinterppwd_

#define OPS_GetAnalysisModel ops_getanalysismodel_
#define OPS_GetNumEigen ops_getnumeigen_
#define OPS_GetStaticIntegrator ops_getstaticintegrator_
#define OPS_builtModel ops_builtmodel_
#define OPS_GetDomain ops_getdomain_

#include <OPS_Globals.h>
#include <logging/Logging.h>
#ifdef __cplusplus
#include <map>
#include <vector>
class AnalysisModel;
class LinearSOE;

class UniaxialMaterial;
class NDMaterial;
class SectionForceDeformation;
class CrdTransf;
class Domain;
class FE_Datastore;
//
extern "C" int  OPS_SetDoubleListsOutput(std::vector<std::vector<double>>& data);
extern "C" int  OPS_SetDoubleDictOutput(std::map<const char*, double>& data);
extern "C" int  OPS_SetDoubleDictListOutput(std::map<const char*, std::vector<double>>& data);

extern NDMaterial* OPS_GetNDMaterial(int matTag);
extern SectionForceDeformation* OPS_GetSectionForceDeformation(int secTag);
extern CrdTransf* OPS_GetCrdTransf(int crdTag);
// extern FrictionModel* OPS_GetFrictionModel(int frnTag);

extern FE_Datastore* OPS_GetFEDatastore();

extern "C" bool* OPS_builtModel(void);
int OPS_numIter();

extern "C" {
#endif // __cplusplus

int         OPS_GetNDM();
int         OPS_GetNDF();
int         OPS_Error(const char* errorMessage, int length);
int         OPS_GetNumRemainingInputArgs();
int         OPS_ResetCurrentInputArg(int cArg);
int         OPS_GetIntInput(int* numData, int* data);
int         OPS_GetDoubleInput(int* numData, double* data);
int         OPS_GetDoubleListInput(int* size, Vector * data);
int         OPS_EvalDoubleStringExpression(const char* theExpression, double& current_val);
const char* OPS_GetString(); // does a strcpy
const char* OPS_GetStringFromAll(char* buffer, int len); // does a strcpy
int         OPS_SetString(const char* str);
int         OPS_GetStringCopy(char** cArray); // returns a new copy
const char* OPS_GetInterpPWD();

//
// extern "C" int       OPS_ResetInput(ClientData clientData, Tcl_Interp * interp, int cArg, int mArg, TCL_Char * *argv, Domain * domain, TclModelBuilder * builder);
// extern "C" int       OPS_ResetInputNoBuilder(ClientData clientData, Tcl_Interp * interp, int cArg, int mArg, TCL_Char * *argv, Domain * domain);
// extern "C" int       OPS_GetString(char *cArray, int sizeArray); // does a strcpy


#ifdef __cplusplus
}
#endif // __cplusplus

#endif // _eleAPI

#include <api/runtimeAPI.h>
