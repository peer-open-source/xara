/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** With a lot of additions by                                         **
**   Boris Jeremic    (jeremic@ucdavis.edu)                           **
**   Zaohui Yang      (zhyang@ucdavis.edu)                            **
**   Zhao Cheng       (zcheng@ucdavis.edu)                            **
**                                                                    **
** ****************************************************************** */
//
#include <tcl.h>
#include <string.h>
#include <Logging.h>
#include <Parsing.h>
#include <ModelRegistry.h>

#include <PressureDependentElastic3D.h>

#include <PlaneStressMaterial.h>
#include <PlateFiberMaterial.h>
#include <BeamFiberMaterial.h>

#include <PressureIndependMultiYield.h>
#include <PressureDependMultiYield.h>
#include <PressureDependMultiYield02.h>
#include <FluidSolidPorousMaterial.h>
// #include <isotropy.h>

// #include <Template3Dep.h>
// #include <NewTemplate3Dep.h>
// #include <FiniteDeformationElastic3D.h>
// #include <FiniteDeformationEP3D.h>


Tcl_CmdProc TclCommand_newPlasticMaterial;
Tcl_CmdProc TclCommand_newElasticMaterial;
Tcl_CmdProc TclCommand_newConcreteMaterial;
// Tcl_CmdProc TclCommand_newIsotropicMaterial;



int
TclCommand_addMaterial(ClientData clientData, Tcl_Interp* interp, 
                        Tcl_Size argc, TCL_Char** const argv)
{
  static
  std::unordered_map<std::string, Tcl_CmdProc*> MaterialLibrary = {
    {"ElasticIsotropic",          TclCommand_newElasticMaterial},
    {"Elastic",                   TclCommand_newElasticMaterial},
    {"Isotropic",                 TclCommand_newElasticMaterial},
    {"J2",                        TclCommand_newPlasticMaterial},
    {"J2Simplified",              TclCommand_newPlasticMaterial},
    {"J2BeamFiber",               TclCommand_newPlasticMaterial},
    {"GeneralizedJ2",             TclCommand_newPlasticMaterial},
    {"NonlinearJ2",               TclCommand_newPlasticMaterial},

    {"PlasticDamageConcrete",     TclCommand_newConcreteMaterial},
    {"FariaPlasticDamage",        TclCommand_newConcreteMaterial},
  };


  if (argc < 2) {
    opserr << OpenSees::PromptValueError
           << "missing argument type"
           << "\n";
    return TCL_ERROR;
  }

  auto cmd = MaterialLibrary.find(std::string(argv[1]));
  if (cmd != MaterialLibrary.end())
    return (*cmd->second)(clientData, interp, argc, &argv[0]);

  opserr << OpenSees::PromptValueError
         << "unknown material type "
         << argv[1]
         << "\n";
  return TCL_ERROR;

}




#if 0
Template3Dep* TclModelBuilder_addTemplate3Dep(ClientData clientData, Tcl_Interp* interp, int argc,
                                              TCL_Char** argv, TclModelBuilder* theTclBuilder,
                                              int eleArgStart);

NewTemplate3Dep* TclModelBuilder_addNewTemplate3Dep(ClientData clientData, Tcl_Interp* interp,
                                                    int argc, TCL_Char** argv,
                                                    TclModelBuilder* theTclBuilder,
                                                    int eleArgStart);

FiniteDeformationElastic3D*
TclModelBuilder_addFiniteDeformationElastic3D(ClientData clientData, Tcl_Interp* interp, int argc,
                                              TCL_Char** argv, TclModelBuilder* theTclBuilder,
                                              int eleArgStart);

FiniteDeformationEP3D* TclModelBuilder_addFiniteDeformationEP3D(ClientData clientData,
                                                                Tcl_Interp* interp, int argc,
                                                                TCL_Char** argv,
                                                                TclModelBuilder* theTclBuilder,
                                                                int eleArgStart);

NDMaterial* TclModelBuilder_addFeapMaterial(ClientData clientData, Tcl_Interp* interp, int argc,
                                            TCL_Char** argv, TclModelBuilder* theTclBuilder);




template <typename MatType>
int
TclCommand_newMinMaxND(ClientData clientData, Tcl_Interp* interp, int argc, const char**const argv)
{
    if (argc < 4) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: uniaxialMaterial MinMax tag? matTag?";
      opserr << " <-min min?> <-max max?>" << "\n";
      return TCL_ERROR;
    }

    int tag, matTag;
    
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid uniaxialMaterial MinMax tag" << "\n";
      return TCL_ERROR;               
    }

    if (Tcl_GetInt(interp, argv[3], &matTag) != TCL_OK) {
      opserr << "WARNING invalid component tag\n";
      opserr << "uniaxialMaterial MinMax: " << tag << "\n";
      return TCL_ERROR;
    }

    // Search for min and max strains
    double epsmin = NEG_INF_STRAIN;
    double epsmax = POS_INF_STRAIN;
      
    for (int j = 4; j < argc; j++) {
      if (strcmp(argv[j],"-min") == 0) {
        if ((j+1) >= argc || Tcl_GetDouble (interp, argv[j+1], &epsmin) != TCL_OK) {
          opserr << "WARNING invalid min\n";
          opserr << "uniaxialMaterial MinMax: " << tag << "\n";
          return TCL_ERROR;
        }
        j++;
      }
      if (strcmp(argv[j],"-max") == 0) {
        if ((j+1) >= argc || Tcl_GetDouble (interp, argv[j+1], &epsmax) != TCL_OK) {
          opserr << "WARNING invalid max\n";
          opserr << "uniaxialMaterial MinMax: " << tag << "\n";
          return TCL_ERROR;
        }
        j++;
      }
    }
      
    UniaxialMaterial *theMat = theTclBuilder->getUniaxialMaterial(matTag);

    if (theMat == 0) {
      opserr << "WARNING component material does not exist\n";
      opserr << "Component material: " << matTag; 
      opserr << "\nuniaxialMaterial MinMax: " << tag << "\n";
      return TCL_ERROR;
    }

    // Parsing was successful, allocate the material
    theMaterial = new MinMaxNDMaterial(tag, *theMat, epsmin, epsmax);    
}
#endif