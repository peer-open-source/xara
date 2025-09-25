
#include <Parsing.h>
#include <Logging.h>
#include <string.h>
#include <vector>

#include <elementAPI.h>
#include "IGASurfacePatch.h"
#include "IGAKLShell.h"
#include "IGAKLShell_BendingStrip.h"

// static void *OPS_IGASurfacePatch(void);

extern "C" int OPS_ResetInputNoBuilder(ClientData clientData, Tcl_Interp *interp, int cArg,
                          int mArg, TCL_Char ** const argv, Domain *domain);

static void
PrintSyntax()
{
    opserr << "IGA Patch $tag $nodeStartTag $P $Q $noPtsX $noPtsY \n\
        -type KLShell \n\
        -nonLinearGeometry [0 or 1] \n\
        -planeStressMatTags [list of tags] \n\
        -gFact $gx $gy $gz \n\
        -theta [list of thetas] \n\
        -thickness [list of layer thicknesses] \n\
        -uKnot [list of uKnots] \n\
        -vKnot [list of vKnots] \n\
        -controlPts [list coordinates of control points, u-direction first] \n\
        " << endln;
}

static
void* OPS_ADD_RUNTIME_VPV(OPS_IGASurfacePatch)
{

    int tag = 0;
    int P = 0;
    int Q = 0;

    int noPtsX = 0;
    int noPtsY = 0;

    int numdata = 1;

    ShellType shtype = ShellType::KLShell;

    int sectionTag = 0;
    int nodeStartTag = 0;
    int elementStartTag = 1;

    int nonLinearGeometry = 1;
    Vector gFact(3);



    opserr << "Creating IGA Patch:" << endln;
    if (OPS_GetIntInput(&numdata, &tag) < 0){
      PrintSyntax();
      return 0;  
    } 
    opserr << "tag = " << tag << endln;
    if (OPS_GetIntInput(&numdata, &nodeStartTag) < 0){
      PrintSyntax();
      return 0;  
    } 
    opserr << "nodeStartTag = " << nodeStartTag << endln;
    if (OPS_GetIntInput(&numdata, &P) < 0){
      PrintSyntax();
      return 0;  
    } 
    opserr << "P = " << P << endln;
    if (OPS_GetIntInput(&numdata, &Q) < 0){
      PrintSyntax();
      return 0;  
    } 
    opserr << "Q = " << Q << endln;
    if (OPS_GetIntInput(&numdata, &noPtsX) < 0){
      PrintSyntax();
      return 0;  
    } 
    opserr << "noPtsX = " << noPtsX << endln;
    if (OPS_GetIntInput(&numdata, &noPtsY) < 0){
      PrintSyntax();
      return 0;  
    } 
    opserr << "noPtsY = " << noPtsY << endln;
    
    // opserr << "OPS_IGASurfacePatch OPS_GetNumRemainingInputArgs() = " << OPS_GetNumRemainingInputArgs() << endln;

    // get other inputs
    std::vector<double> theta_stdVector, 
                        thickness_stdVector, 
                        uKnot_stdVector, 
                        vKnot_stdVector, 
                        controlPts_stdVector;
    std::vector<NDMaterial*>  materials; //  matTags_stdVector;
    // int loc = 5; // Vector (knotVector) start position
    // int loc = 7;
    int loc = 8;
    while (OPS_GetNumRemainingInputArgs() > 0) {

        // next arg
        const char* arg = OPS_GetString();
        loc++;

        // opserr << "arg = " << arg << endln;

        // check arg
        if (strcmp(arg, "-type") == 0) { //matTags
            // const char* type = OPS_GetString();
            // loc++;
            const char* typestring = OPS_GetString();
            opserr << "type = " << typestring << endln;
            if (strcmp(typestring, "KLShell") == 0)
                shtype = ShellType::KLShell;
            else if (strcmp(typestring, "KLShell_BendingStrip") == 0)
                shtype = ShellType::KLShell_BendingStrip;
            else
                opserr << "IGASurfacePatch - Unknown shell of type " << typestring << endln;

            loc++;
        }
        else if (strcmp(arg, "-planeStressMatTags") == 0) { //matTags
            while (OPS_GetNumRemainingInputArgs() > 0) {
                int val;
                if (OPS_GetIntInput(&numdata, &val) < 0) {
                    OPS_ResetCurrentInputArg(loc);
                    break;
                }

//              matTags_stdVector.push_back(val);

                NDMaterial* mat = nullptr;
                if ((mat = OPS_getNDMaterial(val)) {
                    materials.push_back(mat);
                } else {
                    opserr << "Failed to get material with tag " << val << "\n";
                    return nullptr;
                }

                loc++;
            }

        }
        else if (strcmp(arg, "-theta") == 0) { //theta
            while (OPS_GetNumRemainingInputArgs() > 0) {
                double val;
                if (OPS_GetDoubleInput(&numdata, &val) < 0) {
                    OPS_ResetCurrentInputArg(loc);
                    break;
                }
                theta_stdVector.push_back(val);
                loc++;
            }

        }
        else if (strcmp(arg, "-thickness") == 0) { //thickness
            while (OPS_GetNumRemainingInputArgs() > 0) {
                double val;
                if (OPS_GetDoubleInput(&numdata, &val) < 0) {
                    OPS_ResetCurrentInputArg(loc);
                    break;
                }
                thickness_stdVector.push_back(val);
                loc++;
            }

        }
        else if (strcmp(arg, "-uKnot") == 0) { //uKnot
            while (OPS_GetNumRemainingInputArgs() > 0) {
                double val;
                if (OPS_GetDoubleInput(&numdata, &val) < 0) {
                    OPS_ResetCurrentInputArg(loc);
                    break;
                }
                uKnot_stdVector.push_back(val);
                loc++;
            }

        }
        else if (strcmp(arg, "-vKnot") == 0) { //vKnot
            while (OPS_GetNumRemainingInputArgs() > 0) {
                double val;
                if (OPS_GetDoubleInput(&numdata, &val) < 0) {
                    OPS_ResetCurrentInputArg(loc);
                    break;
                }
                vKnot_stdVector.push_back(val);
                loc++;
            }

        }
        else if (strcmp(arg, "-controlPts") == 0) { //controlPts
            while (OPS_GetNumRemainingInputArgs() > 0) {
                double val;
                if (OPS_GetDoubleInput(&numdata, &val) < 0) {
                    OPS_ResetCurrentInputArg(loc);
                    break;
                }
                controlPts_stdVector.push_back(val);
                loc++;
            }
        }
        // else if (strcmp(arg, "-type") == 0) {
        //     const char* typestring = OPS_GetString();
        //     if (strcmp(typestring, "KLShell") == 0)
        //     {
        //         shtype = ShellType::KLShell;
        //     }
        //     else if (strcmp(typestring, "KLShell_BendingStrip") == 0)
        //     {
        //         shtype = ShellType::KLShell_BendingStrip;
        //     }
        //     else
        //     {
        //         opserr << "IGASurfacePatch - Unknown shell of type " << typestring << endln;
        //     }
        //     loc++;
        // }
        else if (strcmp(arg, "-sectionTag") == 0) {
            numdata = 1;

            if (OPS_GetIntInput(&numdata, &sectionTag) < 0) return 0;
            loc++;
        }
        else if (strcmp(arg, "-nodeStartTag") == 0) {
            numdata = 1;

            if (OPS_GetIntInput(&numdata, &nodeStartTag) < 0) return 0;
            loc++;
        }
        else if (strcmp(arg, "-elementStartTag") == 0) {
            numdata = 1;

            if (OPS_GetIntInput(&numdata, &elementStartTag) < 0) return 0;
            loc++;
        }
        else if (strcmp(arg, "-nonLinearGeometry") == 0) {
            numdata = 1;

            if (OPS_GetIntInput(&numdata, &nonLinearGeometry) < 0) return 0;
            loc++;
        }
        else if (strcmp(arg, "-gFact") == 0) {
            int i = 0;
            while (OPS_GetNumRemainingInputArgs() > 0) {
                double val;
                if (OPS_GetDoubleInput(&numdata, &val) < 0) {
                    OPS_ResetCurrentInputArg(loc);
                    break;
                }
                gFact(i) = val;
                i++;
                loc++;
            }
        }


    }

//  ID matTags(&matTags_stdVector[0], (int)(matTags_stdVector.size()));
    Vector theta(&theta_stdVector[0], (int)(theta_stdVector.size()));
    Vector thickness(&thickness_stdVector[0], (int)(thickness_stdVector.size()));
    Vector uKnot(&uKnot_stdVector[0], (int)uKnot_stdVector.size());
    Vector vKnot(&vKnot_stdVector[0], (int)vKnot_stdVector.size());


    int controlPts_size = controlPts_stdVector.size();
    int M = controlPts_size / 4;
    int N = 4;

    // Matrix controlPts(&controlPts_stdVector[0], M, N);
    Matrix controlPts(&controlPts_stdVector[0], N, M);

    IGASurfacePatch* patch = new IGASurfacePatch(tag, nodeStartTag, P, Q, noPtsX, noPtsY, nonLinearGeometry, gFact, materials, theta, thickness, uKnot, vKnot, controlPts, shtype);


    Domain* theDomain = OPS_GetDomain();
    theDomain->addElement(patch);

    // return 0;
    return patch;
}

int
TclCommand_IGA(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  G3_Runtime *rt = G3_getRuntime(interp);

  if (argc < 3) {
    opserr << "IGA_Command : IGA <cmd> <args...>" << endln;
    return -1; 
  }

  if (strcmp(argv[1],"Patch") == 0) {
    // OPS_ResetInput(clientData, interp, 2, argc, argv, theTclDomain, theTclBuilder);
    OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, 0);
    OPS_IGASurfacePatch(rt, argc, argv);
  }
  
  return 0; 
}


