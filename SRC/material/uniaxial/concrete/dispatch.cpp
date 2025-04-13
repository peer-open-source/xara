#include <tcl.h>
#include <OPS_Globals.h>
#include <BasicModelBuilder.h>
#include <Concrete04.h>
#include <Concrete06.h>
#include <Concrete07.h>

int 
TclCommand_newUniaxialConcrete04(ClientData cd, Tcl_Interp* interp, int argc, TCL_Char ** const argv)
{
    UniaxialMaterial *theMaterial = nullptr;
    int tag;
    double fpc, epsc0, ft, epscu, Ec0, etu, beta;

    if (argc != 10 && argc != 9 && argc != 7) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: uniaxialMaterial Concrete04 tag? fpc? epsc0? epscu? "
                "Ec0? <ft? etu? <beta?> >"
             << endln;
      return TCL_ERROR;
    }


    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid uniaxialMaterial Concrete04 tag" << endln;
      return TCL_ERROR;
    }

    // Read required Concrete04 material parameters

    if (Tcl_GetDouble(interp, argv[3], &fpc) != TCL_OK) {
      opserr << "WARNING invalid fpc\n";
      opserr << "Concrete04 material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &epsc0) != TCL_OK) {
      opserr << "WARNING invalid epsc0\n";
      opserr << "Concrete04 material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &epscu) != TCL_OK) {
      opserr << "WARNING invalid epscu\n";
      opserr << "Concrete04 material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[6], &Ec0) != TCL_OK) {
      opserr << "WARNING invalid Ec0\n";
      opserr << "Concrete04 material: " << tag << endln;
      return TCL_ERROR;
    }
    if (argc == 9 || argc == 10) {
      if (Tcl_GetDouble(interp, argv[7], &ft) != TCL_OK) {
        opserr << "WARNING invalid ft\n";
        opserr << "Concrete04 material: " << tag << endln;
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[8], &etu) != TCL_OK) {
        opserr << "WARNING invalid etu\n";
        opserr << "Concrete04 material: " << tag << endln;
        return TCL_ERROR;
      }
    }
    if (argc == 10) {
      if (Tcl_GetDouble(interp, argv[9], &beta) != TCL_OK) {
        opserr << "WARNING invalid beta\n";
        opserr << "Concrete04 material: " << tag << endln;
        return TCL_ERROR;
      }
    }

    // Parsing was successful, allocate the material
    if (argc == 10) {
      theMaterial =
          new Concrete04(tag, fpc, epsc0, epscu, Ec0, ft, etu, beta);

    } else if (argc == 9) {
      theMaterial = new Concrete04(tag, fpc, epsc0, epscu, Ec0, ft, etu);

    } else if (argc == 7) {
      theMaterial = new Concrete04(tag, fpc, epsc0, epscu, Ec0);
    }

    return ((BasicModelBuilder*)cd)->addTaggedObject<UniaxialMaterial>(*theMaterial);
//  return theMaterial;
}


int
TclCommand_newUniaxialConcrete06(ClientData cd, Tcl_Interp* interp, int argc, TCL_Char ** const argv)
{
  UniaxialMaterial *theMaterial = nullptr;
      if (argc < 12) {
        opserr << "WARNING insufficient arguments\n";
        opserr
            << "Want: uniaxialMaterial Concrete06 tag? fc? eo? r? k? alphaC? "
               "fcr? ecr? b? alphaT?"
            << endln;
        return TCL_ERROR;
      }

      int tag;

      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
        opserr << "WARNING invalid uniaxialMaterial Concrete06 tag" << endln;
        return TCL_ERROR;
      }

      // Read required Concrete01 material parameters
      double fc, eo, r, k, fcr, ecr, b, alphaC, alphaT;

      if (Tcl_GetDouble(interp, argv[3], &fc) != TCL_OK) {
        opserr << "WARNING invalid fc\n";
        opserr << "Concrete06 material: " << tag << endln;
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[4], &eo) != TCL_OK) {
        opserr << "WARNING invalid eo\n";
        opserr << "Concrete06 material: " << tag << endln;
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[5], &r) != TCL_OK) {
        opserr << "WARNING invalid r\n";
        opserr << "Concrete06 material: " << tag << endln;
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[6], &k) != TCL_OK) {
        opserr << "WARNING invalid k\n";
        opserr << "Concrete06 material: " << tag << endln;
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[7], &alphaC) != TCL_OK) {
        opserr << "WARNING invalid alphaC\n";
        opserr << "Concrete06 material: " << tag << endln;
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[8], &fcr) != TCL_OK) {
        opserr << "WARNING invalid fcr\n";
        opserr << "Concrete06 material: " << tag << endln;
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[9], &ecr) != TCL_OK) {
        opserr << "WARNING invalid ecr\n";
        opserr << "Concrete06 material: " << tag << endln;
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[10], &b) != TCL_OK) {
        opserr << "WARNING invalid b\n";
        opserr << "Concrete06 material: " << tag << endln;
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[11], &alphaT) != TCL_OK) {
        opserr << "WARNING invalid alphaT\n";
        opserr << "Concrete06 material: " << tag << endln;
        return TCL_ERROR;
      }

      // Parsing was successful, allocate the material
      theMaterial =
          new Concrete06(tag, fc, eo, r, k, alphaC, fcr, ecr, b, alphaT);

    return ((BasicModelBuilder*)cd)->addTaggedObject<UniaxialMaterial>(*theMaterial);
//  return theMaterial;
}

int
TclCommand_newUniaxialConcrete07(ClientData cd, Tcl_Interp* interp, int argc, TCL_Char ** const argv)
{
    UniaxialMaterial *theMaterial = 0;
      // Check to see if there are enough arquements
      if (argc < 11) {
        opserr << "WARNING: Insufficient arguments\n";
        opserr << "Want: uniaxialMaterial Concrete07 tag? fpc? epsc0? Ec? fpt? "
                  "epst0? xcrp? xcrn? r?\n";
        return TCL_ERROR;
      }

      int tag;

      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
        opserr << "WARNING: Invalid uniaxial Concrete07 tag\n";
        return TCL_ERROR;
      }

      // Read in the faluves required for the model
      double fpc, epsc0, Ec, fpt, epst0, xcrp, xcrn, r;

      if (Tcl_GetDouble(interp, argv[3], &fpc) != TCL_OK) {
        opserr << "WARNING: Invalid peak compression stress\n";
        opserr << "uniaxialMaterial Concrete07: " << tag << endln;
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[4], &epsc0) != TCL_OK) {
        opserr << "WARNING: Invalid peak compression strain\n";
        opserr << "uniaxialMaterial Concrete07: " << tag << endln;
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[5], &Ec) != TCL_OK) {
        opserr << "WARNING: Invalid Young's Modulus\n";
        opserr << "uniaxialMaterial Concrete07: " << tag << endln;
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[6], &fpt) != TCL_OK) {
        opserr << "WARNING: Invalid peak tension stress\n";
        opserr << "uniaxialMaterial Concrete07: " << tag << endln;
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[7], &epst0) != TCL_OK) {
        opserr << "WARNING: Invalid peak tension strain\n";
        opserr << "uniaxialMaterial Concrete07: " << tag << endln;
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[8], &xcrp) != TCL_OK) {
        opserr
            << "WARNING: Invalid critical nondimensional strain in tension\n";
        opserr << "uniaxialMaterial Concrete07: " << tag << endln;
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[9], &xcrn) != TCL_OK) {
        opserr << "WARNING: Invalid critical nondimensional strain in "
                  "compression\n";
        opserr << "uniaxialMaterial Concrete07: " << tag << endln;
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[10], &r) != TCL_OK) {
        opserr << "WARNING: Invalid value for r\n";
        opserr << "uniaxialMaterial Concrete07: " << tag << endln;
      }

      // Parsing was successful, allocate the material
      theMaterial =
          new Concrete07(tag, fpc, epsc0, Ec, fpt, epst0, xcrp, xcrn, r);

    return ((BasicModelBuilder*)cd)->addTaggedObject<UniaxialMaterial>(*theMaterial);
//  return theMaterial;
}   
