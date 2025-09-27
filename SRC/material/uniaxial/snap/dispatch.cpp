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

// $Revision: 1.5 $
// $Date: 2005-01-14 23:45:58 $
// $Source:
// /usr/local/cvs/OpenSees/SRC/material/uniaxial/snap/TclSnapMaterialCommand.cpp,v
// $

// Written: Arash Altoontash, Gregory Deierlein,
// Created: Feb 2002
// Modified: Arash June 2004
//
// Description: This file contains the implementation of the
// TclBasicBuilder_addSnapMaterial() function.
#include <Logging.h>
#include <Parsing.h>
#include <Pinching.h>
#include <Clough.h>
#include <CloughHenry.h>
#include <PinchingDamage.h>
#include <CloughDamage.h>
#include <Bilinear.h>
#include <DamageModel.h>
#include <ModelRegistry.h>

#include <Vector.h>
#include <string.h>


UniaxialMaterial *
TclBasicBuilder_addSnapMaterial(ClientData clientData, Tcl_Interp *interp,
                                int argc, TCL_Char ** const argv)
{

  ModelRegistry* builder = (ModelRegistry*)clientData;
  
  if (argc < 3) {
    opserr << "WARNING insufficient number of arguments for the Snap material "
              "model" << OpenSees::SignalMessageEnd;
    return 0;
  }

  int tag;
  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << "WARNING invalid uniaxialMaterial tag" << OpenSees::SignalMessageEnd;
    return 0;
  }

  UniaxialMaterial *theMaterial = 0;

  if (strcmp(argv[1], "Bilinear") == 0) {
    if (argc < 15) {
      opserr << "WARNING insufficient arguments" << OpenSees::SignalMessageEnd;
      return 0;
    }

    Vector input(12);
    double temp;

    for (int i = 3, j = 0; j < 12; i++, j++) {
      if (Tcl_GetDouble(interp, argv[i], &temp) != TCL_OK) {
        opserr << "WARNING invalid input, data " << i << OpenSees::SignalMessageEnd;
        return 0;
      }
      input(j) = temp;
    }

    DamageModel *strength;
    if ((int)input(9) == 0) {
      strength = NULL;
    } else {
      strength = builder->getTypedObject<DamageModel>((int)input(9));
      if (strength == nullptr) {
        return nullptr;
      }
    }

    DamageModel *stiffness;
    if ((int)input(10) == 0) {
      stiffness = NULL;
    } else {
      stiffness = builder->getTypedObject<DamageModel>((int)input(10));
      if (stiffness == 0) {
        return nullptr;
      }
    }

    DamageModel *capping;
    if ((int)input(11) == 0) {
      capping = NULL;
    } else {
      capping = builder->getTypedObject<DamageModel>((int)input(11));
      if (capping == 0) {
        return nullptr;
      }
    }

    theMaterial = new Bilinear(tag, input, strength, stiffness, capping);
  }

  else if ((strcmp(argv[1], "Clough") == 0) ||
           (strcmp(argv[1], "clough") == 0) ||
           (strcmp(argv[1], "CloughHenry") == 0)) {

    if (argc < 19) {
      opserr << "WARNING insufficient arguments" << OpenSees::SignalMessageEnd;
      return 0;
    }

    Vector input(16);
    double temp;

    for (int i = 3, j = 0; j < 16; i++, j++) {
      if (Tcl_GetDouble(interp, argv[i], &temp) != TCL_OK) {
        opserr << "WARNING invalid input, data " << i << OpenSees::SignalMessageEnd;
        return 0;
      }
      input(j) = temp;
    }

    if ((strcmp(argv[1], "Clough") == 0) || (strcmp(argv[1], "clough") == 0))
      theMaterial = new Clough(tag, input);
    else
      theMaterial = new CloughHenry(tag, input);
  }

  else if (strcmp(argv[1], "Clough_Damage") == 0 ||
           strcmp(argv[1], "CloughDamage") == 0 ||
           strcmp(argv[1], "Clough_Damage") == 0 ||
           strcmp(argv[1], "CloughDamage") == 0) {
    if (argc < 15) {
      opserr << "WARNING insufficient arguments"
             << OpenSees::SignalMessageEnd;
      return 0;
    }

    Vector input(12);
    double temp;

    for (int i = 3, j = 0; j < 12; i++, j++) {
      if (Tcl_GetDouble(interp, argv[i], &temp) != TCL_OK) {
        opserr << "WARNING invalid input, data " << i << OpenSees::SignalMessageEnd;
        return 0;
      }
      input(j) = temp;
    }

    DamageModel *strength;
    if ((int)input(8) == 0) {
      strength = NULL;
    } else {
      strength = builder->getTypedObject<DamageModel>((int)input(8));

      if (strength == 0) {
        opserr << OpenSees::PromptModelError
               << "damage model for strength deterioration not found"
               << OpenSees::SignalMessageEnd;
        return 0;
      }
    }

    DamageModel *stiffness;
    if ((int)input(9) == 0) {
      stiffness = NULL;
    } else {
      stiffness = builder->getTypedObject<DamageModel>((int)input(9));

      if (stiffness == 0) {
        opserr << OpenSees::PromptModelError
            << "damage model for stiffness deterioration not found"
            << OpenSees::SignalMessageEnd;
        return 0;
      }
    }

    DamageModel *accelerated;
    if ((int)input(10) == 0) {
      accelerated = NULL;
    } else {
      accelerated = builder->getTypedObject<DamageModel>((int)input(10));

      if (accelerated == 0) {
        opserr << "WARNING damage model for accelerated stiffness "
                  "deterioration not found"
                << OpenSees::SignalMessageEnd;
        return 0;
      }
    }

    DamageModel *capping;
    if ((int)input(11) == 0) {
      capping = NULL;
    } else {
      capping = builder->getTypedObject<DamageModel>((int)input(11));
      if (capping == 0) {
        return 0;
      }
    }
    theMaterial =
        new CloughDamage(tag, input, strength, stiffness, accelerated, capping);
  }

  else if (strcmp(argv[1], "Pinching") == 0 ||
           strcmp(argv[1], "pinching") == 0) {
    if (argc < 22) {
      opserr << "WARNING insufficient arguments\n";
      return 0;
    }

    Vector input(19);
    double temp;

    for (int i = 3, j = 0; j < 19; i++, j++) {
      if (Tcl_GetDouble(interp, argv[i], &temp) != TCL_OK) {
        opserr << "WARNING invalid input, data " << i << OpenSees::SignalMessageEnd;
        return 0;
      }
      input(j) = temp;
    }
    theMaterial = new Pinching(tag, input);
  }

  else if (strcmp(argv[1], "Pinching_Damage") == 0 ||
           strcmp(argv[1], "pinching_Damage") == 0 ||
           strcmp(argv[1], "PinchingDamage") == 0 ||
           strcmp(argv[1], "pinchingDamage") == 0) {
    if (argc < 18) {
      opserr << "WARNING insufficient arguments" << OpenSees::SignalMessageEnd;
      return 0;
    }

    Vector input(15);
    double temp;

    for (int i = 3, j = 0; j < 15; i++, j++) {
      if (Tcl_GetDouble(interp, argv[i], &temp) != TCL_OK) {
        opserr << "WARNING invalid input, data " << i << OpenSees::SignalMessageEnd;
        return 0;
      }
      input(j) = temp;
    }

    DamageModel *strength;
    if ((int)input(11) == 0) {
      strength = NULL;
    } else {
      strength = builder->getTypedObject<DamageModel>((int)input(11));

      if (strength == 0) {
        opserr << "WARNING damage model for strength deterioration not found" << OpenSees::SignalMessageEnd;
        return 0;
      }
    }

    DamageModel *stiffness;
    if ((int)input(12) == 0) {
      stiffness = NULL;
    } else {
      stiffness = builder->getTypedObject<DamageModel>((int)input(12));
      if (stiffness == 0) {
        opserr
            << "WARNING damage model for stiffness deterioration not found" << OpenSees::SignalMessageEnd;
        return 0;
      }
    }

    DamageModel *accelerated;
    if ((int)input(13) == 0) {
      accelerated = NULL;
    } else {
      accelerated = builder->getTypedObject<DamageModel>((int)input(13));
      if (accelerated == 0) {
        opserr << "WARNING damage model for accelerated stiffness "
                  "deterioration not found" << OpenSees::SignalMessageEnd;
        return 0;
      }
    }

    DamageModel *capping;
    if ((int)input(14) == 0) {
      capping = NULL;
    } else {
      capping = builder->getTypedObject<DamageModel>((int)input(14));
      if (capping == nullptr) {
        return 0;
      }
    }
    theMaterial = new PinchingDamage(tag, input, strength, stiffness,
                                     accelerated, capping);
  }

  return theMaterial;
}
