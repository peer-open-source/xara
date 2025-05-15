/* ****************************************************************** **
**    Opensees - Open System for Earthquake Engineering Simulation    **
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
// Written: NM (nmitra@u.washington.edu)
// Created: April 2002
// Revised: August 2004
//
// Description: This file contains the implementation of the commands used
// to add beam column joint to a model.
// Update: Optional User interfaces added in.
//
#include <tcl.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <Domain.h>


#include <BasicModelBuilder.h>
#include <BeamColumnJoint2d.h>
#include <BeamColumnJoint3d.h>
#include <Information.h>
#include <ElementResponse.h>

#include <Node.h>
#include <UniaxialMaterial.h>


int
TclBasicBuilder_addBeamColumnJoint(ClientData clientData, Tcl_Interp *interp,
                                   int argc, TCL_Char ** const argv)
{
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;
  Domain* domain = builder->getDomain();

  constexpr static int eleArgStart = 1;

  int NDM = builder->getNDM(); // dimension of the structure (1d, 2d, or 3d)
  int NDF = builder->getNDF(); // number of degrees of freedom per node

  if ((NDM == 2 && NDF == 3) || (NDM == 3 && NDF == 6)) {

    // check no of arguments
    if ((argc - eleArgStart) != 19 && (argc - eleArgStart) != 21) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: element beamColumnJoint eleTag? node1? node2? node3? "
                "node4? matTag1? matTag2? matTag3?\n";
      opserr << "matTag4? matTag5? matTag6? matTag7? matTag8? matTag9? "
                "matTag10? matTag11? matTag12? matTag13?\n";
      opserr << "<ElementHeightFactor? ElementWidthFactor?>\n";

      return TCL_ERROR;
    }

    int id, nd1, nd2, nd3, nd4, matId1, matId2, matId3, matId4, matId5, matId6,
        matId7, matId8, matId9, matId10;
    int matId11, matId12, matId13;
    double hgtfac, wdtfac;

    UniaxialMaterial *theMaterial1 = 0;
    UniaxialMaterial *theMaterial2 = 0;
    UniaxialMaterial *theMaterial3 = 0;
    UniaxialMaterial *theMaterial4 = 0;
    UniaxialMaterial *theMaterial5 = 0;
    UniaxialMaterial *theMaterial6 = 0;
    UniaxialMaterial *theMaterial7 = 0;
    UniaxialMaterial *theMaterial8 = 0;
    UniaxialMaterial *theMaterial9 = 0;
    UniaxialMaterial *theMaterial10 = 0;
    UniaxialMaterial *theMaterial11 = 0;
    UniaxialMaterial *theMaterial12 = 0;
    UniaxialMaterial *theMaterial13 = 0;

    if (Tcl_GetInt(interp, argv[1 + eleArgStart], &id) != TCL_OK) {
      opserr << "WARNING invalid beamColumnJoint eleTag" << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[2 + eleArgStart], &nd1) != TCL_OK) {
      opserr << "WARNING invalid Node 1\n";
      opserr << "beamColumnJoint Element: " << id << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[3 + eleArgStart], &nd2) != TCL_OK) {
      opserr << "WARNING invalid Node 2\n";
      opserr << "beamColumnJoint Element: " << id << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[4 + eleArgStart], &nd3) != TCL_OK) {
      opserr << "WARNING invalid Node 3\n";
      opserr << "beamColumnJoint Element: " << id << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[5 + eleArgStart], &nd4) != TCL_OK) {
      opserr << "WARNING invalid Node 4\n";
      opserr << "beamColumnJoint Element: " << id << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[6 + eleArgStart], &matId1) != TCL_OK) {
      opserr << "WARNING invalid Material Tag 1\n";
      opserr << "beamColumnJoint Element: " << id << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[7 + eleArgStart], &matId2) != TCL_OK) {
      opserr << "WARNING invalid Material Tag 2\n";
      opserr << "beamColumnJoint Element: " << id << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[8 + eleArgStart], &matId3) != TCL_OK) {
      opserr << "WARNING invalid Material Tag 3\n";
      opserr << "beamColumnJoint Element: " << id << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[9 + eleArgStart], &matId4) != TCL_OK) {
      opserr << "WARNING invalid Material Tag 4\n";
      opserr << "beamColumnJoint Element: " << id << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[10 + eleArgStart], &matId5) != TCL_OK) {
      opserr << "WARNING invalid Material Tag 5\n";
      opserr << "beamColumnJoint Element: " << id << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[11 + eleArgStart], &matId6) != TCL_OK) {
      opserr << "WARNING invalid Material Tag 6\n";
      opserr << "beamColumnJoint Element: " << id << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[12 + eleArgStart], &matId7) != TCL_OK) {
      opserr << "WARNING invalid Material Tag 7\n";
      opserr << "beamColumnJoint Element: " << id << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[13 + eleArgStart], &matId8) != TCL_OK) {
      opserr << "WARNING invalid Material Tag 8\n";
      opserr << "beamColumnJoint Element: " << id << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[14 + eleArgStart], &matId9) != TCL_OK) {
      opserr << "WARNING invalid Material Tag 9\n";
      opserr << "beamColumnJoint Element: " << id << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[15 + eleArgStart], &matId10) != TCL_OK) {
      opserr << "WARNING invalid Material Tag 10\n";
      opserr << "beamColumnJoint Element: " << id << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[16 + eleArgStart], &matId11) != TCL_OK) {
      opserr << "WARNING invalid Material Tag 11\n";
      opserr << "beamColumnJoint Element: " << id << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[17 + eleArgStart], &matId12) != TCL_OK) {
      opserr << "WARNING invalid Material Tag 12\n";
      opserr << "beamColumnJoint Element: " << id << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[18 + eleArgStart], &matId13) != TCL_OK) {
      opserr << "WARNING invalid Material Tag 13\n";
      opserr << "beamColumnJoint Element: " << id << endln;
      return TCL_ERROR;
    }

    if ((argc - eleArgStart) == 21) {
      if (Tcl_GetDouble(interp, argv[19 + eleArgStart], &hgtfac) != TCL_OK) {
        opserr << "WARNING invalid factor for height\n";
        opserr << "beamColumnJoint Element: " << id << endln;
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[20 + eleArgStart], &wdtfac) != TCL_OK) {
        opserr << "WARNING invalid factor for width\n";
        opserr << "beamColumnJoint Element: " << id << endln;
        return TCL_ERROR;
      }
    }

    if (matId1 != 0) {
      theMaterial1 = builder->getTypedObject<UniaxialMaterial>(matId1);

      if (theMaterial1 == 0) {
        opserr << "WARNING material not found\n";
        opserr << "Material: " << matId1;
        opserr << "\nbeamColumnJoint element: " << id << endln;
        return TCL_ERROR;
      }
    } else
      theMaterial1 = 0;

    if (matId2 != 0) {
      theMaterial2 = builder->getTypedObject<UniaxialMaterial>(matId2);

      if (theMaterial2 == 0) {
        opserr << "WARNING material not found\n";
        opserr << "Material: " << matId2;
        opserr << "\nbeamColumnJoint element: " << id << endln;
        return TCL_ERROR;
      }
    } else
      theMaterial2 = 0;

    if (matId3 != 0) {
      theMaterial3 = builder->getTypedObject<UniaxialMaterial>(matId3);

      if (theMaterial3 == 0) {
        opserr << "WARNING material not found\n";
        opserr << "Material: " << matId3;
        opserr << "\nbeamColumnJoint element: " << id << endln;
        return TCL_ERROR;
      }
    } else
      theMaterial3 = 0;

    if (matId4 != 0) {
      theMaterial4 = builder->getTypedObject<UniaxialMaterial>(matId4);

      if (theMaterial4 == 0) {
        opserr << "WARNING material not found\n";
        opserr << "Material: " << matId4;
        opserr << "\nbeamColumnJoint element: " << id << endln;
        return TCL_ERROR;
      }
    } else
      theMaterial4 = 0;

    if (matId5 != 0) {
      theMaterial5 = builder->getTypedObject<UniaxialMaterial>(matId5);

      if (theMaterial5 == 0) {
        opserr << "WARNING material not found\n";
        opserr << "Material: " << matId5;
        opserr << "\nbeamColumnJoint element: " << id << endln;
        return TCL_ERROR;
      }
    } else
      theMaterial5 = 0;

    if (matId6 != 0) {
      theMaterial6 = builder->getTypedObject<UniaxialMaterial>(matId6);

      if (theMaterial6 == 0) {
        opserr << "WARNING material not found\n";
        opserr << "Material: " << matId6;
        opserr << "\nbeamColumnJoint element: " << id << endln;
        return TCL_ERROR;
      }
    } else
      theMaterial6 = 0;

    if (matId7 != 0) {
      theMaterial7 = builder->getTypedObject<UniaxialMaterial>(matId7);

      if (theMaterial7 == 0) {
        opserr << "WARNING material not found\n";
        opserr << "Material: " << matId7;
        opserr << "\nbeamColumnJoint element: " << id << endln;
        return TCL_ERROR;
      }
    } else
      theMaterial7 = 0;

    if (matId8 != 0) {
      theMaterial8 = builder->getTypedObject<UniaxialMaterial>(matId8);

      if (theMaterial8 == 0) {
        opserr << "WARNING material not found\n";
        opserr << "Material: " << matId8;
        opserr << "\nbeamColumnJoint element: " << id << endln;
        return TCL_ERROR;
      }
    } else
      theMaterial8 = 0;

    if (matId9 != 0) {
      theMaterial9 = builder->getTypedObject<UniaxialMaterial>(matId9);

      if (theMaterial9 == 0) {
        opserr << "WARNING material not found\n";
        opserr << "Material: " << matId9;
        opserr << "\nbeamColumnJoint element: " << id << endln;
        return TCL_ERROR;
      }
    } else
      theMaterial9 = 0;

    if (matId10 != 0) {
      theMaterial10 = builder->getTypedObject<UniaxialMaterial>(matId10);

      if (theMaterial10 == 0) {
        opserr << "WARNING material not found\n";
        opserr << "Material: " << matId10;
        opserr << "\nbeamColumnJoint element: " << id << endln;
        return TCL_ERROR;
      }
    } else
      theMaterial10 = 0;

    if (matId11 != 0) {
      theMaterial11 = builder->getTypedObject<UniaxialMaterial>(matId11);

      if (theMaterial11 == 0) {
        opserr << "WARNING material not found\n";
        opserr << "Material: " << matId11;
        opserr << "\nbeamColumnJoint element: " << id << endln;
        return TCL_ERROR;
      }
    } else
      theMaterial11 = 0;

    if (matId12 != 0) {
      theMaterial12 = builder->getTypedObject<UniaxialMaterial>(matId12);

      if (theMaterial12 == 0) {
        opserr << "WARNING material not found\n";
        opserr << "Material: " << matId12;
        opserr << "\nbeamColumnJoint element: " << id << endln;
        return TCL_ERROR;
      }
    } else
      theMaterial12 = 0;

    if (matId13 != 0) {
      theMaterial13 = builder->getTypedObject<UniaxialMaterial>(matId13);

      if (theMaterial13 == 0) {
        opserr << "WARNING material not found\n";
        opserr << "Material: " << matId13;
        opserr << "\nbeamColumnJoint element: " << id << endln;
        return TCL_ERROR;
      }
    } else
      theMaterial13 = 0;

    // create Beam Column Joint Element and add it to Domain
    Element *theBeamColumnJoint = 0;

    if (NDM == 2) {
      if ((argc - eleArgStart) == 19) {
        theBeamColumnJoint = new BeamColumnJoint2d(
            id, nd1, nd2, nd3, nd4, *theMaterial1, *theMaterial2, *theMaterial3,
            *theMaterial4, *theMaterial5, *theMaterial6, *theMaterial7,
            *theMaterial8, *theMaterial9, *theMaterial10, *theMaterial11,
            *theMaterial12, *theMaterial13);
      } else if ((argc - eleArgStart) == 21) {
        theBeamColumnJoint = new BeamColumnJoint2d(
            id, nd1, nd2, nd3, nd4, *theMaterial1, *theMaterial2, *theMaterial3,
            *theMaterial4, *theMaterial5, *theMaterial6, *theMaterial7,
            *theMaterial8, *theMaterial9, *theMaterial10, *theMaterial11,
            *theMaterial12, *theMaterial13, hgtfac, wdtfac);
      }
    } else if (NDM == 3) {
      if ((argc - eleArgStart) == 19) {
        theBeamColumnJoint = new BeamColumnJoint3d(
            id, nd1, nd2, nd3, nd4, *theMaterial1, *theMaterial2, *theMaterial3,
            *theMaterial4, *theMaterial5, *theMaterial6, *theMaterial7,
            *theMaterial8, *theMaterial9, *theMaterial10, *theMaterial11,
            *theMaterial12, *theMaterial13);
      } else if ((argc - eleArgStart) == 21) {
        theBeamColumnJoint = new BeamColumnJoint3d(
            id, nd1, nd2, nd3, nd4, *theMaterial1, *theMaterial2, *theMaterial3,
            *theMaterial4, *theMaterial5, *theMaterial6, *theMaterial7,
            *theMaterial8, *theMaterial9, *theMaterial10, *theMaterial11,
            *theMaterial12, *theMaterial13, hgtfac, wdtfac);
      }
    }

    if (theBeamColumnJoint == 0) {
      opserr << "WARNING ran out of memory creating elements\n";
      opserr << "beamColumnJoint element: " << id << endln;
      return TCL_ERROR;
    }

    if (domain->addElement(theBeamColumnJoint) == false) {
      opserr << "WARNING could not add element to the domain\n";
      opserr << "beamColumnJoint element: " << id << endln;
      delete theBeamColumnJoint;
      return TCL_ERROR;
    }

  } else {
    opserr << "WARNING NDM = " << NDM << " and NDF = " << NDF
           << " is incompatible with available joint elements";
    return TCL_ERROR;
  }

  // the element successfully created and added to the domain
  return TCL_OK;
}
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
// Written: Arash Altoontash, Gregory Deierlein	Created: 04/01
// Revision:
//				AA		02/03
//
// Description: This file contains the implementation of the
// TclBasicBuilder_addJoint2D() command.
//
#include <vector>
#include <stdlib.h>
#include <string.h>
#include <Domain.h>

#include <Joint2D.h>
#include <tcl.h>
#include <DamageModel.h>
#include <UniaxialMaterial.h>


int
TclBasicBuilder_addJoint2D(ClientData clientData, Tcl_Interp *interp, int argc,
                           TCL_Char ** const argv)
{
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;
  Domain* domain = builder->getDomain();
  // check the number of arguments is correct
  int argStart = 2;

  if ((argc - argStart) != 8 && (argc - argStart) != 10 &&
      (argc - argStart) != 12 && (argc - argStart) != 18) {
    opserr << "WARNING incorrect number of arguments\n";
    opserr << "Want:\n";
    opserr
        << "element Joint2D Tag? NodI? NodJ? NodK? NodL? NodC? MatC? LrgDsp?\n";
    opserr << "or:\n";
    opserr << "element Joint2D Tag? NodI? NodJ? NodK? NodL? NodC? MatC? "
              "LrgDsp? -damage DmgTag?\n";
    opserr << "or:\n";
    opserr << "element Joint2D Tag? NodI? NodJ? NodK? NodL? NodC? MatI? MatJ? "
              "MatK? MatL? MatC? LrgDsp?\n";
    opserr << "or:\n";
    opserr << "element Joint2D Tag? NodI? NodJ? NodK? NodL? NodC? MatI? MatJ? "
              "MatK? MatL? MatC? LrgDsp? -damage DmgI DmgJ DmgK DmgL DmgC\n";
    return TCL_ERROR;
  }

  // get the id and end nodes
  int Joint2DId, iNode, jNode, kNode, lNode;
  if (Tcl_GetInt(interp, argv[argStart], &Joint2DId) != TCL_OK) {
    opserr << "WARNING invalid Joint2D eleTag" << endln;
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[1 + argStart], &iNode) != TCL_OK) {
    opserr << "WARNING invalid iNode\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[2 + argStart], &jNode) != TCL_OK) {
    opserr << "WARNING invalid jNode\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[3 + argStart], &kNode) != TCL_OK) {
    opserr << "WARNING invalid kNode\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[4 + argStart], &lNode) != TCL_OK) {
    opserr << "WARNING invalid lNode\n";
    return TCL_ERROR;
  }

  // Get the center node
  int CenterNodeTag;
  if (Tcl_GetInt(interp, argv[5 + argStart], &CenterNodeTag) != TCL_OK) {
    opserr << "WARNING invalid tag for center node\n";
    return TCL_ERROR;
  }

  // check domain for existence of internal node tag
  Node *CenterNode = domain->getNode(CenterNodeTag);
  if (CenterNode != 0) {
    opserr
        << "WARNING node tag specified for the center node already exists.\n";
    opserr << "Use a new node tag.\n";
    return TCL_ERROR;
  }

  UniaxialMaterial *MatI = nullptr;
  UniaxialMaterial *MatJ = nullptr;
  UniaxialMaterial *MatK = nullptr;
  UniaxialMaterial *MatL = nullptr;
  UniaxialMaterial *PanelMaterial = nullptr;
  Joint2D *theJoint2D;
  int LargeDisp;

  // Decide to use which constructor, based on the number of arguments
  if ((argc - argStart) == 8 || (argc - argStart) == 12) {

    // Using Joint2D constructor without damage

    if ((argc - argStart) == 8) {
      int PanelMatId;
      if (Tcl_GetInt(interp, argv[6 + argStart], &PanelMatId) != TCL_OK) {
        opserr << "WARNING invalid matID\n";
        return TCL_ERROR;
      }

      if (Tcl_GetInt(interp, argv[7 + argStart], &LargeDisp) != TCL_OK) {
        // use 0 as default
        LargeDisp = 0;
      }

      PanelMaterial = builder->getTypedObject<UniaxialMaterial>(PanelMatId);
      if (PanelMaterial == nullptr)
        return TCL_ERROR;

    }

    else // if ( (argc-argStart) == 12  )
    {
      int MatIid;
      if (Tcl_GetInt(interp, argv[6 + argStart], &MatIid) != TCL_OK) {
        opserr << "WARNING invalid material ID for spring I\n";
        return TCL_ERROR;
      }

      if (MatIid != 0) {
        MatI = builder->getTypedObject<UniaxialMaterial>(MatIid);
        if (MatI == nullptr)
          return TCL_ERROR;

      } else
        MatI = nullptr;

      int MatJid;
      if (Tcl_GetInt(interp, argv[7 + argStart], &MatJid) != TCL_OK) {
        opserr << "WARNING invalid material ID for spring J\n";
        return TCL_ERROR;
      }

      if (MatJid != 0) {
        MatJ = builder->getTypedObject<UniaxialMaterial>(MatJid);
        if (MatJ == nullptr)
          return TCL_ERROR;

      } else
        MatJ = nullptr;

      int MatKid;
      if (Tcl_GetInt(interp, argv[8 + argStart], &MatKid) != TCL_OK) {
        opserr << "WARNING invalid material ID for spring K\n";

        return TCL_ERROR;
      }
      if (MatKid != 0) {
        MatK = builder->getTypedObject<UniaxialMaterial>(MatKid);
        if (MatK == nullptr)
          return TCL_ERROR;

      } else
        MatK = nullptr;

      int MatLid;
      if (Tcl_GetInt(interp, argv[9 + argStart], &MatLid) != TCL_OK) {
        opserr << "WARNING invalid material ID for spring L\n";
        return TCL_ERROR;
      }
      if (MatLid != 0) {
        MatL = builder->getTypedObject<UniaxialMaterial>(MatLid);

        if (MatL == nullptr) {
          opserr << "WARNING material not found\n";
          opserr << "Material: " << MatLid;
          opserr << "\nJoint2D element: " << Joint2DId << endln;
          return TCL_ERROR;
        }
      } else
        MatL = nullptr;

      int PanelMatId;
      if (Tcl_GetInt(interp, argv[10 + argStart], &PanelMatId) != TCL_OK) {
        opserr << "WARNING invalid matID\n";
        return TCL_ERROR;
      }
      PanelMaterial = builder->getTypedObject<UniaxialMaterial>(PanelMatId);

      if (PanelMaterial == 0) {
        opserr << "WARNING material not found\n";
        opserr << "Material: " << PanelMatId;
        return TCL_ERROR;
      }

      if (Tcl_GetInt(interp, argv[11 + argStart], &LargeDisp) != TCL_OK) {
        // use 0 as default
        LargeDisp = 0;
      }
    }

    UniaxialMaterial *springModels[5] = {MatI, MatJ, MatK, MatL, PanelMaterial};
    theJoint2D =
        new Joint2D(Joint2DId, iNode, jNode, kNode, lNode, CenterNodeTag,
                    springModels, domain, LargeDisp);


    if (domain->addElement(theJoint2D) == false) {
      opserr << "WARNING could not add element to the domain\n";
      opserr << "Joint2D element: " << Joint2DId << endln;
      delete theJoint2D;
      return TCL_ERROR;
    }

    // if get here we have successfully created the element and added it to the
    // domain
    return TCL_OK;
  }

  else if ((argc - argStart) == 10 || (argc - argStart) == 18) {
    // Using Joint2D constructor with damage
    DamageModel *DmgI = nullptr;
    DamageModel *DmgJ = nullptr;
    DamageModel *DmgK = nullptr;
    DamageModel *DmgL = nullptr;
    DamageModel *PanelDamage = nullptr;

    if ((argc - argStart) == 10) {
      int PanelMatId;
      if (Tcl_GetInt(interp, argv[6 + argStart], &PanelMatId) != TCL_OK) {
        opserr << "WARNING invalid matID\n";
        return TCL_ERROR;
      }

      if (Tcl_GetInt(interp, argv[7 + argStart], &LargeDisp) != TCL_OK) {
        // use 0 as default
        LargeDisp = 0;
      }

      PanelMaterial = builder->getTypedObject<UniaxialMaterial>(PanelMatId);

      if (PanelMaterial == 0) {
        opserr << "WARNING material not found\n";
        opserr << "Material: " << PanelMatId;
        return TCL_ERROR;
      }

      if (strcmp(argv[8 + argStart], "-damage") != 0 &&
          strcmp(argv[8 + argStart], "-Damage") != 0) {
        opserr << "WARNING incorrect command line\n";
        return TCL_ERROR;
      }

      int PanelDamageId;
      if (Tcl_GetInt(interp, argv[9 + argStart], &PanelDamageId) != TCL_OK) {
        opserr << "WARNING invalid damageID\n";
        return TCL_ERROR;
      }

      DamageModel *PanelDamage = builder->getTypedObject<DamageModel>(PanelDamageId);

      if (PanelDamage == 0) {
        opserr << "WARNING damage model not found\n";
        return TCL_ERROR;
      }
    }

    else // if ( (argc-argStart) == 18  )
    {
      int MatIid;
      if (Tcl_GetInt(interp, argv[6 + argStart], &MatIid) != TCL_OK) {
        opserr << "WARNING invalid material ID for spring I\n";
        return TCL_ERROR;
      }

      if (MatIid != 0) {
        MatI = builder->getTypedObject<UniaxialMaterial>(MatIid);

        if (MatI == nullptr) {
          opserr << "WARNING material not found\n";
          return TCL_ERROR;
        }
      } else
        MatI = nullptr;

      int MatJid;
      if (Tcl_GetInt(interp, argv[7 + argStart], &MatJid) != TCL_OK) {
        opserr << "WARNING invalid material ID for spring J\n";
        return TCL_ERROR;
      }

      if (MatJid != 0) {
        MatJ = builder->getTypedObject<UniaxialMaterial>(MatJid);

        if (MatJ == nullptr) {
          opserr << "WARNING material not found\n";
          opserr << "Material: " << MatJid;
          return TCL_ERROR;
        }
      } else
        MatJ = nullptr;

      int MatKid;
      if (Tcl_GetInt(interp, argv[8 + argStart], &MatKid) != TCL_OK) {
        opserr << "WARNING invalid material ID for spring K\n";

        return TCL_ERROR;
      }
      if (MatKid != 0) {
        MatK = builder->getTypedObject<UniaxialMaterial>(MatKid);

        if (MatK == nullptr) {
          opserr << "WARNING material not found\n";
          opserr << "Material: " << MatKid;
          opserr << "\nJoint2D element: " << Joint2DId << endln;
          return TCL_ERROR;
        }
      } else
        MatK = nullptr;

      int MatLid;
      if (Tcl_GetInt(interp, argv[9 + argStart], &MatLid) != TCL_OK) {
        opserr << "WARNING invalid material ID for spring L\n";
        return TCL_ERROR;
      }
      if (MatLid != 0) {
        MatL = builder->getTypedObject<UniaxialMaterial>(MatLid);

        if (MatL == nullptr) {
          opserr << "WARNING material not found\n";
          opserr << "Material: " << MatLid;
          opserr << "\nJoint2D element: " << Joint2DId << endln;
          return TCL_ERROR;
        }
      } else
        MatL = nullptr;

      int PanelMatId;
      if (Tcl_GetInt(interp, argv[10 + argStart], &PanelMatId) != TCL_OK) {
        opserr << "WARNING invalid matID\n";
        return TCL_ERROR;
      }
      PanelMaterial = builder->getTypedObject<UniaxialMaterial>(PanelMatId);

      if (PanelMaterial == 0) {
        opserr << "WARNING material not found\n";
        opserr << "Material: " << PanelMatId;
        return TCL_ERROR;
      }

      if (Tcl_GetInt(interp, argv[11 + argStart], &LargeDisp) != TCL_OK) {
        // use 0 as default
        LargeDisp = 0;
      }

      if (strcmp(argv[12 + argStart], "-damage") != 0 &&
          strcmp(argv[12 + argStart], "-Damage") != 0) {
        opserr << "WARNING incorrect command line\n";
        return TCL_ERROR;
      }

      int DmgIid;
      if (Tcl_GetInt(interp, argv[13 + argStart], &DmgIid) != TCL_OK) {
        opserr << "WARNING invalid damage model ID for spring I\n";
        return TCL_ERROR;
      }

      if (DmgIid != 0 && MatI != 0) {
        DmgI = builder->getTypedObject<DamageModel>(DmgIid);

        if (DmgI == nullptr) {
          opserr << "WARNING damage model not found\n";
          opserr << "Damage Model: " << DmgIid;
          opserr << "\nJoint2D element: " << Joint2DId << endln;
          return TCL_ERROR;
        }
      } else
        DmgI = nullptr;

      int DmgJid;
      if (Tcl_GetInt(interp, argv[14 + argStart], &DmgJid) != TCL_OK) {
        opserr << "WARNING invalid damage model ID for spring J\n";
        return TCL_ERROR;
      }

      if (DmgJid != 0 && MatJ != 0) {
        DmgJ = builder->getTypedObject<DamageModel>(DmgJid);

        if (DmgJ == nullptr) {
          opserr << "WARNING damage model not found\n";
          opserr << "Damage Model: " << DmgJid;
          opserr << "\nJoint2D element: " << Joint2DId << endln;
          return TCL_ERROR;
        }
      } else
        DmgJ = nullptr;

      int DmgKid;
      if (Tcl_GetInt(interp, argv[15 + argStart], &DmgKid) != TCL_OK) {
        opserr << "WARNING invalid damage model ID for spring K\n";
        return TCL_ERROR;
      }

      if (DmgKid != 0 && MatK != 0) {
        DmgK = builder->getTypedObject<DamageModel>(DmgKid);

        if (DmgK == nullptr) {
          opserr << "WARNING damage model not found\n";
          opserr << "Damage Model: " << DmgKid;
          opserr << "\nJoint2D element: " << Joint2DId << endln;
          return TCL_ERROR;
        }
      } else
        DmgK = nullptr;

      int DmgLid;
      if (Tcl_GetInt(interp, argv[16 + argStart], &DmgLid) != TCL_OK) {
        opserr << "WARNING invalid damage model ID for spring L\n";
        return TCL_ERROR;
      }

      if (DmgLid != 0 && MatL != 0) {
        DmgL = builder->getTypedObject<DamageModel>(DmgLid);

        if (DmgL == nullptr) {
          opserr << "WARNING damage model not found\n";
          opserr << "Damage Model: " << DmgLid;
          opserr << "\nJoint2D element: " << Joint2DId << endln;
          return TCL_ERROR;
        }
      } else
        DmgL = nullptr;

      int PanelDmgId;
      if (Tcl_GetInt(interp, argv[17 + argStart], &PanelDmgId) != TCL_OK) {
        opserr << "WARNING invalid panel DmgID\n";
        return TCL_ERROR;
      }

      if (PanelDmgId != 0 && PanelMaterial != 0) {
        PanelDamage = builder->getTypedObject<DamageModel>(PanelDmgId);

        if (PanelDamage == nullptr) {
          opserr << "WARNING damage model not found\n";
          opserr << "Damage Model: " << PanelDmgId;
          opserr << "\nJoint2D element: " << Joint2DId << endln;
          return TCL_ERROR;
        }
      } else
        DmgL = nullptr;
    }

    // Create the new material
    DamageModel *damageModels[5] = {DmgI, DmgJ, DmgK, DmgL, PanelDamage};
    UniaxialMaterial *springModels[5] = {MatI, MatJ, MatK, MatL, PanelMaterial};
    theJoint2D =
        new Joint2D(Joint2DId, iNode, jNode, kNode, lNode, CenterNodeTag,
                    springModels, domain, LargeDisp, damageModels);


    if (domain->addElement(theJoint2D) == false) {
      opserr << "WARNING could not add element to the domain\n";
      opserr << "Joint2D element: " << Joint2DId << endln;
      delete theJoint2D;
      return TCL_ERROR;
    }

    return TCL_OK;

  } else {
    return TCL_ERROR;
  }
}
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

// $Revision: 1.3 $
// $Date: 2004-09-01 04:01:27 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/joint/TclJoint3dCommand.cpp,v $

// Written: Arash Altoontash, Gregory Deierlein,	Created: 04/01
// Revision:
//				AAA		02/03
//
// Description: This file contains the implementation of the
// TclBasicBuilder_addJoint3D() command.
//
#include <Joint3D.h>
#include <UniaxialMaterial.h>


int
TclBasicBuilder_addJoint3D(ClientData clientData, Tcl_Interp *interp, int argc,
                           TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder *builder = (BasicModelBuilder*)clientData;
  Domain* domain = builder->getDomain();


  if (builder->getNDM() != 3 || builder->getNDF() != 6) {
    opserr << "WARNING -- model dimensions and/or nodal DOF not compatible "
              "with Joint3D element\n";
    return TCL_ERROR;
  }

  // check the number of arguments is correct
  int argStart = 2;

  if ((argc - argStart) != 12 && (argc - argStart) != 16) {
    opserr << "WARNING incorrect number of arguments\n";
    opserr << "Want:\n";
    opserr << "element Joint3D Tag? NodI? NodJ? NodK? NodL? NodM? NodN? NodC? "
              "MatX? MatY? MatZ? LrgDsp?\n";
    opserr << "or:\n";
    opserr << "element Joint3D Tag? NodI? NodJ? NodK? NodL? NodM? NodN? NodC? "
              "MatX? MatY? MatZ? LrgDsp? -damage DmgX DmgY DmgZ\n";
    return TCL_ERROR;
  }

  // get the id and end nodes
  int Joint3DId, iNode, jNode, kNode, lNode, mNode, nNode;
  if (Tcl_GetInt(interp, argv[argStart], &Joint3DId) != TCL_OK) {
    opserr << "WARNING invalid Joint3D eleTag" << endln;
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[1 + argStart], &iNode) != TCL_OK) {
    opserr << "WARNING invalid iNode\n";
    opserr << "Joint3D element: " << Joint3DId << endln;
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[2 + argStart], &jNode) != TCL_OK) {
    opserr << "WARNING invalid jNode\n";
    opserr << "Joint3D element: " << Joint3DId << endln;
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[3 + argStart], &kNode) != TCL_OK) {
    opserr << "WARNING invalid kNode\n";
    opserr << "Joint3D element: " << Joint3DId << endln;
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[4 + argStart], &lNode) != TCL_OK) {
    opserr << "WARNING invalid lNode\n";
    opserr << "Joint3D element: " << Joint3DId << endln;
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[5 + argStart], &mNode) != TCL_OK) {
    opserr << "WARNING invalid mNode\n";
    opserr << "Joint3D element: " << Joint3DId << endln;
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[6 + argStart], &nNode) != TCL_OK) {
    opserr << "WARNING invalid nNode\n";
    opserr << "Joint3D element: " << Joint3DId << endln;
    return TCL_ERROR;
  }
  // Get the center node
  int CenterNodeTag;
  if (Tcl_GetInt(interp, argv[7 + argStart], &CenterNodeTag) != TCL_OK) {
    opserr << "WARNING invalid tag for center node\n";
    opserr << "Joint3D element: " << Joint3DId << endln;
    return TCL_ERROR;
  }

  // check domain for existence of internal node tag
  Node *CenterNode = domain->getNode(CenterNodeTag);
  if (CenterNode != 0) {
    opserr
        << "WARNING node tag specified for the center node already exists.\n";
    opserr << "Use a new node tag.\n";
    opserr << "Joint3D element: " << Joint3DId << endln;
    return TCL_ERROR;
  }

  UniaxialMaterial *MatX = NULL;
  int MatXid;
  if (Tcl_GetInt(interp, argv[8 + argStart], &MatXid) != TCL_OK) {
    opserr << "WARNING invalid material ID for spring X\n";
    opserr << "Joint3D element: " << Joint3DId << endln;
    return TCL_ERROR;
  }

  MatX = builder->getTypedObject<UniaxialMaterial>(MatXid);
  if (MatX == NULL) {
    opserr << "WARNING material not found\n";
    opserr << "Material: " << MatXid;
    opserr << "\nJoint3D element: " << Joint3DId << endln;
    return TCL_ERROR;
  }

  UniaxialMaterial *MatY = NULL;
  int MatYid;
  if (Tcl_GetInt(interp, argv[9 + argStart], &MatYid) != TCL_OK) {
    opserr << "WARNING invalid material ID for spring Y\n";
    opserr << "Joint3D element: " << Joint3DId << endln;
    return TCL_ERROR;
  }

  MatY = builder->getTypedObject<UniaxialMaterial>(MatYid);
  if (MatY == NULL) {
    opserr << "WARNING material not found\n";
    opserr << "Material: " << MatYid;
    opserr << "\nJoint3D element: " << Joint3DId << endln;
    return TCL_ERROR;
  }

  UniaxialMaterial *MatZ = NULL;
  int MatZid;
  if (Tcl_GetInt(interp, argv[10 + argStart], &MatZid) != TCL_OK) {
    opserr << "WARNING invalid material ID for spring Z\n";
    opserr << "Joint3D element: " << Joint3DId << endln;
    return TCL_ERROR;
  }

  MatZ = builder->getTypedObject<UniaxialMaterial>(MatZid);
  if (MatZ == NULL) {
    opserr << "WARNING material not found\n";
    opserr << "Material: " << MatZid;
    opserr << "\nJoint3D element: " << Joint3DId << endln;
    return TCL_ERROR;
  }

  int LargeDisp;
  if (Tcl_GetInt(interp, argv[11 + argStart], &LargeDisp) != TCL_OK) {
    // use 0 as default
    LargeDisp = 0;
  }

  Joint3D *theJoint3D;
  // Decide to use which constructor, based on the number of arguments
  if ((argc - argStart) == 12) {

    // Using Joint3D constructor without damage
    UniaxialMaterial *springModels[3] = {MatX, MatY, MatZ};
    theJoint3D =
        new Joint3D(Joint3DId, iNode, jNode, kNode, lNode, mNode, nNode,
                    CenterNodeTag, springModels, domain, LargeDisp);

    if (theJoint3D == 0) {
      opserr << "WARNING ran out of memory creating element\n";
      opserr << "Joint3D element: " << Joint3DId << endln;
      return TCL_ERROR;
    }

    if (domain->addElement(theJoint3D) == false) {
      opserr << "WARNING could not add element to the domain\n";
      opserr << "Joint3D element: " << Joint3DId << endln;
      delete theJoint3D;
      return TCL_ERROR;
    }

    // if get here we have successfully created the element and added it to the
    // domain
    return TCL_OK;
  }

  else // if ( (argc-argStart) == 16  )
  {
    // Using Joint3D constructor with damage
    // not implemented in this version
    return TCL_ERROR;
  }
  return TCL_ERROR;
}
