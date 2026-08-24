//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, OpenSees/Xara Developers
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
//

//
// NODES
//
// domain/node.cpp
Tcl_CmdProc XaraCmd_nodeCoord;
Tcl_CmdProc XaraCmd_nodeDOFs;
Tcl_CmdProc XaraCmd_nodeMass;
Tcl_CmdProc XaraCmd_nodePressure;
Tcl_CmdProc XaraCmd_nodeDisp;
Tcl_CmdProc XaraCmd_nodeVel;
Tcl_CmdProc XaraCmd_nodeAccel;
Tcl_CmdProc XaraCmd_nodeReaction;
Tcl_CmdProc XaraCmd_nodeUnbalance;
Tcl_CmdProc XaraCmd_nodeEigenvector;
Tcl_CmdProc XaraCmd_nodeRotation;
Tcl_CmdProc XaraCmd_nodeResponse;
// setters
Tcl_CmdProc XaraCmd_setNodeCoord;
Tcl_CmdProc XaraCmd_setNodeDisp;
Tcl_CmdProc XaraCmd_setNodeVel;
Tcl_CmdProc XaraCmd_setNodeAccel;
Tcl_CmdProc XaraCmd_setNodePressure;
// other
Tcl_CmdProc XaraCmd_nodeBounds;
Tcl_CmdProc XaraCmd_findID;
Tcl_CmdProc XaraCmd_getNodeTags;


// domain/region.cpp
Tcl_CmdProc TclCommand_addMeshRegion;

// domain/recorder.cpp
Tcl_CmdProc OPS_recorderValue;


// domain/section.cpp
Tcl_CmdProc sectionForce;
Tcl_CmdProc sectionDeformation;
Tcl_CmdProc sectionStiffness;
Tcl_CmdProc sectionFlexibility;
Tcl_CmdProc sectionLocation;
Tcl_CmdProc sectionWeight;
Tcl_CmdProc sectionTag;
Tcl_CmdProc sectionDisplacement;



namespace OpenSees {
namespace DomainCommands {
  // domain.cpp
  Tcl_ObjCmdProc removeObject;
  Tcl_ObjCmdProc fixedNodes;
  Tcl_ObjCmdProc constrainedNodes;
  Tcl_ObjCmdProc fixedDOFs;
  Tcl_ObjCmdProc constrainedDOFs;
  Tcl_ObjCmdProc domainChange;
  Tcl_CmdProc    retainedDOFs;
  Tcl_CmdProc    updateElementDomain;

  // domain/element.cpp
  Tcl_CmdProc addElementRayleigh;
  Tcl_CmdProc getEleTags;
  Tcl_CmdProc getNumElements;
  Tcl_CmdProc getEleClassTags;
  Tcl_CmdProc eleNodes;
  Tcl_CmdProc eleType;
  Tcl_CmdProc eleForce;
  Tcl_CmdProc localForce;
  Tcl_CmdProc eleDynamicalForce;
  Tcl_CmdProc eleResponse;

  // modal.cpp
  Tcl_CmdProc modalProperties;
}
}

// parameter.cpp
Tcl_CmdProc getParamTags;
Tcl_CmdProc TclCommand_parameter;
Tcl_CmdProc getParamValue;
Tcl_CmdProc TclCommand_setParameter;


//


// sensitivity.cpp
Tcl_CmdProc computeGradients;
Tcl_CmdProc sensNodeDisp;
Tcl_CmdProc sensLambda; // Abbas
Tcl_CmdProc sensNodeVel;
Tcl_CmdProc sensNodeAccel;
Tcl_CmdProc sensNodePressure;
Tcl_CmdProc sensSectionForce;
Tcl_CmdProc TclCommand_sensitivityAlgorithm;
// Tcl_CmdProc sensitivityIntegrator;


// Tcl_CmdProc startTimer;
// Tcl_CmdProc stopTimer;
Tcl_CmdProc XaraCmd_getTime;
Tcl_CmdProc XaraCmd_setTime;

Tcl_CmdProc XaraCmd_rayleigh;

Tcl_CmdProc modalDamping;

Tcl_CmdProc modalDampingQ;

Tcl_CmdProc basicDeformation;

Tcl_CmdProc basicForce;

Tcl_CmdProc basicStiffness;

// added: Chris McGann, U.Washington for initial state analysis of nDMaterials
Tcl_CmdProc InitialStateAnalysis;



Tcl_CmdProc setLoadConst;

Tcl_CmdProc setCreep;

Tcl_CmdProc XaraCmd_getLoadFactor;

Tcl_CmdProc printModelGID;

Tcl_CmdProc TclAddRecorder;

Tcl_CmdProc addAlgoRecorder;

Tcl_CmdProc addDatabase;

Tcl_CmdProc playbackRecorders;

Tcl_CmdProc playbackAlgorithmRecorders;

Tcl_CmdProc groundExcitation;


Tcl_CmdProc getEleLoadData;
Tcl_CmdProc getEleLoadClassTags;
Tcl_CmdProc getEleLoadTags;





//

Tcl_CmdProc XaraCmd_reactions;

Tcl_CmdProc retainedNodes;