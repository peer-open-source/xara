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
// Description: This file contains functions that are responsible
// for orchestrating an analysis.
//
#include <tcl.h>
#include <assert.h>
#include <Parsing.h>
#include <runtimeAPI.h>
#include <Logging.h>
#include <StandardStream.h>
#include <FileStream.h>

#include <Matrix.h>
#include <Domain.h> // for modal damping
#include <AnalysisModel.h>
#include <ModelRegistry.h>

#include "BasicAnalysisBuilder.h"

#include <EigenSOE.h>
#include <LinearSOE.h>
// for printA
#include <FullGenLinLapackSolver.h>
#include <FullGenLinSOE.h>
#include <GimmeMCK.h>

#include <LoadControl.h>
#include <EquiSolnAlgo.h>

#include <TransientIntegrator.h>
#include <StaticIntegrator.h>

// constraint handlers
#include <PlainHandler.h>
#include <AutoConstraintHandler.h>
#include <PenaltyConstraintHandler.h>
#include <LagrangeConstraintHandler.h>
#include <TransformationConstraintHandler.h>

// numberers
#include <numberer/DOF_Numberer.h>
#include <numberer/PlainNumberer.h>
#include "analysis.h"


// extern int OPS_ResponseSpectrumAnalysis(G3_Runtime*);
extern "C" int OPS_ResetInputNoBuilder(ClientData clientData,
                                       Tcl_Interp *interp, int cArg, int mArg,
                                       TCL_Char ** const argv, Domain *domain);

Tcl_CmdProc TclCommand_clearAnalysis;
Tcl_CmdProc TclCommand_setNumberer;
namespace OpenSees {
Tcl_CmdProc responseSpectrumAnalysis;
}


//
// Add commands to the interpreter that take the AnalysisBuilder as clientData.
//
int
G3_AddTclAnalysisAPI(Tcl_Interp *interp, ModelRegistry& context)
{
  BasicAnalysisBuilder *analysis = new BasicAnalysisBuilder(context);
  Tcl_CreateCommand(interp, "wipeAnalysis", &wipeAnalysis, analysis, nullptr);
  Tcl_CreateCommand(interp, "_clearAnalysis", &TclCommand_clearAnalysis, analysis, nullptr);

  Tcl_CreateCommand(interp, "numberer",   TclCommand_setNumberer, analysis, nullptr);

  Tcl_CreateCommand(interp, "responseSpectrumAnalysis", &OpenSees::responseSpectrumAnalysis, nullptr, nullptr);


  static int ncmd = sizeof(tcl_analysis_cmds)/sizeof(char_cmd);
  for (int i = 0; i < ncmd; ++i)
    Tcl_CreateCommand(interp, 
        tcl_analysis_cmds[i].name, 
        tcl_analysis_cmds[i].func, 
        (ClientData)analysis, nullptr);

  return TCL_OK;
}

//
// command invoked to build an Analysis object
//
static int
specifyAnalysis(ClientData clientData, Tcl_Interp *interp, Tcl_Size argc,
                TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder*)clientData;

  if (argc < 2) {
    opserr << OpenSees::PromptValueError 
           << "need to specify an analysis type (Static, Transient)\n";
    return TCL_ERROR;
  }

  int argi = 1;

  if (strcmp(argv[argi], "-linear") == 0) {
    if (argc < 3) {
      opserr << OpenSees::PromptValueError << "need to specify an analysis type (Static, Transient)\n";
      return TCL_ERROR;
    }
    Tcl_Eval(interp, "algorithm Linear\n"
                     "test FixedNumIter 1\n"
    );
    argi++;
  }
  if (strcmp(argv[argi], "Static") == 0) {
    builder->setStaticAnalysis();
    return TCL_OK;

  } else if (strcmp(argv[argi], "Transient") == 0) {
    builder->setTransientAnalysis();
    return TCL_OK;
  }

  else if (((strcmp(argv[1], "VariableTimeStepTransient") == 0) ||
          (strcmp(argv[1], "TransientWithVariableTimeStep") == 0) ||
          (strcmp(argv[1], "VariableTransient") == 0))) {
    opserr << "Unimplemented\n";
    return TCL_ERROR;

  }
  else {
    opserr << OpenSees::PromptValueError << "Analysis type '" << argv[1]
      << "' does not exists (Static or Transient only). \n";
    return TCL_ERROR;
  }

  return TCL_OK;
}

//
// Command invoked to build the model, i.e. to invoke analyze()
// on the Analysis object
//
static int
analyzeModel(ClientData clientData, Tcl_Interp *interp, Tcl_Size argc,
             TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder*)clientData;

  int result = 0;
  int commit = BasicAnalysisBuilder::Increment
             | BasicAnalysisBuilder::Iterate 
             | BasicAnalysisBuilder::Commit;
  
  for (int i=2; i<argc; i++) {
    if (strcmp(argv[i], "-operation") == 0) {
      if (argc < i+2) {
        opserr << OpenSees::PromptValueError << "operation key requires argument\n";
        return TCL_ERROR;
      }
      i++;
      if (strcmp(argv[i], "commit") == 0) {
        commit = BasicAnalysisBuilder::Commit;
      }
      else if (strcmp(argv[i], "increment") == 0) {
        commit = BasicAnalysisBuilder::Increment;
      }
      else if (strcmp(argv[i], "iteration") == 0) {
        commit = BasicAnalysisBuilder::Iterate;
      }
    }
  }

  switch (builder->CurrentAnalysisFlag) {
    case BasicAnalysisBuilder::STATIC_ANALYSIS: {
      int numIncr;
      if (argc < 2) {
        opserr << OpenSees::PromptValueError << "static analysis: analysis numIncr?\n";
        return TCL_ERROR;
      }

      if (Tcl_GetInt(interp, argv[1], &numIncr) != TCL_OK)
        return TCL_ERROR;

      result = builder->analyze(numIncr, 0.0, commit);
      break;
    }
    case BasicAnalysisBuilder::TRANSIENT_ANALYSIS: {
      double dT;
      int numIncr;
      if (argc < 3) {
        opserr << OpenSees::PromptValueError << "transient analysis: analysis numIncr? deltaT?\n";
        return TCL_ERROR;
      }
      if (Tcl_GetInt(interp, argv[1], &numIncr) != TCL_OK)
        return TCL_ERROR;
      if (Tcl_GetDouble(interp, argv[2], &dT) != TCL_OK)
        return TCL_ERROR;

      if (argc == 6) {
        int Jd;
        double dtMin, dtMax;
        if (Tcl_GetDouble(interp, argv[3], &dtMin) != TCL_OK)
          return TCL_ERROR;
        if (Tcl_GetDouble(interp, argv[4], &dtMax) != TCL_OK)
          return TCL_ERROR;
        if (Tcl_GetInt(interp, argv[5], &Jd) != TCL_OK)
          return TCL_ERROR;

        result = builder->analyzeVariable(numIncr, dT, dtMin, dtMax, Jd);

      } else {
        result = builder->analyze(numIncr, dT, commit);
      }
      break;
    }
    default:
      opserr << OpenSees::PromptValueError << "No Analysis type has been specified \n";
      return TCL_ERROR;
  }

  Tcl_SetObjResult(interp, Tcl_NewIntObj(result));

  return TCL_OK;
}


static int
initializeAnalysis(ClientData clientData, Tcl_Interp *interp, Tcl_Size argc,
                   TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder*)clientData;
  if (builder->initialize() < 0)
    return TCL_ERROR;
  return TCL_OK;
}



static int
eigenAnalysis(ClientData clientData,
              Tcl_Interp *interp, 
              Tcl_Size argc,
              TCL_Char ** const argv)
{

  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder*)clientData;

  Domain *domain = builder->getDomain();

  // make sure at least one other argument to contain type of system
  if (argc < 2) {
    opserr << OpenSees::PromptValueError
           << "eigen <type> numModes?\n";
    return TCL_ERROR;
  }

  bool generalizedAlgo = true; 
      // 0 - frequency/generalized (default),
      // 1 - standard,
      // 2 - buckling

  int typeSolver = EigenSOE_TAGS_ArpackSOE;
  double shift = 0.0;
  bool findSmallest = true;
  int  numEigen = -1;

  // Check type of eigenvalue analysis
  int loc = 1;
  while (loc < argc) {
    if ((strcmp(argv[loc], "frequency") == 0) ||
        (strcmp(argv[loc], "-frequency") == 0) ||
        (strcmp(argv[loc], "generalized") == 0) ||
        (strcmp(argv[loc], "-generalized") == 0))
      generalizedAlgo = true;

    else if ((strcmp(argv[loc], "standard") == 0) ||
             (strcmp(argv[loc], "-standard") == 0)) {
      generalizedAlgo = false;
      typeSolver = EigenSOE_TAGS_SymBandEigenSOE;
    }

    else if ((strcmp(argv[loc], "-findLargest") == 0))
      findSmallest = false;

    else if ((strcmp(argv[loc], "genBandArpack") == 0) ||
             (strcmp(argv[loc], "-genBandArpack") == 0) ||
             (strcmp(argv[loc], "genBandArpackEigen") == 0) ||
             (strcmp(argv[loc], "-genBandArpackEigen") == 0))
      typeSolver = EigenSOE_TAGS_ArpackSOE;

    else if ((strcmp(argv[loc], "symmBandLapack") == 0) ||
             (strcmp(argv[loc], "-symmBandLapack") == 0) ||
             (strcmp(argv[loc], "symmBandLapackEigen") == 0) ||
             (strcmp(argv[loc], "-symmBandLapackEigen") == 0))
      typeSolver = EigenSOE_TAGS_SymBandEigenSOE;

    else if ((strcmp(argv[loc], "fullGenLapack") == 0) ||
             (strcmp(argv[loc], "-fullGenLapack") == 0) ||
             (strcmp(argv[loc], "fullGenLapackEigen") == 0) ||
             (strcmp(argv[loc], "-fullGenLapackEigen") == 0))
      typeSolver = EigenSOE_TAGS_FullGenEigenSOE;

    else if (numEigen == -1) {
      if ((Tcl_GetInt(interp, argv[loc], &numEigen) != TCL_OK) || (numEigen < 0)) {
        opserr << OpenSees::PromptValueError
              << "invalid number of modes " << argv[loc]
              << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }
    }

    else {
      opserr << "Unknown option: " << argv[loc] 
             << OpenSees::SignalMessageEnd;
    }

    loc++;
  }

  if (numEigen < 0) {
    opserr << OpenSees::PromptValueError
           << "eigen command requires number of modes to be specified"
           << OpenSees::SignalMessageEnd;
  }

  //
  // create a transient analysis if no analysis exists
  // 
  builder->newEigenAnalysis(typeSolver, shift);

  int result = builder->eigen(numEigen,generalizedAlgo,findSmallest);


  if (result == 0) {
    Tcl_Obj* eig_values = Tcl_NewListObj(numEigen, nullptr);
    const Vector &eigenvalues = domain->getEigenvalues();
    for (int i = 0; i < numEigen; ++i) {
      Tcl_ListObjAppendElement(interp, eig_values, Tcl_NewDoubleObj(eigenvalues[i]));
    }
    Tcl_SetObjResult(interp, eig_values);
    return TCL_OK;
  }
  else {
    opserr << OpenSees::PromptValueError
           << "Failed to perform eigenvalue analysis"
           << OpenSees::SignalMessageEnd;
    
    return TCL_ERROR;
  }
}



// TODO: Move this to commands/modeling/damping.cpp? ...but it uses and
// AnalysisBuilder
static int
modalDamping(ClientData clientData, Tcl_Interp *interp, Tcl_Size argc,
             TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder*)clientData;
  int numEigen = builder->getNumEigen();

  if (argc < 2) {
    opserr
        << OpenSees::PromptValueError << argv[0] << " ?factor - not enough arguments to command\n";
    return TCL_ERROR;
  }

  int numModes = argc - 1;

  if (numEigen == 0) {
    opserr << G3_WARN_PROMPT 
           << "- " << argv[0] << " - eigen command should be called first\n";

    numEigen = numModes;
    builder->newEigenAnalysis(EigenSOE_TAGS_ArpackSOE, 0.0);
    builder->eigen(numModes, true, true);
    // return TCL_ERROR;
  }

  /* 
   * "quick" modal damping adds modal damping forces to the right-hand side,
   * but does not add modal damping terms to the dynamic tangent.
   *
   * see https://portwooddigital.com/2022/11/08/quick-and-dirty-modal-damping/
   */
  bool do_tangent = true;
  if (strcmp(argv[0], "modalDampingQ") == 0)
    do_tangent = false;

  double factor = 0;
  Vector modalDampingValues(numEigen);

  if (numModes != 1 && numModes != numEigen) {
    // TODO: Just call eigen again?
    opserr << OpenSees::PromptValueError 
           << "same number of damping factors as modes must be "
              "specified\n";
//  opserr << " same damping ratio will be applied to all\n";
    return TCL_ERROR;
  }

  //
  // read in values and set factors
  //
  if (numModes == numEigen) {

    // read in all factors one at a time
    for (int i = 0; i < numEigen; ++i) {
      if (Tcl_GetDouble(interp, argv[1 + i], &factor) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "could not read factor at position "
               << i << "\n";
        return TCL_ERROR;
      }
      modalDampingValues[i] = factor;
    }

  } else {
    //  read in one & set all factors to that value
    if (Tcl_GetDouble(interp, argv[1], &factor) != TCL_OK) {
      opserr << OpenSees::PromptValueError 
             << "rayleigh alphaM? betaK? betaK0? betaKc? - could not "
                "read betaK? \n";
      return TCL_ERROR;
    }

    for (int i = 0; i < numEigen; ++i)
      modalDampingValues[i] = factor;
  }

  // set factors in domain
  Domain *theDomain = builder->getDomain();
  assert(theDomain != nullptr);

  theDomain->setModalDampingFactors(&modalDampingValues, do_tangent);
  return TCL_OK;
}

static int
resetModel(ClientData clientData, Tcl_Interp *interp, Tcl_Size argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder*)clientData;

  Domain* domain = builder->getDomain();
  assert(domain != nullptr);
  domain->revertToStart();

  TransientIntegrator *theTransientIntegrator =  builder->getTransientIntegrator();
  if (theTransientIntegrator != nullptr) {
    theTransientIntegrator->revertToStart();
  }
  return TCL_OK;
}


int
printIntegrator(ClientData clientData, Tcl_Interp *interp, Tcl_Size argc,
                TCL_Char ** const argv, OPS_Stream &output)
{
  assert(clientData != nullptr);
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder*)clientData;

  TransientIntegrator *theTransientIntegrator =  builder->getTransientIntegrator();
  StaticIntegrator *the_static_integrator = builder->getStaticIntegrator();

  int eleArg = 0;
  if (the_static_integrator == nullptr && theTransientIntegrator == nullptr)
    return TCL_OK;

  // if 'print <filename> Algorithm flag' get the flag
  int flag = 0;
  if (argc > 2) {
    if (Tcl_GetInt(interp, argv[eleArg], &flag) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "print algorithm failed to get integer flag: \n";
      opserr << argv[eleArg] << "\n";
      return TCL_ERROR;
    }
  }

  if (the_static_integrator != nullptr)
    the_static_integrator->Print(output, flag);
  else
    theTransientIntegrator->Print(output, flag);
  return TCL_OK;
}

static int
printA(ClientData clientData, Tcl_Interp *interp, Tcl_Size argc, TCL_Char ** const argv)
{
  // printA <filename> - m <double> -c <double> -k <double>
  assert(clientData != nullptr);
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder*)clientData;

  int res = 0;

  enum class Format {
    None
  } format = Format::None;

  FileStream outputFile;
  OPS_Stream *output = &opserr;
  LinearSOE  *oldSOE = builder->getLinearSOE();


  // Cant allocate the Solver on stack because theSOE is going to 
  // delete it
  FullGenLinSOE theSOE(*new FullGenLinLapackSolver());

  builder->set(&theSOE, false);
  // Construct a graph and pass
  // it to the SOE. Otherwise, getA() returns null
  builder->domainChanged();


  bool ret = false;
  double m = 0.0, c = 0.0, k = 0.0;
  bool do_mck = false;
  int currentArg = 1;
  while (currentArg < argc) {
    if ((strcmp(argv[currentArg], "file") == 0) ||
        (strcmp(argv[currentArg], "-file") == 0)) {
      currentArg++;
      if (currentArg == argc) {
        opserr << G3_WARN_PROMPT << "-file missing argument\n";
        return TCL_ERROR;
      }

      if (outputFile.setFile(argv[currentArg]) != 0) {
        opserr << OpenSees::PromptValueError
               << "failed to open file: "
               << argv[currentArg] << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }
      output = &outputFile;
    }
    else if ((strcmp(argv[currentArg], "ret") == 0) ||
             (strcmp(argv[currentArg], "-ret") == 0)) {
      ret = true;
    }

    else if ((strcmp(argv[currentArg], "-m") == 0)) {
      currentArg++;
      if (Tcl_GetDouble(interp, argv[currentArg], &m) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to read float following flag -m\n";
        return TCL_ERROR;
      }
      do_mck = true;
    }

    else if ((strcmp(argv[currentArg], "-c") == 0)) {
      currentArg++;
      if (Tcl_GetDouble(interp, argv[currentArg], &c) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to read float following flag -c\n";
        return TCL_ERROR;
      }
      do_mck = true;
    }

    else if ((strcmp(argv[currentArg], "-k") == 0)) {
      currentArg++;
      if (Tcl_GetDouble(interp, argv[currentArg], &k) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to read float following flag -k\n";
        return TCL_ERROR;
      }
      do_mck = true;
    }
    currentArg++;
  }
  
  //
  // Form the tangent
  //
  TransientIntegrator *oldint = nullptr;

  // construct integrator here so that it is not
  // destructed when the `if` scope ends
  GimmeMCK integrator(m, c, k, 0.0);
  if (do_mck) {
    oldint = builder->getTransientIntegrator();
    builder->set(integrator, false);
    integrator.formTangent(0);
    integrator.revertToLastStep();
  }

  else if (builder->getStaticIntegrator() != nullptr) {
    builder->getStaticIntegrator()->formTangent();
    builder->getStaticIntegrator()->revertToLastStep();
  }

  else if (builder->getTransientIntegrator() != nullptr) {
    builder->getTransientIntegrator()->formTangent(0);
    builder->getTransientIntegrator()->revertToLastStep();
  }
  else {
    opserr << OpenSees::PromptValueError 
           << "No integrator has been set, cannot form tangent"
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  const Matrix *A = theSOE.getA();
  if (A == nullptr) {
    opserr << OpenSees::PromptValueError 
           << "Could not get matrix from linear system\n";
    return TCL_ERROR;
  }

  if (format == Format::None) {
    if (ret) {
      int n = A->noRows();
      int m = A->noCols();
      if (n*m == 0) {
        opserr << OpenSees::PromptValueError 
               << "linear system is empty, got n=" << n << " m= " << m << "\n";
        return TCL_ERROR;
      }

      // Create an empty list with space preallocated for
      // n*m elements. This is not formally documented, but
      // it is mentioned here 
      //   https://wiki.tcl-lang.org/page/Tcl_NewListObj
      //
      // and evident from the source code here:
      //   https://github.com/enthought/tcl/blob/master/generic/tclListObj.c
      //
      Tcl_Obj* list = Tcl_NewListObj(n*m, nullptr);


      for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; j++)
          Tcl_ListObjAppendElement(interp, list, Tcl_NewDoubleObj((*A)(i, j)));
      }
      Tcl_SetObjResult(interp, list);
    }
    else {
      *output << *A;
      outputFile.close();
    }
  }

  builder->getDomain()->revertToLastCommit();
  // put the original SOE back.
  if (oldSOE != nullptr)
    builder->set(oldSOE, true);

  if (oldint != nullptr)
    builder->set(*oldint, true);

  return res;
}

static int
printB(ClientData clientData, Tcl_Interp *interp, Tcl_Size argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder*)clientData;

  int res = 0;

  FileStream outputFile;
  OPS_Stream *output = &opserr;

  bool ret = false;
  int currentArg = 1;
  while (currentArg < argc) {
    if ((strcmp(argv[currentArg], "file") == 0) ||
        (strcmp(argv[currentArg], "-file") == 0)) {
      currentArg++;
      if (currentArg == argc) {
        opserr << G3_WARN_PROMPT << "-file missing argument\n";
        return TCL_ERROR;
      }

      if (outputFile.setFile(argv[currentArg]) != 0) {
        opserr << "print <filename> .. - failed to open file: "
               << argv[currentArg] << "\n";
        return TCL_ERROR;
      }
      output = &outputFile;
    } else if ((strcmp(argv[currentArg], "ret") == 0) ||
               (strcmp(argv[currentArg], "-ret") == 0)) {
      ret = true;
    }
    currentArg++;
  }

  LinearSOE  *theSOE = builder->getLinearSOE();
  if (theSOE != nullptr) {

    // TODO
    builder->formUnbalance();

    if (theSOE->getNumEqn() == 0) {
      opserr << OpenSees::PromptValueError << "System of equations is empty\n";
      return TCL_ERROR;
    }

    const Vector &b = theSOE->getB();

    if (ret) {
      const int size = b.Size();
      Tcl_Obj* list = Tcl_NewListObj(size, nullptr);
      for (int i = 0; i < size; ++i)
        Tcl_ListObjAppendElement(interp, list, Tcl_NewDoubleObj(b[i]));

      Tcl_SetObjResult(interp, list);
    } else {
      *output << b;
      outputFile.close();
    }
  }

  return res;
}

// This is removed from the Tcl_Interp in model.cpp
extern int
TclCommand_clearAnalysis(ClientData cd, Tcl_Interp *interp, Tcl_Size argc, TCL_Char ** const argv)
{

  if (cd != nullptr) {
    BasicAnalysisBuilder *builder = (BasicAnalysisBuilder *)cd;
    builder->wipe();
    delete builder;

    static int ncmd = sizeof(tcl_analysis_cmds)/sizeof(char_cmd);
    for (int i = 0; i < ncmd; ++i)
      Tcl_DeleteCommand(interp, tcl_analysis_cmds[i].name);

    Tcl_CreateCommand(interp, "wipeAnalysis",  &wipeAnalysis, nullptr, nullptr);
    Tcl_CreateCommand(interp, "_clearAnalysis", &TclCommand_clearAnalysis, nullptr, nullptr);
  }

  return TCL_OK;
}

static int
wipeAnalysis(ClientData cd, Tcl_Interp *interp, Tcl_Size argc, TCL_Char ** const argv)
{
  if (cd != nullptr) {
    BasicAnalysisBuilder *builder = (BasicAnalysisBuilder *)cd;
    builder->wipe();
  }
  return TCL_OK;
}


//
// command invoked to allow the ConstraintHandler object to be built
//
static int
specifyConstraintHandler(ClientData clientData, Tcl_Interp *interp, Tcl_Size argc,
                         TCL_Char ** const argv)
{
  
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder*)clientData;

  // make sure at least one other argument to contain type name
  if (argc < 2) {
    opserr << OpenSees::PromptValueError << "need to specify a constraint type \n";
    return TCL_ERROR;
  }

  //
  // Create the handler
  //

  ConstraintHandler *theHandler = nullptr;
  // check argv[1] for type of handler and create the object
  if (strcmp(argv[1], "Plain") == 0)
    theHandler = new PlainHandler();

  else if (strcmp(argv[1], "Transformation") == 0) {
    theHandler = new TransformationConstraintHandler();
  }

  else if (strcmp(argv[1], "Penalty") == 0) {
    // handler Penalty alpha1 alpha2
    if (argc < 4) {
      opserr << OpenSees::PromptValueError
             << "need to specify alpha\n";
      return TCL_ERROR;
    }
    double alpha1, alpha2;
    if (Tcl_GetDouble(interp, argv[2], &alpha1) != TCL_OK)
      return TCL_ERROR;
    if (Tcl_GetDouble(interp, argv[3], &alpha2) != TCL_OK)
      return TCL_ERROR;
    theHandler = new PenaltyConstraintHandler(alpha1, alpha2);
  }

  else if (strcmp(argv[1], "Lagrange") == 0) {
    double alpha1 = 1.0;
    double alpha2 = 1.0;
    if (argc == 4) {
      if (Tcl_GetDouble(interp, argv[2], &alpha1) != TCL_OK)
        return TCL_ERROR;
      if (Tcl_GetDouble(interp, argv[3], &alpha2) != TCL_OK)
        return TCL_ERROR;
    }
    theHandler = new LagrangeConstraintHandler(alpha1, alpha2);
  }

  else if (strcmp(argv[1], "Auto") == 0) {
    bool verbose = false;
    bool auto_penalty = true;
    double auto_penalty_oom = 3.0;
    double user_penalty = 0.0;
    bool auto_penalty_done = false;
    bool user_penalty_done = false;
    for (int i=2; i<argc; i++) {
      if (strcmp(argv[i], "-verbose") == 0)
        verbose = true;
      else if (strcmp(argv[i], "-autoPenalty") == 0) {
        if (user_penalty_done) {
          opserr << OpenSees::PromptValueError << "cannot use with userPenalty\n";
          return TCL_ERROR;
        }
        if (argc < i+1) {
          opserr << OpenSees::PromptValueError << "autoPenalty needs a value\n";
          return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[i+1], &auto_penalty_oom) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "autoPenalty needs a numeric value\n";
          return TCL_ERROR;
        }
        i++;
        auto_penalty = true;
        auto_penalty_done = true;
      }
      else if (strcmp(argv[i], "-userPenalty") == 0) {
        if (argc < i+1) {
          opserr << OpenSees::PromptValueError << "userPenalty needs a value\n";
          return TCL_ERROR;
        }
        if (auto_penalty_done) {
          opserr << OpenSees::PromptValueError << "cannot use userPenalty with autoPenalty\n";
          return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[i+1], &user_penalty) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "userPenalty needs a numeric value\n";
          return TCL_ERROR;
        }
        i++;
        auto_penalty = false;
        user_penalty_done = true;
      }
      else {
        opserr << OpenSees::PromptValueError << "unknown option " << argv[i] << "\n";
        return TCL_ERROR;
      }
    }
    theHandler = new AutoConstraintHandler(verbose, auto_penalty, auto_penalty_oom, user_penalty);
  }

  else {
    opserr << OpenSees::PromptValueError << "ConstraintHandler type '" << argv[1]
      << "' does not exists \n\t(Plain, Penalty, Lagrange, Transformation, Auto) only\n";
    return TCL_ERROR;
  }

  builder->set(theHandler);
  return TCL_OK;
}
