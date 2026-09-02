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
#include <Parsing.h>
//
// NOTE: wipe is not added to the table on purpose, it cannot be added to
// and removed from the interpreter as simply as the other commands.
static Tcl_CmdProc wipeAnalysis;
//
static Tcl_CmdProc XaraCmd_analysis;
static Tcl_CmdProc XaraCmd_analyze;
static Tcl_CmdProc XaraCmd_eigen;
static Tcl_CmdProc XaraCmd_printA;
static Tcl_CmdProc XaraCmd_printB;
static Tcl_CmdProc XaraCmd_applyA;
static Tcl_CmdProc XaraCmd_solveA;
static Tcl_CmdProc initializeAnalysis;
static Tcl_CmdProc resetModel;
static Tcl_CmdProc XaraCmd_constraints;
static Tcl_CmdProc XaraCmd_modalDamping;

Tcl_CmdProc TclCommand_clearAnalysis;
extern Tcl_CmdProc XaraCmd_numberer;
extern Tcl_CmdProc XaraCmd_number;

Tcl_CmdProc XaraCmd_responseSpectrumAnalysis;

// commands/analysis/integrator.cpp
extern Tcl_CmdProc XaraCmd_integrator;

// commands/analysis/solver.cpp
extern Tcl_CmdProc XaraCmd_system;
extern Tcl_CmdProc XaraCmd_systemSize;

// commands/analysis/algorithm.cpp
extern Tcl_CmdProc XaraCmd_algorithm;
extern Tcl_CmdProc XaraCmd_numIter;
extern Tcl_CmdProc XaraCmd_accelCPU;
extern Tcl_CmdProc XaraCmd_totalCPU;
extern Tcl_CmdProc XaraCmd_solveCPU;
extern Tcl_CmdProc XaraCmd_numFact;

// from commands/analysis/ctest.cpp
extern Tcl_CmdProc XaraCmd_test;
extern Tcl_CmdProc XaraCmd_testNorms;
extern Tcl_CmdProc XaraCmd_testIter;
extern Tcl_CmdProc XaraCmd_algorithmRecorder;

// from commands/analysis/sensitivity.cpp
extern Tcl_CmdProc TclCommand_sensitivityAlgorithm;
extern Tcl_CmdProc TclCommand_sensLambda;

Tcl_CmdProc XaraCmd_printEigenMatrices;

struct char_cmd {
  const char* name;
  Tcl_CmdProc*  func;
} const tcl_analysis_cmds[] =  {
    // {"printEigenMatrices", &XaraCmd_printEigenMatrices},

    {"system",              &XaraCmd_system},
    {"systemSize",          &XaraCmd_systemSize},

    {"test",                &XaraCmd_test},
    {"testIter",            &XaraCmd_testIter},
    {"testNorms",           &XaraCmd_testNorms},
    {"integrator",          &XaraCmd_integrator},
    {"constraints",         &XaraCmd_constraints},

    {"eigen",               &XaraCmd_eigen},
    {"analysis",            &XaraCmd_analysis},

    {"analyze",             &XaraCmd_analyze},
    {"initialize",          &initializeAnalysis},
    {"modalDamping",        &XaraCmd_modalDamping},
    {"modalDampingQ",       &XaraCmd_modalDamping},
    {"printA",              &XaraCmd_printA},
    {"printB",              &XaraCmd_printB},
    // for testing 
    {"applyA",              &XaraCmd_applyA},
    {"solveA",              &XaraCmd_solveA},

    {"reset",               &resetModel},

  // From algorithm.cpp
    {"algorithm",           &XaraCmd_algorithm},
    {"numIter",             &XaraCmd_numIter},
    {"numFact",             &XaraCmd_numFact},
    {"accelCPU",            &XaraCmd_accelCPU},
    {"totalCPU",            &XaraCmd_totalCPU},
    {"solveCPU",            &XaraCmd_solveCPU},
  // recorder.cpp
    {"algorithmRecorder",   &XaraCmd_algorithmRecorder},

  // sensitivity
    {"sensitivityAlgorithm", TclCommand_sensitivityAlgorithm},
    {"sensLambda",           TclCommand_sensLambda},
};

