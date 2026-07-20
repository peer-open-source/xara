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
static Tcl_CmdProc initializeAnalysis;
static Tcl_CmdProc resetModel;
static Tcl_CmdProc XaraCmd_constraints;
static Tcl_CmdProc XaraCmd_modalDamping;

Tcl_CmdProc TclCommand_clearAnalysis;
extern Tcl_CmdProc XaraCmd_numberer;
extern Tcl_CmdProc XaraCmd_number;

namespace OpenSees {
Tcl_CmdProc responseSpectrumAnalysis;
}

// commands/analysis/integrator.cpp
extern Tcl_CmdProc XaraCmd_integrator;

// commands/analysis/solver.cpp
extern Tcl_CmdProc XaraCmd_system;
extern Tcl_CmdProc XaraCmd_systemSize;

// commands/analysis/algorithm.cpp
extern Tcl_CmdProc XaraCmd_algorithm;
extern Tcl_CmdProc XaraCmd_numIter;
extern Tcl_CmdProc TclCommand_accelCPU;
extern Tcl_CmdProc TclCommand_totalCPU;
extern Tcl_CmdProc TclCommand_solveCPU;
extern Tcl_CmdProc TclCommand_numFact;

// from commands/analysis/ctest.cpp
extern Tcl_CmdProc specifyCTest;
extern Tcl_CmdProc getCTestNorms;
extern Tcl_CmdProc getCTestIter;
extern Tcl_CmdProc TclCommand_algorithmRecorder;

// from commands/analysis/sensitivity.cpp
extern Tcl_CmdProc TclCommand_sensitivityAlgorithm;
extern Tcl_CmdProc TclCommand_sensLambda;

struct char_cmd {
  const char* name;
  Tcl_CmdProc*  func;
} const tcl_analysis_cmds[] =  {
    {"system",              &XaraCmd_system},
    {"systemSize",          &XaraCmd_systemSize},

    {"test",                &specifyCTest},
    {"testIter",            &getCTestIter},
    {"testNorms",           &getCTestNorms},
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
    {"reset",               &resetModel},

  // From algorithm.cpp
    {"algorithm",           &XaraCmd_algorithm},
    {"numIter",             &XaraCmd_numIter},
    {"numFact",             &TclCommand_numFact},
    {"accelCPU",            &TclCommand_accelCPU},
    {"totalCPU",            &TclCommand_totalCPU},
    {"solveCPU",            &TclCommand_solveCPU},
  // recorder.cpp
    {"algorithmRecorder",   &TclCommand_algorithmRecorder},

  // sensitivity
    {"sensitivityAlgorithm", TclCommand_sensitivityAlgorithm},
    {"sensLambda",           TclCommand_sensLambda},
};

