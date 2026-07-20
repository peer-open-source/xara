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
#pragma once

#include <vector>
#include <runtime/interpreter/Interpreter.h>
#include <MPI_MachineBroker.h>

class XaraClassBroker;
class Channel;

#define MODEL_CHANNELS

namespace Xara {
// struct Interpreter;

class ProcessContext {
public:
  ProcessContext();
  ~ProcessContext();
#if defined(XARA_ENABLE_MPI)
  int getProcessID() {return theMachine.getPID();}
  int getProcessCount() {return theMachine.getNP();}
#else 
  int getProcessID() {return 0;}
  int getProcessCount() {return 1;}
#endif

  Channel** getChannels() const {return theChannels;}
  int getNumChannels() const {return numChannels;}

  int setup(Interpreter&);
  int clean(Interpreter&) {return 0;}

private:
  XaraClassBroker *m_obroker;
#if defined(XARA_ENABLE_MPI)
  MPI_MachineBroker theMachine;
#endif
  Channel** theChannels;
  int numChannels;
};

} // namespace 
