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

class TclPackageClassBroker;
class Channel;

#define MODEL_CHANNELS

namespace OpenSees {
// struct Interpreter;

class ProcessContext {
public:
  ProcessContext();
  ~ProcessContext();

  int getProcessID() {return theMachine.getPID();}
  int getProcessCount() {return theMachine.getNP();}

  Channel** getChannels() const {return theChannels;}
  int getNumChannels() const {return numChannels;}

  int setup(Interpreter&);
  int clean(Interpreter&) {return 0;}

private:
  TclPackageClassBroker *m_obroker;
  MPI_MachineBroker theMachine;
  Channel** theChannels;
  int numChannels;
};

} // namespace 
