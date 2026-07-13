#pragma once

#include <vector> 
#include <Channel.h>
#include <MPI_MachineBroker.h>
#include <TclPackageClassBroker.h>

class AnalysisChannels {
public:
  AnalysisChannels();

  int getProcessID() {return theMachine.getPID();}

  Channel** getChannels() {return theChannels;}
  int getNumChannels() {return numChannels;}

private:
  TclPackageClassBroker m_obroker;
  MPI_MachineBroker theMachine;
  Channel** theChannels;
  int numChannels;
};