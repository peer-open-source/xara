#include "AnalysisChannels.h"

AnalysisChannels::AnalysisChannels()
: m_obroker()
, theMachine(&m_obroker, 0, nullptr)
{
  int OPS_rank = theMachine.getPID();
  int OPS_np = theMachine.getNP();

  if (OPS_rank == 0) {
    theChannels = new Channel *[OPS_np-1];
    numChannels = OPS_np-1;


    for (int j=0; j<OPS_np-1; j++) {
      Channel *otherChannel = theMachine.getRemoteProcess();
      theChannels[j] = otherChannel;
    //   otherChannel->sendID(0,0,data);
    //   otherChannel->sendMsg(0,0,msgChar);
    }
  }
  else {
    theChannels = new Channel *[1];
    numChannels = 1;

    Channel *myChannel = theMachine.getMyChannel();
    theChannels[0] = myChannel;
  }
}