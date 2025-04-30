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
// Description: This file contains the class definition for NodeRecorderRMS.
// A NodeRecorderRMS is used to record the envelop of specified dof responses 
// at a collection of nodes over an analysis. (between commitTag of 0 and
// last commitTag).
//
// Written: fmk 
//
#include <NodeRecorderRMS.h>
#include <Domain.h>
#include <Node.h>
#include <NodeData.h>
#include <NodeIter.h>
#include <Vector.h>
#include <ID.h>
#include <Matrix.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <MeshRegion.h>
#include <TimeSeries.h>

#include <StandardStream.h>
#include <DataFileStream.h>
#include <DataFileStreamAdd.h>
#include <XmlFileStream.h>
#include <BinaryFileStream.h>
#include <TCP_Stream.h>

#include <string.h>
#include <stdlib.h>
#include <math.h>


NodeRecorderRMS::NodeRecorderRMS()
:Recorder(RECORDER_TAGS_NodeRecorderRMS),
 theDofs(0), theNodalTags(0), theNodes(0),
 currentData(0), runningTotal(0), count(0),
 theDomain(0), theHandler(0),
 deltaT(0.0), relDeltaTTol(0.00001), nextTimeStampToRecord(0.0),
 initializationDone(false), 
 numValidNodes(0), addColumnInfo(0), theTimeSeries(0), timeSeriesValues(0),
 dataFlag(NodeData::Disp), dataIndex(-1)
{

}

NodeRecorderRMS::NodeRecorderRMS(const ID &dofs, 
                                 const ID *nodes, 
                                 NodeData _dataFlag,
                                 int _dataIndex,
                                 Domain &theDom,
                                 OPS_Stream &theOutputHandler,
                                 double dT, 
                                 double rTolDt,
                                 TimeSeries **theSeries)
:Recorder(RECORDER_TAGS_NodeRecorderRMS),
 theDofs(0), theNodalTags(0), theNodes(0),
 currentData(0), runningTotal(0), count(0),
 theDomain(&theDom), theHandler(&theOutputHandler),
 deltaT(dT), relDeltaTTol(rTolDt), nextTimeStampToRecord(0.0),
 initializationDone(false), numValidNodes(0), 
 addColumnInfo(0), theTimeSeries(theSeries), timeSeriesValues(0),
 dataFlag(_dataFlag), dataIndex(_dataIndex)
{
  // verify dof are valid 
  int numDOF = dofs.Size();
  theDofs = new ID(0, numDOF);

  int count = 0;
  for (int i=0; i<numDOF; i++) {
    int dof = dofs(i);
    if (dof >= 0) {
      (*theDofs)[count] = dof;
      count++;
    } else {
      opserr << "NodeRecorderRMS::NodeRecorderRMS - invalid dof  " << dof;
      opserr << " will be ignored\n";
    }
  }

  // 
  // create memory to hold nodal ID's (need parallel)
  //

  if (nodes != 0) {
    int numNode = nodes->Size();
    if (numNode != 0) {
      theNodalTags = new ID(*nodes);
      if (theNodalTags == 0 || theNodalTags->Size() != nodes->Size()) {
	opserr << "NodeRecorder::NodeRecorder - out of memory\n";
      }
    }
  } 

  if (theTimeSeries != nullptr) {
    timeSeriesValues = new double [numDOF];
    for (int i=0; i<numDOF; i++)
      timeSeriesValues[i] = 0.0;
  }

  if (dataFlag == NodeData::Reaction
   || dataFlag == NodeData::ReactionInclInertia
   || dataFlag == NodeData::ReactionInclRayleigh) {
    theHandler->setAddCommon(1);
  }
}


NodeRecorderRMS::~NodeRecorderRMS()
{
  //
  // write the data
  //
  if (theHandler != 0 && runningTotal != 0) {
    
    theHandler->tag("Data"); // Data
    
    int numResponse = runningTotal->Size();
    
    if (runningTotal != 0) {
      for (int j=0; j<numResponse; j++)
	if (count != 0) {
	  double value = (*runningTotal)(j);
	  (*runningTotal)(j) = sqrt(value/(count));
	}
    }

    theHandler->write(*runningTotal);
      
    theHandler->endTag(); // Data
  }

  //
  // clean up the memory
  //

  int numDOF = theDofs->Size();

  if (theDofs != 0)
    delete theDofs;

  if (theNodalTags != 0)
    delete theNodalTags;

  if (theHandler != 0)
    delete theHandler;

  if (currentData != 0)
    delete currentData;

  if (runningTotal != 0)
    delete runningTotal;

  if (theNodes != 0)
    delete [] theNodes;

  if (theTimeSeries != 0) {
    for (int i=0; i<numDOF; i++)
      delete theTimeSeries[i];
    delete [] theTimeSeries;
  }

  if (timeSeriesValues != 0)
    delete [] timeSeriesValues;  
}

int 
NodeRecorderRMS::record(int commitTag, double timeStamp)
{
  if (theDomain == 0 || theDofs == 0) {
    return 0;
  }


  if (theHandler == 0) {
    opserr << "NodeRecorderRMS::record() - no DataOutputHandler has been set\n";
    return -1;
  }


  if (initializationDone != true) 
    if (this->initialize() != 0) {
      opserr << "NodeRecorderRMS::record() - failed in initialize()\n";
      return -1;
    }


  int numDOF = theDofs->Size();

  // where relDeltaTTol is the maximum reliable ratio between analysis time step and deltaT
  // and provides tolerance for floating point precision (see floating-point-tolerance-for-recorder-time-step.md)
    if (deltaT == 0.0 || timeStamp - nextTimeStampToRecord >= -deltaT * relDeltaTTol) {
    
    if (deltaT != 0.0) 
      nextTimeStampToRecord = timeStamp + deltaT;

    double timeSeriesTerm = 0.0;

    if (theTimeSeries != 0) {
      for (int i=0; i<numDOF; i++) {
	if (theTimeSeries[i] != 0) 
	  timeSeriesValues[i] = theTimeSeries[i]->getFactor(timeStamp);
      }
    }

    //
    // if need nodal reactions get the domain to calculate them
    // before we iterate over the nodes
    //

    if (dataFlag == NodeData::Reaction)
      theDomain->calculateNodalReactions(0);
    else if (dataFlag == NodeData::ReactionInclInertia)
      theDomain->calculateNodalReactions(1);
    if (dataFlag == NodeData::ReactionInclRayleigh)
      theDomain->calculateNodalReactions(2);
 
    for (int i=0; i<numValidNodes; i++) {
      int cnt = i*numDOF;

      if (dataFlag == NodeData::DisplNorm)
	cnt = i;

      Node *theNode = theNodes[i];

      if (dataFlag == NodeData::DisplTrial) {
	const Vector &response = theNode->getTrialDisp();
	for (int j=0; j<numDOF; j++) {
	  int dof = (*theDofs)(j);

	  if (theTimeSeries != 0) {
	    timeSeriesTerm = timeSeriesValues[j];
	  }

	  if (response.Size() > dof) {
	    (*currentData)(cnt) = response(dof) + timeSeriesTerm;
	  }else 
	    (*currentData)(cnt) = 0.0 + timeSeriesTerm;
	  
	  cnt++;
	}

      } else if (dataFlag == NodeData::DisplNorm) {
	const Vector &response = theNode->getTrialDisp();
	double sum = 0.0;
	for (int j=0; j<numDOF; j++) {
	  int dof = (*theDofs)(j);

	  if (theTimeSeries != 0) {
	    timeSeriesTerm = timeSeriesValues[j];
	  }

	  if (response.Size() > dof) {
	    sum += (response(dof) + timeSeriesTerm) * (response(dof) + timeSeriesTerm);    
	  } else 
	    sum += timeSeriesTerm * timeSeriesTerm;    
	}

	(*currentData)(cnt) = sqrt(sum);
	cnt++;

      } else if (dataFlag == NodeData::VelocTrial) {
	const Vector &response = theNode->getTrialVel();
	for (int j=0; j<numDOF; j++) {

	  if (theTimeSeries != 0) {
	    timeSeriesTerm = timeSeriesValues[j];
	  }

	  int dof = (*theDofs)(j);
	  if (response.Size() > dof) {
	    (*currentData)(cnt) = response(dof) + timeSeriesTerm;
	  } else 
	    (*currentData)(cnt) = 0.0 + timeSeriesTerm;
	  
	  cnt++;
	}
      } else if (dataFlag == NodeData::AccelTrial) {
	const Vector &response = theNode->getTrialAccel();
	for (int j=0; j<numDOF; j++) {

	  if (theTimeSeries != 0) {
	    timeSeriesTerm = timeSeriesValues[j];
	  }

	  int dof = (*theDofs)(j);
	  if (response.Size() > dof) {
	    (*currentData)(cnt) = response(dof) + timeSeriesTerm;
	  } else 
	    (*currentData)(cnt) = 0.0 + timeSeriesTerm;
	  
	  cnt++;
	}
      } else if (dataFlag == NodeData::IncrDisp) {
	const Vector &response = theNode->getIncrDisp();
	for (int j=0; j<numDOF; j++) {
	  int dof = (*theDofs)(j);
	  if (response.Size() > dof) {
	    (*currentData)(cnt) = response(dof);
	  } else 
	    (*currentData)(cnt) = 0.0;
	  
	  cnt++;
	}
      } else if (dataFlag == NodeData::IncrDeltaDisp) {
	const Vector &response = theNode->getIncrDeltaDisp();
	for (int j=0; j<numDOF; j++) {
	  int dof = (*theDofs)(j);
	  if (response.Size() > dof) {
	    (*currentData)(cnt) = response(dof);
	  } else 
	    (*currentData)(cnt) = 0.0;
	  
	  cnt++;
	}
      } else if (dataFlag == NodeData::UnbalancedLoad) {
	const Vector &theResponse = theNode->getUnbalancedLoad();
	for (int j=0; j<numDOF; j++) {
	  int dof = (*theDofs)(j);
	  if (theResponse.Size() > dof) {
	    (*currentData)(cnt) = theResponse(dof);
	  } else 
	    (*currentData)(cnt) = 0.0;
	  
	  cnt++;
	}

      } else if (dataFlag == NodeData::UnbalanceInclInertia) {
	const Vector &theResponse = theNode->getUnbalancedLoadIncInertia();
	for (int j=0; j<numDOF; j++) {
	  int dof = (*theDofs)(j);
	  if (theResponse.Size() > dof) {
	    (*currentData)(cnt) = theResponse(dof);
	  } else 
	    (*currentData)(cnt) = 0.0;
	  
	  cnt++;
	}

      } else if (dataFlag == NodeData::Reaction 
              || dataFlag == NodeData::ReactionInclInertia 
              || dataFlag == NodeData::ReactionInclRayleigh) {
	const Vector &theResponse = theNode->getReaction();
	for (int j=0; j<numDOF; j++) {
	  int dof = (*theDofs)(j);
	  if (theResponse.Size() > dof) {
	    (*currentData)(cnt) = theResponse(dof);
	  } else 
	    (*currentData)(cnt) = 0.0;
	  
	  cnt++;
	}

      } else if (dataFlag == NodeData::EigenVector) {
	int mode = dataIndex;
	int column = mode - 1;
	const Matrix &theEigenvectors = theNode->getEigenvectors();
	if (theEigenvectors.noCols() > column) {
	  int noRows = theEigenvectors.noRows();
	  for (int j=0; j<numDOF; j++) {
	    int dof = (*theDofs)(j);
	    if (noRows > dof) {
	      (*currentData)(cnt) = theEigenvectors(dof,column);
	    } else 
	      (*currentData)(cnt) = 0.0;
	    cnt++;		
	  }
	} else {
	  for (int j=0; j<numDOF; j++) {
	    (*currentData)(cnt) = 0.0;
	    cnt++;		
	  }
	}
      }
    }
  }
 
 // add data contribution to runningTotal
  count++;
  int sizeData = currentData->Size();
  for (int i=0; i<sizeData; i++) {
    double value = (*currentData)(i);
    (*runningTotal)(i) += value*value;;
  }
  return 0;
}


int
NodeRecorderRMS::restart(void)
{
  runningTotal->Zero();
  count = 0;
  return 0;
}


int 
NodeRecorderRMS::setDomain(Domain &theDom)
{
  theDomain = &theDom;
  initializationDone = false;  
  return 0;
}


int
NodeRecorderRMS::sendSelf(int commitTag, Channel &theChannel)
{
  addColumnInfo = 1;

  int numDOF = theDofs->Size();

  if (theChannel.isDatastore() == 1) {
    opserr << "NodeRecorderRMS::sendSelf() - does not send data to a datastore\n";
    return -1;
  }

  initializationDone = false;
  static ID idData(7); 
  idData.Zero();

  if (theDofs != 0)
    idData(0) = numDOF;
  if (theNodalTags != 0)
    idData(1) = theNodalTags->Size();
  if (theHandler != 0) {
    idData(2) = theHandler->getClassTag();
  }

  idData(3) = (int)dataFlag;


  idData(5) = this->getTag();

  if (theTimeSeries == 0)
    idData[6] = 0;
  else
    idData[6] = 1;

  if (theChannel.sendID(0, commitTag, idData) < 0) {
    opserr << "NodeRecorderRMS::sendSelf() - failed to send idData\n";
    return -1;
  }

  if (theDofs != 0) 
    if (theChannel.sendID(0, commitTag, *theDofs) < 0) {
      opserr << "NodeRecorderRMS::sendSelf() - failed to send dof id's\n";
      return -1;
    }

  if (theNodalTags != 0)
    if (theChannel.sendID(0, commitTag, *theNodalTags) < 0) {
      opserr << "NodeRecorderRMS::sendSelf() - failed to send nodal tags\n";
      return -1;
    }

  static Vector data(2);
  data(0) = deltaT;
  data(1) = nextTimeStampToRecord;
  data(2) = relDeltaTTol;
  if (theChannel.sendVector(0, commitTag, data) < 0) {
    opserr << "NodeRecorderRMS::sendSelf() - failed to send data\n";
    return -1;
  }

  if (theHandler->sendSelf(commitTag, theChannel) < 0) {
    opserr << "NodeRecorderRMS::sendSelf() - failed to send the DataOutputHandler\n";
    return -1;
  }

  if (theTimeSeries != 0) {
    ID timeSeriesTags(numDOF);
    for (int i=0; i<numDOF; i++) {
      if (theTimeSeries[i] != 0) {
	timeSeriesTags[i] = theTimeSeries[i]->getClassTag();
      } else
	timeSeriesTags[i] = -1;
    }
    if (theChannel.sendID(0, commitTag, timeSeriesTags) < 0) {
      opserr << "NodeRecorderRMS::sendSelf() - failed to send time series tags\n";
      return -1;
    }    
    for (int i=0; i<numDOF; i++) {
      if (theTimeSeries[i] != 0) {	
	if (theTimeSeries[i]->sendSelf(commitTag, theChannel) < 0) {
	  opserr << "NodeRecorderRMS::sendSelf() - time series failed in send\n";
	  return -1;

	}
      }
    }
  }

  return 0;
}

int 
NodeRecorderRMS::recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker)
{
  addColumnInfo = 1;

  if (theChannel.isDatastore() == 1) {
    opserr << "NodeRecorderRMS::sendSelf() - does not send data to a datastore\n";
    return -1;
  }

  static ID idData(7); 
  if (theChannel.recvID(0, commitTag, idData) < 0) {
    opserr << "NodeRecorderRMS::recvSelf() - failed to send idData\n";
    return -1;
  }

  int numDOFs = idData(0);
  int numNodes = idData(1);

  dataFlag = (NodeData)idData(3);

  this->setTag(idData(5));



  //
  // get the DOF ID data
  //

  if (theDofs == 0 || theDofs->Size() != numDOFs) {
    if (theDofs != 0)
      delete theDofs;

    if (numDOFs != 0) {
      theDofs = new ID(numDOFs);
      if (theDofs == 0 || theDofs->Size() != numDOFs) {
	opserr << "NodeRecorderRMS::recvSelf() - out of memory\n";
	return -1;
      }	
    }
  }
  if (theDofs != 0)
    if (theChannel.recvID(0, commitTag, *theDofs) < 0) {
      opserr << "NodeRecorderRMS::recvSelf() - failed to recv dof data\n";
      return -1;
    } 

  //
  // get the NODAL tag data
  //

  if (theNodalTags == 0 || theNodalTags->Size() != numNodes) {
    if (theNodalTags != 0)
      delete theNodalTags;

    if (numNodes != 0) {
      theNodalTags = new ID(numNodes);
      if (theNodalTags == 0 || theNodalTags->Size() != numNodes) {
	opserr << "NodeRecorderRMS::recvSelf() - out of memory\n";
	return -1;
      }	
    }
  }
  if (theNodalTags != 0)
    if (theChannel.recvID(0, commitTag, *theNodalTags) < 0) {
      opserr << "NodeRecorderRMS::recvSelf() - failed to recv dof data\n";
      return -1;
    } 


  static Vector data(2);
  if (theChannel.recvVector(0, commitTag, data) < 0) {
    opserr << "NodeRecorderRMS::sendSelf() - failed to receive data\n";
    return -1;
  }
  deltaT = data(0);
  nextTimeStampToRecord = data(1);
  relDeltaTTol = data(2);
 
  if (theHandler != 0)
    delete theHandler;

  theHandler = theBroker.getPtrNewStream(idData(2));
  if (theHandler == 0) {
    opserr << "NodeRecorderRMS::sendSelf() - failed to get a data output handler\n";
    return -1;
  }

  if (theHandler->recvSelf(commitTag, theChannel, theBroker) < 0) {
    opserr << "NodeRecorderRMS::sendSelf() - failed to send the DataOutputHandler\n";
    return -1;
  }


  if (idData[6] == 1) {
    theTimeSeries = new TimeSeries *[numDOFs];
    ID timeSeriesTags(numDOFs);
    if (theChannel.recvID(0, commitTag, timeSeriesTags) < 0) {
      opserr << "NodeRecorderRMS::recvSelf() - failed to recv time series tags\n";
      return -1;
    }    
    for (int i=0; i<numDOFs; i++) {
      if (timeSeriesTags[i] == -1)
	theTimeSeries[i] = 0;
      else {
	theTimeSeries[i] = theBroker.getNewTimeSeries(timeSeriesTags(i));
	if (theTimeSeries[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	  opserr << "NodeRecorderRMS::recvSelf() - time series failed in recv\n";
	  return -1;
	}
      }
    }
  }

  return 0;
}


int
NodeRecorderRMS::initialize(void)
{
  if (theDofs == 0 || theDomain == 0) {
    opserr << "NodeRecorderRMS::initialize() - either nodes, dofs or domain has not been set\n";
    return -1;
  }

  //
  // create & set nodal array pointer
  //

  if (theNodes != 0) 
    delete [] theNodes;
  
  numValidNodes = 0;

  if (theNodalTags != 0) {
    int numNode = theNodalTags->Size();
    theNodes = new Node *[numNode];
    if (theNodes == 0) {
      opserr << "NodeRecorderRMS::domainChanged - out of memory\n";
      return -1;
    }

    for (int i=0; i<numNode; i++) {
      int nodeTag = (*theNodalTags)(i);
      Node *theNode = theDomain->getNode(nodeTag);
      if (theNode != 0) {
	theNodes[numValidNodes] = theNode;
	numValidNodes++;
      } 
    }
  } else {

    int numNodes = theDomain->getNumNodes();
    if (numNodes != 0) {
      theNodes = new Node *[numNodes];
      
      if (theNodes == 0) {
	opserr << "NodeRecorder::domainChanged - out of memory\n";
	return -1;
      }
      NodeIter &theDomainNodes = theDomain->getNodes();
      Node *theNode;
      numValidNodes = 0;
      while (((theNode = theDomainNodes()) != 0) && (numValidNodes < numNodes)) {
	theNodes[numValidNodes] = theNode;
	numValidNodes++;
      }
    } else
      numValidNodes = 0;
  }

  //
  // need to create the data description, i.e. what each column of data is
  //

  //
  // need to create the data description, i.e. what each column of data is
  //

  char outputData[128];
  char dataType[128];

  if (dataFlag == NodeData::DisplTrial) {
    strcpy(dataType,"D");
  } else if (dataFlag == NodeData::VelocTrial) {
    strcpy(dataType,"V");
  } else if (dataFlag == NodeData::AccelTrial) {
    strcpy(dataType,"A");
  } else if (dataFlag == NodeData::IncrDisp) {
    strcpy(dataType,"dD");
  } else if (dataFlag == NodeData::IncrDeltaDisp) {
    strcpy(dataType,"ddD");
  } else if (dataFlag == NodeData::UnbalancedLoad) {
    strcpy(dataType,"U");
  } else if (dataFlag == NodeData::UnbalanceInclInertia) {
    strcpy(dataType,"U");
  } else if (dataFlag == NodeData::Reaction) {
    strcpy(dataType,"R");
  } else if (dataFlag == NodeData::ReactionInclInertia) {
    strcpy(dataType,"R");
  } else if (dataFlag == NodeData::DisplNorm) {
    strcpy(dataType,"|D|");
  } else if (dataFlag == NodeData::EigenVector) {
    sprintf(dataType,"E%d", dataIndex);
  } else
    strcpy(dataType,"Unknown");

  //
  // resize the output matrix
  //
  int numDOF = theDofs->Size();
  int numValidResponse = numValidNodes*numDOF;

  if (dataFlag == NodeData::DisplNorm)
    numValidResponse = numValidNodes;  

  currentData = new Vector(numValidResponse);
  runningTotal = new Vector(numValidResponse);
  runningTotal->Zero();

  ID dataOrder(numValidResponse);
  ID xmlOrder(numValidNodes);

  if (theNodalTags != nullptr && addColumnInfo == 1) {

    int count = 0;
    int nodeCount = 0;

    int numNode = theNodalTags->Size();
    for (int i=0; i<numNode; i++) {
      int nodeTag = (*theNodalTags)(i);
      Node *theNode = theDomain->getNode(nodeTag);
      if (theNode != 0) {
	xmlOrder(nodeCount++) = i+1;
	for (int j=0; j<numDOF; j++)
	  dataOrder(count++) = i+1;
      }
    }

    theHandler->setOrder(xmlOrder);
  }

  for (int i=0; i<numValidNodes; i++) {
    int nodeTag = theNodes[i]->getTag();

    theHandler->tag("NodeOutput");
    theHandler->attr("nodeTag", nodeTag);

    for (int j=0; j<theDofs->Size(); j++) {
      
      sprintf(outputData, "%s%d", dataType, j+1);
      theHandler->tag("ResponseType",outputData);
    }

    theHandler->endTag();
  }

  if (theNodalTags != 0 && addColumnInfo == 1) {
    theHandler->setOrder(dataOrder);
  }

  initializationDone = true;

  return 0;
}
//added by SAJalali
double NodeRecorderRMS::getRecordedValue(int clmnId, int rowOffset, bool reset)
{
	double res = 0;
	if (!initializationDone)
		return res;
	if (clmnId >= runningTotal->Size() && count == 0)
	  return res;
	if (count != 0)
	  res = sqrt(((*runningTotal)(clmnId)*(*runningTotal)(clmnId))/count);
	if (reset)
	  count = 0;
	return res;
}

int NodeRecorderRMS::flush(void) {
  if (theHandler != 0) {
    return theHandler->flush();
  }
  return 0;
}
