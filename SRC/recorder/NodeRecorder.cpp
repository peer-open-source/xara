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
// Description: This file contains the class definition for NodeRecorder.
// A NodeRecorder is used to record the specified dof responses 
// at a collection of nodes over an analysis. (between commitTag of 0 and
// last commitTag).
//
// Written: fmk 
//
#include <NodeRecorder.h>
#include <Domain.h>
#include <Parameter.h>
#include <Node.h>
#include <NodeIter.h>
#include <Vector.h>
#include <ID.h>
#include <Matrix.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Pressure_Constraint.h>
#include <MeshRegion.h>
#include <TimeSeries.h>

#include <StandardStream.h>
#include <DataFileStream.h>
#include <DataFileStreamAdd.h>
#include <XmlFileStream.h>
#include <BinaryFileStream.h>
#include <TCP_Stream.h>

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

NodeRecorder::NodeRecorder()
: Recorder(RECORDER_TAGS_NodeRecorder),
  theDofs(0), theNodalTags(0), theNodes(0), response(0),
  theDomain(0), theOutputHandler(0),
  echoTimeFlag(true), dataFlag(NodeData::DisplTrial), 
  deltaT(0), nextTimeStampToRecord(0.0), 
  gradIndex(-1),
  initializationDone(false), numValidNodes(0), 
  addColumnInfo(0), theTimeSeries(nullptr), timeSeriesValues(0)
{

}


NodeRecorder::NodeRecorder(const ID &dofs, 
                           const ID *nodes, 
                           int pgradIndex,
                           NodeData _dataFlag,
                           int _dataIndex,
                           Domain &theDom,
                           OPS_Stream &theOutput,
                           double dT,
                           double _relDeltaTTol,
                           bool timeFlag,
                           TimeSeries **theSeries)
:Recorder(RECORDER_TAGS_NodeRecorder),
 theDofs(0), theNodalTags(0), theNodes(0), response(0),
 theDomain(&theDom), theOutputHandler(&theOutput),
 echoTimeFlag(timeFlag), dataFlag(_dataFlag), dataIndex(_dataIndex),
 deltaT(dT), relDeltaTTol(_relDeltaTTol), nextTimeStampToRecord(0.0), 
 gradIndex(pgradIndex), 
 initializationDone(false), numValidNodes(0), addColumnInfo(0), 
 theTimeSeries(theSeries), timeSeriesValues(0)
{

  //
  // store copy of dof's to be recorder, verifying dof are valid, i.e. >= 0
  //

  int numDOF = dofs.Size();

  if (numDOF != 0) {

    theDofs = new ID(numDOF);

    int count = 0;
    int i;
    for (i=0; i<numDOF; i++) {
      int dof = dofs(i);
      if (dof >= 0) {
	(*theDofs)[count] = dof;
	count++;
      } else {
	opserr << "NodeRecorder::NodeRecorder - invalid dof  " << dof;
	opserr << " will be ignored\n";
      }
    }
  }

  //
  // create memory to hold nodal ID's
  //
  if (nodes != nullptr) {
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

  //
  // set the data flag used as a switch to get the response in a record
  //
  if (dataFlag == NodeData::Reaction
   || dataFlag == NodeData::ReactionInclInertia
   || dataFlag == NodeData::ReactionInclRayleigh) {
    if (timeFlag == true)
      theOutputHandler->setAddCommon(2);
    else
      theOutputHandler->setAddCommon(1);
  }
}


NodeRecorder::~NodeRecorder()
{
  if (theOutputHandler != nullptr) {
    theOutputHandler->endTag(); // Data
    delete theOutputHandler;
  }

  if (theDofs != nullptr) {
    int numDOF = theDofs->Size();
    delete theDofs;

    if (theTimeSeries != nullptr) {
      for (int i=0; i<numDOF; i++) {
        if (theTimeSeries[i] != nullptr)
          delete theTimeSeries[i];
      }
      delete [] theTimeSeries;
    }
  } else {
    assert(theTimeSeries == nullptr);
  }
  
  if (timeSeriesValues != nullptr) 
    delete [] timeSeriesValues;

  if (theNodalTags != nullptr)
    delete theNodalTags;

  if (theNodes != nullptr)
    delete [] theNodes;
}


int
NodeRecorder::record(int commitTag, double timeStamp)
{
  if (theDomain == nullptr || theDofs == nullptr)
    return 0;

  if (theOutputHandler == nullptr) {
    opserr << "NodeRecorder::record() - no DataOutputHandler has been set\n";
    return -1;
  }

  if (initializationDone == false) {
    if (this->initialize() != 0) {
      opserr << "NodeRecorder::record() - failed in initialize()\n";
      return -1;
    }
  }

  int numDOF = theDofs->Size();

  // where relDeltaTTol is the maximum reliable ratio between analysis time step and deltaT
  // to provide adequate tolerance for floating point precision (default=1.0e-5)
    if (deltaT == 0.0 || timeStamp - nextTimeStampToRecord >= -deltaT * relDeltaTTol) {

    if (deltaT != 0.0)
      nextTimeStampToRecord = timeStamp + deltaT;

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

    //
    // add time information if requested
    //
    int timeOffset = 0;
    if (echoTimeFlag == true) {
      timeOffset = 1;
      response(0) = timeStamp;
    }

    double timeSeriesTerm = 0.0;

    if (theTimeSeries != nullptr) {
      for (int i=0; i<numDOF; i++) { 
	if (theTimeSeries[i] != nullptr) 
	  timeSeriesValues[i] = theTimeSeries[i]->getFactor(timeStamp);
      }
    }

    //
    // now we go get the responses from the nodes & place them in disp vector
    //
    if (dataFlag != NodeData::Empty) { // != 10

      for (int i=0; i<numValidNodes; i++) {

	int cnt = i*numDOF + timeOffset; 
	if (dataFlag == NodeData::DisplNorm 
         || dataFlag == NodeData::Pressure)
	  cnt = i + timeOffset;

	Node *theNode = theNodes[i];
	if (dataFlag == NodeData::DisplTrial) { // == 0

	  // AddingSensitivity:BEGIN ///////////////////////////////////
	  if (gradIndex < 0) {
	    const Vector &theResponse = theNode->getTrialDisp();
	    for (int j=0; j<numDOF; j++) {

	      if (theTimeSeries != nullptr)
		timeSeriesTerm = timeSeriesValues[j];

	      int dof = (*theDofs)(j);
	      if (theResponse.Size() > dof)
		response(cnt) = theResponse(dof)  + timeSeriesTerm;
	      else
		response(cnt) = 0.0 + timeSeriesTerm;

	      cnt++;
	    }
	  }
	  else {
	    for (int j=0; j<numDOF; j++) {
	      int dof = (*theDofs)(j);
	      response(cnt) = theNode->getDispSensitivity(dof+1, gradIndex);
	      cnt++;
	    }
	  } 
	  // AddingSensitivity:END /////////////////////////////////////
	} else if (dataFlag == NodeData::VelocTrial) {

	  const Vector &theResponse = theNode->getTrialVel();
	  for (int j=0; j<numDOF; j++) {

	    if (theTimeSeries != nullptr) {
	      timeSeriesTerm = timeSeriesValues[j];
	    }

	    int dof = (*theDofs)(j);
	    if (theResponse.Size() > dof) {
	      response(cnt) = theResponse(dof) + timeSeriesTerm;
	    } else
	      response(cnt) = 0.0 + timeSeriesTerm;

	    cnt++;
	  }


	} else if (dataFlag == NodeData::DisplNorm) {
	  const Vector &theResponse = theNode->getTrialDisp();
	  double sum = 0.0;

	  for (int j=0; j<numDOF; j++) {

	    if (theTimeSeries != nullptr)
	      timeSeriesTerm = timeSeriesValues[j];

	    int dof = (*theDofs)(j);
	    if (theResponse.Size() > dof) {
	      sum += (theResponse(dof) + timeSeriesTerm) * (theResponse(dof) + timeSeriesTerm);
	    } else
	      sum += timeSeriesTerm * timeSeriesTerm;
	  }

	  response(cnt) = sqrt(sum);
	  cnt++;

	} else if (dataFlag == NodeData::Pressure) {

	    if (theTimeSeries != nullptr)
              timeSeriesTerm = timeSeriesValues[0];
	    
	    // get node pressure
	    double pressure = timeSeriesTerm;
	    Pressure_Constraint* pc = theDomain->getPressure_Constraint(theNode->getTag());
	    if (pc != nullptr) {
	      pressure += pc->getPressure();
	    }

	    response(cnt) = pressure;
	    cnt++;

	} else if (dataFlag == NodeData::AccelTrial) {

	  const Vector &theResponse = theNode->getTrialAccel();
	  for (int j=0; j<numDOF; j++) {

	    if (theTimeSeries != nullptr)
	      timeSeriesTerm = timeSeriesValues[j];

	    int dof = (*theDofs)(j);
	    if (theResponse.Size() > dof) {
	      response(cnt) = theResponse(dof) + timeSeriesTerm;
	    } else
	      response(cnt) = 0.0 + timeSeriesTerm;

	    cnt++;
	  }

	} else if (dataFlag == NodeData::IncrDisp) {
	  const Vector &theResponse = theNode->getIncrDisp();
	  for (int j=0; j<numDOF; j++) {

	    if (theTimeSeries != nullptr)
	      timeSeriesTerm = timeSeriesValues[j];

	    int dof = (*theDofs)(j);
	    if (theResponse.Size() > dof) {
	      response(cnt) = theResponse(dof);
	    } else
	      response(cnt) = 0.0;

	    cnt++;
	  }
	} else if (dataFlag == NodeData::IncrDeltaDisp) {
	  const Vector &theResponse = theNode->getIncrDeltaDisp();
	  for (int j=0; j<numDOF; j++) {
	    int dof = (*theDofs)(j);
	    if (theResponse.Size() > dof) {
	      response(cnt) = theResponse(dof);
	    } else
	      response(cnt) = 0.0;

	    cnt++;
	  }
	} else if (dataFlag == NodeData::UnbalancedLoad) {
	  const Vector &theResponse = theNode->getUnbalancedLoad();
	  for (int j=0; j<numDOF; j++) {
	    int dof = (*theDofs)(j);
	    if (theResponse.Size() > dof) {
	      response(cnt) = theResponse(dof);
	    } else
	      response(cnt) = 0.0;

	    cnt++;
	  }
	  
	} else if (dataFlag == NodeData::UnbalanceInclInertia) {
	  const Vector &theResponse = theNode->getUnbalancedLoadIncInertia();
	  for (int j=0; j<numDOF; j++) {
	    int dof = (*theDofs)(j);
	    if (theResponse.Size() > dof) {
	      response(cnt) = theResponse(dof);
	    } else
	      response(cnt) = 0.0;

	    cnt++;
	  }

	} else if (dataFlag == NodeData::RayleighForces) {
	  const Vector *theResponse = theNode->getResponse(NodeData::RayleighForces);
	  for (int j=0; j<numDOF; j++) {
	    int dof = (*theDofs)(j);
	    if (theResponse->Size() > dof) {
	      response(cnt) = (*theResponse)(dof);
	    } else
	      response(cnt) = 0.0;
	    cnt++;
	  }

	} else if (dataFlag == NodeData::Reaction 
                || dataFlag == NodeData::ReactionInclInertia 
                || dataFlag == NodeData::ReactionInclRayleigh) {
	  const Vector &theResponse = theNode->getReaction();
	  for (int j=0; j<numDOF; j++) {
	    int dof = (*theDofs)(j);
	    if (theResponse.Size() > dof) {
	      response(cnt) = theResponse(dof);
	    } else
	      response(cnt) = 0.0;
	    cnt++;
	  }
	  
	} else if (dataFlag == NodeData::EigenVector){// (10 <= dataFlag  && dataFlag < 1000) {
	  int mode = dataIndex; // dataFlag - 10;
	  int column = mode - 1;

	  const Matrix &theEigenvectors = theNode->getEigenvectors();
	  if (theEigenvectors.noCols() > column) {
	    int noRows = theEigenvectors.noRows();
	    for (int j=0; j<numDOF; j++) {
	      int dof = (*theDofs)(j);
	      if (noRows > dof) {
		response(cnt) = theEigenvectors(dof,column);
	      } else
		response(cnt) = 0.0;
	      cnt++;
	    }
	  }
	  
	} else if (dataFlag == NodeData::DisplSensitivity) { // (dataFlag  >= 1000 && dataFlag < 2000) {
	  int grad = dataIndex; // dataFlag - 1000;
	  
	  for (int j=0; j<numDOF; j++) {
	    int dof = (*theDofs)(j);
	    dof += 1; // Terje uses 1 through DOF for the dof indexing
	    response(cnt) = theNode->getDispSensitivity(dof, grad);
	    cnt++;
	  }
	  
	} else if (dataFlag == NodeData::VelocSensitivity) { // (dataFlag  >= 2000 && dataFlag < 3000) {
	  int grad = dataIndex; // dataFlag - 2000;
	  
	  for (int j=0; j<numDOF; j++) {
	    int dof = (*theDofs)(j);
	    dof += 1; // Terje uses 1 through DOF for the dof indexing
	    response(cnt) = theNode->getVelSensitivity(dof, grad);
	    cnt++;
	  }

	} else if (dataFlag == NodeData::AccelSensitivity) { // (dataFlag  >= 3000) {
	  int grad = dataIndex; // dataFlag - 3000;
	  
	  for (int j=0; j<numDOF; j++) {
	    int dof = (*theDofs)(j);
	    dof += 1; // Terje uses 1 through DOF for the dof indexing
	    response(cnt) = theNode->getAccSensitivity(dof, grad);
	    cnt++;
	  }
	}

	else {
	  // unknown response
	  for (int j=0; j<numDOF; j++) {
	    response(cnt) = 0.0;
	  }
	}
      }

      // insert the data into the database
      theOutputHandler->write(response);

    } else { // output all eigenvalues

      Node *theNode = theNodes[0];
      const Matrix &theEigenvectors = theNode->getEigenvectors();
      int numValidModes = theEigenvectors.noCols();

      for (int mode=0; mode<numValidModes; mode++) {

	for (int i=0; i<numValidNodes; i++) {
	  int cnt = i*numDOF + timeOffset;
	  theNode = theNodes[i];
	  int column = mode;

	  const Matrix &theEigenvectors = theNode->getEigenvectors();
	  if (theEigenvectors.noCols() > column) {
	    int noRows = theEigenvectors.noRows();
	    for (int j=0; j<numDOF; j++) {
	      int dof = (*theDofs)(j);
	      if (noRows > dof) {
		response(cnt) = theEigenvectors(dof,column);
	      } else
		response(cnt) = 0.0;
	      cnt++;
	    }
	  }
	}
	theOutputHandler->write(response);
      }
    }
  }

  return 0;
}


int
NodeRecorder::setDomain(Domain &theDom)
{
  theDomain = &theDom;
  return 0;
}


int
NodeRecorder::sendSelf(int commitTag, Channel &theChannel)
{
  addColumnInfo = 1;

  if (theChannel.isDatastore() == 1) {
    opserr << "NodeRecorder::sendSelf() - does not send data to a datastore\n";
    return -1;
  }

  initializationDone = false;

  int numDOF = theDofs->Size();

  static ID idData(8);
  idData.Zero();
  if (theDofs != 0)
    idData(0) = numDOF;
  if (theNodalTags != 0)
    idData(1) = theNodalTags->Size();
  if (theOutputHandler != 0) {
    idData(2) = theOutputHandler->getClassTag();
  }

  if (echoTimeFlag == true)
    idData(3) = 1;
  else
    idData(3) = 0;

  idData(4) = (int)dataFlag;
  idData(5) = gradIndex;

  idData(6) = this->getTag();
  if (theTimeSeries == nullptr)
    idData[7] = 0;
  else
    idData[7] = 1;


  if (theChannel.sendID(0, commitTag, idData) < 0) {
    opserr << "NodeRecorder::sendSelf() - failed to send idData\n";
    return -1;
  }

  if (theDofs != 0)
    if (theChannel.sendID(0, commitTag, *theDofs) < 0) {
      opserr << "NodeRecorder::sendSelf() - failed to send dof id's\n";
      return -1;
    }

  if (theNodalTags != 0)
    if (theChannel.sendID(0, commitTag, *theNodalTags) < 0) {
      opserr << "NodeRecorder::sendSelf() - failed to send nodal tags\n";
      return -1;
    }

  static Vector data(3);
  data(0) = deltaT;
  data(1) = nextTimeStampToRecord;
  data(2) = relDeltaTTol;
  if (theChannel.sendVector(0, commitTag, data) < 0) {
    opserr << "NodeRecorder::sendSelf() - failed to send data\n";
    return -1;
  }

  if (theOutputHandler->sendSelf(commitTag, theChannel) < 0) {
    opserr << "NodeRecorder::sendSelf() - failed to send the DataOutputHandler\n";
    return -1;
  }


  if (theTimeSeries != nullptr) {
    ID timeSeriesTags(numDOF);
    for (int i=0; i<numDOF; i++) {
      if (theTimeSeries[i] != 0) {
	timeSeriesTags[i] = theTimeSeries[i]->getClassTag();
      } else
	timeSeriesTags[i] = -1;
    }
    if (theChannel.sendID(0, commitTag, timeSeriesTags) < 0) {
      opserr << "EnvelopeNodeRecorder::sendSelf() - failed to send time series tags\n";
      return -1;
    }
    for (int i=0; i<numDOF; i++) {
      if (theTimeSeries[i] != 0) {
	if (theTimeSeries[i]->sendSelf(commitTag, theChannel) < 0) {
	  opserr << "EnvelopeNodeRecorder::sendSelf() - time series failed in send\n";
	  return -1;

	}
      }
    }
  }


  return 0;
}


int
NodeRecorder::recvSelf(int commitTag, Channel &theChannel,
		       FEM_ObjectBroker &theBroker)
{
  addColumnInfo = 1;

  if (theChannel.isDatastore() == 1) {
    opserr << "NodeRecorder::sendSelf() - does not send data to a datastore\n";
    return -1;
  }

  static ID idData(8);
  if (theChannel.recvID(0, commitTag, idData) < 0) {
    opserr << "NodeRecorder::recvSelf() - failed to send idData\n";
    return -1;
  }


  int numDOFs = idData(0);
  int numNodes = idData(1);

  this->setTag(idData(6));

  if (idData(3) == 1)
    echoTimeFlag = true;
  else
    echoTimeFlag = false;

  dataFlag = (NodeData)idData(4);
  gradIndex = idData(5);

  //
  // get the DOF ID data
  //

  if (theDofs == 0 || theDofs->Size() != numDOFs) {
    if (theDofs != 0)
      delete theDofs;

    if (numDOFs != 0) {
      theDofs = new ID(numDOFs);
      if (theDofs == 0 || theDofs->Size() != numDOFs) {
	opserr << "NodeRecorder::recvSelf() - out of memory\n";
	return -1;
      }
    }
  }
  if (theDofs != 0)
    if (theChannel.recvID(0, commitTag, *theDofs) < 0) {
      opserr << "NodeRecorder::recvSelf() - failed to recv dof data\n";
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
	opserr << "NodeRecorder::recvSelf() - out of memory\n";
	return -1;
      }
    }
  }
  if (theNodalTags != 0)
    if (theChannel.recvID(0, commitTag, *theNodalTags) < 0) {
      opserr << "NodeRecorder::recvSelf() - failed to recv dof data\n";
      return -1;
    }


  static Vector data(3);
  if (theChannel.recvVector(0, commitTag, data) < 0) {
    opserr << "NodeRecorder::sendSelf() - failed to receive data\n";
    return -1;
  }
  deltaT = data(0);
  nextTimeStampToRecord = data(1);
  relDeltaTTol = data(2);


  if (theOutputHandler != 0)
    delete theOutputHandler;

  theOutputHandler = theBroker.getPtrNewStream(idData(2));
  if (theOutputHandler == 0) {
    opserr << "NodeRecorder::sendSelf() - failed to get a data output handler\n";
    return -1;
  }

  if (theOutputHandler->recvSelf(commitTag, theChannel, theBroker) < 0) {
    opserr << "NodeRecorder::sendSelf() - failed to send the DataOutputHandler\n";
    return -1;
  }

  if (idData[7] == 1) {

    timeSeriesValues = new double [numDOFs];
    for (int i=0; i<numDOFs; i++)
      timeSeriesValues[i] = 0.0;

    theTimeSeries = new TimeSeries *[numDOFs];
    ID timeSeriesTags(numDOFs);
    if (theChannel.recvID(0, commitTag, timeSeriesTags) < 0) {
      opserr << "EnvelopeNodeRecorder::recvSelf() - failed to recv time series tags\n";
      return -1;
    }    
    for (int i=0; i<numDOFs; i++) {
      if (timeSeriesTags[i] == -1)
	theTimeSeries[i] = 0;
      else {
	theTimeSeries[i] = theBroker.getNewTimeSeries(timeSeriesTags(i));
	if (theTimeSeries[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	  opserr << "EnvelopeNodeRecorder::recvSelf() - time series failed in recv\n";
	  return -1;
	} 
	  }
	}
  }
  
  return 0;
}


int
NodeRecorder::domainChanged(void)
{
  return 0;
}


int
NodeRecorder::initialize(void)
{
  if (theDofs == 0 || theDomain == 0) {
    opserr << "NodeRecorder::initialize() - either nodes, dofs or domain has not been set\n";
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
      opserr << "NodeRecorder::domainChanged - out of memory\n";
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
  }

  //
  // resize the response vector
  //
  int timeOffset = 0;
  if (echoTimeFlag == true)
    timeOffset = 1;

  int numValidResponse = numValidNodes*theDofs->Size() + timeOffset;
  if (dataFlag == NodeData::DisplNorm || dataFlag == NodeData::Pressure) {
      numValidResponse = numValidNodes + timeOffset;
  }

  response.resize(numValidResponse);
  response.Zero();

  ID orderResponse(numValidResponse);

  //
  // need to create the data description, i.e. what each column of data is
  //
  
  char outputData[256];
  char dataType[256];

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

  int numDOF = theDofs->Size();
  // 
  // write out info to handler if parallel execution
  //
  ID xmlOrder(numValidNodes);

  if (echoTimeFlag == true)  
    xmlOrder.resize(numValidNodes+1);

  if (theNodalTags != 0 && addColumnInfo == 1) {

    int numNode = theNodalTags->Size();
    int count = 0;
    int nodeCount = 0;

    if (echoTimeFlag == true)  {
      orderResponse(count++) = 0;
      xmlOrder(nodeCount++) = 0;
    }
    
    for (int i=0; i<numNode; i++) {
      int nodeTag = (*theNodalTags)(i);
      Node *theNode = theDomain->getNode(nodeTag);
      if (theNode != 0) {
        xmlOrder(nodeCount++) = i+1;
        for (int j=0; j<numDOF; j++)
          orderResponse(count++) = i+1;
      }
    }

    theOutputHandler->setOrder(xmlOrder);
  }

  char nodeCrdData[20];
  sprintf(nodeCrdData,"coord");

  if (echoTimeFlag == true) {
    if (theNodalTags != 0 && addColumnInfo == 1) {
      theOutputHandler->tag("TimeOutput");
      theOutputHandler->tag("ResponseType", "time");
      theOutputHandler->endTag();
    }
  }

  for (int i=0; i<numValidNodes; i++) {
    int nodeTag = theNodes[i]->getTag();
    const Vector &nodeCrd = theNodes[i]->getCrds();
    int numCoord = nodeCrd.Size();


    theOutputHandler->tag("NodeOutput");
    theOutputHandler->attr("nodeTag", nodeTag);

    for (int j=0; j<3; j++) {
      sprintf(nodeCrdData,"coord%d",j+1);
      if (j < numCoord)
        theOutputHandler->attr(nodeCrdData, nodeCrd(j));      
      else
        theOutputHandler->attr(nodeCrdData, 0.0);      
    }

    for (int k=0; k<theDofs->Size(); k++) {
      sprintf(outputData, "%s%d", dataType, k+1);
      theOutputHandler->tag("ResponseType",outputData);
    }

    theOutputHandler->endTag();
  }

  if (theNodalTags != 0 && addColumnInfo == 1) {
    theOutputHandler->setOrder(orderResponse);
  }

  theOutputHandler->tag("Data");
  initializationDone = true;

  return 0;
}

//added by SAJalali
double NodeRecorder::getRecordedValue(int clmnId, int rowOffset, bool reset)
{
  double res = 0;
  if (!initializationDone)
    return res;

  if (clmnId >= response.Size())
    return res;

  res = response(clmnId);
  return res;
}

int NodeRecorder::flush(void) {
  if (theOutputHandler != 0) {
    return theOutputHandler->flush();
  }
  return 0;
}
