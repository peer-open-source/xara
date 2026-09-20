/* ****************************************************************** **
**    OpenSeess - Open System for Earthquake Engineering Simulation   **
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
                                                                        
// $Revision: 1.19 $
// $Date: 2010-09-16 00:07:11 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/subdomain/ActorSubdomain.cpp,v $
                                                                        
#include <ActorSubdomain.h>
#include <FEM_ObjectBroker.h>
#include <Element.h>
#include <Node.h>
#include <SP_Constraint.h>
#include <MP_Constraint.h>
#include <ElementalLoad.h>
#include <NodalLoad.h>
#include <LoadPattern.h>
#include <Matrix.h>
#include <Vector.h>
#include <DomainDecompositionAnalysis.h>
#include <analysis/criteria/ConvergenceTest.h>

#include <EquiSolnAlgo.h>
#include <IncrementalIntegrator.h>
#include <LinearSOE.h>
#include <EigenSOE.h>
#include <Recorder.h>
#include <Parameter.h>
#include <Message.h>

#include <ArrayOfTaggedObjects.h>
#include <ShadowActorSubdomain.h>

// 2 procedurs defined in SP_Constraint.cpp
int SP_Constraint_GetNextTag();
int SP_Constraint_SetNextTag(int);

ActorSubdomain::ActorSubdomain(Channel &theChannel,
			       FEM_ObjectBroker &theBroker)
:Subdomain(0), Actor(theChannel,theBroker,0),
 msgData(4),lastResponse(0)
{
  // does nothing
}
    
ActorSubdomain::~ActorSubdomain()
{
  // does nothing
}


int
ActorSubdomain::run()
{
    return 0;
}



const Vector &
ActorSubdomain::getLastExternalSysResponse()
{
    int numDOF = this->getNumDOF();
    numDOF = this->getNumDOF();

    if (lastResponse == 0)
		lastResponse = new Vector(numDOF);
	else if (lastResponse->Size() != numDOF) {
		delete lastResponse;
		lastResponse = new Vector(numDOF);
    }
    
    if (mapBuilt == false)
      this->buildMap();

    ID &theMap = *map;
    Vector &localResponse = *lastResponse;
    int numberDOF = this->getNumDOF();
    for (int i=0; i<numberDOF; i++)
      (*mappedVect)(theMap(i)) = localResponse(i);

    return *mappedVect;

}

int
ActorSubdomain::update(void)
{
  int res = this->Domain::update();

  res = this->barrierCheck(res);

  return res;
}

int
ActorSubdomain::updateTimeDt(void)
{
  static Vector data(2);

  this->recvVector(data);

  double newTime = data(0);
  double dT = data(1);
  int res = this->Domain::update(newTime, dT);
  return this->barrierCheck(res);
}

int
ActorSubdomain::barrierCheck(int myResult)
{
  static ID data(1);
  data(0) = myResult;
  this->sendID(data);
  this->recvID(data);

  return data(0);
}







