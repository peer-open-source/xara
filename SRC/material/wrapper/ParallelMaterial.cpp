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
                                                                        
// $Revision: 1.12 $
// $Date: 2007-02-02 01:19:30 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ParallelMaterial.cpp,v $
                                                                        
                                                                        
// File: ~/material/ParallelModel.C
//
// Written: fmk 
// Created: 07/98
// Revision: A
//
// Description: This file contains the class definition for 
// ParallelModel. ParallelModel is an aggregation
// of UniaxialMaterial objects all considered acting in parallel.
//
// What: "@(#) ParallelModel.C, revA"

#include <ParallelMaterial.h>
#include <ID.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <stdlib.h>
#include <string.h>
#include <MaterialResponse.h>


ParallelMaterial::ParallelMaterial(
				 int tag, 
				 int num, 
				 UniaxialMaterial **theMaterialModels,
				 Vector *factors)
:UniaxialMaterial(tag, MAT_TAG_ParallelMaterial),
 trialStrain(0.0), trialStrainRate(0.0), numMaterials(num),
 theModels(0), theFactors(0)
{
    // create an array (theModels) to store copies of the MaterialModels
    theModels = new UniaxialMaterial *[num];

    // into the newly created array store a pointer to a copy
    // of the UniaxialMaterial stored in theMaterialModels
    for (int i=0; i<num; i++) {
      theModels[i] = theMaterialModels[i]->getCopy();
    }

    // copy the factors
    if (factors != 0) {
      theFactors = new Vector(*factors);
    }
}



// this constructor is used for a ParallelMaterailModel object that
// needs to be constructed in a remote actor process. recvSelf() needs
// to be called on this object
ParallelMaterial::ParallelMaterial()
:UniaxialMaterial(0,MAT_TAG_ParallelMaterial),
 trialStrain(0.0), trialStrainRate(0.0), numMaterials(0),
 theModels(0), theFactors(0)
{

}


ParallelMaterial::~ParallelMaterial()
{
  // invoke the destructor on each MaterialModel object
  for (int i=0; i<numMaterials; i++)
    delete theModels[i];

  // now we can delete the array
  if (theModels != 0) // just in case blank constructor called and no recvSelf()
    delete [] theModels;

  // clean up the factors
  if (theFactors != 0)
      delete theFactors;
}



int 
ParallelMaterial::setTrialStrain(double strain, double strainRate)
{
  // set the trialStrain and the trialStrain in each of the
  // local MaterialModel objects 
  trialStrain = strain;
  trialStrainRate = strainRate;

  for (int i=0; i<numMaterials; i++)
    theModels[i]->setTrialStrain(strain, strainRate);

  return 0;
}


double 
ParallelMaterial::getStrain()
{
  return trialStrain;
}

double 
ParallelMaterial::getStrainRate()
{
    return trialStrainRate;
}

double 
ParallelMaterial::getStress()
{
  // get the stress = sum of stress in all local MaterialModel objects
  double stress = 0.0;
  if (theFactors == 0) {
      for (int i=0; i<numMaterials; i++)
          stress += theModels[i]->getStress();
  } else {
      for (int i=0; i<numMaterials; i++)
          stress += (*theFactors)(i) * theModels[i]->getStress();
  }

  return stress;
}



double 
ParallelMaterial::getTangent()
{
  // get the tangent = sum of tangents in all local MaterialModel objects    
  double E = 0.0;
  if (theFactors == 0) {
      for (int i=0; i<numMaterials; i++)
          E += theModels[i]->getTangent();
  } else {
      for (int i=0; i<numMaterials; i++)
          E += (*theFactors)(i) * theModels[i]->getTangent();
  }

  return E;
}

double 
ParallelMaterial::getInitialTangent()
{
  // get the tangent = sum of tangents in all local MaterialModel objects    
  double E = 0.0;
  if (theFactors == 0) {
      for (int i=0; i<numMaterials; i++)
          E += theModels[i]->getInitialTangent();
  } else {
      for (int i=0; i<numMaterials; i++)
          E += (*theFactors)(i) * theModels[i]->getInitialTangent();
  }

  return E;
}

double 
ParallelMaterial::getDampTangent()
{
  // get the damp tangent = sum of damp tangents in all local MaterialModel objects    
  double eta = 0.0;
  if (theFactors == 0) {
      for (int i=0; i<numMaterials; i++)
          eta += theModels[i]->getDampTangent();
  } else {
      for (int i=0; i<numMaterials; i++)
          eta += (*theFactors)(i) * theModels[i]->getDampTangent();
  }

  return eta;
}

int 
ParallelMaterial::commitState()
{
  // invoke commitState() on each of local MaterialModel objects
  for (int i=0; i<numMaterials; i++)
    if (theModels[i]->commitState() != 0) {
      opserr << "WARNING ParallelMaterial::commitState() ";
      opserr << "MaterialModel failed to commitState():" ;
      theModels[i]->Print(opserr);
    }
  
  return 0;    
}

int 
ParallelMaterial::revertToLastCommit(void)
{
  // invoke commitState() on each of local MaterialModel objects
  for (int i=0; i<numMaterials; i++)
    if (theModels[i]->revertToLastCommit() != 0) {
      opserr << "WARNING ParallelMaterial::revertToLastCommit() ";
      opserr << "MaterialModel failed to revertToLastCommit():" ;
      theModels[i]->Print(opserr);
    }
  
    return 0;    
}


int 
ParallelMaterial::revertToStart(void)
{
    trialStrain = 0.0;
    trialStrainRate = 0.0;

    // invoke commitState() on each of local MaterialModel objects
    for (int i=0; i<numMaterials; i++)
	if (theModels[i]->revertToStart() != 0) {
	    opserr << "WARNING ParallelMaterial::revertToStart() ";
	    opserr << "MaterialModel failed to revertToStart():" ;
	    theModels[i]->Print(opserr);
	}
    
    return 0;    
}



UniaxialMaterial *
ParallelMaterial::getCopy()
{
  ParallelMaterial *theCopy = new 
    ParallelMaterial(this->getTag(), numMaterials, theModels, theFactors);
  
  theCopy->trialStrain = trialStrain;
  theCopy->trialStrainRate = trialStrainRate;
  
  return theCopy;
}


int 
ParallelMaterial::sendSelf(int cTag, Channel &theChannel)
{
    int res = 0;

    static ID data(3);

    // send ID of size 3 so no possible conflict with classTags ID
    int dbTag = this->getDbTag();
    data(0) = this->getTag();
    data(1) = numMaterials;
    data(2) = 0;
    if (theFactors != 0)
        data(2) = 1;

    res = theChannel.sendID(dbTag, cTag, data);
    if (res < 0) {
      opserr << "ParallelMaterial::sendSelf() - failed to send data\n";
      return res;
    }

    // send the factors if not null
    if (theFactors != 0) {
        res = theChannel.sendVector(dbTag, cTag, *theFactors);
        if (res < 0) {
            opserr << "ParallelMaterial::sendSelf() - failed to send factors\n";
            return res;
        }
    }

    // now create an ID containing the class tags and dbTags of all
    // the MaterialModel objects in this ParallelMaterial
    // then send each of the MaterialModel objects
    ID classTags(numMaterials*2);
    for (int i=0; i<numMaterials; i++) {
	classTags(i) = theModels[i]->getClassTag();
	int matDbTag = theModels[i]->getDbTag();
	if (matDbTag == 0) {
	  matDbTag  = theChannel.getDbTag();
	  if (matDbTag != 0)
	    theModels[i]->setDbTag(matDbTag);
	}
	classTags(i+numMaterials) = matDbTag;
    }

    res = theChannel.sendID(dbTag, cTag, classTags);
    if (res < 0) {
      opserr << "ParallelMaterial::sendSelf() - failed to send classTags\n";
      return res;
    }

    for (int j=0; j<numMaterials; j++)
	theModels[j]->sendSelf(cTag, theChannel);
    
    return 0;
}

int 
ParallelMaterial::recvSelf(int cTag, Channel &theChannel, 
				FEM_ObjectBroker &theBroker)
{
    int res = 0;
    static ID data(3);
    int dbTag = this->getDbTag();

    res = theChannel.recvID(dbTag, cTag, data);
    if (res < 0) {
      opserr << "ParallelMaterial::recvSelf() - failed to receive data\n";
      return res;
    }

    this->setTag(int(data(0)));
    int numMaterialsSent = int(data(1));
    if (numMaterials != numMaterialsSent) { 
      numMaterials = numMaterialsSent;
      if (theModels != 0) {
	for (int i=0; i<numMaterials; i++)
	  delete theModels[i];

	delete [] theModels;
      }

      theModels = new UniaxialMaterial *[numMaterials];      
      if (theModels == 0) {
	opserr << "FATAL ParallelMaterial::recvSelf() - ran out of memory";
	opserr << " for array of size: " << numMaterials << "\n";
	return -2;
      }
      for (int i=0; i<numMaterials; i++)
	theModels[i] = 0;
    }

    if (data(2) == 1) {
        theFactors = new Vector(numMaterials);
        res = theChannel.recvVector(dbTag, cTag, *theFactors);
        if (res < 0) {
            opserr << "ParallelMaterial::recvSelf() - failed to receive factors\n";
        return res;
        }
    }

    // create and receive an ID for the classTags and dbTags of the local 
    // MaterialModel objects
    ID classTags(numMaterials*2);
    res = theChannel.recvID(dbTag, cTag, classTags);
    if (res < 0) {
      opserr << "ParallelMaterial::recvSelf() - failed to receive classTags\n";
      return res;
    }

    // now for each of the MaterialModel objects, create a new object
    // and invoke recvSelf() on it
    for (int i=0; i<numMaterials; i++) {
      int matClassTag = classTags(i);
      if (theModels[i] == nullptr || theModels[i]->getClassTag() != matClassTag) {
        if (theModels[i] != nullptr)
          delete theModels[i];

        UniaxialMaterial *theMaterialModel = 
            theBroker.getNewUniaxialMaterial(matClassTag);

        if (theMaterialModel != 0) {
            theModels[i] = theMaterialModel;
            theMaterialModel->setDbTag(classTags(i+numMaterials));
        }
        else {
            opserr << "FATAL ParallelMaterial::recvSelf() ";
            opserr << " could not get a UniaxialMaterial \n";
            exit(-1);
        }
      }
      theModels[i]->recvSelf(cTag, theChannel, theBroker);
    }
    return 0;
}

void 
ParallelMaterial::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
        s << "ParallelMaterial tag: " << this->getTag() << endln;
        for (int i = 0; i < numMaterials; i++) {
            s << " ";
            theModels[i]->Print(s, flag);
        }
        if (theFactors != 0)
            opserr << " Factors: " << *theFactors;
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << OPS_PRINT_JSON_MATE_INDENT << "{";
        s << "\"name\": \"" << this->getTag() << "\", ";
        s << "\"type\": \"ParallelMaterial\", ";
        s << "\"materials\": [";
        for (int i = 0; i < numMaterials-1; i++)
            s << "\"" << theModels[i]->getTag() << "\", ";
        s << "\"" << theModels[numMaterials-1]->getTag() << "\"]}";
        if (theFactors != 0) {
            s << "\"factors\": [";
            for (int i = 0; i < numMaterials-1; i++)
                s << (*theFactors)(i) << ", ";
            s << (*theFactors)(numMaterials-1) << "]}";
        }
        return;
    }
}

Response*
ParallelMaterial::setResponse(const char **argv, int argc, OPS_Stream &theOutput)
{
  Response *theResponse = 0;

  theOutput.tag("UniaxialMaterialOutput");
  theOutput.attr("matType", this->getClassType());
  theOutput.attr("matTag", this->getTag());

  // stress
  if (strcmp(argv[0],"stress") == 0) {
    theOutput.tag("ResponseType", "sigma11");
    theResponse =  new MaterialResponse(this, 1, this->getStress());
  }  
  // tangent
  else if (strcmp(argv[0],"tangent") == 0) {
    theOutput.tag("ResponseType", "C11");
    theResponse =  new MaterialResponse(this, 2, this->getTangent());
  }

  // strain
  else if (strcmp(argv[0],"strain") == 0) {
    theOutput.tag("ResponseType", "eps11");
    theResponse =  new MaterialResponse(this, 3, this->getStrain());
  }

  // strain
  else if ((strcmp(argv[0],"stressStrain") == 0) || 
	   (strcmp(argv[0],"stressANDstrain") == 0)) {
    theOutput.tag("ResponseType", "sig11");
    theOutput.tag("ResponseType", "eps11");
    theResponse =  new MaterialResponse(this, 4, Vector(2));
  }

  else if (strcmp(argv[0],"stresses") == 0) {
    for (int i=0; i<numMaterials; i++) {
      theOutput.tag("UniaxialMaterialOutput");
      theOutput.attr("matType", this->getClassType());
      theOutput.attr("matTag", this->getTag());
      theOutput.tag("ResponseType", "sigma11");
      theOutput.endTag();
    }
    theResponse =  new MaterialResponse(this, 100, Vector(numMaterials));
  }

  else if (strcmp(argv[0],"material") == 0 ||
	   strcmp(argv[0],"component") == 0) {
    if (argc > 1) {
      int matNum = atoi(argv[1]) - 1;
      if (matNum >= 0 && matNum < numMaterials)
	theResponse =  theModels[matNum]->setResponse(&argv[2], argc-2, theOutput);
    }
  }

  theOutput.endTag();
  return theResponse;
}

int
ParallelMaterial::getResponse(int responseID, Information &info)
{
  Vector stresses(numMaterials);
  int i;

  switch (responseID) {
  case 100:
    for (i = 0; i < numMaterials; i++)
      stresses(i) = theModels[i]->getStress();
    return info.setVector(stresses);

  default:
    return this->UniaxialMaterial::getResponse(responseID, info);
  }
}
