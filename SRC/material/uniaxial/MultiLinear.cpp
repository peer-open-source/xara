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

// Written: fmk 
// Created: 05/12
// Revision: A
//
// Description: This file contains the class implementation for multilinear material. 
//
// What: "@(#) MultiLinear.C, revA"

#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>

#include <MultiLinear.h>
#include <Vector.h>
#include <Matrix.h>
#include <ID.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>

#include <OPS_Globals.h>
#include <elementAPI.h>

void * OPS_ADD_RUNTIME_VPV(OPS_MultiLinear)
{
    // Pointer to a uniaxial material that will be returned
    UniaxialMaterial* theMaterial = 0;
    
    if (OPS_GetNumRemainingInputArgs() < 5) {
        opserr << "Invalid #args,  want: uniaxialMaterial MultiLinear tag? e1 s1 e2 s2 ... " << endln;
        return 0;
    }
    
    int iData[1];
    int numData = 1;
    if (OPS_GetIntInput(&numData, iData) != 0) {
        opserr << "WARNING invalid tag or soilType uniaxialMaterial MultiLinearMaterial" << endln;
        return 0;
    }
    
    numData = OPS_GetNumRemainingInputArgs();
    
    int numSlope = numData / 2;
    double* dData = new double[numData];
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
        opserr << "Invalid pyData data for material uniaxial MultiLinear " << iData[0] << endln;
        return 0;
    }
    
    Vector e(numSlope);
    Vector s(numSlope);
    for (int i = 0; i < numSlope; i++) {
        e(i) = dData[2 * i];
        s(i) = dData[2 * i + 1];
    }
    
    // Parsing was successful, allocate the material
    theMaterial = new MultiLinear(iData[0], s, e);
    if (theMaterial == 0) {
        opserr << "WARNING could not create uniaxialMaterial of type MultiLinear\n";
        return 0;
    }
    
    return theMaterial;
}


MultiLinear::MultiLinear(int tag, const Vector& s, const Vector& e)
    :UniaxialMaterial(tag, MAT_TAG_MultiLinear),
    numSlope(0)
{
    numSlope = e.Size();
    data.resize(numSlope, 6);
    
    for (int i = 1; i < numSlope; i++) {
        if (e(i) < e(i - 1)) {
            opserr << "ERROR: MultiLinear strain_i+1 < strain_i\n";
            exit(-1);
        }
    }
    
    data(0, 0) = -e(0);      // neg yield strain
    data(0, 1) = e(0);       // pos yield strain
    data(0, 2) = -s(0);      // neg yield stress
    data(0, 3) = s(0);       // pos yield stress
    data(0, 4) = s(0) / e(0);  // slope
    data(0, 5) = e(0);       // dist - (0-1)/2
    
    for (int i = 1; i < numSlope; i++) {
        data(i, 0) = -e(i);
        data(i, 1) = e(i);
        data(i, 2) = -s(i);
        data(i, 3) = s(i);
        data(i, 4) = (s(i) - s(i - 1)) / (e(i) - e(i - 1));
        data(i, 5) = e(i) - e(i - 1);
    }
    
    tStrain = 0.0;
    tStress = 0.0;
    tTangent = data(0, 4);
    
    cStrain = 0.0;
    cStress = 0.0;
    cTangent = tTangent;
    
    tSlope = 0;
}


MultiLinear::MultiLinear()
    :UniaxialMaterial(0, MAT_TAG_MultiLinear),
    numSlope(0)
{
    tStrain = 0;
    tStress = 0;
    tTangent = 0;
    cStrain = 0;
    cStress = 0;
    cTangent = 0;
    tSlope = 0;
}


MultiLinear::~MultiLinear()
{
    // does nothing
}

int
MultiLinear::setTrialStrain(double strain, double strainRate)
{
    if (fabs(tStrain - strain) < DBL_EPSILON)
        return 0;
    
    tStrain = strain;
    tSlope = 0;
    
    if (tStrain >= data(0, 0) && tStrain <= data(0, 1)) { // elastic
        tSlope = 0;
        tStress = data(0, 2) + (tStrain - data(0, 0)) * data(0, 4);
        tTangent = data(0, 4);
    }
    else if (tStrain < data(0, 0)) { // search neg of data
        tSlope = 1;
        while (tSlope < numSlope && tStrain < data(tSlope, 0))
            tSlope++;
        if (tSlope == numSlope)
            tSlope = numSlope - 1;
        tStress = data(tSlope, 2) + (tStrain - data(tSlope, 0)) * data(tSlope, 4);
        tTangent = data(tSlope, 4);
    }
    else { // search pos side of data
        tSlope = 1;
        while (tSlope < numSlope && tStrain > data(tSlope, 1))
            tSlope++;
        if (tSlope == numSlope)
            tSlope = numSlope - 1;
        tStress = data(tSlope, 3) + (tStrain - data(tSlope, 1)) * data(tSlope, 4);
        tTangent = data(tSlope, 4);
    }
    
    return 0;
}

double
MultiLinear::getStrain(void)
{
    return tStrain;

}

double
MultiLinear::getStress(void)
{
    return tStress;
}

double
MultiLinear::getTangent(void)
{
    return tTangent;
}

int
MultiLinear::commitState(void)
{
    // if yielded we need to reset the values               
    if (tSlope != 0) { // yielded
        
        if (tStrain > data(0, 1)) {  // positive yield direction
            
            // set elastic bounds
            data(0, 1) = tStrain;
            data(0, 3) = tStress;
            data(0, 0) = tStrain - 2 * data(0, 5);
            data(0, 2) = tStress - 2 * data(0, 5) * data(0, 4);
            
            // reset bounds for all those pts before yield
            //  - pos & neg affected
            for (int i = 1; i < tSlope; i++) {
                data(i, 1) = tStrain;
                data(i, 3) = tStress;
                data(i, 0) = data(i - 1, 0) - 2 * data(i, 5);
                data(i, 2) = data(i - 1, 2) - 2 * data(i, 5) * data(i, 4);
            }
            
            // reset bounds for all those pts after
            //  - neg affected 
            data(tSlope, 0) = data(tSlope - 1, 0) - 2 * data(tSlope, 5)
                            + data(tSlope, 1) - data(tSlope - 1, 1);
            data(tSlope, 2) = data(tSlope - 1, 2)
                            + (data(tSlope, 0) - data(tSlope - 1, 0)) * data(tSlope, 4);
            
            for (int i = tSlope + 1; i < numSlope; i++) {
                data(i, 0) = data(i - 1, 0) - 2 * data(i, 5)
                           + data(i, 1) - data(i - 1, 1);
                data(i, 2) = data(i - 1, 2)
                           + (data(i, 0) - data(i - 1, 0)) * data(i, 4);
            }
        }
        else {  // neg yield direction

            // set elastic bounds
            data(0, 0) = tStrain;
            data(0, 2) = tStress;
            data(0, 1) = tStrain + 2 * data(0, 5);
            data(0, 3) = tStress + 2 * data(0, 5) * data(0, 4);
            
            // reset bounds for all those pts before yield slope
            //  - pos & neg affected 
            for (int i = 1; i < tSlope; i++) {
                data(i, 0) = tStrain;
                data(i, 2) = tStress;
                data(i, 1) = data(i - 1, 1) + 2 * data(i, 5);
                data(i, 3) = data(i - 1, 3) + 2 * data(i, 5) * data(i, 4);
            }
            
            // reset bounds for all those pts after
            //  - pos pts affected 
            data(tSlope, 1) = data(tSlope - 1, 1) + 2 * data(tSlope, 5)
                            + data(tSlope, 0) - data(tSlope - 1, 0);
            data(tSlope, 3) = data(tSlope - 1, 3)
                            + (data(tSlope, 1) - data(tSlope - 1, 1)) * data(tSlope, 4);
            
            for (int i = tSlope + 1; i < numSlope; i++) {
                data(i, 1) = data(i - 1, 1) + 2 * data(i, 5)
                           + data(i, 0) - data(i - 1, 0);
                data(i, 3) = data(i - 1, 3)
                           + (data(i, 1) - data(i - 1, 1)) * data(i, 4);
            }
        }
    }
    
    cStress = tStress;
    cStrain = tStrain;
    cTangent = tTangent;
    
    return 0;
}

int
MultiLinear::revertToLastCommit(void)
{
    tStrain = cStrain;
    tTangent = cTangent;
    tStress = cStress;
    
    return 0;
}

int
MultiLinear::revertToStart(void)
{
    data(0, 1) = data(0, 5);
    data(0, 3) = data(0, 5) * data(0, 4);
    data(0, 0) = -data(0, 1);
    data(0, 2) = -data(0, 3);
    
    for (int i = 1; i < numSlope; i++) {
        data(i, 1) = data(i - 1, 1) + data(i, 5);
        data(i, 3) = data(i - 1, 3) + data(i, 5) * data(i, 4);
        data(i, 0) = -data(i, 1);
        data(i, 2) = -data(i, 3);
    }
    
    tStrain = cStrain = 0.0;
    tStress = cStress = 0.0;
    tTangent = data(0, 4);
    cTangent = tTangent;
    return 0;
}

UniaxialMaterial*
MultiLinear::getCopy(void)
{
    MultiLinear* theCopy =
        new MultiLinear();
    theCopy->data = this->data;
    theCopy->numSlope = this->numSlope;
    
    theCopy->tSlope = this->tSlope;
    theCopy->tStress = this->tStress;
    theCopy->tStrain = this->tStrain;
    theCopy->tTangent = this->tTangent;
    theCopy->cStress = this->cStress;
    theCopy->cStrain = this->cStrain;
    theCopy->cTangent = this->cTangent;
    
    return theCopy;
}

int
MultiLinear::sendSelf(int cTag, Channel& theChannel)
{
    int res = 0;
    static ID iData(2);
    iData(0) = this->getTag();
    iData(1) = numSlope;
    
    res = theChannel.sendID(this->getDbTag(), cTag, iData);
    if (res < 0) {
        opserr << "ElasticMaterial::sendSelf() - failed to send data\n";
        return res;
    }
    
    res = theChannel.sendMatrix(this->getDbTag(), cTag, data);
    if (res < 0)
        opserr << "ElasticMaterial::sendSelf() - failed to send data\n";
    
    return res;
}

int
MultiLinear::recvSelf(int cTag, Channel& theChannel,
    FEM_ObjectBroker& theBroker)
{
    int res = 0;
    static ID iData(2);
    res = theChannel.recvID(this->getDbTag(), cTag, iData);
    if (res < 0) {
        opserr << "ElasticMaterial::recvSelf() - failed to recv data\n";
        return res;
    }
    
    this->setTag(iData(0));
    numSlope = iData(1);
    
    data.resize(numSlope, 6);
    
    res = theChannel.recvMatrix(this->getDbTag(), cTag, data);
    if (res < 0)
        opserr << "ElasticMaterial::recvSelf() - failed to recv data\n";
    
    return res;
}

void
MultiLinear::Print(OPS_Stream& s, int flag)
{
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << OPS_PRINT_JSON_MATE_INDENT << "{";
        s << "\"type\": \"MultiLinear\", ";
        s << "\"tag\": " << this->getTag();
        s << "}";
    } else {
        s << "MultiLinear tag: " << this->getTag() << endln;
        s << "  stress: " << tStress << " tangent: " << tTangent << endln;
        s << "tSlope: " << tSlope << "numSlope: " << numSlope << endln;
        s << data;
    }
}

// AddingSensitivity:BEGIN ///////////////////////////////////
int
MultiLinear::setParameter(const char **argv, int argc, Parameter &param)
{
    // trying to make this a little more general for arbitrary numbers of points
    // expecting parameters named stress1, strain1, stress3, etc.
    int dindx = -1;

    if (argc < 2) {
      return -1;
    }
    
    dindx = atoi(argv[1]);
    if (dindx < 1 || dindx > numSlope) {
      return -1;
    }
    
    if (strcmp(argv[0],"stress") == 0) {
        param.setValue( data(dindx-1,3) );
        return param.addObject(100+dindx, this);
    }
    if (strcmp(argv[0],"strain") == 0) {
        param.setValue( data(dindx-1,1) );
        return param.addObject(200+dindx, this);
    }

    return -1;
}


int
MultiLinear::updateParameter(int parameterID, Information &info)
{
    int dindx = 0;
    double sprev = 0;
    double eprev = 0;
    
    if (parameterID == -1)
        return -1;
    else if (parameterID > 100 && parameterID <= 100+numSlope) {
        // stress terms
        dindx = parameterID-100;
        data(dindx-1,2) = -1.0*info.theDouble;
        data(dindx-1,3) = info.theDouble;
    }
    else if (parameterID > 200 && parameterID <= 200+numSlope) {
        // strain terms
        dindx = parameterID-200;
        data(dindx-1,0) = -1.0*info.theDouble;
        data(dindx-1,1) = info.theDouble;
    }
    else {
        return -1;
    }
    
    // update slopes
    if (dindx > 1) {
        sprev = data(dindx-2,3);
        eprev = data(dindx-2,1);
    }
    data(dindx-1,4) = (data(dindx-1,3) - sprev) / (data(dindx-1,1) - eprev);
    data(dindx-1,5) = data(dindx-1,1) - eprev;
    
    return 0;
}
