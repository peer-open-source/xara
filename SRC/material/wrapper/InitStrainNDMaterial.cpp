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
// Written: Massimo Petracca - ASDEA Software, Italy
// Created: 03/2024
//
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "InitStrainNDMaterial.h"
#include <Matrix.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <string.h>
#include <Parameter.h>

namespace {

    // a base index to try not to shadow the adpted material's parameters...
    constexpr int param_base_index = 111000;

}


InitStrainNDMaterial::InitStrainNDMaterial(int tag, NDMaterial& material, const Vector& eps0)
    : NDMaterial(tag, ND_TAG_InitStrainNDMaterial)
    , theMaterial(nullptr)
    , epsInit(eps0)
{
    // get copy of the main material
    theMaterial = material.getCopy("ThreeDimensional");
    if (theMaterial == nullptr) {
        opserr << "InitStrainNDMaterial::InitStrainNDMaterial -- failed to get copy of material (a 3D material is required)\n";
        exit(-1);
    }

    // make sure the input strain vector is of correct size
    assert(epsInit.Size() == 6);
}

InitStrainNDMaterial::InitStrainNDMaterial(int tag, NDMaterial& material, double eps0)
    : NDMaterial(tag, ND_TAG_InitStrainNDMaterial)
    , theMaterial(0)
{
    // get copy of the main material
    theMaterial = material.getCopy("ThreeDimensional");
    if (theMaterial == nullptr) {
        opserr << "InitStrainNDMaterial::InitStrainNDMaterial -- failed to get copy of material (a 3D material is required)\n";
        exit(-1);
    }

    // initialize epsInit
    epsInit.resize(6);
    epsInit.Zero();
    for (int i = 0; i < 3; ++i)
        epsInit(i) = eps0;
}

InitStrainNDMaterial::InitStrainNDMaterial()
    : NDMaterial(0, ND_TAG_InitStrainNDMaterial)
    , theMaterial(nullptr)
    , epsInit(6)
{

}

InitStrainNDMaterial::~InitStrainNDMaterial()
{
    if (theMaterial)
        delete theMaterial;
}

int
InitStrainNDMaterial::setTrialStrain(const Vector& strain)
{
    static Vector total_strain(6);
    total_strain = strain;
    total_strain.addVector(1.0, epsInit, 1.0);
    return theMaterial->setTrialStrain(total_strain);
}

int
InitStrainNDMaterial::setTrialStrain(const Vector& strain, const Vector& /*strainRate*/)
{
    return setTrialStrain(strain);
}

int
InitStrainNDMaterial::setTrialStrainIncr(const Vector& strain)
{
    static Vector strain_from_ele(6);
    strain_from_ele = theMaterial->getStrain();
    strain_from_ele.addVector(1.0, epsInit, -1.0);
    strain_from_ele.addVector(1.0, strain, 1.0);
    return setTrialStrain(strain_from_ele);
}

int
InitStrainNDMaterial::setTrialStrainIncr(const Vector& strain, const Vector& /*strainRate*/)
{
    return setTrialStrainIncr(strain);
}

const Vector&
InitStrainNDMaterial::getStress()
{
    return theMaterial->getStress();
}

const Matrix&
InitStrainNDMaterial::getTangent()
{
    return theMaterial->getTangent();
}

const Matrix&
InitStrainNDMaterial::getInitialTangent()
{
    return theMaterial->getInitialTangent();
}

const Vector&
InitStrainNDMaterial::getStrain()
{
    return theMaterial->getStrain();
}

int
InitStrainNDMaterial::commitState()
{
    return theMaterial->commitState();
}

int
InitStrainNDMaterial::revertToLastCommit()
{
    return theMaterial->revertToLastCommit();
}

int
InitStrainNDMaterial::revertToStart()
{
    return theMaterial->revertToStart();
}

double
InitStrainNDMaterial::getRho()
{
    return theMaterial->getRho();
}

NDMaterial*
InitStrainNDMaterial::getCopy()
{
    InitStrainNDMaterial* theCopy = new InitStrainNDMaterial(getTag(), *theMaterial, epsInit);
    return theCopy;
}

NDMaterial *
InitStrainNDMaterial::getCopy(const char *type)
{
    if (strcmp(type, "ThreeDimensional") == 0)
        return getCopy();
    return NDMaterial::getCopy(type);
}

const char*
InitStrainNDMaterial::getType() const
{
    return theMaterial->getType();
}

int InitStrainNDMaterial::getOrder() const
{
    return 6;
}

int
InitStrainNDMaterial::sendSelf(int cTag, Channel& theChannel)
{
    if (theMaterial == nullptr) {
        opserr << "InitStrainNDMaterial::sendSelf() - theMaterial is null, nothing to send" << endln;
        return -1;
    }

    int dbTag = this->getDbTag();

    static ID dataID(3);
    dataID(0) = this->getTag();
    dataID(1) = theMaterial->getClassTag();
    int matDbTag = theMaterial->getDbTag();
    if (matDbTag == 0) {
        matDbTag = theChannel.getDbTag();
        theMaterial->setDbTag(matDbTag);
    }
    dataID(2) = matDbTag;
    if (theChannel.sendID(dbTag, cTag, dataID) < 0) {
        opserr << "InitStrainNDMaterial::sendSelf() - failed to send the ID\n";
        return -1;
    }

    if (theChannel.sendVector(dbTag, cTag, epsInit) < 0) {
        opserr << "InitStrainNDMaterial::sendSelf() - failed to send epsInit\n";
        return -2;
    }

    if (theMaterial->sendSelf(cTag, theChannel) < 0) {
        opserr << "InitStrainNDMaterial::sendSelf() - failed to send the Material\n";
        return -3;
    }

    return 0;
}

int
InitStrainNDMaterial::recvSelf(int cTag, Channel& theChannel,
    FEM_ObjectBroker& theBroker)
{
    int dbTag = this->getDbTag();

    static ID dataID(3);
    if (theChannel.recvID(dbTag, cTag, dataID) < 0) {
        opserr << "InitStrainNDMaterial::recvSelf() - failed to get the ID\n";
        return -1;
    }
    setTag(dataID(0));

    // as no way to change material, don't have to check classTag of the material 
    if (theMaterial == 0) {
        int matClassTag = dataID(1);
        theMaterial = theBroker.getNewNDMaterial(matClassTag);
        if (theMaterial == 0) {
            opserr << "InitStrainNDMaterial::recvSelf() - failed to create Material with classTag "
                << matClassTag << endln;
            return -2;
        }
    }
    theMaterial->setDbTag(dataID(2));

    epsInit.resize(6);
    if (theChannel.recvVector(dbTag, cTag, epsInit) < 0) {
        opserr << "InitStrainNDMaterial::recvSelf() - failed to get the epsInit vector\n";
        return -3;
    }

    if (theMaterial->recvSelf(cTag, theChannel, theBroker) < 0) {
        opserr << "InitStrainNDMaterial::recvSelf() - failed to get the Material\n";
        return -4;
    }

    return 0;
}

void
InitStrainNDMaterial::Print(OPS_Stream& s, int flag)
{

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << OPS_PRINT_JSON_MATE_INDENT << "{";
        s << "\"name\": \"" << this->getTag() << "\", ";
        s << "\"type\": \"InitialStrain\", ";
        if (theMaterial)
          s << "\"Material\": " << theMaterial->getTag() << ", ";
        else
          s << "\"Material\": " << "null" << ", ";
        s << "\"initial_strain\": " << epsInit;
        s <<  "}";
        return;
    }
    else {
        s << "InitStrainNDMaterial tag: " << this->getTag() << endln;
        s << "\tMaterial: " << theMaterial->getTag() << endln;
        s << "\tinitital strain: " << epsInit << endln;
    }
}

int
InitStrainNDMaterial::setParameter(const char** argv, int argc, Parameter& param)
{
    if (argc > 0) {
        if (strcmp(argv[0], "eps0") == 0) {
            // eps0 is assumed to impose with a single scalar a volumetric stress = I*eps0
            param.setValue(epsInit(0));
            return param.addObject(param_base_index, this);
        }
        else if (strcmp(argv[0], "eps0_11") == 0) {
            param.setValue(epsInit(0));
            return param.addObject(param_base_index + 1, this);
        }
        else if (strcmp(argv[0], "eps0_22") == 0) {
            param.setValue(epsInit(1));
            return param.addObject(param_base_index + 2, this);
        }
        else if (strcmp(argv[0], "eps0_33") == 0) {
            param.setValue(epsInit(2));
            return param.addObject(param_base_index + 3, this);
        }
        else if (strcmp(argv[0], "eps0_12") == 0) {
            param.setValue(epsInit(3));
            return param.addObject(param_base_index + 4, this);
        }
        else if (strcmp(argv[0], "eps0_23") == 0) {
            param.setValue(epsInit(4));
            return param.addObject(param_base_index + 5, this);
        }
        else if (strcmp(argv[0], "eps0_13") == 0) {
            param.setValue(epsInit(5));
            return param.addObject(param_base_index + 6, this);
        }
    }
    return theMaterial->setParameter(argv, argc, param);
}

int InitStrainNDMaterial::updateParameter(int parameterID, Information& info)
{
    switch (parameterID) {
    case param_base_index:
        epsInit(0) = epsInit(1) = epsInit(2) = info.theDouble;
        return 0;
    case param_base_index + 1:
        epsInit(0) = info.theDouble;
        return 0;
    case param_base_index + 2:
        epsInit(1) = info.theDouble;
        return 0;
    case param_base_index + 3:
        epsInit(2) = info.theDouble;
        return 0;
    case param_base_index + 4:
        epsInit(3) = info.theDouble;
        return 0;
    case param_base_index + 5:
        epsInit(4) = info.theDouble;
        return 0;
    case param_base_index + 6:
        epsInit(5) = info.theDouble;
        return 0;
    default:
        return -1;
    }
}

Response* InitStrainNDMaterial::setResponse(const char** argv, int argc, OPS_Stream& output)
{
    return theMaterial->setResponse(argv, argc, output);
}

const Vector&
InitStrainNDMaterial::getStressSensitivity(int gradIndex, bool conditional)
{
    return theMaterial->getStressSensitivity(gradIndex, conditional);
}

int
InitStrainNDMaterial::commitSensitivity(const Vector& depsdh, int gradIndex, int numGrads)
{
    return theMaterial->commitSensitivity(depsdh, gradIndex, numGrads);
}
