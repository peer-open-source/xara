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

#include <ShellFiberSection.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <SensitiveResponse.h>

typedef SensitiveResponse<SectionForceDeformation> SectionResponse;

#include <cassert>
#include <cstdlib>
#include <cstring>

#define SEC_TAG_ShellFiberSection 0

Vector ShellFiberSection::stressResultant(8);
Matrix ShellFiberSection::tangent(8, 8);
ID ShellFiberSection::array(8);

ShellFiberSection::ShellFiberSection()
  : SectionForceDeformation(0, SEC_TAG_ShellFiberSection),
    fibers(),
    h(0.0),
    strainResultant(8)
{

}

ShellFiberSection::ShellFiberSection(int tag,
                                     const FiberVector& inputFibers)
  : SectionForceDeformation(tag, SEC_TAG_ShellFiberSection),
    fibers(),
    h(0.0),
    strainResultant(8)
{
  assert(!inputFibers.empty());

  fibers.reserve(inputFibers.size());
  for (const auto& [z, w, material] : inputFibers)
    copyFiber(z, w, material);

  assert(!fibers.empty());
  assert(h > 0.0);
}

ShellFiberSection::ShellFiberSection(int tag,
                                     int nLayers,
                                     double* thickness,
                                     NDMaterial** materials)
  : SectionForceDeformation(tag, SEC_TAG_ShellFiberSection),
    fibers(),
    h(0.0),
    strainResultant(8)
{
  assert(nLayers > 0);
  assert(thickness != nullptr);
  assert(materials != nullptr);

  double totalThickness = 0.0;
  for (int i = 0; i < nLayers; ++i) {
    assert(thickness[i] > 0.0);
    totalThickness += thickness[i];
  }

  assert(totalThickness > 0.0);

  fibers.reserve(static_cast<std::size_t>(nLayers));

  double z0 = -0.5 * totalThickness;
  for (int i = 0; i < nLayers; ++i) {
    const double w = thickness[i];
    const double z = z0 + 0.5 * w;

    copyFiber(z, w, materials[i]);
    z0 += w;
  }

  assert(h > 0.0);
}


ShellFiberSection::ShellFiberSection(const ShellFiberSection& other)
  : SectionForceDeformation(other.getTag(), SEC_TAG_ShellFiberSection),
    fibers(),
    h(0.0),
    strainResultant(8)
{
  strainResultant = other.strainResultant;

  fibers.reserve(other.fibers.size());
  for (const auto& [z, w, material] : other.fibers)
    copyFiber(z, w, material);
}


ShellFiberSection::~ShellFiberSection()
{
  clearFibers();
}

void
ShellFiberSection::clearFibers()
{
  for (auto& fiber : fibers) {
    delete std::get<2>(fiber);
    std::get<2>(fiber) = nullptr;
  }

  fibers.clear();
  h = 0.0;
}

void
ShellFiberSection::copyFiber(double z, double w, NDMaterial* material)
{
  assert(w > 0.0);
  assert(material != nullptr);

  NDMaterial* copy = material->getCopy("PlateFiber");
  assert(copy != nullptr);

  fibers.emplace_back(z, w, copy);
  h += w;
}

int
ShellFiberSection::getNumFibers() const
{
  return static_cast<int>(fibers.size());
}

SectionForceDeformation*
ShellFiberSection::getCopy()
{
  return new ShellFiberSection(*this);
}

int
ShellFiberSection::getOrder() const
{
  return 8;
}

const ID&
ShellFiberSection::getType()
{
  static bool initialized = false;
  if (!initialized) {
    array(0)    = SECTION_RESPONSE_FXX;
    array(1)    = SECTION_RESPONSE_FYY;
    array(2)    = SECTION_RESPONSE_FXY;
    array(3)    = SECTION_RESPONSE_MXX;
    array(4)    = SECTION_RESPONSE_MYY;
    array(5)    = SECTION_RESPONSE_MXY;
    array(6)    = SECTION_RESPONSE_VXZ;
    array(7)    = SECTION_RESPONSE_VYZ;
    initialized = true;
  }
  return array;
}

int
ShellFiberSection::commitState()
{
  int success = 0;

  for (auto& [z, w, material] : fibers)
    success += material->commitState();

  return success;
}

int
ShellFiberSection::revertToLastCommit()
{
  int success = 0;

  for (auto& [z, w, material] : fibers)
    success += material->revertToLastCommit();

  return success;
}

int
ShellFiberSection::revertToStart()
{
  int success = 0;

  for (auto& [z, w, material] : fibers)
    success += material->revertToStart();

  return success;
}

double
ShellFiberSection::getRho()
{
  double rhoH = 0.0;

  for (const auto& [z, w, material] : fibers)
    rhoH += material->getRho() * w;

  return rhoH;
}

Response*
ShellFiberSection::setResponse(const char** argv, int argc, OPS_Stream& output)
{
  Response* theResponse = nullptr;

  if (strcmp(argv[0], "fiber") == 0 || strcmp(argv[0], "Fiber") == 0) {
    if (argc < 3) {
      opserr << "ShellFiberSection::setResponse() - need to specify more data\n";
      return nullptr;
    }

    const int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= getNumFibers()) {
      const auto& [z, w, material] = fibers[static_cast<std::size_t>(pointNum - 1)];

      output.tag("FiberOutput");
      output.attr("number", pointNum);
      output.attr("zLoc", z);
      output.attr("weight", w);

      theResponse = material->setResponse(&argv[2], argc - 2, output);

      output.endTag();
    }
  }

  if (theResponse == nullptr)
    return SectionForceDeformation::setResponse(argv, argc, output);

  return theResponse;
}

int
ShellFiberSection::getResponse(int responseID, Information& secInfo)
{
  return SectionForceDeformation::getResponse(responseID, secInfo);
}

int
ShellFiberSection::setTrialSectionDeformation(const Vector& strainResultant_from_element)
{
  this->strainResultant = strainResultant_from_element;

  static Vector strain(5);

  int success = 0;

  for (auto& [z, w, material] : fibers) {
    strain(0) = strainResultant(0) - z * strainResultant(3);
    strain(1) = strainResultant(1) - z * strainResultant(4);
    strain(2) = strainResultant(2) - z * strainResultant(5);
    strain(3) = strainResultant(6);
    strain(4) = strainResultant(7);

    success += material->setTrialStrain(strain);
  }

  return success;
}

const Vector&
ShellFiberSection::getSectionDeformation()
{
  return this->strainResultant;
}

const Vector&
ShellFiberSection::getStressResultant()
{
  stressResultant.Zero();

  for (const auto& [z, w, material] : fibers) {
    const Vector& stress = material->getStress();

    stressResultant(0) += stress(0) * w;
    stressResultant(1) += stress(1) * w;
    stressResultant(2) += stress(2) * w;

    stressResultant(3) += z * stress(0) * w;
    stressResultant(4) += z * stress(1) * w;
    stressResultant(5) += z * stress(2) * w;

    stressResultant(6) += stress(3) * w;
    stressResultant(7) += stress(4) * w;
  }

  return this->stressResultant;
}

const Matrix&
ShellFiberSection::getSectionTangent()
{
  static constexpr int stressIndex[8] = {0, 1, 2, 0, 1, 2, 3, 4};
  static constexpr int strainIndex[8] = {0, 1, 2, 0, 1, 2, 3, 4};

  tangent.Zero();

  for (const auto& [z, w, material] : fibers) {
    const Matrix& dd = material->getTangent();

    const double stressFactor[8] = {1.0, 1.0, 1.0, z, z, z, 1.0, 1.0};
    const double strainFactor[8] = {1.0, 1.0, 1.0, -z, -z, -z, 1.0, 1.0};

    for (int i = 0; i < 8; ++i) {
      const int ii = stressIndex[i];
      const double si = stressFactor[i];

      for (int j = 0; j < 8; ++j) {
        const int jj = strainIndex[j];
        tangent(i, j) += si * dd(ii, jj) * strainFactor[j] * w;
      }
    }
  }

  return this->tangent;
}


void
ShellFiberSection::Print(OPS_Stream& s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_SECTION || flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
    s << "ShellFiber Section tag: " << this->getTag() << endln;
    s << "Total thickness h = " << h << endln;

    for (std::size_t i = 0; i < fibers.size(); ++i) {
      const auto& [z, w, material] = fibers[i];

      s << "Fiber " << i + 1
        << ", z = " << z
        << ", weight = " << w << endln;
      material->Print(s, flag);
      s << endln;
    }
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_MATE_INDENT << "{";
    s << "\"name\": \"" << this->getTag() << "\", ";
    s << "\"type\": \"ShellFiberSection\", ";
    s << "\"totalThickness\": " << h << ", ";
    s << "\"fibers\": [\n";

    for (std::size_t i = 0; i < fibers.size(); ++i) {
      const auto& [z, w, material] = fibers[i];

      s << OPS_PRINT_JSON_MATE_INDENT << "  {\"fiber\": " << i + 1 << ", ";
      s << "\"zLoc\": " << z << ", ";
      s << "\"weight\": " << w << ", ";
      s << "\"material\": \"" << material->getTag() << "\"";

      if (i + 1 < fibers.size())
        s << "},\n";
      else
        s << "}\n";
    }

    s << OPS_PRINT_JSON_MATE_INDENT << "]}";
  }
}

int
ShellFiberSection::sendSelf(int commitTag, Channel& theChannel)
{
  int res = 0;
  const int dataTag = this->getDbTag();
  const int nFibers = getNumFibers();

  ID iData(2);
  iData(0) = this->getTag();
  iData(1) = nFibers;

  res += theChannel.sendID(dataTag, commitTag, iData);
  if (res < 0) {
    opserr << "WARNING ShellFiberSection::sendSelf() - " << this->getTag()
           << " failed to send data" << endln;
    return res;
  }

  if (nFibers == 0)
    return res;

  Vector vecData(2 * nFibers);
  for (int i = 0; i < nFibers; ++i) {
    const auto& [z, w, material] = fibers[static_cast<std::size_t>(i)];
    vecData(i) = z;
    vecData(i + nFibers) = w;
  }

  res += theChannel.sendVector(dataTag, commitTag, vecData);
  if (res < 0) {
    opserr << "WARNING ShellFiberSection::sendSelf() - " << this->getTag()
           << " failed to send fiber data" << endln;
    return res;
  }

  ID idData(2 * nFibers);
  for (int i = 0; i < nFibers; ++i) {
    auto& [z, w, material] = fibers[static_cast<std::size_t>(i)];

    idData(i) = material->getClassTag();

    int matDbTag = material->getDbTag();
    if (matDbTag == 0) {
      matDbTag = theChannel.getDbTag();
      if (matDbTag != 0)
        material->setDbTag(matDbTag);
    }

    idData(i + nFibers) = matDbTag;
  }

  res += theChannel.sendID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING ShellFiberSection::sendSelf() - " << this->getTag()
           << " failed to send material IDs" << endln;
    return res;
  }

  for (auto& [z, w, material] : fibers) {
    res += material->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING ShellFiberSection::sendSelf() - " << this->getTag()
             << " failed to send its material" << endln;
      return res;
    }
  }

  return res;
}

int
ShellFiberSection::recvSelf(int commitTag,
                            Channel& theChannel,
                            FEM_ObjectBroker& theBroker)
{
  int res = 0;
  const int dataTag = this->getDbTag();

  ID iData(2);
  res += theChannel.recvID(dataTag, commitTag, iData);
  if (res < 0) {
    opserr << "WARNING ShellFiberSection::recvSelf() - " << this->getTag()
           << " failed to receive data" << endln;
    return res;
  }

  this->setTag(iData(0));
  const int nFibers = iData(1);

  if (nFibers == 0) {
    clearFibers();
    return res;
  }

  Vector vecData(2 * nFibers);
  res += theChannel.recvVector(dataTag, commitTag, vecData);
  if (res < 0) {
    opserr << "WARNING ShellFiberSection::recvSelf() - " << this->getTag()
           << " failed to receive fiber data" << endln;
    return res;
  }

  ID idData(2 * nFibers);
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING ShellFiberSection::recvSelf() - " << this->getTag()
           << " failed to receive material IDs" << endln;
    return res;
  }

  FiberVector newFibers;
  newFibers.reserve(static_cast<std::size_t>(nFibers));
  double newH = 0.0;

  for (int i = 0; i < nFibers; ++i) {
    const double z = vecData(i);
    const double w = vecData(i + nFibers);
    const int matClassTag = idData(i);

    NDMaterial* material = theBroker.getNewNDMaterial(matClassTag);
    if (material == nullptr) {
      opserr << "ShellFiberSection::recvSelf() - broker could not create NDMaterial of class type "
             << matClassTag << endln;

      for (auto& fiber : newFibers)
        delete std::get<2>(fiber);

      return -1;
    }

    material->setDbTag(idData(i + nFibers));

    res += material->recvSelf(commitTag, theChannel, theBroker);
    if (res < 0) {
      opserr << "ShellFiberSection::recvSelf() - material " << i
             << " failed to receive itself" << endln;

      delete material;
      for (auto& fiber : newFibers)
        delete std::get<2>(fiber);

      return res;
    }

    newFibers.emplace_back(z, w, material);
    newH += w;
  }

  clearFibers();
  fibers.swap(newFibers);
  h = newH;

  return res;
}
