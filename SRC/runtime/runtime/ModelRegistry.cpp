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
// A ModelRegistry stores intermediate "reference" objects like
// materials and sections that are used
// to construct other objects like elements.
//
// Written: cmp
//
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <initializer_list>
#include <string>
#include <unordered_map>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Domain.h>

#include <ModelRegistry.h>
#include <Logging.h>
#include <Parsing.h> // For TCL_OK/ERROR

using OpenSees::LoadCase;

ModelRegistry::ModelRegistry(Domain &domain,
                             int NDM, int NDF, 
                             Rotations::Parameters rotation_type)
  : ndm(NDM), ndf(NDF),
    rotation_type(rotation_type),
    section_builder_is_set(false),
    theDomain(&domain),
    tclEnclosingPattern(nullptr),
    next_node_load(0)
{

}

ModelRegistry::~ModelRegistry()
{

  for (auto& [part, val] : m_registry ) {
    for (auto& [tag, obj] : val)
      delete obj;
  }

  // set the pointers to 0
  theDomain = nullptr;
  tclEnclosingPattern = nullptr;
}


int
ModelRegistry::buildFE_Model()
{
  return 0;
}


int
ModelRegistry::getNDM() const
{
  return ndm;
}

int
ModelRegistry::getNDF() const
{
  return ndf;
}


int
ModelRegistry::incrNodalLoadTag()
{
  return ++next_node_load;
}

int
ModelRegistry::decrNodalLoadTag()
{
  return --next_node_load;
}

int
ModelRegistry::getNodalLoadTag() 
{
  return next_node_load;
}

int
ModelRegistry::incrElemLoadTag()
{
  return ++next_elem_load;
}

int
ModelRegistry::decrElemLoadTag()
{
  return --next_elem_load;
}

int
ModelRegistry::getElemLoadTag()
{
  return next_elem_load++;
}

int
ModelRegistry::addSP_Constraint(int axisDirn, double axisValue, const ID &fixityCodes, double tol)
{
  return theDomain->addSP_Constraint(axisDirn, axisValue, fixityCodes, tol);
}


int
ModelRegistry::getCurrentSectionBuilder(int& tag)
{
  if (section_builder_is_set) {
    tag = current_section_builder;
    return  0;
  } else
    return -1;
}

void 
ModelRegistry::setCurrentSectionBuilder(int tag)
{
  section_builder_is_set   = true;
  current_section_builder  = tag;
}

Domain *
ModelRegistry::getDomain() const 
{
  return theDomain;
}


int 
ModelRegistry::printRegistry(const char *partition, OPS_Stream& stream, int flag) const 
{
  int count = 0;
  auto iter = m_registry.find(partition);
  if (iter == m_registry.end()) {
    return count;
  }

  for (auto const& [key, val] : iter->second) {
    if (count != 0)
      stream << ",\n";

    val->Print(stream, flag);
    count++;
  }

  return count;
}

void* 
ModelRegistry::getRegistryObject(const char* type, const char* specialize, int tag, int flags) const
{
  std::string partition = std::string{type};
  if (specialize)
    partition += std::string{specialize};

  auto iter = m_registry.find(partition);
  if (iter == m_registry.end()) {
    if (flags == 0)
      opserr << "No objects of type \"" << partition.c_str()
             << "\" have been created.\n";
    return nullptr;
  }

  auto iter_objs = iter->second.find(tag) ;
  if (iter_objs == iter->second.end()) {
    if (flags == 0)
      opserr << "No object with tag \"" << tag << "\" of type \"" 
             << partition.c_str() << "\"\n";
    return nullptr;
  }

  return (void*)iter_objs->second;
}

int
ModelRegistry::addRegistryObject(const char* type, const char* specialize, int tag, void *obj)
{
  std::string partition = std::string{type};
  if (specialize)
    partition += std::string{specialize};

  // check for clobbering an existing object
  auto iter = m_registry.find(partition);
  if (iter != m_registry.end()) {
    auto iter_objs = iter->second.find(tag) ;
    if (iter_objs != iter->second.end()) {
      opserr << OpenSees::PromptValueError 
             << "An object with tag \"" << tag << "\" of type \"" 
             << partition.c_str() << "\" already exists.\n";
      return TCL_ERROR;
    }
  }

  m_registry[partition][tag] = (TaggedObject*)obj;
  return TCL_OK;
}

int
ModelRegistry::findFreeTag(const char* partition, int& tag) const
{
  tag = 0;
  const auto iter = m_registry.find(std::string{partition});
  // If we dont have a table with partition name, no objects
  // have been created and tag = 0 works; return success.
  if (iter == m_registry.end())
    return 0;


  // Otherwise, find something larger than all existing tags
  const std::unordered_map<int, TaggedObject*>& table = iter->second;
  for (auto const& [key, val] : table)
    if (key > tag)
      tag = key + 1;

  return 0;
}

int
ModelRegistry::removeRegistryObject(const char* partition, int tag, int flags) 
{
  const auto iter = m_registry.find(std::string{partition});
  if (iter == m_registry.end()) {
    if (flags == 0)
      opserr << "No objects of type \"" << partition
             << "\" have been created.\n";
    return -1;
  }
  std::unordered_map<int, TaggedObject*>& table = iter->second;
  if (table.find(tag) != table.end()) {
      table.erase(tag);
      return 0;
  }
  return -1;
}


