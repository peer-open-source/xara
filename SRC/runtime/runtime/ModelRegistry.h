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
// Description: This file contains the class definition for
// ModelRegistry. A ModelRegistry aims to be a threadsafe
// alternative to the TclBasicBuilder class. 
//
// Written: cmp
// Created: 10/21
//
#pragma once

#include <typeinfo>
#include <string>
#include <vector>
#include <unordered_map>
#include <LoadCase.h>
#include <TaggedObject.h>
#include <StaticPattern.h>
#include <MultiSupportPattern.h>
#include <Rotations.h>
#include "ProcessContext.h"

class LoadPattern;
class StaticPattern;
class MultiSupportPattern;
class OPS_Stream;
class ID;
class Domain;
class InterpreterResponse;
class UniaxialMaterial;
class NDMaterial;

using Xara::ProcessContext;

class ModelRegistry {
public:

  ModelRegistry(Domain &domain, int ndm, int ndf, Rotations::Parameters);
  ~ModelRegistry();


  int getNDM() const;
  int getNDF() const;
  void setDimension(int ndm, int ndf) {this->ndm = ndm; this->ndf = ndf;}
  Domain *getDomain() const;
#ifdef MODEL_CHANNELS
  ProcessContext& getParallelContext() {return m_channels;}
#endif
  //
  // Managing tagged objects
  //
  template<class T> int addTypedObject(int tag, T* obj) {
    auto status = addRegistryObject(typeid(T).name(), nullptr, tag, obj);
    if constexpr (std::is_same<T, NDMaterial>::value) {
      if (status == TCL_OK)
        obj->setDomain(theDomain);
    }
    return status;
  }

  template<class T, const char* specialize=nullptr> int addTaggedObject(T& obj) {
    int tag = obj.getTag();
    auto status =  addRegistryObject(typeid(T).name(), specialize, tag, &obj);

    //
    if constexpr (std::is_same<T, NDMaterial>::value) {
      if (status == TCL_OK)
        obj.setDomain(theDomain);
    }
    // else if constexpr (std::is_same<T, UniaxialMaterial>::value) {
    //   if (status == 0)
    //     obj.setDomain(theDomain);
    // }
    return status;
  }

  constexpr static int SilentLookup = 1;
  template <class T>
  int printRegistry(OPS_Stream& stream, int flag) const {
    return printRegistry(typeid(T).name(), stream, flag);
  }

  template<class T, const char* specialize=nullptr> T* 
  getTypedObject(int tag, int flags=0) const {
    return (T*)getRegistryObject(typeid(T).name(), specialize, tag, flags);
  }

  template<class T> int 
  removeObject(int tag, int flags=0) {
    return removeRegistryObject(typeid(T).name(), tag, flags);
  }

  template <class T> int findFreeTag(int &tag) const {
    return findFreeTag(typeid(T).name(), tag);
  }


  int  getCurrentSectionBuilder(int&);
  void setCurrentSectionBuilder(int);

  int addResponse(InterpreterResponse* response) {
    m_responses.push_back(response);
    return m_responses.size() - 1;
  }

  InterpreterResponse* getResponse(int index) {
    if (index < 0 || index >= std::ssize(m_responses))
      return nullptr;
    return m_responses[index];
  }

  OpenSees::LoadCase& getLoadCase();
  int setLoadCase(std::string& name);
  int newLoadCase(std::string& name);

  template<class T>
  T* getCurrentPattern() {
    if constexpr (std::is_same<T, StaticPattern>::value) {
      return static_pattern;
    }
    else if constexpr (std::is_same<T, MultiSupportPattern>::value) {
      return multi_pattern;
    }
    else if constexpr (std::is_same<T, LoadPattern>::value) {
      return tclEnclosingPattern;
    }
  }

  template <class T>
  int setCurrentPattern(T* pattern) {
    if constexpr (std::is_same<T, StaticPattern>::value) {
      static_pattern = pattern;
    }
    else if constexpr (std::is_same<T, MultiSupportPattern>::value) {
      multi_pattern = pattern;
    }
    tclEnclosingPattern = static_cast<LoadPattern*>(pattern);
    return 0;
  }

  Rotations::Parameters getRotationType() const {return rotation_type;}

  int incrNodalLoadTag();
  int decrNodalLoadTag();
  int getNodalLoadTag();
  int incrElemLoadTag();
  int decrElemLoadTag();
  int getElemLoadTag();

  int addSP_Constraint(int axisDirn, 
         double axisValue, 
         const ID &fixityCodes, 
         double tol=1e-10);

//
private:
  int   addRegistryObject(const char*, const char*, int tag, void* obj); 
  void* getRegistryObject(const char*, const char*, int tag, int flags) const;
  int   removeRegistryObject(const char*, int tag, int flags);
  int   findFreeTag(const char*, int& tag) const;
  int   printRegistry(const char *, OPS_Stream& stream, int flag) const ;


  int ndm; // space dimension of the mesh
  int ndf; // number of degrees of freedom per node

  Domain *theDomain     = nullptr;

  int next_node_load          = 0;
  int next_elem_load          = 0;

  // previously extern variables
  LoadPattern *tclEnclosingPattern = nullptr;
  StaticPattern* static_pattern       = nullptr;
  MultiSupportPattern* multi_pattern  = nullptr;

  bool  section_builder_is_set   = false;
  int   current_section_builder  = 0;

  Rotations::Parameters rotation_type;

// OBJECT CONTAINERS
  std::unordered_map<std::string, std::unordered_map<int, TaggedObject*>> m_registry;
  std::unordered_map<std::string, OpenSees::LoadCase> m_cases;
  std::vector<InterpreterResponse*> m_responses;

  // Parallel
#ifdef MODEL_CHANNELS
  ProcessContext m_channels;
#endif
};

