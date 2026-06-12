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
#pragma once
#include "FiniteElement.h"

template <int nen, int ndm, int ndf>
void 
FiniteElement<nen,ndm,ndf>::setDomain(Domain *theDomain)
{
  if (theDomain == nullptr) {
    for (int i=0; i<nen; i++)
      theNodes[i] = nullptr;
    return;
  }

  for (int i=0; i<nen; i++) {
    theNodes[i] = theDomain->getNode(connectedExternalNodes(i));
    if (theNodes[i] == nullptr) {
      opserr << "FiniteElement::setDomain  tag: " 
              << this->getTag() << " -- Node " 
              << connectedExternalNodes(i) << " does not exist\n";
      return;
    }

    if (theNodes[i]->getNumberDOF() != ndf) {
      opserr << "FiniteElement::setDomain  tag: " << this->getTag() << " -- Node " << connectedExternalNodes(i) 
              << " has incorrect number of DOF\n";
      opserr << " " << theNodes[i]->getNumberDOF() << " should be " << ndf << endln;
      return;
    }
  }

  if (theDomain != nullptr)
    this->Element::link(*theDomain);

  if (this->setState(State::Init) != 0)
    return;

//    if (this->setState(State::Pres) != 0)
//      return;
}