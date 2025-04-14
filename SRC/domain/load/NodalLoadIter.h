//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Written: cmp, Spring 2025
//

// File: ~/domain/loadcase/NodalLoadIter.h
//
// Written: fmk 
// Created: Fri Sep 20 15:27:47: 1996
// Revision: A
//

#ifndef NodalLoadIter_h
#define NodalLoadIter_h

#include <TaggedIterator.hpp>

class NodalLoad;

class NodalLoadIter: public TaggedIterator<NodalLoad> {
  public:
  using TaggedIterator<NodalLoad>::TaggedIterator;
};

#endif

