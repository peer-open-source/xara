//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Description: This file implements the SensitiveResponse template which
// creates Response classes from objects implementing getResponse().
//
// Written: CMP 
// Created: March 2024
//
#ifndef SensitiveResponse_h
#define SensitiveResponse_h
#include "GenericResponse.h"

template <typename T>
class SensitiveResponse : public GenericResponse<T> {

public:
  using GenericResponse<T>::GenericResponse;

  int getResponseSensitivity(int grad) {
    return this->theObject.getResponseSensitivity(this->responseID, grad, this->myInfo);
  }

};
#endif
