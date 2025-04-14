//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// This class supercedes the original MaterialResponse class, so that
// it is only a thin wrapper around the GenericResponse template.
// Eventually this class should be removed altogether, and only
// GenericResponse should be used. The prior implementation of 
// MaterialResponse was the only thing that depended on having a generic
// Material base class, which was otherwise useless. Migrating to 
// GenericResponse eliminates entirely the need for Material.
//
// Written: cmp
// Created: Oct 2024
//
#ifndef MaterialResponse_h
#define MaterialResponse_h

#include <Response.h>
#include <Information.h>
#include <GenericResponse.h>

class Material;

class ID;
class Vector;
class Matrix;

template <typename T>
class MaterialResponse : public GenericResponse<T> {
  public:
  using GenericResponse<T>::GenericResponse;
};
template <typename T, typename ...B>
MaterialResponse(T*, B...)->MaterialResponse<T>;

#endif // include guard

