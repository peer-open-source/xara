//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Description: This file implements the GenericResponse template which
// creates Response classes from objects implementing getResponse().
//
// Written: cmp
// Created: March 2024
//
#ifndef GenericResponse_h
#define GenericResponse_h

#include <Response.h>

class ID;
class Vector;
class Matrix;

template <typename T>
class GenericResponse : public Response
{
 public:
  GenericResponse(T &obj, int id)                    : Response(),    theObject(obj), responseID(id) {};
  GenericResponse(T &obj, int id, int val)           : Response(val), theObject(obj), responseID(id) {};
  GenericResponse(T &obj, int id, double val)        : Response(val), theObject(obj), responseID(id) {};
  GenericResponse(T &obj, int id, const ID &val)     : Response(val), theObject(obj), responseID(id) {};
  GenericResponse(T &obj, int id, const Vector &val) : Response(val), theObject(obj), responseID(id) {};
  GenericResponse(T &obj, int id, const Matrix &val) : Response(val), theObject(obj), responseID(id) {};

  GenericResponse(T *obj, int id)                    : Response(),    theObject(*obj), responseID(id) {};
  GenericResponse(T *obj, int id, int val)           : Response(val), theObject(*obj), responseID(id) {};
  GenericResponse(T *obj, int id, double val)        : Response(val), theObject(*obj), responseID(id) {};
  GenericResponse(T *obj, int id, const ID &val)     : Response(val), theObject(*obj), responseID(id) {};
  GenericResponse(T *obj, int id, const Vector &val) : Response(val), theObject(*obj), responseID(id) {};
  GenericResponse(T *obj, int id, const Matrix &val) : Response(val), theObject(*obj), responseID(id) {};

  ~GenericResponse() {};

  int
  getResponse() 
  {
    return theObject.getResponse(responseID, myInfo);
  }

protected:
  T & theObject;
  const int responseID;
};

#endif
