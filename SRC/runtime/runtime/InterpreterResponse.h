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
// Description: The InterpreterResponse class is used to send a response to the
// interpreter. It is a thin wrapper around the abstract Response class
//
// Created: 06/2026
#pragma once

class Response;

class InterpreterResponse {
public:
  enum class Type {
    Vector
  };
  InterpreterResponse(Response* response, Type type, int size, double* data)
    : response(response), type(type), size(size), data(data) {}

  ~InterpreterResponse();

public:
  const int size;
  const Type type;
  const double* data;
  Response* const response;
};
