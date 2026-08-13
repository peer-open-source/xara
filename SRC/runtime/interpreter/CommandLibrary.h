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
#include <Parsing.h>

#include <string>
#include <algorithm>
#include <unordered_map>

namespace Xara::Internal {

std::string toLower(const std::string & s);

bool equalsIgnoreCase(const std::string & lhs, const std::string & rhs );

class CaseInsensitive
{
public:
  inline size_t operator( ) (const std::string & s ) const
  {  
    static std::hash<std::string> hf;
    return hf( toLower( s ) );
  }
  
  inline bool operator( ) (const std::string & lhs, const std::string & rhs ) const
  {
    return equalsIgnoreCase( lhs, rhs );
  }
};

} // namespace Xara::Internal

using CommandLibrary = std::unordered_map<std::string, Tcl_CmdProc *, 
                                          Xara::Internal::CaseInsensitive, 
                                          Xara::Internal::CaseInsensitive> ;

template <typename T>
using CaseInsensitiveMap = std::unordered_map<std::string, T, 
                                          Xara::Internal::CaseInsensitive, 
                                          Xara::Internal::CaseInsensitive> ;
