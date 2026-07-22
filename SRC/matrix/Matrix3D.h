//===----------------------------------------------------------------------===//
//
//                                   xara
//
//===----------------------------------------------------------------------===//
//                              https://xara.so
//===----------------------------------------------------------------------===//
//
// Claudio Perez
//
#ifndef Matrix3D_H
#define Matrix3D_H

#include "MatrixND.h"
#include <type_traits>

namespace OpenSees {

  using Matrix3D = MatrixND<3,3,double>;

  // Layout
  static_assert(std::is_standard_layout<Matrix3D>::value, "MatrixND is not standard layout.");
  static_assert(std::is_aggregate<Matrix3D>::value, "MatrixND is not an aggregate type.");

  // Triviality
  static_assert(std::is_trivially_copyable<Matrix3D>::value, "MatrixND is not trivially copyable.");
  static_assert(std::is_trivially_move_assignable<Matrix3D>::value, "MatrixND is not trivially move assignable.");
  static_assert(std::is_trivial<Matrix3D>::value, "MatrixND is not trivial.");

  // Nothrow
  static_assert(std::is_nothrow_default_constructible<Matrix3D>::value, "MatrixND is not nothrow default constructible.");
  static_assert(std::is_nothrow_constructible<Matrix3D>::value, "MatrixND is not nothrow constructible.");
  static_assert(std::is_nothrow_move_assignable<Matrix3D>::value, "MatrixND is not nothrow move assignable.");
}

#endif // Matrix3D_H
