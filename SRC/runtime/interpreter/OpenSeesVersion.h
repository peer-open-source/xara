#include "Interpreter.h"

namespace Xara {

enum class OpenSeesVersion : int {
  O3 = 300,
  X1 = 3000,
  XaraLatest = 9999
};

OpenSeesVersion GetCompatibilityVersion(Tcl_Interp *interp);

} // namespace Xara