#include <tcl.h>
#include <stdlib.h>
#include <Parsing.h>

OpenSeesVersion 
GetCompatibilityVersion(Tcl_Interp *interp)
{
  OpenSeesVersion version = OpenSeesVersion::XaraLatest;
  const char *version_str = nullptr;
  if (getenv("XARA_COMPATIBILITY_VERSION") != nullptr) {
    version_str = getenv("XARA_COMPATIBILITY_VERSION");
    if (strcmp(version_str, "O3") == 0)
      version = OpenSeesVersion::O3;
    else {
      opserr << "WARNING: unknown compatibility version '" << version_str << "'.\n";
    }
    return version;
  }

  
  // 
  return version;
}
