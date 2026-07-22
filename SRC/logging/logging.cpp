#include <cstdarg>
#include <tcl.h>
#include <string.h>
class G3_Runtime;


#include <StandardStream.h>
#include <FileStream.h>
#include <DummyStream.h>


#include "Logging.h"


namespace Xara::Logging::Internal {

StandardStream sserr;
DummyStream    ssnul;
OPS_Stream *opserrPtr = &sserr;
OPS_Stream *opsdbgPtr = &ssnul;
OPS_Stream *opslogPtr = &ssnul;
OPS_Stream *opswrnPtr = &sserr;

// namespace Internal {
  const char * WarnPromptColor   = RED "WARNING " COLOR_RESET;
  const char * WarnPromptNoColor = "WARNING ";

  const char * ErrorPromptColor   = BRED "ERROR " COLOR_RESET;
  const char * ErrorPromptNoColor = "ERROR ";

  const char * DebugPromptColor   = GRN "DEBUG " COLOR_RESET;
  const char * DebugPromptNoColor = "DEBUG ";

  const char * AnalysisIterateColor    = BLU "   ITERATE" COLOR_RESET " :: ";
  const char * AnalysisIterateNoColor  =     "   ITERATE"             " :: ";

  const char * AnalysisFailureColor    = RED "   FAILURE" COLOR_RESET " :: ";
  const char * AnalysisFailureNoColor  =     "   FAILURE"             " :: ";

  const char * AnalysisSuccessColor    = GRN "   SUCCESS" COLOR_RESET " :: ";
  const char * AnalysisSuccessNoColor  =     "   SUCCESS"             " :: ";

} // namespace Xara::Logging::Internal 

namespace OpenSees {
  using namespace Xara::Logging;
  // Default to no color
  const char * SignalMessageEnd      = "\n";
  const char * PromptParseError      = Internal::ErrorPromptNoColor;
  const char * PromptValueError      = PromptParseError;
  const char * PromptModelError      = PromptParseError;
  const char * SignalWarning         = Internal::WarnPromptNoColor;

  const char * PromptDomainFailure   = Internal::AnalysisFailureNoColor;
  const char * PromptAnalysisFailure = Internal::AnalysisFailureNoColor;
  const char * PromptAnalysisSuccess = Internal::AnalysisSuccessNoColor;
  const char * PromptAnalysisIterate = Internal::AnalysisIterateNoColor;

} // namespace OpenSees

const char * G3_WARN_PROMPT  = OpenSees::Internal::WarnPromptNoColor;
const char * G3_ERROR_PROMPT = OpenSees::Internal::ErrorPromptNoColor;
const char * G3_DEBUG_PROMPT = OpenSees::Internal::DebugPromptNoColor;

int
G3_SetStreamLevel(int stream, bool on)
{
  using namespace Xara::Logging::Internal;
  OPS_Stream **theStream;
  switch (stream) {
    case G3_LevelError: theStream = &opserrPtr; break;
    case G3_LevelDebug: theStream = &opsdbgPtr; break;
    case G3_LevelWarn : theStream = &opswrnPtr; break;
    default:
      return -1;
  }

  if (on) {
    *theStream = &sserr;
  } else {
    *theStream = &ssnul;
  }
  return 0;
}

int
G3_SetStreamColor(G3_Runtime* rt, int strm, int flag)
{
  using namespace Xara::Logging::Internal;
  if (flag == 1) {
    G3_WARN_PROMPT                  = WarnPromptColor;
    OpenSees::SignalWarning         = WarnPromptColor;
    G3_DEBUG_PROMPT                 = DebugPromptColor;
    OpenSees::PromptParseError      = ErrorPromptColor;
    OpenSees::PromptValueError      = ErrorPromptColor;
    OpenSees::PromptModelError      = ErrorPromptColor;
    OpenSees::PromptAnalysisFailure = AnalysisFailureColor;
    OpenSees::PromptAnalysisSuccess = AnalysisSuccessColor;
    OpenSees::PromptAnalysisIterate = AnalysisIterateColor;

  } else if (flag == 0) {
    G3_WARN_PROMPT             = WarnPromptNoColor;
    OpenSees::SignalWarning    = WarnPromptNoColor;
    G3_DEBUG_PROMPT            = DebugPromptNoColor;
    OpenSees::PromptParseError = ErrorPromptNoColor;
    OpenSees::PromptAnalysisFailure = AnalysisFailureNoColor;
    OpenSees::PromptAnalysisSuccess = AnalysisSuccessNoColor;
    OpenSees::PromptAnalysisIterate = AnalysisIterateNoColor;
  }

  return 0;
}
