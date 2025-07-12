#include <cstdarg>
#include <tcl.h>
#include <string.h>
#include <elementAPI.h>
class G3_Runtime;


#include <StandardStream.h>
#include <FileStream.h>
#include <DummyStream.h>
StandardStream sserr;
DummyStream    ssnul;
OPS_Stream *opserrPtr = &sserr;
OPS_Stream *opsdbgPtr = &ssnul;
OPS_Stream *opslogPtr = &ssnul;
OPS_Stream *opswrnPtr = &sserr;
OPS_Stream *opsmrdPtr = &sserr;


#include "G3_Logging.h"


namespace OpenSees {

namespace Internal {
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

} // namespace OpenSees::Internal

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
}

const char * G3_WARN_PROMPT  = OpenSees::Internal::WarnPromptNoColor;
const char * G3_ERROR_PROMPT = OpenSees::Internal::ErrorPromptNoColor;
const char * G3_DEBUG_PROMPT = OpenSees::Internal::DebugPromptNoColor;

int
G3_SetStreamLevel(int stream, bool on)
{
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

int G3_SetStreamColor(G3_Runtime* rt, int strm, int flag)
{
  if (flag == 1) {
    G3_WARN_PROMPT                  = OpenSees::Internal::WarnPromptColor;
    OpenSees::SignalWarning         = OpenSees::Internal::WarnPromptColor;
    G3_DEBUG_PROMPT                 = OpenSees::Internal::DebugPromptColor;
    OpenSees::PromptParseError      = OpenSees::Internal::ErrorPromptColor;
    OpenSees::PromptValueError      = OpenSees::Internal::ErrorPromptColor;
    OpenSees::PromptModelError      = OpenSees::Internal::ErrorPromptColor;
    OpenSees::PromptAnalysisFailure = OpenSees::Internal::AnalysisFailureColor;
    OpenSees::PromptAnalysisSuccess = OpenSees::Internal::AnalysisSuccessColor;
    OpenSees::PromptAnalysisIterate = OpenSees::Internal::AnalysisIterateColor;

  } else if (flag == 0) {
    G3_WARN_PROMPT             = OpenSees::Internal::WarnPromptNoColor;
    OpenSees::SignalWarning    = OpenSees::Internal::WarnPromptNoColor;
    G3_DEBUG_PROMPT            = OpenSees::Internal::DebugPromptNoColor;
    OpenSees::PromptParseError = OpenSees::Internal::ErrorPromptNoColor;
    OpenSees::PromptAnalysisFailure = OpenSees::Internal::AnalysisFailureNoColor;
    OpenSees::PromptAnalysisSuccess = OpenSees::Internal::AnalysisSuccessNoColor;
    OpenSees::PromptAnalysisIterate = OpenSees::Internal::AnalysisIterateNoColor;
  }

  return 0;
}

