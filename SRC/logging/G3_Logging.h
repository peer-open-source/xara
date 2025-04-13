#pragma once
#include "AnsiColors.h"
#include <OPS_Stream.h>
#define opserr (*opserrPtr)
#define opsdbg (*opsdbgPtr)
#define opswrn ((*opswrnPtr) << G3_WARN_PROMPT)
#define opsmrd opserr

class G3_Runtime;
extern OPS_Stream *opserrPtr;
extern OPS_Stream *opswrnPtr;
extern OPS_Stream *opsdbgPtr;
extern const char *G3_WARN_PROMPT;
extern const char *G3_ERROR_PROMPT;
extern const char *G3_DEBUG_PROMPT;

namespace OpenSees {
extern const char * PromptParseError;
extern const char * PromptValueError;
extern const char * PromptAnalysisFailure;
extern const char * PromptAnalysisSuccess;
extern const char * PromptAnalysisIterate;
};


enum G3_Stream {
  G3_StdOut, G3_StdIn, G3_StdErr, G3_Null
};

enum G3_StreamLevel {
  G3_LevelError, 
  G3_LevelDebug, 
  G3_LevelLog, 
  G3_LevelWarn
};

int G3_SetStreamColor(G3_Runtime* rt, int strm, int flag);
int G3_SetStreamLevel(int stream, bool on);
