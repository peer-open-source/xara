#pragma once

#  include <handler/OPS_Stream.h>
#ifndef opserr
   extern OPS_Stream *opserrPtr;
#  define opserr (*opserrPtr)
#  define endln "\n"
#endif

extern OPS_Stream *opslogPtr;
#  define opslog (*opslogPtr)
#include "G3_Logging.h"

#include <logging/AnsiColors.h>
#define LOG_TEST ":: "
#define LOG_ITERATE BLU "   ITERATE" COLOR_RESET " :: "
#define LOG_FAILURE RED "   FAILURE" COLOR_RESET " :: "
#define LOG_SUCCESS GRN "   SUCCESS" COLOR_RESET " :: "
#define LOG_CONTINUE "\n              "
