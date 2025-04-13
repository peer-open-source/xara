#pragma once

#ifndef opserr
#  include <handler/OPS_Stream.h>
   extern OPS_Stream *opserrPtr;
#  define opserr (*opserrPtr)
#  define endln "\n"
#endif

#include "G3_Logging.h"

#include <logging/AnsiColors.h>
#define LOG_TEST ":: "
#define LOG_ITERATE BLU "   ITERATE" COLOR_RESET " :: "
#define LOG_FAILURE RED "   FAILURE" COLOR_RESET " :: "
#define LOG_SUCCESS GRN "   SUCCESS" COLOR_RESET " :: "
#define LOG_CONTINUE "\n              "
