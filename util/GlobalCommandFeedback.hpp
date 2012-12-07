#ifndef GLOBAL_COMMAND_FEEDBACK_HPP
#define GLOBAL_COMMAND_FEEDBACK_HPP

#include "CommandFeedback.hpp"
#include "ManagedPtrVector.hpp"
#include "CubitUtilConfigure.h"

namespace GlobalCommandFeedback
{
  CUBIT_UTIL_EXPORT void create();
  CUBIT_UTIL_EXPORT void destroy();
  CUBIT_UTIL_EXPORT ManagedPtrVector<CommandFeedback>& instance();
}

#endif
