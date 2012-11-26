//-------------------------------------------------------------------------
// Filename      : CubitProcess.hpp
//
// Purpose       : This file contains utility functions to call other programs
//
// Creator       : Clinton Stimpson
//
// Date          : 07/29/10
//-------------------------------------------------------------------------

#include "CubitUtilConfigure.h"

#include <string>
#include <vector>

#ifndef CubitProcess_hpp
#define CubitProcess_hpp

#ifdef WIN32
#include <windows.h>
typedef PROCESS_INFORMATION PidType;
#else
typedef pid_t PidType;
#endif

class CUBIT_UTIL_EXPORT CubitProcess
{
public:
  
  // execute a program with some arguments
  // returns the pid the child program, or 0 for the process id if there was an error starting the program
  // for windows, a .exe is automatically added to the app name
  // this function returns and doesn't wait for the child process to exit
  static PidType start(const std::string& app, const std::vector<std::string>& args, bool hide=false);
  
  // wait for a child process to exit
  // returns the exit code of the program, or -1 if there was an error starting the program
  // 0 is returned if the programed exited normally and without errors
  // 1 - 127 is returned if the programed exited normally and with errors
  // 128 and greater is returned if the programed exited abnormally (such as segfault, abort, etc..)
  static int wait(PidType pid);

  // execute a program with some arguments and returns the exit code of the program
  // basically calls  start()/wait()
  // returns -1 if there was an error starting the program
  // 0 is returned if the programed exited normally and without errors
  // 1 - 127 is returned if the programed exited normally and with errors
  // 128 and greater is returned if the programed exited abnormally (such as segfault, abort, etc..)
  static int execute(const std::string& app, const std::vector<std::string>& args, bool hide=false);


  // given an executable file, return the absolute path to it
  // the input can already be absolute or relative
  static std::string find_executable(const std::string& exe);
};

#endif

