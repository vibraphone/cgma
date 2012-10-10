//-------------------------------------------------------------------------
// Filename      : AppUtil.hpp 
//
// Purpose       : This file represents the Cubit application itself.
//
// Special Notes : 
//
// Creator       : Darryl Melander
//
// Date          : 06/08/98
//
// Owner         : Darryl Melander
//-------------------------------------------------------------------------

#ifndef APP_UTIL_HPP
#define APP_UTIL_HPP

#include "CubitDefines.h"
#include "CubitString.hpp"
#include "CubitUtilConfigure.h"

#ifndef WIN32
#include <sys/resource.h>
#endif

class ProgressTool;
class CubitString;

class CUBIT_UTIL_EXPORT AppUtil
{
public:
   static AppUtil* instance();
     //- Access to the application object

   static void delete_instance();
        
   void report_resource_usage() const;
     //- Prints out information about the session
     
   bool get_terminal_size( int& rows, int& cols );
   //- Get the size of the terminal, if known.

#ifdef CAT
   int days_to_expire(int year, int month, int day) const;
     //- Returns number of days until the app expires
#endif
   ~AppUtil();
   

  // Signal handing code.  The signal handler provided by
  // AppUtil sets the flag cubit_intr to CUBIT_TRUE when
  // an interrupt (SIGINT) is detected.  See the comments
  // with the declaration of cubit_intr in AppUtil.cpp
  // for more information.

  // returns whether the interupt flag has been set
  CubitBoolean interrupt();
  void set_interrupt(CubitBoolean);

  // clears the interrupt flag
  void clear_interrupt();
  
  static CubitBoolean catch_interrupt() { return catching_sigint_; }
  //- Check if the signal handler provided in AppUtil
  //- is being used.
  
  static void catch_interrupt( CubitBoolean yesno );
  //- yesno:
  //-  CUBIT_TRUE:  Set signal handler for SIGINT to the
  //-               handler provided by AppUtil.
  //-  CUBIT_FALSE: Set the signal hander for SIGINT to
  //-               the system default.
  
  void startup(int /*argc*/, char ** /*argv*/);
   //-  Contains startup code for cubit

  int shutdown();
   //- Contains shutdown code for cubit

  void apputil_getrusage(struct rusage &r_usage) const;
    //- fill the r_usage with data from getrusage; implements special code to
    //- find memory usage when you don't get that from the system call getrusage
   
   void set_terminal_size( int rows, int cols );
   //- This function sets the values returned from
   //- get_terminal_size().  It has no effect on
   //- the actual terminal size.

  ProgressTool *progress_tool();
  void progress_tool(ProgressTool* pTool);
   // pass in a pointer to a progress tool  - the lifetime of the tool will be
   // controlled by AppUtil

  CubitString get_cubit_dir();
    // get/set the cubit dir variable

  static void initialize_settings();
    // initialize settings

private:

  static CubitBoolean catching_sigint_;
   //- Signal handler registered for SIGINT?
   
   int term_width, term_height;
   
   static AppUtil* instance_;
   CubitBoolean mAppStarted;

  CubitString cubitDir;
    //- directory of the cubit executable

   ProgressTool *mProgressTool;

   AppUtil();
};

// returns the directory where cubit binaries reside, 
// or more correctly, the path of the library/executable that contains this AppUtil code.
inline CubitString AppUtil::get_cubit_dir() 
{
  return cubitDir;
}

#endif

