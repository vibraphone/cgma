#include "CubitApp_ccapi.h"
#include "CubitApp.hpp"

/* CubitApp* */ 
void *CubitApp_instance() 
{
  return CubitApp::instance();
}

  //- Access to the application object

void CubitApp_delete_instance()
{
  CubitApp::delete_instance();
}
        
void CubitApp_report_resource_usage()
{
  CubitApp::instance()->report_resource_usage();
}
  //- Prints out information about the session

int CubitApp_days_to_expire(int year, int month, int day)
{
  return CubitApp::instance()->days_to_expire(year, month, day);
}
  //- Returns number of days until the app expires


void gl_cleanup() {}
