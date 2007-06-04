#ifndef CUBIT_APP_CCAPI
#define CUBIT_APP_CCAPI

#ifdef __cplusplus
extern "C" {
#endif
  
/* CubitApp* */ 
void *CubitApp_instance();
    /* Access to the application object */

void CubitApp_delete_instance();
        
void CubitApp_report_resource_usage();
    /* Prints out information about the session */

int CubitApp_days_to_expire(int year, int month, int day);
    /* Returns number of days until the app expires */

  void gl_cleanup();
    /* This function is currently called from CubitApp, even though it
       is defined in getline.cpp, which is outside cgm.
    */

#ifdef __cplusplus
}
#endif

#endif

