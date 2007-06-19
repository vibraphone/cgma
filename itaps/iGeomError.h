#ifndef I_GEOM_ERROR_H
#define I_GEOM_ERROR_H

#ifdef __cplusplus
extern "C" {
#endif

void iGeom_clearLastError();

void iGeom_setLastError( int error_type, const char* description = 0 );

void iGeom_getLastError( int& error_type_out, 
                         char* description_buffer,
                         int description_buffer_length );

int iGeom_getLastErrorType();

const char* iGeom_getLastErrorDesc();

#ifdef __cplusplus
 } // extern "C"
#endif

#endif
