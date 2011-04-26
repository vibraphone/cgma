#ifndef I_GEOM_ERROR_H
#define I_GEOM_ERROR_H

#ifdef __cplusplus
extern "C" {
#endif

void iGeom_clearLastError();

void iGeom_setLastError( int error_type, const char* description = 0 );

int iGeom_getLastErrorType();

void iGeom_getLastErrorDesc( char* description_buffer,
                             int description_buffer_length );

#ifdef __cplusplus
 } // extern "C"
#endif

#endif
