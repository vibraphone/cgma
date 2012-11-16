#ifndef CUBIT_COMPAT_HPP
#define CUBIT_COMPAT_HPP

#include "CubitDefines.h"

class CubitString;
class RefEntity;
template <typename T> class DLIList;

CubitStatus
CubitCompat_import_solid_model( const char* file_name,
                                const char* file_type,
                                const char* logfile_name = NULL,
                                CubitBoolean heal_step = CUBIT_TRUE,
                                CubitBoolean import_bodies = CUBIT_TRUE,
                                CubitBoolean import_surfaces = CUBIT_TRUE,
                                CubitBoolean import_curves = CUBIT_TRUE,
                                CubitBoolean import_vertices = CUBIT_TRUE,
                                CubitBoolean free_surfaces = CUBIT_TRUE,
				DLIList<RefEntity*> *imported_entities = NULL);

CubitStatus 
CubitCompat_export_solid_model( DLIList<RefEntity*>& ref_entity_list,
                                const char* filename,
                                const char * filetype,
                                int &num_ents_exported,
                                const CubitString &cubit_version,
                                const char* logfile_name = NULL );

#endif /* CUBIT_COMPAT_HPP */

