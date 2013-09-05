#include "CubitCompat.h"
#include "CubitCompat.hpp"
#include "GeometryQueryTool.hpp"
#include <string.h>

#define CUBIT_COMPAT_FT_ELIF(TYPE) \
  else if (!strcmp(file_type,#TYPE)) \
    return TYPE ## _TYPE

static Model_File_Type
CubitCompat_file_type( const char* file_type )
{
  if (!file_type)
    return MFT_NOT_DEFINED;
  CUBIT_COMPAT_FT_ELIF(ACIS);
  CUBIT_COMPAT_FT_ELIF(ACIS_SAT);
  CUBIT_COMPAT_FT_ELIF(ACIS_SAB);
  CUBIT_COMPAT_FT_ELIF(ACIS_DEBUG);
  CUBIT_COMPAT_FT_ELIF(IGES);
  CUBIT_COMPAT_FT_ELIF(CATIA);
  CUBIT_COMPAT_FT_ELIF(STEP);
  CUBIT_COMPAT_FT_ELIF(PROE);
  CUBIT_COMPAT_FT_ELIF(GRANITE);
  CUBIT_COMPAT_FT_ELIF(GRANITE_G);
  CUBIT_COMPAT_FT_ELIF(GRANITE_SAT);
  CUBIT_COMPAT_FT_ELIF(GRANITE_PROE_PART);
  CUBIT_COMPAT_FT_ELIF(GRANITE_PROE_ASM);
  CUBIT_COMPAT_FT_ELIF(GRANITE_NEUTRAL);
  CUBIT_COMPAT_FT_ELIF(NCGM);
  CUBIT_COMPAT_FT_ELIF(CATIA_NCGM);
  CUBIT_COMPAT_FT_ELIF(CATPART);
  CUBIT_COMPAT_FT_ELIF(CATPRODUCT);
  CUBIT_COMPAT_FT_ELIF(FACET);
  CUBIT_COMPAT_FT_ELIF(SOLIDWORKS);
  CUBIT_COMPAT_FT_ELIF(OCC);
  else
    return MFT_NOT_DEFINED;
}

CubitStatus
CubitCompat_import_solid_model( const char* file_name,
                                const char* file_type,
                                const char* logfile_name,
                                CubitBoolean heal_step,
                                CubitBoolean import_bodies,
                                CubitBoolean import_surfaces,
                                CubitBoolean import_curves,
                                CubitBoolean import_vertices,
                                CubitBoolean free_surfaces,
				DLIList<RefEntity*> *imported_entities)
{
  const bool print_results = false;
  const bool merge_globally = false;
  const bool no_assembly_level_features = false;
  ModelImportOptions CubitCompat_opts = { heal_step,
                                        print_results,
                                        import_bodies,
                                        import_surfaces,
                                        import_curves,
                                        import_vertices,
                                        free_surfaces,
                                        merge_globally,
                                        no_assembly_level_features,
                                        logfile_name ? logfile_name : "" };

  return GeometryQueryTool::instance()->import_solid_model( 
           file_name,
           CubitCompat_file_type(file_type),
           CubitCompat_opts,
           imported_entities );
}

CubitStatus 
CubitCompat_export_solid_model( DLIList<RefEntity*>& ref_entity_list,
                                const char* filename,
                                const char * filetype,
                                int &num_ents_exported,
                                const CubitString &cubit_version,
                                const char* logfile_name )
{
  ModelExportOptions CubitCompat_opts = {1, logfile_name ? logfile_name : "" }; 

  return GeometryQueryTool::instance()->export_solid_model(
           ref_entity_list,
           filename,
           CubitCompat_file_type(filetype),
           num_ents_exported,
           cubit_version,
           CubitCompat_opts );
}
