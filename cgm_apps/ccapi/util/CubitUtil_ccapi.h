#ifndef CUBIT_UTIL_CCAPI_H
#define CUBIT_UTIL_CCAPI_H

#include <string.h>
#include "CubitDefines.h"

#ifdef __cplusplus
extern "C" {
#endif
  
  void CubitUtil_convert_string_to_lowercase(char *string);

  int  CubitUtil_strcmp_case_insensitive(const char *s1, const char *s2);

  int  CubitUtil_strncmp_case_insensitive(const char *s1, const char *s2,
                                       int n);
  void CubitUtil_list_entity_ids_1( const char *pre_string, 
                                 /* const DLCubitEntityList & */ void ***entity_list, int *entity_list_size,
                               int width, const char *post_string,
                               int sort, int unique,
                               int tab, const char *sep_string,
                               const char *post_string_none);
  
  void CubitUtil_list_entity_ids_2( const char *pre_string, 
                                 /* DLIList<int> & */ int **id_list, int *id_list_size,
                               int width, 
                               const char *post_string,
                               int sort, int unique,
                               int tab_len, const char *sep_string,
                               const char *post_string_none);
  
  /* CubitString */ char *CubitUtil_get_entity_ids_str( const char *pre_string,
                                                       /* DLIList<int> & */ int **int_list, int *int_list_size,
                                                     int width, const char *post_string,
                                                     int sort, int unique,
                                                     int left_tab, const char *sep_string,
                                                     const char *post_string_none);

  void CubitUtil_process_entity_ids( int method,
                                    /* CubitString & */ char **ret_str,
                                  const char *pre_string, 
                                    /* DLIList<int> & */ int **id_list, int *id_list_size,
                                  int max_len, 
                                  const char *post_string,
                                  int sort, int unique,
                                  int tab_len, const char *sep_string,
                                  const char* post_string_none);
  
  int CubitUtil_int_len( int num );
    /* - Finds the number of spaces required to print an integer number */

  void CubitUtil_set_file_ptr( FILE* file_ptr );

  void CubitUtil_reset_file_ptr();

  FILE* CubitUtil_get_file_ptr();
    /* - Used to optionally have list_entity_ids or get_entity_ids_str output dumped */
    /* - to a file as well. */

  int CubitUtil_compare( const char* a, const char* b );

  enum CubitSense CubitUtil_opposite_sense(enum CubitSense sense);
    /* - return the sense opposite from sense */

#ifdef __cplusplus
}
#endif

#endif
