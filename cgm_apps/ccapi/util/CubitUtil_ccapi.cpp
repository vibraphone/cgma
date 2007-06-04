#include "CubitUtil_ccapi.h"

#include "CubitUtil.hpp"
#include "CubitString.hpp"
#include "DLIList.hpp"
#include "DLCubitEntityList.hpp"

#include "copy_defines.h"

void CubitUtil_convert_string_to_lowercase(char *string)
{
  CubitUtil::convert_string_to_lowercase(string);
}

int  CubitUtil_strcmp_case_insensitive(const char *s1, const char *s2)
{
  return CubitUtil::strcmp_case_insensitive(s1, s2);
}

int  CubitUtil_strncmp_case_insensitive(const char *s1, const char *s2,
                                               int n)
{
  return CubitUtil::strncmp_case_insensitive(s1, s2, n);
}

void CubitUtil_list_entity_ids_1( const char *pre_string, 
                                         /* const DLCubitEntityList & */ void ***entity_list, int *entity_list_size,
                                       int width, const char *post_string,
                                       int sort, int unique,
                                       int tab, const char *sep_string,
                                       const char *post_string_none)
{
  DLCubitEntityList temp_entity_list;
  COPY_ARRAY_TO_LIST(*entity_list, *entity_list_size, temp_entity_list);
  
  CubitUtil::list_entity_ids(pre_string, temp_entity_list,
                             width, post_string, sort, unique, tab,
                             sep_string, post_string_none);
}
  
void CubitUtil_list_entity_ids_2( const char *pre_string, 
                                         /* DLIList<int> & */ int **id_list, int *id_list_size,
                                       int width, 
                                       const char *post_string,
                                       int sort, int unique,
                                       int tab_len, const char *sep_string,
                                       const char *post_string_none)
{
  DLIList<int> temp_id_list;
  COPY_ARRAY_TO_ILIST(*id_list, *id_list_size, temp_id_list);
  
  CubitUtil::list_entity_ids(pre_string, temp_id_list,
                             width, post_string, sort, unique, tab_len,
                             sep_string, post_string_none);
}
  
/* CubitString */ char *CubitUtil_get_entity_ids_str( const char *pre_string,
                                                               /* DLIList<int> & */ int **int_list, int *int_list_size,
                                                             int width, const char *post_string,
                                                             int sort, int unique,
                                                             int left_tab, const char *sep_string,
                                                             const char *post_string_none)
{
  DLIList<int> temp_int_list;
  COPY_ARRAY_TO_ILIST(*int_list, *int_list_size, temp_int_list);

  CubitString temp_string(CubitUtil::get_entity_ids_str(pre_string, temp_int_list,
                                                        width, post_string, sort, unique, left_tab,
                                                        sep_string, post_string_none));

  char *temp_char_string = new char[temp_string.length()];
  strcpy(temp_char_string, temp_string.c_str());
  
  return temp_char_string;
}

void CubitUtil_process_entity_ids( int method,
                                            /* CubitString & */ char **ret_str,
                                          const char *pre_string, 
                                            /* DLIList<int> & */ int **id_list, int *id_list_size,
                                          int max_len, 
                                          const char *post_string,
                                          int sort, int unique,
                                          int tab_len, const char *sep_string,
                                          const char* post_string_none)
{
  DLIList<int> temp_id_list;
  COPY_ARRAY_TO_ILIST(*id_list, *id_list_size, temp_id_list);

  CubitString temp_ret_str;

  CubitUtil::process_entity_ids(method, temp_ret_str, pre_string, temp_id_list,
                                max_len, post_string, sort, unique,
                                tab_len, sep_string,
                                post_string_none);

  *ret_str = new char[temp_ret_str.length()];
  strcpy(*ret_str, temp_ret_str.c_str());
}
  
int CubitUtil_int_len( int num )
{
  return CubitUtil::int_len(num);
}

void CubitUtil_set_file_ptr( FILE* file_ptr )
{
  CubitUtil::set_file_ptr(file_ptr);
}

void CubitUtil_reset_file_ptr()
{
  CubitUtil::reset_file_ptr();
}

FILE* CubitUtil_get_file_ptr()
{
  return CubitUtil::get_file_ptr();
}

int CubitUtil_compare( const char* a, const char* b )
{
  return CubitUtil::compare(a, b);
}

enum CubitSense CubitUtil_opposite_sense(enum CubitSense sense)
{
  return CubitUtil::opposite_sense(sense);
}
