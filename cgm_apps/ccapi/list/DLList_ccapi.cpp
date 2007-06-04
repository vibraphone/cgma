#include "DLList_ccapi.h"

#include "DLList.hpp"
#include "CubitDefines.h"

#include "copy_defines.h"

  void DLList_compress(void ***this_list, int *this_list_size)
{
  int i;
  int num_nulls = 0;
  
  for (i = 0; i < *this_list_size; i++) {
    if ((*this_list)[i] == NULL) num_nulls++;
    else (*this_list)[i-num_nulls] = (*this_list)[i];
  }

  *this_list_size -= num_nulls;    
}	

  
  void DLList_reverse(void ***this_list, int *this_list_size)
{
  int i;
  void *temp;
  
  for (i = 0; i < (*this_list_size)/2; i++) {
    temp = (*this_list)[i];
    (*this_list)[i] = (*this_list)[*this_list_size-i-1];
    (*this_list)[*this_list_size-i-1] = temp;
  }
}	

  
    /* DLList& operator= */ void DLList_equals(void ***this_list, int *this_list_size,
                                                 /* const DLList& */ void ***from_list,
                                               int *from_list_size)
{
  COPY_ARRAY(*this_list, *this_list_size, *from_list, *from_list_size);
}	

  
  void DLList_merge_unique(void ***this_list, int *this_list_size,
                             /* DLList& */ void ***merge_list,
                           int *merge_list_size,
                           int merge_list_unique)
{
  int new_size = *this_list_size;
  int *loop_end = (merge_list_unique ? this_list_size :
                   &new_size);
  
  int i, j;
  int found;
  for (i = 0; i < *merge_list_size; i++) {
    found = 0;
    for (j = 0; j < *loop_end; j++) {
      if ((*this_list)[j] == (*merge_list)[i]) {
        found = 1;
        break;
      }
    }
    if (!found) (*this_list)[new_size++] = (*merge_list)[i];
  }

  *this_list_size = new_size;
}

  void DLList_intersect(void ***this_list, int *this_list_size,
                          /* DLList& */ void ***merge_list,
                        int *merge_list_size)
{
  int i, j;
  int found;
  int did_something = 0;
  
  for (i = 0; i < *merge_list_size; i++) {
    found = 0;
    for (j = 0; j < *this_list_size; j++) {
      if ((*this_list)[j] == (*merge_list)[i]) {
        found = 1;
        break;
      }
    }

    if (!found) {
      (*this_list)[j] = NULL;
      did_something = 1;
    }
  }

  if (did_something) DLList_compress(this_list, this_list_size);
}	

  
  int DLList_append_unique(void ***this_list, int *this_list_size,
                           void* new_item)
{
  int i;
  int found = 0;
  
  for (i = 0; i < *this_list_size; i++)
    if ((*this_list)[i] == new_item) {
      found = 1;
      break;
    }

  if (!found) {
    free(*this_list);
    (*this_list_size)++;
    *this_list = (void **) malloc(*this_list_size*sizeof(void *));
    (*this_list)[*this_list_size-1] = new_item;
    return CUBIT_TRUE;
  }
  else
    return CUBIT_FALSE;
}	
