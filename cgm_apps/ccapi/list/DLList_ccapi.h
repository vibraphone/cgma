#ifndef DLLIST_CCAPI_H
#define DLLIST_CCAPI_H

#ifdef __cplusplus
extern "C" {
#endif
  
void DLList_compress(void ***this_list, int *this_list_size);
    /* - Remove all null (0) pointers from list. Must be called after */
    /* - {nullify_link()} is called before calling any list modification  */
    /* - function except {move_to()} or {nullify_link()} */
  
void DLList_reverse(void ***this_list, int *this_list_size);
    /* - Reverse the items in the list. */
  
  /* DLList& operator= */ void DLList_equals(void ***this_list, int *this_list_size,
                                             void ***to_list, int *to_list_size,
                                               /* const DLList& */ void ***from_list,
                                             int *from_list_size);
    /* - Create a copy of a list. */
  
void DLList_merge_unique(void ***this_list, int *this_list_size,
                           /* DLList& */ void ***merge_list,
                         int *merge_list_size,
                         int merge_list_unique);
    /* - Merges the contents of the list, merge_list, with those of "this" */
    /* - list, ensuring that items that are being added do not already */
    /* - exist in "this" list. If the items are known to appear in the */
    /* - merge_list only once, then it can be done faster. */

void DLList_intersect(void ***this_list, int *this_list_size,
                        /* DLList& */ void ***merge_list,
                      int *merge_list_size);
    /* - Merges the contents of the list merge_list, with those of "this" */
    /* - list, but only the items that are common to both lists are kept */
    /* - in "this" list.  This is a rather costly process I think. */
    /* - if merge_list is the same as the calling list, it just returns */
  
int DLList_append_unique(void ***this_list, int *this_list_size,
                         void* new_item);
    /* - Appends the new item to the list, if it doesn't already exist */
    /* - in the list. In either case, index is not changed. */
    /* - Return CUBIT_TRUE if the item was added, else CUBIT_FALSE. */

#ifdef __cplusplus
}
#endif
  
#endif
