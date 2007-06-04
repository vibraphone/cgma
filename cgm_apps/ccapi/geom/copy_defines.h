#ifndef COPY_DEFINES_HPP
#define COPY_DEFINES_HPP

#define COPY_LIST_TO_ARRAY(from, to, to_size) \
{  \
  if (from.size() && to == NULL)  \
    to = (void **) malloc(from.size()); \
  else if (from.size() && to && to_size < from.size()) { \
    free(to);  \
    to = (void **) malloc(from.size()); \
  } \
  if (from.size()) { \
    assert(to); \
    from.copy_to(to); \
  } \
  to_size = from.size(); \
}

#define COPY_ARRAY_TO_LIST(from, from_size, to) \
  if (from && from_size > 0) to.copy_from(from, from_size); 

#define COPY_LIST_TO_STRUCTARRAY(from, to, to_size, to_type) \
{  \
  if (from.size() && to == NULL)  \
    to = (to_type *) malloc(from.size()*sizeof(to[0])); \
  else if (from.size() && to && to_size < from.size()) { \
    free(to);  \
    to = (to_type *) malloc(from.size()*sizeof(to[0])); \
  } \
  if (from.size()) { \
    for (int STRUCT_SIZE = 0; STRUCT_SIZE < from.size(); STRUCT_SIZE++) \
      (to)[STRUCT_SIZE] = *from.get_and_step(); \
  } \
  to_size = from.size(); \
}

#define COPY_STRUCTARRAY_TO_LIST(from, from_size, to, to_type) \
{ \
  if (from && from_size > 0) { \
    for (int STRUCT_SIZE = 0; STRUCT_SIZE < from_size; STRUCT_SIZE++) { \
      to_type *NEW_PTR = new to_type(from[STRUCT_SIZE]); \
      (to).append(NEW_PTR); \
    } \
  } \
} 

#define DELETE_STRUCTLIST(list) \
  {for (int STRUCT_SIZE = (list).size(); STRUCT_SIZE > 0; STRUCT_SIZE--) \
      delete (list).get_and_step(); \
   (list).clean_out(); \
  }

#define COPY_ILIST_TO_ARRAY(from, to, to_size, to_type) \
{  \
  if (from.size() && to == NULL)  \
    to = (void *) malloc(from.size()*sizeof(to_type)); \
  else if (from.size() && to && to_size < from.size()) { \
    free(to);  \
    to = (void *) malloc(from.size()*sizeof(to_type)); \
  } \
  if (from.size()) { \
    assert(to); \
    from.copy_to(to); \
  } \
  to_size = from.size(); \
}

#define COPY_ARRAY_TO_ILIST(from, from_size, to) \
  if (from && from_size > 0) to.copy_from(from, from_size); 

#define COPY_ARRAY(from, from_size, to, to_size) \
  if (from_size > 0 && to == NULL) \
    to = (void **) malloc(from_size);\
  else if (from_size > 0 && to && to_size < from_size) { \
    free(to); \
    to = (void **) malloc(from_size); \
  } \
  memcpy(to, from, from_size*sizeof(void*));

#define GTI GeometryTool::instance()

#endif
