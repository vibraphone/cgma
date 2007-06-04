#ifndef VG_DEFINES_H
#define VG_DEFINES_H

enum {
  SUBCOMP_PARTITION_LAYER = 126,
  COMPOSITE_LAYER = 127,
  SUPERCOMP_PARTITION_LAYER = 128
};

  /* used for printing debug info - remove leading
     digits from type_info::name() on linux */
#ifdef __GNUC__
#  include <ctype.h>
  inline const char* fix_type_name(const char* name)
   { while( isdigit(*name) ) name++; return name; }
#else
#  define fix_type_name(A) A
#endif

#endif
