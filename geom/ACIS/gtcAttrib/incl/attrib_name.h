//  attrib_gtc_name.hpp
//   Author: Arlo L. Ames, Sandia National Laboratories
//   declares the attribute_gtc_name structure; built by the book according to
//   acis documentation

#if !defined( ATTRIB_GTC_NAME_CLASS )
#define ATTRIB_GTC_NAME_CLASS

#include "attrib_gtc.h"
#include "attrib_gtc_export.h"

extern DECL_GTCATTRIB int ATTRIB_GTC_NAME_TYPE;
#define ATTRIB_GTC_NAME_LEVEL (ATTRIB_GTC_LEVEL + 1)

class DECL_GTCATTRIB ATTRIB_GTC_NAME: public ATTRIB_GTC {

private:
   char* name_data;
   char* option_data;

public:

  ATTRIB_GTC_NAME( ENTITY* = NULL, const char* = NULL, const char* = NULL);

  char* name() const { return name_data; }
  char* option() const { return option_data; }
  void set_name(const char*);
  void set_option(const char*);

// Functions called to aid attribute migration during modeling
// (boolean) operations.

  virtual void split_owner( ENTITY * );
//- The owner of this attribute is about to be split in two - the
//- argument is the new piece.  

  virtual void merge_owner( ENTITY *, // "other entity"
                             logical  // deleting_owner
                          );
//- The owner of this attribute is about to be merged with the 
//- given entity.  The logical argument is TRUE if the owner is to 
//- be deleted in the merge.


   ATTRIB_FUNCTIONS( ATTRIB_GTC_NAME, GTCATTRIB )
};

DECL_GTCATTRIB void set_name(ENTITY*, const char* const);
DECL_GTCATTRIB void set_option(ENTITY*, const char* const);

DECL_GTCATTRIB char* get_name(ENTITY*);
DECL_GTCATTRIB char* get_option(ENTITY*);


#endif

