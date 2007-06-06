/* attrib_edge_down_len.h
   Author: Deepesh Kholwadwala, Sandia National Laboratories
   declares the attribute for storing edge length information.  
   Currently the length is stored only for downward pointing edges
   as required.
*/

#if !defined( ATTRIB_LENGTH_CLASS )
#define ATTRIB_PARAMETRIC_LENGTH_CLASS

class EDGE ;


#include "attrib_gtc.h"
#include "attrib_gtc_export.h"

extern DECL_GTCATTRIB int ATTRIB_PARAMETRIC_LENGTH_TYPE;
#define ATTRIB_PARAMETRIC_LENGTH_LEVEL (ATTRIB_GTC_LEVEL + 1)

class DECL_GTCATTRIB ATTRIB_PARAMETRIC_LENGTH: public ATTRIB_GTC {
   double length_data;
   double to_bottom_length_data;

public:

ATTRIB_PARAMETRIC_LENGTH( ENTITY* = NULL, double = 0.0);

double parametric_length() const { return length_data; }
void set_parametric_length(double);
double to_bottom_length() const {return to_bottom_length_data;}
void set_to_bottom_length(double);

// Functions called to aid attribute migration during modeling
// operations.

virtual void split_owner( ENTITY * );

// The owner of this attribute is about to be split in two - the
// argument is the new piece.  The owner of this attribute is about to
// be merged with the given entity.  The logical argument is TRUE if the
// owner is to be deleted in the merge.

   virtual void merge_owner( ENTITY *, // "other entity"
                             logical   // deleting_owner
                           );

   ATTRIB_FUNCTIONS( ATTRIB_PARAMETRIC_LENGTH, GTCATTRIB )
};

DECL_GTCATTRIB double get_parametric_length(EDGE*);
DECL_GTCATTRIB void set_parametric_length(EDGE*, double);

DECL_GTCATTRIB double get_to_bottom_length(EDGE*);
DECL_GTCATTRIB void set_to_bottom_length(EDGE*, double);
#endif

