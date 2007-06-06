/* attrib_edge_down_len.C
   Author: Deepesh Kholwadwala, Sandia National Laboratories
   declares the attribute for storing edge length information.  
   Currently the length is stored only for downward pointing edges
   as required.
*/
// Implementation of the length attribute


#include "attrib_edge_parametric_length.h"

#include "alltop.hxx"
#include "datamsc.hxx"
#include "api.hxx"

//Define macros for this attribute and its parent, to simplify later
// stuff and make it suitable for making into macros.

#define THIS() ATTRIB_PARAMETRIC_LENGTH
#define PARENT() ATTRIB_GTC
//#define PARENT_LIB NONE
//#define THIS_LIB NONE
#define PARENT_LIB GTCATTRIB
#define THIS_LIB GTCATTRIB


// Identifier used externally to identify a particular entity type.
// This is used in the save/restore system for translating to/from
// external file format.

#define ATTRIB_PARAMETRIC_LENGTH_NAME "parametric_length"

// Implement the standard attribute functions.

ATTCOPY_DEF( "parametric_length_attribute" )

LOSE_DEF

DTOR_DEF

DEBUG_DEF

SAVE_DEF
   write_real(length_data);
RESTORE_DEF
   length_data = read_real();
COPY_DEF
   length_data = from->length_data ;
   to_bottom_length_data = from->to_bottom_length_data ;

SCAN_DEF
   // (no specific pointer data)

FIX_POINTER_DEF
   // (no specific pointer data)

TERMINATE_DEF

// make a length attribute

ATTRIB_PARAMETRIC_LENGTH::ATTRIB_PARAMETRIC_LENGTH(
   ENTITY* owner,
   double length_
) : ATTRIB_GTC(owner)
{
   // Initialize members
   set_parametric_length(length_);
}

// Set the member data.

void
ATTRIB_PARAMETRIC_LENGTH::set_parametric_length(
   double length_
)
{
   backup();
   length_data = length_;
}

// Virtual function called when an owner entity is being split.

void
ATTRIB_PARAMETRIC_LENGTH::set_to_bottom_length( double length_)
{
   backup();
   to_bottom_length_data = length_;
}

void
ATTRIB_PARAMETRIC_LENGTH::split_owner(
   ENTITY* new_ent
)
{
}
// Virtual function called when two entities are to be merged.

void
ATTRIB_PARAMETRIC_LENGTH::merge_owner(
   ENTITY* other_ent,
   logical delete_owner
)
{
}

double get_parametric_length(EDGE* edge)
{
  double get_parametric_length_result = 0.0;
  ATTRIB_PARAMETRIC_LENGTH *parametric_length_att = 
    (ATTRIB_PARAMETRIC_LENGTH *) find_attrib(edge, ATTRIB_GTC_TYPE, ATTRIB_PARAMETRIC_LENGTH_TYPE);
  if (parametric_length_att)  get_parametric_length_result =  parametric_length_att->parametric_length();
  return get_parametric_length_result;
}


void set_parametric_length(EDGE* edge, double down_length)
{
API_BEGIN
   ATTRIB_PARAMETRIC_LENGTH *parametric_length_att = (ATTRIB_PARAMETRIC_LENGTH *) find_attrib(edge, ATTRIB_GTC_TYPE, ATTRIB_PARAMETRIC_LENGTH_TYPE);
   if (!parametric_length_att)
   {
     parametric_length_att = new ATTRIB_PARAMETRIC_LENGTH(edge, down_length);
   }
   else
      parametric_length_att->set_parametric_length(down_length);

API_END
 
}

double get_to_bottom_length(EDGE* edge)
{
  double get_to_bottom_length_result = 0.0;
    ATTRIB_PARAMETRIC_LENGTH *parametric_length_att = 
    (ATTRIB_PARAMETRIC_LENGTH *) find_attrib(edge, ATTRIB_GTC_TYPE, ATTRIB_PARAMETRIC_LENGTH_TYPE);
  if (parametric_length_att)
    get_to_bottom_length_result = parametric_length_att->to_bottom_length();

  return get_to_bottom_length_result;
}

void set_to_bottom_length(EDGE* edge, double down_length)
{
API_BEGIN
   ATTRIB_PARAMETRIC_LENGTH *parametric_length_att = (ATTRIB_PARAMETRIC_LENGTH *) find_attrib(edge, ATTRIB_GTC_TYPE, ATTRIB_PARAMETRIC_LENGTH_TYPE);
   if (!parametric_length_att)
   {
     parametric_length_att = new ATTRIB_PARAMETRIC_LENGTH(edge, 0.0);
   }
   parametric_length_att->set_to_bottom_length(down_length);

API_END
 
}

