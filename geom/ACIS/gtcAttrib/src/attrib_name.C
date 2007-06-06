// attrib_gtc_name
// Author: Arlo L. Ames, Sandia National Laboratories
// Implementation of the ATTRIB_GTC_NAME attribute, a child of attrib_gtc

#include "attrib_name.h"


#include <stdio.h>
#include <string.h>
#include <memory.h>

// ACIS Includes
#if ACIS_VER < 1100
#include "kernel/kernapi/api/api.hxx"
#include "kernel/kerndata/attrib/at_tsl.hxx"
#include "kernel/kerndata/data/datamsc.hxx"
#else
#include "api.hxx"
#include "at_tsl.hxx"
#include "datamsc.hxx"
#endif



//Define macros for this attribute and its parent, to simplify later
// stuff and make it suitable for making into macros.

#define THIS() ATTRIB_GTC_NAME
#define THIS_LIB GTCATTRIB
#define PARENT() ATTRIB_GTC
#define PARENT_LIB GTCATTRIB


// Identifier used externally to identify a particular entity type.
// This is used in the save/restore system for translating to/from
// external file format.

#define ATTRIB_GTC_NAME_NAME "name"

// Implement the standard attribute functions.

ATTCOPY_DEF( "name_attribute" )
   rollback->name_data = NULL ;
   rollback->option_data = NULL ;

LOSE_DEF

DTOR_DEF
   if (name_data != NULL)
      delete [] name_data ;
   if (option_data != NULL)
      delete [] option_data ;
   name_data = NULL ;
   option_data = NULL ;

DEBUG_DEF

SAVE_DEF
   if (NULL != name_data)
      write_string(name_data); //Save specific data
   else 
      write_string("") ;

   if (NULL != option_data && 
       strcmp(option_data, "gto_backpointer_claro_interface"))
      write_string(option_data);
   else
      write_string("") ;

RESTORE_DEF
   {
      char buffer[80];
      read_string(buffer);
      set_name(buffer);
      read_string(buffer);
      set_option(buffer);
   }

COPY_DEF
   if (from != this)
   {
      set_name(from->name());
      set_option(from->option());
   }

SCAN_DEF
   // (no specific pointer data)

FIX_POINTER_DEF
   // (no specific pointer data)

TERMINATE_DEF

// make a name attribute

ATTRIB_GTC_NAME::ATTRIB_GTC_NAME(
   ENTITY* owner,
   const char* const nam,
   const char* const opt
   ) : ATTRIB_GTC(owner), name_data(NULL), option_data(NULL)
{
   //Initialize members
   set_name(nam);
   set_option(opt);
}

// Set the member data.

void
ATTRIB_GTC_NAME::set_name(const char* const nam
)
{
   backup();
   if (!nam)
   {
      if (!name_data)
         name_data = new char[2] ;
      strcpy(name_data, "");
   }
   else
   {
      if (name_data) 
         delete [] name_data;
      name_data = new char[strlen(nam)+1];
      if (name_data)
         strcpy(name_data, nam);
   }
}

void
ATTRIB_GTC_NAME::set_option(const char* const opt
)
{
   backup();
   if (option_data) 
      delete [] option_data;
   if (!opt)
   {
      option_data = NULL;
      return;
   }
   else
   {
      option_data = new char[strlen(opt)+1];
      if (option_data) 
         strcpy(option_data, opt);
   }
}

// Virtual function called when an owner entity is being split.

void
ATTRIB_GTC_NAME::split_owner(ENTITY* new_ent)
{
   // Duplicate the attribute on the new entity
  API_BEGIN;
  new ATTRIB_GTC_NAME(new_ent, name_data, option_data);
  API_END;
}

// Virtual function called when two entities are to be merged.

void
ATTRIB_GTC_NAME::merge_owner(ENTITY* other_ent, logical delete_owner)
{
   // If the owner of this attribute is goint to be deleted, and there is
   // no attribute attached to the other entity, the we transfer ourself
   // to that other entity.

   if (delete_owner) {
      ATTRIB* other_att = find_attrib(
         other_ent,
         ATTRIB_GTC_TYPE,
         ATTRIB_GTC_NAME_TYPE
      );
      if (other_att == NULL ) {
         // No name on other entity, so transfer ourself
         move(other_ent);
      }
   }
}


void set_name(ENTITY* entity, const char* const name)
{

   // Search for an existing name attribute and set it.

   ATTRIB_GTC_NAME *name_att = (ATTRIB_GTC_NAME *) find_attrib(entity, 
							       ATTRIB_GTC_TYPE, 
							       ATTRIB_GTC_NAME_TYPE);
   while (name_att != NULL && 
          !strcmp(name_att->option(), "gto_backpointer_claro_interface"))
      name_att = (ATTRIB_GTC_NAME *) find_next_attrib(name_att, 
                                                      ATTRIB_GTC_TYPE, 
                                                      ATTRIB_GTC_NAME_TYPE);
   if (name_att != NULL)
   {
      //cout << "already had a name\n" << endl;
      //name_att->backup(); already done in next line
      name_att->set_name(name);
   }
   else
   {
      //cout << "adding a name\n" << endl;
      // Create a new attribute if none found
     API_BEGIN;
     new ATTRIB_GTC_NAME(entity, (char *)name, "");
     API_END;
   }
}

void set_option(ENTITY* entity, const char* const opt)
{
   // Search for an existing name attribute and set it.

   ATTRIB_GTC_NAME *name_att = (ATTRIB_GTC_NAME *) find_attrib(entity, 
							       ATTRIB_GTC_TYPE, 
							       ATTRIB_GTC_NAME_TYPE);

   while (name_att != NULL && 
          !strcmp(name_att->option(), "gto_backpointer_claro_interface"))
      name_att = (ATTRIB_GTC_NAME *) find_next_attrib(name_att, 
                                                      ATTRIB_GTC_TYPE, 
                                                      ATTRIB_GTC_NAME_TYPE);

   if (name_att != NULL)
   {
      //cout <<"already had a name\n" << endl;
      //name_att->backup(); already done in next line
      name_att->set_option(opt);
   }
   else
   {
      //cout <<"adding a name\n" << endl;
      // Create a new attribute if none found
     API_BEGIN;
     new ATTRIB_GTC_NAME(entity, "", (char *)opt);
     API_END;
   }
}
   

char* get_name(ENTITY* entity)
{
  char *get_name_result = NULL;
  ATTRIB_GTC_NAME *name_att = 
     (ATTRIB_GTC_NAME *) find_attrib(entity,
                                     ATTRIB_GTC_TYPE,
                                     ATTRIB_GTC_NAME_TYPE);

   while (name_att != NULL && 
          !strcmp(name_att->option(), "gto_backpointer_claro_interface"))
      name_att = (ATTRIB_GTC_NAME *) find_next_attrib(name_att, 
                                                      ATTRIB_GTC_TYPE, 
                                                      ATTRIB_GTC_NAME_TYPE);
  if (NULL != name_att) 
     get_name_result = name_att->name();

  return get_name_result;
}


char* get_option(ENTITY* entity)
{
  char *get_option_result = NULL;
  ATTRIB_GTC_NAME *option_att = 
     (ATTRIB_GTC_NAME *) find_attrib(entity, 
                                     ATTRIB_GTC_TYPE,
                                     ATTRIB_GTC_NAME_TYPE);

  if (option_att) 
     get_option_result = option_att->option();

  return get_option_result;
}

