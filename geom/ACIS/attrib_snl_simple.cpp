//- Class:       attrib_snl_simple
//- Owner:       Greg Nielson
//- Description: The implementation of the attrib_snl_simple ACIS
//-              attribute class.
//- Checked by:
//- Version: $Id:

#include <stdio.h>
#include <memory.h>
#include <string.h>

// ACIS Includes
#if CUBIT_ACIS_VERSION < 1100
#include "kernel/kernapi/api/api.hxx"
#include "kernel/kerndata/attrib/at_tsl.hxx"
#include "kernel/kerndata/data/datamsc.hxx"
#else
#include "api.hxx"	
#include "at_tsl.hxx"
#include "datamsc.hxx"
#include "acistype.hxx"
#include "edge.hxx"
#endif

#include "attrib_snl_simple.hpp"
#include "attrib_cubit_owner.hpp"
#include "CubitString.hpp"
#include "CubitSimpleAttrib.hpp"
#include "CubitAttribUser.hpp"
#include "CastTo.hpp"
#include "TopologyEntity.hpp"

// Buffer size
const int MAX_STRING_LENGTH = 256;

//Define macros for this attribute and its parent, to simplify later
// stuff and make it suitable for making into macros.

#define THIS() ATTRIB_SNL_SIMPLE
#define THIS_LIB NONE
#define PARENT() ATTRIB_SNL
#define PARENT_LIB NONE

// Identifier used externally to identify a particular entity type.
// This is used in the save/restore system for translating to/from
// external file format.

#define ATTRIB_SNL_SIMPLE_NAME "simple"

// Implement the standard attribute functions.

ATTCOPY_DEF( "snl_simple_attribute" )
  
  // Implementation of ATTRIB_SNL_SIMPLE::fixup_copy() 
  //  - For rollback ACIS does memcpy of attributes.  Fix things
  //    after memcpy so that a proper deep copy is done.
  if (data)
    data->inc_use_count();

 
  
LOSE_DEF

  // Implementation of ATTRIB_SNL_SIMPLE::lose()  



DTOR_DEF

  // Destructor implementation.
  if (data)
    data->dec_use_count();



DEBUG_DEF

  // Debug imeplementation
  


SAVE_DEF

  // Implement save to file
  
  write_string("NEW_SIMPLE_ATTRIB"); 
  
  if (!data)
  {
    write_int(0);
    write_int(0);
    write_int(0);
  }
  else
  {
    data->save();
  }



RESTORE_DEF
   
   assert(!data);

   char name[128];
   read_string(name);
   
   if(strcmp(name,"NEW_SIMPLE_ATTRIB")==0)
     data = AttribData::read_new_attrib();
   else
     data = AttribData::read_old_attrib(name);



COPY_DEF

  if (data)
    data->dec_use_count();
  
  data = from->data;

  if (data)
    data->inc_use_count();


SCAN_DEF
   // (no specific pointer data)



FIX_POINTER_DEF
   // (no specific pointer data)



TERMINATE_DEF




// Virtual function called when an owner entity is being split.

// ****THIS NEEDS TO BE CHANGED TO CALL THE CUBIT ATTRIBUTE
// ****DEFINED FUNCTIONS!!!!!!!!!!
void ATTRIB_SNL_SIMPLE::split_owner(ENTITY* new_ent)
{
  // If we have split an entity that is a hidden
  // entity of a composite transfer the attribute
  // to the new entity.
  if(!strcmp(data->name().c_str(), "COMPOSITE_GEOM") ||
     !strcmp(data->name().c_str(), "IMPRINT_PREEXISTING"))
  {
    API_BEGIN;
    new ATTRIB_SNL_SIMPLE(new_ent, *this);
    API_END;
    return;
  }

    // handle from attrib_cubit_owner if there is one
  ATTRIB_CUBIT_OWNER *owner_att =
      (ATTRIB_CUBIT_OWNER *)find_attrib(entity(),
                                        ATTRIB_SNL_TYPE, 
                                        ATTRIB_CUBIT_OWNER_TYPE);


  if (owner_att != NULL)
  {
    return;
  }

    // for now, make a new simple attrib for the new entity, but
    // only if this is a name attribute or an assembly attribute
  
  if(data &&
     (!strcmp(data->name().c_str(), "ENTITY_NAME") ||
      !strcmp(data->name().c_str(), "CA_ASSEMBLY_DATA")))
  {
    API_BEGIN;
    new ATTRIB_SNL_SIMPLE(new_ent, *this);
    API_END;
  }
}

void ATTRIB_SNL_SIMPLE::copy_owner( ENTITY *copy_ent )
{
  if(data)
  {
    if( !strcmp(data->name().c_str(), "ENTITY_ID") )
      return;
    
    API_BEGIN;
    new ATTRIB_SNL_SIMPLE(copy_ent, *this);
    API_END;
  }
}


// Virtual function called when two entities are to be merged.

// ****THIS NEEDS TO BE CHANGED TO CALL THE CUBIT ATTRIBUTE
// ****DEFINED FUNCTIONS!!!!!!!!!!
void ATTRIB_SNL_SIMPLE::merge_owner(ENTITY* other_ent, logical delete_owner)
{
  if(!strcmp(data->name().c_str(), "COMPOSITE_GEOM"))
  {
    if(is_EDGE(entity()) && is_EDGE(other_ent))
    {
      EDGE *this_edge = (EDGE*)entity();
      EDGE *other_edge = (EDGE*)other_ent;
      if((this_edge->start() == other_edge->start() &&
          this_edge->end() == other_edge->end()) ||
         (this_edge->start() == other_edge->end() &&
          this_edge->end() == other_edge->start()))
      {
        // If we get to this point we most likely have the
        // case where a geometric operation has cut this
        // body right along the hidden entity of a composite.
        // In this case we will just blow away the composite
        // by removing the attribute.
        unhook();
        lose();
      }
    }
    return;
  }

  if (delete_owner) return;

    // handle from attrib_cubit_owner if there is one
  ATTRIB_CUBIT_OWNER *owner_att =
      (ATTRIB_CUBIT_OWNER *)find_attrib(entity(),
                                        ATTRIB_SNL_TYPE, 
                                        ATTRIB_CUBIT_OWNER_TYPE);
  if (owner_att != NULL) return;
  
  if (data &&
      (!strcmp(data->name().c_str(),"ENTITY_NAME") ||
       !strcmp(data->name().c_str(), "CA_ASSEMBLY_DATA")))
  {
   // If the owner of this attribute is going to be deleted, then we
   // transfer ourself to that other entity.
    if (delete_owner) {
      move(other_ent);
    }
  }
}


void ATTRIB_SNL_SIMPLE::trans_owner(SPAtransf const&)
{
}


ATTRIB_SNL_SIMPLE::AttribData::AttribData( int num_strings, CubitString* strings,
                                           int num_ints, int* ints,
                                           int num_reals, double* reals ) 
  : useCount(1), 
    numStrings(num_strings), 
    numInts(num_ints), 
    numReals(num_reals),
    stringData(strings), 
    intData(ints), 
    realData(reals)
{ }


ATTRIB_SNL_SIMPLE::AttribData::AttribData(CubitSimpleAttrib* csa) 
  : useCount(1),
    numStrings(csa->string_data_list()->size()),
    numInts(csa->int_data_list()->size()),
    numReals(csa->double_data_list()->size()),
    stringData(0), intData(0), realData(0)
{
  int i;
  
  if (numStrings)

  {
    stringData = new CubitString[numStrings];
    csa->string_data_list()->reset();
    for (i = 0; i < numStrings; i++)
      stringData[i] = *csa->string_data_list()->next(i);
  }
  
  if (numReals)
  {
    realData = new double[numReals];
    csa->double_data_list()->reset();
    for (i = 0; i < numReals; i++)
      realData[i] = *csa->double_data_list()->next(i);
  }
  
  if (numInts)
  {
    intData = new int[numInts];
    csa->int_data_list()->reset();
    for (i = 0; i < numInts; i++)
      intData[i] = *csa->int_data_list()->next(i);
  }
}

ATTRIB_SNL_SIMPLE::AttribData::~AttribData()
{
  delete [] stringData;
  delete [] intData;
  delete [] realData;
}

CubitBoolean ATTRIB_SNL_SIMPLE::AttribData::equivalent(
                                 CubitSimpleAttrib *csa_ptr) const
{
  int i;
  
  if( csa_ptr->string_data_list()->size() != numStrings ||
      csa_ptr->int_data_list()->size() != numInts ||
      csa_ptr->double_data_list()->size() != numReals )
    return CUBIT_FALSE;
  
  csa_ptr->string_data_list()->reset();
  for( i = 0; i < numStrings; i++ )
    if( stringData[i] != *csa_ptr->string_data_list()->next(i) )
      return CUBIT_FALSE;
  
  csa_ptr->int_data_list()->reset();
  for( i = 0; i < numInts; i++ )
    if( intData[i] != *csa_ptr->int_data_list()->next(i) )
      return CUBIT_FALSE;
    
  csa_ptr->double_data_list()->reset();
  for( i = 0; i < numReals; i++ )
    if( realData[i] != *csa_ptr->double_data_list()->next(i) )
      return CUBIT_FALSE;

    // if we've gotten here, we're equivalent
  return CUBIT_TRUE;
}

CubitSimpleAttrib* ATTRIB_SNL_SIMPLE::AttribData::make_CSA() const
{
  int i;

  DLIList<CubitString*> string_list(numStrings);
  DLIList<int> int_list(numInts);
  DLIList<double> double_list(numReals);
  
  for( i = 0; i < numStrings; i++ )
    string_list.append( stringData + i );
  
  for( i = 0; i < numInts; i++ )
    int_list.append( *(intData + i) );
  
  for( i = 0; i < numReals; i++ )
    double_list.append( *(realData + i) );

  return new CubitSimpleAttrib( &string_list, 
                                numReals ? &double_list : 0, 
                                numInts ? &int_list : 0 );
}

void ATTRIB_SNL_SIMPLE::AttribData::dec_use_count()
{
  useCount--;
  if (useCount == 0)
    delete this;
}

void ATTRIB_SNL_SIMPLE::AttribData::save() const
{
  int i;
  write_int(numStrings);
  for (i = 0; i < numStrings; i++)
  {
    if( stringData[i].length() >= (unsigned)MAX_STRING_LENGTH )
    {
      PRINT_ERROR("Error saving ATTRIB_SNL_SIMPLE : "
                  "string length exceeds %d characters\n",
                  MAX_STRING_LENGTH - 1);
      write_string(stringData[i].substr(1,MAX_STRING_LENGTH-1).c_str());
    }
    else
    {
      write_string(stringData[i].c_str());
    }
  }

  write_int(numReals);
  for (i = 0; i < numReals; i++)
    write_real( realData[i] );

  write_int(numInts);
  for (i = 0; i < numInts; i++)
    write_int( intData[i] );
}


ATTRIB_SNL_SIMPLE::AttribData* ATTRIB_SNL_SIMPLE::AttribData::read_new_attrib()
{
   char buffer[MAX_STRING_LENGTH];
   int i, num_strings,  num_ints, num_reals;
   CubitString* strings = 0;
   int* ints = 0;
   double* reals = 0;

     // first do strings
   num_strings = read_int();
   if (num_strings)
   {
     strings = new CubitString[num_strings];
     for(i = 0; i < num_strings; i++)
     {
       read_string(buffer);
       strings[i] = buffer;
     }
   }

   num_reals = read_int();
   if (num_reals)
   {
     reals = new double[num_reals];
     for(i = 0; i < num_reals; i++)
       reals[i] = read_real();
   }

   num_ints = read_int();
   if (num_ints)
   {
     ints = new int[num_ints];
     for(i = 0; i < num_ints; i++)
       ints[i] = read_int();
   }
   
   return new AttribData(num_strings, strings, num_ints, ints, num_reals, reals);
}


ATTRIB_SNL_SIMPLE::AttribData* 
ATTRIB_SNL_SIMPLE::AttribData::read_old_attrib( const char* name )
{
     // read old-format attribute
     //  - {name, string1, string2, double, int}

   char buffer[MAX_STRING_LENGTH];
   CubitString* strings = new CubitString[3];
   int* ints = new int[1];
   double* reals = new double[1];

   strings[0] = name;
   read_string(buffer);
   strings[1] = buffer;
   read_string(buffer);
   strings[2] = buffer;
   reals[0] = read_real();
   ints[0] = read_int();

   return new AttribData(3, strings, 1, ints, 1, reals);
}
