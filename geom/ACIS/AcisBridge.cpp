//-------------------------------------------------------------------------
// Filename      : AcisBridge.cpp
//
// Purpose       : Many functions are identical for each Acis-specific
//                 TopologyBridge.  This class implements those functions.
//
// Creator       : Darryl Melander
//
// Creation Date : 03/22/99
//
// Owner         : Darryl Melander
//-------------------------------------------------------------------------

// ************* BEGIN ACIS INCLUDES *************
#if CUBIT_ACIS_VERSION < 1100
#include "kernel/kernapi/api/api.hxx"
#else
#include "api.hxx"
#endif
// ************* END ACIS INCLUDES ***************

// ************* BEGIN CUBIT INCLUDES *************
#include "AcisBridge.hpp"
#include "TopologyBridge.hpp"
#include "GeometryQueryTool.hpp"
#include "AcisQueryEngine.hpp"
#include "CubitSimpleAttrib.hpp"
#include "attrib_cubit_owner.hpp"
#include "attrib_snl_simple.hpp"

#include "attrib_name.h"
#include "attrib_edge_parametric_length.h"
#include "CastTo.hpp"
// ************* END CUBIT INCLUDES ***************

unsigned AcisBridge::bridgeCount = 0;

AcisBridge::AcisBridge(ENTITY* entity) 
    : myENTITY(NULL)
{
  ENTITY_ptr(entity);
  bridgeCount++;
}

AcisBridge::~AcisBridge()
{
  ENTITY_ptr(NULL);
  bridgeCount--;
}

void AcisBridge::ENTITY_ptr(ENTITY* entity)
{
    // If we're not making a change, don't do anything
  if (entity == myENTITY)
    return;
  
    // If we are unhooking an entity...
  if (myENTITY)
  {
    ATTRIB_CUBIT_OWNER::remove_cubit_owner(myENTITY);
  }
  
    // Hook up new entity
  myENTITY = entity;
  if (myENTITY)
  {
    ATTRIB_CUBIT_OWNER::set_cubit_owner(myENTITY, this);
  }
}

void AcisBridge::append_simple_attribute_virt(CubitSimpleAttrib* attrib_ptr)
{
    // Append this name to the underlying ACIS BODY object.
  API_BEGIN;
  new ATTRIB_SNL_SIMPLE ( myENTITY, attrib_ptr );
  API_END;
}

void AcisBridge::remove_simple_attribute_virt(CubitSimpleAttrib* attrib_ptr)
{
  ATTRIB_SNL_SIMPLE *next =
    (ATTRIB_SNL_SIMPLE *) find_attrib(myENTITY,
                                      ATTRIB_SNL_TYPE,
                                      ATTRIB_SNL_SIMPLE_TYPE);
  
  attrib_ptr->string_data_list()->reset();
#ifdef BOYD16
  CubitString name = attrib_ptr->string_data_list()->size() ?
                     *attrib_ptr->string_data_list()->get() : 0;
#endif
  while( next != NULL )
  {
    ATTRIB_SNL_SIMPLE *attribute = next;
    next = (ATTRIB_SNL_SIMPLE *)
        find_next_attrib (attribute,
                          ATTRIB_SNL_TYPE,
                          ATTRIB_SNL_SIMPLE_TYPE);
    if( attribute->equivalent(attrib_ptr) )
    {
      API_BEGIN;
      attribute->unhook();
      API_END;
    }
  }

  CubitBoolean loop_continue = CUBIT_TRUE;
  attrib_ptr->string_data_list()->reset();
  if (attrib_ptr->string_data_list()->size() &&
      strcmp(attrib_ptr->string_data_list()->get()->c_str(),"ENTITY_NAME")==0)
  {
    loop_continue = CUBIT_TRUE;
    ATTRIB_GTC_NAME *attribute_gtc =
      (ATTRIB_GTC_NAME *) find_attrib(myENTITY,
                                      ATTRIB_GTC_TYPE,
                                      ATTRIB_GTC_NAME_TYPE);
    if (attribute_gtc != NULL)
    {
      for(;loop_continue && (attribute_gtc != NULL);
          attribute_gtc = (ATTRIB_GTC_NAME *)
            find_next_attrib (attribute_gtc,
                              ATTRIB_GTC_TYPE,
                              ATTRIB_GTC_NAME_TYPE))
      {
        attrib_ptr->string_data_list()->reset();
        // skip the type stored as the first string in the list
        if((strcmp(attribute_gtc->name(),
                   (attrib_ptr->string_data_list()->step_and_get())->c_str()) == 0) &&
           (strcmp(attribute_gtc->option(),
                   (attrib_ptr->string_data_list()->step_and_get())->c_str()) == 0))
        {
          attribute_gtc->unhook();
          loop_continue = CUBIT_FALSE;
        }
      }
    }

    ATTRIB_PARAMETRIC_LENGTH *attribute_len =
      (ATTRIB_PARAMETRIC_LENGTH *) find_attrib(myENTITY,
                                      ATTRIB_GTC_TYPE,
                                      ATTRIB_PARAMETRIC_LENGTH_TYPE);
    if (attribute_len != NULL)
    {
      for(;loop_continue && (attribute_len != NULL);
          attribute_len = (ATTRIB_PARAMETRIC_LENGTH *)
            find_next_attrib (attribute_len,
                              ATTRIB_GTC_TYPE,
                              ATTRIB_PARAMETRIC_LENGTH_TYPE))
      {
        attrib_ptr->string_data_list()->reset();

        attribute_len->unhook();
        loop_continue = CUBIT_FALSE;

      }
    }
  }
}

void AcisBridge::remove_all_simple_attribute_virt()
{
  ATTRIB_SNL_SIMPLE *attribute =
    (ATTRIB_SNL_SIMPLE *) find_attrib(myENTITY,
                                      ATTRIB_SNL_TYPE,
                                      ATTRIB_SNL_SIMPLE_TYPE);
  for(;attribute != NULL;
      attribute = (ATTRIB_SNL_SIMPLE *)
        find_attrib (myENTITY,
                     ATTRIB_SNL_TYPE,
                     ATTRIB_SNL_SIMPLE_TYPE))
  {
    API_BEGIN;
    attribute->unhook();
    API_END;
  }
  ATTRIB_GTC_NAME *attribute_gtc =
    (ATTRIB_GTC_NAME *) find_attrib(myENTITY,
                                    ATTRIB_GTC_TYPE,
                                    ATTRIB_GTC_NAME_TYPE);
  for(;attribute_gtc != NULL;
      attribute_gtc = (ATTRIB_GTC_NAME *)
        find_attrib (myENTITY,
                     ATTRIB_GTC_TYPE,
                     ATTRIB_GTC_NAME_TYPE))
  {
    API_BEGIN;
    attribute_gtc->unhook();
    API_END;
  }

  ATTRIB_PARAMETRIC_LENGTH *attribute_len =
    (ATTRIB_PARAMETRIC_LENGTH *) find_attrib(myENTITY,
                                    ATTRIB_GTC_TYPE,
                                    ATTRIB_PARAMETRIC_LENGTH_TYPE);
  for(;attribute_len != NULL;
      attribute_len = (ATTRIB_PARAMETRIC_LENGTH *)
        find_attrib (myENTITY,
                     ATTRIB_GTC_TYPE,
                     ATTRIB_PARAMETRIC_LENGTH_TYPE))
  {
    API_BEGIN;
    attribute_len->unhook();
    API_END;
  }
}

CubitStatus AcisBridge::get_simple_attribute( ENTITY *entity,
                                              DLIList<CubitSimpleAttrib*>& attrib_list)
{
  ATTRIB_SNL_SIMPLE *attribute =
    (ATTRIB_SNL_SIMPLE *) find_attrib(entity, 
                                      ATTRIB_SNL_TYPE,
                                      ATTRIB_SNL_SIMPLE_TYPE);
  for(;attribute != NULL;attribute = (ATTRIB_SNL_SIMPLE *)
        find_next_attrib(attribute,
                         ATTRIB_SNL_TYPE,
                         ATTRIB_SNL_SIMPLE_TYPE))
  {
    attrib_list.append( attribute->get_CSA() );
  }
  
    // This is to take care of lingering GTC name attributes.
  ATTRIB_GTC_NAME *attribute_gtc =
    (ATTRIB_GTC_NAME *) find_attrib(entity,
                                    ATTRIB_GTC_TYPE,
                                    ATTRIB_GTC_NAME_TYPE);
  for(;attribute_gtc != NULL;attribute_gtc = (ATTRIB_GTC_NAME *)
        find_next_attrib(attribute_gtc,
                         ATTRIB_GTC_TYPE,
                         ATTRIB_GTC_NAME_TYPE))
  {
     CubitString name_and_option(attribute_gtc->name());
     if (attribute_gtc->option()[0] != '\0') {
       name_and_option += ".";
       name_and_option += attribute_gtc->option();
     }
    CubitSimpleAttrib* cubit_simple_attrib_ptr;
    cubit_simple_attrib_ptr = new CubitSimpleAttrib("ENTITY_NAME",
						     name_and_option.c_str());


    attrib_list.append_unique(cubit_simple_attrib_ptr);


    //attrib_list.reset();
    //CubitBoolean unique = CUBIT_TRUE;
    //for(int j = 0; j < attrib_list.size() && unique; j++){
    //  CubitSimpleAttrib* cs_attrib_ptr =
    //    attrib_list.get_and_step();
    //  if (CubitSimpleAttrib::equivalent(cs_attrib_ptr, cubit_simple_attrib_ptr)){
    //    unique = CUBIT_FALSE;
    //  }
    //}
    //if (unique)
    //{
    //  AcisBridge *bridge = ATTRIB_CUBIT_OWNER::cubit_owner(entity);
    //  TopologyBridge *tb = CAST_TO(bridge, TopologyBridge);
    //  assert(tb != 0);
    //  tb->append_simple_attribute(cubit_simple_attrib_ptr);
    //  attrib_list.append_unique(cubit_simple_attrib_ptr);
    //}
    //else
    //{
    //  delete cubit_simple_attrib_ptr;
    //}
  }

    // This is to take care of lingering GTC parametric length attributes.
  ATTRIB_PARAMETRIC_LENGTH *attribute_parametric_len =
    (ATTRIB_PARAMETRIC_LENGTH *) find_attrib(entity,
                                    ATTRIB_GTC_TYPE,
                                    ATTRIB_PARAMETRIC_LENGTH_TYPE);
  for(;attribute_parametric_len != NULL;attribute_parametric_len = (ATTRIB_PARAMETRIC_LENGTH *)
        find_next_attrib(attribute_parametric_len,
                         ATTRIB_GTC_TYPE,
                         ATTRIB_PARAMETRIC_LENGTH_TYPE))
  {
    CubitSimpleAttrib* cubit_simple_attrib_ptr;
    cubit_simple_attrib_ptr = new CubitSimpleAttrib("MESH_RELATIVE_LENGTH", 
													 "","",0,
                                                     attribute_parametric_len->parametric_length());

    attrib_list.append_unique(cubit_simple_attrib_ptr);
  }

  return CUBIT_SUCCESS;
}

CubitStatus AcisBridge::get_simple_attribute( ENTITY *entity, const CubitString& name,
                                              DLIList<CubitSimpleAttrib*>& attrib_list)
{
  ATTRIB_SNL_SIMPLE *attribute =
    (ATTRIB_SNL_SIMPLE *) find_attrib(entity, 
                                      ATTRIB_SNL_TYPE,
                                      ATTRIB_SNL_SIMPLE_TYPE);
  for(;attribute != NULL;attribute = (ATTRIB_SNL_SIMPLE *)
        find_next_attrib(attribute,
                         ATTRIB_SNL_TYPE,
                         ATTRIB_SNL_SIMPLE_TYPE))
  {
    if( attribute->attribute_name() == name )
      attrib_list.append( attribute->get_CSA() );
  }
  return CUBIT_SUCCESS;
}

AcisQueryEngine *AcisBridge::get_acis_query_engine() const
{
  AcisQueryEngine *aqe = AcisQueryEngine::instance();

    // if this function gets called, it needs to return non-NULL
  assert(aqe != 0);

  return aqe;
}

