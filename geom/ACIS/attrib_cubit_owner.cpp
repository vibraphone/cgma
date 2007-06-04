//-------------------------------------------------------------------------
// Filename      : attrib_cubit_owner.hpp
//
// Purpose       : This attribute represents a pointer from an
//                 ACIS ENTITY to a Cubit AcisBridge
//
// Creator       : Greg Neilson
//
// Owner         : Tim Tautges
//-------------------------------------------------------------------------

// ********** BEGIN STANDARD INCLUDES         **********
#include <stdio.h>
#include <memory.h>
// ********** END STANDARD INCLUDES           **********

// ********** BEGIN ACIS INCLUDES             **********
#if CUBIT_ACIS_VERSION < 1100
#include "kernel/kernapi/api/api.hxx"
#include "kernel/kerndata/data/datamsc.hxx"
#include "kernel/kerndata/geom/transfrm.hxx"
#include "kernel/kerndata/bulletin/bulletin.hxx"
#else
#include "api.hxx"	
#include "datamsc.hxx"
#include "transfrm.hxx"
#include "bulletin.hxx"
#endif
// ********** END ACIS INCLUDES               **********

// ********** BEGIN CUBIT INCLUDES            **********
#include "attrib_cubit_owner.hpp"
#include "AcisBridge.hpp"
#include "GeometryQueryTool.hpp"
#include "GeometryQueryEngine.hpp"
#include "AcisQueryEngine.hpp"
#include "TopologyEntity.hpp"
#include "TopologyBridge.hpp"
#include "CubitMessage.hpp"
#include "CubitAttribUser.hpp"
#include "CastTo.hpp"
#include "attrib_snl_simple.hpp"
#include "DLIList.hpp"
#include "RefEntity.hpp"
// ********** END CUBIT INCLUDES              **********


// ********** BEGIN MACRO DEFINITIONS         **********
#define THIS() ATTRIB_CUBIT_OWNER
#define THIS_LIB NONE
#define PARENT() ATTRIB_SNL
#define PARENT_LIB NONE

#define ATTRIB_CUBIT_OWNER_NAME "cubit_owner"

ATTRIB_DEF( "cubit_owner_attribute")

SAVE_DEF
// Don't save

RESTORE_DEF
// Don't restore

COPY_DEF
   set_cubit_owner ( from->cubit_owner() );
//    PRINT_INFO("COPY_DEF called on attrib_cubit_owner.\n");

SCAN_DEF
// (no specific pointer data)

FIX_POINTER_DEF
// (no specific pointer data)

TERMINATE_DEF
// Don't do anything special
// ********** END MACRO DEFINITIONS           **********

bool ATTRIB_CUBIT_OWNER::copyingAttribs = true;

// make a cubit_owner attribute
ATTRIB_CUBIT_OWNER::ATTRIB_CUBIT_OWNER(ENTITY* owner,
                                       AcisBridge *cubit_owner)
    : ATTRIB_SNL(owner), cubitOwnerData(NULL)
{
  set_cubit_owner(cubit_owner);
}

logical ATTRIB_CUBIT_OWNER::savable() const
  { return FALSE; }

AcisBridge* ATTRIB_CUBIT_OWNER::cubit_owner(ENTITY *ENTITY_ptr)
{
  ATTRIB_CUBIT_OWNER *attribute = NULL;
    // Get the attribute if there is one
   attribute = 
      (ATTRIB_CUBIT_OWNER *) find_attrib(ENTITY_ptr, 
                                         ATTRIB_SNL_TYPE,
                                         ATTRIB_CUBIT_OWNER_TYPE);
  
    // If there is an attribute and an owner, return it, else return NULL.
  if (attribute == NULL)
    return NULL;
  else
    return attribute->cubit_owner();
}

void ATTRIB_CUBIT_OWNER::cubit_owner(ENTITY_LIST &entity_list,
                                     DLIList<AcisBridge*> &tb_list)
{
  for (int i = 0; i < entity_list.count(); i++)
  {
    ENTITY *entity = entity_list[i];
    AcisBridge *tb = cubit_owner(entity);
    if (tb != NULL) tb_list.append(tb);
  }
}

void ATTRIB_CUBIT_OWNER::cubit_owner(ENTITY_LIST &entity_list,
                                     DLIList<TopologyBridge*> &tb_list)
{
  for (int i = 0; i < entity_list.count(); i++)
  {
    ENTITY *entity = entity_list[i];
    AcisBridge *tb = cubit_owner(entity);
    if (tb != NULL)
      tb_list.append(CAST_TO(tb, TopologyBridge));
  }
}

TopologyEntity* ATTRIB_CUBIT_OWNER::get_topology_entity(ENTITY* ENTITY_ptr)
{
    // Get the attribute if there is one
  ATTRIB_CUBIT_OWNER *attribute =
    (ATTRIB_CUBIT_OWNER *) find_attrib(ENTITY_ptr, 
                                       ATTRIB_SNL_TYPE,
                                       ATTRIB_CUBIT_OWNER_TYPE);
  
    // If there is an attribute and an owner, return it, else return NULL.
  if (attribute == NULL || attribute->cubit_owner() == NULL)
    return NULL;
  else
  {
    TopologyBridge* tb = CAST_TO(attribute->cubit_owner(), TopologyBridge);
    return tb->topology_entity();
  }
}

// set the member data.
void ATTRIB_CUBIT_OWNER::set_cubit_owner (AcisBridge* new_cubit_owner)
{
  if (new_cubit_owner != NULL)
    backup();
  cubitOwnerData = new_cubit_owner;
}

void ATTRIB_CUBIT_OWNER::set_cubit_owner(ENTITY* acis_entity,
                                         AcisBridge *cubit_entity)
{
    // Some checks
  if (cubit_entity == NULL)
  {
    PRINT_ERROR("Trying to set the CUBIT_OWNER attrribute of "
                "an ACIS ENTITY to be NULL in set_cubit_owner\n");
    return;
  }
  
  if (acis_entity == NULL)
  {
    PRINT_ERROR("Trying to set the CUBIT_OWNER attribute of a "
                " NULL ACIS ENTITY in set_cubit_owner\n");
    return;
  }
  
    // Search for the attribute
  ATTRIB_CUBIT_OWNER *owner_att =
    (ATTRIB_CUBIT_OWNER *)find_attrib(acis_entity,
                                      ATTRIB_SNL_TYPE, 
                                      ATTRIB_CUBIT_OWNER_TYPE);
  
    // If found, reset it's value
  if (owner_att != NULL)
  {
    owner_att->set_cubit_owner(cubit_entity);
  }
  
    // Otherwise, create a new attribute for acis_entity
  else
  {
    API_BEGIN;
    ATTRIB_CUBIT_OWNER *owner = new ATTRIB_CUBIT_OWNER(acis_entity, cubit_entity);
    owner->history()->set_max_states_to_keep( 1 ); 
    API_END;
  }
}

void
ATTRIB_CUBIT_OWNER::remove_cubit_owner(ENTITY* ENTITY_ptr,
                                       CubitBoolean children_too)
{
  ENTITY_LIST entities;
  entities.add(ENTITY_ptr);
  if (children_too == CUBIT_TRUE) {
      // need to do something special to get the aqe, since this is a
      // static function
    AcisQueryEngine::instance()
      ->get_child_ENTITYs(ENTITY_ptr, entities, CUBIT_TRUE);
  }

  for (int i = 0; i < entities.count(); i++)
  {
    ENTITY *this_ENTITY = entities[i];
    assert(this_ENTITY != NULL);

      // Find an attribute of type ATTRIB_CUBIT_OWNER_TYPE
    ATTRIB_CUBIT_OWNER* ATTRIB_ptr = 
      (ATTRIB_CUBIT_OWNER *) find_attrib(ENTITY_ptr, 
                                         ATTRIB_SNL_TYPE,
                                         ATTRIB_CUBIT_OWNER_TYPE);
    
      // If the attribute is found, remove it. In the following code
      // block, a do-while structure is used to ensure that all such
      // attributes are removed. At present, at most 1 CUBIT_OWNER 
      // attribute is allowed, per ACIS ENTITY.
    while (ATTRIB_ptr != NULL)
    {
          // "unhook" this attribute from the ENTITY (i.e., remove it from
          // the list of attributes associated with this ENTITY
        API_BEGIN;
        ATTRIB_ptr->unhook();
        API_END;
        
          // Find the next CUBIT_OWNER attribute attached to this ENTITY
        ATTRIB_ptr =
          (ATTRIB_CUBIT_OWNER *) find_attrib(ENTITY_ptr, 
                                             ATTRIB_SNL_TYPE,
                                             ATTRIB_CUBIT_OWNER_TYPE);
    }
  }
}

// Virtual function called when an owner entity is being split,
// such as from a subtract operation.
// Also called when an ENTITY is copied in part, such as when
// a FACE is copied without also copying a LUMP to attach it to.
void ATTRIB_CUBIT_OWNER::split_owner(ENTITY *entity)
{
//   PRINT_INFO("split_owner called for attrib_cubit_owner on ");
  
    // must get simple attributes for new entity from CubitAttribUser
  if (!cubitOwnerData)
    return;
  
  TopologyBridge* tb = CAST_TO(cubitOwnerData, TopologyBridge);
//  CubitEntity* cep = CAST_TO(tb->topology_entity(), CubitEntity);
//   if (cep)
//     PRINT_INFO("%s %d\n",
//                cep->class_name(),
//                cep->id());
//   else
//     PRINT_INFO("unknown ENTITY.\n");
  CubitAttribUser *cau = CAST_TO(tb->topology_entity(), CubitAttribUser);
  
  if (cau != NULL && cau->num_cubit_attrib() != 0)
  {
    DLIList<CubitSimpleAttrib*> csa_list;
    cau->split_owner(csa_list);
    
      // now put those csa's on the new entity, if there are any
    if (csa_list.size() != 0)
    {
      for (int i = csa_list.size(); i > 0; i--)
      {
        CubitSimpleAttrib *csa_ptr = csa_list.get_and_step();
        assert(csa_ptr != NULL);
        API_BEGIN;
        new ATTRIB_SNL_SIMPLE ( entity, csa_ptr );
        API_END;
      }
    }
  }
  
    // in the case of a split, lose this attribute
  unhook();
  lose();
}

void ATTRIB_CUBIT_OWNER::set_copyable( bool copying_attribs )
{
  copyingAttribs = copying_attribs;   
}


logical ATTRIB_CUBIT_OWNER::copyable() const
{
  return copyingAttribs;
}


// Virtual function called when two entities are to be merged,
// as may happen with 'unite'
void ATTRIB_CUBIT_OWNER::merge_owner(ENTITY *other_entity,
                                     logical delete_this)
{
//   PRINT_INFO("merge_owner called for attrib_cubit_owner.\n");
    // Handle merging via the survivor
  if (delete_this)
    return;
  
    // Get CAU
  if (!cubitOwnerData)
    return;
  TopologyBridge* tb = CAST_TO(cubitOwnerData, TopologyBridge);
  CubitAttribUser *cau =
    (tb ? CAST_TO(tb->topology_entity(), CubitAttribUser) :
     NULL);
  if (cau == NULL || cau->num_cubit_attrib() == 0)
    return;
  
    // get the owner attribute on other_entity
  ATTRIB_CUBIT_OWNER *other_owner_att =
    (ATTRIB_CUBIT_OWNER *)find_attrib(other_entity,
                                      ATTRIB_SNL_TYPE, 
                                      ATTRIB_CUBIT_OWNER_TYPE);
  RefEntity *ref_entity = NULL;
  if (other_owner_att && other_owner_att->cubit_owner()) {
    tb = CAST_TO(other_owner_att->cubit_owner(), TopologyBridge);
    ref_entity = CAST_TO(tb->topology_entity(), RefEntity);
    assert(ref_entity != NULL);
  }
  
  cau->merge_owner(ref_entity);
}

// Virtual function called when the owner is being translated
void ATTRIB_CUBIT_OWNER::trans_owner(SPAtransf const&transform )
{
//   PRINT_INFO("trans_owner called for attrib_cubit_owner.\n");
  if (!cubitOwnerData)
    return;
  TopologyBridge* tb = CAST_TO(cubitOwnerData, TopologyBridge);
  if (tb == NULL)
    return;
  CubitAttribUser *cau = CAST_TO(tb->topology_entity(), CubitAttribUser);
  if (cau == NULL || cau->num_cubit_attrib() == 0)
    return;
  
    // dummy variables representing transform
  CubitVector translate_vec(transform.translation().x(),
                            transform.translation().y(),
                            transform.translation().z());
  
  CubitVector matrow1(transform.affine().element(0, 0),
                      transform.affine().element(0, 1),
                      transform.affine().element(0, 2));
    
  CubitVector matrow2(transform.affine().element(1, 0),
                      transform.affine().element(1, 1),
                      transform.affine().element(1, 2));
    
  CubitVector matrow3(transform.affine().element(2, 0),
                      transform.affine().element(2, 1),
                      transform.affine().element(2, 2));
    
  double scale_factor = transform.scaling();
    
  cau->transf_owner(matrow1, matrow2, matrow3, translate_vec, scale_factor);
}

//Notifies this ATTRIB that its owning ENTITY is being 
//replaced with a tolerant ENTITY. 
void ATTRIB_CUBIT_OWNER::to_tolerant_owner( ENTITY *tol_ent )
{
  //does tolerant entity already have an attribute on it?  
  ATTRIB_CUBIT_OWNER *tol_ent_owner_att =
    (ATTRIB_CUBIT_OWNER *)find_attrib(tol_ent,
                                      ATTRIB_SNL_TYPE, 
                                      ATTRIB_CUBIT_OWNER_TYPE);
  //if not...append this one to it
  if( tol_ent_owner_att == NULL )
  {
    if( cubitOwnerData )
      set_cubit_owner(tol_ent, cubitOwnerData );
  }
}

