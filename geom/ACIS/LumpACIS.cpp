//-------------------------------------------------------------------------
// Filename      : LumpACIS.cpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/02/96
//
// Owner         : Malcolm J. Panthaki
//-------------------------------------------------------------------------

// ********** BEGIN STANDARD INCLUDES      **********
#include <assert.h>
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN ACIS INCLUDES          **********
#if CUBIT_ACIS_VERSION < 1100
#include "kernel/kernapi/api/kernapi.hxx"
#include "kernel/kerndata/data/datamsc.hxx"
#include "constrct/kernapi/api/cstrapi.hxx"
#include "kernel/kernutil/tensor/tensor.hxx"
#include "kernel/kerndata/top/lump.hxx"
#include "kernel/kerndata/top/body.hxx"
#include "kernel/kerndata/lists/lists.hxx"
#include "intersct/kernapi/api/intrapi.hxx"
#include "baseutil/vector/position.hxx"
#include "baseutil/vector/unitvec.hxx"
#include "baseutil/vector/box.hxx"

#else
#include "kernapi.hxx"
#include "datamsc.hxx"
#include "cstrapi.hxx"
#include "tensor.hxx"
#include "lump.hxx"
#include "body.hxx"
#include "lists.hxx"
#include "intrapi.hxx"
#include "position.hxx"
#include "unitvec.hxx"
#include "box.hxx"
#include "mprop.hxx"
#endif
// ********** END ACIS INCLUDES            **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "attrib_cubit_owner.hpp"
#include "attrib_snl_simple.hpp"

#include "AcisQueryEngine.hpp"
#include "LumpACIS.hpp"
#include "CastTo.hpp"
#include "CubitString.hpp"
#include "BodySM.hpp"
#include "ShellSM.hpp"
#include "RefFace.hpp"
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN FORWARD DECLARATIONS   **********
class RefVolume;
// ********** END FORWARD DECLARATIONS     **********

// ********** BEGIN STATIC DECLARATIONS    **********
// ********** END STATIC DECLARATIONS      **********

// ********** BEGIN PUBLIC FUNCTIONS       **********

//-------------------------------------------------------------------------
// Purpose       : The constructor with a pointer to a LUMP. 
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/02/96
//-------------------------------------------------------------------------
LumpACIS::LumpACIS(LUMP* LUMPPtr)
    : AcisBridge(LUMPPtr)
{
    // Calculate a bounding box if there isn't one already
    //get_acis_query_engine()->bounding_box(LUMPPtr);
}

//-------------------------------------------------------------------------
// Purpose       : The destructor
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 09/06/96
//-------------------------------------------------------------------------
LumpACIS::~LumpACIS()
{}

//-------------------------------------------------------------------------
// Purpose       : Return a pointer to the LUMP associated with
//                 the object.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 09/06/96
//-------------------------------------------------------------------------
LUMP* LumpACIS::get_LUMP_ptr() const
{
  return (LUMP*)ENTITY_ptr();
}

void LumpACIS::set_LUMP_ptr(LUMP* LUMP_ptr)
{
  ENTITY_ptr(LUMP_ptr);
}

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to append a
//                 attribute to the GE. The name is attached to the 
//                 underlying solid model entity this one points to.
//
//
// Special Notes : 
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 11/21/96
//-------------------------------------------------------------------------
void LumpACIS::append_simple_attribute_virt(CubitSimpleAttrib* csattrib_ptr)
{
  AcisBridge::append_simple_attribute_virt(csattrib_ptr);
}

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to remove a simple 
//                 attribute attached to this geometry entity. The name is 
//                 removed from the underlying BODY this points to.
//
// Special Notes : 
//
// Creator       : David R. White
//
// Creation Date : 03/18/97
//-------------------------------------------------------------------------
void LumpACIS::remove_simple_attribute_virt(CubitSimpleAttrib* csattrib_ptr)
{
  AcisBridge::remove_simple_attribute_virt(csattrib_ptr);
}

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to remove all simple 
//                 attributes attached to this geometry entity.  Also
//                 removes lingering GTC attributes.
//
//
// Special Notes : 
//
// Creator       : Greg Nielson
//
// Creation Date : 07/10/98
//-------------------------------------------------------------------------
void LumpACIS::remove_all_simple_attribute_virt()
{
  AcisBridge::remove_all_simple_attribute_virt();
}

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to get the  
//                 attributes attached to this geometry entity. The name is 
//                 attached to the underlying BODY this points to.
//
// Special Notes : 
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 01/23/97
//-------------------------------------------------------------------------
CubitStatus LumpACIS::get_simple_attribute(DLIList<CubitSimpleAttrib*>&
                                           cubit_simple_attrib_list)
{
  return AcisBridge::get_simple_attribute(cubit_simple_attrib_list);
}
CubitStatus LumpACIS::get_simple_attribute(const CubitString& name,
                                       DLIList<CubitSimpleAttrib*>& list)
  { return AcisBridge::get_simple_attribute(name, list); }

//-------------------------------------------------------------------------
// Purpose       : Get the bounding box of the object.
//
// Special Notes :
//
// Creator       : jihong Ma
//
// Creation Date : 10/23/96
//-------------------------------------------------------------------------
CubitBox LumpACIS::bounding_box() const 
{
   // Calculate a bounding box if there isn't one already and return
   // the result.
    SPAbox ACIS_box = get_acis_query_engine()->bounding_box(get_LUMP_ptr());
   
   // Convert to a CubitBox and return it
   return get_acis_query_engine()->bounding_box(ACIS_box);
}

//-------------------------------------------------------------------------
// Purpose       : Get geometry modeling engine: AcisQueryEngine
//
// Special Notes :
//
// Creator       : jihong Ma
//
// Creation Date : 10/22/96
//-------------------------------------------------------------------------
GeometryQueryEngine* 
                 LumpACIS::get_geometry_query_engine() const
{
   return get_acis_query_engine();   
}                 

//-------------------------------------------------------------------------
// Purpose       : Returns the volume of the Lump
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 12/19/96
//-------------------------------------------------------------------------
double LumpACIS::measure()
{
  double volume = 0.0;
    // Create a temporary ACIS BODY using a temporary copy of this LUMP
    // and ask ACIS to compute the volume of that BODY.  The reason this
    // needs to be done is that ACIS does not provide a function to
    // compute the volume of a LUMP. This is wrapped with an API_NOP macro
    // so that ACIS rolls back to destroy the new temporary objects
    // created during this operation.
  API_NOP_BEGIN; // This "prevents" any changes from occurring to the BODY
    // upto the corresponding API_NOP_END
  
    // Make a BODY from this LUMP. Note that this modifies the input LUMP
    // which is one of the reasons for having the API_NOP macros.
  get_LUMP_ptr()->set_next(( LUMP*)NULL );
  BODY* BODYPtr = new BODY(get_LUMP_ptr()); 
  
  mass_props volume_props;
  mass_props_options volume_props_options;
  volume_props_options.set_level( VOLUME_ONLY );

  api_body_mass_props( BODYPtr, volume_props, &volume_props_options );
  volume = volume_props.get_volume();

  API_NOP_END;
  
  return volume;
}


CubitStatus LumpACIS::mass_properties( CubitVector &centroid, double &volume )
{
    // Create a temporary ACIS BODY using a temporary copy of this LUMP
    // and ask ACIS to compute the volume of that BODY.  The reason this
    // needs to be done is that ACIS does not provide a function to
    // compute the volume of a LUMP. This is wrapped with an API_NOP macro
    // so that ACIS rolls back to destroy the new temporary objects
    // created during this operation.
  API_NOP_BEGIN; // This "prevents" any changes from occurring to the BODY
    // upto the corresponding API_NOP_END
  
    // Make a BODY from this LUMP. Note that this modifies the input LUMP
    // which is one of the reasons for having the API_NOP macros.
  get_LUMP_ptr()->set_next(( LUMP*)NULL );
  BODY* BODYPtr = new BODY(get_LUMP_ptr()); 

  mass_props volume_props;
  mass_props_options volume_props_options;
  volume_props_options.set_level( VOLUME_AND_CENTROID );

  outcome result = api_body_mass_props( BODYPtr, volume_props, &volume_props_options );
  if( !result.ok() )
    return CUBIT_FAILURE;

  volume = volume_props.get_volume();
  SPAposition tmp_centroid = volume_props.get_centroid();
  centroid.set( tmp_centroid.x(), tmp_centroid.y(), tmp_centroid.z() );

  API_NOP_END;

  return CUBIT_SUCCESS;
}


//-------------------------------------------------------------------------
// Purpose       : Checks the validity of this lump
//
// Special Notes :
//
// Creator       : David R. White
//
// Creation Date : 11/20/2002
//-------------------------------------------------------------------------
int LumpACIS::validate(const CubitString &entity_name,
                       DLIList <TopologyEntity*> &bad_entities)
{
  ENTITY_LIST entity_list;
  FILE *file_ptr = NULL;
  outcome result = api_check_entity(get_LUMP_ptr(), &entity_list, file_ptr);
  if ( result.ok() && entity_list.count() == 0 )
  {
    return 0;
  }
  else if ( !result.ok() )
  {
    get_acis_query_engine()->ACIS_API_error(result);
  }
  
  if ( entity_list.count() )
  {
    convert_entity_list( entity_list, bad_entities );
    PRINT_ERROR("with '%s' in the ACIS geometry.\n",
                entity_name.c_str());
  }
  else
    return 0;

  return 1;
}

void LumpACIS::convert_entity_list( ENTITY_LIST &ent_list,
                                    DLIList <TopologyEntity*> &topo_entity_list)
{
  TopologyEntity* te_ptr = NULL;
  int i;
  for (i=0;i<ent_list.count();i++)
  {
    te_ptr = ATTRIB_CUBIT_OWNER::get_topology_entity(ent_list[i]);
    if ( te_ptr != NULL )
      topo_entity_list.append(te_ptr);
  }
}


/*
void LumpACIS::bodysms(DLIList<BodySM*> &bodies) 
{
  AcisBridge::bodysms(bodies);
}

void LumpACIS::lumps(DLIList<Lump*> &lumps)
{
  AcisBridge::lumps(lumps);
}

void LumpACIS::shellsms(DLIList<ShellSM*> &shellsms)
{
  AcisBridge::shellsms(shellsms);
}

void LumpACIS::surfaces(DLIList<Surface*> &surfaces)
{
  AcisBridge::surfaces(surfaces);
}

void LumpACIS::loopsms(DLIList<LoopSM*> &loopsms)
{
  AcisBridge::loopsms(loopsms);
}

void LumpACIS::curves(DLIList<Curve*> &curves)
{
  AcisBridge::curves(curves);
}

void LumpACIS::coedgesms(DLIList<CoEdgeSM*> &coedgesms)
{
  AcisBridge::coedgesms(coedgesms);
}

void LumpACIS::points(DLIList<Point*> &points)
{
  AcisBridge::points(points);
}
*/


void LumpACIS::get_parents_virt( DLIList<TopologyBridge*>& parents )
{
  BODY* body = get_LUMP_ptr()->body();
  if( body )
  {
    AcisBridge* acis_bridge = ATTRIB_CUBIT_OWNER::cubit_owner( body );
    TopologyBridge* bridge = dynamic_cast<TopologyBridge*>(acis_bridge);
    if (bridge)
      parents.append( bridge );
  }
}

void LumpACIS::get_children_virt( DLIList<TopologyBridge*>& children )
{
  ENTITY_LIST entities;
  api_get_shells( get_LUMP_ptr(), entities );
  ATTRIB_CUBIT_OWNER::cubit_owner( entities, children );
}


// ********** END PUBLIC FUNCTIONS         **********

// ********** BEGIN PROTECTED FUNCTIONS    **********
// ********** END PROTECTED FUNCTIONS      **********

// ********** BEGIN PRIVATE FUNCTIONS      **********
// ********** END PRIVATE FUNCTIONS        **********

// ********** BEGIN HELPER CLASSES         **********
// ********** END HELPER CLASSES           **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

// ********** BEGIN STATIC FUNCTIONS       **********
// ********** END STATIC FUNCTIONS         **********
