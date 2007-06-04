//-------------------------------------------------------------------------
// Filename      : BodyACIS.cpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/06/96
//
// Owner         : Malcolm J. Panthaki
//-------------------------------------------------------------------------

// ********** BEGIN STANDARD INCLUDES      **********
//#include <stddef.h>
//#include <iostream.h>
#include <assert.h>
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN ACIS INCLUDES          **********
#if CUBIT_ACIS_VERSION < 1100
#include "constrct/kernapi/api/cstrapi.hxx"
#include "kernel/kerndata/top/body.hxx"
#include "kernel/kerndata/top/lump.hxx"
#include "kernel/kerndata/top/shell.hxx"
#include "kernel/kerndata/top/face.hxx"
#include "kernel/kerndata/transent/transent.hxx"
#include "kernel/kerndata/geom/transfrm.hxx"
#include "kernel/kernapi/api/api.hxx"
#include "kernel/kernapi/api/kernapi.hxx"
#include "kernel/kerndata/lists/lists.hxx"
#include "intersct/kernapi/api/intrapi.hxx"
#include "baseutil/vector/vector.hxx"
#include "baseutil/vector/transf.hxx"
#include "baseutil/vector/box.hxx"
#include "operator/kernapi/api/operapi.hxx"
#include "kernel/kernutil/tensor/tensor.hxx"
#include "intersct/kernapi/api/ptcont.hxx"

#else
#include "cstrapi.hxx"
#include "body.hxx"
#include "lump.hxx"
#include "shell.hxx"
#include "face.hxx"
#include "transent.hxx"
#include "transfrm.hxx"
#include "api.hxx"
#include "kernapi.hxx"
#include "lists.hxx"
#include "intrapi.hxx"
#include "vector.hxx"
#include "transf.hxx"
#include "box.hxx"
#include "operapi.hxx"
#include "warp_api.hxx"
#include "tensor.hxx"
#include "ptcont.hxx"
#include "insanity_list.hxx"
#include "err_ent.hxx"
#endif
// ********** END ACIS INCLUDES            **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "attrib_cubit_owner.hpp"
#include "attrib_snl_simple.hpp"
#include "CubitDefines.h"
#include "CubitString.hpp"
#include "CastTo.hpp"
#include "CubitTransformMatrix.hpp"

#include "BodyACIS.hpp"
#include "BodySM.hpp"
#include "Lump.hpp"
#include "Body.hpp"

#include "GeometryQueryTool.hpp"
#include "AcisQueryEngine.hpp"
#include "CubitSimpleAttrib.hpp"

#include "CubitUtil.hpp"
#include "RefVertex.hpp"
#include "RefEdge.hpp"
#include "CoEdge.hpp"
#include "Loop.hpp"
#include "RefFace.hpp"
#include "Shell.hpp"
#include "RefVolume.hpp"

// ********** END CUBIT INCLUDES           **********

// ********** BEGIN STATIC DECLARATIONS    **********
// ********** END STATIC DECLARATIONS      **********

// ********** BEGIN PUBLIC FUNCTIONS       **********

//-------------------------------------------------------------------------
// Purpose       : A constructor with a pointer to a ACIS BODY.
//
// Special Notes :
//
// Creator       : Xuechen Liu
//
// Creation Date : 08/06/96
//-------------------------------------------------------------------------
BodyACIS::BodyACIS(BODY* BODYPtr)
    : AcisBridge(BODYPtr)
{
    // Calculate a bounding box if there isn't one already
    get_acis_query_engine()->bounding_box(BODYPtr);
}

//-------------------------------------------------------------------------
// Purpose       : The destructor.
//
// Special Notes :
//
// Creator       : Raikanta Sahu
//
// Creation Date : 09/06/96
//-------------------------------------------------------------------------
BodyACIS::~BodyACIS()
{
}

//-------------------------------------------------------------------------
// Purpose       : This function returns a void pointer that points to the
//                 "appropriate" portion of this object.  The appropriate
//                 portion is determined by the input EntityType variable.
//                 Returns NULL if the input type and the type of "this"
//                 are not related by inheritance.
//
// Special Notes : 
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 10/29/96
//-------------------------------------------------------------------------
BODY* BodyACIS::get_BODY_ptr(Body *body_ptr, CubitBoolean verbose)
{
  BodySM* OSMEPtr =
    body_ptr->get_body_sm_ptr();
  BodyACIS* bodyACISPtr = CAST_TO(OSMEPtr, BodyACIS) ;
  if(!bodyACISPtr)
  {
    if( verbose )
      PRINT_ERROR("Referenced body is not ACIS.\n"
                  "Possible incompatible geometry engine for this operation.\n");

    return (BODY *)NULL;
  }
  
  return bodyACISPtr->get_BODY_ptr();
}

BODY* BodyACIS::get_BODY_ptr(const Body *body_ptr)
{
  BodySM* OSMEPtr =
    body_ptr->get_body_sm_ptr();
  BodyACIS* bodyACISPtr = CAST_TO(OSMEPtr, BodyACIS) ;
  assert(bodyACISPtr != NULL);
  
  return bodyACISPtr->get_BODY_ptr();
}


//-------------------------------------------------------------------------
// Purpose       : Get geometry modeling engine: AcisQueryEngine
//
// Special Notes :
//
// Creator       : Stephen J. Verzi
//
// Creation Date : 02/26/97
//-------------------------------------------------------------------------
GeometryQueryEngine* BodyACIS::get_geometry_query_engine() const
{
  return get_acis_query_engine();   
}                 

//-------------------------------------------------------------------------
// Purpose       : Check if this body is a sheet body
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 08/18/03
//-------------------------------------------------------------------------
bool BodyACIS::is_sheet_body() const 
{
  BODY* BODY_ptr = get_BODY_ptr();
  if (!BODY_ptr)
  {
    PRINT_ERROR("BodyACIS without BODY at %s:%d\n", __FILE__, __LINE__ );
    assert(BODY_ptr != NULL);
    return true;
  }
    
  return is_sheet_body( BODY_ptr );
}

bool BodyACIS::is_sheet_body(BODY* BODY_ptr) 
{
  for( LUMP* LUMP_ptr = BODY_ptr->lump();
       LUMP_ptr != NULL; 
       LUMP_ptr = LUMP_ptr->next() )
  {
        
      // For each LUMP, traverse its SHELLs
    for( SHELL* SHELL_ptr = LUMP_ptr->shell();
         SHELL_ptr != NULL;
         SHELL_ptr = SHELL_ptr->next() )
    {

        // Get the FACEs of this SHELL
      for( FACE* FACE_ptr = SHELL_ptr->first_face();
           FACE_ptr != NULL;
           FACE_ptr = FACE_ptr->next_face() )
      {
        if ( FACE_ptr->sides() != DOUBLE_SIDED ||
             FACE_ptr->cont()  != BOTH_OUTSIDE )
          return false;
      }
    }
  }
  
  return true;
}

  

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to append a
//                 attribute to the OSME. The name is attached to the 
//                 underlying solid model entity this one points to.
//
//
// Special Notes : 
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 11/21/96
//-------------------------------------------------------------------------
void BodyACIS::append_simple_attribute_virt(CubitSimpleAttrib* csattrib_ptr)
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
void BodyACIS::remove_simple_attribute_virt(CubitSimpleAttrib* csattrib_ptr)
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
void BodyACIS::remove_all_simple_attribute_virt()
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
CubitStatus BodyACIS::get_simple_attribute(DLIList<CubitSimpleAttrib*>&
                                           cubit_simple_attrib_list)
{
  return AcisBridge::get_simple_attribute(cubit_simple_attrib_list);
}
CubitStatus BodyACIS::get_simple_attribute(const CubitString& name,
                                       DLIList<CubitSimpleAttrib*>& list)
  { return AcisBridge::get_simple_attribute(name, list); }


//-------------------------------------------------------------------------
// Purpose       : return the ACIS body transformation matrix
//
// Special Notes : 
//
// Creator       : sjowen
//
// Creation Date : 4/18/01
//-------------------------------------------------------------------------
CubitStatus BodyACIS::get_transforms( CubitTransformMatrix &tfm )
{
   // Get the current transformation matrix of the BODY
  if ( get_BODY_ptr()->transform() == NULL )
  {
    PRINT_INFO("Sorry, transform data was lost.\n");
    return CUBIT_FAILURE;
  }
   SPAtransf  transformation = get_BODY_ptr()->transform()->transform() ; //last_transf();
   
   // If there has been any shear, the BODY cannot be restored.
   if (transformation.shear()) 
   {
      PRINT_ERROR("Can't restore body that has been sheared\n");
      return CUBIT_FAILURE;
   }
  
   tfm.set_to_identity();
   SPAmatrix affine = transformation.affine();
   int ii,jj;
   for (ii=0; ii<3; ii++)
   {
     for (jj=0; jj<3; jj++)
     {
       tfm.set( ii, jj, affine.element(ii, jj) );
     }
   }
   
   return CUBIT_SUCCESS;
}


void BodyACIS::set_BODY_ptr(BODY* the_BODY)
{
  ENTITY_ptr(the_BODY);
}

//-------------------------------------------------------------------------
// Purpose       : This function does and api_entity_check for the body.
//
// Special Notes :
//
// Creator       : David White
//
// Creation Date : 3/6/98
//-------------------------------------------------------------------------

int BodyACIS::validate( const CubitString &user_name,
                        DLIList <TopologyEntity*> &bad_entities)
{
  insanity_list *entity_list = NULL;
  outcome result = api_check_entity (get_BODY_ptr(), entity_list );

  if ( result.ok() && entity_list == NULL )
  {
    return 0;
  }
  else if ( !result.ok() )
  {
    get_acis_query_engine()->ACIS_API_error(result);
  }
  
  if ( entity_list )
  {
    convert_entity_list( entity_list, bad_entities );
    show_bad_geom( bad_entities, user_name );
    PRINT_ERROR("with '%s' in the ACIS geometry.\n",
                user_name.c_str());
  }

  else
    return 0;

  return 1;
}
/*
void BodyACIS::bodysms(DLIList<BodySM*> &bodies) 
{
  AcisBridge::bodysms(bodies);
}

void BodyACIS::lumps(DLIList<Lump*> &lumps)
{
  AcisBridge::lumps(lumps);
}

void BodyACIS::shellsms(DLIList<ShellSM*> &shellsms)
{
  AcisBridge::shellsms(shellsms);
}

void BodyACIS::surfaces(DLIList<Surface*> &surfaces)
{
  AcisBridge::surfaces(surfaces);
}

void BodyACIS::loopsms(DLIList<LoopSM*> &loopsms)
{
  AcisBridge::loopsms(loopsms);
}

void BodyACIS::curves(DLIList<Curve*> &curves)
{
  AcisBridge::curves(curves);
}

void BodyACIS::coedgesms(DLIList<CoEdgeSM*> &coedgesms)
{
  AcisBridge::coedgesms(coedgesms);
}

void BodyACIS::points(DLIList<Point*> &points)
{
  AcisBridge::points(points);
}
*/

void BodyACIS::get_parents_virt( DLIList<TopologyBridge*>& )
{
}

void BodyACIS::get_children_virt( DLIList<TopologyBridge*>& children )
{
  ENTITY_LIST entities;
  api_get_lumps( get_BODY_ptr(), entities );
  ATTRIB_CUBIT_OWNER::cubit_owner( entities, children );
}


CubitStatus BodyACIS::mass_properties( CubitVector& result, double& volume )
{
  BODY *body_ACIS = get_BODY_ptr();
  outcome rc;
  SPAbox temp_box = AcisQueryEngine::instance()->bounding_box(body_ACIS);
  double px = temp_box.high().x() - temp_box.low().x();
  double py = temp_box.high().y() - temp_box.low().y();
  double pz = temp_box.high().z() - temp_box.low().z();
  
  SPAposition temp_position(px, py, pz);
  SPAunit_vector temp_vector(0, 0, 1);
  SPAposition cofg_position;
  tensor temp_tensor;
  double p_moments[3];
  SPAunit_vector p_axes[3];
  double real_rel_acc;

//   DELTA_STATE* ds = NULL; 
//   api_note_state(ds); 
  
  rc = api_body_mass_pr ( body_ACIS,  //body to be inspected
                          temp_position, //this and temp_vector locate the plane
                          temp_vector, //the body is projected onto
                          1, //0, calculate all props, 1, don't do inertia
                          .01, //request accuracy of 1%
                          volume, //fill this with volume
                          cofg_position, //fill this with c of g
                          temp_tensor, //fill this with inertia
                          p_moments, //fill with p_moments
                          p_axes, //fill with p_axes
                          real_rel_acc);

//   api_note_state(ds); 
//   if(ds != NULL)
//   {
//     api_change_state(ds); 
//     api_delete_ds(ds);
//   }
  

  if (!rc.ok())
  {
    PRINT_ERROR("Mass Properties calculation failed");
    return CUBIT_FAILURE;
  }
  result.set(cofg_position.x(),cofg_position.y(), cofg_position.z());
  return CUBIT_SUCCESS;
}


CubitPointContainment BodyACIS::point_containment( const CubitVector &point_coords )
{
  outcome rc;
  CubitPointContainment is_point_in;
  
  SPAposition is_in_point(point_coords.x(), point_coords.y(), point_coords.z());

  BODY *body_ACIS = get_BODY_ptr();
  enum point_containment is_in;

  rc = api_point_in_body(is_in_point, body_ACIS, is_in);

  if (!rc.ok())
  {
    PRINT_ERROR("Interior point calculation failed");
    return CUBIT_PNT_UNKNOWN;
  }

  if (is_in == point_unknown)
     is_point_in = CUBIT_PNT_UNKNOWN; // undetermined location of point
  else if (is_in == point_inside)
     is_point_in = CUBIT_PNT_INSIDE;
  else if (is_in == point_boundary)
     is_point_in = CUBIT_PNT_BOUNDARY;
  else
     is_point_in = CUBIT_PNT_OUTSIDE; // point is outside body

  return is_point_in;
}

// ********** END PUBLIC FUNCTIONS         **********

// ********** BEGIN PROTECTED FUNCTIONS    **********
// ********** END PROTECTED FUNCTIONS      **********

// ********** BEGIN PRIVATE FUNCTIONS      **********

void BodyACIS::convert_entity_list( insanity_list *ent_list,
                                    DLIList <TopologyEntity*> &topo_entity_list)
{
  int i;
  if( ent_list->count() )
  {
    insanity_list *tmp_list = ent_list;
    for(; tmp_list; tmp_list = tmp_list->next() )
    {
      insanity_data *error_data = tmp_list->data();
      ENTITY *bad_ent = error_data->get_ent();
      ENTITY_LIST bad_ents;
      if(bad_ent->identity() == ERROR_ENTITY_TYPE )
      {
        ERROR_ENTITY *error_ent = (ERROR_ENTITY*)bad_ent;
        bad_ents.add( error_ent->get_owner(0) ); 
        bad_ents.add( error_ent->get_owner(1) ); 
      }
      else
        bad_ents.add( bad_ent );
      
      bad_ents.init();

      PRINT_INFO("%s: ", error_data->get_message() );

      while( (bad_ent = bad_ents.next()) != NULL )
      {
        TopologyEntity* te_ptr = NULL;
        te_ptr = ATTRIB_CUBIT_OWNER::get_topology_entity( bad_ent );
        if( te_ptr )
        {
          RefEntity *ref_ent = CAST_TO( te_ptr, RefEntity );
          topo_entity_list.append_unique(te_ptr);
          if( ref_ent )
            PRINT_INFO(" %s %d", ref_ent->class_name(), ref_ent->id() );
        }
        else
          PRINT_INFO("%s",bad_ent->type_name() );
      }
      PRINT_INFO("\n");
    }
  }
}

 
CubitStatus BodyACIS::show_bad_geom(DLIList <TopologyEntity*> &ent_list,
                                    const CubitString &user_name)
{
  DLIList<RefVertex*> vertex_list;
  DLIList<RefEdge*> edge_list;
  DLIList<RefEdge*> curve_coedge_list;
  DLIList<CoEdge*> coedge_list;
  DLIList<Loop*> loop_list;
  DLIList<RefFace*> face_list;
  DLIList<Shell*> shell_list;
  DLIList<RefVolume*> volume_list;
  int i;

  // We need an encapsulated class to handle this sort of thing...
  for (i=0;i<ent_list.size();i++)
  {
    TopologyEntity* TE_ptr = ent_list.get_and_step();
    RefVertex* ref_vertex;
    RefEdge* ref_edge;
    CoEdge* co_edge;
    Loop* loop;
    RefFace* face;
    Shell* shell;
    RefVolume* volume;

    if( (ref_vertex = CAST_TO( TE_ptr, RefVertex ) ) != NULL )
      vertex_list.append_unique(ref_vertex);
    else if( (ref_edge = CAST_TO( TE_ptr, RefEdge ) ) != NULL )
      edge_list.append_unique(ref_edge);
    else if( (co_edge = CAST_TO( TE_ptr, CoEdge ) ) != NULL )
      coedge_list.append_unique(co_edge);
    else if( (loop = CAST_TO( TE_ptr, Loop ) ) != NULL )
      loop_list.append_unique(loop);
    else if( (face = CAST_TO( TE_ptr, RefFace ) ) != NULL )
      face_list.append_unique(face);
    else if( (shell = CAST_TO( TE_ptr, Shell ) ) != NULL )
      shell_list.append_unique(shell);
    else if( (volume = CAST_TO( TE_ptr, RefVolume ) ) != NULL )
      volume_list.append_unique(volume);
  }

  char pre[100];

  if( vertex_list.size() )
  {
    sprintf( pre, "Found %d bad vertices in %s. The vertices: ", vertex_list.size(),
      user_name.c_str() );
    DLIList<CubitEntity*> cubit_list;
    CAST_LIST( vertex_list, cubit_list, CubitEntity );
    CubitUtil::list_entity_ids( pre, cubit_list );
  }

  if( edge_list.size() )
  {
    sprintf( pre, "Found %d bad curves in %s. The curves: ", edge_list.size(), 
      user_name.c_str() );
    DLIList<CubitEntity*> cubit_list;
    CAST_LIST( edge_list, cubit_list, CubitEntity );
    CubitUtil::list_entity_ids( pre, cubit_list );
  }

  for( i=0; i<coedge_list.size(); i++ )
  {
    RefEdge* ref_edge_ptr = coedge_list.get_and_step()->get_ref_edge_ptr();
    curve_coedge_list.append_unique( ref_edge_ptr );
  }
  if( coedge_list.size() )
  {
    sprintf( pre, "Found %d bad coedges in %s. The curves: ", coedge_list.size(), 
      user_name.c_str() );
    DLIList<CubitEntity*> cubit_list;
    CAST_LIST( curve_coedge_list, cubit_list, CubitEntity );
    CubitUtil::list_entity_ids( pre, cubit_list );
  }

  // LOOPS
  Loop* loop_ptr;
  DLIList<RefEdge*> curve_loop_list;
  for( i=0; i<curve_loop_list.size(); i++ )
  {
    loop_ptr = loop_list.get_and_step();
    DLIList<RefEdge*> tmp_curve_list;
    loop_ptr->ordered_ref_edges( tmp_curve_list );
    curve_loop_list.merge_unique( tmp_curve_list );
  }
  curve_loop_list.reset();
  if( curve_loop_list.size() )
  {
    sprintf( pre, "Found %d bad loops in %s. The curves: ", loop_list.size(), user_name.c_str() );
    DLIList<CubitEntity*> cubit_list;
    CAST_LIST( curve_loop_list, cubit_list, CubitEntity );
    CubitUtil::list_entity_ids( pre, cubit_list );
  }

  // SURFACES
  if( face_list.size() )
  {
    sprintf( pre, "Found %d bad surfaces in %s. The surfaces: ", face_list.size(), user_name.c_str() );
    DLIList<CubitEntity*> cubit_list;
    CAST_LIST( face_list, cubit_list, CubitEntity );
    CubitUtil::list_entity_ids( pre, cubit_list );
  }

  // SHELLS
  Shell* shell_ptr;
  DLIList<RefFace*> face_shell_list;
  for( i=0; i<shell_list.size(); i++ )
  {
    shell_ptr = shell_list.get_and_step();
    DLIList<RefFace*> tmp_surface_list;
    shell_ptr->ref_faces( tmp_surface_list );
    face_shell_list.merge_unique( tmp_surface_list );
  }
  if( face_shell_list.size() )
  {
    sprintf( pre, "Found %d bad shells in %s. The surfaces: ", shell_list.size(), user_name.c_str() );
    DLIList<CubitEntity*> cubit_list;
    CAST_LIST( face_shell_list, cubit_list, CubitEntity );
    CubitUtil::list_entity_ids( pre, cubit_list );
  }

  // LUMPS
  if( volume_list.size() )
  {
    sprintf( pre, "Found %d bad volumes in %s. The volumes: ", volume_list.size(), user_name.c_str() );
    DLIList<CubitEntity*> cubit_list;
    CAST_LIST( volume_list, cubit_list, CubitEntity );
    CubitUtil::list_entity_ids( pre, cubit_list );
  }
  
  return CUBIT_SUCCESS;
}


// ********** END PRIVATE FUNCTIONS        **********

// ********** BEGIN HELPER CLASSES         **********
// ********** END HELPER CLASSES           **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

// ********** BEGIN STATIC FUNCTIONS       **********
// ********** END STATIC FUNCTIONS         **********


