//-------------------------------------------------------------------------
// Filename      : AcisToolUtil.cpp
//
// Purpose       : Provide old RefEntity-based functions that were removed
//                 from Acis{Query|Modify}Engine for use by misc. Acis*Tool
//                 classes that still have a RefEntity-based interface.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/23/03
//-------------------------------------------------------------------------

#include "DLIList.hpp"

#include "AcisToolUtil.hpp"
#include "AcisQueryEngine.hpp"
#include "AcisModifyEngine.hpp"
#include "GeometryQueryTool.hpp"
#include "GeometryModifyTool.hpp"

#include "Body.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"

#include "BodyACIS.hpp"
#include "SurfaceACIS.hpp"
#include "CurveACIS.hpp"
#include "attrib_cubit_owner.hpp"

Body* AcisToolUtil::get_new_Body( Body* old_body, BODY* old_BODY,
                                  BODY* new_BODY, bool keep_old,
                                  bool topo_check, bool delete_old )
{
  BodySM* body_sm = old_body->get_body_sm_ptr();
  if (!body_sm)
  {
    PRINT_ERROR("Body %d is invalid.  No attached BodySM.\n", old_body->id());
    return 0;
  }
  
  BodySM* new_bodysm = AcisModifyEngine::instance()->
    get_new_Body( body_sm, old_BODY, new_BODY, keep_old, topo_check, delete_old);
  if (!new_bodysm)
    return 0;
  
  DLIList<Body*> input_bodies(1), result_bodies(1);
  DLIList<BodySM*> new_bodies(1);
  new_bodies.append(new_bodysm);
  input_bodies.append(old_body);
  GeometryModifyTool::instance()->finish_sm_op( input_bodies, new_bodies, result_bodies );
  return result_bodies.size() ? result_bodies.get() : 0;
}



CubitStatus AcisToolUtil::get_ACIS_Curves( DLIList<RefEdge*>& edge_list,
                                           DLIList<CurveACIS*>& curve_list )
{
  DLIList<TopologyBridge*> bridge_list(2);
  RefEdge* ref_edge_ptr;
  CurveACIS* acis_curve_ptr;
  
  edge_list.reset();
  for (int i = edge_list.size(); i--; )
  {
    ref_edge_ptr = edge_list.get_and_step();
    bridge_list.clean_out();
    ref_edge_ptr->bridge_manager()->get_bridge_list(bridge_list);
    
    acis_curve_ptr = 0;
    while (bridge_list.size() && !acis_curve_ptr)
      acis_curve_ptr = dynamic_cast<CurveACIS*>(bridge_list.pop());
    
    if (!acis_curve_ptr)
    {
      PRINT_ERROR("Non-ACIS RefEdge %d at %s:%d\n", 
        ref_edge_ptr->id(), __FILE__, __LINE__ );
      return CUBIT_FAILURE;
    }
    
    curve_list.append( acis_curve_ptr );
  }
  
  return CUBIT_SUCCESS;
}



CubitStatus AcisToolUtil::get_ACIS_Surfaces( DLIList<RefFace*>& face_list,
                                             DLIList<SurfaceACIS*>& surf_list )
{
  DLIList<TopologyBridge*> bridge_list(2);
  RefFace* ref_face_ptr;
  SurfaceACIS* acis_surf_ptr;
  
  face_list.reset();
  for (int i = face_list.size(); i--; )
  {
    ref_face_ptr = face_list.get_and_step();
    bridge_list.clean_out();
    ref_face_ptr->bridge_manager()->get_bridge_list(bridge_list);
    
    acis_surf_ptr = 0;
    while (bridge_list.size() && !acis_surf_ptr)
      acis_surf_ptr = dynamic_cast<SurfaceACIS*>(bridge_list.pop());
    
    if (!acis_surf_ptr)
    {
      PRINT_ERROR("Non-ACIS RefFace %d at %s:%d\n", 
        ref_face_ptr->id(), __FILE__, __LINE__ );
      return CUBIT_FAILURE;
    }
    
    surf_list.append( acis_surf_ptr );
  }
  
  return CUBIT_SUCCESS;
}




CubitStatus AcisToolUtil::get_copied_FACES_of_body( 
                                  DLIList<RefFace*>& ref_face_list,
                                  DLIList<FACE*>& FACE_list,
                                  BODY*& copied_body_ptr ) 
{
  DLIList<RefFace*> removed_face_list;
  return get_copied_FACES_of_body( ref_face_list, FACE_list, removed_face_list, copied_body_ptr);
}



CubitStatus AcisToolUtil::get_copied_FACES_of_body( 
                                  DLIList<RefFace*>& ref_face_list,
                                  DLIList<FACE*>& FACE_list,
                                  DLIList<RefFace*>& removed_face_list,
                                  BODY*& copied_body_ptr ) 
{
  DLIList<SurfaceACIS*> surf_list(ref_face_list.size()), removed_surfaces;
  if (!get_ACIS_Surfaces( ref_face_list, surf_list ))
    return CUBIT_FAILURE;
  
  DLIList<SurfaceACIS*> copied_surf_list(surf_list);
  if (!AcisModifyEngine::instance()->get_copied_FACES_of_body(
        copied_surf_list, FACE_list, removed_surfaces, copied_body_ptr ))
    return CUBIT_FAILURE;
  
  copied_surf_list.reset();
  surf_list.reset();
  ref_face_list.reset();
  
  if (copied_surf_list.size())
  {
    for (int i = ref_face_list.size(); i--; )
    {
      SurfaceACIS* surf_ptr = surf_list.get_and_step();
      if (surf_ptr == copied_surf_list.get())
      {
        copied_surf_list.step();
        ref_face_list.step();
      }
      else
      {
        RefFace* dead = ref_face_list.remove();
        removed_face_list.append(dead);
        assert(removed_surfaces.move_to(surf_ptr) && removed_surfaces.extract());
      }
    } 
  }
  else
  {
    removed_face_list += ref_face_list;
    ref_face_list.clean_out();
    removed_surfaces.clean_out();
  }
  assert(removed_surfaces.size() == 0);
  
  AcisBridge* bridge_ptr = ATTRIB_CUBIT_OWNER::cubit_owner(copied_body_ptr);
  BodySM* body_sm_ptr = dynamic_cast<BodySM*>(bridge_ptr);
  GeometryQueryTool::instance()->make_Body(body_sm_ptr);
  
  return CUBIT_SUCCESS;
}



Body* AcisToolUtil::get_body_of_ENTITY(ENTITY* ENTITY_ptr)
{
  BodySM* body_sm = AcisQueryEngine::instance()->get_body_sm_of_ENTITY(ENTITY_ptr);
  if (!body_sm)
    return 0;
  
  TopologyEntity* ent_ptr = body_sm->topology_entity();
  return dynamic_cast<Body*>(ent_ptr);
}



