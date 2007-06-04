//- Class:       AcisSurfaceTool
//- Description: Functions to create surfaces
// 
//- Owner: Steve Storm      
//- Checked by:
//- Version: $Id:

// ********** BEGIN ACIS INCLUDES             **********
#if CUBIT_ACIS_VERSION < 1100
#include "kernel/acis.hxx"
#include "kernel/kernapi/api/api.hxx"
#include "kernel/kernapi/api/kernapi.hxx"
#include "kernel/kerndata/geom/transfrm.hxx"
#include "kernel/kerndata/top/edge.hxx"
#include "kernel/kerndata/top/loop.hxx"
#include "kernel/kerndata/top/coedge.hxx"
#include "kernel/kerndata/top/face.hxx"
#include "kernel/kerndata/top/wire.hxx"
#include "kernel/kerndata/top/body.hxx"
#include "kernel/kerndata/lists/lists.hxx"
#include "kernel/kernutil/law/law.hxx"
#include "kernel/spline/api/spl_api.hxx"
#include "intersct/kernapi/api/intrapi.hxx"
#include "constrct/kernapi/api/cstrapi.hxx"
#include "skin/kernapi/api/skinapi.hxx"
#include "skin/sg_husk/skin/skin.hxx"
#include "offset/kernapi/api/ofstapi.hxx"
#include "cover/kernapi/api/coverapi.hxx"
#include "kernel/geomhusk/bnd_crv.hxx"
#include "kernel/geomhusk/getowner.hxx"
#include "baseutil/vector/interval.hxx"
#include "intersct/kerngeom/curve/curextnd.hxx"
#include "kernel/kerngeom/curve/extend.hxx"
#include "kernel/kerngeom/curve/curdef.hxx"
#include "kernel/kerndata/geom/curve.hxx"
#include "kernel/kerndata/geom/point.hxx"
#include "kernel/kerndata/geom/surface.hxx"
#include "kernel/kerndata/top/shell.hxx"
#include "kernel/kerndata/top/lump.hxx"
#include "kernel/kerndata/top/vertex.hxx"
#include "kernel/kerndata/geom/cnstruct.hxx"

#include "baseutil/vector/position.hxx"
#include "baseutil/vector/unitvec.hxx"
#include "baseutil/logical.h"

#ifdef ACIS_LOCAL_OPS
#include "lop_husk/api/lop_api.hxx"
#endif
#include "boolean/kernapi/api/boolapi.hxx"
#include "kernel/geomhusk/copyent.hxx" // copy single entity

#else
#include "acis.hxx" 
#include "api.hxx" 
#include "kernapi.hxx" 
#include "transfrm.hxx" 
#include "edge.hxx" 
#include "loop.hxx" 
#include "coedge.hxx" 
#include "face.hxx" 
#include "wire.hxx"
#include "body.hxx" 
#include "lists.hxx" 
#include "law.hxx" 
#include "spl_api.hxx" 
#include "intrapi.hxx" 
#include "cstrapi.hxx" 
#include "skinapi.hxx"
#include "skin.hxx"
#include "ofstapi.hxx"
#include "coverapi.hxx"
#include "bnd_crv.hxx"
#include "getowner.hxx"
#include "interval.hxx"
#include "curextnd.hxx"
#include "extend.hxx"
#include "curdef.hxx"
#include "curve.hxx"
#include "point.hxx"
#include "surface.hxx"
#include "shell.hxx"
#include "lump.hxx"
#include "vertex.hxx"
#include "cnstruct.hxx" 
#include "position.hxx"
#include "unitvec.hxx" 
#include "logical.h"

#ifdef ACIS_LOCAL_OPS
  #include "lop_api.hxx"
#endif
#include "boolapi.hxx" 

#include "copyent.hxx" // copy single entity
#endif

// ********** END ACIS INCLUDES               **********

// ********** BEGIN CUBIT INCLUDES            **********
#include "AcisSurfaceTool.hpp"
#include "CubitUtil.hpp"
#include "RefEntity.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"

#include "AcisQueryEngine.hpp"
#include "GeometryQueryTool.hpp"
#include "Body.hpp"
#include "BodyACIS.hpp"
#include "CurveACIS.hpp"
#include "SurfaceACIS.hpp"

#ifdef ACIS_HEALER
  #include "AcisHealerTool.hpp"
#endif

#include "AcisTweakTool.hpp" // REMOVE THIS WHEN FIXED (see create_offset_body)

#include "DLIList.hpp"

// For multiple surfaces
#include "Loop.hpp"
#include "RefVertex.hpp"
#include "CoEdge.hpp"

#include "attrib_cubit_owner.hpp"

#include "AcisModifyEngine.hpp"
#include "GeometryModifyTool.hpp"

// ********** END CUBIT INCLUDES              **********

AcisSurfaceTool* AcisSurfaceTool::instance_ = 0;

// Method: instance
// provides access to the unique model for this execution.
// sets up this instance on first access
AcisSurfaceTool* AcisSurfaceTool::instance()
{
  if (instance_ == 0) {
    instance_ = new AcisSurfaceTool();
  }
  return instance_;
}

AcisSurfaceTool::AcisSurfaceTool() 
{
}

AcisSurfaceTool::~AcisSurfaceTool()
{

}

CubitStatus
AcisSurfaceTool::create_net_surface( DLIList<Curve*>& u_curves_in, 
                                     DLIList<Curve*>& v_curves_in, 
                                     BodySM*& new_body,
                                     double net_tol, 
                                     CubitBoolean heal,
                                     CubitBoolean verbose)
{
  new_body = NULL;
  int num_u = u_curves_in.size();
  if( num_u < 2 ) 
  {
    PRINT_ERROR( "Need at least 2 curves in the u-direction\n" );
    return CUBIT_FAILURE;
  }
  
  int num_v = v_curves_in.size();
  if( num_v < 2  )
  {
    PRINT_ERROR( "Need at least 2 curves in the v-direction\n" );
    return CUBIT_FAILURE;
  }
  
  DLIList<CurveACIS*> u_curves(num_u), v_curves(num_v);
  CAST_LIST( u_curves_in, u_curves, CurveACIS );
  CAST_LIST( v_curves_in, v_curves, CurveACIS );
  if (u_curves.size() != num_u || v_curves.size() != num_v) 
  {
    int diff = num_u - u_curves.size() + num_v - v_curves.size();
    PRINT_ERROR( "%d input Curves are not ACIS Curves.\n", diff);
    return CUBIT_FAILURE;
  }

  int i, j;
  EDGE* copied_EDGE_ptr;
  outcome result;
  
  // Create the wires from the edges 
  BODY** u_body = new BODY*[num_u];
  BODY** v_body = new BODY*[num_v];
  
  CurveACIS* curve_ptr;
  EDGE* EDGE_ptr;
  u_curves.reset();
  int num_bodies_u = 0;
  int num_bodies_v = 0;
  for( i=0; i<num_u; i++ )
  {
    copied_EDGE_ptr = 0;
    u_body[i] = 0;
    curve_ptr = u_curves.get_and_step();
    
    // Get the EDGE ptr from the RefEdge
    EDGE_ptr = curve_ptr->get_EDGE_ptr();
    if( EDGE_ptr == NULL )
    {
      PRINT_ERROR( "Invalid ACIS curve\n");
      break;
    }
    
    result = api_edge( EDGE_ptr, copied_EDGE_ptr );
    if (!result.ok())
    {
      PRINT_ERROR( "Unable to copy curve -- aborting\n" );
      AcisModifyEngine::instance()->get_acis_query_engine()->ACIS_API_error (result);
      break;
    }
    
    result = api_make_ewire( 1, &copied_EDGE_ptr, u_body[i]); 
    
    if( !result.ok() )
    {
      AcisQueryEngine::instance()->ACIS_API_error(result);
      PRINT_ERROR( "unable make ACIS WIRE BODY from curve -- aborting\n");
      break;
    }
    
    // Make sure the wire doesn't self-intersect
    result = api_check_wire_self_inters(u_body[i]);
    
    if( !result.ok() )
    {
      AcisQueryEngine::instance()->ACIS_API_error(result);
      PRINT_ERROR( "curve self-intersects -- aborting\n" );
      break;
    } 
  }
  
  if (i < num_u)  // error, clean up and return
  {
    for (j = 0; j < i; j++)
      api_delent( u_body[j] );
    if (u_body[i])
      api_delent( u_body[i] );
    else if(copied_EDGE_ptr)
      api_delent( copied_EDGE_ptr );
    delete [] u_body;
    delete [] v_body;
    return CUBIT_FAILURE;
  }
  
  v_curves.reset();
  for( i=0; i<num_v; i++ )
  {
    copied_EDGE_ptr = 0;
    v_body[i] = 0;
    curve_ptr = v_curves.get_and_step();

    // Get the EDGE ptr from the RefEdge
    EDGE_ptr = curve_ptr->get_EDGE_ptr();
    if( EDGE_ptr == NULL )
    {
      PRINT_ERROR( "Invalid ACIS curve\n" );
      break;
    }

    result = api_edge( EDGE_ptr, copied_EDGE_ptr );
    if (!result.ok())
    {
      PRINT_ERROR( "Unable to copy curve -- aborting\n" );
      AcisQueryEngine::instance()->ACIS_API_error (result);
      break;
    }
    
    result = api_make_ewire( 1, &copied_EDGE_ptr, v_body[i]); 
    
    if( !result.ok() )
    {
      AcisQueryEngine::instance()->ACIS_API_error(result);
      PRINT_ERROR( "unable make WIRE BODY from curve -- aborting\n" );
      break;
    }
    
    // Make sure the wire doesn't self-intersect
    result = api_check_wire_self_inters(v_body[i]);
    
    if( !result.ok() )
    {
      AcisQueryEngine::instance()->ACIS_API_error(result);
      PRINT_ERROR( "curve self-intersects -- aborting\n");
      break;
    } 
  } 
  
  if (i < num_v)  // error, clean up and return
  {
    for (j = 0; j < num_u; j++)
      api_delent( u_body[j] );
    for (j = 0; j < i; j++)
      api_delent( v_body[j] );
    if (v_body[i])
      api_delent( u_body[i] );
    else if(copied_EDGE_ptr)
      api_delent( copied_EDGE_ptr );
    delete [] u_body;
    delete [] v_body;
    return CUBIT_FAILURE;
  }
  
  BODY *net_body = NULL; 
  
  if(verbose)
     PRINT_INFO( "Creating net surface...\n" );
  result = api_net_wires( num_v, v_body, num_u, u_body, net_body, 
                          TRUE, TRUE, net_tol); 

  // Cleanup memory
  for( i=0; i<num_bodies_u; i++ )
    api_delent( u_body[i] );
  for( i=0; i<num_bodies_v; i++ )
    api_delent( v_body[i] );
  delete [] u_body;
  delete [] v_body;
  
  if( !result.ok() )
  {
    AcisQueryEngine::instance()->ACIS_API_error(result);
    PRINT_ERROR( "unable to create net surface\n" );
    return CUBIT_FAILURE;
  } 

#ifdef ACIS_HEALER
  // Heal the body, if requested
  if( heal == CUBIT_TRUE )
  {
      // Now heal the combined body
    //PRINT_INFO("Healing net surface...\n");
    if( AcisHealerTool::instance()->init_BODY_for_healing( net_body ) == CUBIT_SUCCESS )
    {
      int percent_before, percent_after, number_splines_simplified;
      if( AcisHealerTool::instance()->heal_BODY( net_body, percent_before,
                                                 percent_after, number_splines_simplified ) == CUBIT_FAILURE )
        PRINT_ERROR( "Error healing the net body\n" );
      //else
        //PRINT_INFO( "Healed the net body to %d%% good geometry.\n", percent_after );
      AcisHealerTool::instance()->end_BODY_for_healing( net_body );
    }
  }
#endif
  
  // Make a new body
  new_body = AcisQueryEngine::instance()->populate_topology_bridges(net_body);
  return CUBIT_SUCCESS;
}

CubitStatus
AcisSurfaceTool::create_net_surface( BodySM*& new_body, 
                                     DLIList<DLIList<CubitVector*>*> &vec_lists_u, 
                                     DLIList<DLIList<CubitVector*>*> &vec_lists_v, 
                                     double net_tol, CubitBoolean heal,
                                     CubitBoolean verbose)
{
  new_body = NULL;
  int i, j;

  int num_v = vec_lists_v.size();
  int num_u = vec_lists_u.size();

  outcome result;
  DLIList<CubitVector*>* vec_list_ptr;

  // Create EDGE's from the node lists
  EDGE** v_edge = new EDGE*[vec_lists_v.size()]; 
  EDGE** u_edge = new EDGE*[vec_lists_u.size()]; 

  // Make splines
  //CubitVector coords;
  vec_lists_v.reset();
  int num_splines_v = 0;
  int num_splines_u = 0;
  for( i=0; i<vec_lists_v.size(); i++ )
  {
    vec_list_ptr = vec_lists_v.get_and_step();
    int num_pnts = vec_list_ptr->size();
    SPAposition* pos_array = new SPAposition[num_pnts];
    
    // Do the vector_list of points
    vec_list_ptr->reset();
    for( j=0; j<num_pnts; j++)
    {
      // Set the coordinates of the SPAposition object
      pos_array[j].set_x(vec_list_ptr->get()->x());
      pos_array[j].set_y(vec_list_ptr->get()->y());
      pos_array[j].set_z(vec_list_ptr->get()->z());
      vec_list_ptr->step();
    }

    result = api_curve_spline(num_pnts, pos_array, NULL,NULL, v_edge[i]); 
    if( !result.ok() )
    {
      AcisQueryEngine::instance()->ACIS_API_error(result);
      PRINT_ERROR( "unable to make spline from node list\n" );
      delete [] pos_array;
      for( i=0; i<num_splines_u; i++ )
        api_delent( u_edge[i] );
      for( i=0; i<num_splines_v; i++ )
        api_delent( v_edge[i] );
      delete [] u_edge;
      delete [] v_edge;
      return CUBIT_FAILURE;
    }
    num_splines_u++;
    delete [] pos_array;
  }
  vec_lists_u.reset();
  for( i=0; i<vec_lists_u.size(); i++ )
  {
    vec_list_ptr = vec_lists_u.get_and_step();
    int num_pnts = vec_list_ptr->size();
    SPAposition* pos_array = new SPAposition[num_pnts];
    
    // Do the vector_list of points
    vec_list_ptr->reset();
    for( j=0; j<num_pnts; j++)
    {
      // Set the coordinates of the SPAposition object
      pos_array[j].set_x(vec_list_ptr->get()->x());
      pos_array[j].set_y(vec_list_ptr->get()->y());
      pos_array[j].set_z(vec_list_ptr->get()->z());
      vec_list_ptr->step();
    }

    result = api_curve_spline(num_pnts, pos_array, NULL,NULL, u_edge[i]); 
    if( !result.ok() )
    {
      AcisQueryEngine::instance()->ACIS_API_error(result);
      PRINT_ERROR( "unable to make spline from node list\n" );
      delete [] pos_array;
      for( i=0; i<num_splines_u; i++ )
        api_delent( u_edge[i] );
      for( i=0; i<num_splines_v; i++ )
        api_delent( v_edge[i] );
      delete [] u_edge;
      delete [] v_edge;
      return CUBIT_FAILURE;
    }
    num_splines_v++;
    delete [] pos_array;
  }

  // Create the wires from the edges 
  BODY** u_body = new BODY*[vec_lists_u.size()];
  BODY** v_body = new BODY*[vec_lists_v.size()];
  
  int num_bodies_u = 0;
  int num_bodies_v = 0;
  for( i=0; i<vec_lists_u.size(); i++ )
  {
    result = api_make_ewire( 1, &u_edge[i], u_body[i]); 
    
    if( !result.ok() )
    {
      AcisQueryEngine::instance()->ACIS_API_error(result);
      PRINT_ERROR( "unable make WIRE BODY from curve, aborting\n");
      for( i=0; i<num_bodies_u; i++ )
        api_delent( u_body[i] );
      for( i=0; i<num_bodies_v; i++ )
        api_delent( v_body[i] );
      delete [] u_body;
      delete [] v_body;
      // Deleting the bodies removed some of the edges...
      for( i=num_bodies_u; i<num_splines_u; i++ )
        api_delent( u_edge[i] );
      for( i=num_bodies_v; i<num_splines_v; i++ )
        api_delent( v_edge[i] );
      delete [] u_edge;
      delete [] v_edge;
      return CUBIT_FAILURE;
    }
    
    // Make sure the wire doesn't self-intersect
    result = api_check_wire_self_inters(u_body[i]);
    
    if( !result.ok() )
    {
      AcisQueryEngine::instance()->ACIS_API_error(result);
      PRINT_ERROR( "curve self-intersects, aborting\n" );
      for( i=0; i<num_bodies_u; i++ )
        api_delent( u_body[i] );
      for( i=0; i<num_bodies_v; i++ )
        api_delent( v_body[i] );
      delete [] u_body;
      delete [] v_body;
      // Deleting the bodies removed some of the edges...
      for( i=num_bodies_u; i<num_splines_u; i++ )
        api_delent( u_edge[i] );
      for( i=num_bodies_v; i<num_splines_v; i++ )
        api_delent( v_edge[i] );
      delete [] u_edge;
      delete [] v_edge;
      return CUBIT_FAILURE;
    } 
    num_bodies_u++;
  }

  for( i=0; i<vec_lists_v.size(); i++ )
  {
    result = api_make_ewire( 1, &v_edge[i], v_body[i]);
    
    if( !result.ok() )
    {
      AcisQueryEngine::instance()->ACIS_API_error(result);
      PRINT_ERROR( "unable make WIRE BODY from curve, aborting\n" );
      for( i=0; i<num_bodies_u; i++ )
        api_delent( u_body[i] );
      for( i=0; i<num_bodies_v; i++ )
        api_delent( v_body[i] );
      delete [] u_body;
      delete [] v_body;
      // Deleting the bodies removed some of the edges...
      for( i=num_bodies_u; i<num_splines_u; i++ )
        api_delent( u_edge[i] );
      for( i=num_bodies_v; i<num_splines_v; i++ )
        api_delent( v_edge[i] );
      delete [] u_edge;
      delete [] v_edge;
      return CUBIT_FAILURE;
    }
    
    // Make sure the wire doesn't self-intersect
    result = api_check_wire_self_inters(v_body[i]);
    
    if( !result.ok() )
    {
      AcisQueryEngine::instance()->ACIS_API_error(result);
      PRINT_ERROR( "curve self-intersects, aborting\n" );
      for( i=0; i<num_bodies_u; i++ )
        api_delent( u_body[i] );
      for( i=0; i<num_bodies_v; i++ )
        api_delent( v_body[i] );
      delete [] u_body;
      delete [] v_body;
      // Deleting the bodies removed some of the edges...
      for( i=num_bodies_u; i<num_splines_u; i++ )
        api_delent( u_edge[i] );
      for( i=num_bodies_v; i<num_splines_v; i++ )
        api_delent( v_edge[i] );
      delete [] u_edge;
      delete [] v_edge;
      return CUBIT_FAILURE;
    } 
    num_bodies_v++;
  } 

  BODY *net_body = NULL; 
 
  // Attempt to create the net
  if(verbose)
     PRINT_INFO( "Creating net surface...\n" );
  result = api_net_wires( num_v, v_body, num_u, u_body,  
                          net_body, TRUE, TRUE, net_tol ); 

  for( i=0; i<num_bodies_u; i++ )
    api_delent( u_body[i] );
  for( i=0; i<num_bodies_v; i++ )
    api_delent( v_body[i] );
  delete [] u_edge;
  delete [] v_edge;
  delete [] u_body;
  delete [] v_body;
  
  if( !result.ok() )
  {
    AcisQueryEngine::instance()->ACIS_API_error(result);
    PRINT_ERROR("Unable to create net surface approximating surfaces\n");
    return CUBIT_FAILURE;
  } 

#ifdef ACIS_HEALER
  // Heal the body, if requested
  if( heal == CUBIT_TRUE )
  {
      // Now heal the combined body
    //PRINT_INFO("Healing net surface...\n");
    if( AcisHealerTool::instance()->init_BODY_for_healing( net_body ) == CUBIT_SUCCESS )
    {
      int percent_before, percent_after, number_splines_simplified;
      if( AcisHealerTool::instance()->heal_BODY( net_body, percent_before,
                                                 percent_after, number_splines_simplified ) == CUBIT_FAILURE )
        PRINT_ERROR( "Error healing the net body\n" );
      //else
        //PRINT_INFO( "Healed the net body to %d%% good geometry.\n", percent_after );
      AcisHealerTool::instance()->end_BODY_for_healing( net_body );
    }
  }
#endif
  
  // Make a new body
  new_body = AcisQueryEngine::instance()->populate_topology_bridges(net_body);
  return CUBIT_SUCCESS;
}

CubitStatus
AcisSurfaceTool::create_offset_surface( Surface* input_surface_ptr, 
                                        BodySM*& new_body, 
                                        double offset_distance )
{
  SurfaceACIS* surface_ptr = dynamic_cast<SurfaceACIS*>(input_surface_ptr);
  if (!surface_ptr)
  {
    PRINT_ERROR("Non-ACIS Surface in AcisSurfaceTool::create_offset_surface\n");
    return CUBIT_FAILURE;
  }

  outcome result;
  FACE* FACE_ptr;
  FACE* new_FACE_ptr;

  // Get the FACE to offset
  FACE_ptr = surface_ptr->get_FACE_ptr();
  if( FACE_ptr == NULL )
  {
    PRINT_ERROR( "Invalid ACIS surface\n" );
    return CUBIT_FAILURE;
  }

  result = api_offset_face ( FACE_ptr, offset_distance, new_FACE_ptr );

  if( new_FACE_ptr == NULL )
  {
    AcisQueryEngine::instance()->ACIS_API_error(result);
    PRINT_ERROR( "unablet to offset surface\n" );
    return CUBIT_FAILURE;
  }

  // Make a sheet body out of the face
  FACE *face_list[1];
  BODY *sheet_body;
  face_list[0] = new_FACE_ptr;
  result = api_sheet_from_ff(1, face_list, sheet_body);
  if (!result.ok())
  {
    AcisQueryEngine::instance()->ACIS_API_error (result);
    PRINT_ERROR("Couldn't build a body from new offset surface.\n" );
    api_delent( new_FACE_ptr );
    return CUBIT_FAILURE;
  }
  result  = api_body_to_2d(sheet_body);
  if (!result.ok())
  {
    PRINT_ERROR("Couldn't build a body from new offset surface.\n" );
    api_delent( new_FACE_ptr );
    return CUBIT_FAILURE;
  }
  // Make sure we were successful
//  if ( sheet_body->lump() == NULL ||
//    sheet_body->lump()->shell() == NULL ||
//    sheet_body->lump()->shell()->first_face() == NULL )
//  {
//    PRINT_ERROR("Couldn't build a body from new offset surface.\n" );
//    api_delent( new_FACE_ptr );
//    return CUBIT_FAILURE;
//  }

  // Make a new body
  new_body = AcisQueryEngine::instance()->populate_topology_bridges(sheet_body);
  return CUBIT_SUCCESS;
}

CubitStatus
AcisSurfaceTool::create_offset_body( BodySM* body_ptr, 
                                     BodySM*& new_body, 
                                     double offset_distance )
{
#ifndef ACIS_LOCAL_OPS
      PRINT_ERROR( "The ACIS Local Operations Husk is required for offsetting\n"
            "       bodies.  It has not been licensed for this installation.\n" );
      return CUBIT_FAILURE;
#endif

  outcome result;

  // we need an AcisModifyEngine to exist; make sure it does
  assert(AcisModifyEngine::instance() != 0);

  // Temporary workaround (remove after create_offset_body is moved
  // into AcisTweakTool).  Need to call api_initialize_local_ops.
  /*AcisTweakTool* acis_tweak =*/ AcisTweakTool::instance();

  // Get the BODY to offset
  
  // Get the BodyACIS part of the OSMEPtr
  BodyACIS* bodyACISPtr = CAST_TO(body_ptr, BodyACIS) ;
  
  if ( bodyACISPtr == NULL )
  {
    PRINT_ERROR("Body is not an ACIS Body\n");
    return CUBIT_FAILURE;
  }
  
  // If we have a valid BodyACIS, proceed.
  
  // Get the BODY pointer
  BODY* BODY_ptr = bodyACISPtr->get_BODY_ptr() ;
  
  // Make sure there is a valid associated ACIS BODY.
  if (BODY_ptr == NULL)
  {
    PRINT_ERROR("Body has no associated ACIS Body\n");
    return CUBIT_FAILURE;
  }
  
  // Copy the body.
  BODY* copied_BODY_ptr = AcisModifyEngine::instance()->copy_BODY(BODY_ptr, CUBIT_TRUE);
  
  // Now, offset the entire copied body.
  SPAposition box_l(0,0,0);
  SPAposition box_h(0,0,0);
#ifdef ACIS_LOCAL_OPS
  result = api_offset_body( copied_BODY_ptr, offset_distance, box_l, box_h );
#endif
  if( !result.ok() )
  {
    AcisQueryEngine::instance()->ACIS_API_error(result);
    PRINT_ERROR( "Error offsetting body\n");
    api_delent(copied_BODY_ptr);
    return CUBIT_FAILURE;
  }
  
  // Make a new body
  new_body = AcisQueryEngine::instance()->populate_topology_bridges(copied_BODY_ptr);
  return CUBIT_SUCCESS;
}

CubitStatus
AcisSurfaceTool::create_skin_surface( DLIList<Curve*>& curves,
                                      BodySM*& new_body )
{
  new_body = NULL;
  int i;
  EDGE* copied_EDGE_ptr;
  outcome result;

    // check the first against the last to see if we should be closed
  curves.reset();
  bool closed = false;
  if (curves.get() == curves.prev()) {
    closed = true;
    curves.pop();
  }
  
  // Create the wires from the edges 
  BODY** wires = new BODY*[curves.size()];
  
  Curve* curve_ptr;
  EDGE* EDGE_ptr;
  curves.reset();
  int num_wires = 0;
  for( i=0; i<curves.size(); i++ )
  {
    curve_ptr = curves.get_and_step();
    CurveACIS* acis_curve_ptr = dynamic_cast<CurveACIS*>(curve_ptr);
    
    // Get the EDGE ptr from the RefEdge
    if( !acis_curve_ptr || !(EDGE_ptr = acis_curve_ptr->get_EDGE_ptr()) )
    {
      PRINT_ERROR( "Curve is not an ACIS curve\n" );
      for( i=0; i<num_wires; i++ )
        api_delent( wires[i] );
      delete [] wires;
      return CUBIT_FAILURE;
    }
    
    result = api_edge( EDGE_ptr, copied_EDGE_ptr );
    if (!result.ok())
    {
      for( i=0; i<num_wires; i++ )
        api_delent( wires[i] );
      delete [] wires;
      PRINT_ERROR( "Unable to copy curve -- aborting\n" );
      AcisQueryEngine::instance()->ACIS_API_error (result);
      return CUBIT_FAILURE;
    }
    
    result = api_make_ewire( 1, &copied_EDGE_ptr, wires[i]); 
    
    if( !result.ok() )
    {
      AcisQueryEngine::instance()->ACIS_API_error(result);
      PRINT_ERROR( "unable make ACIS WIRE BODY from curve -- aborting\n" );
      for( i=0; i<num_wires; i++ )
        api_delent( wires[i] );
      delete [] wires;
      return CUBIT_FAILURE;
    }
    
    // Make sure the wire doesn't self-intersect
    result = api_check_wire_self_inters(wires[i]);
    
    if( !result.ok() )
    {
      AcisQueryEngine::instance()->ACIS_API_error(result);
      PRINT_ERROR( "curve self-intersects, aborting\n" );
      for( i=0; i<num_wires; i++ )
        api_delent( wires[i] );
      delete [] wires;
      return CUBIT_FAILURE;
    }
    num_wires++;
  }

  BODY *skin_body = NULL; 
  
  PRINT_INFO( "Creating skin surface...\n" );
  skin_options tmp_skin_options;
  tmp_skin_options.set_solid(0);
  if (closed) {
    tmp_skin_options.set_closed(2);
    tmp_skin_options.set_periodic(TRUE);
    tmp_skin_options.set_solid(1);
  }
  
  result = api_skin_wires( num_wires, wires, skin_body, &tmp_skin_options ); 

  // Cleanup memory
  for( i=0; i<num_wires; i++ )
    api_delent( wires[i] );
  delete [] wires;
  
  if( !result.ok() )
  {
    AcisQueryEngine::instance()->ACIS_API_error(result);
    PRINT_ERROR( "unable to create skin surface\n" );
    return CUBIT_FAILURE;
  }

#ifdef ACIS_HEALER
  // Now heal the skin body (found problems if we don't)
  if( AcisHealerTool::instance()->init_BODY_for_healing( skin_body ) == CUBIT_SUCCESS )
  {
     int percent_before, percent_after, number_splines_simplified;
     if( AcisHealerTool::instance()->heal_BODY( skin_body, percent_before,
        percent_after, number_splines_simplified ) == CUBIT_FAILURE )
        PRINT_ERROR( "Error healing the skin body\n" );
     AcisHealerTool::instance()->end_BODY_for_healing( skin_body );

     // Okay, also copy the body (strange but true...we need to copy the body too or
     // there may be problems).  Hmmm...not a good situation but for now just do it
     // (SRS - 9-21-2000)
     BODY* new_body_ptr = NULL;
     
     outcome result = api_copy_body(skin_body, new_body_ptr);
     if (!result.ok())
     {
        AcisQueryEngine::instance()->ACIS_API_error(result);
        new_body_ptr = NULL;
     }
     else
     {
        result = api_delent( skin_body );
        if (!result.ok())
        {
           PRINT_ERROR("Problems deleting an intermediate skin body..this may use up computer memory.\n");
           AcisQueryEngine::instance()->ACIS_API_error(result);
        }
        skin_body = new_body_ptr;
     }
  }
  else
     PRINT_WARNING( "Unable to heal the skin body\n" );
#endif
  
  // Make a new body
  new_body = AcisQueryEngine::instance()->populate_topology_bridges(skin_body);
  return CUBIT_SUCCESS;
}





CubitStatus
AcisSurfaceTool::loft_surfaces( Surface *face1, const double &takeoff1,
                                Surface *face2, const double &takeoff2,
                                BodySM*& new_body,
                                CubitBoolean arc_length_option,
                                CubitBoolean twist_option,
                                CubitBoolean align_direction,
                                CubitBoolean perpendicular,
                                CubitBoolean simplify_option)
{

  new_body = NULL;

  SurfaceACIS* surf1;
  SurfaceACIS* surf2;
  FACE* FACE_ptr1;
  FACE* FACE_ptr2;
  ENTITY* ENTITY_ptr;
  outcome result;
  
  surf1 = dynamic_cast<SurfaceACIS*>(face1);
  surf2 = dynamic_cast<SurfaceACIS*>(face2);
  if (!surf1 || !surf2)
  {
    PRINT_ERROR( "Surface is not an ACIS surface -- aborting\n" );
    return CUBIT_FAILURE;
  }
  
  FACE_ptr1 = surf1->get_FACE_ptr();
  FACE_ptr2 = surf2->get_FACE_ptr();
  if (!FACE_ptr1 || !FACE_ptr2)
  {
    PRINT_ERROR( "Invalid ACIS Surface -- aborting.\n");
    return CUBIT_FAILURE;
  }
  
  api_copy_entity_contents( (ENTITY *)FACE_ptr1, ENTITY_ptr);
  FACE_ptr1 = (FACE *)ENTITY_ptr;

  api_copy_entity_contents( (ENTITY *)FACE_ptr2, ENTITY_ptr);
  FACE_ptr2 = (FACE *)ENTITY_ptr;


  BODY *loft_body = NULL; 

    // build the two lists of coedges
  Loft_Connected_Coedge_List loft_struc[3];
  
  int i;
  FACE *temp_FACE = FACE_ptr1;
  for (i = 0; i <= 1; i++) {
    DLIList<COEDGE*> co_edges;
    LOOP *temp_LOOP = temp_FACE->loop();
    while (temp_LOOP) {
      
      COEDGE* temp_COEDGE = temp_LOOP->start();
      do {
        co_edges.append(temp_COEDGE);
        temp_COEDGE = temp_COEDGE->next();
      }
      while (temp_COEDGE && temp_COEDGE != temp_LOOP->start());
      
      temp_LOOP = temp_LOOP->next();
    }
    
    loft_struc[i].coedge_list = new COEDGE*[co_edges.size()];
    int j;
    co_edges.reset();
    for (j = 0; j < co_edges.size(); j++) 
      loft_struc[i].coedge_list[j] = co_edges.get_and_step();
    
    loft_struc[i].n_list = co_edges.size();
    if (i == 1) loft_struc[i].coedge_orient = REVERSED;
    else loft_struc[i].coedge_orient = FORWARD;
    
    loft_struc[i].cross_tg_attr = (i == 0 ? takeoff1 : takeoff2);
    loft_struc[i].law_list = NULL;
  }

  loft_struc[2] = loft_struc[1];
    
  PRINT_INFO( "Creating loft surface...\n" );
  CubitBoolean closed = twist_option;
  CubitBoolean solid = perpendicular;
  result = api_loft_coedges( 2, loft_struc, 
                             loft_body, 
                             arc_length_option, CUBIT_FALSE,
                             align_direction, CUBIT_FALSE,
                             simplify_option,
                             closed, solid); 

  delete [] loft_struc[0].coedge_list;
  delete [] loft_struc[1].coedge_list;
  
  if( !result.ok() )
  {
    AcisQueryEngine::instance()->ACIS_API_error(result);
    PRINT_ERROR( "unable to create loft surface\n" );
    return CUBIT_FAILURE;
  }
  
  // Make a new body
  new_body = AcisQueryEngine::instance()->populate_topology_bridges(loft_body);
  return CUBIT_SUCCESS;
}




CubitStatus
AcisSurfaceTool::create_surface( DLIList<CubitVector*>& vec_list,
                                 BodySM *&new_body,
                                 Surface *surface_ptr,
                                 CubitBoolean project_points )

{

  int i;
  if (surface_ptr)
  {
    // Check the project_points option and do the necessary checks or projections.
    if (project_points)
    {
      // Create a new list of points that are projected to the surface
      vec_list.reset();
      DLIList<CubitVector*> new_vec_list;
      CubitVector *vec_ptr, new_vec;
      for( i=0; i<vec_list.size(); i++ )
      {
        vec_ptr = vec_list.get_and_step();
        surface_ptr->closest_point( *vec_ptr, &new_vec );
        new_vec_list.append( new CubitVector(new_vec) );
      }

      // Make the surface in the solid modeller
      // Currently only allowed in ACIS
      CubitStatus status = create_surface( new_vec_list, new_body, surface_ptr );
      for( i=0; i<new_vec_list.size(); i++ )
      {
        vec_ptr = new_vec_list.get_and_step();
        delete vec_ptr;
      }
      return status;
    }
    
    else
    {
      // Make sure the points lie on the surface
      //CubitVector new_vec;
      vec_list.reset();
      CubitVector loc_on_surf;
      CubitVector *vec_ptr;
      for( i=0; i<vec_list.size(); i++ )
      {
        vec_ptr = vec_list.get_and_step();
        surface_ptr->closest_point( *vec_ptr, &loc_on_surf );
   
        if (!vec_ptr->within_tolerance(loc_on_surf, GEOMETRY_RESABS))
        {
          PRINT_ERROR("all locations must lie on Surface\n" );
          return CUBIT_FAILURE;
        }
      }
      // Make the surface in the solid modeller
      create_surface( vec_list, new_body, surface_ptr );
    }

  }
  return create_surface( vec_list, new_body );//In this case, the ref_face_ptr is NULL
}


CubitStatus 
AcisSurfaceTool::create_surface( DLIList<CubitVector*>& vec_list,
                                 BodySM *&new_body,
                                 Surface *surface_ptr )
{
  int i, j;
  outcome result;
  FACE *FACE_ptr = NULL;
  SurfaceACIS* acis_surf_ptr = dynamic_cast<SurfaceACIS*>(surface_ptr);

  if( surface_ptr )
  {
    if( !acis_surf_ptr || !(FACE_ptr = acis_surf_ptr->get_FACE_ptr()) )
    {
      PRINT_ERROR( "Surface is not an ACIS surface -- aborting\n" );
      return CUBIT_FAILURE;
    }
  }

  // Create curves connecting the points
  DLIList<EDGE*> EDGE_list;
  EDGE *EDGE_ptr;
  if( acis_surf_ptr )
  {
    vec_list.reset();
    CubitVector *start = vec_list.get_and_step();
    CubitVector *end;
    for( i=0; i<vec_list.size(); i++ )
    {
      end = vec_list.get_and_step();

      if( start->within_tolerance( *end, GEOMETRY_RESABS) )
      {
         PRINT_ERROR( "Attempt to create a line between coincident points at %f, %f, %f\n",
            start->x(), start->y(), start->z() );
         for( j=0; j<EDGE_list.size(); j++ )
        {
          EDGE_ptr = EDGE_list.get_and_step();
          api_delent( EDGE_ptr );
        }
        return CUBIT_FAILURE;
      }
      
      CubitVector surf_norm;
      acis_surf_ptr->closest_point( *start, 0, &surf_norm );
      
      CubitVector plane_norm;
      
      CubitVector edge0 = *end - *start;
      CubitVector edge1 = surf_norm;
      plane_norm = edge0*edge1;
      if ( plane_norm.length() < CUBIT_RESABS )
      {
        //try using the normal at point 2.
        acis_surf_ptr->closest_point( *end, 0, &surf_norm );
        edge1 = surf_norm;
        edge0 = *start - *end;
        plane_norm = edge0*edge1;
        if ( plane_norm.length() < CUBIT_RESABS )
        {
          PRINT_ERROR( "unable to create line on surface from (%f, %f, %f) to (%f, %f, %f)\n",
             start->x(), start->y(), start->z(), end->x(), end->y(), end->z() );
          for( j=0; j<EDGE_list.size(); j++ )
          {
             EDGE_ptr = EDGE_list.get_and_step();
             api_delent( EDGE_ptr );
          }
          return CUBIT_FAILURE;
        }
      }
      plane_norm.normalize();
      
      // Create vertices at the ends
      VERTEX *start_VERTEX_ptr = AcisModifyEngine::instance()->make_VERTEX( *start );
      VERTEX *end_VERTEX_ptr = AcisModifyEngine::instance()->make_VERTEX( *end );

      if( !start_VERTEX_ptr || !end_VERTEX_ptr )
      {
        if (start_VERTEX_ptr)
          api_delent( start_VERTEX_ptr );
        for( j=0; j<EDGE_list.size(); j++ )
        {
          EDGE_ptr = EDGE_list.get_and_step();
          api_delent( EDGE_ptr );
        }
        return CUBIT_FAILURE;
      }

      EDGE_ptr = AcisModifyEngine::instance()->make_surface_EDGE( start_VERTEX_ptr, end_VERTEX_ptr,
        FACE_ptr, plane_norm );

      if( !EDGE_ptr )
      {
        api_delent( (ENTITY *)start_VERTEX_ptr );
        api_delent( (ENTITY *)end_VERTEX_ptr );
        for( j=0; j<EDGE_list.size(); j++ )
        {
          EDGE_ptr = EDGE_list.get_and_step();
          api_delent( EDGE_ptr );
        }
        return CUBIT_FAILURE;
      }
      else
        EDGE_list.append( EDGE_ptr );

      start = end;
    }
  }
  else
  {
    vec_list.reset();
    CubitVector *coord1 = vec_list.get_and_step();
    CubitVector *coord2 = NULL;
    SPAposition start( coord1->x(), coord1->y(), coord1->z() );
    SPAposition end;
    for( i=0; i<vec_list.size(); i++ )
    {
      coord2 = vec_list.get_and_step();
      end.set_x( coord2->x() );
      end.set_y( coord2->y() );
      end.set_z( coord2->z() );

      if( coord1->within_tolerance( *coord2, GEOMETRY_RESABS ) )
      {
         PRINT_ERROR( "Attempt to create a line between coincident points at (%f, %f, %f)\n",
            start.x(), start.y(), start.z() );
         for( j=0; j<EDGE_list.size(); j++ )
        {
          EDGE_ptr = EDGE_list.get_and_step();
          api_delent( EDGE_ptr );
        }
        return CUBIT_FAILURE;
      }

      result = api_mk_ed_line ( start, end, EDGE_ptr );
      if( !result.ok() )
      {
        AcisQueryEngine::instance()->ACIS_API_error(result);
        PRINT_ERROR( "unable to create line from (%f, %f, %f) to (%f, %f, %f)\n",
          start.x(), start.y(), start.z(), end.x(), end.y(), end.z() );
        for( j=0; j<EDGE_list.size(); j++ )
        {
          EDGE_ptr = EDGE_list.get_and_step();
          api_delent( EDGE_ptr );
        }
        return CUBIT_FAILURE;
      }
      EDGE_list.append( EDGE_ptr );

      start = end;
      coord1 = coord2;
    }
  }

  // Make a WIRE from the EDGES
  int edge_count = vec_list.size();
  EDGE** EDGEs = new EDGE* [edge_count];
  EDGE_list.reset();
  for( i=0; i<edge_count; i++ )
    EDGEs[i] = EDGE_list.get_and_step();

  BODY* wire_BODY = NULL;
  result = api_make_ewire(edge_count, EDGEs, wire_BODY) ;
  delete [] EDGEs;
  if ( !result.ok() || wire_BODY == NULL )
  {
    AcisQueryEngine::instance()->ACIS_API_error( result );
    if (wire_BODY)
      api_delent(wire_BODY);
    else
    {
      for( j=0; j<EDGE_list.size(); j++ )
      {
        EDGE_ptr = EDGE_list.get_and_step();
        api_delent( EDGE_ptr );
      }
    }
    return CUBIT_FAILURE;
  }
  
    // Check to make sure the wire is closed
  result = api_closed_wire( wire_BODY );
  if( !result.ok() )
  {
    AcisQueryEngine::instance()->ACIS_API_error( result );
    PRINT_ERROR( "ACIS reports that a closed loop was not found\n" );
    api_delent( wire_BODY ); // This deletes the underlying EDGEs
    return CUBIT_FAILURE;   
  }
  
    // Use the WIRE to make a FACE.
    // Note that the call to api_cover_wires creates not
    // only the FACE but the entire topological data
    // structure above FACE (till BODY).
  ENTITY_LIST acis_FACE_list;
  if ( FACE_ptr != NULL )
  {
    result = api_cover_wires( wire_BODY,
                              FACE_ptr->geometry()->equation(),
                              acis_FACE_list );
    
      // Now check to see if we covered the wire successfully
    if ( !result.ok() || acis_FACE_list.count() == 0 )
    {
      AcisQueryEngine::instance()->ACIS_API_error( result );
      PRINT_ERROR("Error in constructing surface on top of another"
                  "surface, the edges must fit the surface...\n");
      api_delent( wire_BODY );
      return CUBIT_FAILURE;
    }
  }
  else
  {
    result = api_cover_wires ( wire_BODY,
                               *(surface*)NULL_REF,
                               acis_FACE_list );
    
    if ( !result.ok() || acis_FACE_list.count() == 0 )
    {
      AcisQueryEngine::instance()->ACIS_API_error( result );
      api_delent( wire_BODY );
      return CUBIT_FAILURE;
    }
  }
  if (acis_FACE_list.count() != 1)
  {
    PRINT_ERROR("ACIS api_cover_wires function created more than one FACE.\n"
                "       At this time, we cannot deal with this.\n");
    api_delent( wire_BODY );
    return CUBIT_FAILURE;
  }
  
    //make this a 2d body.
  result = api_body_to_2d( wire_BODY );
  if ( !result.ok() )
  {
    AcisQueryEngine::instance()->ACIS_API_error( result );
    PRINT_ERROR("ACIS api_body_to_2d function failed.\n");
    api_delent( wire_BODY );
    return CUBIT_FAILURE;
  }
  FACE_ptr = (FACE*) acis_FACE_list.next();
  if ( FACE_ptr == NULL || FACE_ptr->geometry() == NULL )
  {
    PRINT_ERROR("a surface could not be constructed from the resultant polygon.\n");
    api_delent( wire_BODY );
    return CUBIT_FAILURE;
  }

  Surface *new_surface_ptr = AcisQueryEngine::instance()->populate_topology_bridges( FACE_ptr );
  new_body = AcisModifyEngine::instance()->make_BodySM( new_surface_ptr );
 
  return CUBIT_SUCCESS;
}

static SPAtransf const &get_edge_trans(EDGE *edge)
{
	// Find the body

	ENTITY *entity=NULL;
	if( edge->coedge()==NULL)
	{
		return *(SPAtransf *)NULL_REF;
	}
	else if(edge->coedge()->loop()!=NULL)
	{
		entity=edge->coedge()->loop()->face()->shell()->lump()->body();
	}
	else if(edge->coedge()->wire()!=NULL)
	{
		entity=edge->coedge()->wire()->body();
		if(!entity)
		{
			entity=edge->coedge()->wire()->shell()->lump()->body();
		}
	}
	else
	{
		return *(SPAtransf *)NULL_REF;
	}

	// Get the transform
	
	if( ((BODY *)entity)->transform() != NULL )
	{
		return ((BODY *)entity)->transform()->transform();
	}
	else
	{
		return *(SPAtransf *)NULL_REF;
	}
}

CubitStatus 
AcisSurfaceTool::create_weld_surface( CubitVector &root, 
                                      Surface *surface1_in, double leg1, 
                                      Surface *surface2_in, double leg2, 
                                      BodySM *&new_body )
{
  CubitVector root1, root2, dir1, dir2, pnt1, pnt2;
  
  SurfaceACIS* surface1 = dynamic_cast<SurfaceACIS*>(surface1_in);
  SurfaceACIS* surface2 = dynamic_cast<SurfaceACIS*>(surface2_in);
  if (!surface1 || !surface2)
  {
    PRINT_ERROR("Surface is not an ACIS Surface.\n");
    return CUBIT_FAILURE;
  }
  
  // First make sure root point is shared between both surfaces
  surface1->closest_point( root, &root1, &dir1 );
  surface2->closest_point( root, &root2, &dir2 );
  if( root1.within_tolerance( root2, GEOMETRY_RESABS) == CUBIT_FALSE )
  {
    PRINT_ERROR( "Root point must lie at common location on both Surfaces for weld creation\n" );         
    return CUBIT_FAILURE;
  }

  root.next_point( dir1, leg1, pnt1 );
  root.next_point( dir2, leg2, pnt2 );

  surface2->closest_point( CubitVector(pnt1), &pnt1 );
  surface1->closest_point( CubitVector(pnt2), &pnt2 );

  // Create a curve from root to pnt1, on ref_face_ptr1
  EDGE *EDGE_ptr1 = NULL;
  if( create_EDGE_on_Surface( root, pnt1, surface2, EDGE_ptr1 ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;
  double EDGE_ptr1_length = EDGE_ptr1->length();
  if( fabs(EDGE_ptr1_length - leg1) > GEOMETRY_RESABS )
     extend_weld_EDGE( EDGE_ptr1, leg1 );

  // Create a curve from root to pnt1, on ref_face_ptr2
  EDGE *EDGE_ptr2 = NULL;
  if( create_EDGE_on_Surface( root, pnt2, surface1, EDGE_ptr2 ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;
  double EDGE_ptr2_length = EDGE_ptr2->length();
  if( fabs(EDGE_ptr2_length - leg2) > GEOMETRY_RESABS )
     extend_weld_EDGE( EDGE_ptr2, leg2 );

  // Create a curve from the end of EDGE_ptr1 to the end of EDGE_ptr2
  EDGE *EDGE_ptr3;
  SPAposition start = EDGE_ptr1->end_pos();
  SPAposition end = EDGE_ptr2->end_pos();
  outcome result = api_mk_ed_line( start, end, EDGE_ptr3 );
  if ( !result.ok() )
  {
    AcisQueryEngine::instance()->ACIS_API_error( result );
    PRINT_ERROR("ACIS api_mk_ed_line function failed.\n");
    api_delent( EDGE_ptr1 );
    api_delent( EDGE_ptr2 );
    return CUBIT_FAILURE;
  }

  // Create a surface from these 3 curves

  // Make a WIRE from the EDGES
  EDGE** EDGEs = new EDGE* [3];
  EDGEs[0] = EDGE_ptr1;
  EDGEs[1] = EDGE_ptr2;
  EDGEs[2] = EDGE_ptr3;

  BODY* wire_BODY = NULL;
  result = api_make_ewire( 3, EDGEs, wire_BODY );
  delete [] EDGEs;
  if ( !result.ok() || wire_BODY == NULL )
  {
    AcisQueryEngine::instance()->ACIS_API_error( result );
    if (wire_BODY)
      api_delent(wire_BODY);
    else
    {
       api_delent( EDGE_ptr1 );
       api_delent( EDGE_ptr2 );
       api_delent( EDGE_ptr3 );
    }
    return CUBIT_FAILURE;
  }
  
  // Check to make sure the wire is closed
  result = api_closed_wire( wire_BODY );
  if( !result.ok() )
  {
    AcisQueryEngine::instance()->ACIS_API_error( result );
    PRINT_ERROR( "ACIS reports that a closed loop was not found\n" );
    api_delent( wire_BODY ); // This deletes the underlying EDGEs
    return CUBIT_FAILURE;   
  }
  
  // Use the WIRE to make a FACE.
  // Note that the call to api_cover_wires creates not
  // only the FACE but the entire topological data
  // structure above FACE (till BODY).
  ENTITY_LIST acis_FACE_list;
  result = api_cover_wires ( wire_BODY,
    *(surface*)NULL_REF,
    acis_FACE_list );
  
  if ( !result.ok() || acis_FACE_list.count() == 0 )
  {
    AcisQueryEngine::instance()->ACIS_API_error( result );
    api_delent( wire_BODY );
    return CUBIT_FAILURE;
  }
  if (acis_FACE_list.count() != 1)
  {
    PRINT_ERROR("ACIS api_cover_wires function created more than one FACE.\n"
                "       At this time, we cannot deal with this.\n");
    api_delent( wire_BODY );
    return CUBIT_FAILURE;
  }
  
  //make this a 2d body.
  result = api_body_to_2d( wire_BODY );
  if ( !result.ok() )
  {
    AcisQueryEngine::instance()->ACIS_API_error( result );
    PRINT_ERROR("ACIS api_body_to_2d function failed.\n");
    api_delent( wire_BODY );
    return CUBIT_FAILURE;
  }
  FACE *FACE_ptr = (FACE*) acis_FACE_list.next();
  if ( FACE_ptr == NULL || FACE_ptr->geometry() == NULL )
  {
    PRINT_ERROR("a surface could not be constructed from the resultant polygon.\n");
    api_delent( wire_BODY );
    return CUBIT_FAILURE;
  }

  Surface *new_surface_ptr = AcisQueryEngine::instance()->populate_topology_bridges( FACE_ptr );
  new_body = AcisModifyEngine::instance()->make_BodySM( new_surface_ptr );
  return CUBIT_SUCCESS;
}

CubitStatus
AcisSurfaceTool::create_EDGE_on_Surface( const CubitVector &start_in, 
                                         const CubitVector &end_in,
                                         SurfaceACIS* surface_ptr, 
                                         EDGE *&EDGE_ptr )
{
  FACE *FACE_ptr = surface_ptr->get_FACE_ptr();
  if( FACE_ptr == NULL )
  {
    PRINT_ERROR( "Surface is not an ACIS surface -- aborting\n" );
    return CUBIT_FAILURE;
  }

  CubitVector start, end, surf_norm;
  surface_ptr->closest_point( start_in, &start, &surf_norm );
  surface_ptr->closest_point( end_in, &end );
  
  if( start.within_tolerance( end, GEOMETRY_RESABS) )
  {
    PRINT_ERROR( "Attempt to create a curve on Surface between coincident points at (%f, %f, %f)\n",
      start.x(), start.y(), start.z() );         
    return CUBIT_FAILURE;
  }
  
  CubitVector plane_norm;
  
  CubitVector edge0 = end - start;
  CubitVector edge1 = surf_norm;
  plane_norm = edge0*edge1;
  if ( plane_norm.length() < CUBIT_RESABS )
  {
    //try using the normal at point 2.
    surface_ptr->closest_point( end, 0, &surf_norm );
    edge1 = surf_norm;
    edge0 = start - end;
    plane_norm = edge0*edge1;
    if ( plane_norm.length() < CUBIT_RESABS )
    {
      PRINT_ERROR( "unable to create line on Surface from (%f, %f, %f) to (%f, %f, %f)\n",
        start.x(), start.y(), start.z(), end.x(), end.y(), end.z() );
      return CUBIT_FAILURE;
    }
  }
  plane_norm.normalize();
  
  // Create vertices at the ends
  VERTEX *start_VERTEX_ptr = AcisModifyEngine::instance()->make_VERTEX( start );
  VERTEX *end_VERTEX_ptr = AcisModifyEngine::instance()->make_VERTEX( end );
  
  if( !start_VERTEX_ptr || !end_VERTEX_ptr )
  {
    if (start_VERTEX_ptr)
      api_delent(start_VERTEX_ptr);
    return CUBIT_FAILURE;
  }
  
  EDGE_ptr = AcisModifyEngine::instance()->make_surface_EDGE( start_VERTEX_ptr, end_VERTEX_ptr,
    FACE_ptr, plane_norm );
  
  if( !EDGE_ptr )
  {
    api_delent( (ENTITY *)start_VERTEX_ptr );
    api_delent( (ENTITY *)end_VERTEX_ptr );
    return CUBIT_FAILURE;
  }

  return CUBIT_SUCCESS;
}

void
AcisSurfaceTool::extend_weld_EDGE( EDGE *&EDGE_ptr, double desired_length )
{
  double length = EDGE_ptr->length();
  
  // Assume start is at root of weld
  double start = EDGE_ptr->start_param();
  double end = EDGE_ptr->end_param();
  
  // Assume parameter space is well behaved
  double ratio = desired_length/length;
  
  double new_end = start + ratio*(end-start);
  
  SPAtransf const &etrans=get_edge_trans( EDGE_ptr );
  curve *cur = EDGE_ptr->geometry()->trans_curve( etrans );
  SPAinterval newint( start, new_end );
  //SPAinterval test1=
    extend_curve( *cur, newint );
  
  // get ends of the edge
  SPAposition pt0=cur->eval_position( start );
  SPAposition pt1=cur->eval_position( new_end );
  
  // create the vertices
  VERTEX *svert= ACIS_NEW VERTEX( ACIS_NEW APOINT(pt0));
  VERTEX *evert= ACIS_NEW VERTEX( ACIS_NEW APOINT(pt1));	
  
  // create the curve
  CURVE *the_curve=make_curve(*cur);
  ACIS_DELETE cur;
  
  // make the edge
  EDGE *new_edge= ACIS_NEW EDGE(svert,evert,the_curve,FORWARD);	
  
  // set the param range of the edge
  SPAinterval range(start,new_end);
  new_edge->set_param_range(range);
  
  api_delent( EDGE_ptr );
  EDGE_ptr = new_edge;
}

CubitStatus
AcisSurfaceTool::loft_surfaces_to_body( 
                       Surface *face1, const double &takeoff1,
                       Surface *face2, const double &takeoff2,
                       BodySM*& new_body,
                       CubitBoolean arc_length_option,
                       CubitBoolean twist_option,
                       CubitBoolean align_direction,
                       CubitBoolean perpendicular,
                       CubitBoolean simplify_option)
{
   
   new_body = NULL;

   if (face1 == face2)
   {
      PRINT_ERROR("Can't loft surface to itself.\n");
	  return CUBIT_FAILURE;
   }

  SurfaceACIS* surf1 = dynamic_cast<SurfaceACIS*>(face1);
  SurfaceACIS* surf2 = dynamic_cast<SurfaceACIS*>(face2);
  FACE* FACE_ptr1;
  FACE* FACE_ptr2;
  ENTITY* ENTITY_ptr;
  outcome result;
  
  if( !surf1 || !surf2 )
  {
    PRINT_ERROR( "Surface is not an ACIS surface -- aborting\n" );
    return CUBIT_FAILURE;
  }
  
  FACE_ptr1 = surf1->get_FACE_ptr();
  FACE_ptr2 = surf2->get_FACE_ptr();
  if (!FACE_ptr1 || !FACE_ptr2)
  {
    PRINT_ERROR("Invalid ACIS Surface -- aborting.\n");
    return CUBIT_FAILURE;
  }
  
  api_copy_entity_contents( (ENTITY *)FACE_ptr1, ENTITY_ptr);
  FACE_ptr1 = (FACE *)ENTITY_ptr;

  api_copy_entity_contents( (ENTITY *)FACE_ptr2, ENTITY_ptr);
  FACE_ptr2 = (FACE *)ENTITY_ptr;

  //make both FACEs single sided before loft...bug in ACIS.
  //may be able to take this out later on.
  //FACE_ptr1->set_sides( SINGLE_SIDED );
  //FACE_ptr2->set_sides( SINGLE_SIDED );

  BODY *loft_body = NULL; 

    
  PRINT_INFO( "Creating new body from lofted surfaces...\n" );
  //CubitBoolean closed = twist_option;
  //CubitBoolean solid = perpendicular;
  
   result = api_loft_faces(FACE_ptr1, 
								takeoff1, 
								FACE_ptr2, 
								takeoff2, 
								loft_body, 
								arc_length_option, 
								twist_option, 
								align_direction,  
								perpendicular,
								simplify_option);
  
  if( !result.ok() )
  {
    AcisQueryEngine::instance()->ACIS_API_error(result);
    if (loft_body)
    {
      api_delent(loft_body);
    }
    else
    {
      api_delent(FACE_ptr1);
      api_delent(FACE_ptr2);
    }
    PRINT_ERROR( "Unable to create loft body\n" );
    return CUBIT_FAILURE;
  }
  
  // Make a new body
  new_body = AcisQueryEngine::instance()->populate_topology_bridges(loft_body);
  return CUBIT_SUCCESS;

}

TopologyBridge* AcisSurfaceTool::get_acis_bridge( TopologyEntity* topo_ptr )
  { return topo_ptr->bridge_manager()->topology_bridge( AcisQueryEngine::instance() ); }
SurfaceACIS* AcisSurfaceTool::get_acis_surface( RefFace* face_ptr )
  { return dynamic_cast<SurfaceACIS*>(get_acis_bridge(face_ptr)); }
CurveACIS* AcisSurfaceTool::get_acis_curve( RefEdge* edge_ptr )
  { return dynamic_cast<CurveACIS*>(get_acis_bridge(edge_ptr)); }
BodyACIS* AcisSurfaceTool::get_acis_body( Body* body_ptr )
  { return dynamic_cast<BodyACIS*>(get_acis_bridge(body_ptr)); }

CubitStatus
AcisSurfaceTool::create_net_surface( DLIList<RefEdge*>& u_curves, 
                                     DLIList<RefEdge*>& v_curves, 
                                     Body*& new_body,
                                     double net_tol, CubitBoolean heal,
                                     CubitBoolean verbose)
{
  DLIList<Curve*> u_list(u_curves.size());
  DLIList<Curve*> v_list(v_curves.size());
  BodySM* this_bodysm;
  CubitStatus result;
  int i;
  
  new_body = 0;
  
  u_curves.reset();
  for (i = u_curves.size(); i--; )
  {
    RefEdge* edge_ptr = u_curves.get_and_step();
    CurveACIS* curve_ptr = get_acis_curve(edge_ptr);
    if (curve_ptr)
      u_list.append(curve_ptr);
    else
      PRINT_ERROR("Curve %d is not an ACIS Curve.\n", edge_ptr->id());
  }
  
  v_curves.reset();
  for (i = v_curves.size(); i--; )
  {
    RefEdge* edge_ptr = v_curves.get_and_step();
    CurveACIS* curve_ptr = get_acis_curve(edge_ptr);
    if (curve_ptr)
      v_list.append(curve_ptr);
    else
      PRINT_ERROR("Curve %d is not an ACIS Curve.\n", edge_ptr->id());
  }
  
  if (u_list.size() != u_curves.size() || v_list.size() != v_curves.size())
    return CUBIT_FAILURE;
  
  result = create_net_surface( u_list, v_list, this_bodysm, net_tol, heal, verbose );
  if (!result)
    return result;
  
  // Make a new body
  new_body = GeometryQueryTool::instance()->make_Body(this_bodysm);
  DLIList<RefFace*> new_face_list;
  new_body->ref_faces( new_face_list );
  
  if (verbose)
     PRINT_INFO( "Created new \"net\" surface %d (in body %d)\n", new_face_list.get()->id(),
    new_body->id() );
  
  return CUBIT_SUCCESS;
}


CubitStatus
AcisSurfaceTool::create_net_surface( Body*& new_body, 
                                     DLIList<DLIList<CubitVector*>*> &vec_lists_u, 
                                     DLIList<DLIList<CubitVector*>*> &vec_lists_v, 
                                     double net_tol, CubitBoolean heal,
                                     CubitBoolean verbose)
{
  BodySM* this_bodysm = 0;
  CubitStatus result;
  
  new_body = 0;

  result = create_net_surface( this_bodysm, vec_lists_u, 
                               vec_lists_v, net_tol, heal, verbose );
  if (!result)
    return result;
  
  new_body = GeometryQueryTool::instance()->make_Body(this_bodysm);
  DLIList<RefFace*> new_face_list;
  new_body->ref_faces( new_face_list );
  
  if(verbose)
     PRINT_INFO( "Created new \"net\" surface %d (in body %d)\n", new_face_list.get()->id(),
    new_body->id() );
  
  return CUBIT_SUCCESS;
}

CubitStatus
AcisSurfaceTool::create_offset_surface( RefFace* ref_face_ptr, Body*& new_body, 
                                        double offset_distance )
{
  CubitStatus result;
  BodySM* this_bodysm;
  SurfaceACIS* surf_ptr;
  
  new_body = 0;
  
  surf_ptr = get_acis_surface(ref_face_ptr);
  if (!surf_ptr) 
  {
    PRINT_ERROR("Surface %d is not an ACIS Surface.\n", ref_face_ptr->id());
    return CUBIT_FAILURE;
  }
  
  result = create_offset_surface( surf_ptr, this_bodysm, offset_distance );
  if (!result)
    return result;
    

  // Make a new body
  Body *body_ptr = GeometryQueryTool::instance()->make_Body(this_bodysm);
  DLIList<RefFace*> new_face_list;
  body_ptr->ref_faces( new_face_list );
  
  PRINT_INFO( "Created new \"offset\" surface %d (in body %d)\n", new_face_list.get()->id(),
    body_ptr->id() );
  
  new_body = body_ptr;

  return CUBIT_SUCCESS;
}

CubitStatus
AcisSurfaceTool::create_offset_body( Body* body_ptr, Body*& new_body, 
                                    double offset_distance )
{
  BodySM* this_bodysm;
  CubitStatus result;
  BodyACIS* bodysm_ptr;
  
  new_body = 0;
  
  bodysm_ptr = get_acis_body(body_ptr);
  if (!bodysm_ptr)
  {
    PRINT_ERROR("Body %d is not an ACIS Body.\n", body_ptr->id());
    return CUBIT_FAILURE;
  }
  
  result = create_offset_body( bodysm_ptr, this_bodysm, offset_distance );
  if (!result)
    return CUBIT_FAILURE;
  
  // Make a new body
  new_body = GeometryQueryTool::instance()->make_Body(this_bodysm);
  
  PRINT_INFO( "Created new body %d offset from body %d\n", new_body->id(),
    body_ptr->id() );
  
  return CUBIT_SUCCESS;
}

CubitStatus
AcisSurfaceTool::create_skin_surface( DLIList<RefEdge*>& ref_edge_list,
                                     Body*& new_body )
{
  CubitStatus result;
  BodySM* this_bodysm;
  DLIList<Curve*> curve_list(ref_edge_list.size());
  
  new_body = 0;
  
  ref_edge_list.reset();
  for (int i = ref_edge_list.size(); i--; )
  {
    RefEdge* edge = ref_edge_list.get_and_step();
    CurveACIS* curve = get_acis_curve(edge);
    if (!curve)
      PRINT_ERROR("Curve %d is not an ACIS Curve.\n", edge->id());
    else
      curve_list.append(curve);
  }
  
  if (curve_list.size() != ref_edge_list.size())
    return CUBIT_FAILURE;


  result = create_skin_surface( curve_list, this_bodysm );
  if (!result)
    return result;
  
  // Make a new body
  new_body = GeometryQueryTool::instance()->make_Body(this_bodysm);
  DLIList<RefFace*> new_face_list;
  new_body->ref_faces( new_face_list );
  
  PRINT_INFO( "Created new \"skin\" surface %d (in body %d)\n", new_face_list.get()->id(),
    new_body->id() );
  
  return CUBIT_SUCCESS;
}





CubitStatus
AcisSurfaceTool::loft_surfaces( RefFace *face1, const double &takeoff1,
                                RefFace *face2, const double &takeoff2,
                                Body*& new_body,
                                CubitBoolean arc_length_option,
                                CubitBoolean twist_option,
                                CubitBoolean align_direction,
                                CubitBoolean perpendicular,
                                CubitBoolean simplify_option)
{
  CubitStatus result;
  BodySM* this_bodysm;
  SurfaceACIS* surf1 = get_acis_surface( face1 );
  SurfaceACIS* surf2 = get_acis_surface( face2 );

  new_body = 0;

  if (!surf1)
    PRINT_ERROR("Surface %d is not an ACIS Surface.\n", face1->id());
  if (!surf2)
    PRINT_ERROR("Surface %d is not an ACIS Surface.\n", face2->id());
  if (!surf1 || !surf2)
    return CUBIT_FAILURE;

  result = loft_surfaces( surf1, takeoff1, surf2, takeoff2, this_bodysm,
                          arc_length_option, twist_option, align_direction,
                          perpendicular, simplify_option );
  if (!result)
    return result;
  
  // Make a new body
  new_body = GeometryQueryTool::instance()->make_Body(this_bodysm);
  
  PRINT_INFO( "Created new loft body %d\n", new_body->id() );
  
  return CUBIT_SUCCESS;

}




CubitStatus
AcisSurfaceTool::create_surface( DLIList<CubitVector*>& vec_list,
				 Body *&new_body,
				 RefFace *ref_face_ptr,
				 CubitBoolean project_points )

{
  CubitStatus result;
  BodySM* this_bodysm;
  SurfaceACIS* surface_ptr = 0;
  
  new_body = 0;

  if (ref_face_ptr)
  {
    surface_ptr = get_acis_surface( ref_face_ptr );
    if (!surface_ptr)
    {
      PRINT_ERROR("Surface %d is not an ACIS Surface.\n", ref_face_ptr->id());
      return CUBIT_FAILURE;
    }
  }


  result = create_surface( vec_list, this_bodysm, surface_ptr, project_points );
  if (!result)
    return result;
  
  
  new_body = GeometryQueryTool::instance()->make_Body( this_bodysm );
 
  return CUBIT_SUCCESS;
}


CubitStatus 
AcisSurfaceTool::create_weld_surface( CubitVector &root, RefFace *ref_face_ptr1, 
    double leg1, RefFace *ref_face_ptr2, double leg2, Body *&new_body )
{
  CubitStatus result;
  BodySM* this_bodysm;
  SurfaceACIS* surface_ptr1 = get_acis_surface( ref_face_ptr1 );
  SurfaceACIS* surface_ptr2 = get_acis_surface( ref_face_ptr2 );
  
  new_body = 0;

  if (!surface_ptr1)
    PRINT_ERROR("Surface %d is not an ACIS Surface.\n", ref_face_ptr1->id());
  if (!surface_ptr2)
    PRINT_ERROR("Surface %d is not an ACIS Surface.\n", ref_face_ptr2->id());
  if (!surface_ptr1 || !surface_ptr2)
    return CUBIT_FAILURE;

  result = create_weld_surface( root, surface_ptr1, leg1, surface_ptr2, leg2, this_bodysm );
  if (!result)  
    return result;

  new_body = GeometryQueryTool::instance()->make_Body( this_bodysm );
  DLIList<RefFace*> new_face_list;
  new_body->ref_faces( new_face_list );

  PRINT_INFO( "Created new \"weld\" surface %d (in body %d)\n", new_face_list.get()->id(),
    new_body->id() );
  
  return CUBIT_SUCCESS;
}


CubitStatus
AcisSurfaceTool::loft_surfaces_to_body( RefFace *face1, const double &takeoff1,
                                        RefFace *face2, const double &takeoff2,
                                        Body*& new_body,
                                        CubitBoolean arc_length_option,
                                        CubitBoolean twist_option,
                                        CubitBoolean align_direction,
                                        CubitBoolean perpendicular,
                                        CubitBoolean simplify_option)
{
  CubitStatus result;
  BodySM* this_bodysm;
  SurfaceACIS* surf1 = get_acis_surface( face1 );
  SurfaceACIS* surf2 = get_acis_surface( face2 );

  new_body = 0;

  if (!surf1)
    PRINT_ERROR("Surface %d is not an ACIS Surface.\n", face1->id());
  if (!surf2)
    PRINT_ERROR("Surface %d is not an ACIS Surface.\n", face2->id());
  if (!surf1 || !surf2)
    return CUBIT_FAILURE;

  result = loft_surfaces_to_body( surf1, takeoff1, 
                                  surf2, takeoff2, 
                                  this_bodysm,
                                  arc_length_option, 
                                  twist_option, 
                                  align_direction,
                                  perpendicular, 
                                  simplify_option );
  if (!result)
    return result;
  
  
  // Make a new body
  new_body = GeometryQueryTool::instance()->make_Body(this_bodysm);
  
  PRINT_INFO( "Created new loft body %d\n", new_body->id() );
  
  return CUBIT_SUCCESS;

}
