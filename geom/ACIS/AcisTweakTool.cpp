//- Class:       AcisTweakTool
//- Description: Functions to tweak ACIS geometry
//
//- Owner: Steve Storm
//- Checked by:
//- Version: $Id:

// ********** BEGIN ACIS INCLUDES             **********
#ifdef UNIX_PORT
#include "undefwin_.h"
#endif
#if CUBIT_ACIS_VERSION < 1100
#include "kernel/acis.hxx"
#include "kernel/kernapi/api/api.hxx"
#include "kernel/kernapi/api/kernapi.hxx"
#include "kernel/kerndata/geom/point.hxx"
#include "kernel/kerndata/geom/transfrm.hxx"
#include "kernel/kerndata/top/edge.hxx"
#include "kernel/kerndata/top/body.hxx"
#include "kernel/kerndata/top/face.hxx"
#include "kernel/kerndata/top/loop.hxx"
#include "kernel/kerndata/top/wire.hxx"
#include "kernel/kernutil/law/law.hxx"
#include "kernel/kerndata/lists/lists.hxx"
#include "kernel/spline/api/spl_api.hxx"
#include "cover/kernapi/api/coverapi.hxx"
#include "blend/kernapi/api/blendapi.hxx"
#include "intersct/kernapi/api/intrapi.hxx"
#include "constrct/kernapi/api/cstrapi.hxx"
#include "skin/kernapi/api/skinapi.hxx"
#include "offset/kernapi/api/ofstapi.hxx"
#include "boolean/kernapi/api/boolapi.hxx"
#include "baseutil/vector/position.hxx"
#include "baseutil/vector/unitvec.hxx"
#include "baseutil/logical.h"

#ifdef ACIS_LOCAL_OPS
  #include "lop_husk/api/lop_api.hxx"
  #include "lopt_husk/api/loptapi.hxx"
  #include "rem_husk/api/rem_api.hxx"
#endif

#include "kernel/geomhusk/copyent.hxx" // copy single entity

#else
#include "acis.hxx"
#include "api.hxx"
#include "kernapi.hxx"
#include "transfrm.hxx"
#include "vertex.hxx"
#include "edge.hxx"
#include "coedge.hxx"
#include "loop.hxx"
#include "face.hxx"
#include "body.hxx"
#include "wire.hxx"
#include "law.hxx"
#include "lists.hxx"
#include "spl_api.hxx"
#include "intrapi.hxx"
#include "cstrapi.hxx"
#include "skinapi.hxx"
#include "ofstapi.hxx"
#include "boolapi.hxx"
#include "position.hxx"
#include "unitvec.hxx"
#include "logical.h"
#include "eulerapi.hxx" // api_separate_body
#include "shl_api.hxx"  // api_sheet_thicken
#include "blendapi.hxx" // api_chamfer_edges, api_blend_edges
#include "ablapi.hxx"   // api_blend_edges_pos_rad
#include "ga_api.hxx"
#include "at_ent.hxx" // to be deleted
#include "at_int.hxx"
#include "point.hxx"
#include "getbox.hxx" //
#include "utils.hxx" // sg_body_to_1d
#include "heal_api.hxx"
#include "rbi_api.hxx"

#if CUBIT_ACIS_VERSION >= 1600
#include "lop_opts.hxx"
#endif

#ifdef ACIS_LOCAL_OPS
  #include "lop_api.hxx"
  #include "loptapi.hxx"
  #include "rem_api.hxx"
#endif

#include "copyent.hxx" // copy_single_entity
#include "acistype.hxx" // is_closed_solid_body
#endif

#ifdef UNIX_PORT
#include "defwin_.h"
#endif
// ********** END ACIS INCLUDES               **********

// ********** BEGIN CUBIT INCLUDES            **********
#include "AcisTweakTool.hpp"

#include "AcisQueryEngine.hpp"
#include "AcisModifyEngine.hpp"
#include "GeometryQueryTool.hpp"
#include "GeometryModifyTool.hpp"
#include "BodyACIS.hpp"
#include "Body.hpp"
#include "SurfaceACIS.hpp"
#include "LoopACIS.hpp"
#include "CurveACIS.hpp"
#include "PointACIS.hpp"

#include "DLIList.hpp"

#include "attrib_cubit_owner.hpp"

#include "CubitUtil.hpp"

#include "GfxDebug.hpp"
#include "GfxPreview.hpp"
#include "AcisDrawTool.hpp"
#include "CubitPlane.hpp"
#include "GMem.hpp"

// ********** END CUBIT INCLUDES              **********

AcisTweakTool* AcisTweakTool::instance_ = 0;

// Method: instance
// provides access to the unique model for this execution.
// sets up this instance on first access
AcisTweakTool* AcisTweakTool::instance()
{
  if (instance_ == 0) {
    instance_ = new AcisTweakTool();
  }
  return instance_;
}

AcisTweakTool::AcisTweakTool()
{
#ifdef ACIS_LOCAL_OPS
   api_initialize_local_ops();
   api_initialize_lopt_ops();
#endif
   api_initialize_blending();

   // For convenience, setup pointers to AcisQueryEngine and AcisModifyEngine
   AQE = AcisQueryEngine::instance();
   AME = AcisModifyEngine::instance();
}

AcisTweakTool::~AcisTweakTool()
{
#ifdef ACIS_LOCAL_OPS
  api_terminate_local_ops();
  api_terminate_lopt_ops();
#endif
  api_terminate_blending();
}

//=============================================================================
// Function   : tweak_chamfer
// Member Type: PUBLIC
// Description: Chamfer curves on solid bodies.  The left and right offsets are
//              with respect to the curve direction.  If the given right offset
//              is negative, the left offset is used.  Users can preview to
//              clarify the meaning of left and right.
// Author     : Eric Nielsen
// Date       :
//=============================================================================
CubitStatus
AcisTweakTool::tweak_chamfer( DLIList<Curve*> &curve_list,
                              double left_offset,
                              DLIList<BodySM*> &new_bodysm_list,
                              double right_offset,
                              CubitBoolean keep_old_body,
                              CubitBoolean preview )
{
  int i;
  outcome result;

  BodySM *body_ptr;

  BODY *BODY_ptr;
  BODY *copied_BODY_ptr;

  int delete_attribs =
    (GeometryModifyTool::instance()->get_new_ids() || keep_old_body);

  if( right_offset <= 0.0 )
    right_offset = left_offset;

  // Copy the incoming ref_edge_list since we will be pulling
  // edges out of it.
  DLIList<CurveACIS*> copied_ref_edge_list(curve_list.size());
  CAST_LIST( curve_list, copied_ref_edge_list, CurveACIS );
  if (curve_list.size() != copied_ref_edge_list.size())
  {
    PRINT_ERROR("Non-ACIS Curve found\n");
    return CUBIT_FAILURE;
  }

  copied_ref_edge_list.reset();
  while( copied_ref_edge_list.size() )
  {
    DLIList<EDGE*> EDGE_list;
    if( AME->get_copied_EDGES_of_body( copied_ref_edge_list, EDGE_list,
      copied_BODY_ptr ) == CUBIT_FAILURE )
    {
      // Return success if any bodies were created
      return new_bodysm_list.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
    }

    // Get original Body and BODY
    body_ptr = AQE->get_body_sm_of_ENTITY( copied_BODY_ptr );
    BODY_ptr = dynamic_cast<BodyACIS*>(body_ptr)->get_BODY_ptr();

    // Now, blend the edges on this body
    if( chamfer_edges( EDGE_list, left_offset, right_offset ) == CUBIT_FAILURE )
    {
      api_delent(copied_BODY_ptr);
      // Return success if any bodies were created
      return new_bodysm_list.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
    }

    if( !preview )
    {
      // If we've made it this far, the copied_BODY has been
      // modified and we can update it in CUBIT

      // Now cleanout the owner attributes from the copied BODY, if required
      if( delete_attribs )
        AQE->remove_cubit_owner_attrib_in_BODY( copied_BODY_ptr );

      BodySM* new_body_ptr = AME->get_new_Body( body_ptr, BODY_ptr,
        copied_BODY_ptr, keep_old_body );

      if( new_body_ptr )
        new_bodysm_list.append( new_body_ptr );
    }
    else
    {
      GfxPreview::clear(); // should previews always be cleared?

      ENTITY_LIST face_list;
      api_get_faces( copied_BODY_ptr, face_list);
      AcisBridge *ab_face_ptr = NULL;
      for( i=0; i<face_list.count(); i++ )
      {
        ab_face_ptr = ATTRIB_CUBIT_OWNER::cubit_owner( face_list[i] );
        if( !ab_face_ptr )
        {
          // Draw this face
          AcisDrawTool::instance()->draw_FACE( (FACE*)face_list[i], CUBIT_BLUE );
        }
      }

      api_delent(copied_BODY_ptr);
    }
  }

  if( preview )
    GfxPreview::flush();

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : tweak_chamfer
// Member Type: PUBLIC
// Description: Chamfer vertices on solid or sheet bodies.  On a solid body
//              there can be up to 3 offsets; on a sheet body up to 2 offsets.
//              The offsets are in the direction of the supplied edges.  If
//              multiple vertices are supplied, only one offset value is
//              allowed and the edges are not used.
// Author     : Steve Storm
// Date       : 03/28/05
//=============================================================================
CubitStatus
AcisTweakTool::tweak_chamfer( DLIList<Point*> &point_list,
                              double offset1,
                              DLIList<BodySM*> &new_bodysm_list,
                              Curve *edge1,
                              double offset2,
                              Curve *edge2,
                              double offset3,
                              Curve *edge3,
                              CubitBoolean keep_old_body,
                              CubitBoolean preview )
{
  // Sort out vertices between sheet and solid bodies
  DLIList<Point*> solid_points, sheet_points;
  if( sort_points_by_body_type( point_list, solid_points, sheet_points ) ==
    CUBIT_FAILURE )
    return CUBIT_FAILURE;

  // Do simple forms
  if( edge1 == NULL || offset2 <= 0.0 )
  {
    if( tweak_chamfer_solid( solid_points, offset1, new_bodysm_list,
      keep_old_body, preview )== CUBIT_FAILURE )
      return CUBIT_FAILURE;
    return tweak_chamfer_fillet_sheet( sheet_points, offset1, 1,
      new_bodysm_list, keep_old_body, preview );
  }

  if( solid_points.size() > 1 || sheet_points.size() > 1 )
  {
    PRINT_ERROR( "cannot chamfer multiple vertices with a variable radius.\n" );
    return CUBIT_FAILURE;
  }

  if( solid_points.size() )
  {
    Point *point_ptr = solid_points.get();
    BodySM *bodysm_ptr = NULL;

    if( tweak_chamfer_solid( point_ptr, offset1, edge1, offset2, edge2,
      offset3, edge3, bodysm_ptr, keep_old_body, preview ) == CUBIT_FAILURE )
      return CUBIT_FAILURE;

    new_bodysm_list.append( bodysm_ptr );
    return CUBIT_SUCCESS;
  }

  if( sheet_points.size() )
  {
    Point *point_ptr = sheet_points.get();
    BodySM *bodysm_ptr = NULL;

    if( tweak_chamfer_sheet( point_ptr, offset1, edge1, offset2, edge2,
      bodysm_ptr, keep_old_body, preview ) == CUBIT_FAILURE )
      return CUBIT_FAILURE;

    new_bodysm_list.append( bodysm_ptr );
    return CUBIT_SUCCESS;
  }

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : tweak_fillet
// Member Type: PUBLIC
// Description: Create a round fillet (or blend) at the given curves on solid
//              bodies.
// Author     : Eric Nielsen
// Date       :
//=============================================================================
CubitStatus
AcisTweakTool::tweak_fillet( DLIList<Curve*> &curve_list,
                             double radius,
                             DLIList<BodySM*> &new_bodysm_list,
                             bool keep_old_body,
                             CubitBoolean preview )
{
  int i;
  outcome result;

  BodySM *body_ptr;

  BODY *BODY_ptr;
  BODY *copied_BODY_ptr;

  int delete_attribs =
    (GeometryModifyTool::instance()->get_new_ids() || keep_old_body);

  // Copy the incoming ref_edge_list since we will be pulling
  // edges out of it.
  DLIList<CurveACIS*> copied_ref_edge_list(curve_list.size());
  CAST_LIST( curve_list, copied_ref_edge_list, CurveACIS );
  if (curve_list.size() != copied_ref_edge_list.size())
  {
    PRINT_ERROR("Non-ACIS Curve encountered\n");
    return CUBIT_FAILURE;
  }

  copied_ref_edge_list.reset();
  while( copied_ref_edge_list.size() )
  {
    DLIList<EDGE*> EDGE_list;
    if( AME->get_copied_EDGES_of_body( copied_ref_edge_list, EDGE_list,
      copied_BODY_ptr ) == CUBIT_FAILURE )
    {
      // Return success if any bodies were created
      return new_bodysm_list.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
    }

    // Get original Body and BODY
    body_ptr = AQE->get_body_sm_of_ENTITY( copied_BODY_ptr );
    BODY_ptr = dynamic_cast<BodyACIS*>(body_ptr)->get_BODY_ptr();

    if( blend_edges( EDGE_list, radius ) == CUBIT_FAILURE )
    {
      AQE->ACIS_API_error(result);
      api_delent(copied_BODY_ptr);
      // Return success if any bodies were created
      return new_bodysm_list.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
    }

    if( !preview )
    {
      // If we've made it this far, the copied_BODY has been
      // modified and we can update it in CUBIT

      // Now cleanout the owner attributes from the copied BODY, if required
      if( delete_attribs )
        AQE->remove_cubit_owner_attrib_in_BODY(copied_BODY_ptr);

      BodySM* new_body_ptr = AME->get_new_Body( body_ptr, BODY_ptr, copied_BODY_ptr,
        keep_old_body );

      if (new_body_ptr)
        new_bodysm_list.append( new_body_ptr );
    }
    else
    {
      GfxPreview::clear();

      ENTITY_LIST face_list;
      api_get_faces( copied_BODY_ptr, face_list);
      AcisBridge *ab_face_ptr = NULL;
      for( i=0; i<face_list.count(); i++ )
      {
        ab_face_ptr = ATTRIB_CUBIT_OWNER::cubit_owner( face_list[i] );
        if( !ab_face_ptr )
        {
          // Draw this face
          AcisDrawTool::instance()->draw_FACE( (FACE*)face_list[i], CUBIT_BLUE );
        }
      }

      api_delent(copied_BODY_ptr);
    }
  }

  if( preview )
    GfxPreview::flush();

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : tweak_fillet
// Member Type: PUBLIC
// Description: Create a round fillet (or blend) at the given curves on a solid
//              body.  The fillet has a variable radius from the start to the
//              end of the curve.
// Author     : Steve Storm
// Date       : 03/28/05
//=============================================================================
CubitStatus AcisTweakTool::tweak_fillet( Curve *curve_ptr,
                                         double start_radius,
                                         double end_radius,
                                         BodySM *&new_bodysm_ptr,
                                         CubitBoolean keep_old_body,
                                         CubitBoolean preview )
{
  outcome result;

  BodySM *body_ptr;

  BODY *BODY_ptr;
  BODY *copied_BODY_ptr;

  int delete_attribs =
    (GeometryModifyTool::instance()->get_new_ids() || keep_old_body);

  DLIList<CurveACIS*> copied_curve_list(1);
  CurveACIS *curve_acis_ptr = CAST_TO( curve_ptr, CurveACIS );
  if( curve_acis_ptr == NULL )
  {
    PRINT_ERROR("Non-ACIS Curve encountered\n");
    return CUBIT_FAILURE;
  }
  copied_curve_list.append( curve_acis_ptr );

  DLIList<EDGE*> EDGE_list;
  if( AME->get_copied_EDGES_of_body( copied_curve_list, EDGE_list,
    copied_BODY_ptr ) == CUBIT_FAILURE )
  {
    return CUBIT_FAILURE;
  }

  // Get original Body and BODY
  body_ptr = AQE->get_body_sm_of_ENTITY( copied_BODY_ptr );
  BODY_ptr = dynamic_cast<BodyACIS*>(body_ptr)->get_BODY_ptr();

  EDGE *EDGE_ptr = EDGE_list.get();
  ENTITY_LIST ENTITY_list;
  ENTITY_list.add( EDGE_ptr );

  int num_fixes = 2;
  SPAposition pos_array[2];
  pos_array[0] = EDGE_ptr->start_pos();
  pos_array[1] = EDGE_ptr->end_pos();

  double fix_radii[2] = { start_radius, end_radius };

  // Setup tweak attributes so we can preserve Cubit owners
  DLIList<FACE*> pre_FACE_list;
  DLIList<EDGE*> pre_EDGE_list;
  DLIList<VERTEX*> pre_VERTEX_list;
  DLIList<AcisBridge*> ab_FACE_list, ab_EDGE_list, ab_VERTEX_list;
  assign_tweak_attribs( copied_BODY_ptr, "tweak", pre_FACE_list, ab_FACE_list,
    pre_EDGE_list, ab_EDGE_list, pre_VERTEX_list, ab_VERTEX_list );

  // Now, blend the edges on this body
  result = api_blend_edges_pos_rad( ENTITY_list, num_fixes, pos_array, fix_radii );

  if( !result.ok() )
  {
    AQE->ACIS_API_error(result);
    api_delent(copied_BODY_ptr);
    return CUBIT_FAILURE;
  }

  // Replace Cubit owners
  reassign_cubit_owners_from_tweak_attribs( copied_BODY_ptr, "tweak", 
    pre_FACE_list, ab_FACE_list, pre_EDGE_list, ab_EDGE_list, pre_VERTEX_list,
    ab_VERTEX_list );

  // Remove tweak attributes
  remove_named_attribs( copied_BODY_ptr, "tweak" );

  if( !preview )
  {
    // If we've made it this far, the copied_BODY has been
    // modified and we can update it in CUBIT

    // Now cleanout the owner attributes from the copied BODY, if required
    if( delete_attribs )
      AQE->remove_cubit_owner_attrib_in_BODY(copied_BODY_ptr);

    new_bodysm_ptr = AME->get_new_Body( body_ptr, BODY_ptr, copied_BODY_ptr,
      keep_old_body );
  }
  else
  {
    GfxPreview::clear();

    ENTITY_LIST face_list;
    api_get_faces( copied_BODY_ptr, face_list);
    AcisBridge *ab_face_ptr = NULL;
    int i;
    for( i=0; i<face_list.count(); i++ )
    {
      ab_face_ptr = ATTRIB_CUBIT_OWNER::cubit_owner( face_list[i] );
      if( !ab_face_ptr )
      {
        // Draw this face
        AcisDrawTool::instance()->draw_FACE( (FACE*)face_list[i], CUBIT_BLUE );
      }
    }

    api_delent(copied_BODY_ptr);

    GfxPreview::flush();
  }

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : tweak_fillet
// Member Type: PUBLIC
// Description: Create a round fillet (or blend) at the given vertices on sheet
//              bodies.
// Author     : Steve Storm
// Date       : 03/28/05
//=============================================================================
CubitStatus AcisTweakTool::tweak_fillet( DLIList<Point*> &point_list,
                                         double radius,
                                         DLIList<BodySM*> &new_bodysm_list,
                                         CubitBoolean keep_old_body,
                                         CubitBoolean preview )
{
  return tweak_chamfer_fillet_sheet( point_list, radius, 2, new_bodysm_list,
    keep_old_body, preview );
}

//=============================================================================
// Function   : tweak_move
// Member Type: PUBLIC
// Description: Tweak specified faces of a volume or volumes along a vector.
// Author     : Steve Storm
// Date       :
//=============================================================================
CubitStatus AcisTweakTool::tweak_move( DLIList<Surface*> &surface_list,
                                       const CubitVector &delta,
                                       DLIList<BodySM*> &new_bodysm_list,
                                       CubitBoolean keep_old_body,
                                       CubitBoolean preview )
{
#ifndef ACIS_LOCAL_OPS
  PRINT_ERROR( "The ACIS Local Operations Husk is required for moving\n"
    "       surfaces.  It has not been licensed for this installation.\n" );
  return CUBIT_FAILURE;
#endif

  int i;
  int num_faces = 0;

  BodySM* body_ptr;

  BODY* BODY_ptr;
  BODY* copied_BODY_ptr;
  outcome result;

  // Setup translation transform
  SPAtransf tr;
  SPAvector vec(delta.x(), delta.y(), delta.z() );
  tr = translate_transf( vec );

  bool delete_attribs =
    (GeometryModifyTool::instance()->get_new_ids() || keep_old_body);

  // Copy the incoming surface_list since we will be pulling
  // surfaces out of it.
  DLIList<SurfaceACIS*> copied_surface_list( surface_list.size() );
  if (!get_ACIS_surfaces( surface_list, copied_surface_list ))
    return CUBIT_FAILURE;

  copied_surface_list.reset();
  while( copied_surface_list.size() )
  {
    DLIList<FACE*> FACE_list;
    if( AME->get_copied_FACES_of_body( copied_surface_list, FACE_list,
      copied_BODY_ptr ) == CUBIT_FAILURE )
    {
      // Return success if any bodies were created
      return new_bodysm_list.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
    }

    // Get original Body and BODY
    body_ptr = AQE->get_body_sm_of_ENTITY( copied_BODY_ptr );
    BODY_ptr = dynamic_cast<BodyACIS*>(body_ptr)->get_BODY_ptr();

    // Store moved FACEs for preview
    DLIList<AcisBridge*> preview_ab_list;
    if( preview ) get_owner_list( FACE_list, preview_ab_list );

    // Now, move the surfaces on this body
    num_faces = FACE_list.size();
    FACE** move_FACES = new FACE*[num_faces];
    for( i=0; i<num_faces; i++ )
      move_FACES[i] = FACE_list.get_and_step();

    SPAposition box_l(0,0,0);
    SPAposition box_h(0,0,0);

    // Setup tweak attributes so we can preserve Cubit owners
    DLIList<FACE*> pre_FACE_list;
    DLIList<EDGE*> pre_EDGE_list;
    DLIList<VERTEX*> pre_VERTEX_list;
    DLIList<AcisBridge*> ab_FACE_list, ab_EDGE_list, ab_VERTEX_list;
    // This operation sometimes destroys owner attributes
    AcisBridge *ab_body_ptr;
    AcisBridge *ab_lump_ptr;
    AcisBridge *ab_shell_ptr;
    if( !delete_attribs )
    {
      assign_tweak_attribs( copied_BODY_ptr, "tweak", pre_FACE_list, ab_FACE_list,
         pre_EDGE_list, ab_EDGE_list, pre_VERTEX_list, ab_VERTEX_list );
      get_owner_attribs( copied_BODY_ptr, ab_body_ptr, ab_lump_ptr, ab_shell_ptr );
    }

    result = api_move_faces( num_faces, move_FACES, tr, box_l, box_h );

    delete [] move_FACES;

    // Replace Cubit owners
    if( !delete_attribs )
    {
      reassign_cubit_owners_from_tweak_attribs( copied_BODY_ptr, "tweak",
        pre_FACE_list, ab_FACE_list, pre_EDGE_list, ab_EDGE_list, pre_VERTEX_list,
        ab_VERTEX_list );
      reset_owner_attribs( copied_BODY_ptr, ab_body_ptr, ab_lump_ptr, ab_shell_ptr );
    
      // Remove tweak attributes
      remove_named_attribs( copied_BODY_ptr, "tweak" );
    }

    if( !result.ok() )
    {
      AQE->ACIS_API_error(result);
      api_delent(copied_BODY_ptr);
      // Return success if any bodies were created
      return new_bodysm_list.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
    }

    if( !preview )
    {
      // Now cleanout the owner attributes from the copied BODY, if required
      if( delete_attribs )
        AQE->remove_cubit_owner_attrib_in_BODY(copied_BODY_ptr);

      // If we've made it this far, the copied_BODY has been
      // modified and we can update it in CUBIT
      BodySM* new_body_ptr = 
        AME->get_new_Body( body_ptr, BODY_ptr, copied_BODY_ptr, keep_old_body );

      if( new_body_ptr )
        new_bodysm_list.append( new_body_ptr );
    }
    else
    {
      GfxPreview::clear();

      // Draw preview EDGEs
      draw_tweak_preview_omt( copied_BODY_ptr, CUBIT_TRUE, &preview_ab_list );
      api_delent( copied_BODY_ptr );
    }
  }

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : tweak_move
// Member Type: PUBLIC
// Description: Tweak specified curves of a sheet body along a vector.
// Author     : Steve Storm
// Date       : 03/28/05
//=============================================================================
CubitStatus AcisTweakTool::tweak_move( DLIList<Curve*> &input_curve_list,
                                       const CubitVector &delta,
                                       DLIList<BodySM*> &new_bodysm_list,
                                       CubitBoolean keep_old_bodies,
                                       CubitBoolean preview )
{
#ifndef ACIS_LOCAL_OPS
  PRINT_ERROR( "The ACIS Local Operations Husk is required for tweaking\n"
    "       curves.  It has not been licensed for this installation.\n" );
  return CUBIT_FAILURE;
#endif

  int i, j;
  outcome result;

  // Setup translation transform
  SPAtransf tr;
  SPAvector vec(delta.x(), delta.y(), delta.z() );
  tr = translate_transf( vec );

  bool delete_attribs =
    (GeometryModifyTool::instance()->get_new_ids() || keep_old_bodies);

  // Get ACIS EDGEs
  DLIList<EDGE*> input_EDGE_list;
  if( get_EDGEs( input_curve_list, input_EDGE_list )== CUBIT_FAILURE )
    return CUBIT_FAILURE;

  input_EDGE_list.reset();
  while( input_EDGE_list.size() )
  {
    BODY *BODY_ptr;
    DLIList<EDGE*> removed_EDGE_list;
    DLIList<BODY*> thickened_BODY_list;
    DLIList<DLIList<EDGE*>*> output_EDGE_lists;
    DLIList<DLIList<FACE*>*> output_FACE_lists;
    DLIList<DLIList<FACE*>*> conjugate_FACE_lists;
    if( get_thickened_BODIES_of_EDGES( "tweak", input_EDGE_list,
      removed_EDGE_list, BODY_ptr, thickened_BODY_list, output_EDGE_lists,
      output_FACE_lists, conjugate_FACE_lists ) == CUBIT_FAILURE )
    {
      // Return success if any bodies were created
      return new_bodysm_list.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
    }

    BodySM *bodysm_ptr = AQE->get_body_sm_of_ENTITY( BODY_ptr );

    // Remove the surfaces on the output BODIEs
    // First copy the BODY
    BODY *copied_BODY_ptr = AME->copy_BODY(BODY_ptr, CUBIT_FALSE);

    // Store moved EDGEs for preview
    DLIList<AcisBridge*> preview_ab_list;
    if( preview ) get_owner_list( removed_EDGE_list, preview_ab_list );

    int num_faces;
    BODY *thickened_BODY_ptr;
    DLIList<FACE*> *output_FACE_list_ptr;
    DLIList<FACE*> *conjugate_FACE_list_ptr;
    thickened_BODY_list.reset();
    output_EDGE_lists.reset();
    output_FACE_lists.reset();
    conjugate_FACE_lists.reset();
    for( i=thickened_BODY_list.size(); i--; )
    {
      thickened_BODY_ptr = thickened_BODY_list.get_and_step();
      output_FACE_list_ptr = output_FACE_lists.get_and_step();
      conjugate_FACE_list_ptr = conjugate_FACE_lists.get_and_step();

      // Note - the conjugate FACE list can have duplicate entries in
      // it, as thicken doesn't always create a separate surface for
      // each curve.
      DLIList<FACE*> tweak_FACE_list;
      conjugate_FACE_list_ptr->reset();
      for( j=conjugate_FACE_list_ptr->size(); j--; )
        tweak_FACE_list.append_unique( conjugate_FACE_list_ptr->get_and_step() );

      // Keep track of the FACEs of interest, via their Cubit owners, which
      // will survive throughout the operation
      DLIList<AcisBridge*> owner_list;
      get_owner_list( *output_FACE_list_ptr, owner_list );

      num_faces = tweak_FACE_list.size();

      FACE** FACES = new FACE*[num_faces];
      output_FACE_list_ptr->reset();
      for( j=0; j<num_faces; j++ )
        FACES[j] = tweak_FACE_list.get_and_step();

      SPAposition box_l(0,0,0);
      SPAposition box_h(0,0,0);

      result = api_move_faces( num_faces, FACES, tr, box_l, box_h );

      delete [] FACES;

      if( !result.ok() )
      {
        AQE->ACIS_API_error(result);
        while( output_FACE_lists.size() ) delete output_FACE_lists.pop();
        while( conjugate_FACE_lists.size() ) delete conjugate_FACE_lists.pop();
        while( thickened_BODY_list.size() )
          api_delent( thickened_BODY_list.pop() );
        api_delent( copied_BODY_ptr );
        // Return success if any bodies were created
        return new_bodysm_list.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
      }

      // Get bodies ready to swap new surfaces for old
      if( prep_for_surface_swap( thickened_BODY_ptr, copied_BODY_ptr,
        owner_list ) == CUBIT_FAILURE )
      {
        while( output_EDGE_lists.size() ) delete output_EDGE_lists.pop();
        while( output_FACE_lists.size() ) delete output_FACE_lists.pop();
        while( conjugate_FACE_lists.size() ) delete conjugate_FACE_lists.pop();
        while( thickened_BODY_list.size() )
          api_delent( thickened_BODY_list.pop() );
        api_delent( copied_BODY_ptr );
        // Return success if any bodies were created
        return new_bodysm_list.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
      }
    } // End loop on thickened (separated) BODIES

    // Free memory
    while( output_EDGE_lists.size() ) delete output_EDGE_lists.pop();
    while( output_FACE_lists.size() ) delete output_FACE_lists.pop();
    while( conjugate_FACE_lists.size() ) delete conjugate_FACE_lists.pop();

    // Unite the thickened BODIEs back into the copied_BODY_ptr
    BODY *master;
    if( unite_BODIES( copied_BODY_ptr, thickened_BODY_list, master ) == CUBIT_FAILURE )
    {
      // If failure, the unite_BODIES function cleaned up the memory
      // Return success if any bodies were created
      return new_bodysm_list.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
    }

    // This BODY is done (whew!)

    if( !preview )
    {
      // Now cleanout the owner attributes from the copied BODY, if required
      if( delete_attribs )
        AQE->remove_cubit_owner_attrib_in_BODY( master );

      BodySM *new_body = AME->get_new_Body( bodysm_ptr, BODY_ptr, master,
        keep_old_bodies );

      if( new_body )
        new_bodysm_list.append( new_body );
    }
    else
    {
      GfxPreview::clear();

      // Draw preview EDGEs
      draw_tweak_preview_omt( copied_BODY_ptr, CUBIT_TRUE, &preview_ab_list );
      api_delent( master );
    }
  }

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : tweak_offset
// Member Type: PUBLIC
// Description: Tweak specified faces of a volume or volumes by offsetting
//              those faces by the offset distance.
// Author     : Steve Storm
// Date       :
//=============================================================================
CubitStatus
AcisTweakTool::tweak_offset( DLIList<Surface*> &surface_list,
                             double offset_distance,
                             DLIList<BodySM*> &new_bodysm_list,
                             CubitBoolean keep_old_body,
                             CubitBoolean preview )
{
#ifndef ACIS_LOCAL_OPS
   PRINT_ERROR( "The ACIS Local Operations Husk is required for offsetting\n"
      "       surfaces.  It has not been licensed for this installation.\n" );
   return CUBIT_FAILURE;
#endif

   int i;
   outcome result;

   // Copy the incoming surface_list since we will be pulling surfaces out of
   // it.

   BodySM *body_ptr;

   BODY *BODY_ptr;
   BODY *copied_BODY_ptr;

   bool delete_attribs =
      (GeometryModifyTool::instance()->get_new_ids() || keep_old_body);

   // Copy the incoming surface_list since we will be pulling surfaces out of
   // it.
   DLIList<SurfaceACIS*> copied_surface_list( surface_list.size() );
   if (!get_ACIS_surfaces( surface_list, copied_surface_list ))
     return CUBIT_FAILURE;

   copied_surface_list.reset();
   while( copied_surface_list.size() )
   {
     DLIList<FACE*> FACE_list;
     if( AME->get_copied_FACES_of_body( copied_surface_list, FACE_list,
       copied_BODY_ptr ) == CUBIT_FAILURE )
     {
       // Return success if any bodies were created
       return new_bodysm_list.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
     }

     // Get original Body and BODY
     body_ptr = AQE->get_body_sm_of_ENTITY( copied_BODY_ptr );
     BODY_ptr = dynamic_cast<BodyACIS*>(body_ptr)->get_BODY_ptr();

     // Store offset FACEs for preview
     DLIList<AcisBridge*> preview_ab_list;
     if( preview ) get_owner_list( FACE_list, preview_ab_list );

     // Now, offset the surfaces on this body
     int num_offsets = FACE_list.size();
     FACE** offset_FACES = new FACE*[num_offsets];
     for( i=0; i<num_offsets; i++ )
       offset_FACES[i] = FACE_list.get_and_step();

     SPAposition box_l(0,0,0);
     SPAposition box_h(0,0,0);

     //This has been added as a temporary fix for bug #5559.
     //Lump attrib should persist across offset operation, but does not
     //ENTITY_LIST lump_list;
     //api_get_lumps(copied_BODY_ptr, lump_list);
     //AcisBridge *lump_owner_before = ATTRIB_CUBIT_OWNER::cubit_owner(lump_list[0]);

     // Setup tweak attributes so we can preserve Cubit owners
     DLIList<FACE*> pre_FACE_list;
     DLIList<EDGE*> pre_EDGE_list;
     DLIList<VERTEX*> pre_VERTEX_list;
     DLIList<AcisBridge*> ab_FACE_list, ab_EDGE_list, ab_VERTEX_list;
     // This operation sometimes destroys owner attributes
     AcisBridge *ab_body_ptr;
     AcisBridge *ab_lump_ptr;
     AcisBridge *ab_shell_ptr;
     if( !delete_attribs )
     {
       assign_tweak_attribs( copied_BODY_ptr, "tweak", pre_FACE_list, ab_FACE_list,
          pre_EDGE_list, ab_EDGE_list, pre_VERTEX_list, ab_VERTEX_list );
       get_owner_attribs( copied_BODY_ptr, ab_body_ptr, ab_lump_ptr, ab_shell_ptr );
     }

     result = api_offset_faces( num_offsets, offset_FACES, offset_distance,
                                box_l, box_h );
     delete [] offset_FACES;

     if( !result.ok() )
     {
       AQE->ACIS_API_error(result);
       api_delent(copied_BODY_ptr);
       // Return success if any bodies were created
       return new_bodysm_list.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
     }

     if( !preview )
     {
       // Now cleanout the owner attributes from the copied BODY, if required
       if( delete_attribs )
         AQE->remove_cubit_owner_attrib_in_BODY( copied_BODY_ptr );
       else
       {
         /*
         //continued temporary fix for bug #5559
         lump_list.clear();
         api_get_lumps(copied_BODY_ptr, lump_list);
         AcisBridge *lump_owner_after = ATTRIB_CUBIT_OWNER::cubit_owner(lump_list[0]);
         if( lump_owner_after == NULL )
           ATTRIB_CUBIT_OWNER::set_cubit_owner( lump_list[0], lump_owner_before );
           */

         // Replace Cubit owners
         if( !delete_attribs )
         {
           reassign_cubit_owners_from_tweak_attribs( copied_BODY_ptr, "tweak",
             pre_FACE_list, ab_FACE_list, pre_EDGE_list, ab_EDGE_list, pre_VERTEX_list,
             ab_VERTEX_list );
           reset_owner_attribs( copied_BODY_ptr, ab_body_ptr, ab_lump_ptr, ab_shell_ptr );
         
           // Remove tweak attributes
           remove_named_attribs( copied_BODY_ptr, "tweak" );
         }
       }

       BodySM *new_body_ptr = AME->get_new_Body( body_ptr, BODY_ptr,
         copied_BODY_ptr, keep_old_body );

       if( new_body_ptr )
         new_bodysm_list.append( new_body_ptr );
     }
     else
     {
       GfxPreview::clear();

       // Draw preview EDGEs
       draw_tweak_preview_omt( copied_BODY_ptr, CUBIT_TRUE, &preview_ab_list );
       api_delent( copied_BODY_ptr );
     }
   }

   return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : tweak_offset
// Member Type: PUBLIC
// Description: Tweak specified curves of a sheet body or bodies by offsetting
//              those curves by the offset distance.
// Author     :
// Date       :
//=============================================================================
CubitStatus AcisTweakTool::tweak_offset( DLIList<Curve*> &input_curve_list,
                                         double offset_distance,
                                         DLIList<BodySM*> &new_bodysm_list,
                                         CubitBoolean keep_old_bodies,
                                         CubitBoolean preview )
{
#ifndef ACIS_LOCAL_OPS
  PRINT_ERROR( "The ACIS Local Operations Husk is required for tweaking\n"
    "       curves.  It has not been licensed for this installation.\n" );
  return CUBIT_FAILURE;
#endif

  int i, j;
  outcome result;

  bool delete_attribs =
    (GeometryModifyTool::instance()->get_new_ids() || keep_old_bodies);

  // Get ACIS EDGEs
  DLIList<EDGE*> input_EDGE_list;
  if( get_EDGEs( input_curve_list, input_EDGE_list ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  input_EDGE_list.reset();
  while( input_EDGE_list.size() )
  {
    BODY *BODY_ptr;
    DLIList<EDGE*> removed_EDGE_list;
    DLIList<BODY*> thickened_BODY_list;
    DLIList<DLIList<EDGE*>*> output_EDGE_lists;
    DLIList<DLIList<FACE*>*> output_FACE_lists;
    DLIList<DLIList<FACE*>*> conjugate_FACE_lists;
    if( get_thickened_BODIES_of_EDGES( "tweak", input_EDGE_list,
      removed_EDGE_list, BODY_ptr, thickened_BODY_list, output_EDGE_lists,
      output_FACE_lists, conjugate_FACE_lists ) == CUBIT_FAILURE )
    {
      // Return success if any bodies were created
      return new_bodysm_list.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
    }

    // Get the BodySM from the BODY
    BodySM *bodysm_ptr = AQE->get_body_sm_of_ENTITY( BODY_ptr );

    // Remove the surfaces on the output BODIEs

    // Copy the BODY with all owner attributes intact
    BODY *copied_BODY_ptr = AME->copy_BODY( BODY_ptr, CUBIT_FALSE );

    // Store offset EdGEs for preview
    DLIList<AcisBridge*> preview_ab_list;
    if( preview ) get_owner_list( removed_EDGE_list, preview_ab_list );

    int num_faces;
    BODY *thickened_BODY_ptr;
    DLIList<FACE*> *output_FACE_list_ptr;
    DLIList<FACE*> *conjugate_FACE_list_ptr;
    thickened_BODY_list.reset();
    output_EDGE_lists.reset();
    output_FACE_lists.reset();
    conjugate_FACE_lists.reset();
    for( i=thickened_BODY_list.size(); i--; )
    {
      thickened_BODY_ptr = thickened_BODY_list.get_and_step();
      output_FACE_list_ptr = output_FACE_lists.get_and_step();
      conjugate_FACE_list_ptr = conjugate_FACE_lists.get_and_step();

      // Note - the conjugate FACE list can have duplicate entries in
      // it, as thicken doesn't always create a separate surface for
      // each curve.
      DLIList<FACE*> tweak_FACE_list;
      conjugate_FACE_list_ptr->reset();
      for( j=conjugate_FACE_list_ptr->size(); j--; )
        tweak_FACE_list.append_unique( conjugate_FACE_list_ptr->get_and_step() );

      // Keep track of the FACEs of interest, via their Cubit owners, which
      // will survive throughout the operation
      DLIList<AcisBridge*> owner_list;
      get_owner_list( *output_FACE_list_ptr, owner_list );

      num_faces = tweak_FACE_list.size();

      FACE** FACES = new FACE*[num_faces];
      output_FACE_list_ptr->reset();
      for( j=0; j<num_faces; j++ )
        FACES[j] = tweak_FACE_list.get_and_step();

      SPAposition box_l(0,0,0);
      SPAposition box_h(0,0,0);

#ifdef ACIS_LOCAL_OPS
      result = api_offset_faces( num_faces, FACES, offset_distance,
       box_l, box_h );
#endif

      delete [] FACES;

      if( !result.ok() )
      {
        AQE->ACIS_API_error(result);
        while( output_EDGE_lists.size() ) delete output_EDGE_lists.pop();
        while( output_FACE_lists.size() ) delete output_FACE_lists.pop();
        while( conjugate_FACE_lists.size() ) delete conjugate_FACE_lists.pop();
        while( thickened_BODY_list.size() )
          api_delent( thickened_BODY_list.pop() );
        api_delent( copied_BODY_ptr );
        // Return success if any bodies were created
        return new_bodysm_list.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
      }

      // Get bodies ready to swap new surfaces for old
      if( prep_for_surface_swap( thickened_BODY_ptr, copied_BODY_ptr,
        owner_list ) == CUBIT_FAILURE )
      {
        while( output_EDGE_lists.size() ) delete output_EDGE_lists.pop();
        while( output_FACE_lists.size() ) delete output_FACE_lists.pop();
        while( conjugate_FACE_lists.size() ) delete conjugate_FACE_lists.pop();
        while( thickened_BODY_list.size() )
          api_delent( thickened_BODY_list.pop() );
        api_delent( copied_BODY_ptr );
        // Return success if any bodies were created
        return new_bodysm_list.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
      }
    } // End loop on thickened (separated) BODIES

    // Free memory
    while( output_EDGE_lists.size() ) delete output_EDGE_lists.pop();
    while( output_FACE_lists.size() ) delete output_FACE_lists.pop();
    while( conjugate_FACE_lists.size() ) delete conjugate_FACE_lists.pop();

    // Unite the thickened BODIEs back into the copied_input_BODY_ptr
    BODY *master;
    if( unite_BODIES( copied_BODY_ptr, thickened_BODY_list, master )
      == CUBIT_FAILURE )
    {
      // If failure, the unite_BODIES function cleaned up the memory
      // Return success if any bodies were created
      return new_bodysm_list.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
    }

    // This BODY is done (whew!)

    if( !preview )
    {
      // Now cleanout the owner attributes from the copied BODY, if required
      if( delete_attribs )
        AQE->remove_cubit_owner_attrib_in_BODY( master );

      BodySM *new_body = AME->get_new_Body( bodysm_ptr, BODY_ptr, master,
        keep_old_bodies );

      if( new_body )
        new_bodysm_list.append( new_body );
    }
    else
    {
      GfxPreview::clear();

      // Draw preview EDGEs
      draw_tweak_preview_omt( master, CUBIT_TRUE, &preview_ab_list );
      api_delent( master );
    }
  }

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : tweak_remove
// Member Type: PUBLIC
// Description: Function to remove surfaces from a body and then extend the
//              remaining surfaces to fill the gap or hole.
// Author     : Steve Storm
// Date       : 03/28/05
//=============================================================================
CubitStatus AcisTweakTool::tweak_remove( DLIList<Surface*> &surface_list,
                                         DLIList<BodySM*> &new_bodysm_list,
                                         CubitBoolean extend_adjoining,
                                         CubitBoolean keep_surface,
                                         CubitBoolean keep_old_body,
                                         CubitBoolean preview )
{
#ifndef ACIS_LOCAL_OPS
   if( extend_adjoining )
   {
      PRINT_ERROR( "The ACIS Local Operations Husk is required for extending\n"
            "       adjoining surfaces.  It has not been licensed for this\n"
            "       installation.\n" );
      return CUBIT_FAILURE;
   }
#endif

   int i;
   outcome result;

   // Copy the incoming surface_list since we will be pulling
   // surfaces out of it.
   DLIList<SurfaceACIS*> copied_surface_list(surface_list.size());
   if (!get_ACIS_surfaces( surface_list, copied_surface_list ))
     return CUBIT_FAILURE;

   BodySM *body_ptr;

   BODY *BODY_ptr;
   BODY *copied_BODY_ptr;
   FACE *FACE_ptr;
   ENTITY *copied_entity_ptr;
   DLIList<FACE*> kept_FACE_list, remove_FACE_list;
   DLIList<SurfaceACIS*> body_remove_face_list;

   bool delete_attribs =
     (GeometryModifyTool::instance()->get_new_ids() || keep_old_body);

   // Remove sometimes destroys the owner attributes, so we keep track of them
   AcisBridge *ab_body_ptr;
   AcisBridge *ab_lump_ptr;
   AcisBridge *ab_shell_ptr;

   copied_surface_list.reset();
   while( copied_surface_list.size() )
   {
     remove_FACE_list.clean_out();
     body_remove_face_list.clean_out();
     if( AME->get_copied_FACES_of_body( copied_surface_list, remove_FACE_list,
       body_remove_face_list, copied_BODY_ptr ) == CUBIT_FAILURE )
     {
       // Return success if any bodies were created
       return new_bodysm_list.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
     }

     // Get original Body and BODY
     body_ptr = AQE->get_body_sm_of_ENTITY( copied_BODY_ptr );
     BODY_ptr = dynamic_cast<BodyACIS*>(body_ptr)->get_BODY_ptr();

     // Now cleanout the owner attributes from the copied BODY, if required
     if( delete_attribs )
       AQE->remove_cubit_owner_attrib_in_BODY(copied_BODY_ptr);

     // Keep a copy of the surfaces if required
     if( keep_surface && !preview )
     {
       // Make a copy of the FACES to hold onto
       for( i=0; i<remove_FACE_list.size(); i++ )
       {
         FACE_ptr = remove_FACE_list.get_and_step();
         copy_single_entity((ENTITY*)FACE_ptr, copied_entity_ptr);
         kept_FACE_list.append( (FACE*)copied_entity_ptr );
       }
     }

     // Finally, remove the faces from this body
     if( extend_adjoining )
     {
       FACE** FACES = new FACE*[remove_FACE_list.size()];
       for( i=0; i<remove_FACE_list.size(); i++ )
         FACES[i] = remove_FACE_list.get_and_step();

       // Remove sometimes destroys owner attributes
       // Setup tweak attributes so we can preserve Cubit owners
       DLIList<FACE*> pre_FACE_list;
       DLIList<EDGE*> pre_EDGE_list;
       DLIList<VERTEX*> pre_VERTEX_list;
       DLIList<AcisBridge*> ab_FACE_list, ab_EDGE_list, ab_VERTEX_list;
       if( !delete_attribs )
       {
         get_owner_attribs( copied_BODY_ptr, ab_body_ptr, ab_lump_ptr, ab_shell_ptr );
         assign_tweak_attribs( copied_BODY_ptr, "tweak", pre_FACE_list, ab_FACE_list,
           pre_EDGE_list, ab_EDGE_list, pre_VERTEX_list, ab_VERTEX_list );
       }

       SPAposition box_l(0,0,0);
       SPAposition box_h(0,0,0);

#ifdef ACIS_LOCAL_OPS

       // Add tweak_preview attributes to adjoining FACEs to remove if needed
       if( preview )
         tag_tweak_remove_FACEs_for_preview( copied_BODY_ptr, remove_FACE_list );

#if CUBIT_ACIS_VERSION < 1600
       // When removing FACEs, we don't want ACIS to regularize the volume.
       // Here we set the EDGEs not owned by the FACEs being removed as not
       // mergeable so they persist through the removal operation.

       // Get all the EDGEs on the BODY
       ENTITY_LIST body_edges;
       api_get_edges( copied_BODY_ptr, body_edges );

       // Get all the EDGEs on the FACE to remove
       ENTITY_LIST surf_edges;
       for( i=0; i<remove_FACE_list.size(); i++ )
         api_get_edges( remove_FACE_list.get_and_step(), surf_edges );

       // Remove all EDGEs of the FACE to remove
       body_edges.remove( surf_edges );

       // Set the EDGEs as not mergeable
       api_set_no_merge_attrib( body_edges );
#endif

       if( !preview )
         PRINT_INFO( "Removing %d surface(s) from volume...\n", remove_FACE_list.size() );
       else
         PRINT_INFO( "Previewing removal of %d surface(s) from volume...\n", remove_FACE_list.size() );

       result = api_remove_faces( remove_FACE_list.size(), FACES, box_l, box_h );

#if CUBIT_ACIS_VERSION < 1600
       // Remove the merge attributes
       body_edges.clear();
       api_get_edges( copied_BODY_ptr, body_edges );
       api_remove_no_merge_attrib( body_edges );
#endif

       delete [] FACES;

       if( !result.ok() )
       {
         PRINT_INFO("result was bad....\n");
         AQE->ACIS_API_error(result);
         api_delent(copied_BODY_ptr);
         // Return success if any bodies were created
         return new_bodysm_list.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
       }
       else if( !preview )
         PRINT_INFO( "Successfully removed surfaces from volume\n" );
       else
         PRINT_INFO( "Successfully previewed removal of surfaces from volume\n" );


       if( !delete_attribs )
       {
         reset_owner_attribs( copied_BODY_ptr, ab_body_ptr, ab_lump_ptr, ab_shell_ptr );
         
         // Replace Cubit owners
         reassign_cubit_owners_from_tweak_attribs( copied_BODY_ptr, "tweak",
           pre_FACE_list, ab_FACE_list, pre_EDGE_list, ab_EDGE_list, pre_VERTEX_list,
           ab_VERTEX_list );

         // Remove tweak attributes
         remove_named_attribs( copied_BODY_ptr, "tweak" );
       }
#endif
     }
     else
     {
       // Just unhook the faces (if previewing do this anyway to check for
       // errors)
       for( i=0; i<remove_FACE_list.size(); i++ )
       {
        
         //make sure we're not trying to unhook the last face
         ENTITY_LIST tmp_face_list;
         api_get_faces( copied_BODY_ptr, tmp_face_list);

         if( tmp_face_list.count() == 1 )
         {
           AcisBridge *tmp_bridge = ATTRIB_CUBIT_OWNER::cubit_owner( copied_BODY_ptr );
           if( tmp_bridge )
           {
             BodyACIS *tmp_acis_body = CAST_TO( tmp_bridge, BodyACIS );
             Body *tmp_body = CAST_TO(tmp_acis_body->topology_entity(), Body ); 
             PRINT_WARNING("Will not remove ALL faces from body %d\n", tmp_body->id() ); 
           }
           continue;
         }

         BODY *new_BODY;
         FACE_ptr = remove_FACE_list.get_and_step();

         result = api_unhook_face( FACE_ptr, new_BODY );

         if( !result.ok() )
         {
           AQE->ACIS_API_error(result);
           api_delent(copied_BODY_ptr);
           // Return success if any bodies were created
           return new_bodysm_list.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
         }

         if( preview )
         {
           GfxPreview::clear();

           AcisDrawTool::instance()->draw_ENTITY( new_BODY, CUBIT_BLUE,
           CUBIT_FALSE, CUBIT_TRUE );
         }

         // This function puts the pulled-off FACE into a sheet
         // body.  No need to keep this - we copied the surfaces
         // before if the keep_surface flag is on.
         api_delent( new_BODY );
       }
     }

     // If we are keeping the pulled-off surfaces, just make free RefFaces from
     // them.
     if( keep_surface && !preview )
     {
       for( i=0; i<kept_FACE_list.size(); i++ )
       {
         FACE* kept_FACE = kept_FACE_list.get_and_step();
         FACE *face_list[1];
         face_list[0] = kept_FACE;
         BODY *sheet_BODY_ptr = NULL;
         result = api_sheet_from_ff( 1, face_list, sheet_BODY_ptr );
         api_body_to_2d( sheet_BODY_ptr );

         if( !result.ok() )
         {
           PRINT_ERROR("Problem with 'keepsurface' option\n" );
         }
         else
         {
           BodySM *body_ptr = AQE->populate_topology_bridges( sheet_BODY_ptr );
           DLIList<Surface*> tmp_surfs;
           body_ptr->surfaces( tmp_surfs );
           Surface *surface_ptr = tmp_surfs.get();
           GeometryModifyTool::instance()->make_Body( surface_ptr );
         }
       }
     }

     // If we've made it this far, the copied_BODY has been modified and we can
     // update it in CUBIT
     if( !preview )
     {
       BodySM *new_body = AME->get_new_Body( body_ptr, BODY_ptr,
         copied_BODY_ptr, keep_old_body );

       if( new_body )
         new_bodysm_list.append( new_body );
     }
     else if( extend_adjoining ) // Only case we need to preview
     {
       GfxPreview::clear();

       draw_tweak_preview_tagged_FACEs( copied_BODY_ptr, CUBIT_TRUE );
       api_delent( copied_BODY_ptr );
     }
   }

   return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : tweak_remove
// Member Type: PUBLIC
// Description: Function to remove curves from a sheet body and then extend the
//              remaining curves or fill the gap or hole.
// Author     : Steve Storm
// Date       : 03/28/05
//=============================================================================
CubitStatus
AcisTweakTool::tweak_remove( DLIList<Curve*> &input_curve_list,
                             DLIList<BodySM*> &new_bodysm_list,
                             CubitBoolean keep_old_bodies,
                             CubitBoolean preview )
{
  int i, j, k;
  outcome result;

  bool delete_attribs =
    (GeometryModifyTool::instance()->get_new_ids() || keep_old_bodies);

  // We have a very efficient way to remove internal loops, so if all curves
  // being removed on a given body are internal loops, we do that with a
  // special case first.  We will loop on each set of EDGEs from a common
  // BODY.
  DLIList<BODY*> BODY_list;
  DLIList<DLIList<EDGE*>*> BODY_EDGE_lists;
  if( get_EDGES_by_BODY( input_curve_list, BODY_list, BODY_EDGE_lists )
    == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  BODY *BODY_ptr;
  BodySM *bodysm_ptr;
  DLIList<EDGE*> *BODY_EDGE_list_ptr;
  BODY_list.reset();
  BODY_EDGE_lists.reset();
  for( i=BODY_EDGE_lists.size(); i--; )
  {
    BODY_ptr = BODY_list.get_and_step();
    BODY_EDGE_list_ptr = BODY_EDGE_lists.get_and_step();

    // Get the BodySM from the BODY
    bodysm_ptr = AQE->get_body_sm_of_ENTITY( BODY_ptr );

    // Copy the BODY with all owner attributes intact
    BODY *copied_BODY_ptr = AME->copy_BODY( BODY_ptr, CUBIT_FALSE );

    DLIList<LOOP*> LOOP_list;
    if( all_complete_internal_loops( *BODY_EDGE_list_ptr, LOOP_list )
      == CUBIT_SUCCESS )
    {
      // Efficiently remove these LOOPS

      // Get the LOOPs to remove on the copied BODY
      DLIList<AcisBridge*> owner_list;
      get_owner_list( LOOP_list, owner_list );
      DLIList<LOOP*> copied_LOOP_list;
      get_corresponding_LOOP_list( owner_list, copied_BODY_ptr, copied_LOOP_list );

      // Tag EDGEs for preview
      if( preview )
      {
        LOOP *LOOP_ptr;
        for( j=copied_LOOP_list.size(); j--; )
        {
          LOOP_ptr = copied_LOOP_list.get_and_step();
          DLIList<EDGE*> LOOP_EDGE_list;
          AQE->get_EDGEs( LOOP_ptr, LOOP_EDGE_list );
          tag_tweak_remove_FACEs_for_preview( copied_BODY_ptr, LOOP_EDGE_list );
        }
      }

      if( remove_LOOPs( copied_LOOP_list ) == CUBIT_SUCCESS )
      {
        if( !preview )
        {
          // Now cleanout the owner attributes from the copied BODY, if required
          if( delete_attribs )
            AQE->remove_cubit_owner_attrib_in_BODY(copied_BODY_ptr);

          BodySM *new_body = AME->get_new_Body( bodysm_ptr, BODY_ptr,
            copied_BODY_ptr, keep_old_bodies );

          if( new_body )
            new_bodysm_list.append( new_body );
        }
        else
        {
          GfxPreview::clear();

          draw_tweak_preview_tagged_FACEs( copied_BODY_ptr, CUBIT_TRUE );
          api_delent( copied_BODY_ptr );
        }
      }
      else
      {
        // Cleanup memory
        while( BODY_EDGE_lists.size() ) delete BODY_EDGE_lists.pop();
        api_delent( copied_BODY_ptr );
        // Return success if any bodies were created
        return new_bodysm_list.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
      }
    }
    // For this BODY, all curves to remove are not complete loops
    else
    {
      // Remove sometimes destroys the owner attributes
      AcisBridge *ab_body_ptr;
      AcisBridge *ab_lump_ptr;
      AcisBridge *ab_shell_ptr;

      BODY *BODY_ptr2;
      DLIList<EDGE*> removed_EDGE_list;
      DLIList<BODY*> thickened_BODY_list;
      DLIList<DLIList<EDGE*>*> output_EDGE_lists;
      DLIList<DLIList<FACE*>*> output_FACE_lists;
      DLIList<DLIList<FACE*>*> conjugate_FACE_lists;
      if( get_thickened_BODIES_of_EDGES( "tweak", *BODY_EDGE_list_ptr, 
        removed_EDGE_list, BODY_ptr2, thickened_BODY_list, output_EDGE_lists,
        output_FACE_lists, conjugate_FACE_lists ) == CUBIT_FAILURE )
      {
        // Cleanup memory
        while( BODY_EDGE_lists.size() ) delete BODY_EDGE_lists.pop();
        api_delent( copied_BODY_ptr );
        // Return success if any bodies were created
        return new_bodysm_list.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
      }

      if( BODY_ptr != BODY_ptr2 )
      {
        PRINT_ERROR( "Internal error - please report\n" );
        // Cleanup memory
        while( BODY_EDGE_lists.size() ) delete BODY_EDGE_lists.pop();
        api_delent( copied_BODY_ptr );
        // Return success if any bodies were created
        return new_bodysm_list.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
      }

      // Tag EDGEs for preview
      if( preview )
      {
        thickened_BODY_list.reset();
        output_EDGE_lists.reset();
        for( j=output_EDGE_lists.size(); j--; )
        {
          DLIList<EDGE*> *EDGE_list_ptr = output_EDGE_lists.get_and_step();
          BODY *tmp_BODY_ptr = thickened_BODY_list.get_and_step();
          tag_tweak_remove_FACEs_for_preview( tmp_BODY_ptr, *EDGE_list_ptr );
        }
      }

      // Remove the surfaces on the output BODIEs

      int num_faces;
      BODY *thickened_BODY_ptr;
      DLIList<FACE*> *output_FACE_list_ptr;
      DLIList<FACE*> *conjugate_FACE_list_ptr;
      thickened_BODY_list.reset();
      output_EDGE_lists.reset();
      output_FACE_lists.reset();
      conjugate_FACE_lists.reset();
      for( j=thickened_BODY_list.size(); j--; )
      {
        thickened_BODY_ptr = thickened_BODY_list.get_and_step();
        output_FACE_list_ptr = output_FACE_lists.get_and_step();
        conjugate_FACE_list_ptr = conjugate_FACE_lists.get_and_step();

        // Note - the conjugate FACE list can have duplicate entries in
        // it, as thicken doesn't always create a separate surface for
        // each curve.
        DLIList<FACE*> tweak_FACE_list;
        conjugate_FACE_list_ptr->reset();
        for( k=conjugate_FACE_list_ptr->size(); k--; )
          tweak_FACE_list.append_unique( conjugate_FACE_list_ptr->get_and_step() );

        // Keep track of the FACEs of interest, via their Cubit owners, which
        // will survive throughout the operation
        DLIList<AcisBridge*> owner_list;
        get_owner_list( *output_FACE_list_ptr, owner_list );

        num_faces = tweak_FACE_list.size();

        FACE** FACES = new FACE*[num_faces];
        output_FACE_list_ptr->reset();
        for( k=0; k<num_faces; k++ )
          FACES[k] = tweak_FACE_list.get_and_step();

        // Remove sometimes destroys owner attributes
        get_owner_attribs( thickened_BODY_ptr, ab_body_ptr, ab_lump_ptr, ab_shell_ptr );

        SPAposition box_l(0,0,0);
        SPAposition box_h(0,0,0);

#ifdef ACIS_LOCAL_OPS
        result = api_remove_faces( num_faces, FACES, box_l, box_h );
#endif

        delete [] FACES;

        if( !result.ok() )
        {
          AQE->ACIS_API_error(result);
          while( BODY_EDGE_lists.size() ) delete BODY_EDGE_lists.pop();
          while( output_EDGE_lists.size() ) delete output_EDGE_lists.pop();
          while( output_FACE_lists.size() ) delete output_FACE_lists.pop();
          while( conjugate_FACE_lists.size() ) delete conjugate_FACE_lists.pop();
          while( thickened_BODY_list.size() )
            api_delent( thickened_BODY_list.pop() );
          api_delent( copied_BODY_ptr );
          // Return success if any bodies were created
          return new_bodysm_list.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
        }

        // Reset the owner attributes
        reset_owner_attribs( thickened_BODY_ptr, ab_body_ptr, ab_lump_ptr, ab_shell_ptr );

        // Get bodies ready to swap new surfaces for old
        if( prep_for_surface_swap( thickened_BODY_ptr, copied_BODY_ptr,
          owner_list ) == CUBIT_FAILURE )
        {
          while( BODY_EDGE_lists.size() ) delete BODY_EDGE_lists.pop();
          while( output_EDGE_lists.size() ) delete output_EDGE_lists.pop();
          while( output_FACE_lists.size() ) delete output_FACE_lists.pop();
          while( conjugate_FACE_lists.size() ) delete conjugate_FACE_lists.pop();
          while( thickened_BODY_list.size() )
            api_delent( thickened_BODY_list.pop() );
          api_delent( copied_BODY_ptr );
          // Return success if any bodies were created
          return new_bodysm_list.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
        }
      } // End loop on thickened (separated) BODIES

      // Free memory
      while( output_EDGE_lists.size() ) delete output_EDGE_lists.pop();
      while( output_FACE_lists.size() ) delete output_FACE_lists.pop();
      while( conjugate_FACE_lists.size() ) delete conjugate_FACE_lists.pop();

      // Unite the thickened BODIEs back into the copied_BODY_ptr
      BODY *master;
      if( unite_BODIES( copied_BODY_ptr, thickened_BODY_list, master ) == CUBIT_FAILURE )
      {
        // If failure, the unite_BODIES function cleaned up the
        // thickened_BODY_list memory
        while( BODY_EDGE_lists.size() ) delete BODY_EDGE_lists.pop();
        // Return success if any bodies were created
        return new_bodysm_list.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
      }

      // This BODY is done (whew!)

      if( !preview )
      {
        // Now cleanout the owner attributes from the copied BODY, if required
        if( delete_attribs )
          AQE->remove_cubit_owner_attrib_in_BODY( master );

        BodySM *new_body = AME->get_new_Body( bodysm_ptr, BODY_ptr, master,
          keep_old_bodies );

        if( new_body )
          new_bodysm_list.append( new_body );
      }
      else
      {
        GfxPreview::clear();

        draw_tweak_preview_tagged_FACEs( master, CUBIT_TRUE );
        api_delent( master );
      }
    }
  }

  // Free memory
  while( BODY_EDGE_lists.size() ) delete BODY_EDGE_lists.pop();

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : tweak_target
// Member Type: PUBLIC
// Description: Tweak specified faces of a volume or volumes up to a set of
//              target surfaces.
// Author     : Steve Storm
// Date       : 10/27/05
//=============================================================================
CubitStatus AcisTweakTool::tweak_target( DLIList<Surface*> &input_surface_list,
                                         DLIList<Surface*> &target_surf_list,
                                         DLIList<BodySM*> &new_bodysm_list,
                                         CubitBoolean reverse_flg,
                                         CubitBoolean keep_old_body,
                                         CubitBoolean preview )
{
#ifndef ACIS_LOCAL_OPS
  PRINT_ERROR( "The ACIS Local Operations Husk is required for tweaking\n"
    "       surfaces to a target.  It has not been licensed for this installation.\n" );
  return CUBIT_FAILURE;
#endif

  if( target_surf_list.size() == 1 )
  {
    return tweak_target_single( input_surface_list, target_surf_list.get(),
      new_bodysm_list, reverse_flg, keep_old_body, preview );
  }
  else
  {
    return tweak_target_multiple( input_surface_list, target_surf_list,
      new_bodysm_list, reverse_flg, keep_old_body, preview );
  }
}

//=============================================================================
// Function   : tweak_target_single
// Member Type: PRIVATE
// Description: Tweak specified faces of a volume or volumes up to a single
//              target surface.
// Author     : Steve Storm
// Date       :
//=============================================================================
CubitStatus AcisTweakTool::tweak_target_single(
                                         DLIList<Surface*> &input_surface_list,
                                         Surface *target_surf_ptr,
                                         DLIList<BodySM*> &new_bodysm_list,
                                         CubitBoolean reverse_flg,
                                         CubitBoolean keep_old_body,
                                         CubitBoolean preview )
{
  int i;
  BODY* copied_BODY_ptr;
  BODY* BODY_ptr;
  BodySM* body_ptr;

  bool delete_attribs =
    (GeometryModifyTool::instance()->get_new_ids() || keep_old_body);

  // Copy the incoming input_surface_list since we will be pulling
  // surfaces out of it.
  DLIList<SurfaceACIS*> surface_list( input_surface_list.size() );
  if (!get_ACIS_surfaces( input_surface_list, surface_list ))
    return CUBIT_FAILURE;

  // Make sure that the none of the surfaces are sheet bodies (2-sided)
  for( i=surface_list.size(); i--; )
  {
    if( surface_list.get_and_step()->get_FACE_ptr()->sides() == DOUBLE_SIDED )
    {
      PRINT_ERROR("Cannot tweak surfaces that do not belong to a closed volume\n");
      return CUBIT_FAILURE;
    }
  }

  // Get target FACE
  Surface *target_surface_ptr = target_surf_ptr;
  SurfaceACIS *target_surf_acis = dynamic_cast<SurfaceACIS*>( target_surf_ptr );
  if (!target_surf_acis)
  {
    PRINT_ERROR("Target cannot be a non-ACIS Surface.\n");
    return CUBIT_FAILURE;
  }

  // Copy the target face (in some cases, the target face can be modified -
  // if for example it is tweaked to itself, which is sometimes done to
  // remove unwanted topology or to fix an error).  This preserves the target
  // in that case.
  ENTITY *copied_entity_ptr;
  outcome result = api_copy_entity_contents( target_surf_acis->get_FACE_ptr(),
    copied_entity_ptr );
  if( !result.ok() )
  {
    AQE->ACIS_API_error(result);
    return CUBIT_FAILURE;
  }
  FACE *target_FACE_ptr = (FACE *)copied_entity_ptr;

  // This acis api doesn't work exactly as I'd like.  The
  // only thing that makes sense using it, it seems, is to
  // tweak multiple surfaces to one surface.
  DLIList<FACE*> tweak_FACE_list;
  DLIList<SurfaceACIS*> tweak_surface_list;
  surface_list.reset();
  while( surface_list.size() )
  {
    // Remove a group of surfaces from surface_list that belong to the same
    // BODY.  Copy the BODY and put the copied surfaces in tweak_FACE_list.
    tweak_FACE_list.clean_out();
    tweak_surface_list.clean_out();
    if( AME->get_copied_FACES_of_body( surface_list, tweak_FACE_list,
      tweak_surface_list, copied_BODY_ptr ) == CUBIT_FAILURE )
    {
      api_delent( target_FACE_ptr );
      // Return success if any bodies were created
      return new_bodysm_list.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
    }

    // Store source FACEs for preview
    DLIList<AcisBridge*> preview_ab_list;
    if( preview )
    {
      get_owner_list( tweak_FACE_list, preview_ab_list );
      preview_ab_list.append_unique( 
        ATTRIB_CUBIT_OWNER::cubit_owner( target_surf_acis->get_FACE_ptr() ) );
    }

    // Do the tweak
    if( tweak_FACEs_to_target( tweak_FACE_list, target_FACE_ptr, reverse_flg )
      == CUBIT_FAILURE )
    {
      api_delent( copied_BODY_ptr );
      api_delent( target_FACE_ptr );
      // Return success if any bodies were created
      return new_bodysm_list.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
    }

    // Get original Body and BODY
    body_ptr = AQE->get_body_sm_of_ENTITY( copied_BODY_ptr );
    BODY_ptr = dynamic_cast<BodyACIS*>(body_ptr)->get_BODY_ptr();

    if( !preview )
    {
      // Now cleanout the owner attributes from the copied BODY, if required
      if( delete_attribs )
        AQE->remove_cubit_owner_attrib_in_BODY( copied_BODY_ptr );

      BodySM *new_body = AME->get_new_Body( body_ptr, BODY_ptr,
        copied_BODY_ptr, keep_old_body );

      if( new_body )
        new_bodysm_list.append( new_body );
    }
    else
    {
      GfxPreview::clear();

      // Note the preview_ab_list is used for a special case in the preview -
      // when the target is a FACE in the BODY itself.
      draw_tweak_preview_omt( copied_BODY_ptr, CUBIT_TRUE, &preview_ab_list );
      
      api_delent( copied_BODY_ptr );
    }
  }

  api_delent( target_FACE_ptr );

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : tweak_target_multiple
// Member Type: PRIVATE
// Description: Tweak specified faces of a volume or volumes up to a set of
//              target surfaces.
// Author     : Steve Storm
// Date       : 10/27/05
//=============================================================================
CubitStatus AcisTweakTool::tweak_target_multiple(
                                         DLIList<Surface*> &input_surface_list,
                                         DLIList<Surface*> &target_surf_list,
                                         DLIList<BodySM*> &new_bodysm_list,
                                         CubitBoolean reverse_flg,
                                         CubitBoolean keep_old_body,
                                         CubitBoolean preview )
{
  bool delete_attribs =
    (GeometryModifyTool::instance()->get_new_ids() || keep_old_body);

  // Copy the incoming input_surface_list since we will be pulling surfaces out
  // of it.
  DLIList<SurfaceACIS*> surface_list( input_surface_list.size() );
  if (!get_ACIS_surfaces( input_surface_list, surface_list ))
    return CUBIT_FAILURE;

  // Make sure that the none of the surfaces are sheet bodies (2-sided)
  int i;
  for( i=surface_list.size(); i--; )
  {
    if( surface_list.get_and_step()->get_FACE_ptr()->sides() == DOUBLE_SIDED )
    {
      PRINT_ERROR("Cannot tweak surfaces that do not belong to a solid volume\n");
      return CUBIT_FAILURE;
    }
  }

  // Note: target surfaces don't necessarily have to be from same BODY

  // Get target FACEs
  DLIList<FACE*> target_FACE_list;
  if( get_FACEs( target_surf_list, target_FACE_list ) == CUBIT_FAILURE )
  {
    PRINT_ERROR( "Targets must be ACIS surfaces\n" );
    return CUBIT_FAILURE;
  }

  // Get extended target surfaces (in a single sheet BODY)
  BODY *ext_target_BODY_ptr;
  if( create_extended_sheet( target_FACE_list, ext_target_BODY_ptr )
    == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  if( DEBUG_FLAG(168) )
  {
    BODY *tmp_BODY = AME->copy_BODY( ext_target_BODY_ptr, CUBIT_TRUE );
    BodySM *this_bodysm = AcisQueryEngine::instance()->
      populate_topology_bridges( tmp_BODY );
    new_bodysm_list.append( this_bodysm );
  }

  // Tweak the sources to the targets on a per body basis
  DLIList<FACE*> source_FACE_list;
  BODY *copied_BODY_ptr;
  BODY *BODY_ptr;
  BodySM *body_ptr;
  surface_list.reset();
  while( surface_list.size() )
  {
    // Remove a group of surfaces from surface_list that belong to the same
    // BODY. Copy the BODY and put the copied sources in source_FACE_list.
    source_FACE_list.clean_out();
    if( AME->get_copied_FACES_of_body( surface_list, source_FACE_list,
      copied_BODY_ptr ) == CUBIT_FAILURE )
      break;

    // Get original Body and BODY
    body_ptr = AQE->get_body_sm_of_ENTITY( copied_BODY_ptr );
    BODY_ptr = dynamic_cast<BodyACIS*>(body_ptr)->get_BODY_ptr();

    // Store source and target FACEs for preview (necessary especially if
    // tweaking surfaces on a body to surfaces on the same body)
    DLIList<AcisBridge*> preview_ab_list;
    if( preview )
    {
      get_owner_list( source_FACE_list, preview_ab_list );
      get_owner_list( target_FACE_list, preview_ab_list );
    }

    // Tweak the source FACEs up to the target sheet body
    if( tweak_target_multiple( source_FACE_list, target_FACE_list,
      ext_target_BODY_ptr, new_bodysm_list, reverse_flg ) == CUBIT_FAILURE )
    {
      api_delent( ext_target_BODY_ptr );
      api_delent( copied_BODY_ptr );
      // Return success if any bodies were created
      return new_bodysm_list.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
    }

    if( !preview )
    {
      // Now cleanout the owner attributes from the copied BODY, if required
      if( delete_attribs )
        AQE->remove_cubit_owner_attrib_in_BODY( copied_BODY_ptr );

      BodySM *new_body = AME->get_new_Body( body_ptr, BODY_ptr,
        copied_BODY_ptr, keep_old_body );

      if( new_body )
        new_bodysm_list.append( new_body );
    }
    else
    {
      GfxPreview::clear();

      // Draw preview EDGEs
      draw_tweak_preview_omt( copied_BODY_ptr, CUBIT_TRUE, &preview_ab_list );
      api_delent(copied_BODY_ptr);
    }
  }

  api_delent( ext_target_BODY_ptr );

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : custom_att_merge_method
// Member Type: Global
// Description: Used to mark EDGEs that are merged for tweak_target_multiple
//              method.  This is a custom merge_owner method for ACIS
//              generic attributes.
// Author     : Steve Storm
// Date       : 3/18/06
//=============================================================================
void
custom_att_merge_method( 
              ATTRIB_GEN_NAME *att,  // attribute whose owner is being merged 
              ENTITY *ENTITY_ptr2,   // other entity involved in merge 
              logical merge_method ) // attribute's owner will be lost?
{
  ENTITY *ENTITY_ptr1 = att->entity();

  if( IS_ENTITY_TYPE( ENTITY_ptr1, EDGE ) )
  {
    api_add_generic_named_attribute( ENTITY_ptr1, "tweak_merged", 1, SplitCopy );
    api_add_generic_named_attribute( ENTITY_ptr2, "tweak_merged", 1, SplitCopy );
  }

  return;
}

//=============================================================================
// Function   : tweak_target_multiple
// Member Type: PRIVATE
// Description: Tweak specified faces of a volume up to a sheet body
// Author     : Steve Storm
// Date       : 10/27/05
//=============================================================================
CubitStatus AcisTweakTool::tweak_target_multiple( DLIList<FACE*> &source_FACE_list,
                                                  DLIList<FACE*> &target_FACE_list,
                                                  BODY *ext_target_BODY_ptr,
                                                  DLIList<BodySM*> &debug_BodySM_list,
                                                  CubitBoolean reverse_flg )
{
  // Summary of method:
  //  1. Check for special case - extended target could contain only one
  //     surface - if so, use single tweak method.
  //  2. Assign attributes to side FACEs (FACEs alongside the sources) so that
  //     we can get back to these FACEs after the tweak.  Since we are cutting
  //     the body with an extended target, the cutting sheet could cross parts
  //     of the body we don't want to cut - if it crosses any FACEs without
  //     these attributes we throw away that portion of the cutting sheet.
  //  3. Copy the targets into a single sheet BODY.  We will make the surface
  //     normals consistent in this sheet BODY.  Our method relies on having
  //     consistent normals on the cutting sheet.  Since when extending the
  //     target surfaces, cylinders may become full cylinders (i.e., extend
  //     a planar FACE with one filleted surface adjoining - you will get a
  //     plane connected to a cylinder), we will ultimately have to remove the
  //     portions of the cylinder we don't want (remove the "back side" of the
  //     cylinder).  By computing the average normal of the original target
  //     surfaces, we will have a good idea as to which side of the extended
  //     cylinders to keep.
  //  4. Get weighted average FACE normal of source and target surfaces.
  //  5. Create a plane just past the targets aligned with the target normal.
  //     Tweak the source surfaces to this plane (the body is now tweaked past
  //     the targets).  This is the tweaked BODY.
  //  6. Make normals in the extended target sheet consistent with the normals
  //     in #3, since we will need to compare to the extended target sheet
  //     normals.
  //  7. Intersect the tweaked BODY (#5) and the extended target sheet - this
  //     makes a cutting sheet.  Note this cutting sheet may cross parts of the
  //     BODY where we don't want it to cut.
  //  8. Imprint so as to score the cutting sheet.
  //  9. Throw away individual cutting sheet surfaces that don't cross the
  //     FACEs we want to cut (the side FACEs tracked in #2).
  // 10. Since the extended target could contain complete cylinders, torii and
  //     spheres, AND it could cut through other portions of our body, we need
  //     to do some special handling before just cutting the body with it.  Our
  //     approach is to trim the extended target body to the desired cutting
  //     sheet.
  //     a) Remove periodic surfaces near nonmanifold edges.
  //     b) Align the normals of the cutting sheet with the extended target
  //        sheet - again, we need to keep the normals consistent to the
  //        original targets.
  //     c) Remove unwanted periodic FACEs on the back of the cutting sheet.
  //        Since we diligently kept the normals consistent, we can simply
  //        remove periodic surfaces with normals aligned with the negative
  //        of the target normal (from #4).
  // 11. Make sure cutting sheet normal is pointing in opposite direction of
  //     the source normal.  Note we can tweak forwards or backwards - no
  //     special handling is needed.
  // 12. Chopoff the body.

  // Check for special case - if ext_target_BODY_ptr only contains one
  // FACE, tweak to it and exit.
  DLIList<FACE*> ext_FACE_list;
  AQE->get_FACEs( ext_target_BODY_ptr, ext_FACE_list );
  if( ext_FACE_list.size() == 1 )
  {
    FACE *target_FACE_ptr = ext_FACE_list.get();
    if( tweak_FACEs_to_target( source_FACE_list, target_FACE_ptr, reverse_flg )
      == CUBIT_FAILURE )
      return CUBIT_FAILURE;

    return CUBIT_SUCCESS;
  }

  // Get the BODY the sources are in
  FACE *FACE_ptr = source_FACE_list.get();
  BODY *copied_BODY_ptr = AQE->get_BODY_of_ENTITY( FACE_ptr );

  // This operation sometimes destroys owner attributes
  AcisBridge *ab_body_ptr;
  AcisBridge *ab_lump_ptr;
  AcisBridge *ab_shell_ptr;
  get_owner_attribs( copied_BODY_ptr, ab_body_ptr, ab_lump_ptr, ab_shell_ptr );

  // Assign attributes to side FACEs for later use
  int i, j, k;
  outcome result;
  ATTRIB_GEN_NAME *ret_att;
  for( i=source_FACE_list.size(); i--; )
  {
    FACE_ptr = source_FACE_list.get_and_step();

    DLIList<EDGE*> EDGE_list;
    AQE->get_EDGEs( FACE_ptr, EDGE_list );

    for( j=EDGE_list.size(); j--; )
    {
      EDGE *EDGE_ptr = EDGE_list.get_and_step();
      DLIList<FACE*> attached_FACEs;
      AQE->get_FACEs( EDGE_ptr, attached_FACEs );

      for( k=attached_FACEs.size(); k--; )
      {
        FACE_ptr = attached_FACEs.get_and_step();
        result = api_find_named_attribute( FACE_ptr, "tweak", ret_att );
        if( !ret_att )
        {
          result = api_add_generic_named_attribute( FACE_ptr, "tweak", 1, SplitCopy );
          if( !result.ok() )
          {
            AQE->ACIS_API_error(result);
            return CUBIT_FAILURE;
          }
        }
      }
    }
  }

  // Remove named attributes from the sources
  for( i=source_FACE_list.size(); i--; )
  {
    FACE_ptr = source_FACE_list.get_and_step();
    api_remove_generic_named_attribute(FACE_ptr, "tweak");
  }

  // Store the sources in a sheet - needed for distance calculation and may be
  // needed later for self-intersecting BODY issue.
  BODY *source_sheet;
  if( copy_FACEs_into_sheet( source_FACE_list, source_sheet ) == CUBIT_FAILURE )
  {
    PRINT_ERROR( "Unable to tweak to targets\n" );
    remove_named_attribs( copied_BODY_ptr, "tweak" );
    return CUBIT_FAILURE;
  }

  // Make a target sheet with copies of original target FACEs in it.  We will
  // extract the target normal from this BODY.
  BODY *target_sheet;
  if( copy_FACEs_into_sheet( target_FACE_list, target_sheet ) == CUBIT_FAILURE )
  {
    api_delent( source_sheet );
    remove_named_attribs( copied_BODY_ptr, "tweak" );
    return CUBIT_FAILURE;
  }

  // Make the surface normals consistent.  Note user could pick a variety of
  // target surfaces from different bodies, without consistent surface normals.
  if( make_surf_normals_consistent( target_sheet )== CUBIT_FAILURE )
  {
    PRINT_ERROR( "Unable to make target surface normals consistent - invalid target\n" );
    api_delent( source_sheet );
    api_delent( target_sheet );
    remove_named_attribs( copied_BODY_ptr, "tweak" );
    return CUBIT_FAILURE;
  }

  // FYI: Property of dot product:
  //   * angle between vectors acute (<90) if dot product > 0
  //   * angle between vectors obtuse (>90) if dot product < 0
  //   * angle between vectors 90 deg if dot product = 0

  // Get average direction of source normals.
  CubitVector source_norm = weighted_average_FACE_normal( source_FACE_list );
  if( source_norm.length() < 1e-12 )
  {
    PRINT_ERROR( "Unable to tweak source surfaces that form a loop\n" );
    api_delent( source_sheet );
    api_delent( target_sheet );
    remove_named_attribs( copied_BODY_ptr, "tweak" );
    return CUBIT_FAILURE;
  }

  // Get average target normal
  DLIList<FACE*> target_sheet_FACE_list;
  AQE->get_FACEs( target_sheet, target_sheet_FACE_list );
  CubitVector target_norm = weighted_average_FACE_normal( target_sheet_FACE_list );
  if( target_norm.length() < 1e-12 )
  {
    PRINT_ERROR( "Unable to tweak to target surfaces that form a loop\n" );
    api_delent( source_sheet );
    api_delent( target_sheet );
    remove_named_attribs( copied_BODY_ptr, "tweak" );
    return CUBIT_FAILURE;
  }

  // Get a plane through one of the source's vertices (with normal aligned with
  // target normal (but pointing in direction of source normal)
  source_FACE_list.reset();
  DLIList<VERTEX*> temp_VERTEX_list;
  AQE->get_VERTICEs( source_FACE_list.get(), temp_VERTEX_list );
  VERTEX *VERTEX_ptr = temp_VERTEX_list.get();
  SPAposition acis_coords = (VERTEX_ptr->geometry())->coords();
  CubitVector coords( acis_coords.x(), acis_coords.y(), acis_coords.z() );
  CubitPlane source_plane( target_norm, coords );

  // Have plane's normal point in the general direction of source_norm
  if( source_norm % source_plane.normal() < 0.0 ) // Dot product
  {
    // Obtuse angle - reverse the plane
    CubitVector neg_norm = -source_plane.normal();
    source_plane.normal( neg_norm );
    source_plane.coefficient( -source_plane.coefficient() );
  }

  // Offset the plane past the targets.  Do this in two steps - in some
  // cases, our plane is far from the sources and this can cause problems
  // (if there is a cylindrical surface, the API error "no solution for an
  // edge" has occurred).

  // Find offset distance for first tweak - just past the sources (this 
  // will most likely make one planar source)
  double offset_dist;
  if( extrema_pln_BODY( source_plane, source_sheet, offset_dist )
    == CUBIT_FAILURE )
  {
    api_delent( source_sheet );
    api_delent( target_sheet );
    remove_named_attribs( copied_BODY_ptr, "tweak" );
    return CUBIT_FAILURE;
  }

  BODY *planar_target_BODY;
  if( create_offset_planar_body( source_plane, offset_dist,
    planar_target_BODY ) == CUBIT_FAILURE )
  {
    api_delent( source_sheet );
    api_delent( target_sheet );
    remove_named_attribs( copied_BODY_ptr, "tweak" );
    return CUBIT_FAILURE;
  }

  if( DEBUG_FLAG(168) )
  {
    BODY *tmp_BODY = AME->copy_BODY( planar_target_BODY, CUBIT_TRUE );
    BodySM *this_bodysm = AQE->populate_topology_bridges( tmp_BODY );
    debug_BodySM_list.append( this_bodysm );
  }

  DLIList<FACE*> planar_FACE_list;;
  AQE->get_FACEs( planar_target_BODY, planar_FACE_list );
  FACE *planar_target_FACE = planar_FACE_list.get();

  // Mark the sources and other FACEs
  for( i=source_FACE_list.size(); i--; )
  {
    FACE_ptr = source_FACE_list.get_and_step();
    api_add_generic_named_attribute( FACE_ptr, "tweak_source", 1, SplitCopy );
  }
  DLIList<FACE*> copied_BODY_FACE_list;
  AQE->get_FACEs( copied_BODY_ptr, copied_BODY_FACE_list );
  for( i=copied_BODY_FACE_list.size(); i--; )
  {
    FACE_ptr = copied_BODY_FACE_list.get_and_step();
    api_find_named_attribute( FACE_ptr, "tweak_source", ret_att );
    if( !ret_att )
      api_add_generic_named_attribute( FACE_ptr, "tweak_other", 1 );
  }

  // Tweak the sources to the first planar target
  CubitBoolean skip_self_int_check = CUBIT_TRUE;

  if( tweak_FACEs_to_target( source_FACE_list, planar_target_FACE, 
                             CUBIT_FALSE, skip_self_int_check ) == CUBIT_FAILURE )
  {
    api_delent( source_sheet );
    api_delent( target_sheet );
    api_delent( planar_target_BODY );
    remove_named_attribs( copied_BODY_ptr, "tweak" );
    return CUBIT_FAILURE;
  }

  // Delete the first planar target
  api_delent( planar_target_BODY );

  // Get new source list
  source_FACE_list.clean_out();
  copied_BODY_FACE_list.clean_out();
  AQE->get_FACEs( copied_BODY_ptr, copied_BODY_FACE_list );
  for( i=copied_BODY_FACE_list.size(); i--; )
  {
    FACE_ptr = copied_BODY_FACE_list.get_and_step();
    api_find_named_attribute( FACE_ptr, "tweak_source", ret_att );
    if( ret_att )
    {
      source_FACE_list.append( FACE_ptr );
      api_remove_generic_named_attribute(FACE_ptr, "tweak_source");
    }
    else
    {
      api_find_named_attribute( FACE_ptr, "tweak_other", ret_att );
      if( ret_att )
        api_remove_generic_named_attribute(FACE_ptr, "tweak_other");
      else
      {
        // No attribute - this must be the new source
        source_FACE_list.append( FACE_ptr );
      }
    }
  }

  // Get a plane for the second tweak
  target_FACE_list.reset();
  temp_VERTEX_list.clean_out();
  AQE->get_VERTICEs( target_FACE_list.get(), temp_VERTEX_list );
  VERTEX_ptr = temp_VERTEX_list.get();
  acis_coords = (VERTEX_ptr->geometry())->coords();
  coords.set( acis_coords.x(), acis_coords.y(), acis_coords.z() );
  CubitPlane target_plane( target_norm, coords );

  // Have plane's normal point in the general direction of source_norm
  if( source_norm % target_plane.normal() < 0.0 ) // Dot product
  {
    // Obtuse angle - reverse the plane
    CubitVector neg_norm = -target_plane.normal();
    target_plane.normal( neg_norm );
    target_plane.coefficient( -target_plane.coefficient() );
  }

  // Find offset distance for second tweak - past the targets.  Calculate the
  // offset distance as the maximum distance from the plane to any point on the
  // extended target body, in the direction of the plane's normal.
  if( extrema_pln_BODY( target_plane, ext_target_BODY_ptr, offset_dist )
    == CUBIT_FAILURE )
  {
    api_delent( source_sheet );
    api_delent( target_sheet );
    remove_named_attribs( copied_BODY_ptr, "tweak" );
    return CUBIT_FAILURE;
  }

  // Create a target FACE offset from the plane (arbitrarily 10.0 X 10.0
  // in size)
  if( create_offset_planar_body( target_plane, offset_dist,
    planar_target_BODY ) == CUBIT_FAILURE )
  {
    api_delent( source_sheet );
    api_delent( target_sheet );
    remove_named_attribs( copied_BODY_ptr, "tweak" );
    return CUBIT_FAILURE;
  }

  if( DEBUG_FLAG(168) )
  {
    BODY *tmp_BODY = AME->copy_BODY( planar_target_BODY, CUBIT_TRUE );
    BodySM *this_bodysm = AQE->populate_topology_bridges( tmp_BODY );
    debug_BodySM_list.append( this_bodysm );
  }

  planar_FACE_list.clean_out();
  AQE->get_FACEs( planar_target_BODY, planar_FACE_list );
  planar_target_FACE = planar_FACE_list.get();

  // Tweak the sources to the planar target
  if( tweak_FACEs_to_target( source_FACE_list, planar_target_FACE, 
                             CUBIT_FALSE, skip_self_int_check ) == CUBIT_FAILURE )
  {
    api_delent( source_sheet );
    api_delent( target_sheet );
    api_delent( planar_target_BODY );
    remove_named_attribs( copied_BODY_ptr, "tweak" );
    return CUBIT_FAILURE;
  }

  // Free the plane since we are done with it
  api_delent( planar_target_BODY );

  // Copy the ext_target_BODY_ptr so as not to modify it
  BODY *ext_target_BODY_ptr2 = AME->copy_BODY( ext_target_BODY_ptr, CUBIT_TRUE );

  // Heal the body for good measure
  heal_BODY( ext_target_BODY_ptr2 );

  // Make the shell normals in the ext_target_BODY_ptr2 consistent with the
  // target sheet.
  FACE *seed_FACE_ptr = get_seed_FACE( target_sheet_FACE_list );

  DLIList<FACE*> ext_target_BODY_FACE_list;
  AQE->get_FACEs( ext_target_BODY_ptr2, ext_target_BODY_FACE_list );

  FACE *overlap_FACE_ptr = find_overlap_FACE( seed_FACE_ptr,
    ext_target_BODY_FACE_list );
  if( !overlap_FACE_ptr )
  {
    PRINT_ERROR( "Problem creating extended target to tweak to\n" );
    if( DEBUG_FLAG(168) )
    {
      BODY *tmp_BODY = AME->copy_BODY( ext_target_BODY_ptr2, CUBIT_TRUE );
      BodySM *this_bodysm = AQE->populate_topology_bridges( tmp_BODY );
      debug_BodySM_list.append( this_bodysm );
    }
    api_delent( source_sheet );
    api_delent( target_sheet );
    api_delent( ext_target_BODY_ptr2 );
    remove_named_attribs( copied_BODY_ptr, "tweak" );
    return CUBIT_FAILURE;
  }

  CubitVector norm_loc;
  CubitVector ext_norm = FACE_normal( overlap_FACE_ptr, norm_loc );
  CubitVector seed_norm = FACE_normal( seed_FACE_ptr, norm_loc, CUBIT_FALSE );

  if( ext_norm % seed_norm < 0.0 )
    api_reverse_face( overlap_FACE_ptr );

  if( make_surf_normals_consistent( ext_target_BODY_ptr2, overlap_FACE_ptr )
    == CUBIT_FAILURE )
  {
    PRINT_ERROR( "Unable to make extended target shell normals consistent - invalid target\n" );
    if( DEBUG_FLAG(168) )
    {
      BODY *tmp_BODY = AME->copy_BODY( ext_target_BODY_ptr2, CUBIT_TRUE );
      BodySM *this_bodysm = AQE->populate_topology_bridges( tmp_BODY );
      debug_BodySM_list.append( this_bodysm );
    }
    api_delent( source_sheet );
    api_delent( target_sheet );
    api_delent( ext_target_BODY_ptr2 );
    remove_named_attribs( copied_BODY_ptr, "tweak" );
    return CUBIT_FAILURE;
  }

  if( DEBUG_FLAG(168) )
  {
    BODY *tmp_BODY = AME->copy_BODY( ext_target_BODY_ptr2, CUBIT_TRUE );
    BodySM *this_bodysm = AQE->populate_topology_bridges( tmp_BODY );
    debug_BodySM_list.append( this_bodysm );
  }

  // We are done with the target_sheet
  api_delent( target_sheet );

  // Intersect the tweaked body and the ext_target_BODY_ptr2.  We need to then
  // determine which loops on the resulting body to keep as the cutting
  // sheet.  Note a NONREG_INTERSECTION sometimes corrupts the cutting_sheet,
  // so we do a regularized intersection then imprint afterwards.
  BODY *cutting_sheet = NULL;
  result = api_boolean( copied_BODY_ptr, ext_target_BODY_ptr2, INTERSECTION,
    NDBOOL_KEEP_BOTH, cutting_sheet );
  if( !result.ok() )
  {
    // A likely problem is that the BODY intersects with itself - we need to
    // handle this case.  To fix we will chop off the tweaked part of the BODY
    // using the original sources; then our booleans will work.  We then have
    // to reassemble the BODY.
    
    BODY *copied_BODY_ptr2 = AME->copy_BODY( copied_BODY_ptr, CUBIT_FALSE );

    // Chopoff the copied_BODY_ptr with cutting_sheet
    if( chop_off_with_sheet( copied_BODY_ptr2, source_sheet ) == CUBIT_FAILURE )
    {
      AQE->ACIS_API_error(result);
      api_delent( source_sheet );
      api_delent( ext_target_BODY_ptr2 );
      api_delent( copied_BODY_ptr2 );
      remove_named_attribs( copied_BODY_ptr, "tweak" );
      return CUBIT_FAILURE;
    }

    // Chopoff the copied_BODY_ptr with cutting_sheet the other way
    api_reverse_body( source_sheet );
    if( chop_off_with_sheet( copied_BODY_ptr, source_sheet ) == CUBIT_FAILURE )
    {
      AQE->ACIS_API_error(result);
      api_delent( source_sheet );
      api_delent( ext_target_BODY_ptr2 );
      api_delent( copied_BODY_ptr2 );
      remove_named_attribs( copied_BODY_ptr, "tweak" );
      return CUBIT_FAILURE;
    }

    // We are done with the source_sheet
    api_delent( source_sheet );

    // Note copied_BODY_ptr2 is the cutoff tweaked piece - perform steps on
    // this piece.  The good thing is we won't have to worry about looking for
    // FACEs crossed by the extended target that we don't want to affect.

    result = api_boolean( copied_BODY_ptr2, ext_target_BODY_ptr2, INTERSECTION,
      NDBOOL_KEEP_BOTH, cutting_sheet );
    if( !result.ok() )
    {
      api_delent( ext_target_BODY_ptr2 );
      api_delent( copied_BODY_ptr2 );
      remove_named_attribs( copied_BODY_ptr, "tweak" );
      return CUBIT_FAILURE;
    }

    // Prepare the cutting sheet
    if( prep_cutting_sheet( cutting_sheet, copied_BODY_ptr2, ext_target_BODY_ptr2,
      source_norm, target_norm, CUBIT_FALSE ) == CUBIT_FAILURE )
    {
      if( DEBUG_FLAG(168) )
      {
        BodySM *this_bodysm = AQE->populate_topology_bridges( cutting_sheet );
        debug_BodySM_list.append( this_bodysm );
      }
      else
        api_delent( cutting_sheet );
      api_delent( ext_target_BODY_ptr2 );
      api_delent( copied_BODY_ptr2 );
      remove_named_attribs( copied_BODY_ptr, "tweak" );
      return CUBIT_FAILURE;
    }

    // We are done with the "tweak" attributes (used in prep_cutting_sheet)
    remove_named_attribs( copied_BODY_ptr, "tweak" );

    // We are now done with the copied extended target sheet
    api_delent( ext_target_BODY_ptr2 );

    // Reverse if forced to
    if( reverse_flg )
      api_reverse_body( cutting_sheet );

    if( DEBUG_FLAG(168) )
    {
      BODY *tmp_BODY = AME->copy_BODY( cutting_sheet, CUBIT_TRUE );
      BodySM *this_bodysm = AQE->populate_topology_bridges( tmp_BODY );
      debug_BodySM_list.append( this_bodysm );
    }

    // Apply a "tweak2" named attribute as part of a multi-step process we use
    // to determine which EDGEs to regularize out when we reassemble the BODY.  
    // The "tweak2" attribute will not exist on the EDGEs that result from the
    // chopoff operation with the cutting sheet.  We will use a custom function
    // that gets called during boolean merge operations.
    DLIList<EDGE*> BODY2_EDGE_list;
    AQE->get_EDGEs( copied_BODY_ptr2, BODY2_EDGE_list );
    EDGE *EDGE_ptr;
    for( i=BODY2_EDGE_list.size(); i--; )
    {
      EDGE_ptr = BODY2_EDGE_list.get_and_step();
      api_add_generic_named_attribute( EDGE_ptr, "tweak2", 1, SplitLose,
        MergeCustom );
    }

    // Chopoff the copied_BODY_ptr2 with cutting_sheet
    if( chop_off_with_sheet( copied_BODY_ptr2, cutting_sheet ) == CUBIT_FAILURE )
    {
      api_delent( copied_BODY_ptr2 );
      api_delent( cutting_sheet );
      return CUBIT_FAILURE;
    }

    // We are done with the cutting sheet
    api_delent( cutting_sheet );

    // This will add a "tweak_merged" attribute to any EDGE containing a
    // "tweak2" attribute that is merged with another EDGE during the following
    // unite.  We can then clean these EDGEs out.  Note the EDGEs that were
    // sliced at the cutting sheet location won't have "tweak2" attributes on
    // them - thus this method allows us to efficiently find the EDGEs we need
    // to remove.
    set_merge_method( "tweak2", custom_att_merge_method );

    // Unite the copied_BODY_ptr2 with copied_BODY_ptr (this deletes
    // copied_BODY_ptr2)
    result = api_boolean( copied_BODY_ptr2, copied_BODY_ptr, NONREG_UNION );
    if( !result.ok() )
    {
      if( copied_BODY_ptr2 ) api_delent( copied_BODY_ptr2 );
      return CUBIT_FAILURE;
    }

    // The unite can leave double-sided FACEs within the BODY - we need to
    // to remove those.
    DLIList<FACE*>  body_FACE_list;
    AQE->get_FACEs( copied_BODY_ptr, body_FACE_list );
    for( i=body_FACE_list.size(); i--; )
    {
      FACE_ptr = body_FACE_list.get_and_step();

      // If this is an internal FACE remove it
      if( FACE_ptr->sides() == DOUBLE_SIDED &&
          FACE_ptr->cont()==BOTH_INSIDE )
      {
        BODY *unhooked_BODY_ptr = NULL;
        api_unhook_face( FACE_ptr, unhooked_BODY_ptr );
        if( unhooked_BODY_ptr )
          api_delent( unhooked_BODY_ptr );
      }
    }

    // Remove the merged EDGEs
    DLIList<EDGE*> merged_EDGE_list;
    DLIList<EDGE*> body_EDGE_list;
    AQE->get_EDGEs( copied_BODY_ptr, body_EDGE_list );
    for( i=body_EDGE_list.size(); i--; )
    {
      EDGE_ptr = body_EDGE_list.get_and_step();
      api_find_named_attribute( EDGE_ptr, "tweak_merged", ret_att );
      if( ret_att )
      {
        merged_EDGE_list.append( EDGE_ptr );
        api_remove_generic_named_attribute(EDGE_ptr, "tweak_merged");
      }
    }

    for( i=merged_EDGE_list.size(); i--; )
    {
      EDGE_ptr = merged_EDGE_list.get_and_step();

      // Make sure EDGE was not already cleaned out
      DLIList<EDGE*> tmp_EDGE_list;
      AQE->get_EDGEs( copied_BODY_ptr, tmp_EDGE_list );

      if( tmp_EDGE_list.is_in_list( EDGE_ptr ) )
        api_clean_entity( EDGE_ptr );
    }

    // Cleanout any left-over named attributes
    remove_named_attribs( copied_BODY_ptr, "tweak2" );

    // Reset the owner attributes
    reset_owner_attribs( copied_BODY_ptr, ab_body_ptr, ab_lump_ptr, ab_shell_ptr );

    return CUBIT_SUCCESS;
  } // End of workaround for self-intersecting case

  // Since the intersection was successful the source_sheet is not needed
  api_delent( source_sheet );

  // Prepare the cutting sheet
  if( prep_cutting_sheet( cutting_sheet, copied_BODY_ptr, ext_target_BODY_ptr2,
    source_norm, target_norm ) == CUBIT_FAILURE )
  {
    if( DEBUG_FLAG(168) )
    {
      BodySM *this_bodysm = AQE->populate_topology_bridges( cutting_sheet );
      debug_BodySM_list.append( this_bodysm );
    }
    else
      api_delent( cutting_sheet );
    api_delent( ext_target_BODY_ptr2 );
    remove_named_attribs( copied_BODY_ptr, "tweak" );
    return CUBIT_FAILURE;
  }

  // We are done with the "tweak" attributes (used in prep_cutting_sheet)
  remove_named_attribs( copied_BODY_ptr, "tweak" );

  // We are now done with the copied extended target sheet
  api_delent( ext_target_BODY_ptr2 );

  // Reverse if forced to
  if( reverse_flg )
    api_reverse_body( cutting_sheet );

  if( DEBUG_FLAG(168) )
  {
    BODY *tmp_BODY = AME->copy_BODY( cutting_sheet, CUBIT_TRUE );
    BodySM *this_bodysm = AQE->populate_topology_bridges( tmp_BODY );
    debug_BodySM_list.append( this_bodysm );
  }

  // Chopoff the copied_BODY_ptr with cutting_sheet
  if( chop_off_with_sheet( copied_BODY_ptr, cutting_sheet ) == CUBIT_FAILURE )
  {
    api_delent( cutting_sheet );
    return CUBIT_FAILURE;
  }

  // Delete the cutting sheet
  api_delent( cutting_sheet );

  // Reset the owner attributes
  reset_owner_attribs( copied_BODY_ptr, ab_body_ptr, ab_lump_ptr, ab_shell_ptr );

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : FACE_surrounded
// Member Type: PRIVATE
// Description: Determine if the given FACE is entirely surrounded by the
//              other FACEs. The given FACE can exist in FACE_list_in.
// Author     : Steve Storm
// Date       : 01/05/06
//=============================================================================
CubitBoolean
AcisTweakTool::FACE_surrounded( FACE *ref_FACE_ptr,
                                DLIList<FACE*> &FACE_list_in )
{
  // Get outer loop curves in ref_FACE_ptr
  DLIList<EDGE*> ref_FACE_edges;
  get_outer_EDGEs( ref_FACE_ptr, ref_FACE_edges );

  // Loop through ref_FACE_edges, get all attached FACEs, check if any are in
  // FACE_list_in.

  int i, j;
  int num_found = 0;
  EDGE *ref_EDGE_ptr;
  for( i=ref_FACE_edges.size(); i--; )
  {
    ref_EDGE_ptr = ref_FACE_edges.get_and_step();

    DLIList<FACE *> attached_FACE_list;
    AQE->get_FACEs( ref_EDGE_ptr, attached_FACE_list );

    if( attached_FACE_list.move_to( ref_FACE_ptr ) == CUBIT_TRUE )
      attached_FACE_list.extract();

    FACE *FACE_ptr;
    for( j=attached_FACE_list.size(); j--; )
    {
      FACE_ptr = attached_FACE_list.get_and_step();
      if( FACE_list_in.is_in_list( FACE_ptr ) == CUBIT_TRUE )
      {
        num_found++;
        break;
      }
    }
  }

  if( num_found == ref_FACE_edges.size() )
    return CUBIT_TRUE;

  return CUBIT_FALSE;
}

//=============================================================================
// Function   : get_outer_EDGEs
// Member Type: PRIVATE
// Description: Get EDGEs from outer LOOP of FACE
// Author     : Steve Storm
// Date       : 01/03/06
//=============================================================================
CubitStatus
AcisTweakTool::get_outer_EDGEs( FACE *FACE_ptr, DLIList<EDGE*> &EDGE_list )
{
  // Loop through LOOPs
  outcome result;
  logical external;

  LOOP *LOOP_ptr = FACE_ptr->loop();
  while( LOOP_ptr != NULL )
  {
    result = api_loop_external( LOOP_ptr, &external );
    if( !result.ok() )
    {
      AQE->ACIS_API_error(result);
      return CUBIT_FAILURE;
    }

    if( external == TRUE )
    {
      AQE->get_EDGEs( LOOP_ptr, EDGE_list );
      return CUBIT_SUCCESS;
    }

    LOOP_ptr = LOOP_ptr->next();
  }

  return CUBIT_FAILURE;
}

//=============================================================================
// Function   : tweak_target
// Member Type: PUBLIC
// Description: Tweak specified edges of a surface or set of surfaces (in sheet
//              bodies) up to a single target surface or a set of target
//              surfaces.
// Author     : Steve Storm
// Date       : 03/28/05
//=============================================================================
CubitStatus
AcisTweakTool::tweak_target( DLIList<Curve*> &input_curve_list,
                             DLIList<Surface*> &target_surf_list,
                             DLIList<BodySM*> &new_bodysm_list,
                             CubitBoolean reverse_flg,
                             CubitBoolean keep_old_bodies,
                             CubitBoolean preview )
{
#ifndef ACIS_LOCAL_OPS
  PRINT_ERROR( "The ACIS Local Operations Husk is required for tweaking\n"
    "       surfaces to a target.  It has not been licensed for this installation.\n" );
  return CUBIT_FAILURE;
#endif

  if( target_surf_list.size() == 1 )
  {
    return tweak_target_single( input_curve_list, target_surf_list.get(),
      new_bodysm_list, reverse_flg, keep_old_bodies, preview );
  }
  else
  {
    return tweak_target_multiple( input_curve_list, target_surf_list,
      new_bodysm_list, reverse_flg, keep_old_bodies, preview );
  }
}

//=============================================================================
// Function   : tweak_target_single
// Member Type: PRIVATE
// Description: Tweak specified edges of a surface or set of surfaces (in sheet
//              bodies) up to a single target surface.
// Author     : Steve Storm
// Date       : 03/28/05
//=============================================================================
CubitStatus
AcisTweakTool::tweak_target_single( DLIList<Curve*> &input_curve_list,
                                    Surface *target_surf_ptr,
                                    DLIList<BodySM*> &new_bodysm_list,
                                    CubitBoolean reverse_flg,
                                    CubitBoolean keep_old_bodies,
                                    CubitBoolean preview )
{
#ifndef ACIS_LOCAL_OPS
  PRINT_ERROR( "The ACIS Local Operations Husk is required for tweaking\n"
    "       curves.  It has not been licensed for this installation.\n" );
  return CUBIT_FAILURE;
#endif

  // Get target FACE
  SurfaceACIS* target_surface = dynamic_cast<SurfaceACIS*>(target_surf_ptr);
  if( !target_surface )
  {
    PRINT_ERROR("Cannot tweak to non-ACIS target surface.\n");
    return CUBIT_FAILURE;
  }

  // Copy the target face (in some cases, the target face can be modified -
  // if for example it is tweaked to itself, which is sometimes done to
  // remove unwanted topology or to fix an error).  This preserves the target
  // in that case.
  ENTITY *copied_entity_ptr;
  outcome result = api_copy_entity_contents( target_surface->get_FACE_ptr(),
    copied_entity_ptr );
  if( !result.ok() )
  {
    AcisQueryEngine::instance()->ACIS_API_error(result);
    return CUBIT_FAILURE;
  }
  FACE *target_FACE_ptr = (FACE *)copied_entity_ptr;

  // Get source EDGEs
  DLIList<EDGE*> input_EDGE_list( input_curve_list.size() );
  if( get_EDGEs( input_curve_list, input_EDGE_list ) == CUBIT_FAILURE )
  {
    api_delent( target_FACE_ptr );
    return CUBIT_FAILURE;
  }

  // Call private function to do the work
  if( tweak_target_single( input_EDGE_list, target_FACE_ptr, new_bodysm_list,
    reverse_flg, keep_old_bodies, preview ) == CUBIT_FAILURE )
  {
    api_delent( target_FACE_ptr );
    return CUBIT_FAILURE;
  }

  api_delent( target_FACE_ptr );

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : tweak_target_multiple
// Member Type: PRIVATE
// Description: Tweak specified edges of a surface or set of surfaces (in sheet
//              bodies) up to a set of target surfaces.
// Author     : Steve Storm
// Date       : 11/03/05
//=============================================================================
CubitStatus
AcisTweakTool::tweak_target_multiple( DLIList<Curve*> &input_curve_list,
                                      DLIList<Surface*> &target_surf_list,
                                      DLIList<BodySM*> &new_bodysm_list,
                                      CubitBoolean reverse_flg,
                                      CubitBoolean keep_old_bodies,
                                      CubitBoolean preview )
{
#ifndef ACIS_LOCAL_OPS
  PRINT_ERROR( "The ACIS Local Operations Husk is required for extending\n"
    "       surfaces.  It has not been licensed for this installation.\n" );
  return CUBIT_FAILURE;
#endif

  // Get source curves
  DLIList<EDGE*> input_EDGE_list( input_curve_list.size() );
  if( get_EDGEs( input_curve_list, input_EDGE_list ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  // Get target FACES
  DLIList<FACE*> target_FACE_list( target_surf_list.size() );
  if( get_FACEs( target_surf_list, target_FACE_list ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  // Call private function to do the work
  if( tweak_target_multiple( input_EDGE_list, target_FACE_list,
    new_bodysm_list, reverse_flg, keep_old_bodies, preview ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : tweak_target
// Member Type: PUBLIC
// Description: Tweak specified edges of a sheet body or bodies up to a set of
//              target curves that are part of a sheet body.  The target is a
//              surface or surfaces created by thickening the owning surfaces
//              of the target curves.
// Author     : Steve Storm
// Date       : 03/28/05
//=============================================================================
CubitStatus
AcisTweakTool::tweak_target( DLIList<Curve*> &input_curve_list,
                             DLIList<Curve*> &target_curve_list,
                             DLIList<BodySM*> &new_bodysm_list,
                             CubitBoolean reverse_flg,
                             CubitBoolean keep_old_bodies,
                             CubitBoolean preview )
{
#ifndef ACIS_LOCAL_OPS
  PRINT_ERROR( "The ACIS Local Operations Husk is required for extending\n"
    "       surfaces.  It has not been licensed for this installation.\n" );
  return CUBIT_FAILURE;
#endif

  if( target_curve_list.size() == 1 )
  {
    return tweak_target_single( input_curve_list, target_curve_list.get(),
      new_bodysm_list, reverse_flg, keep_old_bodies, preview );
  }
  else
  {
    return tweak_target_multiple( input_curve_list, target_curve_list,
      new_bodysm_list, reverse_flg, keep_old_bodies, preview );
  }
}

//=============================================================================
// Function   : tweak_target_single
// Member Type: PRIVATE
// Description: Tweak specified edges of a sheet body or bodies up to a single
//              target curve that is part of a sheet body.  The target is a
//              surface created by thickening the owning surface of the target
//              curve.
// Author     : Steve Storm
// Date       : 03/28/05
//=============================================================================
CubitStatus
AcisTweakTool::tweak_target_single( DLIList<Curve*> &curve_list,
                                    Curve *target_curve_ptr,
                                    DLIList<BodySM*> &new_bodysm_list,
                                    CubitBoolean reverse_flg,
                                    CubitBoolean keep_old_bodies,
                                    CubitBoolean preview )
{
#ifndef ACIS_LOCAL_OPS
  PRINT_ERROR( "The ACIS Local Operations Husk is required for extending\n"
    "       surfaces.  It has not been licensed for this installation.\n" );
  return CUBIT_FAILURE;
#endif

  // We need to build the target FACE by grabbing the associated surfaces and
  // thickening them.  For convenience, use the get_thickened_BODIES_of_EDGES
  // function to do this.
  EDGE *target_EDGE_ptr = AcisQueryEngine::instance()->
    get_EDGE( target_curve_ptr );
  if( target_EDGE_ptr == NULL )
  {
    PRINT_ERROR( "Target must be an ACIS curve.\n" );
    return CUBIT_FAILURE;
  }
  BODY *BODY_ptr;
  DLIList<EDGE*> removed_EDGE_list;
  DLIList<EDGE*> target_EDGE_list( 1 );
  target_EDGE_list.append( target_EDGE_ptr );
  DLIList<BODY*> thickened_BODY_list;
  DLIList<DLIList<EDGE*>*> output_EDGE_lists;
  DLIList<DLIList<FACE*>*> output_FACE_lists;
  DLIList<DLIList<FACE*>*> conjugate_FACE_lists;
  if( get_thickened_BODIES_of_EDGES( "tweak to", target_EDGE_list,
    removed_EDGE_list, BODY_ptr, thickened_BODY_list, output_EDGE_lists,
    output_FACE_lists, conjugate_FACE_lists ) == CUBIT_FAILURE )
  {
    PRINT_ERROR( "unable to setup target -- aborting\n" );
    return CUBIT_FAILURE;
  }

  FACE *target_FACE_ptr = conjugate_FACE_lists.get()->get();

  // Clean up memory allocated for the output lists
  while( output_EDGE_lists.size() ) delete output_EDGE_lists.pop();
  while( output_FACE_lists.size() ) delete output_FACE_lists.pop();
  while( conjugate_FACE_lists.size() ) delete conjugate_FACE_lists.pop();

  if( target_FACE_ptr == NULL )
  {
    while( thickened_BODY_list.size() ) api_delent( thickened_BODY_list.pop() );
    PRINT_ERROR( "unable to setup target -- aborting\n" );
    return CUBIT_FAILURE;
  }

  // Needed for preview - target might be on BODY we are tweaking and we need
  // to check for that in the preview.
  DLIList<AcisBridge*> preview_ab_list;
  preview_ab_list.append( ATTRIB_CUBIT_OWNER::cubit_owner(target_EDGE_ptr) );

  // Call private function to do the work
  DLIList<EDGE*> EDGE_list( curve_list.size() );
  if( get_EDGEs( curve_list, EDGE_list ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  CubitStatus status;
  status = tweak_target_single( EDGE_list, target_FACE_ptr, new_bodysm_list,
    reverse_flg, keep_old_bodies, preview, &preview_ab_list );

  // Delete thickened BODY
  while( thickened_BODY_list.size() ) api_delent( thickened_BODY_list.pop() );

  return status;
}

//=============================================================================
// Function   : tweak_target_multiple
// Member Type: PRIVATE
// Description: Tweak specified edges of a sheet body or bodies up to a set
//              of target curves that are part of a sheet body or bodies.  The
//              target is a set of surfaces created by thickening the owning
//              surfaces of the target curves.
// Author     : Steve Storm
// Date       : 02/13/06
//=============================================================================
CubitStatus
AcisTweakTool::tweak_target_multiple( DLIList<Curve*> &input_curve_list,
                                      DLIList<Curve*> &target_curve_list,
                                      DLIList<BodySM*> &new_bodysm_list,
                                      CubitBoolean reverse_flg,
                                      CubitBoolean keep_old_bodies,
                                      CubitBoolean preview )
{
#ifndef ACIS_LOCAL_OPS
  PRINT_ERROR( "The ACIS Local Operations Husk is required for extending\n"
    "       surfaces.  It has not been licensed for this installation.\n" );
  return CUBIT_FAILURE;
#endif

  // Get source curves
  DLIList<EDGE*> input_EDGE_list( input_curve_list.size() );
  if( get_EDGEs( input_curve_list, input_EDGE_list ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  // Get target FACES
  DLIList<EDGE*> target_EDGE_list;
  if( get_EDGEs( target_curve_list, target_EDGE_list ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  // Build list of FACEs which could be from multiple BODIES
  DLIList<FACE*> target_FACE_list;
  DLIList<BODY*> thickened_BODY_list;
  while( target_EDGE_list.size() )
  {
    BODY *BODY_ptr;
    DLIList<EDGE*> removed_EDGE_list;
    DLIList<BODY*> tmp_thickened_BODY_list;
    DLIList<DLIList<EDGE*>*> output_EDGE_lists;
    DLIList<DLIList<FACE*>*> output_FACE_lists;
    DLIList<DLIList<FACE*>*> conjugate_FACE_lists;
    if( get_thickened_BODIES_of_EDGES( "tweak to", target_EDGE_list,
      removed_EDGE_list, BODY_ptr, tmp_thickened_BODY_list, output_EDGE_lists,
      output_FACE_lists, conjugate_FACE_lists ) == CUBIT_FAILURE )
    {
      PRINT_ERROR( "unable to setup target -- aborting\n" );
      return CUBIT_FAILURE;
    }

    // The conjugate FACEs are what we are after
    int i, j;
    conjugate_FACE_lists.reset();
    for( i=conjugate_FACE_lists.size(); i--; )
    {
      DLIList<FACE*> *FACE_list_ptr = conjugate_FACE_lists.get_and_step();

      FACE *FACE_ptr;
      FACE_list_ptr->reset();
      for( j=FACE_list_ptr->size(); j--; )
      {
        FACE_ptr = FACE_list_ptr->get_and_step();
        target_FACE_list.append( FACE_ptr );
      }
    }

    // Clean up memory allocated for the output lists
    while( output_EDGE_lists.size() ) delete output_EDGE_lists.pop();
    while( output_EDGE_lists.size() ) delete output_EDGE_lists.pop();
    while( conjugate_FACE_lists.size() ) delete conjugate_FACE_lists.pop();

    thickened_BODY_list += tmp_thickened_BODY_list;
  }

  // Needed for preview - targets might be on BODY we are tweaking and we need
  // to check for that in the preview.
  DLIList<AcisBridge*> preview_ab_list;
  if( preview ) get_owner_list( target_EDGE_list, preview_ab_list );

  // Call private function to do the work
  if( tweak_target_multiple( input_EDGE_list, target_FACE_list,
    new_bodysm_list, reverse_flg, keep_old_bodies, preview, &preview_ab_list )
    == CUBIT_FAILURE )
  {
    while( thickened_BODY_list.size() ) api_delent( thickened_BODY_list.pop() );
    return CUBIT_FAILURE;
  }

  // Delete the thickened target BODIES
  while( thickened_BODY_list.size() ) api_delent( thickened_BODY_list.pop() );

  return CUBIT_SUCCESS;
}

CubitStatus AcisTweakTool::get_ACIS_surfaces( DLIList<Surface*>& surf_list,
                                              DLIList<SurfaceACIS*>& acis_list )
{
  int i;
  surf_list.reset();
  for( i=surf_list.size(); i--; )
  {
    Surface* surf = surf_list.get_and_step();
    if (SurfaceACIS* surf_acis = dynamic_cast<SurfaceACIS*>(surf))
      acis_list.append(surf_acis);
  }

  if (surf_list.size() != acis_list.size())
  {
    PRINT_ERROR("%d non-ACIS Surfaces.\n", surf_list.size() - acis_list.size());
    return CUBIT_FAILURE;
  }

  return CUBIT_SUCCESS;
}

CubitStatus AcisTweakTool::get_FACEs( DLIList<Surface*> &surface_list,
                                      DLIList<FACE*> &FACE_list )
{
  int i;
  FACE *FACE_ptr;
  surface_list.reset();
  for( i=surface_list.size(); i--; )
  {
    Surface* surface_ptr = surface_list.get_and_step();
    FACE_ptr = AcisQueryEngine::instance()->get_FACE( surface_ptr );
    if( FACE_ptr ) FACE_list.append( FACE_ptr );
  }

  if( surface_list.size() != FACE_list.size() )
  {
    PRINT_ERROR("%d non-ACIS surfaces.\n", surface_list.size() - FACE_list.size());
    return CUBIT_FAILURE;
  }

  return CUBIT_SUCCESS;
}

CubitStatus AcisTweakTool::get_EDGEs( DLIList<Curve*> &curve_list,
                                      DLIList<EDGE*> &EDGE_list )
{
  int i;
  EDGE *EDGE_ptr;
  curve_list.reset();
  for( i=curve_list.size(); i--; )
  {
    Curve* curve_ptr = curve_list.get_and_step();
    EDGE_ptr = AcisQueryEngine::instance()->get_EDGE( curve_ptr );
    if( EDGE_ptr ) EDGE_list.append( EDGE_ptr );
  }

  if( curve_list.size() != EDGE_list.size() )
  {
    PRINT_ERROR("%d non-ACIS curves.\n", curve_list.size() - EDGE_list.size());
    return CUBIT_FAILURE;
  }

  return CUBIT_SUCCESS;
}

CubitStatus AcisTweakTool::get_VERTICEs( DLIList<Point*> &point_list,
                                         DLIList<VERTEX*> &VERTEX_list )
{
  int i;
  VERTEX *VERTEX_ptr;
  point_list.reset();
  for( i=point_list.size(); i--; )
  {
    Point* point_ptr = point_list.get_and_step();
    VERTEX_ptr = AcisQueryEngine::instance()->get_VERTEX( point_ptr );
    if( VERTEX_ptr ) VERTEX_list.append( VERTEX_ptr );
  }

  if( point_list.size() != VERTEX_list.size() )
  {
    PRINT_ERROR("%d non-ACIS vertices.\n", point_list.size() - VERTEX_list.size());
    return CUBIT_FAILURE;
  }

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : tweak_target
// Member Type: PRIVATE
// Description: Tweak specified EDGES of a surface or set of surfaces (in
//              sheet bodies) up to a single target FACE.
// Author     : Steve Storm
// Date       : 03/28/05
//=============================================================================
CubitStatus
AcisTweakTool::tweak_target_single( DLIList<EDGE*> &input_EDGE_list,
                                    FACE *target_FACE_ptr,
                                    DLIList<BodySM*> &new_bodysm_list,
                                    CubitBoolean reverse_flg,
                                    CubitBoolean keep_old_bodies,
                                    CubitBoolean preview,
                                    DLIList<AcisBridge*> *t_ab_list_ptr )
{
  assert( input_EDGE_list.size() );

  int i, j;

  bool delete_attribs =
    (GeometryModifyTool::instance()->get_new_ids() || keep_old_bodies);

  // Copy the input EDGE list, since we will be removing EDGEs from it
  DLIList<EDGE*> copied_EDGE_list = input_EDGE_list;

  copied_EDGE_list.reset();
  while( copied_EDGE_list.size() )
  {
    BODY *BODY_ptr;
    DLIList<EDGE*> removed_EDGE_list;
    DLIList<BODY*> thickened_BODY_list;
    DLIList<DLIList<EDGE*>*> output_EDGE_lists;
    DLIList<DLIList<FACE*>*> output_FACE_lists;
    DLIList<DLIList<FACE*>*> conjugate_FACE_lists;
    if( get_thickened_BODIES_of_EDGES( "tweak", copied_EDGE_list, 
      removed_EDGE_list, BODY_ptr, thickened_BODY_list, output_EDGE_lists,
      output_FACE_lists, conjugate_FACE_lists ) == CUBIT_FAILURE )
    {
      // Return success if any bodies were created
      return new_bodysm_list.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
    }

    BodySM *bodysm_ptr = AQE->get_body_sm_of_ENTITY( BODY_ptr );

    // Store source and target EDGEs for preview
    DLIList<AcisBridge*> preview_ab_list;
    if( preview )
    {
      get_owner_list( removed_EDGE_list, preview_ab_list );
      if( t_ab_list_ptr ) preview_ab_list += *t_ab_list_ptr;
    }
   
    // Tweak the surfaces on the output BODIEs, pull the tweaked surfaces off,
    // and replace them in the original BODY
    // First copy the BODY
    BODY *copied_BODY_ptr = AME->copy_BODY(BODY_ptr, CUBIT_FALSE);

    BODY *thickened_BODY_ptr;
    DLIList<FACE*> *output_FACE_list_ptr;
    DLIList<FACE*> *conjugate_FACE_list_ptr;
    thickened_BODY_list.reset();
    output_EDGE_lists.reset();
    output_FACE_lists.reset();
    conjugate_FACE_lists.reset();
    for( i=thickened_BODY_list.size(); i--; )
    {
      thickened_BODY_ptr = thickened_BODY_list.get_and_step();
      output_FACE_list_ptr = output_FACE_lists.get_and_step();
      conjugate_FACE_list_ptr = conjugate_FACE_lists.get_and_step();

      // Note - the conjugate FACE list can have duplicate entries in
      // it, as thicken doesn't always create a separate surface for
      // each curve.
      DLIList<FACE*> tweak_FACE_list;
      conjugate_FACE_list_ptr->reset();
      for( j=conjugate_FACE_list_ptr->size(); j--; )
        tweak_FACE_list.append_unique( conjugate_FACE_list_ptr->get_and_step() );

      // Keep track of the FACEs of interest, via their Cubit owners, which
      // will survive throughout the operation
      DLIList<AcisBridge*> owner_list;
      get_owner_list( *output_FACE_list_ptr, owner_list );

      // Do the tweak
      if( tweak_FACEs_to_target( tweak_FACE_list, target_FACE_ptr, reverse_flg )
        == CUBIT_FAILURE )
      {
        while( output_EDGE_lists.size() ) delete output_EDGE_lists.pop();
        while( output_FACE_lists.size() ) delete output_FACE_lists.pop();
        while( conjugate_FACE_lists.size() ) delete conjugate_FACE_lists.pop();
        while( thickened_BODY_list.size() )
          api_delent( thickened_BODY_list.pop() );
        api_delent( copied_BODY_ptr );
        // Return success if any bodies were created
        return new_bodysm_list.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
      }

      // Get bodies ready to swap new surfaces for old
      if( prep_for_surface_swap( thickened_BODY_ptr, copied_BODY_ptr,
        owner_list ) == CUBIT_FAILURE )
      {
        while( output_EDGE_lists.size() ) delete output_EDGE_lists.pop();
        while( output_FACE_lists.size() ) delete output_FACE_lists.pop();
        while( conjugate_FACE_lists.size() ) delete conjugate_FACE_lists.pop();
        while( thickened_BODY_list.size() )
          api_delent( thickened_BODY_list.pop() );
        api_delent( copied_BODY_ptr );
        // Return success if any bodies were created
        return new_bodysm_list.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
      }

    } // End loop on thickened (separated) BODIES

    // Free memory
    while( output_EDGE_lists.size() ) delete output_EDGE_lists.pop();
    while( output_FACE_lists.size() ) delete output_FACE_lists.pop();
    while( conjugate_FACE_lists.size() ) delete conjugate_FACE_lists.pop();

    // Unite the thickened BODIEs back into the copied_input_BODY_ptr
    BODY *master;
    if( unite_BODIES( copied_BODY_ptr, thickened_BODY_list, master ) == CUBIT_FAILURE )
    {
      // If failure, the unite_BODIES function cleaned up the memory
      return new_bodysm_list.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
    }

    // This BODY is done (whew!)

    if( !preview )
    {
      // Now cleanout the owner attributes from the copied BODY, if required
      if( delete_attribs )
        AQE->remove_cubit_owner_attrib_in_BODY( master );

      BodySM *new_body = AME->get_new_Body( bodysm_ptr, BODY_ptr, master,
        keep_old_bodies );

      if( new_body )
        new_bodysm_list.append( new_body );
    }
    else
    {
      GfxPreview::clear();

      draw_tweak_preview_omt( master, CUBIT_TRUE, &preview_ab_list );
      api_delent( master );
    }
  }

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : tweak_target_multiple
// Member Type: PRIVATE
// Description: Tweak specified EDGES of a surface or set of surfaces (in
//              sheet bodies) up to a multiple target FACEs.
// Author     : Steve Storm
// Date       : 03/28/05
//=============================================================================
CubitStatus
AcisTweakTool::tweak_target_multiple( DLIList<EDGE*> &input_EDGE_list,
                                      DLIList<FACE*> &target_FACE_list,
                                      DLIList<BodySM*> &new_bodysm_list,
                                      CubitBoolean reverse_flg,
                                      CubitBoolean keep_old_bodies,
                                      CubitBoolean preview,
                                      DLIList<AcisBridge*> *preview_ab_list_ptr )
{
  assert( input_EDGE_list.size() );

  int i, j;

  bool delete_attribs =
    (GeometryModifyTool::instance()->get_new_ids() || keep_old_bodies);

  // Get extended target surfaces (in a single sheet BODY)
  BODY *ext_target_BODY_ptr;
  if( create_extended_sheet( target_FACE_list, ext_target_BODY_ptr )
    == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  if( DEBUG_FLAG(168) )
  {
    BODY *tmp_BODY = AME->copy_BODY( ext_target_BODY_ptr, CUBIT_TRUE );
    BodySM *this_bodysm = AQE->populate_topology_bridges( tmp_BODY );
    new_bodysm_list.append( this_bodysm );
  }

  // Copy the input EDGE list, since we will be removing EDGEs from it
  DLIList<EDGE*> copied_EDGE_list = input_EDGE_list;

  copied_EDGE_list.reset();
  while( copied_EDGE_list.size() )
  {
    BODY *BODY_ptr;
    DLIList<EDGE*> removed_EDGE_list;
    DLIList<BODY*> thickened_BODY_list;
    DLIList<DLIList<EDGE*>*> output_EDGE_lists;
    DLIList<DLIList<FACE*>*> output_FACE_lists;
    DLIList<DLIList<FACE*>*> conjugate_FACE_lists;
    if( get_thickened_BODIES_of_EDGES( "tweak", copied_EDGE_list, 
      removed_EDGE_list, BODY_ptr, thickened_BODY_list, output_EDGE_lists,
      output_FACE_lists, conjugate_FACE_lists ) == CUBIT_FAILURE )
    {
      api_delent( ext_target_BODY_ptr );
      // Return success if any bodies were created
      return new_bodysm_list.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
    }

    BodySM *bodysm_ptr = AQE->get_body_sm_of_ENTITY( BODY_ptr );

    // Tweak the surfaces on the output BODIEs, pull the tweaked surfaces off,
    // and replace them in the original BODY
    // First copy the BODY
    BODY *copied_BODY_ptr = AME->copy_BODY(BODY_ptr, CUBIT_FALSE);

    // Store source and target EDGEs for preview
    DLIList<AcisBridge*> preview_ab_list;
    if( preview )
    {
      get_owner_list( removed_EDGE_list, preview_ab_list );
      if( preview_ab_list_ptr ) preview_ab_list += *preview_ab_list_ptr;
    }

    BODY *thickened_BODY_ptr;
    DLIList<FACE*> *output_FACE_list_ptr;
    DLIList<FACE*> *conjugate_FACE_list_ptr;
    thickened_BODY_list.reset();
    output_EDGE_lists.reset();
    output_FACE_lists.reset();
    conjugate_FACE_lists.reset();
    for( i=thickened_BODY_list.size(); i--; )
    {
      thickened_BODY_ptr = thickened_BODY_list.get_and_step();
      output_FACE_list_ptr = output_FACE_lists.get_and_step();
      conjugate_FACE_list_ptr = conjugate_FACE_lists.get_and_step();

      // Note - the conjugate FACE list can have duplicate entries in
      // it, as thicken doesn't always create a separate surface for
      // each curve.
      DLIList<FACE*> tweak_FACE_list;
      conjugate_FACE_list_ptr->reset();
      for( j=conjugate_FACE_list_ptr->size(); j--; )
        tweak_FACE_list.append_unique( conjugate_FACE_list_ptr->get_and_step() );

      // Keep track of the FACEs of interest, via their Cubit owners, which
      // will survive throughout the operation
      DLIList<AcisBridge*> owner_list;
      get_owner_list( *output_FACE_list_ptr, owner_list );

      // Do the tweak
      if( tweak_target_multiple( tweak_FACE_list, target_FACE_list,
        ext_target_BODY_ptr, new_bodysm_list, reverse_flg ) 
        == CUBIT_FAILURE )
      {
        while( output_EDGE_lists.size() ) delete output_EDGE_lists.pop();
        while( output_FACE_lists.size() ) delete output_FACE_lists.pop();
        while( conjugate_FACE_lists.size() ) delete conjugate_FACE_lists.pop();
        while( thickened_BODY_list.size() )
          api_delent( thickened_BODY_list.pop() );
        api_delent( ext_target_BODY_ptr );
        api_delent( copied_BODY_ptr );
        // Return success if any bodies were created
        return new_bodysm_list.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
      }

      // Get bodies ready to swap new surfaces for old
      if( prep_for_surface_swap( thickened_BODY_ptr, copied_BODY_ptr,
        owner_list ) == CUBIT_FAILURE )
      {
        while( output_EDGE_lists.size() ) delete output_EDGE_lists.pop();
        while( output_FACE_lists.size() ) delete output_FACE_lists.pop();
        while( conjugate_FACE_lists.size() ) delete conjugate_FACE_lists.pop();
        while( thickened_BODY_list.size() )
          api_delent( thickened_BODY_list.pop() );
        api_delent( ext_target_BODY_ptr );
        api_delent( copied_BODY_ptr );
        // Return success if any bodies were created
        return new_bodysm_list.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
      }

    } // End loop on thickened (separated) BODIES

    // Free memory
    while( output_EDGE_lists.size() ) delete output_EDGE_lists.pop();
    while( output_FACE_lists.size() ) delete output_FACE_lists.pop();
    while( conjugate_FACE_lists.size() ) delete conjugate_FACE_lists.pop();

    // Unite the thickened BODIEs back into the copied_input_BODY_ptr
    BODY *master;
    if( unite_BODIES( copied_BODY_ptr, thickened_BODY_list, master ) == CUBIT_FAILURE )
    {
      api_delent( ext_target_BODY_ptr );
      // If failure, the unite_BODIES function cleaned up the memory
      return new_bodysm_list.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
    }

    // This BODY is done (whew!)

    if( !preview )
    {
      // Now cleanout the owner attributes from the copied BODY, if required
      if( delete_attribs )
        AQE->remove_cubit_owner_attrib_in_BODY( master );

      BodySM *new_body = AME->get_new_Body( bodysm_ptr, BODY_ptr, master,
        keep_old_bodies );

      if( new_body )
        new_bodysm_list.append( new_body );
    }
    else
    {
      GfxPreview::clear();

      // Preview EDGEs
      draw_tweak_preview_omt( master, CUBIT_TRUE, &preview_ab_list );
      api_delent( master );
    }
  }

  api_delent( ext_target_BODY_ptr );

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : copy_FACES_from_BODY
// Member Type: PRIVATE
// Description: Input FACE list must be in a single BODY.  Copy these FACEs off
//              into a new BODY, while retaining all the Cubit attributes.  The
//              input BODY (the parent BODY of the FACE_list) is NOT modified.
// Author     : Steve Storm
// Date       : 03/28/05
//=============================================================================
CubitStatus
AcisTweakTool::copy_FACES_from_BODY( DLIList<FACE*> &input_FACE_list,
                                     BODY *&copied_BODY_ptr)
{
  // Method - copy the body and remove the other faces from the body.  We tried
  // copying the faces off into a new body using api_unhook_faces, or
  // copy_single_entity, etc., but the cubit attributes on the copied faces
  // always seem to disappear.  This method, which may seem like an odd way
  // to approach the problem, is guaranteed to keep the cubit attributes on
  // the resultant body.  It also guarantees a single body as a result (note
  // the resultant body may need to be separated, as it may have disconnected
  // faces).
  int i;
  FACE *FACE_ptr;
  BODY *input_BODY_ptr;

  // Get owning BODY of the FACE_list
  input_FACE_list.reset();
  FACE_ptr = input_FACE_list.get_and_step();
  input_BODY_ptr = AcisQueryEngine::instance()->get_BODY_of_ENTITY( FACE_ptr );

  // Error check - make sure remaining FACEs are from same BODY
  for( i=1; i<input_FACE_list.size(); i++ )
  {
    FACE_ptr = input_FACE_list.get_and_step();
    if( input_BODY_ptr != AcisQueryEngine::instance()->get_BODY_of_ENTITY( FACE_ptr ) )
      return CUBIT_FAILURE;
  }

  // Get a list of all the associated Surfaces to the input FACEs
  AcisBridge *ab_ptr;
  SurfaceACIS *asurf_ptr;
  DLIList<SurfaceACIS*> input_asurface_list( input_FACE_list.size() );
  input_FACE_list.reset();
  for( i=input_FACE_list.size(); i--; )
  {
    FACE_ptr = input_FACE_list.get_and_step();
    ab_ptr = ATTRIB_CUBIT_OWNER::cubit_owner(FACE_ptr);
    asurf_ptr = CAST_TO( ab_ptr, SurfaceACIS );
    if( asurf_ptr == NULL )
    {
      PRINT_ERROR( "Internal error -- please report.\n" );
      return CUBIT_FAILURE;
    }
    input_asurface_list.append( asurf_ptr );
  }

  // Copy the input BODY, keeping all the cubit attributes
  copied_BODY_ptr = AcisModifyEngine::instance()->copy_BODY(input_BODY_ptr,
    CUBIT_FALSE);

  DLIList<FACE*> remove_FACE_list;

  // Remove the "unused" FACEs from the copied BODY.  Loop on FACEs in
  // copied_BODY_ptr

  DLIList<FACE*> FACE_list;
  AcisQueryEngine::instance()->get_FACEs( copied_BODY_ptr, FACE_list );
  for( i=FACE_list.size(); i--; )
  {
    FACE_ptr = FACE_list.get_and_step();

    ab_ptr = ATTRIB_CUBIT_OWNER::cubit_owner(FACE_ptr);
    asurf_ptr = CAST_TO( ab_ptr, SurfaceACIS );

    // If we don't find a matching asurf_ptr in input_asurface_list remove
    // this FACE
    if( !asurf_ptr || (asurf_ptr && !input_asurface_list.is_in_list( asurf_ptr )) )
      remove_FACE_list.append( FACE_ptr );
  }

  if( remove_FACES_from_BODY( copied_BODY_ptr, remove_FACE_list ) == CUBIT_FAILURE )
  {
    api_delent( copied_BODY_ptr );
    copied_BODY_ptr = NULL;
  }

  return CUBIT_SUCCESS;
}

CubitStatus
AcisTweakTool::remove_FACES_from_BODY( BODY *BODY_ptr,
                                       DLIList<FACE*> &remove_FACE_list )
{
  if( remove_FACE_list.size() == 0 )
    return CUBIT_SUCCESS;

  // I'm not sure why, but this operation can cause cubit owner attributes to
  // disappear from copied_BODY_ptr.  We need to put them back.
  AcisBridge *ab_body_ptr;
  AcisBridge *ab_lump_ptr;
  AcisBridge *ab_shell_ptr;
  get_owner_attribs( BODY_ptr, ab_body_ptr, ab_lump_ptr, ab_shell_ptr );

  int i;
  outcome result;
  //ENTITY_LIST faces_to_remove;
  FACE *FACE_ptr;
  for( i=remove_FACE_list.size(); i--; )
  {
    FACE_ptr = remove_FACE_list.get_and_step();
    //faces_to_remove.add( FACE_ptr );
    result = api_remove_face( FACE_ptr );
    if( !result.ok() )
    {
      AcisQueryEngine::instance()->ACIS_API_error(result);
      return CUBIT_FAILURE;
    }
  }

  // This doesn't work, because all attributes are removed from the input body.
  // I think this would be faster than the method we are using.
  // This must be an acis bug...
  //ENTITY_LIST unhooked_bodies;
  //result = api_unhook_faces( faces_to_remove, false, unhooked_bodies );
  //unhooked_bodies.init();
  //ENTITY *ENTITY_ptr = NULL;
  //while( ENTITY_ptr = unhooked_bodies.next() )
  //  api_delent( ENTITY_ptr );
  //if( !result.ok() )
  //{
  //  AcisQueryEngine::instance()->ACIS_API_error(result);
  //  return CUBIT_FAILURE;
  //}

  // Put cubit owner attributes back
  reset_owner_attribs( BODY_ptr, ab_body_ptr, ab_lump_ptr, ab_shell_ptr );

  // For some reason, removing faces causes the remaining faces to be single sided
  result = api_body_to_2d( BODY_ptr );
  if( !result.ok() )
  {
    AcisQueryEngine::instance()->ACIS_API_error(result);
    return CUBIT_FAILURE;
  }

  return CUBIT_SUCCESS;
}

CubitStatus
AcisTweakTool::remove_FACES_from_BODY_except( BODY *BODY_ptr,
                                              DLIList<FACE*> &keep_FACE_list )
{
  DLIList<FACE*> FACE_list;
  AcisQueryEngine::instance()->get_FACEs( BODY_ptr, FACE_list );

  DLIList<FACE*> remove_FACE_list;
  int i;
  FACE *FACE_ptr;
  for( i=FACE_list.size(); i--; )
  {
    FACE_ptr = FACE_list.get_and_step();
    if( !keep_FACE_list.is_in_list( FACE_ptr ) )
      remove_FACE_list.append( FACE_ptr );
  }

  return remove_FACES_from_BODY( BODY_ptr, remove_FACE_list );
}

CubitStatus
AcisTweakTool::thicken_BODY( BODY *BODY_ptr, double thickness )
{
#if CUBIT_ACIS_VERSION >= 1600
  // ACIS v16 sometimes removes the Cubit owner attributes when thickening -
  // this workaround puts them back.

  // Setup tweak attributes so we can preserve Cubit owners on vertices,
  // edges and faces
  DLIList<FACE*> pre_FACE_list;
  DLIList<EDGE*> pre_EDGE_list;
  DLIList<VERTEX*> pre_VERTEX_list;
  DLIList<AcisBridge*> ab_FACE_list, ab_EDGE_list, ab_VERTEX_list;
  assign_tweak_attribs( BODY_ptr, "ttweak", pre_FACE_list, ab_FACE_list,
    pre_EDGE_list, ab_EDGE_list, pre_VERTEX_list, ab_VERTEX_list );

  // Also preserve on body, lump and shell
  AcisBridge *ab_body_ptr;
  AcisBridge *ab_lump_ptr;
  AcisBridge *ab_shell_ptr;
  get_owner_attribs( BODY_ptr, ab_body_ptr, ab_lump_ptr, ab_shell_ptr );
#endif

  // Keep track of input VERTICEs, EDGEs and FACEs for later use
  DLIList<VERTEX*> input_VERTEX_list;
  AQE->get_VERTICEs( BODY_ptr, input_VERTEX_list );

  DLIList<EDGE*> input_EDGE_list;
  AQE->get_EDGEs( BODY_ptr, input_EDGE_list );

  DLIList<FACE*> input_FACE_list;
  AQE->get_FACEs( BODY_ptr, input_FACE_list );

  // Thicken the body
  SPAposition box_l(0,0,0);
  SPAposition box_h(0,0,0);
  outcome result;
  result = api_sheet_thicken( BODY_ptr, thickness, false, box_l, box_h);
  if (!result.ok())
  {
    AQE->ACIS_API_error(result);
    // Thicken doesn't give any sort of error.  Sometimes it won't thicken
    // because of bad geometry, so just give a hint.
    PRINT_ERROR( "geometry problem encountered - try healing the body first\n" );
    return CUBIT_FAILURE;
  }

  // Thicken copies the cubit owner attributes from the original entities
  // to the other side.  Remove those.
  int i;

  VERTEX *VERTEX_ptr;
  DLIList<VERTEX*> VERTEX_list;
  AQE->get_VERTICEs( BODY_ptr, VERTEX_list );
  for( i=VERTEX_list.size(); i--; )
  {
    VERTEX_ptr = VERTEX_list.get_and_step();
    if( !input_VERTEX_list.is_in_list( VERTEX_ptr ) )
      ATTRIB_CUBIT_OWNER::remove_cubit_owner(VERTEX_ptr);
  }

  EDGE *EDGE_ptr;
  DLIList<EDGE*> EDGE_list;
  AQE->get_EDGEs( BODY_ptr, EDGE_list );
  for( i=EDGE_list.size(); i--; )
  {
    EDGE_ptr = EDGE_list.get_and_step();
    if( !input_EDGE_list.is_in_list( EDGE_ptr ) )
      ATTRIB_CUBIT_OWNER::remove_cubit_owner(EDGE_ptr);
  }

  FACE *FACE_ptr;
  DLIList<FACE*> FACE_list;
  AQE->get_FACEs( BODY_ptr, FACE_list );
  for( i=FACE_list.size(); i--; )
  {
    FACE_ptr = FACE_list.get_and_step();
    if( !input_FACE_list.is_in_list( FACE_ptr ) )
      ATTRIB_CUBIT_OWNER::remove_cubit_owner(FACE_ptr);
  }

#if CUBIT_ACIS_VERSION >= 1600
  // Put cubit owner attributes back on body, lump and shell
  reset_owner_attribs( BODY_ptr, ab_body_ptr, ab_lump_ptr, ab_shell_ptr );

  // Replace Cubit owners on vertices, edges, and faces
  reassign_cubit_owners_from_tweak_attribs( BODY_ptr, "ttweak",
    pre_FACE_list, ab_FACE_list, pre_EDGE_list, ab_EDGE_list, pre_VERTEX_list,
    ab_VERTEX_list );

  // Remove thicken tweak attributes
  remove_named_attribs( BODY_ptr, "ttweak" );
#endif

  return CUBIT_SUCCESS;
}

VERTEX *
AcisTweakTool::find_corresponding_VERTEX( VERTEX *ref_VERTEX_ptr,
                                          DLIList<VERTEX*> &VERTEX_list )
{
  // Get the reference CUBIT owner attribute
  AcisBridge *ab_ptr = ATTRIB_CUBIT_OWNER::cubit_owner(ref_VERTEX_ptr);
  if( !ab_ptr )
    return NULL;

  PointACIS *ref_point_acis_ptr = CAST_TO( ab_ptr, PointACIS );
  if( !ref_point_acis_ptr )
    return NULL;

  // Find the point with the same CUBIT owner attribute in the input VERTEX list
  VERTEX *VERTEX_ptr;
  PointACIS *point_acis_ptr;
  int i;
  for( i=VERTEX_list.size(); i--; )
  {
    VERTEX_ptr = VERTEX_list.get_and_step();

    ab_ptr = ATTRIB_CUBIT_OWNER::cubit_owner(VERTEX_ptr);
    if( !ab_ptr )
      continue;

    point_acis_ptr = CAST_TO( ab_ptr, PointACIS );
    if( !point_acis_ptr )
      continue;

    if( point_acis_ptr == ref_point_acis_ptr )
      return VERTEX_ptr;
  }

  // No match found
  return NULL;
}

EDGE *
AcisTweakTool::find_corresponding_EDGE( EDGE *ref_EDGE_ptr,
                                        DLIList<EDGE*> &EDGE_list )
{
  // Get the reference CUBIT owner attribute
  AcisBridge *ab_ptr = ATTRIB_CUBIT_OWNER::cubit_owner(ref_EDGE_ptr);
  if( !ab_ptr )
    return NULL;

  CurveACIS *ref_curve_acis_ptr = CAST_TO( ab_ptr, CurveACIS );
  if( !ref_curve_acis_ptr )
    return NULL;

  // Find the curve with the same CUBIT owner attribute in the input EDGE list
  EDGE *EDGE_ptr;
  CurveACIS *curve_acis_ptr;
  int i;
  for( i=EDGE_list.size(); i--; )
  {
    EDGE_ptr = EDGE_list.get_and_step();

    ab_ptr = ATTRIB_CUBIT_OWNER::cubit_owner(EDGE_ptr);
    if( !ab_ptr )
      continue;

    curve_acis_ptr = CAST_TO( ab_ptr, CurveACIS );
    if( !curve_acis_ptr )
      continue;

    if( curve_acis_ptr == ref_curve_acis_ptr )
      return EDGE_ptr;
  }

  // No match found
  return NULL;
}


LOOP *
AcisTweakTool::find_corresponding_LOOP( LOOP *ref_LOOP_ptr,
                                        DLIList<LOOP*> &LOOP_list )
{
  // Get the reference CUBIT owner attribute
  AcisBridge *ab_ptr = ATTRIB_CUBIT_OWNER::cubit_owner(ref_LOOP_ptr);
  if( !ab_ptr )
    return NULL;

  LoopACIS *ref_loop_acis_ptr = CAST_TO( ab_ptr, LoopACIS );
  if( !ref_loop_acis_ptr )
    return NULL;

  // Find the loop with the same CUBIT owner attribute in the input LOOP list
  LOOP *LOOP_ptr;
  LoopACIS *loop_acis_ptr;
  int i;
  for( i=LOOP_list.size(); i--; )
  {
    LOOP_ptr = LOOP_list.get_and_step();

    ab_ptr = ATTRIB_CUBIT_OWNER::cubit_owner(LOOP_ptr);
    if( !ab_ptr )
      continue;

    loop_acis_ptr = CAST_TO( ab_ptr, LoopACIS );
    if( !loop_acis_ptr )
      continue;

    if( loop_acis_ptr == ref_loop_acis_ptr )
      return LOOP_ptr;
  }

  // No match found
  return NULL;
}

FACE *
AcisTweakTool::find_corresponding_FACE( FACE *ref_FACE_ptr,
                                        DLIList<FACE*> &FACE_list )
{
  // Get the reference CUBIT owner attribute
  AcisBridge *ab_ptr = ATTRIB_CUBIT_OWNER::cubit_owner(ref_FACE_ptr);
  if( !ab_ptr )
    return NULL;

  SurfaceACIS *ref_surf_acis_ptr = CAST_TO( ab_ptr, SurfaceACIS );
  if( !ref_surf_acis_ptr )
    return NULL;

  // Find the surface with the same CUBIT owner attribute in the input FACE list
  FACE *FACE_ptr;
  SurfaceACIS *surf_acis_ptr;
  int i;
  for( i=FACE_list.size(); i--; )
  {
    FACE_ptr = FACE_list.get_and_step();

    ab_ptr = ATTRIB_CUBIT_OWNER::cubit_owner(FACE_ptr);
    if( !ab_ptr )
      continue;

    surf_acis_ptr = CAST_TO( ab_ptr, SurfaceACIS );
    if( !surf_acis_ptr )
      continue;

    if( surf_acis_ptr == ref_surf_acis_ptr )
      return FACE_ptr;
  }

  // No match found
  return NULL;
}

FACE *
AcisTweakTool::find_corresponding_FACE( AcisBridge *ab_ptr, BODY *BODY_ptr)
{
  int i;
  FACE *FACE_ptr;

  DLIList<FACE*> FACE_list;
  AcisQueryEngine::instance()->get_FACEs( BODY_ptr, FACE_list );

  FACE_list.reset();
  for( i=FACE_list.size(); i--; )
  {
    FACE_ptr = FACE_list.get_and_step();
    if( ATTRIB_CUBIT_OWNER::cubit_owner( FACE_ptr ) == ab_ptr )
      return FACE_ptr;
  }

  return NULL;
}

LOOP *
AcisTweakTool::find_corresponding_LOOP( AcisBridge *ab_ptr, BODY *BODY_ptr)
{
  int i;
  LOOP *LOOP_ptr;

  DLIList<LOOP*> LOOP_list;
  AcisQueryEngine::instance()->get_LOOPs( BODY_ptr, LOOP_list );

  LOOP_list.reset();
  for( i=LOOP_list.size(); i--; )
  {
    LOOP_ptr = LOOP_list.get_and_step();
    if( ATTRIB_CUBIT_OWNER::cubit_owner( LOOP_ptr ) == ab_ptr )
      return LOOP_ptr;
  }

  return NULL;
}

CubitStatus
AcisTweakTool::get_owner_list( DLIList<EDGE*> &EDGE_list,
                               DLIList<AcisBridge*> &owner_list )
{
  int i;
  EDGE *EDGE_ptr;

  EDGE_list.reset();
  for( i=EDGE_list.size(); i--; )
  {
    EDGE_ptr = EDGE_list.get_and_step();
    owner_list.append_unique( ATTRIB_CUBIT_OWNER::cubit_owner( EDGE_ptr ) );
  }

  return CUBIT_SUCCESS;
}

CubitStatus
AcisTweakTool::get_owner_list( DLIList<FACE*> &FACE_list,
                               DLIList<AcisBridge*> &owner_list )
{
  int i;
  FACE *FACE_ptr;

  FACE_list.reset();
  for( i=FACE_list.size(); i--; )
  {
    FACE_ptr = FACE_list.get_and_step();
    owner_list.append_unique( ATTRIB_CUBIT_OWNER::cubit_owner( FACE_ptr ) );
  }

  return CUBIT_SUCCESS;
}

CubitStatus
AcisTweakTool::get_owner_list( DLIList<LOOP*> &LOOP_list,
                               DLIList<AcisBridge*> &owner_list )
{
  int i;
  LOOP *LOOP_ptr;

  LOOP_list.reset();
  for( i=LOOP_list.size(); i--; )
  {
    LOOP_ptr = LOOP_list.get_and_step();
    owner_list.append_unique( ATTRIB_CUBIT_OWNER::cubit_owner( LOOP_ptr ) );
  }

  return CUBIT_SUCCESS;
}

CubitStatus
AcisTweakTool::get_corresponding_FACE_list( DLIList<AcisBridge*> &owner_list,
                                            BODY *BODY_ptr,
                                            DLIList<FACE*> &corresponding_FACE_list )
{
  int i;
  AcisBridge *ab_ptr;

  owner_list.reset();
  for( i=owner_list.size(); i--; )
  {
    ab_ptr = owner_list.get_and_step();
    corresponding_FACE_list.append( find_corresponding_FACE( ab_ptr, BODY_ptr ) );
  }

  return CUBIT_SUCCESS;
}

CubitStatus
AcisTweakTool::get_corresponding_LOOP_list( DLIList<AcisBridge*> &owner_list,
                                            BODY *BODY_ptr,
                                            DLIList<LOOP*> &corresponding_LOOP_list )
{
  int i;
  AcisBridge *ab_ptr;

  owner_list.reset();
  for( i=owner_list.size(); i--; )
  {
    ab_ptr = owner_list.get_and_step();
    corresponding_LOOP_list.append( find_corresponding_LOOP( ab_ptr, BODY_ptr ) );
  }

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : get_thickened_BODIES_of_EDGES
// Member Type: PRIVATE
// Description: Get thickened BODIES from an input list of EDGEs.  This
//              function can be called multiple times on the same list of EDGEs
//              - each time, the EDGEs from the common body that is thickened
//              are removed from the input EDGE list.  The thickened BODIES are
//              copied from the sheet the EDGEs are attached to (multiple
//              BODIES can be returned because the returned thickened BODIES
//              must be nonmanifold).  For each BODY returned, also output a
//              list of EDGES (corresponding to the original input EDGEs), a
//              list of FACEs(corresponding to the original FACEs of the input
//              BODY) and a list of conjugate FACEs (the "side" FACEs of the
//              thickened BODY - these are the FACEs we can tweak).
// Author     : Steve Storm
// Date       : 03/01/05
//=============================================================================
CubitStatus
AcisTweakTool::get_thickened_BODIES_of_EDGES( const char *command_name,
                                              DLIList<EDGE*> &input_EDGE_list,
                                              DLIList<EDGE*> &removed_EDGE_list,
                                              BODY *&common_BODY_ptr,
                                              DLIList<BODY*> &thickened_BODY_list,
                                              DLIList<DLIList<EDGE*>*> &output_EDGE_lists,
                                              DLIList<DLIList<FACE*>*> &output_FACE_lists,
                                              DLIList<DLIList<FACE*>*> &conjugate_FACE_lists,
                                              double thickness )
{
  assert( thickened_BODY_list.size() == 0 );

  int i, j;

  EDGE *EDGE_ptr;
  FACE *FACE_ptr;
  BODY *BODY_ptr;
  outcome result;

  // Pull EDGEs out of the input EDGE list that are from a common BODY
  input_EDGE_list.reset();
  common_BODY_ptr = NULL;
  for( i=input_EDGE_list.size(); i--; )
  {
    EDGE_ptr = input_EDGE_list.get();

    BODY_ptr = AQE->get_BODY_of_ENTITY( EDGE_ptr );

    if( common_BODY_ptr == NULL )
      common_BODY_ptr = BODY_ptr;

    if( common_BODY_ptr == BODY_ptr )
    {
      removed_EDGE_list.append( EDGE_ptr );
      input_EDGE_list.change_to( NULL );
    }

    input_EDGE_list.step();
  }
  input_EDGE_list.remove_all_with_value( NULL );

  // Get all FACEs attached to the EDGEs
  DLIList<FACE*> attached_FACE_list;
  removed_EDGE_list.reset();
  for( i=removed_EDGE_list.size(); i--; )
  {
    EDGE_ptr = removed_EDGE_list.get_and_step();
    DLIList<FACE*> tmp_FACE_list;
    AQE->get_FACEs( EDGE_ptr, tmp_FACE_list );

    if( tmp_FACE_list.size() == 0 )
    {
      PRINT_ERROR( "Cannot %s curves that are free\n",
        command_name );
      return CUBIT_FAILURE;
    }

    if( tmp_FACE_list.size() != 1 )
    {
      PRINT_ERROR( "Can only %s curves attached to one surface\n",
        command_name );
      return CUBIT_FAILURE;
    }

    attached_FACE_list.append( tmp_FACE_list.get() );
  }

  // Make sure that all of the FACEs are on a sheet body
  for( i=attached_FACE_list.size(); i--; )
  {
    FACE_ptr = attached_FACE_list.get_and_step();
    if( FACE_ptr->sides() != DOUBLE_SIDED )
    {
      PRINT_ERROR("Cannot %s curves that are not on sheet bodies\n",
        command_name );
      return CUBIT_FAILURE;
    }
  }

  // Copy these FACEs off into a new body
  // Copy the BODY
  BODY *copied_BODY_ptr;
  if( copy_FACES_from_BODY( attached_FACE_list, copied_BODY_ptr ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  // Make the BODY nonmanifold, as thicken can't handle manifold bodies.  Note
  // the api_unstitch_nonmani api will remove the cubit owner from the resultant
  // BODY, LUMP and SHELL - we manually put them back on, for proper updating
  // of the Cubit entities (luckily it leaves them on the FACEs, CURVEs, etc..
  // Also note the api destroys the input BODY.
  AcisBridge *ab_body_ptr;
  AcisBridge *ab_lump_ptr;
  AcisBridge *ab_shell_ptr;
  get_owner_attribs( copied_BODY_ptr, ab_body_ptr, ab_lump_ptr, ab_shell_ptr );

  BODY* lumps  = NULL;
  BODY* sheets = NULL;
  BODY* lamina = NULL;
  BODY* wires  = NULL;
  result = api_unstitch_nonmani( copied_BODY_ptr, lumps, sheets, lamina, wires );
  if (!result.ok())
  {
    AQE->ACIS_API_error(result);
    AQE->delete_ACIS_BODY( copied_BODY_ptr,CUBIT_TRUE );
    return CUBIT_FAILURE;
  }

  // Make sure cubit owner attributes remain on resultant BODY, LUMP, SHELL
  reset_owner_attribs( sheets, ab_body_ptr, ab_lump_ptr, ab_shell_ptr );

  // Separate bodies (thicken won't work if more than one vol per body)
  BODY **separated_BODY_array = NULL;
  int n_body = 0;

  // Note first BODY will always be the input BODY
  result = api_separate_body( sheets, n_body, separated_BODY_array );

  if (!result.ok())
  {
    AQE->ACIS_API_error(result);
    AQE->delete_ACIS_BODY( sheets, CUBIT_TRUE );
    return CUBIT_FAILURE;
  }

  // These BODIES can now be thickened
  for( i=0; i<n_body; i++ )
    thickened_BODY_list.append(separated_BODY_array[i]);

  // Loop on separated BODIES
  thickened_BODY_list.reset();
  for( i=thickened_BODY_list.size(); i--; )
  {
    BODY_ptr = thickened_BODY_list.get_and_step();

    // Make sure cubit owner attributes remain on resultant BODY, LUMP, SHELL
    reset_owner_attribs( BODY_ptr, ab_body_ptr, ab_lump_ptr, ab_shell_ptr );

    // The FACEs in this BODY can efficiently be added to the output list here
    DLIList<FACE*> *FACE_list_ptr = new DLIList<FACE*>;
    output_FACE_lists.append( FACE_list_ptr );
    AQE->get_FACEs( BODY_ptr, *FACE_list_ptr );

    // Thicken the BODY
    if( thicken_BODY( BODY_ptr, thickness ) == CUBIT_FAILURE )
    {
      // Delete all BODIEs, since we failed
      while( thickened_BODY_list.size() )
        AQE->delete_ACIS_BODY(thickened_BODY_list.pop(),CUBIT_TRUE);
      // Clean up memory allocated for the output lists
      while( output_EDGE_lists.size() ) delete output_EDGE_lists.pop();
      while( output_FACE_lists.size() ) delete output_FACE_lists.pop();
      while( conjugate_FACE_lists.size() ) delete conjugate_FACE_lists.pop();
      return CUBIT_FAILURE;
    }

    // Make sure cubit owner attributes remain on resultant BODY, LUMP, SHELL
    reset_owner_attribs( BODY_ptr, ab_body_ptr, ab_lump_ptr, ab_shell_ptr );

    // Now find EDGEs and conjugate FACEs for the thickened BODY
    DLIList<EDGE*> *EDGE_list_ptr = new DLIList<EDGE*>;
    output_EDGE_lists.append( EDGE_list_ptr );

    DLIList<EDGE*> tmp_EDGE_list;
    AQE->get_EDGEs( BODY_ptr, tmp_EDGE_list );
    EDGE *tmp_EDGE_ptr;
    for( j=removed_EDGE_list.size(); j--; )
    {
      EDGE_ptr = removed_EDGE_list.get_and_step();

      tmp_EDGE_ptr = find_corresponding_EDGE( EDGE_ptr, tmp_EDGE_list );
      if( tmp_EDGE_ptr )
        EDGE_list_ptr->append( tmp_EDGE_ptr );
    }

    if( !EDGE_list_ptr->size() )
    {
      PRINT_ERROR( "Internal error - please report.\n" );
      // Delete all BODIEs, since we failed
      while( thickened_BODY_list.size() )
        AQE->delete_ACIS_BODY(thickened_BODY_list.pop(),CUBIT_TRUE);
      // Clean up memory allocated for the output lists
      while( output_EDGE_lists.size() ) delete output_EDGE_lists.pop();
      while( output_FACE_lists.size() ) delete output_FACE_lists.pop();
      while( conjugate_FACE_lists.size() ) delete conjugate_FACE_lists.pop();
      return CUBIT_FAILURE;
    }

    // Find the "conjugate" FACEs (the FACEs that will be tweaked or removed
    // by the calling function).  These are the FACEs that adjoin the original
    // FACEs - they will be attached to the input EDGEs and won't exist in the
    // FACE_list_ptr.

    DLIList<FACE*> *conjugate_FACE_list_ptr = new DLIList<FACE*>;
    conjugate_FACE_lists.append( conjugate_FACE_list_ptr );

    EDGE_list_ptr->reset();
    for( j=EDGE_list_ptr->size(); j--; )
    {
      EDGE_ptr = EDGE_list_ptr->get_and_step();

      DLIList<FACE*> tmp_FACE_list;
      AQE->get_FACEs( EDGE_ptr, tmp_FACE_list );

      if( tmp_FACE_list.size() != 2 )
      {
        PRINT_ERROR( "Unexpected topology during %s curve function\n", command_name );
        // Delete all BODIEs, since we failed
        while( thickened_BODY_list.size() )
          AQE->delete_ACIS_BODY(thickened_BODY_list.pop(),CUBIT_TRUE);
        // Clean up memory allocated for the output lists
        while( output_EDGE_lists.size() ) delete output_EDGE_lists.pop();
        while( output_FACE_lists.size() ) delete output_FACE_lists.pop();
        while( conjugate_FACE_lists.size() ) delete conjugate_FACE_lists.pop();
        return CUBIT_FAILURE;
      }

      FACE *FACE_ptr1 = tmp_FACE_list.get_and_step();
      FACE *FACE_ptr2 = tmp_FACE_list.get();

      if( FACE_list_ptr->is_in_list( FACE_ptr1 ) )
        conjugate_FACE_list_ptr->append( FACE_ptr2 );
      else
        conjugate_FACE_list_ptr->append( FACE_ptr1 );
    }
  }

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : get_thickened_BODIES_of_VERTICES
// Member Type: PRIVATE
// Description: Get thickened BODIES from an input list of VERTICEs.  This
//              function can be called multiple times on the same list of
//              VERTICEs - each time, the VERTICEs from the common body that is
//              thickened are removed from the input EDGE list.  These
//              thickened BODIES are copied from the sheets the VERTICEs are
//              attached to (multiple BODIES can be returned because the
//              thickened BODIES must be nonmanifold).  For each BODY returned,
//              also output a list of the FACEs (corresponding to the original
//              FACEs of the sheet body), the VERTICEs (corresponding to the
//              original VERTICEs), and a list of EDGEs created by the
//              thickening process (from sweeping the original VERTEX - note it
//              is possible that no EDGE was created from a VERTEX, in which
//              case a NULL value will exist in the list).
// Author     : Steve Storm
// Date       : 02/10/06
//=============================================================================
CubitStatus
AcisTweakTool::get_thickened_BODIES_of_VERTICES( const char *command_name,
                                                 DLIList<VERTEX*> &input_VERTEX_list,
                                                 BODY *&common_BODY_ptr,
                                                 DLIList<BODY*> &thickened_BODY_list,
                                                 DLIList<DLIList<FACE*>*> &output_FACE_lists,
                                                 DLIList<DLIList<VERTEX*>*> &output_VERTEX_lists,
                                                 DLIList<DLIList<EDGE*>*> &output_EDGE_lists,
                                                 double thickness )
{
  assert( thickened_BODY_list.size() == 0 );

  int i, j;

  VERTEX *VERTEX_ptr;
  FACE *FACE_ptr;

  outcome result;

  // Pull VERTICEs out of the input VERTEX list that are from a common BODY
  DLIList<VERTEX*> common_VERTEX_list;
  BODY *BODY_ptr;
  common_BODY_ptr = NULL;
  input_VERTEX_list.reset();
  for( i=input_VERTEX_list.size(); i--; )
  {
    VERTEX_ptr = input_VERTEX_list.get();
    BODY_ptr = AQE->get_BODY_of_ENTITY( VERTEX_ptr );

    if( common_BODY_ptr == NULL )
      common_BODY_ptr = BODY_ptr;

    if( common_BODY_ptr == BODY_ptr )
    {
      common_VERTEX_list.append( VERTEX_ptr );
      input_VERTEX_list.change_to( NULL );
    }

    input_VERTEX_list.step();
  }
  input_VERTEX_list.remove_all_with_value( NULL );

  // Get all FACEs attached to the VERTICEs
  DLIList<FACE*> attached_FACE_list;
  common_VERTEX_list.reset();
  for( i=common_VERTEX_list.size(); i--; )
  {
    VERTEX_ptr = common_VERTEX_list.get_and_step();
    DLIList<FACE*> tmp_FACE_list;
    AQE->get_FACEs( VERTEX_ptr, tmp_FACE_list );

    if( tmp_FACE_list.size() == 0 )
    {
      PRINT_ERROR( "Cannot %s vertices that are free\n",
        command_name );
      return CUBIT_FAILURE;
    }

    if( tmp_FACE_list.size() != 1 )
    {
      PRINT_ERROR( "Can only %s vertices attached to one surface for sheet bodies\n",
        command_name );
      return CUBIT_FAILURE;
    }

    attached_FACE_list.append_unique( tmp_FACE_list.get() );
  }

  // Make sure that all of the FACEs are on a sheet body
  for( i=attached_FACE_list.size(); i--; )
  {
    FACE_ptr = attached_FACE_list.get_and_step();
    if( FACE_ptr->sides() != DOUBLE_SIDED )
    {
      PRINT_ERROR("Cannot %s vertices that are not on sheet bodies\n",
        command_name );
      return CUBIT_FAILURE;
    }
  }

  // Copy these FACEs off into a new body
  // Copy the BODY
  BODY *copied_BODY_ptr;
  if( copy_FACES_from_BODY( attached_FACE_list, copied_BODY_ptr ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  // Make the BODY nonmanifold, as thicken can't handle manifold bodies.  Note
  // the api_unstitch_nonmani api will remove the cubit owner from the resultant
  // BODY, LUMP and SHELL - we manually put them back on, for proper updating
  // of the Cubit entities (luckily it leaves them on the FACEs, CURVEs, etc..
  // Also note the api destroys the input BODY.
  AcisBridge *ab_body_ptr;
  AcisBridge *ab_lump_ptr;
  AcisBridge *ab_shell_ptr;
  get_owner_attribs( copied_BODY_ptr, ab_body_ptr, ab_lump_ptr, ab_shell_ptr );

  BODY* lumps  = NULL;
  BODY* sheets = NULL;
  BODY* lamina = NULL;
  BODY* wires  = NULL;
  result = api_unstitch_nonmani( copied_BODY_ptr, lumps, sheets, lamina, wires );
  if (!result.ok())
  {
    AQE->ACIS_API_error(result);
    api_delent( copied_BODY_ptr );
    return CUBIT_FAILURE;
  }

  // Make sure cubit owner attributes remain on resultant BODY, LUMP, SHELL
  reset_owner_attribs( sheets, ab_body_ptr, ab_lump_ptr, ab_shell_ptr );

  // Separate bodies (thicken won't work if more than one vol per body)
  BODY **separated_BODY_array = NULL;
  int n_body = 0;

  // Note first BODY will always be the input BODY
  result = api_separate_body( sheets, n_body, separated_BODY_array );
  if (!result.ok())
  {
    AQE->ACIS_API_error(result);
    api_delent( sheets );
    return CUBIT_FAILURE;
  }

  // These BODIES can now be thickened
  for( i=0; i<n_body; i++ )
    thickened_BODY_list.append(separated_BODY_array[i]);

  // Loop on separated BODIES
  thickened_BODY_list.reset();
  for( i=thickened_BODY_list.size(); i--; )
  {
    BODY_ptr = thickened_BODY_list.get_and_step();

    // Make sure cubit owner attributes remain on resultant BODY, LUMP, SHELL
    reset_owner_attribs( BODY_ptr, ab_body_ptr, ab_lump_ptr, ab_shell_ptr );

    // Keep track of the original EDGEs in the sheet body
    DLIList<EDGE*> original_EDGE_list;
    AQE->get_EDGEs( BODY_ptr, original_EDGE_list );

    // The FACEs in this BODY can efficiently be added to the output list here
    DLIList<FACE*> *FACE_list_ptr = new DLIList<FACE*>;
    output_FACE_lists.append( FACE_list_ptr );
    AQE->get_FACEs( BODY_ptr, *FACE_list_ptr );

    // Thicken the BODY
    if( thicken_BODY( BODY_ptr, thickness ) == CUBIT_FAILURE )
    {
      // Delete all BODIEs, since we failed
      while( thickened_BODY_list.size() )
        api_delent( thickened_BODY_list.pop() );
      // Clean up memory allocated for the output lists
      while( output_FACE_lists.size() ) delete output_FACE_lists.pop();
      while( output_VERTEX_lists.size() ) delete output_VERTEX_lists.pop();
      while( output_EDGE_lists.size() ) delete output_EDGE_lists.pop();
      return CUBIT_FAILURE;
    }

    // Make sure cubit owner attributes remain on resultant BODY, LUMP, SHELL
    reset_owner_attribs( BODY_ptr, ab_body_ptr, ab_lump_ptr, ab_shell_ptr );

    // Now find VERTICEs and output EDGEs for the thickened BODY
    DLIList<VERTEX*> *VERTEX_list_ptr = new DLIList<VERTEX*>;
    output_VERTEX_lists.append( VERTEX_list_ptr );

    DLIList<VERTEX*> tmp_VERTEX_list;
    AQE->get_VERTICEs( BODY_ptr, tmp_VERTEX_list );
    VERTEX *tmp_VERTEX_ptr;
    for( j=common_VERTEX_list.size(); j--; )
    {
      VERTEX_ptr = common_VERTEX_list.get_and_step();

      tmp_VERTEX_ptr = find_corresponding_VERTEX( VERTEX_ptr, tmp_VERTEX_list );
      if( tmp_VERTEX_ptr )
        VERTEX_list_ptr->append( tmp_VERTEX_ptr );
    }

    if( !VERTEX_list_ptr->size() )
    {
      PRINT_ERROR( "Internal error - please report.\n" );
      // Delete all BODIEs, since we failed
      while( thickened_BODY_list.size() )
        api_delent( thickened_BODY_list.pop() );
      // Clean up memory allocated for the output lists
      while( output_FACE_lists.size() ) delete output_FACE_lists.pop();
      while( output_VERTEX_lists.size() ) delete output_VERTEX_lists.pop();
      while( output_EDGE_lists.size() ) delete output_EDGE_lists.pop();
      return CUBIT_FAILURE;
    }

    // Find the "conjugate" EDGEs (the EDGEs that can be chamfered or filleted
    // by the calling function).  These are the EDGEs that were created by
    // the sweeping the VERTICEs during the thicken operation.

    DLIList<EDGE*> *EDGE_list_ptr = new DLIList<EDGE*>;
    output_EDGE_lists.append( EDGE_list_ptr );

    VERTEX_list_ptr->reset();
    for( j=VERTEX_list_ptr->size(); j--; )
    {
      VERTEX_ptr = VERTEX_list_ptr->get_and_step();

      DLIList<EDGE*> tmp_EDGE_list;
      AQE->get_EDGEs( VERTEX_ptr, tmp_EDGE_list );

      int k;
      EDGE *tmp_EDGE_ptr;
      int found = 0;
      for( k=tmp_EDGE_list.size(); k--; )
      {
        tmp_EDGE_ptr = tmp_EDGE_list.get_and_step();
        if( !original_EDGE_list.is_in_list( tmp_EDGE_ptr ) )
        {
          EDGE_list_ptr->append( tmp_EDGE_ptr );
          found = 1;
          break;
        }
      }
      if( !found )
        EDGE_list_ptr->append( NULL );
    }
  }

  return CUBIT_SUCCESS;
}

CubitStatus
AcisTweakTool::prep_for_surface_swap( BODY *thickened_BODY_ptr,
                                      BODY *copied_BODY_ptr,
                                      DLIList<AcisBridge*> &owner_FACE_list )
{
  // Pull original FACEs off and replace with extended FACEs

  // Remove all but the extended FACEs from thickened_BODY_ptr
  DLIList<FACE*> tweaked_FACE_list;
  get_corresponding_FACE_list( owner_FACE_list, thickened_BODY_ptr, tweaked_FACE_list );
  if( remove_FACES_from_BODY_except( thickened_BODY_ptr, tweaked_FACE_list )
    == CUBIT_FAILURE )
  {
    return CUBIT_FAILURE;
  }

  // Now remove the extended surfaces from the copied_input_BODY_ptr,
  // to be replaced by the output_FACE_list_ptr
  tweaked_FACE_list.clean_out();
  get_corresponding_FACE_list( owner_FACE_list, copied_BODY_ptr, tweaked_FACE_list );
  if( remove_FACES_from_BODY( copied_BODY_ptr, tweaked_FACE_list )
    == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  return CUBIT_SUCCESS;
}

CubitStatus
AcisTweakTool::unite_BODIES( BODY *copied_input_BODY_ptr,
                             DLIList<BODY*> &thickened_BODY_list,
                             BODY *&master )
{
  int i;

  // Unite the thickened BODIEs back into the copied_input_BODY_ptr
  outcome result;
  thickened_BODY_list.reset();
  BODY *BODY_ptr;
  master = copied_input_BODY_ptr;
  for( i=thickened_BODY_list.size(); i--; )
  {
    BODY_ptr = thickened_BODY_list.get();

    // Do the union of the master and the BODY_ptr.
    // If this is successful, the result is master and
    //   BODY_ptr will be deleted
    result = api_boolean( BODY_ptr, master, NONREG_UNION );
    if( !result.ok() || (!master) )
    {
      AQE->ACIS_API_error(result);
      while( thickened_BODY_list.size() )
        AQE->delete_ACIS_BODY(thickened_BODY_list.pop(),CUBIT_TRUE);
      if (master != NULL) AQE->delete_ACIS_BODY(master,CUBIT_TRUE);
      return CUBIT_FAILURE;
    }

    thickened_BODY_list.remove();
  }

  return CUBIT_SUCCESS;
}

CubitStatus
AcisTweakTool::get_owner_attribs( BODY *BODY_ptr, AcisBridge *&ab_body_ptr,
     AcisBridge *&ab_lump_ptr, AcisBridge *&ab_shell_ptr)
{
  ab_body_ptr = ATTRIB_CUBIT_OWNER::cubit_owner( BODY_ptr );

  ENTITY_LIST lump_list;
  api_get_lumps( BODY_ptr, lump_list);
  ab_lump_ptr = NULL;
  if( lump_list.count() == 1 )
    ab_lump_ptr = ATTRIB_CUBIT_OWNER::cubit_owner( lump_list[0] );

  ENTITY_LIST shell_list;
  api_get_shells( BODY_ptr, shell_list);
  ab_shell_ptr = NULL;
  if( shell_list.count() == 1 )
    ab_shell_ptr = ATTRIB_CUBIT_OWNER::cubit_owner( shell_list[0] );

  return CUBIT_SUCCESS;
}

CubitStatus
AcisTweakTool::reset_owner_attribs( BODY *BODY_ptr,
                                    AcisBridge *ab_body_ptr,
                                    AcisBridge *ab_lump_ptr,
                                    AcisBridge *ab_shell_ptr)
{
  // Check BODY
  if( !ATTRIB_CUBIT_OWNER::cubit_owner(BODY_ptr) )
  {
    if( ab_body_ptr )
      ATTRIB_CUBIT_OWNER::set_cubit_owner( BODY_ptr, ab_body_ptr );
  }

  // Check LUMP
  ENTITY_LIST lump_list;
  api_get_lumps( BODY_ptr, lump_list);
  if( lump_list.count() == 1 )
  {
    if( !ATTRIB_CUBIT_OWNER::cubit_owner(lump_list[0]) )
    {
      if( ab_lump_ptr )
        ATTRIB_CUBIT_OWNER::set_cubit_owner( lump_list[0], ab_lump_ptr );
    }
  }

  // Check SHELL
  ENTITY_LIST shell_list;
  api_get_shells( BODY_ptr, shell_list);
  if( shell_list.count() == 1 )
  {
    if( !ATTRIB_CUBIT_OWNER::cubit_owner(shell_list[0]) )
    {
      if( ab_shell_ptr )
        ATTRIB_CUBIT_OWNER::set_cubit_owner( shell_list[0], ab_shell_ptr );
    }
  }

  return CUBIT_SUCCESS;
}

CubitStatus
AcisTweakTool::sort_points_by_body_type( DLIList<Point*> &point_list,
                                         DLIList<Point*> &solid_points,
                                         DLIList<Point*> &sheet_points )
{
  int i, j;
  Point *point_ptr;
  PointACIS *point_acis_ptr;
  VERTEX *VERTEX_ptr;
  FACE *FACE_ptr;
  CubitBoolean sheet_face, solid_face;
  point_list.reset();
  for( i=point_list.size(); i--; )
  {
    point_ptr = point_list.get_and_step();

    point_acis_ptr = dynamic_cast<PointACIS*>(point_ptr);
    if( point_acis_ptr == NULL )
    {
      PRINT_ERROR( "Non-ACIS point encountered.\n" );
      return CUBIT_FAILURE;
    }

    VERTEX_ptr = point_acis_ptr->get_VERTEX_ptr();

    // Retrieve FACEs attached to this VERTEX
    DLIList<FACE*> FACE_list;
    AcisQueryEngine::instance()->get_FACEs( (ENTITY*)VERTEX_ptr, FACE_list );

    if( FACE_list.size() == 0 )
    {
      PRINT_ERROR( "Vertex found not attached to any surfaces.\n" );
      return CUBIT_FAILURE;
    }

    // Check all FACEs
    sheet_face = CUBIT_FALSE;
    solid_face = CUBIT_FALSE;

    for( j=FACE_list.size(); j--; )
    {
      FACE_ptr = FACE_list.get_and_step();

      if( FACE_ptr->sides() == DOUBLE_SIDED )
        sheet_face = CUBIT_TRUE;
      else
        solid_face = CUBIT_TRUE;
    }

    // Error if attached to both solid and sheet
    if( sheet_face == CUBIT_TRUE && solid_face == CUBIT_TRUE )
    {
      PRINT_ERROR( "Encountered vertex attached to both a sheet and a solid.\n" );
      return CUBIT_FAILURE;
    }

    if( sheet_face == CUBIT_TRUE )
      sheet_points.append( point_ptr );
    else
      solid_points.append( point_ptr );
  }

  return CUBIT_SUCCESS;
}

CubitStatus
AcisTweakTool::tweak_chamfer_solid( DLIList<Point*> &point_list,
                                    double radius,
                                    DLIList<BodySM*> &new_bodysm_list,
                                    CubitBoolean keep_old_body,
                                    CubitBoolean preview )
{
  if( point_list.size() == 0 )
    return CUBIT_SUCCESS;

  outcome result;

  BodySM *body_ptr;

  BODY *BODY_ptr;
  BODY *copied_BODY_ptr;

  int delete_attribs =
    (GeometryModifyTool::instance()->get_new_ids() || keep_old_body);

  // Copy the incoming point_list since we will be pulling points out of it.
  DLIList<PointACIS*> copied_point_list(point_list.size());
  CAST_LIST( point_list, copied_point_list, PointACIS );
  if (point_list.size() != copied_point_list.size())
  {
    PRINT_ERROR("Non-ACIS vertices encountered\n");
    return CUBIT_FAILURE;
  }

  copied_point_list.reset();
  while( copied_point_list.size() )
  {
    DLIList<VERTEX*> VERTEX_list;
    if( AME->get_copied_VERTICES_of_body( copied_point_list, VERTEX_list,
      copied_BODY_ptr ) == CUBIT_FAILURE )
    {
      // Return success if any bodies were created
      return new_bodysm_list.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
    }

    // Get original Body and BODY
    body_ptr = AQE->get_body_sm_of_ENTITY( copied_BODY_ptr );
    BODY_ptr = dynamic_cast<BodyACIS*>(body_ptr)->get_BODY_ptr();

    // Now, blend the edges on this body
    if( chamfer_vertices( VERTEX_list, radius ) == CUBIT_FAILURE )
    {
      api_delent(copied_BODY_ptr);
      // Return success if any bodies were created
      return new_bodysm_list.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
    }

    if( !preview )
    {
      // If we've made it this far, the copied_BODY has been
      // modified and we can update it in CUBIT

      // Now cleanout the owner attributes from the copied BODY, if required
      if( delete_attribs )
        AQE->remove_cubit_owner_attrib_in_BODY(copied_BODY_ptr);

      BodySM* new_body_ptr = AME->get_new_Body( body_ptr, BODY_ptr, copied_BODY_ptr,
        keep_old_body );

      if (new_body_ptr)
        new_bodysm_list.append( new_body_ptr );
    }
    else
    {
      GfxPreview::clear();

      ENTITY_LIST face_list;
      api_get_faces( copied_BODY_ptr, face_list);
      AcisBridge *ab_face_ptr = NULL;
      int i;
      for( i=0; i<face_list.count(); i++ )
      {
        ab_face_ptr = ATTRIB_CUBIT_OWNER::cubit_owner( face_list[i] );
        if( !ab_face_ptr )
        {
          // Draw this face
          AcisDrawTool::instance()->draw_FACE( (FACE*)face_list[i], CUBIT_BLUE );
        }
      }

      api_delent(copied_BODY_ptr);
    }
  }

  if( preview )
    GfxPreview::flush();

  return CUBIT_SUCCESS;
}

CubitStatus
AcisTweakTool::tweak_chamfer_fillet_sheet( DLIList<Point*> &input_point_list,
                                           double radius,
                                           int type,
                                           DLIList<BodySM*> &new_bodysm_list,
                                           CubitBoolean keep_old_bodies,
                                           CubitBoolean preview )
{
  if( input_point_list.size() == 0 )
    return CUBIT_SUCCESS;

  int i, j;
  EDGE *EDGE_ptr;
  outcome result;

  bool delete_attribs =
    (GeometryModifyTool::instance()->get_new_ids() || keep_old_bodies);

  // Copy the input point list as we will be removing points from it
  DLIList<Point*> copied_input_point_list = input_point_list;

  char type_name[8];
  if( type == 1 )
    strcpy( type_name, "chamfer" );
  else
    strcpy( type_name, "fillet" );

  // Get VERTICEs to chamfer or fillet
  DLIList<VERTEX*> VERTEX_list;
  if( get_VERTICEs( input_point_list, VERTEX_list ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  while( VERTEX_list.size() )
  {
    BODY *input_BODY_ptr;
    BodySM *input_bodysm_ptr;
    DLIList<BODY*> thickened_BODY_list;
    DLIList<DLIList<FACE*>*> output_FACE_lists;
    DLIList<DLIList<VERTEX*>*> output_VERTEX_lists;
    DLIList<DLIList<EDGE*>*> output_EDGE_lists;

    // Get thickened BODIES of VERTICES from a common BODY
    if( get_thickened_BODIES_of_VERTICES( type_name, VERTEX_list,
      input_BODY_ptr, thickened_BODY_list, output_FACE_lists,
      output_VERTEX_lists, output_EDGE_lists ) == CUBIT_FAILURE )
    {
      // Return success if any bodies were created
        return new_bodysm_list.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
    }

    input_bodysm_ptr = AQE->get_body_sm_of_ENTITY( input_BODY_ptr );

    // Tweak the surfaces on the output BODIEs, pull the tweaked surfaces off,
    // and replace them in the original BODY
    // First copy the BODY
    BODY *copied_input_BODY_ptr = AME->copy_BODY(input_BODY_ptr, CUBIT_FALSE);

    BODY *thickened_BODY_ptr;
    DLIList<FACE*> *output_FACE_list_ptr;
    DLIList<VERTEX*> *output_VERTEX_list_ptr;
    DLIList<EDGE*> *output_EDGE_list_ptr;
    thickened_BODY_list.reset();
    output_FACE_lists.reset();
    output_VERTEX_lists.reset();
    output_EDGE_lists.reset();
    for( i=thickened_BODY_list.size(); i--; )
    {
      thickened_BODY_ptr = thickened_BODY_list.get_and_step();
      output_FACE_list_ptr = output_FACE_lists.get_and_step();
      output_VERTEX_list_ptr = output_VERTEX_lists.get_and_step();
      output_EDGE_list_ptr = output_EDGE_lists.get_and_step();

      // Note - the output EDGE list can have NULL entries in it, as thicken
      // doesn't always create a curve on the side of the thickened body at
      // each vertex.  If there is a NULL entry, it means that there isn't
      // a discernable topology change in the body there anyway, which means
      // we wouldn't be able to chamfer it anyway.  For now just ignore.

      // Keep track of the FACEs of interest, via their Cubit owners, which
      // will survive throughout the operation
      DLIList<AcisBridge*> owner_list;
      get_owner_list( *output_FACE_list_ptr, owner_list );

      // Error check
      int cnt = 0;
      for( j=output_EDGE_list_ptr->size(); j--; )
      {
        EDGE_ptr = output_EDGE_list_ptr->get_and_step();
        if( EDGE_ptr )
          cnt++;
      }

      if( cnt == 0 )
      {
        PRINT_ERROR( "Unable to %s vertices.\n", type_name );
        while( output_FACE_lists.size() ) delete output_FACE_lists.pop();
        while( output_VERTEX_lists.size() ) delete output_VERTEX_lists.pop();
        while( output_EDGE_lists.size() ) delete output_EDGE_lists.pop();
        while( thickened_BODY_list.size() )
          api_delent( thickened_BODY_list.pop() );
        api_delent( copied_input_BODY_ptr );
        // Return success if any bodies were created
        return new_bodysm_list.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
      }

      // Now, blend the edges on this thickened body
      CubitStatus status;
      if( type == 1 )
        status = chamfer_edges( *output_EDGE_list_ptr, radius );
      else
        status = blend_edges( *output_EDGE_list_ptr, radius );


      if( status == CUBIT_FAILURE )
      {
        while( output_FACE_lists.size() ) delete output_FACE_lists.pop();
        while( output_VERTEX_lists.size() ) delete output_VERTEX_lists.pop();
        while( output_EDGE_lists.size() ) delete output_EDGE_lists.pop();
        while( thickened_BODY_list.size() )
          api_delent( thickened_BODY_list.pop() );
        api_delent( copied_input_BODY_ptr );
        // Return success if any bodies were created
        return new_bodysm_list.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
      }

      // Get bodies ready to swap new surfaces for old
      if( prep_for_surface_swap( thickened_BODY_ptr, copied_input_BODY_ptr,
        owner_list ) == CUBIT_FAILURE )
      {
        while( output_FACE_lists.size() ) delete output_FACE_lists.pop();
        while( output_VERTEX_lists.size() ) delete output_VERTEX_lists.pop();
        while( output_EDGE_lists.size() ) delete output_EDGE_lists.pop();
        while( thickened_BODY_list.size() )
          api_delent( thickened_BODY_list.pop() );
        api_delent( copied_input_BODY_ptr );
        // Return success if any bodies were created
        return new_bodysm_list.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
      }

    } // End loop on thickened (separated) BODIES

    // Free memory
    while( output_FACE_lists.size() ) delete output_FACE_lists.pop();
    while( output_VERTEX_lists.size() ) delete output_VERTEX_lists.pop();
    while( output_EDGE_lists.size() ) delete output_EDGE_lists.pop();

    // Unite the thickened BODIEs back into the copied_input_BODY_ptr
    BODY *master;
    if( unite_BODIES( copied_input_BODY_ptr, thickened_BODY_list, master ) == CUBIT_FAILURE )
    {
      // If failure, the unite_BODIES function cleaned up the memory
      // Return success if any bodies were created
      return new_bodysm_list.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
    }

    // This BODY is done (whew!)

    if( !preview )
    {
      // Now cleanout the owner attributes from the copied BODY, if required
      if( delete_attribs )
        AQE->remove_cubit_owner_attrib_in_BODY(master);

      BodySM *new_body = AME->get_new_Body( input_bodysm_ptr, input_BODY_ptr,
        master, keep_old_bodies );

      if( new_body )
        new_bodysm_list.append( new_body );
    }
    else
    {
      GfxPreview::clear();

      ENTITY_LIST edge_list;
      api_get_edges( master, edge_list );
      AcisBridge *ab_edge_ptr = NULL;
      for( i=0; i<edge_list.count(); i++ )
      {
        ab_edge_ptr = ATTRIB_CUBIT_OWNER::cubit_owner( edge_list[i] );
        if( !ab_edge_ptr )
        {
          // Draw this edge
          AcisDrawTool::instance()->draw_EDGE( (EDGE*)edge_list[i], CUBIT_BLUE );
        }
      }

      api_delent(master);
    }
  }

  if( preview )
    GfxPreview::flush();

  return CUBIT_SUCCESS;
}

CubitStatus
AcisTweakTool::tweak_chamfer_solid( Point* point_ptr,
                                    double r1,
                                    Curve *c1,
                                    double r2,
                                    Curve *c2,
                                    double r3,
                                    Curve *c3,
                                    BodySM *&new_bodysm_ptr,
                                    CubitBoolean keep_old_body,
                                    CubitBoolean preview )
{
  outcome result;

  BodySM *body_ptr;
  BODY *BODY_ptr;
  BODY *copied_BODY_ptr;

  int delete_attribs =
    (GeometryModifyTool::instance()->get_new_ids() || keep_old_body);

  // Make sure at least curve 1 supplied
  if( c1 == NULL )
  {
    PRINT_ERROR( "Curve not supplied for radius.\n" );
    return CUBIT_FAILURE;
  }

  // Check for duplicate input curves.  Note c2 and c3 could be NULL, so we
  // don't want to error in that case.
  if( c1==c2 || c1==c3 || (c2 && c2==c3) )
  {
    PRINT_ERROR( "Same curve cannot be specified for multiple radii.\n" );
    return CUBIT_FAILURE;
  }

  PointACIS *point_acis_ptr = CAST_TO( point_ptr, PointACIS );
  if( point_acis_ptr == NULL )
  {
    PRINT_ERROR( "Non-ACIS vertex encountered.\n" );
    return CUBIT_FAILURE;
  }

  DLIList<PointACIS*> copied_point_list(1);
  copied_point_list.append( point_acis_ptr );

  DLIList<VERTEX*> VERTEX_list;
  if( AME->get_copied_VERTICES_of_body( copied_point_list, VERTEX_list,
    copied_BODY_ptr ) == CUBIT_FAILURE )
  {
    return CUBIT_FAILURE;
  }

  // Get original Body and BODY
  body_ptr = AQE->get_body_sm_of_ENTITY( copied_BODY_ptr );
  BODY_ptr = dynamic_cast<BodyACIS*>(body_ptr)->get_BODY_ptr();

  if( VERTEX_list.size() != 1 )
  {
    PRINT_ERROR( "internal error - please report.\n" );
    api_delent(copied_BODY_ptr);
    return CUBIT_FAILURE;
  }

  VERTEX *VERTEX_ptr = VERTEX_list.get();

  // Get EDGEs attached to this vertex
  DLIList<EDGE*> EDGE_list;
  AQE->get_EDGEs( VERTEX_ptr, EDGE_list );
  if( EDGE_list.size() != 3 )
  {
    PRINT_ERROR( "found more than 3 curves attached to vertex.\n" );
    api_delent(copied_BODY_ptr);
    return CUBIT_FAILURE;
  }

  EDGE *ref_EDGE_ptr, *EDGE_ptr1, *EDGE_ptr2, *EDGE_ptr3;
  CurveACIS *curve_acis_ptr = CAST_TO( c1, CurveACIS );
  ref_EDGE_ptr = curve_acis_ptr->get_EDGE_ptr();
  EDGE_ptr1 = find_corresponding_EDGE( ref_EDGE_ptr, EDGE_list );
  if( EDGE_ptr1 == NULL )
  {
    PRINT_ERROR( "Supplied curve not found in volume.\n" );
    api_delent(copied_BODY_ptr);
    return CUBIT_FAILURE;
  }

  if( EDGE_list.move_to( EDGE_ptr1 ) == CUBIT_FALSE )
  {
    PRINT_ERROR( "First supplied curve not attached to vertex.\n" );
    api_delent(copied_BODY_ptr);
    return CUBIT_FAILURE;
  }

  // Remove this EDGE from the list
  EDGE_list.remove();

  // Handle cases where not all curves or radii were supplied
  if( c2 )
  {
    curve_acis_ptr = CAST_TO( c2, CurveACIS );
    ref_EDGE_ptr = curve_acis_ptr->get_EDGE_ptr();
    EDGE_ptr2 = find_corresponding_EDGE( ref_EDGE_ptr, EDGE_list );
    if( EDGE_ptr2 == NULL )
    {
      PRINT_ERROR( "Supplied curve not found in volume.\n" );
      api_delent(copied_BODY_ptr);
      return CUBIT_FAILURE;
    }

    if( EDGE_list.move_to( EDGE_ptr2 ) == CUBIT_FALSE )
    {
      PRINT_ERROR( "Second supplied curve not attached to vertex.\n" );
      api_delent(copied_BODY_ptr);
      return CUBIT_FAILURE;
    }

    EDGE_list.remove();

    if( c3 )
    {
      curve_acis_ptr = CAST_TO( c3, CurveACIS );
      ref_EDGE_ptr = curve_acis_ptr->get_EDGE_ptr();
      EDGE_ptr3 = find_corresponding_EDGE( ref_EDGE_ptr, EDGE_list );
      if( EDGE_ptr3 == NULL )
      {
        PRINT_ERROR( "Supplied curve not found in volume.\n" );
        api_delent(copied_BODY_ptr);
        return CUBIT_FAILURE;
      }

      if( EDGE_list.move_to( EDGE_ptr3 ) == CUBIT_FALSE )
      {
        PRINT_ERROR( "Third supplied curve not attached to vertex.\n" );
        api_delent(copied_BODY_ptr);
        return CUBIT_FAILURE;
      }
    }
    else
      EDGE_ptr3 = EDGE_list.get();
  }
  else
  {
    // Random results for the user - ok
    EDGE_ptr2 = EDGE_list.get_and_step();
    EDGE_ptr3 = EDGE_list.get();
  }

  if( r2 <= 0.0 )
    r2 = r1;
  if( r3 <= 0.0 )
    r3 = r2;

  // Setup tweak attributes so we can preserve Cubit owners
  DLIList<FACE*> pre_FACE_list;
  DLIList<EDGE*> pre_EDGE_list;
  DLIList<VERTEX*> pre_VERTEX_list;
  DLIList<AcisBridge*> ab_FACE_list, ab_EDGE_list, ab_VERTEX_list;
  assign_tweak_attribs( copied_BODY_ptr, "tweak", pre_FACE_list, ab_FACE_list,
    pre_EDGE_list, ab_EDGE_list, pre_VERTEX_list, ab_VERTEX_list );

  // Now, blend the edges on this body
  result = api_chamfer_vertex( VERTEX_ptr, r1, EDGE_ptr1, r2, EDGE_ptr2, r3,
    EDGE_ptr3 );

  if( !result.ok() )
  {
    AQE->ACIS_API_error(result);
    api_delent(copied_BODY_ptr);
    return CUBIT_FAILURE;
  }

  // Replace Cubit owners
  reassign_cubit_owners_from_tweak_attribs( copied_BODY_ptr, "tweak",
    pre_FACE_list, ab_FACE_list, pre_EDGE_list, ab_EDGE_list, pre_VERTEX_list,
    ab_VERTEX_list );

  // Remove tweak attributes
  remove_named_attribs( copied_BODY_ptr, "tweak" );

  if( !preview )
  {

  // If we've made it this far, the copied_BODY has been
  // modified and we can update it in CUBIT

  // Now cleanout the owner attributes from the copied BODY, if required
  if( delete_attribs )
    AQE->remove_cubit_owner_attrib_in_BODY(copied_BODY_ptr);

  new_bodysm_ptr = AME->get_new_Body( body_ptr, BODY_ptr, copied_BODY_ptr,
    keep_old_body );
  }
  else
  {
    GfxPreview::clear();

    ENTITY_LIST face_list;
    api_get_faces( copied_BODY_ptr, face_list );
    AcisBridge *ab_face_ptr = NULL;
    int i;
    for( i=0; i<face_list.count(); i++ )
    {
      ab_face_ptr = ATTRIB_CUBIT_OWNER::cubit_owner( face_list[i] );
      if( !ab_face_ptr )
      {
        // Draw this face
        AcisDrawTool::instance()->draw_FACE( (FACE*)face_list[i], CUBIT_BLUE );
      }
    }

    api_delent(copied_BODY_ptr);

    GfxPreview::flush();
  }

  return CUBIT_SUCCESS;
}

CubitStatus
AcisTweakTool::tweak_chamfer_sheet( Point* point_ptr,
                                    double r1,
                                    Curve *c1,
                                    double r2,
                                    Curve *c2,
                                    BodySM *&new_bodysm_ptr,
                                    CubitBoolean keep_old_body,
                                    CubitBoolean preview )
{
  int i, j;
  outcome result;
  EDGE *EDGE_ptr;

  int delete_attribs =
    (GeometryModifyTool::instance()->get_new_ids() || keep_old_body);

  // Make sure at least curve 1 supplied
  if( c1 == NULL )
  {
    PRINT_ERROR( "Curve not supplied for radius.\n" );
    return CUBIT_FAILURE;
  }

  // Check for duplicate input curves.  Note c2 and c3 could be NULL, so we
  // don't want to error in that case.
  if( c1==c2 )
  {
    PRINT_ERROR( "Same curve cannot be specified for both radii.\n" );
    return CUBIT_FAILURE;
  }

  VERTEX *VERTEX_ptr = AQE->get_VERTEX( point_ptr );
  if( !VERTEX_ptr )
  {
    PRINT_ERROR( "Chamfer vertex must be an ACIS vertex.\n" );
    return CUBIT_FAILURE;
  }

  DLIList<VERTEX*> copied_input_VERTEX_list(1);
  copied_input_VERTEX_list.append( VERTEX_ptr );

  BodySM *input_bodysm_ptr;
  BODY *input_BODY_ptr;
  DLIList<BODY*> thickened_BODY_list;
  DLIList<DLIList<FACE*>*> output_FACE_lists;
  DLIList<DLIList<VERTEX*>*> output_VERTEX_lists;
  DLIList<DLIList<EDGE*>*> output_EDGE_lists;
  if( get_thickened_BODIES_of_VERTICES( "chamfer", copied_input_VERTEX_list,
    input_BODY_ptr, thickened_BODY_list, output_FACE_lists,
    output_VERTEX_lists, output_EDGE_lists ) == CUBIT_FAILURE )
  {
    return CUBIT_FAILURE;
  }

  input_bodysm_ptr = AQE->get_body_sm_of_ENTITY( input_BODY_ptr );

  // Tweak the surfaces on the output BODIEs, pull the tweaked surfaces off,
  // and replace them in the original BODY
  // First copy the BODY
  BODY *copied_input_BODY_ptr = AME->copy_BODY(input_BODY_ptr, CUBIT_FALSE);

  BODY *thickened_BODY_ptr;
  DLIList<FACE*> *output_FACE_list_ptr;
  DLIList<VERTEX*> *output_VERTEX_list_ptr;
  DLIList<EDGE*> *output_EDGE_list_ptr;
  thickened_BODY_list.reset();
  output_FACE_lists.reset();
  output_VERTEX_lists.reset();
  output_EDGE_lists.reset();
  for( i=thickened_BODY_list.size(); i--; )
  {
    thickened_BODY_ptr = thickened_BODY_list.get_and_step();
    output_FACE_list_ptr = output_FACE_lists.get_and_step();
    output_VERTEX_list_ptr = output_VERTEX_lists.get_and_step();
    output_EDGE_list_ptr = output_EDGE_lists.get_and_step();

    // Note - the output EDGE list can have NULL entries in it, as thicken
    // doesn't always create a curve on the side of the thickened body at
    // each vertex.  If there is a NULL entry, it means that there isn't
    // a discernable topology change in the body there anyway, which means
    // we wouldn't be able to chamfer it anyway.  For now just ignore.

    // Keep track of the FACEs of interest, via their Cubit owners, which
    // will survive throughout the operation
    DLIList<AcisBridge*> owner_list;
    get_owner_list( *output_FACE_list_ptr, owner_list );

    // Error check
    int cnt = 0;
    for( j=output_EDGE_list_ptr->size(); j--; )
    {
      EDGE_ptr = output_EDGE_list_ptr->get_and_step();
      if( EDGE_ptr )
        cnt++;
    }

    if( cnt == 0 )
    {
      PRINT_ERROR( "Unable to chamfer the vertex.\n" );
      while( output_FACE_lists.size() ) delete output_FACE_lists.pop();
      while( output_VERTEX_lists.size() ) delete output_VERTEX_lists.pop();
      while( output_EDGE_lists.size() ) delete output_EDGE_lists.pop();
      while( thickened_BODY_list.size() )
        api_delent( thickened_BODY_list.pop() );
      api_delent( copied_input_BODY_ptr );
      return CUBIT_FAILURE;
    }

    // Check supplied curves
    VERTEX *ref_VERTEX_ptr = output_VERTEX_list_ptr->get();
    EDGE *conj_EDGE_ptr = output_EDGE_list_ptr->get();

    // Get EDGEs attached to this vertex
    DLIList<EDGE*> EDGE_list;
    AQE->get_EDGEs( ref_VERTEX_ptr, EDGE_list );
    if( EDGE_list.size() != 3 )
    {
      PRINT_ERROR( "found too many curves attached to vertex.\n" );
      while( output_FACE_lists.size() ) delete output_FACE_lists.pop();
      while( output_VERTEX_lists.size() ) delete output_VERTEX_lists.pop();
      while( output_EDGE_lists.size() ) delete output_EDGE_lists.pop();
      while( thickened_BODY_list.size() )
        api_delent( thickened_BODY_list.pop() );
      api_delent( copied_input_BODY_ptr );
      return CUBIT_FAILURE;
    }

    EDGE *ref_EDGE_ptr, *EDGE_ptr1, *EDGE_ptr2;
    CurveACIS *curve_acis_ptr = CAST_TO( c1, CurveACIS );
    ref_EDGE_ptr = curve_acis_ptr->get_EDGE_ptr();
    EDGE_ptr1 = find_corresponding_EDGE( ref_EDGE_ptr, EDGE_list );
    if( EDGE_ptr1 == NULL )
    {
      PRINT_ERROR( "Supplied curve not found in volume.\n" );
      while( output_FACE_lists.size() ) delete output_FACE_lists.pop();
      while( output_VERTEX_lists.size() ) delete output_VERTEX_lists.pop();
      while( output_EDGE_lists.size() ) delete output_EDGE_lists.pop();
      while( thickened_BODY_list.size() )
        api_delent( thickened_BODY_list.pop() );
      api_delent( copied_input_BODY_ptr );
      return CUBIT_FAILURE;
    }

    if( EDGE_list.move_to( EDGE_ptr1 ) == CUBIT_FALSE )
    {
      PRINT_ERROR( "First curve not attached to vertex.\n" );
      while( output_FACE_lists.size() ) delete output_FACE_lists.pop();
      while( output_VERTEX_lists.size() ) delete output_VERTEX_lists.pop();
      while( output_EDGE_lists.size() ) delete output_EDGE_lists.pop();
      while( thickened_BODY_list.size() )
        api_delent( thickened_BODY_list.pop() );
      api_delent( copied_input_BODY_ptr );
      return CUBIT_FAILURE;
    }

    // Remove this EDGE from the list
    EDGE_list.remove();

    // Handle cases where not all curves or radii were supplied
    if( c2 )
    {
      curve_acis_ptr = CAST_TO( c2, CurveACIS );
      ref_EDGE_ptr = curve_acis_ptr->get_EDGE_ptr();
      EDGE_ptr2 = find_corresponding_EDGE( ref_EDGE_ptr, EDGE_list );
      if( EDGE_ptr2 == NULL )
      {
        PRINT_ERROR( "Second curve not found in volume.\n" );
        while( output_FACE_lists.size() ) delete output_FACE_lists.pop();
        while( output_VERTEX_lists.size() ) delete output_VERTEX_lists.pop();
        while( output_EDGE_lists.size() ) delete output_EDGE_lists.pop();
        while( thickened_BODY_list.size() )
          api_delent( thickened_BODY_list.pop() );
        api_delent( copied_input_BODY_ptr );
        return CUBIT_FAILURE;
      }

      if( EDGE_list.move_to( EDGE_ptr2 ) == CUBIT_FALSE )
      {
        PRINT_ERROR( "Second curve not attached to vertex.\n" );
        while( output_FACE_lists.size() ) delete output_FACE_lists.pop();
        while( output_VERTEX_lists.size() ) delete output_VERTEX_lists.pop();
        while( output_EDGE_lists.size() ) delete output_EDGE_lists.pop();
        while( thickened_BODY_list.size() )
          api_delent( thickened_BODY_list.pop() );
        api_delent( copied_input_BODY_ptr );
        return CUBIT_FAILURE;
      }

      EDGE_list.remove();
    }
    else
    {
      EDGE_ptr2 = EDGE_list.get_and_step();
      if( EDGE_ptr2 == conj_EDGE_ptr )
        EDGE_ptr2 = EDGE_list.get();
    }

    // Determine whether r1 and r2 need to be swapped

    // Are we on the end of the conjugate curve?
    int on_conj_end = 0;
    if( conj_EDGE_ptr->end() == ref_VERTEX_ptr )
      on_conj_end = 1;

    FACE *FACE_ptr = output_FACE_list_ptr->get();
    LOOP *LOOP_ptr = FACE_ptr->loop();

    // Find CoEdge whose end is our vertex
    COEDGE *COEDGE_ptr = LOOP_ptr->start();
    VERTEX *end_VERTEX_ptr;
    while( COEDGE_ptr != NULL )
    {
      end_VERTEX_ptr = COEDGE_ptr->end();
      if( end_VERTEX_ptr == ref_VERTEX_ptr )
        break;

      COEDGE_ptr = COEDGE_ptr->next();
      if( COEDGE_ptr == LOOP_ptr->start() )
        break;
    }

    double tmp_r;
    EDGE_ptr = COEDGE_ptr->edge();

    if( (EDGE_ptr == EDGE_ptr2 && on_conj_end) ||
        (EDGE_ptr == EDGE_ptr1 && !on_conj_end) )
    {
      // Swap
      tmp_r = r1;
      r1 = r2;
      r2 = tmp_r;
    }

    // Now, blend the edges on this thickened body
    if( chamfer_edges( *output_EDGE_list_ptr, r1, r2 ) == CUBIT_FAILURE )
    {
      while( output_FACE_lists.size() ) delete output_FACE_lists.pop();
      while( output_VERTEX_lists.size() ) delete output_VERTEX_lists.pop();
      while( output_EDGE_lists.size() ) delete output_EDGE_lists.pop();
      while( thickened_BODY_list.size() )
        api_delent( thickened_BODY_list.pop() );
      api_delent( copied_input_BODY_ptr );
      return CUBIT_FAILURE;
    }

    // Get bodies ready to swap new surfaces for old
    if( prep_for_surface_swap( thickened_BODY_ptr, copied_input_BODY_ptr,
      owner_list ) == CUBIT_FAILURE )
    {
      while( output_FACE_lists.size() ) delete output_FACE_lists.pop();
      while( output_VERTEX_lists.size() ) delete output_VERTEX_lists.pop();
      while( output_EDGE_lists.size() ) delete output_EDGE_lists.pop();
      while( thickened_BODY_list.size() )
        api_delent( thickened_BODY_list.pop() );
      api_delent( input_BODY_ptr );
      return CUBIT_FAILURE;
    }

  } // End loop on thickened (separated) BODIES

  // Free memory
  while( output_FACE_lists.size() ) delete output_FACE_lists.pop();
  while( output_VERTEX_lists.size() ) delete output_VERTEX_lists.pop();
  while( output_EDGE_lists.size() ) delete output_EDGE_lists.pop();

  // Unite the thickened BODIEs back into the input_BODY_ptr
  BODY *master;
  if( unite_BODIES( input_BODY_ptr, thickened_BODY_list, master ) == CUBIT_FAILURE )
  {
    // If failure, the unite_BODIES function cleaned up the memory
    return CUBIT_FAILURE;
  }

  // This BODY is done (whew!)

  if( !preview )
  {
    // Now cleanout the owner attributes from the copied BODY, if required
    if( delete_attribs )
      AQE->remove_cubit_owner_attrib_in_BODY(master);

    new_bodysm_ptr = AME->get_new_Body( input_bodysm_ptr, input_BODY_ptr,
      master, keep_old_body );
  }
  else
  {
    ENTITY_LIST edge_list;
    api_get_edges( master, edge_list );
    AcisBridge *ab_edge_ptr = NULL;
    for( i=0; i<edge_list.count(); i++ )
    {
      ab_edge_ptr = ATTRIB_CUBIT_OWNER::cubit_owner( edge_list[i] );
      if( !ab_edge_ptr )
      {
        // Draw this edge
        AcisDrawTool::instance()->draw_EDGE( (EDGE*)edge_list[i], CUBIT_BLUE );
      }
    }

    api_delent(master);

    GfxPreview::flush();
  }

  return CUBIT_SUCCESS;
}

CubitStatus
AcisTweakTool::assign_tweak_attribs( BODY *BODY_ptr, const char *att_name,
     DLIList<FACE*> &FACE_list, DLIList<AcisBridge*> &ab_FACE_list,
     DLIList<EDGE*> &EDGE_list, DLIList<AcisBridge*> &ab_EDGE_list,
     DLIList<VERTEX*> &VERTEX_list, DLIList<AcisBridge*> &ab_VERTEX_list )
{
  // Assign named attribs to faces, edges and vertices so we can get back to
  // them after a fillet or chamfer.  These operations split these entities
  // to make room for the new surfaces, thus we can copy the owner att
  // during the split.

  int i;
  outcome result;

  // FACEs
  AQE->get_FACEs( BODY_ptr, FACE_list );
  FACE *FACE_ptr;
  FACE_list.reset();
  for( i=0; i<FACE_list.size(); i++ )
  {
    FACE_ptr = FACE_list.get_and_step();
    result = api_add_generic_named_attribute( FACE_ptr, att_name, i, SplitCopy );
    if( !result.ok() )
    {
      ab_FACE_list.append( NULL );
      continue;
    }
    ab_FACE_list.append( ATTRIB_CUBIT_OWNER::cubit_owner( FACE_ptr ) );
  }


  // EDGEs
  AQE->get_EDGEs( BODY_ptr, EDGE_list );
  EDGE *EDGE_ptr;
  EDGE_list.reset();
  for( i=0; i<EDGE_list.size(); i++ )
  {
    EDGE_ptr = EDGE_list.get_and_step();
    result = api_add_generic_named_attribute( EDGE_ptr, att_name, i, SplitCopy );
    if( !result.ok() )
    {
      ab_EDGE_list.append( NULL );
      continue;
    }
    ab_EDGE_list.append( ATTRIB_CUBIT_OWNER::cubit_owner( EDGE_ptr ) );
  }

  // VERTICEs
  AQE->get_VERTICEs( BODY_ptr, VERTEX_list );
  VERTEX *VERTEX_ptr;
  VERTEX_list.reset();
  for( i=0; i<VERTEX_list.size(); i++ )
  {
    VERTEX_ptr = VERTEX_list.get_and_step();
    result = api_add_generic_named_attribute( VERTEX_ptr, att_name, i, SplitCopy );
    if( !result.ok() )
    {
      ab_VERTEX_list.append( NULL );
      continue;
    }
    ab_VERTEX_list.append( ATTRIB_CUBIT_OWNER::cubit_owner( VERTEX_ptr ) );
  }

  return CUBIT_SUCCESS;
}

CubitStatus
AcisTweakTool::find_corresponding_entity_from_tweak_attrib( BODY *BODY_ptr,
                                                            const char *att_name,
                                                            const int ent_type,
                                                            int input_ent_id,
                                                            ENTITY *&output_ENTITY_ptr )
{
  output_ENTITY_ptr = NULL;
  ATTRIB_GEN_NAME *ret_att;

  outcome result;

  int i, ent_id;
  int found = 0;
  if( ent_type == FACE_TYPE )
  {
    DLIList<FACE*> FACE_list;
    AQE->get_FACEs( BODY_ptr, FACE_list );
    FACE *FACE_ptr;
    FACE_list.reset();
    for( i=FACE_list.size(); i--; )
    {
      FACE_ptr = FACE_list.get_and_step();
      result = api_find_named_attribute( FACE_ptr, att_name, ret_att );

      if( !result.ok() || !ret_att )
        continue;

      ent_id = ((ATTRIB_GEN_INTEGER *) ret_att)->value();

      if( ent_id == input_ent_id )
      {
        if( found == 1 )
        {
          // More than one entity with same att
          output_ENTITY_ptr = NULL;
          return CUBIT_FAILURE;
        }
        output_ENTITY_ptr = (ENTITY *)FACE_ptr;
        found = 1;
        break;
      }
    }
    if( found )
      return CUBIT_SUCCESS;
    else
      return CUBIT_FAILURE;
  }

  else if( ent_type == EDGE_TYPE )
  {
    DLIList<EDGE*> EDGE_list;
    AQE->get_EDGEs( BODY_ptr, EDGE_list );
    EDGE *EDGE_ptr;
    EDGE_list.reset();
    for( i=EDGE_list.size(); i--; )
    {
      EDGE_ptr = EDGE_list.get_and_step();

      result = api_find_named_attribute( EDGE_ptr, att_name, ret_att );

      if( !result.ok() || !ret_att )
        continue;

      ent_id = ((ATTRIB_GEN_INTEGER *) ret_att)->value();

      if( ent_id == input_ent_id )
      {
        if( found == 1 )
        {
          // More than one entity with same att
          output_ENTITY_ptr = NULL;
          return CUBIT_FAILURE;
        }
        output_ENTITY_ptr = (ENTITY *)EDGE_ptr;
        found = 1;
        break;
      }
    }
    if( found )
      return CUBIT_SUCCESS;
    else
      return CUBIT_FAILURE;
  }

  else if( ent_type == VERTEX_TYPE )
  {
    DLIList<VERTEX*> VERTEX_list;
    AQE->get_VERTICEs( BODY_ptr, VERTEX_list );
    VERTEX *VERTEX_ptr;
    VERTEX_list.reset();
    for( i=VERTEX_list.size(); i--; )
    {
      VERTEX_ptr = VERTEX_list.get_and_step();

      result = api_find_named_attribute( VERTEX_ptr, att_name, ret_att );

      if( !result.ok() || !ret_att )
        continue;

      ent_id = ((ATTRIB_GEN_INTEGER *) ret_att)->value();

      if( ent_id == input_ent_id )
      {
        if( found == 1 )
        {
          // More than one entity with same att
          output_ENTITY_ptr = NULL;
          return CUBIT_FAILURE;
        }
        output_ENTITY_ptr = (ENTITY *)VERTEX_ptr;
        found = 1;
        break;
      }
    }
    if( found )
      return CUBIT_SUCCESS;
    else
      return CUBIT_FAILURE;
  }

  return CUBIT_FAILURE;
}

CubitStatus
AcisTweakTool::reassign_cubit_owners_from_tweak_attribs( BODY *BODY_ptr,
     const char *att_name,
     DLIList<FACE*> &FACE_list, DLIList<AcisBridge*> &ab_FACE_list,
     DLIList<EDGE*> &EDGE_list, DLIList<AcisBridge*> &ab_EDGE_list,
     DLIList<VERTEX*> &VERTEX_list, DLIList<AcisBridge*> &ab_VERTEX_list )
{

  int i;
  ENTITY *corr_ENT_ptr;
  AcisBridge *ab_ptr;

  // FACEs
  FACE *FACE_ptr;
  FACE_list.reset();
  ab_FACE_list.reset();
  for( i=0; i<FACE_list.size(); i++ )
  {
    FACE_ptr = FACE_list.get_and_step();
    ab_ptr = ab_FACE_list.get_and_step();
    if( find_corresponding_entity_from_tweak_attrib(
      BODY_ptr, att_name, FACE_TYPE, i, corr_ENT_ptr ) == CUBIT_SUCCESS )
    {
      if( ab_ptr )
        ATTRIB_CUBIT_OWNER::set_cubit_owner( (FACE*)corr_ENT_ptr, ab_ptr );
    }
  }

  // EDGEs
  EDGE *EDGE_ptr;
  EDGE_list.reset();
  ab_EDGE_list.reset();
  for( i=0; i<EDGE_list.size(); i++ )
  {
    EDGE_ptr = EDGE_list.get_and_step();
    ab_ptr = ab_EDGE_list.get_and_step();
    if( find_corresponding_entity_from_tweak_attrib(
      BODY_ptr, att_name, EDGE_TYPE, i, corr_ENT_ptr ) == CUBIT_SUCCESS )
    {
      if( ab_ptr )
        ATTRIB_CUBIT_OWNER::set_cubit_owner( (EDGE*)corr_ENT_ptr, ab_ptr );
    }
  }

  // VERTICEs
  VERTEX *VERTEX_ptr;
  VERTEX_list.reset();
  ab_VERTEX_list.reset();
  for( i=0; i<VERTEX_list.size(); i++ )
  {
    VERTEX_ptr = VERTEX_list.get_and_step();
    ab_ptr = ab_VERTEX_list.get_and_step();
    if( find_corresponding_entity_from_tweak_attrib(
      BODY_ptr, att_name, VERTEX_TYPE, i, corr_ENT_ptr ) == CUBIT_SUCCESS )
    {
      if( ab_ptr )
        ATTRIB_CUBIT_OWNER::set_cubit_owner( (VERTEX*)corr_ENT_ptr, ab_ptr );
    }
  }

  return CUBIT_SUCCESS;
}

CubitStatus
AcisTweakTool::remove_named_attribs( BODY *BODY_ptr, const char *name )
{
  DLIList<FACE*> FACE_list;
  AQE->get_FACEs( BODY_ptr, FACE_list );

  int i;
  FACE *FACE_ptr;
  for( i=FACE_list.size(); i--; )
  {
    FACE_ptr = FACE_list.get_and_step();
    api_remove_generic_named_attribute(FACE_ptr, name);
  }

  // EDGEs
  DLIList<EDGE*> EDGE_list;
  AQE->get_EDGEs( BODY_ptr, EDGE_list );

  EDGE *EDGE_ptr;
  for( i=EDGE_list.size(); i--; )
  {
    EDGE_ptr = EDGE_list.get_and_step();
    api_remove_generic_named_attribute( EDGE_ptr, name );
  }

  // VERTICEs
  DLIList<VERTEX*> VERTEX_list;
  AQE->get_VERTICEs( BODY_ptr, VERTEX_list );

  VERTEX *VERTEX_ptr;
  for( i=VERTEX_list.size(); i--; )
  {
    VERTEX_ptr = VERTEX_list.get_and_step();
    api_remove_generic_named_attribute( VERTEX_ptr, name );
  }

  return CUBIT_SUCCESS;
}

CubitStatus
AcisTweakTool::blend_edges( DLIList<EDGE*> EDGE_list, double radius  )
{
  int i;
  EDGE *EDGE_ptr;
  ENTITY_LIST ENTITY_list;

  EDGE_list.reset();
  for( i=EDGE_list.size(); i--; )
  {
    EDGE_ptr = EDGE_list.get_and_step();
    if( EDGE_ptr )
      ENTITY_list.add( EDGE_ptr );
  }

  if( ENTITY_list.count() == 0 )
    return CUBIT_FAILURE;

  BODY *BODY_ptr = AcisQueryEngine::instance()->get_BODY_of_ENTITY( ENTITY_list[0] );

  // Setup tweak attributes so we can preserve Cubit owners
  DLIList<FACE*> pre_FACE_list;
  DLIList<EDGE*> pre_EDGE_list;
  DLIList<VERTEX*> pre_VERTEX_list;
  DLIList<AcisBridge*> ab_FACE_list, ab_EDGE_list, ab_VERTEX_list;
  assign_tweak_attribs( BODY_ptr, "tweak", pre_FACE_list, ab_FACE_list,
    pre_EDGE_list, ab_EDGE_list, pre_VERTEX_list, ab_VERTEX_list );

  // Now, blend the edges on this body
  outcome result = api_blend_edges( ENTITY_list, radius );

  if( !result.ok() )
  {
    AcisQueryEngine::instance()->ACIS_API_error(result);
    return CUBIT_FAILURE;
  }

  //printout warnings too
  err_mess_type *warnings;
  int nwarn = get_warnings( warnings );
  if( nwarn > 0 )
  {
    PRINT_WARNING("API Warning(s):\n");
    for( i=0; i<nwarn; i++ )
      PRINT_INFO("%s\n", find_err_mess( warnings[i]));
  }

  // Replace Cubit owners
  reassign_cubit_owners_from_tweak_attribs( BODY_ptr, "tweak", pre_FACE_list,
    ab_FACE_list, pre_EDGE_list, ab_EDGE_list, pre_VERTEX_list,
    ab_VERTEX_list );

  // Remove tweak attributes
  remove_named_attribs( BODY_ptr, "tweak" );

  return CUBIT_SUCCESS;
}

CubitStatus
AcisTweakTool::chamfer_edges( DLIList<EDGE*> EDGE_list, double r1, double r2  )
{
  if( r2 <= 0.0 )
    r2 = r1;

  int i;
  EDGE *EDGE_ptr;
  ENTITY_LIST ENTITY_list;

  EDGE_list.reset();
  for( i=EDGE_list.size(); i--; )
  {
    EDGE_ptr = EDGE_list.get_and_step();
    if( EDGE_ptr )
      ENTITY_list.add( EDGE_ptr );
  }

  if( ENTITY_list.count() == 0 )
    return CUBIT_FAILURE;

  BODY *BODY_ptr = AcisQueryEngine::instance()->get_BODY_of_ENTITY( ENTITY_list[0] );

  // Setup tweak attributes so we can preserve Cubit owners
  DLIList<FACE*> pre_FACE_list;
  DLIList<EDGE*> pre_EDGE_list;
  DLIList<VERTEX*> pre_VERTEX_list;
  DLIList<AcisBridge*> ab_FACE_list, ab_EDGE_list, ab_VERTEX_list;
  assign_tweak_attribs( BODY_ptr, "tweak", pre_FACE_list, ab_FACE_list,
    pre_EDGE_list, ab_EDGE_list, pre_VERTEX_list, ab_VERTEX_list );

  // Now, blend the edges on this body
  outcome result = api_chamfer_edges( ENTITY_list, r1, r2 );

  if( !result.ok() )
  {
    AcisQueryEngine::instance()->ACIS_API_error(result);
    return CUBIT_FAILURE;
  }
  else
    PRINT_INFO("Result was okay\n");

  // Replace Cubit owners
  reassign_cubit_owners_from_tweak_attribs( BODY_ptr, "tweak", pre_FACE_list,
    ab_FACE_list, pre_EDGE_list, ab_EDGE_list, pre_VERTEX_list,
    ab_VERTEX_list );

  // Remove tweak attributes
  remove_named_attribs( BODY_ptr, "tweak" );

  return CUBIT_SUCCESS;
}

CubitStatus
AcisTweakTool::chamfer_vertices( DLIList<VERTEX*> VERTEX_list, double radius  )
{
  int i;
  VERTEX *VERTEX_ptr;
  ENTITY_LIST ENTITY_list;

  VERTEX_list.reset();
  for( i=VERTEX_list.size(); i--; )
  {
    VERTEX_ptr = VERTEX_list.get_and_step();
    if( VERTEX_ptr )
      ENTITY_list.add( VERTEX_ptr );
  }

  if( ENTITY_list.count() == 0 )
    return CUBIT_FAILURE;

  BODY *BODY_ptr = AcisQueryEngine::instance()->get_BODY_of_ENTITY( ENTITY_list[0] );

  // Setup tweak attributes so we can preserve Cubit owners
  DLIList<FACE*> pre_FACE_list;
  DLIList<EDGE*> pre_EDGE_list;
  DLIList<VERTEX*> pre_VERTEX_list;
  DLIList<AcisBridge*> ab_FACE_list, ab_EDGE_list, ab_VERTEX_list;
  assign_tweak_attribs( BODY_ptr, "tweak", pre_FACE_list, ab_FACE_list,
    pre_EDGE_list, ab_EDGE_list, pre_VERTEX_list, ab_VERTEX_list );

  // Now, blend the edges on this body
  outcome result = api_chamfer_vertices( ENTITY_list, radius );

  if( !result.ok() )
  {
    AcisQueryEngine::instance()->ACIS_API_error(result);
    return CUBIT_FAILURE;
  }

  // Replace Cubit owners
  reassign_cubit_owners_from_tweak_attribs( BODY_ptr, "tweak", pre_FACE_list,
    ab_FACE_list, pre_EDGE_list, ab_EDGE_list, pre_VERTEX_list,
    ab_VERTEX_list );

  // Remove tweak attributes
  remove_named_attribs( BODY_ptr, "tweak" );

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : get_EDGES_by_BODY
// Member Type: PRIVATE
// Description: Get separate lists of EDGEs from the input Curves by common
//              BODY
// Author     : Steve Storm
// Date       : 12/18/05
//=============================================================================
CubitStatus
AcisTweakTool::get_EDGES_by_BODY( DLIList<Curve*> &input_curve_list,
                                  DLIList<BODY*> &BODY_list,
                                  DLIList<DLIList<EDGE*>*> &BODY_EDGE_lists )
{
  int i;
  Curve *curve_ptr;
  EDGE *EDGE_ptr;
  BODY *BODY_ptr;
  BODY *curr_BODY_ptr;

  // Copy input curve list since we don't want to modify it
  DLIList<Curve*> copied_input_curve_list = input_curve_list;

  // Pull curves out of the input curve list that are from a common BODY
  DLIList<EDGE*> *curr_BODY_EDGE_list_ptr = NULL;
  while( copied_input_curve_list.size() )
  {
    curr_BODY_ptr = NULL;
    copied_input_curve_list.reset();
    for( i=copied_input_curve_list.size(); i--; )
    {
      curve_ptr = copied_input_curve_list.get();
      EDGE_ptr = AQE->get_EDGE( curve_ptr );
      if( EDGE_ptr == NULL )
      {
        PRINT_ERROR( "Non-ACIS curve encountered.\n" );
        BODY_list.clean_out();
        while( BODY_EDGE_lists.size() ) delete BODY_EDGE_lists.pop();
        return CUBIT_FAILURE;
      }

      BODY_ptr = AQE->get_BODY_of_ENTITY( EDGE_ptr );
      if( curr_BODY_ptr == NULL )
      {
        BODY_list.append( BODY_ptr );
        curr_BODY_ptr = BODY_ptr;
        curr_BODY_EDGE_list_ptr = new DLIList<EDGE*>;
        BODY_EDGE_lists.append( curr_BODY_EDGE_list_ptr );
      }

      if( curr_BODY_ptr == BODY_ptr )
      {
        curr_BODY_EDGE_list_ptr->append( EDGE_ptr );
        copied_input_curve_list.change_to( NULL );
      }

      copied_input_curve_list.step();
    }

    copied_input_curve_list.remove_all_with_value( NULL );
    curr_BODY_ptr = NULL;
  }

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : get_EDGES_by_BODY
// Member Type: PRIVATE
// Description: Get separate lists of EDGEs from the input Curves by common
//              BODY
// Author     : Steve Storm
// Date       : 12/18/05
//=============================================================================
CubitStatus
AcisTweakTool::get_EDGES_by_BODY( DLIList<EDGE*> &input_EDGE_list,
                                  DLIList<BODY*> &BODY_list,
                                  DLIList<DLIList<EDGE*>*> &BODY_EDGE_lists )
{
  int i;
  EDGE *EDGE_ptr;
  BODY *BODY_ptr;
  BODY *curr_BODY_ptr;

  // Copy input EDGE list since we don't want to modify it
  DLIList<EDGE*> copied_input_EDGE_list = input_EDGE_list;

  // Pull EDGEs out of the input EDGE list that are from a common BODY
  DLIList<EDGE*> *curr_BODY_EDGE_list_ptr = NULL;
  while( copied_input_EDGE_list.size() )
  {
    curr_BODY_ptr = NULL;
    copied_input_EDGE_list.reset();
    for( i=copied_input_EDGE_list.size(); i--; )
    {
      EDGE_ptr = copied_input_EDGE_list.get();

      BODY_ptr = AQE->get_BODY_of_ENTITY( EDGE_ptr );
      if( curr_BODY_ptr == NULL )
      {
        BODY_list.append( BODY_ptr );
        curr_BODY_ptr = BODY_ptr;
        curr_BODY_EDGE_list_ptr = new DLIList<EDGE*>;
        BODY_EDGE_lists.append( curr_BODY_EDGE_list_ptr );
      }

      if( curr_BODY_ptr == BODY_ptr )
      {
        curr_BODY_EDGE_list_ptr->append( EDGE_ptr );
        copied_input_EDGE_list.change_to( NULL );
      }

      copied_input_EDGE_list.step();
    }

    copied_input_EDGE_list.remove_all_with_value( NULL );
    curr_BODY_ptr = NULL;
  }

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : all_complete_holes
// Member Type: PRIVATE
// Description: Determine if given EDGEs form complete internal holes.  LOOPs
//              must be on a sheet BODY.  Return those LOOPs.  All given EDGEs
//              must be in the same BODY.
// Author     : Steve Storm
// Date       : 12/18/05
//=============================================================================
CubitStatus
AcisTweakTool::all_complete_internal_loops( DLIList<EDGE*> &BODY_EDGE_list,
                                            DLIList<LOOP*> &LOOP_list )
{
  int i;
  EDGE *EDGE_ptr;
  LOOP *LOOP_ptr;

  DLIList<EDGE*> copied_EDGE_list = BODY_EDGE_list;

  while( copied_EDGE_list.size() )
  {
    copied_EDGE_list.reset();
    EDGE_ptr = copied_EDGE_list.get();

    DLIList<LOOP*> EDGE_LOOP_list;
    AQE->get_LOOPs( EDGE_ptr, EDGE_LOOP_list );
    if( EDGE_LOOP_list.size() == 0 )
    {
      LOOP_list.clean_out();
      return CUBIT_FAILURE;
    }

    // Should always be just one LOOP per EDGE
    LOOP_ptr = EDGE_LOOP_list.get();

    loop_type type;
    outcome result = api_loop_type( LOOP_ptr, type );
    if( !result.ok() )
    {
      AQE->ACIS_API_error(result);
      LOOP_list.clean_out();
      return CUBIT_FAILURE;
    }
    if( type != loop_hole )
    {
      LOOP_list.clean_out();
      return CUBIT_FAILURE;
    }

    DLIList<EDGE*> LOOP_EDGE_list;
    AQE->get_EDGEs( LOOP_ptr, LOOP_EDGE_list );

    LOOP_EDGE_list.reset();
    for( i=LOOP_EDGE_list.size(); i--; )
    {
      EDGE_ptr = LOOP_EDGE_list.get_and_step();

      if( copied_EDGE_list.move_to( EDGE_ptr ) == CUBIT_TRUE )
      {
        // Remove from the list
        copied_EDGE_list.change_to( NULL );
      }
      else
      {
        // All EDGEs must be in the input EDGE list
        LOOP_list.clean_out();
        return CUBIT_FAILURE;
      }
    }

    LOOP_list.append( LOOP_ptr );

    copied_EDGE_list.remove_all_with_value( NULL );
  }

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : remove_LOOPs
// Member Type: PRIVATE
// Description: Remove the given internal LOOPs (from a sheet BODY).
// Author     : Steve Storm
// Date       : 12/18/05
//=============================================================================
CubitStatus
AcisTweakTool::remove_LOOPs( DLIList<LOOP*> &LOOP_list )
{
  if( LOOP_list.size() == 0 )
    return CUBIT_SUCCESS;

  int i;
  FACE *FACE_ptr;
  LOOP *LOOP_ptr;

  LOOP_list.reset();
  for( i=LOOP_list.size(); i--; )
  {
    LOOP_ptr = LOOP_list.get_and_step();

    FACE_ptr = LOOP_ptr->face();
    if( FACE_ptr->sides() != DOUBLE_SIDED )
    {
      PRINT_ERROR("Cannot remove loops that are not on sheet bodies\n" );
      return CUBIT_FAILURE;
    }

    LOOP *tmp_LOOP_ptr = FACE_ptr->loop();

    if( tmp_LOOP_ptr == LOOP_ptr )
    {
      // First loop in list is an internal loop we are trying to remove
      tmp_LOOP_ptr = tmp_LOOP_ptr->next();
      FACE_ptr->set_loop( tmp_LOOP_ptr );
    }
    else
    {
      LOOP *prev_LOOP_ptr = NULL;
      LOOP *next_LOOP_ptr = NULL;
      while ( tmp_LOOP_ptr != LOOP_ptr )
      {
        prev_LOOP_ptr = tmp_LOOP_ptr;
        tmp_LOOP_ptr = tmp_LOOP_ptr->next();
      }
      next_LOOP_ptr = LOOP_ptr->next();

      // Now set the next LOOP pointer of the previous to the next
      prev_LOOP_ptr->set_next( next_LOOP_ptr );
    }

    // Now delete this LOOP and all entities below it
    outcome result = api_delent( LOOP_ptr );
    if( !result.ok() )
    {
      AcisQueryEngine::instance()->ACIS_API_error(result);
      PRINT_ERROR( "Unable to remove loop - aborting\n" );
      return CUBIT_FAILURE;
    }
  }

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : remove_holes
// Member Type: PRIVATE
// Description: Remove all holes from the given sheet BODY.
// Author     : Steve Storm
// Date       : 2/12/06
//=============================================================================
CubitStatus
AcisTweakTool::remove_holes( BODY *sheet_BODY_ptr )
{
  // Holes in a sheet body will be internal LOOPs with only 1 coedge on each
  // curve.

  // Get all LOOPs of the sheet BODY
  DLIList<LOOP*> LOOP_list;
  AQE->get_LOOPs( sheet_BODY_ptr, LOOP_list );

  // Put internal LOOPs to remove into a separate list
  DLIList<LOOP*> hole_LOOP_list;
  int i;
  outcome result;
  loop_type type;
  ENTITY_LIST COEDGES;
  LOOP *LOOP_ptr;
  for( i=LOOP_list.size(); i--; )
  {
    LOOP_ptr = LOOP_list.get_and_step();

    result = api_loop_type( LOOP_ptr, type );
    if( !result.ok() )
    {
      AQE->ACIS_API_error(result);
      return CUBIT_FAILURE;
    }
    if( type == loop_hole )
    {
      // Check coedges on edges in loop - if more than one it is shared and we
      // don't want to remove it.
      DLIList<EDGE*> EDGE_list;
      AQE->get_EDGEs( LOOP_ptr, EDGE_list );
      if( EDGE_list.size() )
      {
        EDGE *EDGE_ptr = EDGE_list.get();

        result = api_get_coedges( EDGE_ptr, COEDGES);
        if( !result.ok() )
        {
          AQE->ACIS_API_error(result);
          return CUBIT_FAILURE;
        }

        ENTITY *COEDGE_ptr;
        COEDGES.init();

        int good = 1;

        while( (COEDGE_ptr = COEDGES.next() ) != NULL )
        {
          if( COEDGES.iteration_count() == 1 )
            continue;
           good--;
        }

        COEDGES.clear();

        if( good > 0 )
          hole_LOOP_list.append( LOOP_ptr );
      }
    }
  }

  // Remove the internal LOOPs
  if( remove_LOOPs( hole_LOOP_list ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : extrema_pln_BODY
// Member Type: PRIVATE
// Description: Finds the extrema distance from a plane to an ACIS BODY
//              (perpendicular distance from plane to farthest extent of BODY),
//              on one side of the plane (by default on the front of the plane
//              - direction of plane normal).  If the entire BODY lies on the
//              other side of the plane, the extrema distance is 0.0.
// Author     : Steve Storm
// Date       : 1/15/06
//=============================================================================
CubitStatus AcisTweakTool::extrema_pln_BODY( CubitPlane &plane, BODY *BODY_ptr,
                                             double &extrema_dist,
                                             int back_side )
{
  // Method: find EXTREMA point of BODY in direction of plane norm, check
  // distance to plane.
  ENTITY_LIST ENTITY_list;
  ENTITY_list.add( BODY_ptr );

  CubitVector dir =  back_side ? -(plane.normal()) : plane.normal();

  SPAvector acis_dirs[3];
  acis_dirs[0].set_x( dir.x() );
  acis_dirs[0].set_y( dir.y() );
  acis_dirs[0].set_z( dir.z() );

  SPAposition acis_pos;
	param_info info;
	outcome result = api_entity_extrema( ENTITY_list, 1, acis_dirs, acis_pos,
                                       info );
  if( !result.ok() )
	{
    AcisQueryEngine::instance()->ACIS_API_error( result );
    return CUBIT_FAILURE;
	}

  CubitVector extrema_pnt( acis_pos.x(), acis_pos.y(), acis_pos.z() );

  // Get intersection point on plane
  CubitVector int_pnt;
  int_pnt = plane.intersect( extrema_pnt, dir );

  // Get direction from int_pnt to extrema_pnt
  CubitVector dir2( int_pnt, extrema_pnt  );

  // Note: angle between vectors acute (<90) if dot product > 0.  In
  // this case the vectors would be going the same direction and the extrema
  // point is on the correct side of the plane
  if( dir % dir2 > 0.0 )
    extrema_dist = int_pnt.distance_between( extrema_pnt );
  else
    extrema_dist = 0.0;

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : tweak_FACEs_to_target
// Member Type: PRIVATE
// Description: Tweak the given FACEs (which must all be part of the same BODY)
//              to the given target FACE.
// Author     : Steve Storm
// Date       : 1/15/06
//=============================================================================
CubitStatus
AcisTweakTool::tweak_FACEs_to_target( DLIList<FACE*> &source_FACE_list,
                                      FACE *target_FACE_ptr,
                                      CubitBoolean reverse_flg,
                                      CubitBoolean skip_self_int_check )
{
  int i;
  outcome result;
  int num_faces = source_FACE_list.size();

  // Get normal of tool (target) surface
  CubitVector location;
  CubitVector target_normal_vec = surface_normal( target_FACE_ptr, location );
  if( reverse_flg )
    target_normal_vec = -target_normal_vec;

  // Save the center point of the bounding box of the target face
  SPAposition target_surf_bb_center( location.x(), location.y(), location.z() );

  // Create arrays of input FACEs for tweak API.
  // Note: Apparently the API requires a copy of the tool geometry for each
  //       tweak surface
  FACE** source_FACES = new FACE*[num_faces];
  DLIList<FACE*> target_FACE_list(num_faces);
  SURFACE** target_SURFACES = new SURFACE*[num_faces];
  int* target_sense = new int[num_faces];
  FACE *FACE_ptr;
  source_FACE_list.reset();
  for( i=0; i<num_faces; i++ )
  {
    // Fill source_FACES array
    FACE_ptr = source_FACE_list.get_and_step();
    source_FACES[i] = FACE_ptr;

    // Fill target_SURFACES array
    ENTITY* copied_entity_ptr = 0;
    api_copy_entity_contents(target_FACE_ptr, copied_entity_ptr);
    if( copied_entity_ptr == 0 )
    {
      while( target_FACE_list.size() ) api_delent( target_FACE_list.pop() );
      delete []source_FACES;
      delete []target_SURFACES;
      delete []target_sense;
      PRINT_ERROR( "Invalid ACIS target surface -- aborting\n" );
      return CUBIT_FAILURE;
    }
    target_FACE_list.append( (FACE*)copied_entity_ptr );
    target_SURFACES[i] = ((FACE*)copied_entity_ptr)->geometry();

    SPAposition dummy;
    SPAunit_vector tmp_vec;
    FACE_ptr->geometry()->equation().point_perp( target_surf_bb_center, dummy, tmp_vec );
    CubitVector tweak_normal_vec( tmp_vec.x(), tmp_vec.y(), tmp_vec.z() );

    // Fill target_sense array
    target_sense[i] = ((tweak_normal_vec % target_normal_vec) < 0.0);
  }

  SPAposition box_h(0,0,0);
  SPAposition box_l(0,0,0);

  // This operation sometimes destroys owner attributes
  BODY *copied_BODY_ptr = AQE->get_BODY_of_ENTITY( source_FACE_list.get() );
  AcisBridge *ab_body_ptr;
  AcisBridge *ab_lump_ptr;
  AcisBridge *ab_shell_ptr;
  get_owner_attribs( copied_BODY_ptr, ab_body_ptr, ab_lump_ptr, ab_shell_ptr );

  // Setup tweak attributes so we can preserve Cubit owners
  DLIList<FACE*> pre_FACE_list;
  DLIList<EDGE*> pre_EDGE_list;
  DLIList<VERTEX*> pre_VERTEX_list;
  DLIList<AcisBridge*> ab_FACE_list, ab_EDGE_list, ab_VERTEX_list;
  assign_tweak_attribs( copied_BODY_ptr, "tweak", pre_FACE_list, ab_FACE_list,
    pre_EDGE_list, ab_EDGE_list, pre_VERTEX_list, ab_VERTEX_list );

#if CUBIT_ACIS_VERSION < 1600
  result = api_tweak_faces( num_faces, (FACE**)source_FACES,
                            target_SURFACES, target_sense, box_h, box_l );
#else
  if( skip_self_int_check )
  {
    lop_options lop_opt = lop_options();
    lop_opt.set_repair_self_int(FALSE);
    result = api_tweak_faces( num_faces, (FACE**)source_FACES,
                              target_SURFACES, target_sense, box_h, box_l, &lop_opt );
  }
  else
  {
    result = api_tweak_faces( num_faces, (FACE**)source_FACES,
                              target_SURFACES, target_sense, box_h, box_l );
  }
#endif

  // Reset the owner attributes
  reset_owner_attribs( copied_BODY_ptr, ab_body_ptr, ab_lump_ptr, ab_shell_ptr );

  // Replace Cubit owners
  reassign_cubit_owners_from_tweak_attribs( copied_BODY_ptr, "tweak",
    pre_FACE_list, ab_FACE_list, pre_EDGE_list, ab_EDGE_list, pre_VERTEX_list,
    ab_VERTEX_list );

  remove_named_attribs( copied_BODY_ptr, "tweak" );

  // Freeup memory
  while( target_FACE_list.size() ) api_delent( target_FACE_list.pop() );
  delete []source_FACES;
  delete []target_SURFACES;
  delete []target_sense;

  if( !result.ok() )
  {
    AcisQueryEngine::instance()->ACIS_API_error(result);
    return CUBIT_FAILURE;
  }

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : surface_normal
// Member Type: PRIVATE
// Description: Get normal of underlying surface of given FACE.  The normal is
//              found at the approximate centroid of the FACE.  Note the normal
//              is NOT adjusted for the FACE sense.
// Author     : Steve Storm
// Date       : 1/15/06
//=============================================================================
CubitVector
AcisTweakTool::surface_normal( FACE *FACE_ptr, CubitVector &location,
                               CubitBoolean calc_loc )
{
  SPAunit_vector acis_normal;

  if( calc_loc )
  {
    // Get bounding box of FACE
    CubitBox box = AQE->bounding_box(AQE->bounding_box( FACE_ptr ));

    // Get center point of bounding box
    CubitVector center = box.center();

    // Move to the FACE and get the normal at that location
    SPAposition acis_point ( center.x(), center.y(), center.z() );
    SPAposition acis_point_on_surf;
    (FACE_ptr->geometry()->equation()).point_perp( acis_point,
      acis_point_on_surf, acis_normal );

    location.set( acis_point_on_surf.x(), acis_point_on_surf.y(),
      acis_point_on_surf.z() );
  }
  else
  {
    SPAposition acis_point ( location.x(), location.y(), location.z() );
    SPAposition acis_point_on_surf;
    (FACE_ptr->geometry()->equation()).point_perp( acis_point,
      acis_point_on_surf, acis_normal );
  }

  CubitVector return_vector( acis_normal.x(), acis_normal.y(),
    acis_normal.z() );

  return return_vector;
}

//=============================================================================
// Function   : FACE_normal
// Member Type: PRIVATE
// Description: Get normal of given FACE.  The normal is found at the
//              approximate centroid of the FACE.  Note the normal is adjusted
//              for the FACE sense.
// Author     : Steve Storm
// Date       : 1/28/06
//=============================================================================
CubitVector
AcisTweakTool::FACE_normal( FACE *FACE_ptr, CubitVector &location,
                            CubitBoolean calc_loc )
{
  CubitVector FACE_normal = surface_normal( FACE_ptr, location, calc_loc );

  // Adjust for the FACE sense
  if( FACE_ptr->sense() == REVERSED )
    FACE_normal = -FACE_normal;

  return FACE_normal;
}

//=============================================================================
// Function   : weighted_average_FACE_normal
// Member Type: PRIVATE
// Description: Get weighted average normal of given FACEs.  The normals are
//              weighted by FACE area utilizing the graphics facets.
// Author     : Steve Storm
// Date       : 2/21/06
//=============================================================================
CubitVector
AcisTweakTool::weighted_average_FACE_normal( DLIList<FACE*> &FACE_list )
{
  int i;
  FACE *FACE_ptr;
  CubitVector total_vec;

  for( i=FACE_list.size(); i--; )
  {
    FACE_ptr = FACE_list.get_and_step();
    CubitVector FACE_norm;
    double FACE_weight;
    weighted_average_FACE_normal( FACE_ptr, FACE_norm, FACE_weight );

    total_vec = total_vec += FACE_weight * FACE_norm;
  }

  //PRINT_INFO( "Weighted normal vec = %f, %f, %f\n", total_vec.x(), total_vec.y(),
  //  total_vec.z() );

  total_vec.normalize();
  //PRINT_INFO( "Weighted normal normalized = %f, %f, %f\n", total_vec.x(), total_vec.y(),
  //  total_vec.z() );

  return total_vec;
}

//=============================================================================
// Function   : weighted_average_FACE_normal
// Member Type: PRIVATE
// Description: Get normal (unit vector) and weight (area factor) of a given
//              FACE.  The graphics facets are utilized in the calculation.
// Author     : Steve Storm
// Date       : 2/21/06
//=============================================================================
CubitStatus
AcisTweakTool::weighted_average_FACE_normal( FACE *FACE_ptr,
                                             CubitVector &normal,
                                             double &weight )
{
  int i;
  int num_tris, num_pnts, num_facets;
  GMem g_mem;
  unsigned short norm_tol = 10;
  double dist_tol = -1.0;
  AQE->get_graphics( FACE_ptr, num_tris, num_pnts, num_facets,
    &g_mem, norm_tol, dist_tol );

  if(num_tris < 1)
  {
    // Decrease tolerance and try again (we can get this for small features)
    norm_tol /= 2;
    AQE->get_graphics( FACE_ptr, num_tris, num_pnts, num_facets, &g_mem,
      norm_tol, dist_tol);
  }

  if(num_tris < 1)
  {
    // Lets give up
    PRINT_ERROR( "Unable to find average normal of a surface\n" );
    return CUBIT_FAILURE;
  }

  // Initialize
  weight = 0.0;
  normal.set( 0.0, 0.0, 0.0 );

  // Loop through the triangles
  double tri_weight, A, B, C;
  GPoint p[3];
  GPoint* plist = g_mem.point_list();
  int* facet_list = g_mem.facet_list();
  int c = 0;
  for( i=0; i<num_tris; i++ )
  {
    p[0] = plist[facet_list[++c]];
    p[2] = plist[facet_list[++c]];
    p[1] = plist[facet_list[++c]];
    c++;

    // Get centroid
    CubitVector p1( p[0].x, p[0].y, p[0].z );
    CubitVector p2( p[2].x, p[2].y, p[2].z );
    CubitVector p3( p[1].x, p[1].y, p[1].z );

    CubitVector center = (p1 + p2 + p3)/3.0;

    // Get normal direction at the point closest to the location on the
    // surface.  This method usually returns the normal, but if the nearest
    // point is a singularity (like the apex of a cone), this method still
    // returns an outward direction. The base class definition returns
    // point_normal, which is used by default on simple surfaces.
    SPAposition acis_center( center.x(), center.y(), center.z() );
    SPAunit_vector acis_normal = (FACE_ptr->geometry()->equation()).
      point_outdir( acis_center );

    CubitVector norm( acis_normal.x(), acis_normal.y(), acis_normal.z() );

    // Get triangle area
    A = p1.y() * p2.z() + p1.z() * p3.y() + p2.y() * p3.z() -
      p2.z() * p3.y() - p1.y() * p3.z() - p1.z() * p2.y();

    B = p1.z() * p2.x() + p1.x() * p3.z() + p2.z() * p3.x() -
      p2.x() * p3.z() - p1.z() * p3.x() - p1.x() * p2.z();

    C = p1.x() * p2.y() + p1.y() * p3.x() + p2.x() * p3.y() -
      p2.y() * p3.x() - p1.x() * p3.y() - p1.y() * p2.x();

    //Note: triangle area = 0.5*(sqrt(A*A+B*B+C*C));

    tri_weight = 0.5*(A*A+B*B+C*C);

    normal += tri_weight * norm;

    weight += tri_weight;

  }

  //double FACE_area;
  //double accur_achieved;
  //outcome result = api_ent_area( FACE_ptr, .001, FACE_area, accur_achieved );
  //PRINT_INFO( "\nArea by ACIS = %f\n", FACE_area );
  //PRINT_INFO( "Area by facets = %f\n", weight );

  normal.normalize();

  if( FACE_ptr->sense() == REVERSED )
    normal = -normal;

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : create_offset_planar_body
// Member Type: PRIVATE
// Description: Create a planar ACIS BODY from the given CubitPlane and offset.
//              The offset is in the direction of the plane's normal.  The
//              plane is arbitrarily 10x10x10 in size near the origin.
// Author     : Steve Storm
// Date       : 3/19/06
//=============================================================================
CubitStatus
AcisTweakTool::create_offset_planar_body( CubitPlane &plane, double offset,
                                          BODY *&planar_BODY_ptr )
{
  // Get corners of the surface
  CubitVector normal = plane.normal();
  CubitVector x, y;
  normal.orthogonal_vectors( x, y );
  // There is a bug in the point_on_plane (it needs to check against a small
  // value instead of an exact 0) - for now, do our own calculation.  This
  // allows us to check for an invalid plane anyway...
  //CubitVector p1 = target_plane.point_on_plane();
  CubitVector p1;
  double coeff = plane.coefficient();
  if( fabs( normal.x() ) > 1e-6 )
    p1.set( -coeff/normal.x(), 0.0, 0.0 );
  else if( fabs( normal.y() ) > 1e-6 )
    p1.set( 0.0, -coeff/normal.y(), 0.0 );
  else if( fabs( normal.z() ) > 1e-6 )
    p1.set( 0.0, 0.0, -coeff/normal.z() );
  else
  {
    PRINT_ERROR( "Error constructing offset plane past the targets\n" );
    return CUBIT_FAILURE;
  }

  p1.next_point( normal, offset, p1 );
  CubitVector p2, p3, p4;
  p1.next_point( x, 5.0, p1 );
  p1.next_point( y, 5.0, p1 );
  p1.next_point( -x, 10.0, p2 );
  p2.next_point( -y, 10.0, p3 );
  p3.next_point( x, 10.0, p4 );

  planar_BODY_ptr = AME->make_planar_quad_BODY( p1, p2, p3, p4 );
  if( !planar_BODY_ptr )
  {
    PRINT_ERROR( "Error constructing offset plane past the targets\n" );
    return CUBIT_FAILURE;
  }

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : create_extended_sheet
// Member Type: PRIVATE
// Description: Extend out a set of FACEs and place into a new sheet BODY.
//              Holes are removed and the BODY is healed for consistent FACE
//              normals.  The BODY is also regularized.
// Author     : Steve Storm
// Date       : 1/15/06
//=============================================================================
CubitStatus
AcisTweakTool::create_extended_sheet( DLIList<FACE*> &FACE_list,
                                      BODY *&ext_BODY_ptr )
{
  // Method:
  // Copy FACEs into separate sheet bodies
  // If not entirely surrounded by other FACEs, extend past super bounding box
  // Clip by a large bounding box so they are same size
  // Unite them all together
  // Remove FACEs that don't overlap original FACEs
  // Make sure FACE normals are consistent (we do this by healing)
  // Remove any holes
  // Remove unwanted topology (regularize)

  // It is useful to heal and regularize the surfaces first.  The method we use
  // has problems with spline surfaces - by healing we hopefully simplify the
  // splines to analytics.
  BODY *sheet_body;
  if( copy_FACEs_into_sheet( FACE_list, sheet_body ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  // Get the new FACEs
  DLIList<FACE*> new_FACE_list;
  AQE->get_FACEs( sheet_body, new_FACE_list );

  int i;
  FACE *FACE_ptr;
  FACE *new_FACE_ptr;
  BODY *BODY_ptr;
  DLIList<BODY*> ext_BODY_list;
  new_FACE_list.reset();
  for( i=new_FACE_list.size(); i--; )
  {
    FACE_ptr = new_FACE_list.get_and_step();

    // Do not extend if entirely enclosed by the other surfaces
    CubitBoolean extend = FACE_surrounded( FACE_ptr, new_FACE_list )
      ? CUBIT_FALSE : CUBIT_TRUE;

    new_FACE_ptr = AME->make_FACE( FACE_ptr, extend );
    if( new_FACE_ptr == NULL )
    {
      PRINT_ERROR( "Problem extending surface - aborting\n" );
      api_delent( sheet_body );
      while( ext_BODY_list.size() )
        api_delent( ext_BODY_list.pop() );
      return CUBIT_FAILURE;
    }
    // Get the BODY this new FACE is in
    BODY_ptr = AQE->get_BODY_of_ENTITY( new_FACE_ptr );
    if( BODY_ptr == NULL )
    {
      PRINT_ERROR( "Problem extending surface - aborting\n" );
      api_delent( sheet_body );
      while( ext_BODY_list.size() )
        api_delent( ext_BODY_list.pop() );
      return CUBIT_FAILURE;
    }

    ext_BODY_list.append( BODY_ptr );
  }

  // We are done with the sheet_body
  api_delent( sheet_body );

  // Unite all of the extended BODIES together
  outcome result;
  DLIList<BODY*> BODY_list;
  ext_BODY_list.reset();
  BODY *master = ext_BODY_list.extract();
  for( i=ext_BODY_list.size(); i--; )
  {
    BODY_ptr = ext_BODY_list.extract();

    // Do the union of the master and the BODY_ptr.
    // If this is successful, the result is master and
    //   BODY_ptr will be deleted
    result = api_boolean( BODY_ptr, master, UNION );
    if( !result.ok() || (!master) )
    {
      AQE->ACIS_API_error(result);
      while( ext_BODY_list.size() )
        api_delent( ext_BODY_list.pop() );
      if (master != NULL) api_delent(master);
      return CUBIT_FAILURE;
    }
  }

  // Cut by a bounding box, so all entities are the same extents.  Otherwise,
  // cylinders won't be handled correctly (they may remain full cylinders
  // which is undesireable).
  CubitBox box = GeometryQueryTool::instance()->model_bounding_box();

  // Check the box for zero dimensions, which can definitely occur if we are
  // working with 2D sheet bodies.  If those exist, add a factor in that
  // direction.  Here we pick 0.25 since when we thicken bodies we do it by
  // 0.2 - thus this should be safe to overlap those thickened bodies
  // completely.
  CubitVector box_min = box.minimum();
  CubitVector box_max = box.maximum();

  double f = .25;

  if( box_max.x() - box_min.x() < 1e-6 )
  {
    box_min.set( box_min.x()-f, box_min.y(), box_min.z() );
    box_max.set( box_max.x()+f, box_min.y(), box_min.z() );
  }
  if( box_max.y() - box_min.y() < 1e-6 )
  {
    box_min.set( box_min.x(), box_min.y()-f, box_min.z() );
    box_max.set( box_max.x(), box_max.y()+f, box_max.z() );
  }
  if( box_max.z() - box_min.z() < 1e-6 )
  {
    box_min.set( box_min.x(), box_min.y(), box_min.z()-f );
    box_max.set( box_max.x(), box_max.y(), box_max.z()+f );
  }

  CubitBox temp_box( box_min, box_max );
  box |= temp_box;

  // Scale the box up overall by a factor as well
  box *= 1.25;

  // Create sheet bodies (6) at box sides
  CubitVector c[8];
  box.get_corners( c );

  //     3+----------+2
  //     /|         /|
  //    / |        / |
  //   /  |       /  |      Y
  // 7+----------+6  |      ^
  //  |   |      |   |      |
  //  |  0+------|---+1     |
  //  |  /       |  /       +---->X
  //  | /        | /       /
  //  |/         |/       /
  // 4+----------+5      Z

  SPAposition min( box.minimum().x(), box.minimum().y(), box.minimum().z() );
  SPAposition max( box.maximum().x(), box.maximum().y(), box.maximum().z() );
  SPAbox super_box( min, max );

  // Note - infinite plane is swept in the opposite direction from the normal
  // of the passed in points
  // Create infinite plane object at the FRONT of the box
  BODY_ptr = AME->create_infinite_plane_cutting_tool( c[6], c[5],
    c[4], super_box, false );
  // Subtract this object from the master (BODY_ptr is deleted)
  result = api_boolean( BODY_ptr, master, SUBTRACTION );
  // Do this for the other 5 sides of the box
  // RIGHT SIDE
  BODY_ptr = AME->create_infinite_plane_cutting_tool( c[6], c[2],
    c[1], super_box, false );
  result = api_boolean( BODY_ptr, master, SUBTRACTION );
  // BACK SIDE
  BODY_ptr = AME->create_infinite_plane_cutting_tool( c[1], c[2],
    c[3], super_box, false );
  result = api_boolean( BODY_ptr, master, SUBTRACTION );
  // LEFT SIDE
  BODY_ptr = AME->create_infinite_plane_cutting_tool( c[0], c[3],
    c[7], super_box, false );
  result = api_boolean( BODY_ptr, master, SUBTRACTION );
  // BOTTOM
  BODY_ptr = AME->create_infinite_plane_cutting_tool( c[0], c[4],
    c[5], super_box, false );
  result = api_boolean( BODY_ptr, master, SUBTRACTION );
  // TOP
  BODY_ptr = AME->create_infinite_plane_cutting_tool( c[3], c[2],
    c[6], super_box, false );
  result = api_boolean( BODY_ptr, master, SUBTRACTION );

  // Put original FACEs into a single sheet BODY
  BODY *orig_sheet;
  if( copy_FACEs_into_sheet( FACE_list, orig_sheet )
    == CUBIT_FAILURE )
  {
    api_delent( master );
    return CUBIT_FAILURE;
  }

  // Remove FACEs that don't overlap with the original target FACEs
  double overlap_area;
  DLIList<FACE*> master_FACE_list;
  AQE->get_FACEs( master, master_FACE_list );
  master_FACE_list.reset();
  for( i=master_FACE_list.size(); i--; )
  {
    FACE_ptr = master_FACE_list.get_and_step();

    // Make a copy of the FACE into a new BODY
    FACE *temp_FACE_ptr = AME->make_FACE( FACE_ptr, CUBIT_FALSE );
    if( temp_FACE_ptr == NULL )
    {
      api_delent( master );
      return CUBIT_FAILURE;
    }
    BODY *temp_BODY_ptr = AQE->get_BODY_of_ENTITY( temp_FACE_ptr );

    if( get_overlap_area( temp_BODY_ptr, orig_sheet, overlap_area )
      == CUBIT_FAILURE )
    {
      api_delent( temp_BODY_ptr );
      api_delent( master );
      return CUBIT_FAILURE;
    }

    api_delent( temp_BODY_ptr );

    if( overlap_area > 1e-5 )
      continue;
    else
      api_remove_face( FACE_ptr );
  }

  // Remove any internal holes in the body
  remove_holes( master );

  // Clean (regularize) to remove unwanted topology
  api_clean_entity( master );

  // Make sure FACE normals are consistent (api_body_to_1d could be used for
  // this, then convert back to double sided).  Healing will also work.
  // Ignore return status as this is not fatal.
  heal_BODY( master );

  // Set output BODY
  ext_BODY_ptr = master;

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : prep_cutting_sheet
// Member Type: PRIVATE
// Description: Prepare the cutting sheet for chopping off the BODY in a 
//              tweak_target_multiple operation.
// Author     : Steve Storm
// Date       : 3/18/06
//=============================================================================
CubitStatus
AcisTweakTool::prep_cutting_sheet( BODY *&cutting_sheet,
                                   BODY *tweaked_BODY_ptr,
                                   BODY *ext_target_BODY_ptr,
                                   CubitVector &source_norm,
                                   CubitVector &target_norm,
                                   CubitBoolean check_crossing )
{
  int i, j;
  outcome result;

  // Imprint, so as to score the cutting sheet.  We have to do this with a
  // copied BODY, since the extended target could cross FACEs we don't want
  // to affect.
  BODY *tweaked_BODY_ptr2 = AME->copy_BODY( tweaked_BODY_ptr, CUBIT_TRUE );
  result = api_imprint( cutting_sheet, tweaked_BODY_ptr2 );
  if( !result.ok() )
  {
    AQE->ACIS_API_error(result);
    api_delent( tweaked_BODY_ptr2 );
    return CUBIT_FAILURE;
  }
  api_delent( tweaked_BODY_ptr2 );

  if( check_crossing == CUBIT_TRUE )
  {
    // Throw away individual cutting sheet patches that don't cross the FACEs we
    // want to cut.
    DLIList<FACE*> cut_FACE_list;
    DLIList<FACE*> temp_FACE_list;
    AQE->get_FACEs( tweaked_BODY_ptr, temp_FACE_list );
    FACE *FACE_ptr;
    ATTRIB_GEN_NAME *ret_att;
    for( i=temp_FACE_list.size(); i--; )
    {
      FACE_ptr = temp_FACE_list.get_and_step();
      api_find_named_attribute( FACE_ptr, "tweak", ret_att );
      if( ret_att )
        cut_FACE_list.append( FACE_ptr );
    }
    BODY *cut_FACE_sheet;
    // Note: healing a .2 thick sheet here will collapse the side FACEs out,
    // so we skip healing.
    copy_FACEs_into_sheet( cut_FACE_list, cut_FACE_sheet, CUBIT_FALSE );

    // Intersect this sheet with each cutting sheet. NULL or incomplete wires
    // will indicate an invalid crossing.
    BODY **sheets = NULL;
    int n_body = 0;
    result = api_separate_body( cutting_sheet, n_body, sheets );
    if (!result.ok())
    {
      AQE->ACIS_API_error(result);
      api_delent( cut_FACE_sheet );
      return CUBIT_FAILURE;
    }

    CubitBoolean good_crossing = CUBIT_FALSE;
    for( i=0; i<n_body; i++ )
    {
      BODY* slice = NULL;

      // Create the intersection graph between the given bodies
      result = api_slice( sheets[i], cut_FACE_sheet,
        *(SPAunit_vector *)NULL_REF, slice);

      if (!result.ok() || slice == NULL )
        good_crossing = CUBIT_FALSE;
      else
      {
        api_clean_wire( slice );
        int num = 0;
        BODY **wires = NULL;
        result = api_separate_body( slice, num, wires );
        if (!result.ok())
        {
          AQE->ACIS_API_error(result);
          api_delent( slice );
          for( j=0; j<n_body; j++ )
            if( sheets[j] ) api_delent(sheets[j]);
          api_delent( cut_FACE_sheet );
          return CUBIT_FAILURE;
        }
        if( num == 1 )
          good_crossing = CUBIT_TRUE;
        else
          good_crossing = CUBIT_FALSE;

        for( j=0; j<num; j++ )
          api_delent( wires[j] );
      }

      if( good_crossing == CUBIT_FALSE )
      {
        // Delete this BODY
        api_delent( sheets[i] );
        sheets[i] = NULL;
      }
    }

    // Delete the cut_FACE_sheet since we are done with it
    api_delent( cut_FACE_sheet );

    // Unite the remaining sheets back into cutting_sheet
    BODY *first_BODY_ptr = NULL;
    DLIList<BODY*> unite_BODY_list;
    for( i=0; i<n_body; i++ )
    {
      if( sheets[i] )
      {
        if( !first_BODY_ptr )
          first_BODY_ptr = sheets[i];
        else
          unite_BODY_list.append( sheets[i] );
      }
    }
    if( !first_BODY_ptr )
    {
      PRINT_ERROR( "Problem intersecting extended surface(s) with target\n" );
      return CUBIT_FAILURE;
    }
    if( unite_BODIES( first_BODY_ptr, unite_BODY_list, cutting_sheet )
      == CUBIT_FAILURE )
    {
      for( i=0; i<n_body; i++ )
        if( sheets[i] ) api_delent( sheets[i] );
      return CUBIT_FAILURE;
    }
  }

  // Remove periodic surfaces near nonmanifold EDGEs. Note if any remaining
  // nonmanifold EDGEs exist, this function will return failure (a nonmanifold
  // surface won't work as a cutting sheet).  See picture in the function
  // for an example.
  if( remove_per_nonmanifold_FACEs( cutting_sheet ) == CUBIT_FAILURE )
    return CUBIT_FAILURE;

  // We need to check the normals on the cutting_sheet now against the normals
  // on the ext_target_BODY_ptr, since we can't be sure what the intersection
  // did with the normals.  Unfortunately, we can't use the
  // make_surf_normals_consistent function because we might have disconnected
  // FACEs that are part of the back side of cylinders in the cutting_sheet.
  // If we were to use make_surf_normals_consistent it would align those FACEs
  // normals with the front FACE, and the cutting sheet would be invalid,
  // because the remove_aligned_periodic_FACEs will fail to remove those FACEs.
  // We use a computationally intense method instead.
  if( align_normals( cutting_sheet, ext_target_BODY_ptr ) == CUBIT_FAILURE )
  {
    PRINT_ERROR( "Unable to make target shell normals consistent - invalid target\n" );
    return CUBIT_FAILURE;
  }

  // Now remove unwanted periodic FACEs on the back of the cutting sheet.  We
  // use the target normal as the basis.  Since we worked very hard to keep the
  // normals consisent in relation to the target surfaces the user picked, this
  // will be fairly robust for most cases.
  CubitVector tmp_vec = -target_norm;
  remove_aligned_periodic_FACEs( cutting_sheet, tmp_vec );

  // Make sure there are still some FACEs left!
  DLIList<FACE*> cutting_sheet_FACE_list;
  AQE->get_FACEs( cutting_sheet, cutting_sheet_FACE_list );
  if( cutting_sheet_FACE_list.size() == 0 )
  {
    PRINT_ERROR( "Unable to tweak to targets\n" );
    return CUBIT_FAILURE;
  }

  // Clean the cutting sheet before using it - remove unnecessary topology -
  // this makes booleans more robust.  Ignore result since not fatal.
  result = api_clean_entity( cutting_sheet );

  // Also remove any internal holes
  remove_holes( cutting_sheet );

  // Sheet direction should ALWAYS be pointing in the opposite direction as the
  // source_norm
  cutting_sheet_FACE_list.clean_out();
  AQE->get_FACEs( cutting_sheet, cutting_sheet_FACE_list );
  CubitVector sheet_norm = weighted_average_FACE_normal( cutting_sheet_FACE_list );
  if( sheet_norm.length() < 1e-12 )
  {
    PRINT_ERROR( "Unable to tweak to given surfaces\n" );
    return CUBIT_FAILURE;
  }
  if( sheet_norm % source_norm > 0.0 )
    api_reverse_body( cutting_sheet );

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : chop_off_with_sheet
// Member Type: PRIVATE
// Description: Chop off all the material in a BODY to one side of a sheet
//              body.  The material is removed on the side of the sheet body
//              according to the surface normals of the sheet body (material
//              removed in the opposite direction of the surface normals of the
//              sheet body, which must be consistent).
// Author     : Steve Storm
// Date       : 1/29/06
//=============================================================================
CubitStatus
AcisTweakTool::chop_off_with_sheet( BODY *BODY_ptr, BODY *sheet_BODY_ptr )
{
  outcome result;

  // First test to see if these intersect.  If they don't nothing will be done.
  if ( !AME->BODYs_interfering( BODY_ptr, sheet_BODY_ptr ) )
  {
    // BODIES don't intersect.  Clean up and exit.
    PRINT_WARNING("Cutting sheet does not intersect the original volume.\n");
    return CUBIT_FAILURE;
  }

  BODY *copied_sheet_BODY = AME->copy_BODY( sheet_BODY_ptr );

  // Set all faces to one sided
  result = api_body_to_1d( copied_sheet_BODY, false );
  if( !result.ok() )
  {
    api_delent( copied_sheet_BODY );
    AQE->ACIS_API_error( result );
    return CUBIT_FAILURE;
  }

  // Setup tweak attributes so we can preserve Cubit owners
  DLIList<FACE*> pre_FACE_list;
  DLIList<EDGE*> pre_EDGE_list;
  DLIList<VERTEX*> pre_VERTEX_list;
  DLIList<AcisBridge*> ab_FACE_list, ab_EDGE_list, ab_VERTEX_list;
  assign_tweak_attribs( BODY_ptr, "tweak_own", pre_FACE_list, ab_FACE_list,
    pre_EDGE_list, ab_EDGE_list, pre_VERTEX_list, ab_VERTEX_list );

  // Now we simply have to subtract the copied_sheet_BODY from BODY_ptr
  result = api_boolean( copied_sheet_BODY, BODY_ptr, NONREG_SUBTRACTION );

  // Replace Cubit owners
  reassign_cubit_owners_from_tweak_attribs( BODY_ptr, "tweak_own",
    pre_FACE_list, ab_FACE_list, pre_EDGE_list, ab_EDGE_list, pre_VERTEX_list,
    ab_VERTEX_list );

  // Remove tweak attributes
  remove_named_attribs( BODY_ptr, "tweak_own" );

  if( !result.ok() || BODY_ptr->lump() == NULL ||
      !is_closed_solid_body( BODY_ptr ) )
  {
    CubitBoolean is_error = CUBIT_TRUE;
    if( result.ok() &&
      ( BODY_ptr->lump() == NULL ||
      !is_closed_solid_body( BODY_ptr )))
    {
      // This is where the sheet just intersects with the shell
      // of the outside body.  Should do nothing here.
      is_error = CUBIT_FALSE;
    }
    if( !result.ok() )  //was already deleted if api_subtract was successfull
      api_delent( copied_sheet_BODY );

    if ( is_error )
    {
      if( !result.ok() )
        AQE->ACIS_API_error( result );
      else
        PRINT_ERROR( "problem trimming body\n" );
    }
    else
    {
      // Cutting surfaces graze volume
      PRINT_ERROR( "problem trimming body\n" );
    }
    return CUBIT_FAILURE;
  }

  AQE->clear_bounding_box( BODY_ptr );
  AQE->bounding_box( BODY_ptr );

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : get_overlap_area
// Member Type: PRIVATE
// Description: Determine the area of overlap of two BODIES.  Uses an ACIS
//              intersection boolean so is quite expensive.
// Author     : Steve Storm
// Date       : 2/10/06
//=============================================================================
CubitStatus
AcisTweakTool::get_overlap_area( BODY *BODY_ptr1, BODY *BODY_ptr2,
                                 double &overlap_area, double accuracy )
{
  overlap_area = 0.0;

  BODY *int_sheet;
  outcome result = api_boolean( BODY_ptr1, BODY_ptr2, INTERSECTION,
    NDBOOL_KEEP_BOTH, int_sheet );
  if( !result.ok() )
  {
    AQE->ACIS_API_error(result);
    return CUBIT_FAILURE;
  }

  // Get surface area of int_sheet
  double accur_achieved;
  result = api_ent_area( int_sheet, accuracy, overlap_area, accur_achieved );
  api_delent( int_sheet );
  if( !result.ok() )
  {
    AQE->ACIS_API_error(result);
    return CUBIT_FAILURE;
  }

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : copy_FACEs_into_sheet
// Member Type: PRIVATE
// Description: Copies the given FACEs
// Author     : Steve Storm
// Date       : 2/10/06
//=============================================================================
CubitStatus
AcisTweakTool::copy_FACEs_into_sheet( DLIList<FACE*> &FACE_list,
                                      BODY *&sheet_ptr,
                                      CubitBoolean heal )
{
  DLIList<BODY*> temp_BODY_list;
  FACE *FACE_ptr;
  FACE *temp_FACE_ptr;
  BODY *temp_BODY_ptr;
  int i;
  FACE_list.reset();
  for( i=FACE_list.size(); i--; )
  {
    FACE_ptr = FACE_list.get_and_step();

    temp_FACE_ptr = AME->make_FACE( FACE_ptr, CUBIT_FALSE );
    if( temp_FACE_ptr == NULL )
    {
      while( temp_BODY_list.size() ) api_delent( temp_BODY_list.pop() );
      return CUBIT_FAILURE;
    }

    // Get the BODY this new FACE is in
    temp_BODY_ptr = AQE->get_BODY_of_ENTITY( temp_FACE_ptr );
    if( temp_BODY_ptr == NULL )
    {
      while( temp_BODY_list.size() ) api_delent( temp_BODY_list.pop() );
      return CUBIT_FAILURE;
    }

    temp_BODY_list.append( temp_BODY_ptr );
  }

  // Unite all of the temp BODIES together
  temp_BODY_list.reset();
  outcome result;
  BODY *master = temp_BODY_list.extract();
  for( i=temp_BODY_list.size(); i--; )
  {
    temp_BODY_ptr = temp_BODY_list.extract();

    // Do the union of the master and the BODY_ptr.
    // If this is successful, the result is master and
    //   temp_BODY_ptr will be deleted
    result = api_boolean( temp_BODY_ptr, master, UNION );
    if( !result.ok() || (!master) )
    {
      AQE->ACIS_API_error(result);
      api_delent( temp_BODY_ptr );
      while( temp_BODY_list.size() ) api_delent( temp_BODY_list.pop() );
      if (master != NULL) api_delent(master);
      return CUBIT_FAILURE;
    }
  }

  // Heal the BODY if required
  if( heal )
    heal_BODY( master );

  // Clean (regularize) to remove unwanted topology
  api_clean_entity( master );

  sheet_ptr = master;

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : remove_per_nonmanifold_FACEs
// Member Type: PRIVATE
// Description: Removes cone, sphere or torus FACEs that are attached to a
//              nonmanifold EDGE.  The portion of the periodic that is tangent
//              to the attached FACE is kept.
// Author     : Steve Storm
// Date       : 2/2/06
//=============================================================================
CubitStatus
AcisTweakTool::remove_per_nonmanifold_FACEs( BODY *cutting_tool )
{
  outcome result;

  // List to store FACEs to remove
  DLIList<FACE*> remove_FACE_list;

  // Find nonmanifold edges, remove attached periodic surface on side we
  // don't want to keep.  We for sure want to remove these faces.
  //
  // Example: Diagram is an end view of a full cylinder with a plane attached.
  // The nonmanifold EDGE is on the right side.  In this example we would
  // remove the short portion of the cylinder between the periods.
  //     _____
  //   /       \
  //  /         \
  // |           |
  // |           . <-- nonmanifold EDGE (runs along cylinder length here)
  // |           |
  //  \         /|
  //   \ __.__ / |
  //             |
  //             |
  //             |
  //             |

  // Loop through EDGEs
  DLIList<EDGE*> EDGE_list;
  AQE->get_EDGEs( cutting_tool, EDGE_list );

  int i, j;
  EDGE *EDGE_ptr;
  EDGE_list.reset();
  for( i=EDGE_list.size(); i--; )
  {
    EDGE_ptr = EDGE_list.get_and_step();

    DLIList<FACE*> FACE_list;
    AQE->get_FACEs( EDGE_ptr, FACE_list );

    if( FACE_list.size() < 3 )
      continue;

    // EDGE is nonmanifold - it has 3 FACEs attached to it.  If two of the
    // FACEs are periodic types and one is not, then we will mark one of the
    // periodic types for deletion.

    // We can't handle this case
    if( FACE_list.size() > 3 )
    {
      PRINT_ERROR( "Ill-formed extended target encountered\n" );
      return CUBIT_FAILURE;
    }

    // Find the attached periodic surfaces
    FACE *per_FACE_ptr1 = NULL;
    FACE *per_FACE_ptr2 = NULL;
    FACE *other_FACE_ptr = NULL;
    FACE_list.reset();
    for( j=FACE_list.size(); j--; )
    {
      FACE *FACE_ptr = FACE_list.get_and_step();
      SURFACE *SURFACE_ptr = FACE_ptr->geometry();
      int type = SURFACE_ptr->identity();
      if( type==CONE_TYPE || type==SPHERE_TYPE || type== TORUS_TYPE )
      {
        if( per_FACE_ptr1 == NULL ){
          per_FACE_ptr1 = FACE_ptr;
          continue;
        }

        if( per_FACE_ptr2 == NULL )
        {
          per_FACE_ptr2 = FACE_ptr;
          continue;
        }

        if( per_FACE_ptr1 && per_FACE_ptr2 )
        {
          PRINT_ERROR( "Ill-formed extended target encountered\n" );
          return CUBIT_FAILURE;
        }
      }
      else if( other_FACE_ptr == NULL )
        other_FACE_ptr = FACE_ptr;
    }

    if( !(per_FACE_ptr1 && per_FACE_ptr2 && other_FACE_ptr) )
    {
      PRINT_ERROR( "Ill-formed extended target encountered\n" );
      return CUBIT_FAILURE;
    }

    // Figure out which periodic to remove.  Compare vectors tangent to and
    // pointing away from the FACEs at the middle of the EDGE.

    // Get mid position
    SPAposition acis_mid_pnt = EDGE_ptr->mid_pos();
    CubitVector mid_pnt( acis_mid_pnt.x(), acis_mid_pnt.y(), acis_mid_pnt.z() );

    CubitVector vec1, vec2, vec3;
    tangent_outdir( per_FACE_ptr1, EDGE_ptr, mid_pnt, vec1 );
    tangent_outdir( per_FACE_ptr2, EDGE_ptr, mid_pnt, vec2 );
    tangent_outdir( other_FACE_ptr, EDGE_ptr, mid_pnt, vec3 );

    // FACE that is in same direction as other_FACE gets removed.  If angle
    // between vectors is acute (dot product > 0) flag that FACE
    if( vec1 % vec3 > 0.0 )
      remove_FACE_list.append_unique( per_FACE_ptr1 );
    else if( vec2 % vec3 > 0.0 )
      remove_FACE_list.append_unique( per_FACE_ptr2 );
  }

  // Remove these FACEs
  for( i=remove_FACE_list.size(); i--; )
    api_remove_face( remove_FACE_list.get_and_step() );

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : tangent_outdir
// Member Type: PRIVATE
// Description: Given a FACE, EDGE and position, find the direction that is
//              normal to the EDGE and tangent to the FACE, pointing away from
//              the FACE boundary.
// Author     : Steve Storm
// Date       : 2/18/06
//=============================================================================
CubitStatus
AcisTweakTool::tangent_outdir( FACE *FACE_ptr, EDGE *EDGE_ptr,
                               CubitVector &pos, CubitVector &tangent_outvec )
{
  // Get the normal to the FACE at the position
  SPAposition acis_pos( pos.x(), pos.y(), pos.z() );
  SPAposition acis_pos_on_surf;
  SPAunit_vector acis_face_norm;
  (FACE_ptr->geometry()->equation()).point_perp( acis_pos,
    acis_pos_on_surf, acis_face_norm );

  // Adjust for the FACE sense
  if( FACE_ptr->sense() == REVERSED )
    acis_face_norm = -acis_face_norm;

  // Get a vector tangent to the EDGE, adjusted for the EDGE sense on the FACE
  SPAunit_vector acis_tangent_vec = (EDGE_ptr->geometry()->equation()).
    point_direction( acis_pos );
  if( EDGE_ptr->sense() == REVERSED )
    acis_tangent_vec = -acis_tangent_vec;

  COEDGE *COEDGE_ptr = EDGE_ptr->coedge( FACE_ptr );
  if( COEDGE_ptr->sense() == REVERSED )
    acis_tangent_vec = -acis_tangent_vec;

  // Cross edge tangent with face normal to find desired direction
  CubitVector face_norm( acis_face_norm.x(), acis_face_norm.y(),
    acis_face_norm.z() );
  CubitVector tangent_vec( acis_tangent_vec.x(), acis_tangent_vec.y(),
    acis_tangent_vec.z() );

  tangent_outvec = tangent_vec * face_norm;

  tangent_outvec.normalize();

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : remove_aligned_periodic_FACEs
// Member Type: PRIVATE
// Description: Removes cone, sphere or torus FACEs that have normals aligned
//              with the given basis vector.
// Author     : Steve Storm
// Date       : 2/2/06
//=============================================================================
CubitStatus
AcisTweakTool::remove_aligned_periodic_FACEs( BODY *cutting_tool,
                                              CubitVector &basis_vec )
{
  outcome result;
  int i;
  CubitVector location;
  DLIList<FACE*> FACE_list;
  AQE->get_FACEs( cutting_tool, FACE_list );
  FACE_list.reset();
  for( i=FACE_list.size(); i--; )
  {
    FACE *FACE_ptr = FACE_list.get_and_step();

    SURFACE *SURFACE_ptr = FACE_ptr->geometry();
    int type = SURFACE_ptr->identity();
    if( type==CONE_TYPE || type==SPHERE_TYPE || type== TORUS_TYPE )
    {
      CubitVector FACE_norm = FACE_normal( FACE_ptr, location );

      // If dot product greater than zero, angle is acute
      if( FACE_norm % basis_vec > 0.0 )
      {
        result = api_remove_face( FACE_ptr );
        if( !result.ok() )
        {
          AQE->ACIS_API_error( result );
          return CUBIT_FAILURE;
        }
      }
    }
  }

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : heal_BODY
// Member Type: PRIVATE
// Description: Convenience function to use the ACIS healer to heal a BODY
// Author     : Steve Storm
// Date       : 2/20/06
//=============================================================================
CubitStatus
AcisTweakTool::heal_BODY( BODY *BODY_ptr )
{
  outcome result;
  result = api_hh_init_body_for_healing( BODY_ptr );
  if( !result.ok() )
  {
    AQE->ACIS_API_error( result );
    return CUBIT_FAILURE;
  }

  result = api_hh_auto_heal( BODY_ptr );
  if( !result.ok() )
  {
    AQE->ACIS_API_error( result );
    api_hh_end_body_for_healing( BODY_ptr );
    return CUBIT_FAILURE;
  }

  api_hh_end_body_for_healing( BODY_ptr );

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : align_normals
// Member Type: PRIVATE
// Description: Align the normals of the child's FACEs to the master FACEs.  It
//              is assumed the child is a subset of the master (i.e., the
//              master is an extended BODY).
// Author     : Steve Storm
// Date       : 2/20/06
//=============================================================================
CubitStatus
AcisTweakTool::align_normals( BODY *child, BODY *master )
{
  DLIList<FACE*> child_FACE_list;
  AQE->get_FACEs( child, child_FACE_list );

  DLIList<FACE*> master_FACE_list;
  AQE->get_FACEs( master, master_FACE_list );

  int i;
  CubitVector loc;
  FACE *child_FACE_ptr = NULL;
  for( i=child_FACE_list.size(); i--; )
  {
    child_FACE_ptr = child_FACE_list.get_and_step();

    CubitVector child_norm = FACE_normal( child_FACE_ptr, loc );

    FACE *master_FACE_ptr = find_overlap_FACE( child_FACE_ptr, master_FACE_list );

    if( !master_FACE_ptr )
      return CUBIT_FAILURE;

    // Align the normals
    CubitVector master_norm = FACE_normal( master_FACE_ptr, loc, CUBIT_FALSE );
    if( child_norm % master_norm < 0.0 )
      api_reverse_face( child_FACE_ptr );
  }

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : find_overlap_FACE
// Member Type: PRIVATE
// Description: Find the FACE in FACE_list that overlaps with FACE_ptr.
//              Returns NULL if no overlap FACE found.
// Author     : Steve Storm
// Date       : 2/20/06
//=============================================================================
FACE *
AcisTweakTool::find_overlap_FACE( FACE *FACE_ptr, DLIList<FACE*> &FACE_list )
{
  int i;
  BODY *tmp_BODY_ptr1;
  FACE *face_list[1];
  face_list[0] = FACE_ptr;
  api_sheet_from_ff( 1, face_list, tmp_BODY_ptr1 );
  api_body_to_2d( tmp_BODY_ptr1 );

  FACE *FACE_ptr2;
  FACE *match_FACE_ptr = NULL;
  double overlap_area;
  for( i=FACE_list.size(); i--; )
  {
    FACE_ptr2 = FACE_list.get_and_step();
    if( !FACE_ptr2 )
      continue;
    BODY *tmp_BODY_ptr2;
    face_list[0] = FACE_ptr2;
    api_sheet_from_ff( 1, face_list, tmp_BODY_ptr2 );
    api_body_to_2d( tmp_BODY_ptr2 );
    get_overlap_area( tmp_BODY_ptr1, tmp_BODY_ptr2, overlap_area );
    api_delent( tmp_BODY_ptr2 );
    if( overlap_area > 1e-5 )
    {
      match_FACE_ptr = FACE_ptr2;
      break;
    }
  }

  api_delent( tmp_BODY_ptr1 );

  return match_FACE_ptr;
}

//=============================================================================
// Function   : make_surf_normals_consistent
// Member Type: PRIVATE
// Description: Make the surface normals in the given sheet BODY consistent.
//              The optional seed FACE must be part of the BODY.  If none is
//              given, a random FACE is selected.
// Author     : Steve Storm
// Date       : 2/20/06
//=============================================================================
CubitStatus
AcisTweakTool::make_surf_normals_consistent( BODY *BODY_ptr,
                                             FACE *seed_FACE_ptr )
{
  int i;
  DLIList<FACE*> patch_list;
  DLIList<FACE*> FACE_list;
  AQE->get_FACEs( BODY_ptr, FACE_list );

  if( seed_FACE_ptr == NULL )
    seed_FACE_ptr = get_seed_FACE( FACE_list );

  api_add_generic_named_attribute( seed_FACE_ptr, "tweak_seed", 1, SplitCopy );
  patch_list.append( seed_FACE_ptr );

  int c = 0;

  while( patch_list.size() < FACE_list.size() )
  {
    // Find all neighbors.  Note patch_list will change size within the loop.
    for( i=c; i<patch_list.size(); i++ )
    {
      if( append_neighbors( patch_list ) == CUBIT_FAILURE )
      {
        remove_named_attribs( BODY_ptr, "tweak_seed" );
        return CUBIT_FAILURE;
      }

      patch_list.step();
    }

    // Include disconnected FACE patches
    if( patch_list.size() < FACE_list.size() )
    {
      c = patch_list.size();
      ATTRIB_GEN_NAME *ret_att;
      FACE *FACE_ptr;
      DLIList<FACE*> processed_FACE_list;

      FACE_list.reset();
      for( i=FACE_list.size(); i--; )
      {
        FACE_ptr = FACE_list.get_and_step();

        api_find_named_attribute( FACE_ptr, "tweak_seed", ret_att );
        if( ret_att )
          processed_FACE_list.append( FACE_ptr );
      }

      CubitVector avg_norm = weighted_average_FACE_normal( processed_FACE_list );

      seed_FACE_ptr = get_seed_FACE( FACE_list );
      api_add_generic_named_attribute( seed_FACE_ptr, "tweak_seed", 1, SplitCopy );
      patch_list.append( seed_FACE_ptr );

      CubitVector loc;
      CubitVector seed_norm = FACE_normal( seed_FACE_ptr, loc );

      if( seed_norm % avg_norm < 0.0 )
        api_reverse_face( seed_FACE_ptr );
    }
  }

  remove_named_attribs( BODY_ptr, "tweak_seed" );

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : append_neighbors
// Member Type: PRIVATE
// Description: Helper function for make_surf_normals_consistent
// Author     : Steve Storm
// Date       : 2/20/06
//=============================================================================
CubitStatus
AcisTweakTool::append_neighbors( DLIList<FACE*> &FACE_list )
{
  int i, j;
  ATTRIB_GEN_NAME *ret_att;

  FACE *seed_FACE_ptr = FACE_list.get();

  // Get the edges
  DLIList<EDGE*> EDGE_list;
  AQE->get_EDGEs( seed_FACE_ptr, EDGE_list );

  EDGE *EDGE_ptr;
  FACE *FACE_ptr;
  COEDGE *seed_COEDGE_ptr;
  for( i=0; i<EDGE_list.size(); i++ )
  {
    EDGE_ptr = EDGE_list.get_and_step();

    seed_COEDGE_ptr = EDGE_ptr->coedge( seed_FACE_ptr );

    ENTITY_LIST COEDGES;
    api_get_coedges( EDGE_ptr, COEDGES );
    DLIList<COEDGE*> COEDGE_list;
    ENTITY *ENTITY_ptr;
    COEDGES.init();
    while( (ENTITY_ptr = COEDGES.next() ) != NULL )
      COEDGE_list.append( (COEDGE *)ENTITY_ptr );
    COEDGES.clear();

    if( COEDGE_list.size() == 1 )
      continue;
    if( COEDGE_list.size() == 2 )
    {
      // Normal case
      COEDGE *COEDGE_ptr2 = COEDGE_list.get_and_step();
      if( COEDGE_ptr2 == seed_COEDGE_ptr )
        COEDGE_ptr2 = COEDGE_list.get();

      DLIList<FACE*> temp_FACE_list;
      AQE->get_FACEs( COEDGE_ptr2, temp_FACE_list );
      if( temp_FACE_list.size() > 1 )
      {
        PRINT_ERROR( "Found more than one FACE attached to a COEDGE\n" );
        return CUBIT_FAILURE;
      }
      FACE_ptr = temp_FACE_list.get();

      // Don't consider FACEs that are already in the list
      api_find_named_attribute( FACE_ptr, "tweak_seed", ret_att );
      if( ret_att )
        continue;

      // If the sense's are the same, the FACE needs to be reversed
      if( seed_COEDGE_ptr->sense() == COEDGE_ptr2->sense() )
        api_reverse_face( FACE_ptr );

      FACE_list.append( FACE_ptr );
      api_add_generic_named_attribute( FACE_ptr, "tweak_seed", 1, SplitCopy );
    }
    else if( COEDGE_list.size() == 3 )
    {
      // Nonmanifold.  Check for special case - periodic surface intersecting
      // (usually tangent) another surface.  For now don't allow any other
      // special cases.

      // Find the attached periodic surfaces
      FACE *FACE_ptr2 = NULL;
      int num_per = 0;
      DLIList<FACE*> attached_FACE_list;
      AQE->get_FACEs( EDGE_ptr, attached_FACE_list );
      attached_FACE_list.reset();
      for( j=attached_FACE_list.size(); j--; )
      {
        FACE_ptr = attached_FACE_list.get_and_step();
        SURFACE *SURFACE_ptr = FACE_ptr->geometry();
        int type = SURFACE_ptr->identity();
        if( type==CONE_TYPE || type==SPHERE_TYPE || type== TORUS_TYPE )
          num_per++;

        if( FACE_ptr == seed_FACE_ptr )
          continue;

        api_find_named_attribute( FACE_ptr, "tweak_seed", ret_att );
        if( ret_att )
          continue;

        FACE_ptr2 = FACE_ptr;
      }

      if( FACE_ptr2 == NULL )
        continue;
      else if( attached_FACE_list.size() == 2 && num_per != 1 )
      {
        PRINT_ERROR( "Extended target cannot contain nonmanifold geometry\n" );
        return CUBIT_FAILURE;
      }
      else if( attached_FACE_list.size() == 3 && num_per != 2 )
      {
        PRINT_ERROR( "Extended target cannot contain nonmanifold geometry\n" );
        return CUBIT_FAILURE;
      }

      // Compare normal vectors at the mid position of the EDGE

      // Get mid position
      SPAposition acis_mid_pnt = EDGE_ptr->mid_pos();
      CubitVector mid_pnt( acis_mid_pnt.x(), acis_mid_pnt.y(), acis_mid_pnt.z() );

      CubitVector vec_seed, vec2;
      vec_seed = FACE_normal( seed_FACE_ptr, mid_pnt, CUBIT_FALSE );
      vec2 = FACE_normal( FACE_ptr2, mid_pnt, CUBIT_FALSE );

      // Check angle between FACEs
      if( vec_seed % vec2 < 0.0 )
        api_reverse_face( FACE_ptr2 );

      FACE_list.append( FACE_ptr2 );
      api_add_generic_named_attribute( FACE_ptr2, "tweak_seed", 1, SplitCopy );
    }
    else if( COEDGE_list.size() > 3 )
    {
      PRINT_ERROR( "Extended target cannot contain nonmanifold geometry\n" );
      return CUBIT_FAILURE;
    }
    else
    {
      PRINT_ERROR( "Targets appear to be ill-formed\n" );
      return CUBIT_FAILURE;
    }
  }

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : get_seed_FACE
// Member Type: PRIVATE
// Description: Get a seed FACE from the given list of FACEs.  Only a seed FACE
//              without a named "tweak_seed" attribute on it will be returned.
//              Also, it may not be of advantage, but we try to find a FACE that
//              is not of type CONE, TORUS or SPHERE.  If that is not possible,
//              we just return the first FACE in the list with no "tweak_seed"
//              attribute on it.
// Author     : Steve Storm
// Date       : 2/20/06
//=============================================================================
FACE *
AcisTweakTool::get_seed_FACE( DLIList<FACE*> &FACE_list )
{
  int i;
  FACE *seed_FACE_ptr = NULL;
  ATTRIB_GEN_NAME *ret_att;

  // Attempt to get a non-periodic seed FACE first
  FACE *FACE_ptr;
  FACE_list.reset();
  for( i=FACE_list.size(); i--; )
  {
    FACE_ptr = FACE_list.get_and_step();
    SURFACE *SURFACE_ptr = FACE_ptr->geometry();
    int type = SURFACE_ptr->identity();
    if( type==CONE_TYPE || type==SPHERE_TYPE || type== TORUS_TYPE )
      continue;

    api_find_named_attribute( FACE_ptr, "tweak_seed", ret_att );

    if( ret_att )
      continue;

    seed_FACE_ptr = FACE_ptr;
    break;
  }

  if( seed_FACE_ptr == NULL )
  {
    FACE_list.reset();
    for( i=FACE_list.size(); i--; )
    {
      FACE_ptr = FACE_list.get_and_step();

      api_find_named_attribute( FACE_ptr, "tweak_seed", ret_att );

      if( ret_att )
        continue;

      seed_FACE_ptr = FACE_ptr;
      break;
    }
  }

  return seed_FACE_ptr;
}

//=============================================================================
// Function   : draw_tweak_preview_omt
// Member Type: PRIVATE
// Description: Draw the preview EDGEs for tweak offset, move and target.  The
//              method is to draw:
//              For sheet bodies:
//              -----------------
//               - EDGEs without owners and EDGEs attached to them
//               - EDGEs attached to given (by owner attribute) EDGEs 
//              For solids:
//              -----------
//               - EDGEs on FACEs without owners and EDGEs attached to them
//               - EDGEs attached to given (by owner attribute) FACEs
//                 
// Author     : Steve Storm
// Date       : 2/20/06
//=============================================================================
CubitStatus
AcisTweakTool::draw_tweak_preview_omt( BODY *BODY_ptr,
                                       CubitBoolean flush,
                                       DLIList<AcisBridge*> *ab_list )
{
  DLIList<FACE*> FACE_list;
  AQE->get_FACEs( BODY_ptr, FACE_list );
  if( !FACE_list.size() )
    return CUBIT_FAILURE;

  int i, j, k;
  DLIList<EDGE*> draw_EDGE_list;
  AcisBridge *ab_ptr = NULL;
  EDGE *EDGE_ptr;

  FACE *FACE_ptr = FACE_list.get_and_step();
  if( FACE_ptr->sides() == DOUBLE_SIDED )
  {
    // Sheet BODY
    // Draw EDGEs attached to EDGEs without Cubit owner attributes on them -
    // - this will draw the targets and "side" EDGEs too.
    // Note special case where we need to check if source and/or target still
    // exists - source could be tweaked to an existing EDGE on the same BODY.
    // For some reason the resulting EDGE at the target location inherits the
    // owner attributes of the target EDGE.
    DLIList<EDGE*> EDGE_list;
    AQE->get_EDGEs( BODY_ptr, EDGE_list );

    for( i=EDGE_list.size(); i--; )
    {
      EDGE_ptr = EDGE_list.get_and_step();
      ab_ptr = ATTRIB_CUBIT_OWNER::cubit_owner( EDGE_ptr );
      if( !ab_ptr || (ab_list && ab_list->is_in_list( ab_ptr )) )
      {
        // Get VERTICEs on this EDGE
        DLIList<VERTEX*> VERTEX_list;
        AQE->get_VERTICEs( EDGE_ptr, VERTEX_list );

        VERTEX *VERTEX_ptr;
        for( j=VERTEX_list.size(); j--; )
        {
          VERTEX_ptr = VERTEX_list.get_and_step();
          DLIList<EDGE*> att_EDGE_list;
          AQE->get_EDGEs( VERTEX_ptr, att_EDGE_list );

          for( k=att_EDGE_list.size(); k--; )
          {
            EDGE_ptr = att_EDGE_list.get_and_step();
            draw_EDGE_list.append_unique( EDGE_ptr );
          }
        }
      }
    }
  }
  else
  {
    // Solid BODY
    // Draw EDGEs attached to FACEs without Cubit owner attributes on them -
    // - this will draw the targets and "side" EDGEs too.  
    // Note special case where we need to check if source still exists - source
    // could be tweaked to an existing FACE on the same BODY.  For some reason
    // the resulting FACE at the target location inherits the owner attributes
    // of one of the source FACEs.
    for( i=FACE_list.size(); i--; )
    {
      FACE_ptr = FACE_list.get_and_step();
      ab_ptr = ATTRIB_CUBIT_OWNER::cubit_owner( FACE_ptr );
      if( !ab_ptr || (ab_list && ab_list->is_in_list( ab_ptr )) )
      {
        // Get VERTICEs on this FACE
        DLIList<VERTEX*> VERTEX_list;
        AQE->get_VERTICEs( FACE_ptr, VERTEX_list );

        VERTEX *VERTEX_ptr;
        for( j=VERTEX_list.size(); j--; )
        {
          VERTEX_ptr = VERTEX_list.get_and_step();
          DLIList<EDGE*> EDGE_list;
          AQE->get_EDGEs( VERTEX_ptr, EDGE_list );

          for( k=EDGE_list.size(); k--; )
          {
            EDGE_ptr = EDGE_list.get_and_step();
            draw_EDGE_list.append_unique( EDGE_ptr );
          }
        }
      }
    }
  }

  for( i=draw_EDGE_list.size(); i--; )
  {
    EDGE_ptr = draw_EDGE_list.get_and_step();
    AcisDrawTool::instance()->draw_EDGE( EDGE_ptr, CUBIT_BLUE );
  }

  if( flush )
    GfxPreview::flush();

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : tag_tweak_remove_faces_for_preview
// Member Type: PRIVATE
// Description: Add an attribute ("tweak_preview") to the required FACEs for a 
//              tweak remove.  For the preview we draw the EDGEs on the
//              surviving FACEs adjoining those that are removed.  Input is the
//              copied BODY, along with the FACEs being removed.
// Author     : Steve Storm
// Date       : 2/20/06
//=============================================================================
CubitStatus
AcisTweakTool::tag_tweak_remove_FACEs_for_preview( BODY *BODY_ptr, 
                                             DLIList<FACE*> &remove_FACE_list )
{
  // Add an attribute to adjoining FACEs we are removing.  Then we can
  // draw these FACEs for the preview.  Don't worry about putting
  // attributes on the FACEs we are removing since these FACEs won't exist
  // when we are done anyway.
  int i, j, k;
  FACE *FACE_ptr;
  ATTRIB_GEN_NAME *ret_att;
  for( i=remove_FACE_list.size(); i--; )
  {
    FACE_ptr = remove_FACE_list.get_and_step();
    DLIList<EDGE*> EDGE_list;
    AQE->get_EDGEs( FACE_ptr, EDGE_list );

    EDGE *EDGE_ptr;
    for( j=EDGE_list.size(); j--; )
    {
      EDGE_ptr = EDGE_list.get_and_step();
      DLIList<FACE*> att_FACE_list;
      AQE->get_FACEs( EDGE_ptr, att_FACE_list );

      for( k=att_FACE_list.size(); k--; )
      {
        FACE_ptr = att_FACE_list.get_and_step();

        api_find_named_attribute( FACE_ptr, "tweak_preview", ret_att );
        if( !ret_att )
          api_add_generic_named_attribute( FACE_ptr, "tweak_preview", 1, SplitCopy );
      }
    }
  }

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : tag_tweak_remove_faces_for_preview
// Member Type: PRIVATE
// Description: Add an attribute ("tweak_preview") to the required FACEs for a 
//              tweak remove.  For the preview we draw the EDGEs on the
//              surviving FACEs adjoining those that are removed.  Input is the
//              copied BODY, along with the FACEs being removed.
// Author     : Steve Storm
// Date       : 2/20/06
//=============================================================================
CubitStatus
AcisTweakTool::tag_tweak_remove_FACEs_for_preview( BODY *BODY_ptr, 
                                             DLIList<EDGE*> &remove_EDGE_list )
{
  // Add an attribute to adjoining FACEs we are removing.  Then we can
  // draw these FACEs for the preview.
  int i, j;
  EDGE *EDGE_ptr;
  FACE *FACE_ptr;
  ATTRIB_GEN_NAME *ret_att;
  for( i=remove_EDGE_list.size(); i--; )
  {
    EDGE_ptr = remove_EDGE_list.get_and_step();
    DLIList<FACE*> att_FACE_list;
    AQE->get_FACEs( EDGE_ptr, att_FACE_list );

    for( j=att_FACE_list.size(); j--; )
    {
      FACE_ptr = att_FACE_list.get_and_step();

      api_find_named_attribute( FACE_ptr, "tweak_preview", ret_att );
      if( !ret_att )
        api_add_generic_named_attribute( FACE_ptr, "tweak_preview", 1, SplitCopy );
    }
  }

  return CUBIT_SUCCESS;
}

//=============================================================================
// Function   : draw_tweak_preview_tagged_FACEs
// Member Type: PRIVATE
// Description: Draw EDGEs on FACEs with the "tweak_preview" attribute on them.
// Author     : Steve Storm
// Date       : 2/20/06
//=============================================================================
CubitStatus
AcisTweakTool::draw_tweak_preview_tagged_FACEs( BODY *BODY_ptr, 
                                                CubitBoolean flush )
{
  int i, j;
  FACE *FACE_ptr;
  ATTRIB_GEN_NAME *ret_att;
  DLIList<EDGE*> draw_EDGE_list;
  DLIList<FACE*> FACE_list;
  AQE->get_FACEs( BODY_ptr, FACE_list );
  for( i=FACE_list.size(); i--; )
  {
    FACE_ptr = FACE_list.get_and_step();
    api_find_named_attribute( FACE_ptr, "tweak_preview", ret_att );
    if( ret_att )
    {
      api_remove_generic_named_attribute( FACE_ptr, "tweak_preview" );

      DLIList<EDGE*> EDGE_list;
      AQE->get_EDGEs( FACE_ptr, EDGE_list );

      for( j=EDGE_list.size(); j--; )
        draw_EDGE_list.append_unique( EDGE_list.get_and_step() );
    }
  }

  for( i=draw_EDGE_list.size(); i--; )
    AcisDrawTool::instance()->draw_EDGE( draw_EDGE_list.get_and_step(),
    CUBIT_BLUE );

  if( flush )
    GfxPreview::flush();

  return CUBIT_SUCCESS;
}


