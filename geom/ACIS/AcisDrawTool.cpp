//- Class:       AcisDrawTool
//- Description: Functions to draw ACIS geometry (useful for debugging or
//-              in previewing operations for the user)
// 
//- Owner: Steve Storm      
//- Checked by:
//- Version: $Id:

// ********** BEGIN ACIS INCLUDES             **********
#include "entity.hxx"
#include "body.hxx"
#include "lump.hxx"
#include "shell.hxx"
#include "loop.hxx"
#include "coedge.hxx"
#include "edge.hxx"
#include "vertex.hxx"
#include "kernapi.hxx"
#include "position.hxx"
#include "point.hxx"
#include "SurfaceACIS.hpp"
#include "CurveACIS.hpp"
// ********** END ACIS INCLUDES               **********

// ********** BEGIN CUBIT INCLUDES            **********
#include "AcisDrawTool.hpp"
#include "BodySM.hpp"
#include "Lump.hpp"
#include "Surface.hpp"
#include "Curve.hpp"
#include "Point.hpp"
#include "AcisQueryEngine.hpp"
#include "DLIList.hpp"
#include "GMem.hpp"
#include "GfxPreview.hpp"
// ********** END CUBIT INCLUDES              **********

AcisDrawTool* AcisDrawTool::instance_ = 0;

// Method: instance
// provides access to the unique model for this execution.
// sets up this instance on first access
AcisDrawTool* AcisDrawTool::instance()
{
  if (instance_ == 0) {
    instance_ = new AcisDrawTool();
  }
  return instance_;
}

AcisDrawTool::AcisDrawTool() 
{ 
}

AcisDrawTool::~AcisDrawTool()
{
}

CubitStatus
AcisDrawTool::draw_ENTITY( ENTITY *ENTITY_ptr, int color, CubitBoolean tessellate,
                           CubitBoolean flush )
{
  outcome result;
  AcisQueryEngine *AQE = AcisQueryEngine::instance();

  if( IS_ENTITY_TYPE( ENTITY_ptr, BODY ) ||
      IS_ENTITY_TYPE( ENTITY_ptr, LUMP ) ||
      IS_ENTITY_TYPE( ENTITY_ptr, SHELL ) )
  {
    DLIList<FACE*> FACE_list;
    AQE->get_FACEs( ENTITY_ptr, FACE_list );
    int i;
    for( i=FACE_list.size(); i--; )
      draw_FACE( FACE_list.get_and_step(), color, tessellate );
    if( flush )
      GfxPreview::flush();
    return CUBIT_SUCCESS;
  }
  else if (IS_ENTITY_TYPE( ENTITY_ptr, FACE ) )
  {
    return draw_FACE( (FACE*)ENTITY_ptr, color, tessellate, flush );
  }
  else if( IS_ENTITY_TYPE( ENTITY_ptr, LOOP ) ||
           IS_ENTITY_TYPE( ENTITY_ptr, COEDGE ) )
  {
    DLIList<EDGE*> EDGE_list;
    AQE->get_EDGEs( ENTITY_ptr, EDGE_list );
    int i;
    for( i=EDGE_list.size(); i--; )
      draw_EDGE( EDGE_list.get_and_step(), color );
    if( flush )
      GfxPreview::flush();
    return CUBIT_SUCCESS;
  }
  else if (IS_ENTITY_TYPE( ENTITY_ptr, EDGE ) )
  {
    return draw_EDGE( (EDGE*)ENTITY_ptr, color, flush );
  }
  else if (IS_ENTITY_TYPE( ENTITY_ptr, VERTEX ) )
  {
    return draw_VERTEX( (VERTEX*)ENTITY_ptr, color, flush );
  }
  else
  {
    PRINT_ERROR( "Unsupported entity type specified - cannot draw\n" );
    return CUBIT_FAILURE;
  }

  return CUBIT_SUCCESS;
}

CubitStatus
AcisDrawTool::draw_EDGE( EDGE *EDGE_ptr, int color, CubitBoolean flush )
{
  outcome result;
  AcisQueryEngine *AQE = AcisQueryEngine::instance();

  // Copy so we don't interfere with any possible Cubit attributes
  ENTITY *copied_ENTITY_ptr;
  result = api_copy_entity_contents( EDGE_ptr, copied_ENTITY_ptr );
  if( !result.ok() )
  {
    PRINT_ERROR( "Unable to prepare curve for display\n" );
    return CUBIT_FAILURE;
  }

  Curve *curve_ptr = AQE->populate_topology_bridges( (EDGE *)copied_ENTITY_ptr );

  GMem g_mem;
  int num_points;

  // get the graphics
  CubitStatus stat;
  stat = AQE->get_graphics( curve_ptr, num_points, &g_mem );

  if (stat==CUBIT_FAILURE || num_points == 0)
  {
    PRINT_ERROR("Unable to tessellate a curve for display\n" );
    AQE->delete_solid_model_entities( curve_ptr );
    return CUBIT_FAILURE;
  }
  else
  {
    // Draw the polyline
    GfxPreview::draw_polyline( g_mem.point_list(), g_mem.pointListCount, color );
  }

  AQE->delete_solid_model_entities( curve_ptr );

  if( flush )
    GfxPreview::flush();

  return CUBIT_SUCCESS;
}

CubitStatus 
AcisDrawTool::draw_surface( Surface *surface, int color, CubitBoolean tessellate )
{
  SurfaceACIS *surf_acis = CAST_TO( surface, SurfaceACIS );
  FACE *face1 = surf_acis->get_FACE_ptr();
  draw_FACE( face1, color, tessellate );
  return CUBIT_SUCCESS;
}

CubitStatus 
AcisDrawTool::draw_curve( Curve *curve, int color )
{
  CurveACIS *curve_acis = CAST_TO( curve, CurveACIS );
  EDGE *edge1 = curve_acis->get_EDGE_ptr();
  draw_EDGE( edge1, color );
  return CUBIT_SUCCESS;
}

CubitStatus
AcisDrawTool::draw_FACE( FACE *FACE_ptr, int color, CubitBoolean tessellate, 
                         CubitBoolean flush )
{
  AcisQueryEngine *AQE = AcisQueryEngine::instance();

  if( tessellate )
  {
    int num_tris, num_pnts, num_facets;
    GMem g_mem;
    unsigned short norm_tol = 10;
    double dist_tol = -1.0;
    AQE->get_graphics( FACE_ptr, num_tris, num_pnts, num_facets,
      &g_mem, norm_tol, dist_tol );

    if(num_tris < 1)
    {
      // decrease tolerance and try again (we can get this for small features)
      norm_tol /= 2;
      AQE->get_graphics( FACE_ptr, num_tris, num_pnts, num_facets, &g_mem, 
        norm_tol, dist_tol);
    }

    if(num_tris < 1)
    {
      // lets give up 
      PRINT_ERROR( "Unable to tessellate surface for display\n" );
      return CUBIT_FAILURE;
    }

    // Draw the triangles
    GPoint p[3];
    GPoint* plist = g_mem.point_list();
    int* facet_list = g_mem.facet_list();
    int i, c = 0;
    for( i=0; i<num_tris; i++ )
    {
      p[0] = plist[facet_list[++c]];
      p[2] = plist[facet_list[++c]];
      p[1] = plist[facet_list[++c]]; 
      c++;
      GfxPreview::draw_tri( p, color );
    }
  }
  else
  {
    // Draw curves
    DLIList<EDGE*> EDGE_list;
    AQE->get_EDGEs( FACE_ptr, EDGE_list );
    int i;
    for( i=EDGE_list.size(); i--; )
      draw_EDGE( EDGE_list.get_and_step(), color );
  }

  if( flush )
    GfxPreview::flush();

  return CUBIT_SUCCESS;
}

CubitStatus
AcisDrawTool::draw_VERTEX( VERTEX *VERTEX_ptr, int color, CubitBoolean flush )
{
  const SPAposition &coords = VERTEX_ptr->geometry()->coords();
  GfxPreview::draw_point( coords.x(), coords.y(), coords.z(), color );
  return CUBIT_SUCCESS;
}
