//-------------------------------------------------------------------------
// Filename      : SplitSurfaceVirtual.cpp
//
// Purpose       : Split a surface and create virtual surfaces
//
// Special Notes : 
//
// Creator       : KGM
//
// Creation Date : 30-Jun-2005
//-------------------------------------------------------------------------

#include "SplitSurfaceVirtual.hpp"
#include "SplitSurfaceTool.hpp"
#include "RefVertex.hpp"
#include "RefEdge.hpp"
#include "RefFace.hpp"
#include "Curve.hpp"
#include "Surface.hpp"
#include "GeometryQueryTool.hpp"
#include "PartitionTool.hpp"
#include "CubitMessage.hpp"
#include "CubitUtil.hpp"
#include "DLIList.hpp"
#include "GfxDebug.hpp"
#include "GMem.hpp"

CubitStatus                                
SplitSurfaceVirtual::split_surface_virtual( RefFace *ref_face_ptr,
                                 DLIList<CubitVector*> &locations,
                                 DLIList<DLIList<CubitVector*>*> &vec_lists )
{
  // find the splitting curves. 
  SplitSurfaceTool sst;
  DLIList<Curve*> curve_list;
  CubitStatus err  = sst.calculate_split_curves( ref_face_ptr, locations, 
                                             vec_lists, curve_list );

  GMem gmem;
  int  num_points;
  DLIList<CubitVector*> segments;

  // loop over all the curves 
  int i;
  for (i = 0; i < curve_list.size(); i++)
  {
    Curve* curve_ptr = curve_list.get_and_step();

    // get the curve facets
    err = curve_ptr->get_geometry_query_engine()->
                          get_graphics( curve_ptr, num_points, &gmem );

    // load the graphics points into a CubitVector for insert_curve
    --num_points;  
    int j;
    for (j = 0; j < num_points; j++)
    {
      const GPoint& p = gmem.point_list()[j];
      segments.append( new CubitVector( p.x, p.y, p.z ) );
    }
  }

  // add last point to the end of the list
  const GPoint& p = gmem.point_list()[num_points];
  segments.append( new CubitVector(p.x, p.y, p.z ) );

  // now partition the surface
  DLIList<RefEdge*> new_edges;
  RefFace* new_surf = PartitionTool::instance()->insert_edge( ref_face_ptr, segments, 
                                                                false, new_edges );
  //DLIList<Curve*> new_curves;
  //Surface* new_surf = PartitionEngine::instance().insert_curve( old_surf, segments, new_curves );

  // called routines allocate memory.  Clean up for them.
  while( curve_list.size() ) 
    delete curve_list.pop();

  //while( new_curves.size() ) 
  //  delete new_curves.pop();

  return CUBIT_SUCCESS;
}


//===============================================================================
// Function   : split_surfaces
// Member Type: PUBLIC
// Description: Split a chain of surfaces into one or more pieces
// Author     : Steve Storm (CAT)
// Date       : 01/04
//===============================================================================
CubitStatus 
SplitSurfaceVirtual::split_surfaces_virtual( DLIList<RefFace*> &ref_face_list, 
                                             int num_segs, 
                                             double fraction,
                                             double distance,
                                             RefEdge *from_curve_ptr,
                                             DLIList<RefVertex*> &corner_vertex_list,
                                             DLIList<RefVertex*> &through_vertex_list,
                                             RefEdge *curve_dir_ptr,
                                             CubitBoolean preview_flg,
                                             CubitBoolean create_ref_edges_flg )
{
  // Get parent bodies - all surfs must be from same body
  int i;
  DLIList<Body*> old_body_list;
  RefFace *ref_face_ptr;
  ref_face_list.reset();
  for( i=ref_face_list.size(); i--; )
  {
    ref_face_ptr = ref_face_list.get_and_step();

    DLIList<Body*> body_list;
    ref_face_ptr->bodies( body_list );
    old_body_list.merge_unique( body_list );
  }

  if( old_body_list.size() > 1 )
  {
    PRINT_ERROR( "This operation requires all surfaces to be from the same volume\n" );
    // Note: this restriction could be pretty easily lifted by sorting the 
    //       input lists and calling the SplitSurfaceTool separately for each set of
    //       surfaces on each body.
    return CUBIT_FAILURE;
  }

  //bad geom with no body -- dont try to imprint this...  
  //quick and dirty fix by (aga@cat|1/7/04)
  if( old_body_list.size() < 1 )
  {
    PRINT_ERROR( "A surface is not contained within a parent body.\n"
      "       It cannot be split.\n");
    return CUBIT_FAILURE;
  }
 /* KGM -- need this? 
  // Check for virtual geometry
  if ( contains_intermediate_geometry(ref_face_list) )
  {
    PRINT_ERROR("SPLITTING surfaces containing virtual geometry is not\n"
      "       allowed. Delete virtual geometry on these surfaces\n"
      "       before operation.\n" );
    return CUBIT_FAILURE;
  }
  */

  // Make sure all surfaces are from same geometry engine
  DLIList<RefEntity*> ref_ent_list;
  CAST_LIST_TO_PARENT(ref_face_list, ref_ent_list);
  /* KGM -- need this?
  if ( !same_modify_engine(ref_ent_list, CUBIT_TRUE) )   
  {
    PRINT_ERROR("Performing SPLIT with surfaces containing geometry from\n"
      "different modeling engines is not allowed.\n"
      "Delete uncommon geometry on these surfaces before operation.\n\n");
    return CUBIT_FAILURE;
  }
  */

  // get the splitting curves
  DLIList<DLIList<Curve*>*> curve_lists_list;
  SplitSurfaceTool sst;
  CubitStatus err = sst.calculate_split_curves( ref_face_list, num_segs, fraction, distance,
                                                from_curve_ptr,corner_vertex_list,
                                                through_vertex_list, curve_dir_ptr,
                                                preview_flg, create_ref_edges_flg,
                                                false, curve_lists_list );


  // loop over all the curves 
  GMem gmem;
  int  num_points;
  DLIList<CubitVector*> segments;
  int k;
  for (k=0; k < curve_lists_list.size(); k++)
  {
    DLIList<Curve*> *curve_list = curve_lists_list[k];

    int i;
    for (i = 0; i < curve_list->size(); i++)
    {
      Curve* curve_ptr = curve_list->get_and_step();

      // get the curve facets
      err = curve_ptr->get_geometry_query_engine()->
                            get_graphics( curve_ptr, num_points, &gmem );

      // load the graphics points into a CubitVector for insert_curve
      --num_points;  
      int j;
      for (j = 0; j < num_points; j++)
      {
        const GPoint& p = gmem.point_list()[j];
        segments.append( new CubitVector( p.x, p.y, p.z ) );
      }
    }
  }

  // add last point to the end of the list
  const GPoint& p = gmem.point_list()[num_points];
  segments.append( new CubitVector(p.x, p.y, p.z ) );

  // Get the underlying surface (what if it is virtual does this still work?)
  //Surface *old_surf = ref_face_list[k]->get_surface_ptr(); KGM

  // now partition the surface
  DLIList<RefEdge*> new_edges;
  DLIList<RefFace*> new_faces;
  err = PartitionTool::instance()->insert_edge( ref_face_list, 
                          segments, new_faces, new_edges); 
  //DLIList<Curve*> new_curves;
  //Surface* new_surf = PartitionEngine::instance().insert_curve( old_surf, segments, new_curves );

  // insert_curve allocates memory.  Clean up.
  //while( new_curves.size() ) 
  //  delete new_curves.pop();

  // called routines allocate memory.  Clean up for them.
  for (k=0; k < curve_lists_list.size(); k++)
  {
    DLIList<Curve*> *curve_list = curve_lists_list[k];
    while( curve_list->size() ) 
      delete curve_list->pop();
  }

  return CUBIT_SUCCESS;
}
