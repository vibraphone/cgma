//--------------------------------------------------------------------------
//  Class: ImprintBoundaryTool
//  Description:  Imprints the boundaries of two surfaces.  The boundaries
//                are discritized or faceted, and the then intersected.
//                Virtual geometry is used for the actual imprinting of the
//                topology.
//  Owner: David R. White
//  Creation Date: 4/15/2002
//--------------------------------------------------------------------------
#include "ImprintBoundaryTool.hpp"
#include "ImprintPointData.hpp"
#include "ImprintLineSegment.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include "CoEdge.hpp"
#include "Loop.hpp"
#include "GeometryQueryEngine.hpp"
#include "Curve.hpp"
#include "DLIList.hpp"
#include "CubitVector.hpp"
#include "GMem.hpp"
#include "GfxDebug.hpp"
#include "CubitBox.hpp"
#include "CastTo.hpp"
#include "IntersectionTool.hpp"
#include "PartitionTool.hpp"
//#include "VirtualQueryEngine.hpp"
#include "GeometryQueryTool.hpp"
#include "MergeTool.hpp"
#include "KDDTree.hpp"
#include "AbstractTree.hpp"
#include "FacetModifyEngine.hpp"
#include "Point.hpp"
#include "FacetCurve.hpp"
#include "CubitPointData.hpp"
#include "CubitFacetEdgeData.hpp"
#include "CubitFacetEdge.hpp"
#include "CurveFacetEvalTool.hpp"

//-----------------------------------------------------
// Constructor
//-----------------------------------------------------
ImprintBoundaryTool::ImprintBoundaryTool(RefFace *ref_face_1,
                                         RefFace *ref_face_2,
                                         double tol)
{
  refFacePtr1 = ref_face_1;
  refFacePtr2 = ref_face_2;
  myTolerance = tol;
  modBound1 = CUBIT_FALSE;
  modBound2 = CUBIT_FALSE;
  allocatedPointData = new PointList;
  allocatedLineData =  new SegList;
  allocatedPointLoops =  new PointLoopList;
  allocatedLineLoops =  new SegLoopList;
  allocatedMatchData = new DLIList<ImprintMatchData*>;
  allocatedRefEdge = new DLIList<RefEdge*>;
}
//-----------------------------------------------------
// Destructor
//-----------------------------------------------------
ImprintBoundaryTool::~ImprintBoundaryTool()
{
  int ii;
  for ( ii = allocatedPointData->size(); ii > 0; ii-- )
    delete allocatedPointData->pop();
  delete allocatedPointData;
  for ( ii = allocatedLineData->size(); ii > 0; ii-- )
    delete allocatedLineData->pop();
  delete allocatedLineData;
  for ( ii = allocatedPointLoops->size(); ii > 0; ii-- )
    delete allocatedPointLoops->pop();
  delete allocatedPointLoops;
  for ( ii = allocatedLineLoops->size(); ii > 0; ii-- )
    delete allocatedLineLoops->pop();
  delete allocatedLineLoops;
  for ( ii = allocatedMatchData->size(); ii > 0; ii-- )
    delete allocatedMatchData->pop();
  delete allocatedMatchData;
  for ( ii = allocatedRefEdge->size(); ii > 0; ii-- )
  {
    RefEdge *tmp_edge_ptr = allocatedRefEdge->pop();
    if ( tmp_edge_ptr->num_ref_faces() == 0 )
    {
      RefEntity *tmp_ent = CAST_TO(tmp_edge_ptr, RefEntity);
      GeometryQueryTool::instance()->delete_RefEntity(tmp_ent);
    }
  }

}
//-----------------------------------------------------
// Public Function: imprint
// Description: Does the actual imprinting of the boundaries
//              of the two surfaces passed in. This funcntion
//              is the primary interface for the tool.
//-----------------------------------------------------
CubitStatus ImprintBoundaryTool::imprint(DLIList <RefFace*> &results,
                                         CubitBoolean merge)
{
  CubitBox ref_box1 = refFacePtr1->bounding_box();
  CubitBox ref_box2 = refFacePtr2->bounding_box();

  if ( !ref_box1.overlap(myTolerance, ref_box2 ) )
    return CUBIT_SUCCESS;
  if ( refFacePtr1->num_loops() == 0 ||
	  refFacePtr2->num_loops() == 0 )
  {
    PRINT_WARNING("Virtual imprinting currently can't imprint\n"
                  "surfaces that have no boundary curves (surfaces %d or %d\n"
                  "has no boundary curves.\n", refFacePtr1->id(), refFacePtr2->id());
    return CUBIT_SUCCESS;
  }
  if ( merge )
    PRINT_WARNING("Can't merge while virtual imprinting yet.\n");
  PRINT_DEBUG_129("Imprinting surfaces: %d and %d\n",
                  refFacePtr1->id(),
                  refFacePtr2->id() );
  int debug5 = 0;
  if ( debug5)
  {
    GfxDebug::clear();
    GfxDebug::draw_ref_face_edges(refFacePtr1, CUBIT_GREEN);
    GfxDebug::draw_ref_face_edges(refFacePtr2, CUBIT_YELLOW);
    GfxDebug::flush();
    GfxDebug::mouse_xforms();
    CubitMessage::instance()->debug_flag(129, CUBIT_TRUE);
  }
  int ii, jj;
  PointLoopList boundary_point_loops1;
  PointLoopList boundary_point_loops2;
  CubitStatus stat = get_boundary_points(refFacePtr1,
                                         boundary_point_loops1);
  if ( stat != CUBIT_SUCCESS )
    return stat;

  stat = get_boundary_points(refFacePtr2,
                             boundary_point_loops2);
  if ( stat != CUBIT_SUCCESS )
    return stat;
  boundary_point_loops1.reset();
  boundary_point_loops2.reset();
  if ( DEBUG_FLAG(129) )
  {
      //Drawing the boundary loops.
    GfxDebug::clear();
    GfxDebug::flush();
    PointList *point_list;
    for ( ii = boundary_point_loops1.size(); ii > 0; ii-- )
    {
      point_list = boundary_point_loops1.get_and_step();
      for ( jj = point_list->size(); jj > 0; jj-- )
      {
        ImprintPointData *point = point_list->get_and_step();
        //RefEntity *entity = point->owner();
        draw_point(point);
      }
    }
    if (debug5)
    {
      GfxDebug::mouse_xforms();
      GfxDebug::clear();	
      GfxDebug::draw_ref_face_edges(refFacePtr2);
      GfxDebug::flush();
    }
    for ( ii = boundary_point_loops2.size(); ii > 0; ii-- )
    {
      point_list = boundary_point_loops2.get_and_step();
      for ( jj = point_list->size(); jj > 0; jj-- )
      {
        ImprintPointData *point = point_list->get_and_step();
        //RefEntity *entity = point->owner();
        draw_point(point);
      }
    }
    if (debug5)
    {
      GfxDebug::mouse_xforms();
    }
  }
    //Now intersect the two boundary loops.  Either mark the nodes
    //that are intersecting, or insert new nodes into the boundary loops,
    //that intersect the boundaries.
  stat = imprint_boundaries( boundary_point_loops1,
                             boundary_point_loops2,
                             refFacePtr1, refFacePtr2,
                             results);
  if ( results.size() == 0 )
  {
    if ( modBound1 )
      results.append(refFacePtr1);
    if ( modBound2 )
      results.append(refFacePtr2);
  }
  if ( stat != CUBIT_SUCCESS )
  {
    PRINT_ERROR("Imprinting Surface %d and %d failed\n",
                refFacePtr1->id(), refFacePtr2->id() );
    return CUBIT_FAILURE;
  }
  return CUBIT_SUCCESS;
}
CubitStatus ImprintBoundaryTool::get_boundary_points( RefFace *ref_face,
                                                      PointLoopList &boundary_point_loops )
{
  DLIList<DLIList<CoEdge*>*> co_edge_loops;
  ref_face->co_edge_loops(co_edge_loops);
  int ii, jj;
  DLIList <CoEdge*> *co_edge_list_ptr;
  PointList *new_point_loop_ptr, tmp_point_list;
  RefEdge *ref_edge_ptr;
  CoEdge *co_edge_ptr;
  CubitStatus stat;
  CubitSense sense;
  RefEntity *temp_ref;
  RefVertex *temp_vert;
  
  for ( ii = co_edge_loops.size(); ii > 0; ii-- )
  {
    co_edge_list_ptr = co_edge_loops.get_and_step();
    new_point_loop_ptr = new PointList;
    allocatedPointLoops->append(new_point_loop_ptr);
    for ( jj = co_edge_list_ptr->size(); jj > 0; jj-- )
    {
      co_edge_ptr = co_edge_list_ptr->get_and_step();
      ref_edge_ptr = co_edge_ptr->get_ref_edge_ptr();
      tmp_point_list.clean_out();
      stat = get_curve_facets( ref_edge_ptr, tmp_point_list );
      PRINT_DEBUG_129("curve %d has %d points\n",
                      ref_edge_ptr->id(),
                      tmp_point_list.size());
      if ( stat != CUBIT_SUCCESS ) {
        while (co_edge_loops.size())  
          delete co_edge_loops.pop();
        return CUBIT_FAILURE;
      }
      tmp_point_list.reset();
        //the points are in order from start vertex to end vertex.
        //append them now according to the loop.

        //Assign the points to be part owned by the vertex, rather than the
        //curve.
      if ( ref_edge_ptr->start_vertex() !=
           ref_edge_ptr->end_vertex() )
      {
        ImprintPointData *temp_point = tmp_point_list.get();
        CubitVector v1 = temp_point->coordinates();
        CubitVector v2 = ref_edge_ptr->start_vertex()->coordinates();
        if ( !v1.within_tolerance(v2, myTolerance) )
        {
          PRINT_ERROR("Problem with surface geometry\n"
                      "Check surface %d and %d, especially curve %d\n",
                      refFacePtr1->id(), refFacePtr2->id(), ref_edge_ptr->id());
          while (co_edge_loops.size())  
            delete co_edge_loops.pop();
        }
        temp_vert = ref_edge_ptr->start_vertex();
        temp_ref = CAST_TO(temp_vert, RefEntity);
        temp_point->owner(temp_ref);
        temp_point = tmp_point_list.prev();
        v1 = temp_point->coordinates();
        v2 = ref_edge_ptr->end_vertex()->coordinates();
        if (!v1.within_tolerance(v2, myTolerance))
        {
          PRINT_ERROR("Problem with surface geometry\n"
                      "Check surface %d and %d, especially curve %d\n",
                      refFacePtr1->id(), refFacePtr2->id(), ref_edge_ptr->id());
          while (co_edge_loops.size())  
            delete co_edge_loops.pop();
          return CUBIT_FAILURE;
        }
        temp_vert = ref_edge_ptr->end_vertex();
        temp_ref = CAST_TO(temp_vert, RefEntity);
        temp_point->owner(temp_ref);
      }
      else
      {
        ImprintPointData *temp_point = tmp_point_list.get();
        CubitVector v1 = temp_point->coordinates();
        CubitVector v2 = ref_edge_ptr->start_vertex()->coordinates();
        if ( !v1.about_equal(v2) )
        {
          temp_point = tmp_point_list.prev();
          v1 = temp_point->coordinates();
          v2 = ref_edge_ptr->start_vertex()->coordinates();
          assert(v1.about_equal(v2));
        }
        temp_vert = ref_edge_ptr->end_vertex();
        temp_ref = CAST_TO( temp_vert, RefEntity );
        temp_point->owner(temp_ref);
      }
      tmp_point_list.reset();
      sense = co_edge_ptr->get_sense();
      if ( CUBIT_FORWARD != sense )
        tmp_point_list.reverse();
        //Now take off the last point as it is a duplicate with the
        //other list...
      tmp_point_list.reset();
      if ( co_edge_list_ptr->size() != 1 )
        tmp_point_list.pop();
      (*new_point_loop_ptr) += tmp_point_list;

      if ( num_coedges_on_face(ref_edge_ptr, ref_face) > 1 )
      {
        PRINT_ERROR("Surface %d has a sipe or hard line.\n"
                    "Virtual imprinting does not support that type\n"
                    "of surface.\n", ref_face->id());
          while (co_edge_loops.size())  
            delete co_edge_loops.pop();
        return CUBIT_FAILURE;
      }
    }
    CubitVector curr, prev;
    for ( jj = new_point_loop_ptr->size(); jj > 0; jj-- )
    {
      prev = new_point_loop_ptr->prev()->coordinates();
      curr = new_point_loop_ptr->get_and_step()->coordinates();
      if ( prev.about_equal(curr) )
      {
        PRINT_DEBUG_129("Points within tolerance in boundaryloop.\n");
        int debug = 0;
        if ( debug )
        {
          GfxDebug::draw_point(new_point_loop_ptr->prev(3)->coordinates(), CUBIT_RED);
          GfxDebug::draw_point(new_point_loop_ptr->prev(2)->coordinates(), CUBIT_YELLOW);
          GfxDebug::draw_point(new_point_loop_ptr->prev()->coordinates(), CUBIT_BLUE);
          GfxDebug::draw_point(new_point_loop_ptr->get()->coordinates(), CUBIT_GREEN);
        }
        new_point_loop_ptr->back();
        new_point_loop_ptr->remove();
      }
    }
    boundary_point_loops.append(new_point_loop_ptr);
  }
  CubitVector curr1, curr2;
  int kk, ll;
  PointList *check_list;
  double min_dist = CUBIT_DBL_MAX;
  PointLoopList boundary_point_loops2 = boundary_point_loops;
  for ( ii = boundary_point_loops.size(); ii > 0; ii-- )
  {
    new_point_loop_ptr = boundary_point_loops.get_and_step();
    for ( kk = boundary_point_loops.size(); kk > 0; kk-- )
    {
      check_list = boundary_point_loops2.get_and_step();
      if ( new_point_loop_ptr == check_list )
        continue;
      for ( jj = new_point_loop_ptr->size(); jj > 0; jj-- )
      {
        curr1 = new_point_loop_ptr->get_and_step()->coordinates();
        for ( ll = check_list->size(); ll > 0; ll-- )
        {
          curr2 = check_list->get_and_step()->coordinates();
          if (curr1.about_equal(curr2, myTolerance) )
          {
            double dist = (curr1-curr2).length();
            if ( dist < min_dist )
              min_dist = dist;
          }
        }
      }
    }
  }
  if ( min_dist < CUBIT_DBL_MAX && min_dist < myTolerance  )
  {
    PRINT_INFO("Two loops on Surface %d are within tolerance of each other!\n",
               ref_face->id());
			
    PRINT_INFO("Changing the tolerance to less than %f\n"
               "(the smallest distance between the two loops.)\n",
               min_dist);

    myTolerance = min_dist/4;
    PRINT_INFO("For surface %d and %d tolerance changed to %f\n",
               refFacePtr1->id(), refFacePtr2->id(), myTolerance);
      //return CUBIT_FAILURE;
  }
  
    //clean up the list memory.
  for(ii = co_edge_loops.size(); ii>0; ii-- )
    delete co_edge_loops.pop();
  co_edge_loops.clean_out();
      
  return CUBIT_SUCCESS;
}

CubitStatus ImprintBoundaryTool::get_curve_facets( RefEdge* curve, PointList &segments ) 
{
//  const double COS_ANGLE_TOL =  0.965925826289068312213715; // cos(15)
//  const double COS_ANGLE_TOL =  0.984807753012208020315654; // cos(10)
  //const double COS_ANGLE_TOL =  0.996194698091745545198705; // cos(5)
  GMem curve_graphics;
    //make this tol bigger than myTolerance, to
    //make sure the segments are larger than the tolerance.
  const double dist_tol = 2*myTolerance + .5*myTolerance;// + .05*myTolerance;
  //const double dist_tol_sqr = dist_tol*dist_tol;
  Curve* curve_ptr = curve->get_curve_ptr();
  curve_ptr->get_geometry_query_engine()->get_graphics( curve_ptr, &curve_graphics );
  
  GPoint* gp = curve_graphics.point_list();
  ImprintPointData* last = new ImprintPointData( gp[0].x, gp[0].y, gp[0].z );
  allocatedPointData->append(last);
  last->owner(dynamic_cast<RefEntity*>(curve));
  CubitVector lastv = last->coordinates();
  int num_points = curve_graphics.pointListCount;
  segments.append( last );
  int ii;
  CubitBoolean remove_second_to_end = CUBIT_FALSE;
  for ( ii = 1; ii < num_points; ii++ )
  {
    CubitVector pos(  gp[ii].x, gp[ii].y, gp[ii].z );
    CubitVector step1 = (pos - lastv);
    double len1 = step1.length();
    if( len1 < dist_tol && ii != num_points - 1) 
      continue;
    else if ( len1 < dist_tol && ii == num_points-1 )
    {
      remove_second_to_end = CUBIT_TRUE;
    }
    last = new ImprintPointData( pos );
    if ( DEBUG_FLAG(129) )
      draw_point(last);
    allocatedPointData->append(last);
    last->owner(dynamic_cast<RefEntity*>(curve));
    segments.append( last );
    lastv = last->coordinates();
  }
    // Now check if the segment list is reversed wrt the curve direction.
  segments.reset();
  if ( remove_second_to_end )
  {
    if ( segments.size() == 2 )
    {
      PRINT_DEBUG_129("Tolerance size is small for imprinting\n" 
                      "(curve %d is small)\n",
                      curve->id());
      double leng = curve->measure();
      if ( leng < myTolerance )
      {
        PRINT_ERROR("Tolerance for surfaces %d and %d is too small.\n"
                    "Curve %d is of size %f.  Try making the tolerance\n"
                    "1/4 of that value or (%f).\n", refFacePtr1->id(),
                    refFacePtr2->id(), curve->id(), leng, leng/4.0);
        return CUBIT_FAILURE;
      }
    }
    else
    {
        //Remove the second to last one.  To do
        //this efficiently (don't do remove), pop
        //the last one, then save that and
        //re-add it after poping the second one.
      ImprintPointData *temp = segments.pop();
      segments.pop();
      segments.append(temp);
    }
  }
  segments.reset();
  if( curve->start_vertex() != curve->end_vertex() )
  {
    CubitVector start_vec, end_vec;
    start_vec = curve->start_vertex()->coordinates();
    end_vec = curve->end_vertex()->coordinates();
    CubitVector start_seg = segments.get()->coordinates();
    double dist_1 = (start_seg - start_vec).length_squared();
    double dist_2 = (start_seg - end_vec).length_squared();
    if ( dist_1 > dist_2 )
      segments.reverse();
      //Now make sure that the start and end vertices are "right-on" with
      //the facet representations...
    segments.reset();
    start_seg = segments.get()->coordinates();
    CubitVector end_seg = segments.prev()->coordinates();
      //make sure the start and end positions match up "exactly"
      //with the vertex coordinates of the curves.  The facets sometimes
      //aren't really "right-on"
    segments.get()->set(start_vec);
    segments.prev()->set(end_seg);
  }
  else
  {
    double u1, u2;
    u1 = curve->u_from_position( (segments.next(1)->coordinates()) );
    u2 = curve->u_from_position( (segments.next(2)->coordinates()) );    
    if( (u2 < u1) && (curve->start_param() <= curve->end_param()) )
      segments.reverse();
  }
    //clean up the periodic curve case (last seg may be too small.)
  if ( curve->start_vertex() == curve->end_vertex() )
  {
    segments.reset();
    CubitVector start_v = segments.get()->coordinates();
    CubitVector last_v = segments.prev()->coordinates();
    double dist = (start_v - last_v).length();
    if ( dist < dist_tol )
    {
        //remove the last one.
      segments.pop();
    }
      //now make sure the one vertex matches up "exactly" with the
      //first segement.
    segments.reset();
    CubitVector start_vec = curve->start_vertex()->coordinates();
      //replace start_seg with start_vec.
    segments.get()->set(start_vec);
  }
    
    //Make sure we don't have duplicate points.
  int jj;
  CubitVector curr, prev;
  for ( jj = segments.size(); jj > 0; jj-- )
  {
    prev = segments.prev()->coordinates();
    curr = segments.get_and_step()->coordinates();
    if ( prev.about_equal(curr) )
    {
      PRINT_DEBUG_129("Points on curve %d within tolerance...\n", curve->id());
      segments.back();
      segments.remove();
    }
  }
  if ( segments.size() < 2 )
  {
    PRINT_ERROR("Tolerance size is too small for imprinting\n"
                "(curve %d is small)\n",
                curve->id());
    return CUBIT_FAILURE;
  }
  return CUBIT_SUCCESS;
}
//CubitStatus ImprintBoundaryTool::find_intersections(PointLoopList &boundary_line_loops_1,
//                                                    PointLoopList &boundary_line_loops_2,
//--------------------------------------------------------
// Public Function: imprint_boundaries
// Description: Given the two lists of boundary loops, where
//              the nodes in the loops are assumed to be from either
//              a refedge or a refvertex, intersect them with the other
//              loops.
//--------------------------------------------------------
CubitStatus ImprintBoundaryTool::imprint_boundaries(PointLoopList &boundary_loops_1,
                                                    PointLoopList &boundary_loops_2,
                                                    RefFace *ref_face_1,
                                                    RefFace *ref_face_2,
                                                    DLIList <RefFace*> &results)
{
  int ii, jj;
  SegLoopList boundary_line_loops_1, boundary_line_loops_2;
  
    //Now convert the point lists to line segments.
    //REMEMBER, delete all the data from boundary_line_loops, the lists and the data...
  CubitStatus st = convert_to_lines( boundary_loops_1, boundary_line_loops_1, ref_face_1, CUBIT_TRUE);
  if ( st != CUBIT_SUCCESS )
  {
    return CUBIT_FAILURE;
  }
  st = convert_to_lines( boundary_loops_2, boundary_line_loops_2, ref_face_2, CUBIT_FALSE);
  if ( st != CUBIT_SUCCESS )
  {
    return CUBIT_FAILURE;
  }
  if ( boundary_line_loops_1.size() == 0 || boundary_line_loops_2.size() == 0 )
    return CUBIT_FAILURE;    
  //Do the actual intersecting and imprinting.  Basically splits segments and
    //matches the points on the two boundaries.  Also classifies the connections.
  CubitStatus stat = imprint_segments( boundary_line_loops_1,
                                       boundary_line_loops_2,
                                       boundary_loops_1,
                                       boundary_loops_2);

  if (stat != CUBIT_SUCCESS )
    return CUBIT_FAILURE;
    //Note that the boundary_line_loops are outof date at this point and shouldn't be used. 
    //These lists are also stored in the allocatedLineLoops list and the memory will be
    //cleaned up in the destructor.  Just clean them out to avoid useing the data.
  boundary_line_loops_1.clean_out();
  boundary_line_loops_2.clean_out();
  
    //Mark each node as to which loop its in and its position
    //for efficient searching.
  PointList *point_loop;
  for ( ii = 0; ii < boundary_loops_1.size(); ii++ )
  {
    point_loop = boundary_loops_1.get_and_step();
    int loop_size = point_loop->size();
    for ( jj = 0; jj < loop_size; jj++ )
      point_loop->get_and_step()->set_loop_pos(ii,jj, loop_size);
  }
  for ( ii = 0; ii < boundary_loops_2.size(); ii++ )
  {
    point_loop = boundary_loops_2.get_and_step();
    int loop_size = point_loop->size();
    for ( jj = 0; jj < loop_size; jj++ )
      point_loop->get_and_step()->set_loop_pos(ii,jj, loop_size);
  }

    //Now go through and find the curves that need to be created, and vertices
    //that need to be imprinted for each surface.
  PointLoopList part_segs_1, part_segs_2;
  PointList partition_points_1, partition_points_2;
  
  stat = find_graph_for_surf( boundary_loops_1,
                              boundary_loops_2,
                              ref_face_1,
                              part_segs_1,
                              partition_points_1,
                              CUBIT_TRUE);
  if ( stat != CUBIT_SUCCESS )
    return CUBIT_FAILURE;
  stat = find_graph_for_surf( boundary_loops_2,
                              boundary_loops_1,
                              ref_face_2,
                              part_segs_2,
                              partition_points_2,
                              CUBIT_FALSE);
  if ( stat != CUBIT_SUCCESS )
    return CUBIT_FAILURE;
  
    //Okay now go through and partition the boundaries.
  CubitBoolean mod_bound1, mod_bound2;
  stat = imprint_boundary_vertices( partition_points_1,
                                    mod_bound1);
  if ( stat != CUBIT_SUCCESS )
    return CUBIT_FAILURE;
  stat = imprint_boundary_vertices( partition_points_2,
                                    mod_bound2);
  if ( stat != CUBIT_SUCCESS )
    return CUBIT_FAILURE;
    //Now go through and paritition the surfaces.
  DLIList<RefFace*> results_1, results_2;
  stat = imprint_surface(ref_face_1, part_segs_1, results_1);
  if ( stat != CUBIT_SUCCESS )
    return CUBIT_FAILURE;
  stat = imprint_surface(ref_face_2, part_segs_2, results_2);
  if ( stat != CUBIT_SUCCESS )
    return CUBIT_FAILURE;
  if ( mod_bound1 || results_1.size() )
    modBound1 = CUBIT_TRUE;
  if ( mod_bound2 || results_2.size() )
    modBound2 = CUBIT_TRUE;
  results += results_1;
  results += results_2;
  return CUBIT_SUCCESS;
}

CubitStatus ImprintBoundaryTool::intersect_segments( ImprintLineSegment *seg_1,
                                                     ImprintLineSegment *seg_2,
                                                     IntersectResult &int_result,
                                                     ImprintLineSegment **new_segments)
{
  int ii;
    //initialize the results first.
  for ( ii = 0; ii < 4; ii++)
    new_segments[ii] = NULL;
  ImprintPointData* imp_point_0 = seg_1->get_start();
  ImprintPointData* imp_point_1 = seg_1->get_end();
  ImprintPointData* imp_point_2 = seg_2->get_start();
  ImprintPointData* imp_point_3 = seg_2->get_end();
  int debug4 = 0;
  if ( debug4 )
  {
    draw_point(imp_point_0);
    draw_point(imp_point_1);
    draw_point(imp_point_2);
    draw_point(imp_point_3);
  }
  int debug = 0;
    //First go through and test the end points.
    //This will test for cases L_INTERSECT, and SEGS_EQUAL, and OVERLAP_JOIN.
  MatchType type_0, type_1, type_2, type_3;
  CubitStatus stat1 = match_points(seg_1, seg_2,type_0,
                                   type_1, type_2, type_3);
    //First handle the basic stuff;
  if ( type_0 == MATCH_0_2 && type_1 == MATCH_1_3 &&
       type_2 == MATCH_0_2 && type_3 == MATCH_1_3 )
  {
    imp_point_0->set_matching_point(imp_point_2);
    imp_point_2->set_matching_point(imp_point_0);
    imp_point_1->set_matching_point(imp_point_3);
    imp_point_3->set_matching_point(imp_point_1);
    
      //These segments are equal.
    set_type_for_equal(imp_point_0, imp_point_2,
                       imp_point_1, imp_point_3);
    int_result = SEGS_EQUAL_0_2;
    PRINT_DEBUG_129("Found SEGS_EQUAL_0_2\n");
    return CUBIT_SUCCESS;
  }
  else if ( type_0 == MATCH_0_3 && type_1 == MATCH_1_2 &&
            type_3 == MATCH_0_3 && type_2 == MATCH_1_2 )
  {
    imp_point_0->set_matching_point(imp_point_3);
    imp_point_3->set_matching_point(imp_point_0);
    imp_point_1->set_matching_point(imp_point_2);
    imp_point_2->set_matching_point(imp_point_1);
    
      //These segments are equal.
    set_type_for_equal(imp_point_0, imp_point_3,
                       imp_point_1, imp_point_2);
    int_result = SEGS_EQUAL_0_3;
    PRINT_DEBUG_129("Found SEGS_EQUAL_0_3\n");
    return CUBIT_SUCCESS;
  }
  else if ( type_0 == MATCH_0_2 && type_2 == MATCH_0_2 &&
            (type_1 == NO_MATCH || type_1 == MATCH_1_6 ||
             type_1 == MATCH_1_7) &&
            (type_3 == NO_MATCH || type_3 == MATCH_3_4 ||
             type_3 == MATCH_3_5) )
  {
    imp_point_0->set_matching_point(imp_point_2);
    imp_point_2->set_matching_point(imp_point_0);
    return case_0_2_equal(seg_1, seg_2,
                          imp_point_0, imp_point_1,
                          imp_point_2, imp_point_3,
                          type_1, type_3,
                          int_result, new_segments);
  }
  else if ( type_0 == MATCH_0_3 && type_3 == MATCH_0_3 &&
            (type_1 == NO_MATCH || type_1 == MATCH_1_6 ||
             type_1 == MATCH_1_7) &&
            (type_2 == NO_MATCH || type_2 == MATCH_2_4 ||
             type_2 == MATCH_2_5) )
  {
    
    imp_point_0->set_matching_point(imp_point_3);
    imp_point_3->set_matching_point(imp_point_0);
    return case_0_3_equal(seg_1, seg_2,
                          imp_point_0, imp_point_1,
                          imp_point_2, imp_point_3,
                          type_1, type_2,
                          int_result, new_segments);
  }
  else if ( type_1 == MATCH_1_2 && type_2 == MATCH_1_2 &&
            (type_0 == NO_MATCH || type_0 == MATCH_0_6 ||
             type_0 == MATCH_0_7) &&
            (type_3 == NO_MATCH || type_3 == MATCH_3_4 ||
             type_3 == MATCH_3_5) )
  {
    
    imp_point_1->set_matching_point(imp_point_2);
    imp_point_2->set_matching_point(imp_point_1);
    return case_1_2_equal(seg_1, seg_2,
                          imp_point_0, imp_point_1,
                          imp_point_2, imp_point_3,
                          type_0, type_3,
                          int_result, new_segments);
  }
  else if ( type_1 == MATCH_1_3 && type_3 == MATCH_1_3 &&
            (type_0 == NO_MATCH || type_0 == MATCH_0_6 ||
             type_0 == MATCH_0_7) &&
            (type_2 == NO_MATCH || type_2 == MATCH_2_4 ||
             type_2 == MATCH_2_5) )
  {
    
    imp_point_1->set_matching_point(imp_point_3);
    imp_point_3->set_matching_point(imp_point_1);
    return case_1_3_equal(seg_1, seg_2,
                          imp_point_0, imp_point_1,
                          imp_point_2, imp_point_3,
                          type_0, type_2,
                          int_result, new_segments);
  }
  else if ( type_0 != NO_MATCH || type_1 != NO_MATCH ||
            type_2 != NO_MATCH || type_3 != NO_MATCH )
  {
    if ( type_0 == MATCH_0_2 || type_0 == MATCH_0_3 ||
         type_1 == MATCH_1_2 || type_1 == MATCH_1_3 ||
         type_2 == MATCH_1_2 || type_2 == MATCH_0_2 ||
         type_3 == MATCH_0_3 || type_3 == MATCH_1_3 )
    {
      if (debug )
      {
        draw_point(imp_point_0);
        draw_point(imp_point_1);
        draw_point(imp_point_2);
        draw_point(imp_point_3);
        draw_point(seg_1->get_prev()->get_start());
        draw_point(seg_1->get_next()->get_end());
        draw_point(seg_2->get_prev()->get_start());
        draw_point(seg_2->get_next()->get_end());
        GfxDebug::mouse_xforms();
      }
      if ( type_1 == MATCH_1_2 && type_2 == MATCH_1_2 &&
           imp_point_1->get_matching_point() == imp_point_2 )
      {
        if ( type_0 == MATCH_0_3 && type_3 == MATCH_3_5 )
        {
          int_result = NO_INTERSECT;
          return CUBIT_SUCCESS;
        }
        else if ( type_3 == MATCH_0_3 && type_0 == MATCH_0_7 )
        {
          int_result = NO_INTERSECT;
          return CUBIT_SUCCESS;
        }
        else
        {
          if (( type_0 == MATCH_0_2 ||  type_3 == MATCH_1_3 ) ||
              ( type_0 == MATCH_0_3 && type_3 == MATCH_3_4 )  ||
              ( type_0 == MATCH_0_7 && type_3 == MATCH_0_3 )  ||
              ( type_0 == MATCH_0_6 && type_3 == MATCH_0_3 )  ||
              ( type_0 == MATCH_0_3 && type_3 == MATCH_3_5 ) )
          {
              //This means that 0 and 2 are within tolerance
              //but 1 and 2 are closer.  For now just ignore this...
              //or...
              //This means that 1 and 3 are within tolerance
              //but 1 and 2 are closer.  For now just ignore this...
            imp_point_1->set_matching_point(imp_point_2);
            imp_point_2->set_matching_point(imp_point_1);
            return case_1_2_equal(seg_1, seg_2,
                                  imp_point_0, imp_point_1,
                                  imp_point_2, imp_point_3,
                                  type_0, type_3,
                                  int_result, new_segments);

          }
          else
          {
            PRINT_ERROR("Problems with tolerance. (Need special match code.)\n");
            return CUBIT_FAILURE;
          }
        }
      }
      else if ( type_0 == MATCH_0_3 && type_3 == MATCH_0_3 &&
                imp_point_0->get_matching_point() == imp_point_3 )
      {
        if ( type_1 == MATCH_1_2 && type_2 == MATCH_2_4 )
        {
          int_result = NO_INTERSECT;
          return CUBIT_SUCCESS;
        }
        else if ( type_2 == MATCH_1_2 && type_1 == MATCH_1_7 )
        {
          int_result = NO_INTERSECT;
          return CUBIT_SUCCESS;
        }
        else
        {
          if (( type_1 == MATCH_1_3 ||  type_2 == MATCH_0_2 ) ||
              ( type_1 == MATCH_1_2 && type_2 == MATCH_2_5 ) ||
              ( type_1 == MATCH_1_2 && type_2 == MATCH_2_4 ) ||
              ( type_1 == MATCH_1_6 && type_2 == MATCH_1_2 ) ||
              ( type_1 == MATCH_1_7 && type_2 == MATCH_1_2 ) )

          {
              //This means that 1 and 3 are within tolerance
              //but 0 and 3 are closer.  For now just ignore this...
              //or...
              //This means that 0 and 2 are within tolerance
              //but 0 and 3 are closer.  For now just ignore this...
            imp_point_0->set_matching_point(imp_point_3);
            imp_point_3->set_matching_point(imp_point_0);
            return case_0_3_equal(seg_1, seg_2,
                                  imp_point_0, imp_point_1,
                                  imp_point_2, imp_point_3,
                                  type_1, type_2,
                                  int_result, new_segments);
            
          }
          else
          {
            PRINT_ERROR("Problems with tolerance. (Need special match code.)\n");
            return CUBIT_FAILURE;
          }
        }
      }
      else if ( type_1 == MATCH_1_3 && type_3 == MATCH_1_3 &&
                imp_point_1->get_matching_point() == imp_point_3 )
      {
        if ( type_0 == MATCH_0_2 && type_2 == MATCH_2_5 )
        {
          int_result = NO_INTERSECT;
          return CUBIT_SUCCESS;
        }
        else if ( type_2 == MATCH_0_2 && type_0 == MATCH_0_7 )
        {
          int_result = NO_INTERSECT;
          return CUBIT_SUCCESS;
        }
        else
        {
          if (( type_0 == MATCH_0_3 ||  type_2 == MATCH_1_2 ) ||
              ( type_0 == MATCH_0_2 &&  type_2 == MATCH_2_4 ) ||
              ( type_0 == MATCH_0_2 &&  type_2 == MATCH_2_5 ) ||
              ( type_0 == MATCH_0_6 && type_2 == MATCH_0_2 )  ||
              ( type_0 == MATCH_0_7 && type_2 == MATCH_0_2 ) )

          {
              //This means that 0 and 3 are within tolerance
              //but 1 and 3 are closer.  For now just ignore this...
              //or...
              //This means that 1 and 2 are within tolerance
              //but 1 and 3 are closer.  For now just ignore this...
            imp_point_1->set_matching_point(imp_point_3);
            imp_point_3->set_matching_point(imp_point_1);
            return case_1_3_equal(seg_1, seg_2,
                                  imp_point_0, imp_point_1,
                                  imp_point_2, imp_point_3,
                                  type_0, type_2,
                                  int_result, new_segments);
            
          }
          else
          {
            PRINT_ERROR("Problems with tolerance. (Need special match code.)\n");
            return CUBIT_FAILURE;
          }
        }
      }
      else if ( type_0 == MATCH_0_2 && type_2 == MATCH_0_2 &&
                imp_point_0->get_matching_point() == imp_point_2 )
      {
        if ( type_1 == MATCH_1_3 && type_3 == MATCH_3_4 )
        {
          int_result = NO_INTERSECT;
          return CUBIT_SUCCESS;
        }
        else if ( type_3 == MATCH_1_3 && type_1 == MATCH_1_6 )
        {
          int_result = NO_INTERSECT;
          return CUBIT_SUCCESS;
        }
        else
        {
          if (( type_1 == MATCH_1_2 ||  type_3 == MATCH_0_3 ) ||
              ( type_1 == MATCH_1_3 && type_3 == MATCH_3_5) ||
              ( type_1 == MATCH_1_3 && type_3 == MATCH_3_4) ||
              ( type_1 == MATCH_1_6 && type_3 == MATCH_1_3) ||
              ( type_1 == MATCH_1_7 && type_3 == MATCH_1_3) )

          {
              //This means that 0 and 3 are within tolerance
              //but 0 and 2 are closer.  For now just ignore this...
              //or...
              //This means that 1 and 2 are within tolerance
              //but 0 and 2 are closer.  For now just ignore this...
            imp_point_0->set_matching_point(imp_point_2);
            imp_point_2->set_matching_point(imp_point_0);
            return case_0_2_equal(seg_1, seg_2,
                                  imp_point_0, imp_point_1,
                                  imp_point_2, imp_point_3,
                                  type_1, type_3,
                                  int_result, new_segments);
            
          }
          else
          {
            PRINT_ERROR("Problems with tolerance. (Need special match code.)\n");
            return CUBIT_FAILURE;
          }
        }
      }
      else if ( type_0 == MATCH_0_2 || type_0 == MATCH_0_3 && 
                (type_2 != MATCH_0_2 && type_3 != MATCH_0_3 ) )
      {
        if ( type_2 != MATCH_1_2 && type_3 != MATCH_1_3 &&
             type_1 != MATCH_1_2 && type_1 != MATCH_1_3 )
        {
          int_result = NO_INTERSECT;
          return CUBIT_SUCCESS;
        }
        else
        {
          PRINT_ERROR("Problems with tolerance. (Need special match code.)\n");
          return CUBIT_FAILURE;
        }
      }
      else if ( type_1 == MATCH_1_2 || type_1 == MATCH_1_3 && 
                (type_2 != MATCH_1_2 && type_3 != MATCH_1_3 ) )
      {
        if ( type_2 != MATCH_0_2 && type_3 != MATCH_0_3 &&
             type_0 != MATCH_0_2 && type_1 != MATCH_0_3 )
        {
          int_result = NO_INTERSECT;
          return CUBIT_SUCCESS;
        }
        else
        {
          PRINT_ERROR("Problems with tolerance. (Need special match code.)\n");
          return CUBIT_FAILURE;
        }
      }
      else if ( type_2 == MATCH_0_2 || type_2 == MATCH_1_2 && 
                (type_0 != MATCH_0_2 && type_1 != MATCH_1_2 ) )
      {
        if ( type_0 != MATCH_0_2 && type_1 != MATCH_1_2 &&
             type_3 != MATCH_0_3 && type_3 != MATCH_1_3 )
        {
          int_result = NO_INTERSECT;
          return CUBIT_SUCCESS;
        }
        else
        {
          PRINT_ERROR("Problems with tolerance. (Need special match code.)\n");
          return CUBIT_FAILURE;
        }
      }
      else if ( type_3 == MATCH_0_3 || type_3 == MATCH_1_3 && 
                (type_0 != MATCH_0_3 && type_1 != MATCH_1_3 ) )
      {
        if ( type_0 != MATCH_0_3 && type_1 != MATCH_1_3 &&
             type_2 != MATCH_0_2 && type_2 != MATCH_1_2 )
        {
          int_result = NO_INTERSECT;
          return CUBIT_SUCCESS;
        }
        else
        {
          PRINT_ERROR("Problems with tolerance. (Need special match code.)\n");
          return CUBIT_FAILURE;
        }
      }
      else
      {
        PRINT_ERROR("Problems with tolerance. (Need special match code.)\n");
        return CUBIT_FAILURE;
      }
    }
  }
  
    //At this point we know that none of the end points of the
    //line segements are within tolerance.  The remaining possibilities
    //Are: T_INTERSECT, OVERLAP_PART, OVERLAP_ALL, or CROSS_INTERSECT.
  
  CubitBoolean point_0_on_seg_2 = CUBIT_FALSE;
  CubitBoolean point_1_on_seg_2 = CUBIT_FALSE;
  CubitBoolean point_2_on_seg_1 = CUBIT_FALSE;
  CubitBoolean point_3_on_seg_1 = CUBIT_FALSE;
    //Test to get the correct values for these flags.  From the above
    //flags we can figure out everything except for CROSS_INTERSECT.
  stat1 = determine_on_boundary( seg_1, seg_2, type_0, type_1,
                                 type_2, type_3 );
  if ( stat1 != CUBIT_SUCCESS )
    return stat1;
      
  if ( type_0 == ON_SEG_2 )
    point_0_on_seg_2 = CUBIT_TRUE;
  if ( type_1 == ON_SEG_2 )
    point_1_on_seg_2 = CUBIT_TRUE;
  if ( type_2 == ON_SEG_1 )
    point_2_on_seg_1 = CUBIT_TRUE;
  if ( type_3 == ON_SEG_1 )
    point_3_on_seg_1 = CUBIT_TRUE;
    //Okay, test for T_INTERSECT.
  if ( point_0_on_seg_2 && !point_1_on_seg_2 &&
       !point_2_on_seg_1 && !point_3_on_seg_1 )
  {
    PRINT_DEBUG_129("Found T_INTERSECT_0\n");
    int_result = T_INTERSECT_0;
      //split seg_2
    ImprintLineSegment *new_seg_2_0, *new_seg_0_3;
      //copy imp_point_0 to use it in list 2.
    CubitVector closest_point;
    closest_point_seg( imp_point_0, seg_2, closest_point);
    ImprintPointData *new_point_0 = new ImprintPointData(imp_point_0, closest_point);
    new_point_0->owner(seg_2->owner());
    allocatedPointData->append(new_point_0);
      //set the point type to be so we create a vertex here.
    new_point_0->set_point_type(CREATE_NEW_VERTEX);
    if ( imp_point_0->is_owner_vertex() &&
         imp_point_0->get_point_type() != CREATE_NEW_VERTEX )
      imp_point_0->set_point_type(VERTEX_ON_BOTH_BOUNDARIES);
    else
      imp_point_0->set_point_type(CREATE_NEW_VERTEX);
      
    new_seg_2_0 = new ImprintLineSegment(imp_point_2, new_point_0, seg_2->owner());
    new_seg_0_3 = new ImprintLineSegment(new_point_0, imp_point_3, seg_2->owner());
    allocatedLineData->append(new_seg_2_0);
    allocatedLineData->append(new_seg_0_3);
    new_segments[2] = new_seg_2_0;
    new_segments[3] = new_seg_0_3;
      //update the linked list.
    update_list(new_segments, seg_1, seg_2, CUBIT_FALSE);
    seg_2->set_inactive(CUBIT_TRUE);
    return CUBIT_SUCCESS;
  }
  else if ( point_1_on_seg_2 && !point_0_on_seg_2 &&
            !point_2_on_seg_1 && !point_3_on_seg_1 )
  {
      //split seg_2
    PRINT_DEBUG_129("Found T_INTERSECT_1\n");
    int_result = T_INTERSECT_1;
    ImprintLineSegment *new_seg_2_1, *new_seg_1_3;
      //copy imp_point_1 to use it in list 2.
    CubitVector closest_point;
    closest_point_seg( imp_point_1, seg_2, closest_point);
    ImprintPointData *new_point_1 = new ImprintPointData(imp_point_1, closest_point);
    new_point_1->owner(seg_2->owner());
    allocatedPointData->append(new_point_1);
      //set the point type to be so we create a vertex here.
    new_point_1->set_point_type(CREATE_NEW_VERTEX);
    if ( imp_point_1->is_owner_vertex() &&
         imp_point_1->get_point_type() != CREATE_NEW_VERTEX )
      imp_point_1->set_point_type(VERTEX_ON_BOTH_BOUNDARIES);
    else
      imp_point_1->set_point_type(CREATE_NEW_VERTEX);
    new_seg_2_1 = new ImprintLineSegment(imp_point_2,
                                         new_point_1, seg_2->owner());
    new_seg_1_3 = new ImprintLineSegment(new_point_1,
                                         imp_point_3, seg_2->owner());
    allocatedLineData->append(new_seg_2_1);
    allocatedLineData->append(new_seg_1_3);
    new_segments[2] = new_seg_2_1;
    new_segments[3] = new_seg_1_3;
      //update the linked list
    update_list(new_segments, seg_1, seg_2, CUBIT_FALSE);
    seg_2->set_inactive(CUBIT_TRUE);
    return CUBIT_SUCCESS;
  }
  else if ( point_2_on_seg_1 && !point_0_on_seg_2 &&
            !point_1_on_seg_2 && !point_3_on_seg_1 )
  {
    PRINT_DEBUG_129("Found T_INTERSECT_2\n");
    int_result = T_INTERSECT_2;
      //split seg_1
    ImprintLineSegment *new_seg_0_2, *new_seg_2_1;
      //copy imp_point_2 to use it in list 1.
    CubitVector closest_point;
    closest_point_seg( imp_point_2, seg_1, closest_point);
    ImprintPointData *new_point_2 = new ImprintPointData(imp_point_2, closest_point);
    new_point_2->owner(seg_1->owner());
    allocatedPointData->append(new_point_2);
      //set the point type to be so we create a vertex here.
    new_point_2->set_point_type(CREATE_NEW_VERTEX);
    if ( imp_point_2->is_owner_vertex() &&
         imp_point_2->get_point_type() != CREATE_NEW_VERTEX )
      imp_point_2->set_point_type(VERTEX_ON_BOTH_BOUNDARIES);
    else
      imp_point_2->set_point_type(CREATE_NEW_VERTEX);
    new_seg_0_2 = new ImprintLineSegment(imp_point_0, new_point_2,
                                         seg_1->owner());
    new_seg_2_1 = new ImprintLineSegment(new_point_2, imp_point_1,
                                         seg_1->owner());
    allocatedLineData->append(new_seg_0_2);
    allocatedLineData->append(new_seg_2_1);
    new_segments[0] = new_seg_0_2;
    new_segments[1] = new_seg_2_1;
      //update the linked list
    update_list(new_segments, seg_1, seg_2, CUBIT_TRUE);
    seg_1->set_inactive(CUBIT_TRUE);
    return CUBIT_SUCCESS;
  }
  else if ( point_3_on_seg_1 && !point_0_on_seg_2 &&
            !point_1_on_seg_2 && !point_2_on_seg_1 )
  {
    PRINT_DEBUG_129("Found T_INTERSECT_3\n");
    int_result = T_INTERSECT_3;
      //split seg_1
      //copy imp_point_3 to use it in list 1.
    CubitVector closest_point;
    closest_point_seg( imp_point_3, seg_1, closest_point);
    ImprintPointData *new_point_3 = new ImprintPointData(imp_point_3, closest_point);
    new_point_3->owner(seg_1->owner());
    allocatedPointData->append(new_point_3);
      //set the point type to be so we create a vertex here.
    new_point_3->set_point_type(CREATE_NEW_VERTEX);
    if ( imp_point_3->is_owner_vertex() &&
         imp_point_3->get_point_type() != CREATE_NEW_VERTEX )
      imp_point_3->set_point_type(VERTEX_ON_BOTH_BOUNDARIES);
    else
      imp_point_3->set_point_type(CREATE_NEW_VERTEX);

    ImprintLineSegment *new_seg_0_3, *new_seg_3_1;
    new_seg_0_3 = new ImprintLineSegment(imp_point_0, new_point_3,
                                         seg_1->owner());
    new_seg_3_1 = new ImprintLineSegment(new_point_3, imp_point_1,
                                         seg_1->owner());
    allocatedLineData->append(new_seg_0_3);
    allocatedLineData->append(new_seg_3_1);
    new_segments[0] = new_seg_0_3;
    new_segments[1] = new_seg_3_1;
      //update the linked list
    update_list(new_segments, seg_1, seg_2, CUBIT_TRUE);
    seg_1->set_inactive(CUBIT_TRUE);
    return CUBIT_SUCCESS;
  }
    //Okay, we are done with T intersections, it now must be one
    //of the remaining overlaps (part or all.) or the cross intersection.

    //First test for the part overlaps.
  if ( point_1_on_seg_2 && point_2_on_seg_1 &&
       !point_0_on_seg_2 && !point_3_on_seg_1 )
  {
    PRINT_DEBUG_129("Found OVEVERLAP_PART_0_3\n");
    int_result = OVERLAP_PART_0_3;
    ImprintLineSegment *new_seg_0_2, *new_seg_2_1;
    ImprintLineSegment *new_seg_2_1_v2, *new_seg_1_3;
      //Copy point 2 to use on seg_1
    CubitVector closest_point;
    closest_point_seg( imp_point_2, seg_1, closest_point);
    ImprintPointData *new_point_2 = new ImprintPointData(imp_point_2, closest_point);
    new_point_2->owner(seg_1->owner());
    allocatedPointData->append(new_point_2);
    if ( imp_point_2->is_owner_vertex() )
    {
      new_point_2->set_point_type(CREATE_NEW_VERTEX);
      if ( imp_point_2->get_point_type() != CREATE_NEW_VERTEX )
        imp_point_2->set_point_type(VERTEX_ON_BOTH_BOUNDARIES);
    }
    else
    {
      new_point_2->set_point_type(ON_BOTH_BOUNDARIES);
      imp_point_2->set_point_type(ON_BOTH_BOUNDARIES);
    }
      //split both segs 1 and 2.
    new_seg_0_2 = new ImprintLineSegment(imp_point_0, new_point_2,
                                         seg_1->owner());
    new_seg_2_1 = new ImprintLineSegment(new_point_2, imp_point_1,
                                         seg_1->owner());
    allocatedLineData->append(new_seg_0_2);
    allocatedLineData->append(new_seg_2_1);
      //Copy point 1 to use on seg_2
    closest_point_seg( imp_point_1, seg_2, closest_point);
    ImprintPointData *new_point_1 = new ImprintPointData(imp_point_1, closest_point);
    new_point_1->owner(seg_2->owner());
    allocatedPointData->append(new_point_1);
    if ( imp_point_1->is_owner_vertex() )
    {
      new_point_1->set_point_type(CREATE_NEW_VERTEX);
      if ( imp_point_1->get_point_type() != CREATE_NEW_VERTEX )
        imp_point_1->set_point_type(VERTEX_ON_BOTH_BOUNDARIES);
    }
    else
      new_point_1->set_point_type(ON_BOTH_BOUNDARIES);
    
    new_seg_2_1_v2 = new ImprintLineSegment(imp_point_2, new_point_1,
                                            seg_2->owner());
    new_seg_1_3 = new ImprintLineSegment(new_point_1, imp_point_3,
                                         seg_2->owner());
    allocatedLineData->append(new_seg_2_1_v2);
    allocatedLineData->append(new_seg_1_3);
    new_segments[0] = new_seg_0_2;
    new_segments[1] = new_seg_2_1;
    new_segments[2] = new_seg_2_1_v2;
    new_segments[3] = new_seg_1_3;
      //now update both linked lists.
    update_list(new_segments, seg_1, seg_2, CUBIT_TRUE);
    update_list(new_segments, seg_1, seg_2, CUBIT_FALSE);
    seg_1->set_inactive(CUBIT_TRUE);
    seg_2->set_inactive(CUBIT_TRUE);
    return CUBIT_SUCCESS;
  }
  else if ( point_1_on_seg_2 && !point_2_on_seg_1 &&
            !point_0_on_seg_2 && point_3_on_seg_1 )
  {
    PRINT_DEBUG_129("Found OVEVERLAP_PART_0_2\n");
    int_result = OVERLAP_PART_0_2;
    ImprintLineSegment *new_seg_0_3, *new_seg_3_1;
    ImprintLineSegment *new_seg_2_1, *new_seg_1_3;
      //Copy point 2 to use on seg_1
    CubitVector closest_point;
    closest_point_seg( imp_point_3, seg_1, closest_point);
    ImprintPointData *new_point_3 = new ImprintPointData(imp_point_3, closest_point);
    new_point_3->owner(seg_1->owner());
    allocatedPointData->append(new_point_3);
    int debug1 = 0;
    if ( debug1 )
    {
      GfxDebug::clear();
      draw_seg(seg_1);
      draw_seg(seg_2);
      draw_point(imp_point_0);
      draw_point(imp_point_1);
      draw_point(imp_point_2);
      draw_point(imp_point_3);
      GfxDebug::mouse_xforms();
    }
    if ( imp_point_3->is_owner_vertex() )
    {
      new_point_3->set_point_type(CREATE_NEW_VERTEX);
      if ( imp_point_3->get_point_type() != CREATE_NEW_VERTEX )
        imp_point_3->set_point_type(VERTEX_ON_BOTH_BOUNDARIES);
    }
    else
    {
      new_point_3->set_point_type(ON_BOTH_BOUNDARIES);
      imp_point_3->set_point_type(ON_BOTH_BOUNDARIES);
    }
      //split both segs 1 and 2.
    new_seg_0_3 = new ImprintLineSegment(imp_point_0, new_point_3,
                                         seg_1->owner());
    new_seg_3_1 = new ImprintLineSegment(new_point_3, imp_point_1,
                                         seg_1->owner());
    allocatedLineData->append(new_seg_0_3);
    allocatedLineData->append(new_seg_3_1);
      //Copy point 1 to use on seg_2
    closest_point_seg( imp_point_1, seg_2, closest_point);
    ImprintPointData *new_point_1 = new ImprintPointData(imp_point_1, closest_point);
    new_point_1->owner(seg_2->owner());
    allocatedPointData->append(new_point_1);
    if ( imp_point_1->is_owner_vertex() )
    {
      new_point_1->set_point_type(CREATE_NEW_VERTEX);
      if ( imp_point_1->get_point_type() != CREATE_NEW_VERTEX )
        imp_point_1->set_point_type(VERTEX_ON_BOTH_BOUNDARIES);
    }
    else
      new_point_1->set_point_type(ON_BOTH_BOUNDARIES);
    
    new_seg_2_1 = new ImprintLineSegment(imp_point_2, new_point_1,
                                         seg_2->owner());
    new_seg_1_3 = new ImprintLineSegment(new_point_1, imp_point_3,
                                         seg_2->owner());
    allocatedLineData->append(new_seg_2_1);
    allocatedLineData->append(new_seg_1_3);
    new_segments[0] = new_seg_0_3;
    new_segments[1] = new_seg_3_1;
    new_segments[2] = new_seg_2_1;
    new_segments[3] = new_seg_1_3;
      //now update both linked lists.
    update_list(new_segments, seg_1, seg_2, CUBIT_TRUE);
    update_list(new_segments, seg_1, seg_2, CUBIT_FALSE);
    seg_1->set_inactive(CUBIT_TRUE);
    seg_2->set_inactive(CUBIT_TRUE);
    return CUBIT_SUCCESS;
  }
  else if ( !point_1_on_seg_2 && !point_2_on_seg_1 &&
            point_0_on_seg_2 && point_3_on_seg_1 )
  {
    PRINT_DEBUG_129("Found OVERLAP_PART_1_2\n");
    int_result = OVERLAP_PART_1_2;
    ImprintLineSegment *new_seg_0_3, *new_seg_3_1;
    ImprintLineSegment *new_seg_2_0, *new_seg_0_3_v2;
      //Copy point_3 to use on seg 1.
    CubitVector closest_point;
    closest_point_seg( imp_point_3, seg_1, closest_point);
    ImprintPointData *new_point_3 = new ImprintPointData(imp_point_3, closest_point);
    new_point_3->owner(seg_1->owner());
    allocatedPointData->append(new_point_3);
    if ( imp_point_3->is_owner_vertex() )
    {
      new_point_3->set_point_type(CREATE_NEW_VERTEX);
      if ( imp_point_3->get_point_type() != CREATE_NEW_VERTEX )
        imp_point_3->set_point_type(VERTEX_ON_BOTH_BOUNDARIES);
    }
    else
    {
      new_point_3->set_point_type(ON_BOTH_BOUNDARIES);
      imp_point_3->set_point_type(ON_BOTH_BOUNDARIES);
    }
      //split both segs 1 and 2.
    new_seg_0_3 = new ImprintLineSegment(imp_point_0, new_point_3,
                                         seg_1->owner());
    new_seg_3_1 = new ImprintLineSegment(new_point_3, imp_point_1,
                                         seg_1->owner());
    allocatedLineData->append(new_seg_0_3);
    allocatedLineData->append(new_seg_3_1);
      //Copy point_0 to use on seg 2.
    closest_point_seg( imp_point_0, seg_2, closest_point);
    ImprintPointData *new_point_0 = new ImprintPointData(imp_point_0, closest_point);
    new_point_0->owner(seg_2->owner());
    allocatedPointData->append(new_point_0);
    if ( imp_point_0->is_owner_vertex() )
    {
      new_point_0->set_point_type(CREATE_NEW_VERTEX);
      if ( imp_point_0->get_point_type() != CREATE_NEW_VERTEX )
        imp_point_0->set_point_type(VERTEX_ON_BOTH_BOUNDARIES);
    }
    else
    {
      new_point_0->set_point_type(ON_BOTH_BOUNDARIES);
      imp_point_0->set_point_type(ON_BOTH_BOUNDARIES);
    }
    new_seg_2_0 = new ImprintLineSegment(imp_point_2, new_point_0,
                                         seg_2->owner());
    new_seg_0_3_v2 = new ImprintLineSegment(new_point_0, imp_point_3,
                                            seg_2->owner());
    allocatedLineData->append(new_seg_2_0);
    allocatedLineData->append(new_seg_0_3_v2);

    new_segments[0] = new_seg_0_3;
    new_segments[1] = new_seg_3_1;
    new_segments[2] = new_seg_2_0;
    new_segments[3] = new_seg_0_3_v2;
      //update both lists.
    update_list(new_segments, seg_1, seg_2, CUBIT_TRUE);
    update_list(new_segments, seg_1, seg_2, CUBIT_FALSE);
    seg_1->set_inactive(CUBIT_TRUE);
    seg_2->set_inactive(CUBIT_TRUE);
    return CUBIT_SUCCESS;
  }
  else if ( !point_1_on_seg_2 && point_2_on_seg_1 &&
            point_0_on_seg_2 && !point_3_on_seg_1 )
  {
    PRINT_DEBUG_129("Found OVERLAP_PART_1_3\n");
    int_result = OVERLAP_PART_1_3;
    ImprintLineSegment *new_seg_0_2, *new_seg_2_1;
    ImprintLineSegment *new_seg_2_0, *new_seg_0_3;
      //Copy point_3 to use on seg 1.
    CubitVector closest_point;
    closest_point_seg( imp_point_2, seg_1, closest_point);
    ImprintPointData *new_point_2 = new ImprintPointData(imp_point_2, closest_point);
    new_point_2->owner(seg_1->owner());
    allocatedPointData->append(new_point_2);
    if ( imp_point_2->is_owner_vertex() )
    {
      new_point_2->set_point_type(CREATE_NEW_VERTEX);
      if ( imp_point_2->get_point_type() != CREATE_NEW_VERTEX )
        imp_point_2->set_point_type(VERTEX_ON_BOTH_BOUNDARIES);
    }
    else
    {
      new_point_2->set_point_type(ON_BOTH_BOUNDARIES);
      imp_point_2->set_point_type(ON_BOTH_BOUNDARIES);
    }
      //split both segs 1 and 2.
    new_seg_0_2 = new ImprintLineSegment(imp_point_0, new_point_2,
                                         seg_1->owner());
    new_seg_2_1 = new ImprintLineSegment(new_point_2, imp_point_1,
                                         seg_1->owner());
    allocatedLineData->append(new_seg_0_2);
    allocatedLineData->append(new_seg_2_1);
      //Copy point_0 to use on seg 2.
    closest_point_seg( imp_point_0, seg_2, closest_point);
    ImprintPointData *new_point_0 = new ImprintPointData(imp_point_0, closest_point);
    new_point_0->owner(seg_2->owner());
    allocatedPointData->append(new_point_0);
    if ( imp_point_0->is_owner_vertex() )
    {
      new_point_0->set_point_type(CREATE_NEW_VERTEX);
      if ( imp_point_0->get_point_type() != CREATE_NEW_VERTEX )
        imp_point_0->set_point_type(VERTEX_ON_BOTH_BOUNDARIES);
    }
    else
    {
      new_point_0->set_point_type(ON_BOTH_BOUNDARIES);
      imp_point_0->set_point_type(ON_BOTH_BOUNDARIES);
    }
    new_seg_2_0 = new ImprintLineSegment(imp_point_2, new_point_0,
                                         seg_2->owner());
    new_seg_0_3 = new ImprintLineSegment(new_point_0, imp_point_3,
                                         seg_2->owner());
    allocatedLineData->append(new_seg_2_0);
    allocatedLineData->append(new_seg_0_3);

    new_segments[0] = new_seg_0_2;
    new_segments[1] = new_seg_2_1;
    new_segments[2] = new_seg_2_0;
    new_segments[3] = new_seg_0_3;
      //update both lists.
    update_list(new_segments, seg_1, seg_2, CUBIT_TRUE);
    update_list(new_segments, seg_1, seg_2, CUBIT_FALSE);
    seg_1->set_inactive(CUBIT_TRUE);
    seg_2->set_inactive(CUBIT_TRUE);
    return CUBIT_SUCCESS;
  }

    //Now test for the overlap all conditions.
  else if ( !point_0_on_seg_2 && !point_1_on_seg_2 &&
            point_2_on_seg_1 && point_3_on_seg_1 )
  {
      //First we know that for either of the next cases, we
      //have to copy points 2 and 3, do that first then
      //decide the particular intersection case.
      //copy points 2 and 3 to use on segment 1.
    CubitVector closest_point;
    closest_point_seg( imp_point_2, seg_1, closest_point);
    ImprintPointData *new_point_2 = new ImprintPointData(imp_point_2, closest_point);
    new_point_2->owner(seg_1->owner());
    allocatedPointData->append(new_point_2);
    if ( imp_point_2->is_owner_vertex() )
    {
      new_point_2->set_point_type(CREATE_NEW_VERTEX);
      if ( imp_point_2->get_point_type() != CREATE_NEW_VERTEX )
        imp_point_2->set_point_type(VERTEX_ON_BOTH_BOUNDARIES);
    }
    else
    {
      new_point_2->set_point_type(ON_BOTH_BOUNDARIES);
      imp_point_2->set_point_type(ON_BOTH_BOUNDARIES);
    }
    closest_point_seg( imp_point_3, seg_1, closest_point);
    ImprintPointData *new_point_3 = new ImprintPointData(imp_point_3, closest_point);
    new_point_3->owner(seg_1->owner());
    allocatedPointData->append(new_point_3);
    if ( imp_point_3->is_owner_vertex() )
    {
      new_point_3->set_point_type(CREATE_NEW_VERTEX);
      if ( imp_point_3->get_point_type() != CREATE_NEW_VERTEX )
        imp_point_3->set_point_type(VERTEX_ON_BOTH_BOUNDARIES);
      imp_point_3->set_point_type(VERTEX_ON_BOTH_BOUNDARIES);
    }
    else
    {
      new_point_3->set_point_type(ON_BOTH_BOUNDARIES);
      imp_point_3->set_point_type(ON_BOTH_BOUNDARIES);
    }
      //IMPORTANT: Now breaking up the line segments is special for the
      //overlap all case.  This is the one situation where
      //the one exterior segment gets broken into *3* new segments.
      //Every other instance in this function only breaks each segment
      //into two...
      //To distinguish between OVERLAP_ALL_0_2_3_1 and OVERLAP_ALL_0_3_2_1,
      //test to see which value is closer to point_0 (2 or 3).
    CubitVector point_0 = imp_point_0->coordinates();
    CubitVector point_2 = imp_point_2->coordinates();
    CubitVector point_3 = imp_point_3->coordinates();
    double dist_1_sq = (point_0 - point_2).length_squared();
    double dist_2_sq = (point_0 - point_3).length_squared();
    if ( dist_1_sq < dist_2_sq )
    {
      PRINT_DEBUG_129("Found OVERLAP_ALL_0_2_3_1\n");
      int_result = OVERLAP_ALL_0_2_3_1;
      ImprintLineSegment *new_seg_0_2, *new_seg_2_3, *new_seg_3_1;
        //split seg_1
      new_seg_0_2 = new ImprintLineSegment(imp_point_0, new_point_2,
                                           seg_1->owner());
        //Yes this is really just a copy of seg_2, but do it this
        //way to make the data separate.
      new_seg_2_3 = new ImprintLineSegment(new_point_2, new_point_3,
                                           seg_1->owner());
      new_seg_3_1 = new ImprintLineSegment(new_point_3, imp_point_1,
                                           seg_1->owner());
      allocatedLineData->append(new_seg_0_2);
      allocatedLineData->append(new_seg_2_3);
      allocatedLineData->append(new_seg_3_1);
      new_segments[0] = new_seg_0_2;
      new_segments[1] = new_seg_2_3;
      new_segments[2] = new_seg_3_1;
        //Don't use the update_list function here
        //since this is our special case...

        //update list 1.
        //prev
      new_seg_0_2->set_prev(seg_1->get_prev());
      seg_1->get_prev()->set_next(new_seg_0_2);
        //next
      new_seg_0_2->set_next(new_seg_2_3);
      new_seg_2_3->set_prev(new_seg_0_2);
      new_seg_2_3->set_next(new_seg_3_1);
      new_seg_3_1->set_prev(new_seg_2_3);
      new_seg_3_1->set_next(seg_1->get_next());
      seg_1->get_next()->set_prev(new_seg_3_1);
        //seg_1 will no longer be used
      seg_1->set_inactive(CUBIT_TRUE);
      return CUBIT_SUCCESS;
    }
    else
    {
      PRINT_DEBUG_129("Found OVERLAP_ALL_0_3_2_1\n");
      int_result = OVERLAP_ALL_0_3_2_1;
      ImprintLineSegment *new_seg_0_3, *new_seg_3_2, *new_seg_2_1;

      new_seg_0_3 = new ImprintLineSegment(imp_point_0, new_point_3,
                                           seg_1->owner());
      new_seg_3_2 = new ImprintLineSegment(new_point_3, new_point_2,
                                           seg_1->owner());
      new_seg_2_1 = new ImprintLineSegment(new_point_2, imp_point_1,
                                           seg_1->owner());
      allocatedLineData->append(new_seg_0_3);
      allocatedLineData->append(new_seg_3_2);
      allocatedLineData->append(new_seg_2_1);
      new_segments[0] = new_seg_0_3;
      new_segments[1] = new_seg_3_2;
      new_segments[2] = new_seg_2_1;
        //don't use special function since this is special case.
        //update list 1
        //prev
      new_seg_0_3->set_prev(seg_1->get_prev());
      seg_1->get_prev()->set_next(new_seg_0_3);
        //next
      new_seg_0_3->set_next(new_seg_3_2);
      new_seg_3_2->set_prev(new_seg_0_3);
      new_seg_3_2->set_next(new_seg_2_1);
      new_seg_2_1->set_prev(new_seg_3_2);
      new_seg_2_1->set_next(seg_1->get_next());
      seg_1->get_next()->set_prev(new_seg_2_1);
        //seg_1 will no longer be used
      seg_1->set_inactive(CUBIT_TRUE);
      return CUBIT_SUCCESS;
    }
  }
  else if ( point_0_on_seg_2 && point_1_on_seg_2 &&
            !point_2_on_seg_1 && !point_3_on_seg_1 )
  {
      //First we know that for either of the next cases, we
      //have to copy points 0 and 1, do that first then
      //decide the particular intersection case.
      //copy points 0 and 1 to use on segment 2.
    CubitVector closest_point;
    closest_point_seg( imp_point_0, seg_2, closest_point);
    ImprintPointData *new_point_0 = new ImprintPointData(imp_point_0, closest_point);
    new_point_0->owner(seg_2->owner());
    allocatedPointData->append(new_point_0);
    if ( imp_point_0->is_owner_vertex() )
    {
      new_point_0->set_point_type(CREATE_NEW_VERTEX);
      if ( imp_point_0->get_point_type() != CREATE_NEW_VERTEX )
        imp_point_0->set_point_type(VERTEX_ON_BOTH_BOUNDARIES);
    }
    else
    {
      new_point_0->set_point_type(ON_BOTH_BOUNDARIES);
      imp_point_0->set_point_type(ON_BOTH_BOUNDARIES);
    }
    closest_point_seg( imp_point_1, seg_2, closest_point);
    ImprintPointData *new_point_1 = new ImprintPointData(imp_point_1, closest_point);
    new_point_1->owner(seg_2->owner());
    allocatedPointData->append(new_point_1);
    if ( imp_point_1->is_owner_vertex() )
    {
      new_point_1->set_point_type(CREATE_NEW_VERTEX);
      if ( imp_point_1->get_point_type() != CREATE_NEW_VERTEX )
        imp_point_1->set_point_type(VERTEX_ON_BOTH_BOUNDARIES);
    }
    else
    {
      new_point_1->set_point_type(ON_BOTH_BOUNDARIES);
      imp_point_1->set_point_type(ON_BOTH_BOUNDARIES);
    }
      //IMPORTANT: Now breaking up the line segments is special for the
      //overlap all case.  This is the one situation where
      //the one exterior segment gets broken into *3* new segments.
      //Every other instance in this function only breaks each segment
      //into two...
      //To distinguish between OVERLAP_ALL_2_0_1_3 and OVERLAP_ALL_2_1_0_3,
      //test to see which value is closer to point_2 (0 or 1).
    CubitVector point_0 = imp_point_0->coordinates();
    CubitVector point_1 = imp_point_1->coordinates();
    CubitVector point_2 = imp_point_2->coordinates();
    double dist_1_sq = (point_2 - point_0).length_squared();
    double dist_2_sq = (point_2 - point_1).length_squared();
    if ( dist_1_sq < dist_2_sq )
    {
      PRINT_DEBUG_129("Found OVERLAP_ALL_2_0_1_3\n");
      int_result = OVERLAP_ALL_2_0_1_3;
      ImprintLineSegment *new_seg_2_0, *new_seg_0_1, *new_seg_1_3;

      new_seg_2_0 = new ImprintLineSegment(imp_point_2, new_point_0,
                                           seg_2->owner());
        //Yes this is really just a copy of seg_2, but do it this
        //way to make the data separate.
      new_seg_0_1 = new ImprintLineSegment(new_point_0, new_point_1,
                                           seg_2->owner());
      new_seg_1_3 = new ImprintLineSegment(new_point_1, imp_point_3,
                                           seg_2->owner());
      allocatedLineData->append(new_seg_2_0);
      allocatedLineData->append(new_seg_0_1);
      allocatedLineData->append(new_seg_1_3);
      new_segments[1] = new_seg_2_0;
      new_segments[2] = new_seg_0_1;
      new_segments[3] = new_seg_1_3;
        //Don't use update list function since this is special case.
        //update list 2.
        //prev
      new_seg_2_0->set_prev(seg_2->get_prev());
      seg_2->get_prev()->set_next(new_seg_2_0);
        //next
      new_seg_2_0->set_next(new_seg_0_1);
      new_seg_0_1->set_prev(new_seg_2_0);
      new_seg_0_1->set_next(new_seg_1_3);
      new_seg_1_3->set_prev(new_seg_0_1);
      new_seg_1_3->set_next(seg_2->get_next());
      seg_2->get_next()->set_prev(new_seg_1_3);
        //seg_2 will no longer be used
      seg_2->set_inactive(CUBIT_TRUE);
      return CUBIT_SUCCESS;
    }
    else
    {
      PRINT_DEBUG_129("Found OVERLAP_ALL_2_1_0_3\n");
      int_result = OVERLAP_ALL_2_1_0_3;
      ImprintLineSegment *new_seg_2_1, *new_seg_1_0, *new_seg_0_3;

      new_seg_2_1 = new ImprintLineSegment(imp_point_2, new_point_1,
                                           seg_2->owner());
      new_seg_1_0 = new ImprintLineSegment(new_point_1, new_point_0,
                                           seg_2->owner());
      new_seg_0_3 = new ImprintLineSegment(new_point_0, imp_point_3,
                                           seg_2->owner());
      allocatedLineData->append(new_seg_2_1);
      allocatedLineData->append(new_seg_1_0);
      allocatedLineData->append(new_seg_0_3);
      new_segments[1] = new_seg_2_1;
      new_segments[2] = new_seg_1_0;
      new_segments[3] = new_seg_0_3;
        //Don't use update_list function since this is special case.
        //update list 2.
        //prev
      new_seg_2_1->set_prev(seg_2->get_prev());
      seg_2->get_prev()->set_next(new_seg_2_1);
        //next
      new_seg_2_1->set_next(new_seg_1_0);
      new_seg_1_0->set_prev(new_seg_2_1);
      new_seg_1_0->set_next(new_seg_0_3);
      new_seg_0_3->set_prev(new_seg_1_0);
      new_seg_0_3->set_next(seg_2->get_next());
      seg_2->get_next()->set_prev(new_seg_0_3);
        //seg_2 will no longer be used
      seg_2->set_inactive(CUBIT_TRUE);
      return CUBIT_SUCCESS;
    }
  }
  if ( point_0_on_seg_2 || point_1_on_seg_2 ||
       point_2_on_seg_1 || point_3_on_seg_1 )
  {
    PRINT_ERROR("Intersection problem, case not handled.\n");
    return CUBIT_FAILURE;
  }
    //Okay, finally test for the cross intersection.  If this isnt'
    //within tolerance, then these two segments couldn't
    //posibly intersect!
  IntersectionTool int_tool(GEOMETRY_RESABS); //use GEOMETRY_RESABS
    //rather than myTolerance for calculations...
  CubitVector point_0 = imp_point_0->coordinates();
  CubitVector point_1 = imp_point_1->coordinates();
  CubitVector point_2 = imp_point_2->coordinates();
  CubitVector point_3 = imp_point_3->coordinates();
  CubitVector closest_point_seg_1, closest_point_seg_2;
  double sc, tc;
  CubitStatus stat = int_tool.closest_points_on_segments(point_0, point_1,
                                                         point_2, point_3,
                                                         closest_point_seg_1,
                                                         closest_point_seg_2,
                                                         sc, tc);
  if (stat != CUBIT_SUCCESS )
  {
    PRINT_ERROR("Problems calculation closest points on "
                "segments for boundary imprinting.\n");
    return CUBIT_FAILURE;
  }
    //Make sure the closest points aren't the end points.  If they are
    //and we are within tolerance, it may be that the tolerance is too big
    //cause we shouldn't be at a cross if the closest point is an end point...
  if ( sc > 0. && sc < 1. && tc > 0. && tc < 1. &&
       closest_point_seg_1.within_tolerance(closest_point_seg_2,myTolerance) )
  {
      //okay, these guys actually do intersect!
    PRINT_DEBUG_129("Found CROSS_INTERSECT\n");
    int_result = CROSS_INTERSECT;
    ImprintPointData *new_point_seg_1, *new_point_seg_2;
    new_point_seg_1 = new ImprintPointData(closest_point_seg_1.x(),
                                           closest_point_seg_1.y(),
                                           closest_point_seg_1.z());
    allocatedPointData->append(new_point_seg_1);
    new_point_seg_1->owner(seg_1->owner());
    new_point_seg_1->set_point_type(CREATE_NEW_VERTEX);
    new_point_seg_2 = new ImprintPointData(closest_point_seg_2.x(),
                                           closest_point_seg_2.y(),
                                           closest_point_seg_2.z());
    allocatedPointData->append(new_point_seg_2);
    new_point_seg_2->owner(seg_2->owner());
    new_point_seg_2->set_point_type(CREATE_NEW_VERTEX);

    new_point_seg_2->set_matching_point(new_point_seg_1);
    new_point_seg_1->set_matching_point(new_point_seg_2);
    
    new_segments[0] = new ImprintLineSegment(imp_point_0, new_point_seg_1,
                                             seg_1->owner());
    new_segments[1] = new ImprintLineSegment(new_point_seg_1, imp_point_1,
                                             seg_1->owner());
    new_segments[2] = new ImprintLineSegment(imp_point_2, new_point_seg_2,
                                             seg_2->owner());
    new_segments[3] = new ImprintLineSegment(new_point_seg_2, imp_point_3,
                                             seg_2->owner());
    allocatedLineData->append(new_segments[0]);
    allocatedLineData->append(new_segments[1]);
    allocatedLineData->append(new_segments[2]);
    allocatedLineData->append(new_segments[3]);
      //now update the lists.
      //update list 1.
    update_list(new_segments, seg_1, seg_2, CUBIT_TRUE);
      //update list 2.
    update_list(new_segments, seg_1, seg_2, CUBIT_FALSE);
      //seg_1 and seg_2 will no longer be used
    seg_1->set_inactive(CUBIT_TRUE);
    seg_2->set_inactive(CUBIT_TRUE);
    return CUBIT_SUCCESS;
  }
  else
  {
      //All that for nothing!  Some how I need to find this out
      //earlier to avoid all the comparisions...  Will have
      //to check at least for the intersection though..
    PRINT_DEBUG_129("Found NO_INTERSECT\n");
    int_result = NO_INTERSECT;
    return CUBIT_SUCCESS;
  }
}

void ImprintBoundaryTool::update_list(ImprintLineSegment **new_segments,
                                      ImprintLineSegment *seg_1,
                                      ImprintLineSegment *seg_2,
                                      const CubitBoolean list_1 )
{
  if (list_1)
  {
    new_segments[0]->set_prev(seg_1->get_prev());
    seg_1->get_prev()->set_next(new_segments[0]);
    new_segments[0]->set_next(new_segments[1]);
    new_segments[1]->set_prev(new_segments[0]);
    new_segments[1]->set_next(seg_1->get_next());
    seg_1->get_next()->set_prev(new_segments[1]);
  }
  else
  {
    new_segments[2]->set_prev(seg_2->get_prev());
    seg_2->get_prev()->set_next(new_segments[2]);
    new_segments[2]->set_next(new_segments[3]);
    new_segments[3]->set_prev(new_segments[2]);
    new_segments[3]->set_next(seg_2->get_next());
    seg_2->get_next()->set_prev(new_segments[3]);
  }
}
CubitStatus ImprintBoundaryTool::case_0_2_equal( ImprintLineSegment *seg_1,
                                                 ImprintLineSegment *seg_2,
                                                 ImprintPointData* imp_point_0,
                                                 ImprintPointData* imp_point_1,
                                                 ImprintPointData* imp_point_2,
                                                 ImprintPointData* imp_point_3,
                                                 MatchType type_1, MatchType type_3,
                                                 IntersectResult &int_result,
                                                 ImprintLineSegment **new_segments)
{
    //Okay at this point we know that point_0 and point_2 are
    //equal, and we know that point_1 doesn't match an end_point.
    //We know then that this is either an L or a overlap-join.
    //Test to see if point_1 is on seg_2 or point_3 is on seg_1.
    //If neither then we have an L, other wise an overlap-join.
  CubitBoolean point_1_on_seg_2 =CUBIT_FALSE;
  CubitBoolean point_3_on_seg_1 =CUBIT_FALSE;
  if ( on_interior_segment(imp_point_1, seg_2) )
    point_1_on_seg_2 = CUBIT_TRUE;
  if ( on_interior_segment(imp_point_3, seg_1) )
    point_3_on_seg_1 = CUBIT_TRUE;
  
  if ( !point_1_on_seg_2 && !point_3_on_seg_1 )
  {
    PRINT_DEBUG_129("Found L_INTERSECT_0_2\n");
    int_result = L_INTERSECT_0_2;
    set_type_for_L(imp_point_0, imp_point_2);
    return CUBIT_SUCCESS;
  }
  else if (point_1_on_seg_2 )
  {
    if (type_1 == MATCH_1_6 )
    {
      double dist_1 = (imp_point_1->coordinates() -
                       seg_2->get_prev()->get_start()->coordinates()).length_squared();
      CubitVector point_on_seg_2;
      closest_point_interior_seg(imp_point_1, seg_2, point_on_seg_2);
      double dist_2 = (imp_point_1->coordinates() - point_on_seg_2).length_squared();
      if ( dist_1 < dist_2 )
      {
        PRINT_DEBUG_129("Found L_INTERSECT_0_2\n");
        int_result = L_INTERSECT_0_2;
        set_type_for_L(imp_point_0, imp_point_2);
        return CUBIT_SUCCESS;
      }
    }
    else if (type_1 == MATCH_1_7 )
    {
      double dist_1 = (imp_point_1->coordinates() -
                       seg_2->get_next()->get_end()->coordinates()).length_squared();
      CubitVector point_on_seg_2;
      closest_point_interior_seg(imp_point_1, seg_2, point_on_seg_2);
      double dist_2 = (imp_point_1->coordinates() - point_on_seg_2).length_squared();
      if ( dist_1 < dist_2 )
      {
        PRINT_DEBUG_129("Found L_INTERSECT_0_2\n");
        int_result = L_INTERSECT_0_2;
        set_type_for_L(imp_point_0, imp_point_2);
        return CUBIT_SUCCESS;
      }
    }
  }
  else if (point_3_on_seg_1 )
  {
    if (type_3 == MATCH_3_4 )
    {
      double dist_1 = (imp_point_3->coordinates() -
                       seg_1->get_prev()->get_start()->coordinates()).length_squared();
      CubitVector point_on_seg_1;
      closest_point_interior_seg(imp_point_3, seg_1, point_on_seg_1);
      double dist_2 = (imp_point_3->coordinates() - point_on_seg_1).length_squared();
      if ( dist_1 < dist_2 )
      {
        PRINT_DEBUG_129("Found L_INTERSECT_0_2\n");
        int_result = L_INTERSECT_0_2;
        set_type_for_L(imp_point_0, imp_point_2);
        return CUBIT_SUCCESS;
      }
    }
    else if (type_3 == MATCH_3_5 )
    {
      double dist_1 = (imp_point_3->coordinates() -
                       seg_1->get_next()->get_end()->coordinates()).length_squared();
      CubitVector point_on_seg_1;
      closest_point_interior_seg(imp_point_3, seg_1, point_on_seg_1);
      double dist_2 = (imp_point_3->coordinates() - point_on_seg_1).length_squared();
      if ( dist_1 < dist_2 )
      {
        PRINT_DEBUG_129("Found L_INTERSECT_0_2\n");
        int_result = L_INTERSECT_0_2;
        set_type_for_L(imp_point_0, imp_point_2);
        return CUBIT_SUCCESS;
      }
    }
  }
  if ( point_1_on_seg_2 )
  {
      //This is an overlap_join.
    PRINT_DEBUG_129("Found OVERALP_JOIN_02_1_3\n");
    int_result = OVERLAP_JOIN_02_1_3;
      //Create the new segments, the first segment remains the same so
      //leave positions 0 and 1 empty.
    ImprintLineSegment *new_seg_2_1, *new_seg_1_3;
      //copy point 1 to use on segment 2
    CubitVector closest_point;
    closest_point_seg( imp_point_1, seg_2, closest_point);
    ImprintPointData *new_point_1 = new ImprintPointData(imp_point_1, closest_point);
    new_point_1->owner(seg_2->owner());
    allocatedPointData->append(new_point_1);
    if ( imp_point_1->is_owner_vertex() )
    {
      new_point_1->set_point_type(CREATE_NEW_VERTEX);
      if ( imp_point_1->get_point_type() != CREATE_NEW_VERTEX )
        imp_point_1->set_point_type(VERTEX_ON_BOTH_BOUNDARIES);
    }
    else
    {
      new_point_1->set_point_type(ON_BOTH_BOUNDARIES);
      imp_point_1->set_point_type(ON_BOTH_BOUNDARIES);
    }
    new_seg_2_1 = new ImprintLineSegment(imp_point_2, new_point_1,
                                         seg_2->owner());
    new_seg_1_3 = new ImprintLineSegment(new_point_1, imp_point_3,
                                         seg_2->owner());
    new_segments[2] = new_seg_2_1;
    new_segments[3] = new_seg_1_3;
    allocatedLineData->append(new_segments[2]);
    allocatedLineData->append(new_segments[3]);
      //update the lists.
    update_list(new_segments, seg_1, seg_2, CUBIT_FALSE);
    seg_2->set_inactive(CUBIT_TRUE);
    return CUBIT_SUCCESS;
  }
  else if ( point_3_on_seg_1 )
  {
      //This is an overlap_join.
    PRINT_DEBUG_129("Found OVERALP_JOIN_02_3_1\n");
    int_result = OVERLAP_JOIN_02_3_1;
      //Create the new segments, the first segment remains the same so
      //leave positions 0 and 1 empty.
    ImprintLineSegment *new_seg_0_3, *new_seg_3_1;
      //copy point 3 to use on segment 1
    CubitVector closest_point;
    closest_point_seg( imp_point_3, seg_1, closest_point);
    ImprintPointData *new_point_3 = new ImprintPointData(imp_point_3, closest_point);
    new_point_3->owner(seg_1->owner());
    allocatedPointData->append(new_point_3);
    if ( imp_point_3->is_owner_vertex() )
    {
      new_point_3->set_point_type(CREATE_NEW_VERTEX);
      if ( imp_point_3->get_point_type() != CREATE_NEW_VERTEX )
        imp_point_3->set_point_type(VERTEX_ON_BOTH_BOUNDARIES);
    }
    else
    {
      new_point_3->set_point_type(ON_BOTH_BOUNDARIES);
      imp_point_3->set_point_type(ON_BOTH_BOUNDARIES);
    }

    new_seg_0_3 = new ImprintLineSegment(imp_point_0, new_point_3,
                                         seg_1->owner());
    new_seg_3_1 = new ImprintLineSegment(new_point_3, imp_point_1,
                                         seg_1->owner());
    new_segments[0] = new_seg_0_3;
    new_segments[1] = new_seg_3_1;
    allocatedLineData->append(new_segments[0]);
    allocatedLineData->append(new_segments[1]);
      //update the lists
    update_list(new_segments, seg_1, seg_2, CUBIT_TRUE);
    seg_1->set_inactive(CUBIT_TRUE);
    return CUBIT_SUCCESS;
  }
  else
  {
      //For this case we must have an L condition.
    set_type_for_L(imp_point_0, imp_point_2);
    PRINT_DEBUG_129("Found L_INTERSECT_0_2\n");
    int_result = L_INTERSECT_0_2;
    return CUBIT_SUCCESS;
  }
}
CubitStatus ImprintBoundaryTool::case_0_3_equal( ImprintLineSegment *seg_1,
                                                 ImprintLineSegment *seg_2,
                                                 ImprintPointData* imp_point_0,
                                                 ImprintPointData* imp_point_1,
                                                 ImprintPointData* imp_point_2,
                                                 ImprintPointData* imp_point_3,
                                                 MatchType type_1, MatchType type_2,
                                                 IntersectResult &int_result,
                                                 ImprintLineSegment **new_segments)
{
    //Okay at this point we know that point_0 and point_3 are
    //equal, and we know that point_1 doesn't match an end_point.
    //We know then that this is either an L or a overlap-join.
    //Test to see if point_1 is on seg_2 or point_2 is on seg_1.
    //If neither then we have an L, other wise an overlap-join.
  CubitBoolean point_1_on_seg_2 =CUBIT_FALSE;
  CubitBoolean point_2_on_seg_1 =CUBIT_FALSE;
  if ( on_interior_segment(imp_point_1, seg_2) )
    point_1_on_seg_2 = CUBIT_TRUE;
  if ( on_interior_segment(imp_point_2, seg_1) )
    point_2_on_seg_1 = CUBIT_TRUE;
  
  if ( !point_1_on_seg_2 && !point_2_on_seg_1 )
  {
    PRINT_DEBUG_129("Found L_INTERSECT_0_3\n");
    int_result = L_INTERSECT_0_3;
    set_type_for_L(imp_point_0, imp_point_3);
    return CUBIT_SUCCESS;
  }
  else if (point_1_on_seg_2 )
  {
    if (type_1 == MATCH_1_6 )
    {
      double dist_1 = (imp_point_1->coordinates() -
                       seg_2->get_prev()->get_start()->coordinates()).length_squared();
      CubitVector point_on_seg_2;
      closest_point_interior_seg(imp_point_1, seg_2, point_on_seg_2);
      double dist_2 = (imp_point_1->coordinates() - point_on_seg_2).length_squared();
      if ( dist_1 < dist_2 )
      {
        PRINT_DEBUG_129("Found L_INTERSECT_0_3\n");
        int_result = L_INTERSECT_0_3;
        set_type_for_L(imp_point_0, imp_point_3);
        return CUBIT_SUCCESS;
      }
    }
    else if (type_1 == MATCH_1_7 )
    {
      double dist_1 = (imp_point_1->coordinates() -
                       seg_2->get_next()->get_end()->coordinates()).length_squared();
      CubitVector point_on_seg_2;
      closest_point_interior_seg(imp_point_1, seg_2, point_on_seg_2);
      double dist_2 = (imp_point_1->coordinates() - point_on_seg_2).length_squared();
      if ( dist_1 < dist_2 )
      {
        PRINT_DEBUG_129("Found L_INTERSECT_0_3\n");
        int_result = L_INTERSECT_0_3;
        set_type_for_L(imp_point_0, imp_point_3);
        return CUBIT_SUCCESS;
      }
    }
  }
  else if (point_2_on_seg_1 )
  {
    if (type_2 == MATCH_2_4 )
    {
      double dist_1 = (imp_point_2->coordinates() -
                       seg_1->get_prev()->get_start()->coordinates()).length_squared();
      CubitVector point_on_seg_1;
      closest_point_interior_seg(imp_point_2, seg_1, point_on_seg_1);
      double dist_2 = (imp_point_2->coordinates() - point_on_seg_1).length_squared();
      if ( dist_1 < dist_2 )
      {
        PRINT_DEBUG_129("Found L_INTERSECT_0_3\n");
        int_result = L_INTERSECT_0_3;
        set_type_for_L(imp_point_0, imp_point_3);
        return CUBIT_SUCCESS;
      }
    }
    else if (type_2 == MATCH_2_5 )
    {
      double dist_1 = (imp_point_2->coordinates() -
                       seg_1->get_next()->get_end()->coordinates()).length_squared();
      CubitVector point_on_seg_1;
      closest_point_interior_seg(imp_point_2, seg_1, point_on_seg_1);
      double dist_2 = (imp_point_2->coordinates() - point_on_seg_1).length_squared();
      if ( dist_1 < dist_2 )
      {
        PRINT_DEBUG_129("Found L_INTERSECT_0_3\n");
        int_result = L_INTERSECT_0_3;
        set_type_for_L(imp_point_0, imp_point_3);
        return CUBIT_SUCCESS;
      }
    }
  }

  if ( point_1_on_seg_2 )
  {
      //This is an overlap_join.
    PRINT_DEBUG_129("Found OVERALP_JOIN_03_1_2\n");
    //double dist = (imp_point_1->coordinates() - imp_point_2->coordinates()).length();
    int_result = OVERLAP_JOIN_03_1_2;
      //Create the new segments, the first segment remains the same so
      //leave positions 0 and 1 empty.
    ImprintLineSegment *new_seg_2_1, *new_seg_1_3;
      //copy point 1 to use on segment 2
    CubitVector closest_point;
    closest_point_seg( imp_point_1, seg_2, closest_point);
    ImprintPointData *new_point_1 = new ImprintPointData(imp_point_1, closest_point);
    new_point_1->owner(seg_2->owner());
    allocatedPointData->append(new_point_1);
    if ( imp_point_1->is_owner_vertex() )
    {
      new_point_1->set_point_type(CREATE_NEW_VERTEX);
      if ( imp_point_1->get_point_type() != CREATE_NEW_VERTEX )
        imp_point_1->set_point_type(VERTEX_ON_BOTH_BOUNDARIES);
    }
    else
    {
      new_point_1->set_point_type(ON_BOTH_BOUNDARIES);
      imp_point_1->set_point_type(ON_BOTH_BOUNDARIES);
    }

    new_seg_2_1 = new ImprintLineSegment(imp_point_2, new_point_1,
                                         seg_2->owner());
    new_seg_1_3 = new ImprintLineSegment(new_point_1, imp_point_3,
                                         seg_2->owner());
    new_segments[2] = new_seg_2_1;
    new_segments[3] = new_seg_1_3;
    allocatedLineData->append(new_segments[2]);
    allocatedLineData->append(new_segments[3]);
      //update the list.
    update_list(new_segments, seg_1, seg_2, CUBIT_FALSE);
    seg_2->set_inactive(CUBIT_TRUE);
    return CUBIT_SUCCESS;
  }
  else if ( point_2_on_seg_1 )
  {
      //This is an overlap_join.
    PRINT_DEBUG_129("Found OVERALP_JOIN_03_2_1\n");
    int_result = OVERLAP_JOIN_03_2_1;
      //Create the new segments, the first segment remains the same so
      //leave positions 0 and 1 empty.
    ImprintLineSegment *new_seg_0_2, *new_seg_2_1;
      //copy point 1 to use on segment 2
    CubitVector closest_point;
    closest_point_seg( imp_point_2, seg_1, closest_point);
    ImprintPointData *new_point_2 = new ImprintPointData(imp_point_2, closest_point);
    new_point_2->owner(seg_1->owner());
    allocatedPointData->append(new_point_2);
    if ( imp_point_2->is_owner_vertex() )
    {
      new_point_2->set_point_type(CREATE_NEW_VERTEX);
      if ( imp_point_2->get_point_type() != CREATE_NEW_VERTEX )
        imp_point_2->set_point_type(VERTEX_ON_BOTH_BOUNDARIES);
    }
    else
    {
      new_point_2->set_point_type(ON_BOTH_BOUNDARIES);
      imp_point_2->set_point_type(ON_BOTH_BOUNDARIES);
    }

    new_seg_0_2 = new ImprintLineSegment(imp_point_0, new_point_2,
                                         seg_1->owner());
    new_seg_2_1 = new ImprintLineSegment(new_point_2, imp_point_1,
                                         seg_1->owner());
    new_segments[0] = new_seg_0_2;
    new_segments[1] = new_seg_2_1;
    allocatedLineData->append(new_segments[0]);
    allocatedLineData->append(new_segments[1]);
      //update the list.
    update_list(new_segments, seg_1, seg_2, CUBIT_TRUE);
    seg_1->set_inactive(CUBIT_TRUE);
    return CUBIT_SUCCESS;
  }
  else
  {
      //For this case we must have an L condition.
    set_type_for_L(imp_point_0, imp_point_3);
    PRINT_DEBUG_129("Found L_INTERSECT_0_3\n");
    int_result = L_INTERSECT_0_3;
    return CUBIT_SUCCESS;
  }
}
CubitStatus ImprintBoundaryTool::case_1_2_equal( ImprintLineSegment *seg_1,
                                                 ImprintLineSegment *seg_2,
                                                 ImprintPointData* imp_point_0,
                                                 ImprintPointData* imp_point_1,
                                                 ImprintPointData* imp_point_2,
                                                 ImprintPointData* imp_point_3,
                                                 MatchType type_0, MatchType type_3,
                                                 IntersectResult &int_result,
                                                 ImprintLineSegment **new_segments)
{
    //Okay at this point we know that point_1 and point_2 are
    //equal, and we know that point_0 doesn't match an end_point.
    //We know then that this is either an L or a overlap-join.
    //Test to see if point_0 is on seg_2 or point_3 is on seg_1.
    //If neither then we have an L, other wise an overlap-join.
  CubitBoolean point_0_on_seg_2 =CUBIT_FALSE;
  CubitBoolean point_3_on_seg_1 =CUBIT_FALSE;
  if ( on_interior_segment(imp_point_0, seg_2) )
    point_0_on_seg_2 = CUBIT_TRUE;
  if ( on_interior_segment(imp_point_3, seg_1) )
    point_3_on_seg_1 = CUBIT_TRUE;
  
  if ( !point_0_on_seg_2 && !point_3_on_seg_1 )
  {
    PRINT_DEBUG_129("Found L_INTERSECT_1_2\n");
    int_result = L_INTERSECT_1_2;
    set_type_for_L(imp_point_1, imp_point_2);
    return CUBIT_SUCCESS;
  }
  else if (point_0_on_seg_2 )
  {
    if (type_0 == MATCH_0_6 )
    {
      double dist_1 = (imp_point_0->coordinates() -
                       seg_2->get_prev()->get_start()->coordinates()).length_squared();
      CubitVector point_on_seg_2;
      closest_point_interior_seg(imp_point_0, seg_2, point_on_seg_2);
      double dist_2 = (imp_point_0->coordinates() - point_on_seg_2).length_squared();
      if ( dist_1 < dist_2 )
      {
        PRINT_DEBUG_129("Found L_INTERSECT_1_2\n");
        int_result = L_INTERSECT_1_2;
        set_type_for_L(imp_point_1, imp_point_2);
        return CUBIT_SUCCESS;
      }
    }
    else if (type_0 == MATCH_0_7 )
    {
      double dist_1 = (imp_point_0->coordinates() -
                       seg_2->get_next()->get_end()->coordinates()).length_squared();
      CubitVector point_on_seg_2;
      closest_point_interior_seg(imp_point_0, seg_2, point_on_seg_2);
      double dist_2 = (imp_point_0->coordinates() - point_on_seg_2).length_squared();
      if ( dist_1 < dist_2 )
      {
        PRINT_DEBUG_129("Found L_INTERSECT_1_2\n");
        int_result = L_INTERSECT_1_2;
        set_type_for_L(imp_point_1, imp_point_2);
        return CUBIT_SUCCESS;
      }
    }
  }
  else if (point_3_on_seg_1 )
  {
    if (type_3 == MATCH_3_4 )
    {
      double dist_1 = (imp_point_3->coordinates() -
                       seg_1->get_prev()->get_start()->coordinates()).length_squared();
      CubitVector point_on_seg_1;
      closest_point_interior_seg(imp_point_3, seg_1, point_on_seg_1);
      double dist_2 = (imp_point_3->coordinates() - point_on_seg_1).length_squared();
      if ( dist_1 < dist_2 )
      {
        PRINT_DEBUG_129("Found L_INTERSECT_1_2\n");
        int_result = L_INTERSECT_1_2;
        set_type_for_L(imp_point_1, imp_point_2);
        return CUBIT_SUCCESS;
      }
    }
    else if (type_3 == MATCH_3_5 )
    {
      double dist_1 = (imp_point_3->coordinates() -
                       seg_1->get_next()->get_end()->coordinates()).length_squared();
      CubitVector point_on_seg_1;
      closest_point_interior_seg(imp_point_3, seg_1, point_on_seg_1);
      double dist_2 = (imp_point_3->coordinates() - point_on_seg_1).length_squared();
      if ( dist_1 < dist_2 )
      {
        PRINT_DEBUG_129("Found L_INTERSECT_1_2\n");
        int_result = L_INTERSECT_1_2;
        set_type_for_L(imp_point_1, imp_point_2);
        return CUBIT_SUCCESS;
      }
    }
  }
  if ( point_0_on_seg_2 )
  {
      //This is an overlap_join.
    PRINT_DEBUG_129("Found OVERALP_JOIN_12_0_3\n");
    int_result = OVERLAP_JOIN_12_0_3;
      //Create the new segments, the first segment remains the same so
      //leave positions 0 and 1 empty.
    ImprintLineSegment *new_seg_2_0, *new_seg_0_3;
      //copy point 1 to use on segment 2
    CubitVector closest_point;
    closest_point_seg( imp_point_0, seg_2, closest_point);
    ImprintPointData *new_point_0 = new ImprintPointData(imp_point_0, closest_point);
    new_point_0->owner(seg_2->owner());
    allocatedPointData->append(new_point_0);
    if ( imp_point_0->is_owner_vertex() )
    {
      new_point_0->set_point_type(CREATE_NEW_VERTEX);
      if ( imp_point_0->get_point_type() != CREATE_NEW_VERTEX )
        imp_point_0->set_point_type(VERTEX_ON_BOTH_BOUNDARIES);
    }
    else
    {
      new_point_0->set_point_type(ON_BOTH_BOUNDARIES);
      imp_point_0->set_point_type(ON_BOTH_BOUNDARIES);
    }

    new_seg_2_0 = new ImprintLineSegment(imp_point_2, new_point_0,
                                         seg_2->owner());
    new_seg_0_3 = new ImprintLineSegment(new_point_0, imp_point_3,
                                         seg_2->owner());
    new_segments[2] = new_seg_2_0;
    new_segments[3] = new_seg_0_3;
    allocatedLineData->append(new_segments[2]);
    allocatedLineData->append(new_segments[3]);
      //update the list.
    update_list(new_segments, seg_1, seg_2, CUBIT_FALSE);
    seg_2->set_inactive(CUBIT_TRUE);
    return CUBIT_SUCCESS;
  }
  else if ( point_3_on_seg_1 )
  {
      //This is an overlap_join.
    PRINT_DEBUG_129("Found OVERALP_JOIN_12_3_0\n");
    int_result = OVERLAP_JOIN_12_3_0;
      //Create the new segments, the first segment remains the same so
      //leave positions 0 and 1 empty.
    ImprintLineSegment *new_seg_3_1, *new_seg_0_3;
      //copy point 3 to use on segment 1
    CubitVector closest_point;
    closest_point_seg( imp_point_3, seg_1, closest_point);
    ImprintPointData *new_point_3 = new ImprintPointData(imp_point_3, closest_point);
    new_point_3->owner(seg_1->owner());
    allocatedPointData->append(new_point_3);
    if ( imp_point_3->is_owner_vertex() )
    {
      new_point_3->set_point_type(CREATE_NEW_VERTEX);
      if ( imp_point_3->get_point_type() != CREATE_NEW_VERTEX )
        imp_point_3->set_point_type(VERTEX_ON_BOTH_BOUNDARIES);
    }
    else
    {
      new_point_3->set_point_type(ON_BOTH_BOUNDARIES);
      imp_point_3->set_point_type(ON_BOTH_BOUNDARIES);
    }
    new_seg_0_3 = new ImprintLineSegment(imp_point_0, new_point_3,
                                         seg_1->owner());
    new_seg_3_1 = new ImprintLineSegment(new_point_3, imp_point_1,
                                         seg_1->owner());
    new_segments[0] = new_seg_0_3;
    new_segments[1] = new_seg_3_1;
    allocatedLineData->append(new_segments[0]);
    allocatedLineData->append(new_segments[1]);
      //update the list.
    update_list(new_segments, seg_1, seg_2, CUBIT_TRUE);
    seg_1->set_inactive(CUBIT_TRUE);
    return CUBIT_SUCCESS;
  }
  else
  {
      //For this case we must have an L condition.
    set_type_for_L(imp_point_1, imp_point_2);
    PRINT_DEBUG_129("Found L_INTERSECT_1_2\n");
    int_result = L_INTERSECT_1_2;
    return CUBIT_SUCCESS;
  }
}
CubitStatus ImprintBoundaryTool::determine_on_boundary( ImprintLineSegment *seg_1,
                                                        ImprintLineSegment *seg_2,
                                                        MatchType &type_0, MatchType &type_1,
                                                        MatchType &type_2, MatchType &type_3 )
{
    //Find out if the points are on the boundary.  Also use the matching informatin
    //previously calculated, since some of the nodes may match more closely with
    //other nodes on the boundary.  For instance, if the boundary does a quick
    //reversal, a point may match an existing point better than here.
    //If the points should be on the segment, set the matchtype flags accordingly.
  ImprintPointData* imp_point_0,*imp_point_1,*imp_point_2,*imp_point_3;
  ImprintPointData* imp_point_4,*imp_point_5,*imp_point_6,*imp_point_7;
  imp_point_0 = seg_1->get_start();
  imp_point_1 = seg_1->get_end();
  imp_point_2 = seg_2->get_start();
  imp_point_3 = seg_2->get_end();

  imp_point_4 = seg_1->get_prev()->get_start();
  imp_point_5 = seg_1->get_next()->get_end();
  imp_point_6 = seg_2->get_prev()->get_start();
  imp_point_7 = seg_2->get_next()->get_end();

  CubitBoolean point_0_on_seg_2 =CUBIT_FALSE;
  CubitBoolean point_1_on_seg_2 =CUBIT_FALSE;
  CubitBoolean point_2_on_seg_1 =CUBIT_FALSE;
  CubitBoolean point_3_on_seg_1 =CUBIT_FALSE;

  if ( on_interior_segment(imp_point_0, seg_2) )
    point_0_on_seg_2 = CUBIT_TRUE;
  if ( on_interior_segment(imp_point_1, seg_2) )
    point_1_on_seg_2 = CUBIT_TRUE;
  if ( on_interior_segment(imp_point_2, seg_1) )
    point_2_on_seg_1 = CUBIT_TRUE;
  if ( on_interior_segment(imp_point_3, seg_1) )
    point_3_on_seg_1 = CUBIT_TRUE;

  if ( point_0_on_seg_2 )
  {
    if ( type_0 == MATCH_0_6 )
    {
      double dist_1 = (imp_point_0->coordinates() -
                       seg_2->get_prev()->get_start()->coordinates()).length_squared();
      CubitVector point_on_seg_2;
      closest_point_interior_seg(imp_point_0, seg_2, point_on_seg_2);
      double dist_2 = (imp_point_0->coordinates() - point_on_seg_2).length_squared();
      if ( dist_1 < dist_2 )
        type_0 = NO_MATCH;
      else
        type_0 = ON_SEG_2;
    }
    else if (type_0 == MATCH_0_7 )
    {
      double dist_1 = (imp_point_0->coordinates() -
                       seg_2->get_next()->get_end()->coordinates()).length_squared();
      CubitVector point_on_seg_2;
      closest_point_interior_seg(imp_point_0, seg_2, point_on_seg_2);
      double dist_2 = (imp_point_0->coordinates() - point_on_seg_2).length_squared();
      if ( dist_1 < dist_2 )
        type_0 = NO_MATCH;
      else
        type_0 = ON_SEG_2;
    }
    else 
      type_0 = ON_SEG_2;
  }
    
  if ( point_1_on_seg_2 )
  {
    if ( type_1 == MATCH_1_6 )
    {
      double dist_1 = (imp_point_1->coordinates() -
                       seg_2->get_prev()->get_start()->coordinates()).length_squared();
      CubitVector point_on_seg_2;
      closest_point_interior_seg(imp_point_1, seg_2, point_on_seg_2);
      double dist_2 = (imp_point_1->coordinates() - point_on_seg_2).length_squared();
      if ( dist_1 < dist_2 )
        type_1 = NO_MATCH;
      else
        type_1 = ON_SEG_2;
    }
    else if (type_1 == MATCH_1_7 )
    {
      double dist_1 = (imp_point_1->coordinates() -
                       seg_2->get_next()->get_end()->coordinates()).length_squared();
      CubitVector point_on_seg_2;
      closest_point_interior_seg(imp_point_1, seg_2, point_on_seg_2);
      double dist_2 = (imp_point_1->coordinates() - point_on_seg_2).length_squared();
      if ( dist_1 < dist_2 )
        type_1 = NO_MATCH;
      else
        type_1 = ON_SEG_2;
    }
    else 
      type_1 = ON_SEG_2;
  }
  if ( point_2_on_seg_1 )
  {
    if ( type_2 == MATCH_2_4 )
    {
      double dist_1 = (imp_point_2->coordinates() -
                       seg_1->get_prev()->get_start()->coordinates()).length_squared();
      CubitVector point_on_seg_1;
      closest_point_interior_seg(imp_point_2, seg_1, point_on_seg_1);
      double dist_2 = (imp_point_2->coordinates() - point_on_seg_1).length_squared();
      if ( dist_1 < dist_2 )
        type_2 = NO_MATCH;
      else
        type_2 = ON_SEG_1;
    }
    else if (type_2 == MATCH_2_5 )
    {
      double dist_1 = (imp_point_2->coordinates() -
                       seg_1->get_next()->get_end()->coordinates()).length_squared();
      CubitVector point_on_seg_1;
      closest_point_interior_seg(imp_point_2, seg_1, point_on_seg_1);
      double dist_2 = (imp_point_2->coordinates() - point_on_seg_1).length_squared();
      if ( dist_1 < dist_2 )
        type_2 = NO_MATCH;
      else
        type_2 = ON_SEG_1;
    }
    else 
      type_2 = ON_SEG_1;
  }
  if ( point_3_on_seg_1 )
  {
    if ( type_3 == MATCH_3_4 )
    {
      double dist_1 = (imp_point_3->coordinates() -
                       seg_1->get_prev()->get_start()->coordinates()).length_squared();
      CubitVector point_on_seg_1;
      closest_point_interior_seg(imp_point_3, seg_1, point_on_seg_1);
      double dist_2 = (imp_point_3->coordinates() - point_on_seg_1).length_squared();
      if ( dist_1 < dist_2 )
        type_3 = NO_MATCH;
      else
        type_3 = ON_SEG_1;
    }
    else if (type_3 == MATCH_3_5 )
    {
      double dist_1 = (imp_point_3->coordinates() -
                       seg_1->get_next()->get_end()->coordinates()).length_squared();
      CubitVector point_on_seg_1;
      closest_point_interior_seg(imp_point_3, seg_1, point_on_seg_1);
      double dist_2 = (imp_point_3->coordinates() - point_on_seg_1).length_squared();
      if ( dist_1 < dist_2 )
        type_3 = NO_MATCH;
      else
        type_3 = ON_SEG_1;
    }
    else 
      type_3 = ON_SEG_1;
  }
  return CUBIT_SUCCESS;
}

CubitStatus ImprintBoundaryTool::case_1_3_equal( ImprintLineSegment *seg_1,
                                                 ImprintLineSegment *seg_2,
                                                 ImprintPointData* imp_point_0,
                                                 ImprintPointData* imp_point_1,
                                                 ImprintPointData* imp_point_2,
                                                 ImprintPointData* imp_point_3,
                                                 MatchType type_0, MatchType type_2,
                                                 IntersectResult &int_result,
                                                 ImprintLineSegment **new_segments)
{
    //Okay at this point we know that point_1 and point_3 are
    //equal, and we know that point_0 doesn't match an end_point.
    //We know then that this is either an L or a overlap-join.
    //Test to see if point_0 is on seg_2 or point_2 is on seg_1.
    //If neither then we have an L, other wise an overlap-join.

    //First test to see if they are are on the segments.
  CubitBoolean point_0_on_seg_2 =CUBIT_FALSE;
  CubitBoolean point_2_on_seg_1 =CUBIT_FALSE;
  if ( on_interior_segment(imp_point_0, seg_2) )
    point_0_on_seg_2 = CUBIT_TRUE;
  if ( on_interior_segment(imp_point_2, seg_1) )
    point_2_on_seg_1 = CUBIT_TRUE;
  
  if ( !point_0_on_seg_2 && !point_2_on_seg_1 )
  {
    PRINT_DEBUG_129("Found L_INTERSECT_1_3\n");
    int_result = L_INTERSECT_1_3;
    set_type_for_L(imp_point_1, imp_point_3);
    return CUBIT_SUCCESS;
  }
    //Now resolve if there is a closer match to one of the other
    //points.
  else if (point_0_on_seg_2 )
  {
    if (type_0 == MATCH_0_6 )
    {
      double dist_1 = (imp_point_0->coordinates() -
                       seg_2->get_prev()->get_start()->coordinates()).length_squared();
      CubitVector point_on_seg_2;
      closest_point_interior_seg(imp_point_0, seg_2, point_on_seg_2);
      double dist_2 = (imp_point_0->coordinates() - point_on_seg_2).length_squared();
      if ( dist_1 < dist_2 )
      {
        PRINT_DEBUG_129("Found L_INTERSECT_1_3\n");
        int_result = L_INTERSECT_1_3;
        set_type_for_L(imp_point_1, imp_point_3);
        return CUBIT_SUCCESS;
      }
    }
    else if (type_0 == MATCH_0_7 )
    {
      double dist_1 = (imp_point_0->coordinates() -
                       seg_2->get_next()->get_end()->coordinates()).length_squared();
      CubitVector point_on_seg_2;
      closest_point_interior_seg(imp_point_0, seg_2, point_on_seg_2);
      double dist_2 = (imp_point_0->coordinates() - point_on_seg_2).length_squared();
      if ( dist_1 < dist_2 )
      {
        PRINT_DEBUG_129("Found L_INTERSECT_1_3\n");
        int_result = L_INTERSECT_1_3;
        set_type_for_L(imp_point_1, imp_point_3);
        return CUBIT_SUCCESS;
      }
    }
  }
  else if (point_2_on_seg_1 )
  {
    if (type_2 == MATCH_2_4 )
    {
      double dist_1 = (imp_point_2->coordinates() -
                       seg_1->get_prev()->get_start()->coordinates()).length_squared();
      CubitVector point_on_seg_1;
      closest_point_interior_seg(imp_point_2, seg_1, point_on_seg_1);
      double dist_2 = (imp_point_2->coordinates() - point_on_seg_1).length_squared();
      if ( dist_1 < dist_2 )
      {
        PRINT_DEBUG_129("Found L_INTERSECT_1_3\n");
        int_result = L_INTERSECT_1_3;
        set_type_for_L(imp_point_1, imp_point_3);
        return CUBIT_SUCCESS;
      }
    }
    else if (type_2 == MATCH_2_5 )
    {
      double dist_1 = (imp_point_2->coordinates() -
                       seg_1->get_next()->get_end()->coordinates()).length_squared();
      CubitVector point_on_seg_1;
      closest_point_interior_seg(imp_point_2, seg_1, point_on_seg_1);
      double dist_2 = (imp_point_2->coordinates() - point_on_seg_1).length_squared();
      if ( dist_1 < dist_2 )
      {
        PRINT_DEBUG_129("Found L_INTERSECT_1_3\n");
        int_result = L_INTERSECT_1_3;
        set_type_for_L(imp_point_1, imp_point_3);
        return CUBIT_SUCCESS;
      }
    }
  }

  if ( point_0_on_seg_2 )
  {
      //This is an overlap_join.
    PRINT_DEBUG_129("Found OVERLAP_JOIN_13_0_2\n");
    int_result = OVERLAP_JOIN_13_0_2;
      //Create the new segments, the first segment remains the same so
      //leave positions 0 and 1 empty.
    ImprintLineSegment *new_seg_2_0, *new_seg_0_3;
      //copy point 0 to use on segment 2
    CubitVector closest_point;
    closest_point_seg( imp_point_0, seg_2, closest_point);
    ImprintPointData *new_point_0 = new ImprintPointData(imp_point_0, closest_point);
    new_point_0->owner(seg_2->owner());
    allocatedPointData->append(new_point_0);
    if ( imp_point_0->is_owner_vertex() )
    {
      new_point_0->set_point_type(CREATE_NEW_VERTEX);
      if ( imp_point_0->get_point_type() != CREATE_NEW_VERTEX )
        imp_point_0->set_point_type(VERTEX_ON_BOTH_BOUNDARIES);
    }
    else
    {
      new_point_0->set_point_type(ON_BOTH_BOUNDARIES);
      imp_point_0->set_point_type(ON_BOTH_BOUNDARIES);
    }
    new_seg_2_0 = new ImprintLineSegment(imp_point_2, new_point_0,
                                         seg_2->owner());
    new_seg_0_3 = new ImprintLineSegment(new_point_0, imp_point_3,
                                         seg_2->owner());
    new_segments[2] = new_seg_2_0;
    new_segments[3] = new_seg_0_3;
    allocatedLineData->append(new_segments[2]);
    allocatedLineData->append(new_segments[3]);
      //update the list.
    update_list(new_segments, seg_1, seg_2, CUBIT_FALSE);
    seg_2->set_inactive(CUBIT_TRUE);
    return CUBIT_SUCCESS;
  }
  else if ( point_2_on_seg_1 )
  {
      //This is an overlap_join.
    PRINT_DEBUG_129("Found OVERLAP_JOIN_13_2_0\n");
    int_result = OVERLAP_JOIN_13_2_0;
      //Create the new segments, the first segment remains the same so
      //leave positions 0 and 1 empty.
    ImprintLineSegment *new_seg_0_2, *new_seg_2_1;
      //copy point 2 to use on segment 1
    CubitVector closest_point;
    closest_point_seg( imp_point_2, seg_1, closest_point);
    ImprintPointData *new_point_2 = new ImprintPointData(imp_point_2, closest_point);
    new_point_2->owner(seg_1->owner());
    allocatedPointData->append(new_point_2);
    if ( imp_point_2->is_owner_vertex() )
    {
      new_point_2->set_point_type(CREATE_NEW_VERTEX);
      if ( imp_point_2->get_point_type() != CREATE_NEW_VERTEX )
        imp_point_2->set_point_type(VERTEX_ON_BOTH_BOUNDARIES);
    }
    else
    {
      new_point_2->set_point_type(ON_BOTH_BOUNDARIES);
      imp_point_2->set_point_type(ON_BOTH_BOUNDARIES);
    }
    new_seg_0_2 = new ImprintLineSegment(imp_point_0, new_point_2,
                                         seg_1->owner());
    new_seg_2_1 = new ImprintLineSegment(new_point_2, imp_point_1,
                                         seg_1->owner());
    new_segments[0] = new_seg_0_2;
    new_segments[1] = new_seg_2_1;
    allocatedLineData->append(new_segments[0]);
    allocatedLineData->append(new_segments[1]);
      //update the list.
    update_list(new_segments, seg_1, seg_2, CUBIT_TRUE);
    seg_1->set_inactive(CUBIT_TRUE);
    return CUBIT_SUCCESS;
  }
  else
  {
      //For this case we must have an L condition.
    set_type_for_L(imp_point_1, imp_point_3);
    PRINT_DEBUG_129("Found L_INTERSECT_1_3\n");
    int_result = L_INTERSECT_1_3;
    return CUBIT_SUCCESS;
  }
}
CubitStatus ImprintBoundaryTool::set_type_for_equal(ImprintPointData *pair_1_1,
                                                    ImprintPointData *pair_1_2,
                                                    ImprintPointData *pair_2_1,
                                                    ImprintPointData *pair_2_2)
{
    //These segments are equal.
    //Set the point types accoridingly.
  if ( pair_1_1->is_owner_vertex() &&
       !pair_1_2->is_owner_vertex() )
  {
    pair_1_2->set_point_type(CREATE_NEW_VERTEX);
    if ( pair_1_1->get_point_type() != CREATE_NEW_VERTEX )
      pair_1_1->set_point_type(VERTEX_ON_BOTH_BOUNDARIES);
  }
  else if (!pair_1_1->is_owner_vertex() &&
           pair_1_2->is_owner_vertex() )
  {
    pair_1_1->set_point_type(CREATE_NEW_VERTEX);
    if ( pair_1_2->get_point_type() != CREATE_NEW_VERTEX )
      pair_1_2->set_point_type(VERTEX_ON_BOTH_BOUNDARIES);
  }
  else if ( pair_1_1->is_owner_vertex() &&
            pair_1_2->is_owner_vertex() )
  {
    if ( pair_1_1->get_point_type() != CREATE_NEW_VERTEX )
      pair_1_1->set_point_type(VERTEX_ON_BOTH_BOUNDARIES);
    if (pair_1_2->get_point_type() != CREATE_NEW_VERTEX )
      pair_1_2->set_point_type(VERTEX_ON_BOTH_BOUNDARIES);
  }
  else
  {
    pair_1_1->set_point_type(ON_BOTH_BOUNDARIES);
    pair_1_2->set_point_type(ON_BOTH_BOUNDARIES);
  }
  if ( pair_2_1->is_owner_vertex() &&
       !pair_2_2->is_owner_vertex() )
  {
    pair_2_2->set_point_type(CREATE_NEW_VERTEX);
    if ( pair_2_1->get_point_type() != CREATE_NEW_VERTEX )
      pair_2_1->set_point_type(VERTEX_ON_BOTH_BOUNDARIES);
  }
  else if (!pair_2_1->is_owner_vertex() &&
           pair_2_2->is_owner_vertex() )
  {
    pair_2_1->set_point_type(CREATE_NEW_VERTEX);
    if ( pair_2_2->get_point_type() != CREATE_NEW_VERTEX )
      pair_2_2->set_point_type(VERTEX_ON_BOTH_BOUNDARIES);
  }
  else if ( pair_2_1->is_owner_vertex() &&
            pair_2_2->is_owner_vertex() )
  {
    if ( pair_2_1->get_point_type() != CREATE_NEW_VERTEX )
      pair_2_1->set_point_type(VERTEX_ON_BOTH_BOUNDARIES);
    if ( pair_2_2->get_point_type() != CREATE_NEW_VERTEX )
      pair_2_2->set_point_type(VERTEX_ON_BOTH_BOUNDARIES);
  }
  else
  {
    pair_2_1->set_point_type(ON_BOTH_BOUNDARIES);
    pair_2_2->set_point_type(ON_BOTH_BOUNDARIES);
  }
  return CUBIT_SUCCESS;
}
CubitStatus ImprintBoundaryTool::set_type_for_L( ImprintPointData *imp_point_1,
                                                 ImprintPointData *imp_point_2)
{
  CubitBoolean owner_vert_1 = imp_point_1->is_owner_vertex();
  CubitBoolean owner_vert_2 = imp_point_2->is_owner_vertex();
  if ( owner_vert_1 && !owner_vert_2 )
  {
    imp_point_2->set_point_type(CREATE_NEW_VERTEX);
    if (imp_point_1->get_point_type() != CREATE_NEW_VERTEX )
      imp_point_1->set_point_type(VERTEX_ON_BOTH_BOUNDARIES);
  }
  else if (!owner_vert_1 && owner_vert_2 )
  {
    imp_point_1->set_point_type(CREATE_NEW_VERTEX);
    if (imp_point_2->get_point_type() != CREATE_NEW_VERTEX )
      imp_point_2->set_point_type(VERTEX_ON_BOTH_BOUNDARIES);
  }
  else if ( owner_vert_1 && owner_vert_2 )
  {
      //Only make them verts on both boundaries if
      //they aready all vertices.  Their owner could
      //have been set previously by some other intersection
      //like a cross or a T.  It is important that points
      //that need vertices on them end up being marked
      //as a vertex on the other boundary...
    if(imp_point_1->get_point_type() != CREATE_NEW_VERTEX )
      imp_point_1->set_point_type(VERTEX_ON_BOTH_BOUNDARIES);
    if(imp_point_2->get_point_type() != CREATE_NEW_VERTEX )
      imp_point_2->set_point_type(VERTEX_ON_BOTH_BOUNDARIES);
  }
  else
  {
    imp_point_1->set_point_type(ON_BOTH_BOUNDARIES);
    imp_point_2->set_point_type(ON_BOTH_BOUNDARIES);
  }
  return CUBIT_SUCCESS;
}
CubitBoolean ImprintBoundaryTool::on_interior_segment(ImprintPointData *point,
                                                      ImprintLineSegment *line )
{
    //Find the closest point on the line segement.
  CubitVector closest_point;
  CubitBoolean on_interior = closest_point_interior_seg( point, line, closest_point);
  if (!on_interior )
    return CUBIT_FALSE;
  CubitVector my_point = point->coordinates();
  if ( my_point.within_tolerance( closest_point, myTolerance ) )
    return CUBIT_TRUE;
  else
    return CUBIT_FALSE;
}

CubitBoolean ImprintBoundaryTool::closest_point_interior_seg( ImprintPointData *point,
                                                              ImprintLineSegment *line,
                                                              CubitVector &closest_point )
{
    //Find the closest point between this line and the point.
  double node[3];
  double point_0[3], point_1[3];
  CubitVector point_v = point->coordinates();
  node[0] = point_v.x();
  node[1] = point_v.y();
  node[2] = point_v.z();
  CubitVector start_v = line->get_start()->coordinates();
  CubitVector end_v = line->get_end()->coordinates();
  point_0[0] = start_v.x();
  point_0[1] = start_v.y();
  point_0[2] = start_v.z();
  point_1[0] = end_v.x();
  point_1[1] = end_v.y();
  point_1[2] = end_v.z();
    //use geometry_resabs rather than tolerance to get more accuracy.
  IntersectionTool new_tool(GEOMETRY_RESABS);
  double t;
  double dist = new_tool.distance_point_line(node, point_0,
                                             point_1, t);
  if ( dist < 0.0 )
  {
    return CUBIT_FALSE;
  }
  else
  {
    CubitVector p_0_1 = end_v - start_v;
    closest_point = start_v + t*p_0_1;
  }
  return CUBIT_TRUE;
}
CubitStatus ImprintBoundaryTool::closest_point_seg( ImprintPointData *point,
                                                    ImprintLineSegment *line,
                                                    CubitVector &closest_point )
{
    //Find the closest point between this line and the point.
  double node[3];
  double point_0[3], point_1[3];
  CubitVector point_v = point->coordinates();
  node[0] = point_v.x();
  node[1] = point_v.y();
  node[2] = point_v.z();
  CubitVector start_v = line->get_start()->coordinates();
  CubitVector end_v = line->get_end()->coordinates();
  point_0[0] = start_v.x();
  point_0[1] = start_v.y();
  point_0[2] = start_v.z();
  point_1[0] = end_v.x();
  point_1[1] = end_v.y();
  point_1[2] = end_v.z();
    //use geometry_resabs rather than tolerance to get more accuracy.
  IntersectionTool new_tool(GEOMETRY_RESABS);
  double t;
  double dist = new_tool.distance_point_line(node, point_0,
                                             point_1, t);
  if ( dist < 0.0 )
  {
      //This happens if the point is beyond the line segment, in which
      //case we just want the closer end point.
    if ( t < 0.0 )
    {
      closest_point = start_v;
    }
    else
    {
      closest_point = end_v;
    }
  }
  else
  {
    CubitVector p_0_1 = end_v - start_v;
    closest_point = start_v + t*p_0_1;
  }
  return CUBIT_SUCCESS;
}
CubitStatus ImprintBoundaryTool::closest_point_seg( CubitVector &point_v,
                                                    ImprintLineSegment *line,
                                                    CubitVector &closest_point )
{
    //Find the closest point between this line and the point.
  double node[3];
  double point_0[3], point_1[3];
  node[0] = point_v.x();
  node[1] = point_v.y();
  node[2] = point_v.z();
  CubitVector start_v = line->get_start()->coordinates();
  CubitVector end_v = line->get_end()->coordinates();
  point_0[0] = start_v.x();
  point_0[1] = start_v.y();
  point_0[2] = start_v.z();
  point_1[0] = end_v.x();
  point_1[1] = end_v.y();
  point_1[2] = end_v.z();
    //use geometry_resabs rather than tolerance to get more accuracy.
  IntersectionTool new_tool(GEOMETRY_RESABS);
  double t;
  double dist = new_tool.distance_point_line(node, point_0,
                                             point_1, t);
  if ( dist < 0.0 )
  {
      //This happens if the point is beyond the line segment, in which
      //case we just want the closer end point.
    if ( t < 0.0 )
    {
      closest_point = start_v;
    }
    else
    {
      closest_point = end_v;
    }
  }
  else
  {
    CubitVector p_0_1 = end_v - start_v;
    closest_point = start_v + t*p_0_1;
  }
  return CUBIT_SUCCESS;
}
                                                    
CubitStatus ImprintBoundaryTool::convert_to_lines( PointLoopList &boundary_point_loops,
                                                   SegLoopList &boundary_line_loops,
                                                   RefFace *ref_face,
                                                   const CubitBoolean boundary_1 )
{
    //This does more than just converts the boundary points to line segments.  It also
    //Sets up a doubily linked list through the ImprintLineSegment data structure.
    //Calling codes need to make sure all this data is deleted.
    //This function also sets the PointType for each of the boundary points.  This
    //will depend also on the boundary_1 flag.  If we are doing boundary_1 then,
    //the flag will be true, otherwise it will be false.
  int ii, jj;
  PointList *point_list;
  SegList *segment_list;
  ImprintLineSegment *new_line, *prev = NULL, *next = NULL;
  RefEntity *entity_start, *entity_end, *seg_owner;
  boundary_point_loops.reset();
  PointType vertex_on_boundary, on_boundary;
  if ( boundary_1 )
  {
    vertex_on_boundary = VERTEX_ON_BOUNDARY_1;
    on_boundary = ON_BOUNDARY_1;
  }
  else
  {
    vertex_on_boundary = VERTEX_ON_BOUNDARY_2;
    on_boundary = ON_BOUNDARY_2;
  }
  for ( ii = boundary_point_loops.size(); ii > 0; ii-- )
  {
    point_list = boundary_point_loops.get_and_step();
    segment_list = new SegList;
    prev = NULL;
    next = NULL;
    new_line = NULL;
    seg_owner = NULL;
    for ( jj = point_list->size(); jj > 0; jj-- )
    {
      entity_start = point_list->get()->owner();
      entity_end = point_list->next()->owner();
      seg_owner =  CAST_TO(entity_start, RefEdge);
      if ( seg_owner == NULL )
      {
        seg_owner = CAST_TO(entity_end, RefEdge);
        if ( seg_owner == NULL )
        {
          RefVertex *ref_vert_end = CAST_TO(entity_end, RefVertex);
          RefVertex *ref_vert_start = CAST_TO(entity_start, RefVertex);
          DLIList <RefEdge*> common_edges;
          ref_vert_start->common_ref_edges(ref_vert_end,
                                           common_edges,
                                           ref_face);
          if ( common_edges.size() == 0  )
          {
            PRINT_ERROR("ImprintBoundaryTool::convert_to_lines having problem finding owner of boundary\n"
                        "segment.\n");
            return CUBIT_FAILURE;
          }
          if ( common_edges.size() == 1 )
            seg_owner = common_edges.get();
          else
          {
              //Now we have to decide which of these edges is the one.
              //Lets take a mid point on this segment and check the distances.
              //First find the mid_point of the line segment.
            CubitVector start = point_list->get()->coordinates();
            CubitVector end = point_list->next()->coordinates();
            CubitVector mid = (start+end)/2.0;
            CubitVector closest_point;
            
              //Now test the edges to see which is closest.
            RefEdge *closest_edge = NULL, *curr_ref_edge;
            double min_dist = CUBIT_DBL_MAX, dist;
            int kk;
            for ( kk = common_edges.size(); kk > 0; kk-- )
            {
              curr_ref_edge = common_edges.get_and_step();
              curr_ref_edge->closest_point(mid, closest_point);
              dist = (mid - closest_point).length_squared();
              if ( dist < min_dist )
              {
                min_dist = dist;
                closest_edge = curr_ref_edge;
              }
            }
            if ( closest_edge == NULL )
            {
              PRINT_ERROR("Problems determining segment ownwership.\n"
                          "Internal problem with virtual imprinting.\n");
              boundary_line_loops.clean_out();
			  return CUBIT_FAILURE;
            }
            seg_owner = closest_edge;
          }
        }
      }
      new_line = new ImprintLineSegment( point_list->get(), point_list->next(), seg_owner);
      allocatedLineData->append(new_line);
        //Now for use later on, set the point type.  This will be useful
        //when we get to the intersection parts.
      if ( CAST_TO(entity_start, RefEdge) != NULL )
      {
        new_line->get_start()->set_point_type(on_boundary);
      }
      else {
        assert( CAST_TO(entity_start, RefVertex) != NULL );
        new_line->get_start()->set_point_type(vertex_on_boundary);
      }
      if ( CAST_TO(entity_end, RefEdge) != NULL )
      {
        new_line->get_end()->set_point_type(on_boundary);
      }
      else {
        assert( CAST_TO(entity_end, RefVertex) != NULL );
        new_line->get_end()->set_point_type(vertex_on_boundary);
      }
      point_list->step();
      segment_list->append(new_line);
      if ( prev != NULL )
      {
        prev->set_next(new_line);
        new_line->set_prev(prev);
      }
      prev = new_line;
    }
    if ( prev != NULL )
    {
      segment_list->reset();
      assert(prev == segment_list->prev());
      segment_list->prev()->set_next(segment_list->get());
      segment_list->get()->set_prev(segment_list->prev());
      assert(segment_list->get()->get_next() != NULL );
    }
    if ( segment_list->size() )
    {
      boundary_line_loops.append(segment_list);
      allocatedLineLoops->append(segment_list);
    }
    else
      delete segment_list;
  }
  return CUBIT_SUCCESS;
}
//--------------------------------------------------------------
// Private Function: find_graph_for_surf
// Description:  Given the boundary_line_loops for a surface
//               that have already been intersected against another
//               surface, and the boundary_line_loops for that surface,
//               determine the line segments that are needed
//               to perform the actual imprint.  Also find points
//               That are need to partition or imprint the existing boundaries.
//               These points will also include the points of the segments
//               that touch against the boundary, if necessary.  The
//               RefFace is pretty important as it is used to determine
//               if the points are inside or outside the surface.
//---------------------------------------------------------------
CubitStatus ImprintBoundaryTool::find_graph_for_surf(PointLoopList &boundary_loops_1,
                                                     PointLoopList &boundary_loops_2,
                                                     RefFace *ref_face,
                                                     PointLoopList &part_segs,
                                                     PointList &partition_points,
                                                     CubitBoolean surf_1)
                                                     
{
    //Traverse the boundary_line_loops_2.
    //First search for intersections. If an intersection is found,
    //traverse from there until another intersection is found.
    //If no intersection is found for a loop, test to see if the
    //nodes lie inside the surface.  If they do, then the intire loop
    //should be used for partitioning.
  int ii, jj;
    //ASSUME my boundary is loop_2.
  CubitBoolean start_recording = CUBIT_FALSE;
  CubitBoolean found_start = CUBIT_FALSE;
  
  PointList *loop_1_ptr, *loop_2_ptr, *new_part_line;
  ImprintPointData *curr_point, *next_point;
  for ( ii = boundary_loops_2.size(); ii > 0; ii-- )
  {
    loop_2_ptr = boundary_loops_2.get_and_step();
    start_recording = CUBIT_FALSE;
    new_part_line = NULL;
    found_start = CUBIT_FALSE;
      //find a good place to start.  It messes things
      //up if we start the loop at a non-meeting spot.
    for ( jj = loop_2_ptr->size(); jj > 0; jj-- )
    {
      curr_point = loop_2_ptr->get_and_step();
      if ( curr_point->get_point_type() == CREATE_NEW_VERTEX ||
           curr_point->get_point_type() == VERTEX_ON_BOTH_BOUNDARIES )
      {
        loop_2_ptr->back();
        found_start = CUBIT_TRUE;
        break;
      }
    }
    if ( found_start )
    {
      for ( jj = loop_2_ptr->size(); jj > 0; jj-- )
      {
        curr_point = loop_2_ptr->get_and_step();
        next_point = loop_2_ptr->get();
        if ( DEBUG_FLAG(129) )
          draw_point(curr_point);
        if ( curr_point->get_point_type() == CREATE_NEW_VERTEX ||
             curr_point->get_point_type() == VERTEX_ON_BOTH_BOUNDARIES )
        {
          CubitStatus stat;
            //This section of the code just got too big.  So it may change
            //the loop_2_ptr position, the loop counter jj and other things.
          stat = point_intersects_case( curr_point, next_point,
                                        ref_face,
                                        part_segs,  new_part_line,
                                        start_recording, surf_1 );
        }
        else
        {
          if (start_recording)
          {
            if ( !new_part_line )
            {
              PRINT_ERROR("Bad Logic ImprintBoundaryTool::find_graph_for_surf.\n");
              assert(new_part_line != NULL);
              return CUBIT_FAILURE;
            }
            CubitBox surf_box = ref_face->bounding_box();
            CubitBox curr_box(curr_point->coordinates());
            if( surf_box.overlap(myTolerance, curr_box ) &&
                on_surface(curr_point, ref_face) )
            {
              new_part_line->append(curr_point);
              if ( jj == 1 && 
                   (next_point->get_point_type() == CREATE_NEW_VERTEX ||
                    next_point->get_point_type() == VERTEX_ON_BOTH_BOUNDARIES) )
              {
                  //Do the last node since that is where we have to go to inorder
                  //to stop the partitioning.
                jj++;
              }
            }
            else
            {
                //remove new_part_line from the part_segs list and
                //break out of this loop.
              part_segs.pop();
              start_recording = CUBIT_FALSE;
              break;
            }
          }
        }
      }
      if ( start_recording )
      {
          //This is bad news.  We never should finish the loop and still
          //be recording...
        PRINT_ERROR("Bad Logic ImprintBoundaryTool::find_graph_for_surf.\n");
        assert(!start_recording);
        return CUBIT_FAILURE;
      }
    }
    else
    {
        //Okay this is an interior hole to the ref_face.  Just add the
        //entire loop.
      int mm;
      CubitBox surf_box = ref_face->bounding_box();		
      int debug = 0;
      if ( debug )
      {
        GfxDebug::clear();
        GfxDebug::draw_box(surf_box, CUBIT_RED);
        GfxDebug::flush();
      }
      CubitBox curr_box;
        //Getting the first point is important.  Make
        //sure it is a vertex.
      CubitBoolean start = CUBIT_FALSE;
      for ( mm = loop_2_ptr->size(); mm > 0; mm-- )
      {
        curr_point = loop_2_ptr->get_and_step();
        curr_box.reset(curr_point->coordinates());
        if ( debug )
        {
          draw_point(curr_point);
        }
        if ( surf_box.overlap(myTolerance,curr_box) )
        {
          if ( on_surface(curr_point, ref_face) &&
               curr_point->is_owner_vertex() )
          {
            loop_2_ptr->back();
            start = CUBIT_TRUE;
            break;
          }
        }
      }
      if ( !start )
      {
          //this loop has no affect on the surface. (no points
          //are on the boundary and are vertices...)
        continue;
      }
      else
      {
        new_part_line = new PointList;
        allocatedPointLoops->append(new_part_line);
          //Test one the nodes to see if they are interior to the surface.
        loop_2_ptr->get()->set_start_partition();
        loop_2_ptr->get()->set_end_partition();
        part_segs.append(new_part_line);
        CubitBoolean valid_loop = CUBIT_TRUE;
        for ( mm = loop_2_ptr->size(); mm > 0; mm-- )
        {
          curr_point = loop_2_ptr->get_and_step();
          curr_box.reset(curr_point->coordinates());
          if ( surf_box.overlap(myTolerance, curr_box) )
          {
            if ( on_surface(curr_point, ref_face) )
            {
              new_part_line->append(curr_point);
            }
            else
              valid_loop = CUBIT_FALSE;
          }
          else
            valid_loop = CUBIT_FALSE;
        }
        if ( !valid_loop )
        {
            //remove the part_line.  The
            //entire loop must go into the surface.
          part_segs.pop();
        }
      }
        //Otherwise just ignore this loop it has no affect on the ref_face.
    }
  }
    //Okay now go through and find the points that are needed to create as imprints on the boundary.
  for ( ii = boundary_loops_1.size(); ii > 0; ii-- )
  {
    loop_1_ptr = boundary_loops_1.get_and_step();
    for ( jj = loop_1_ptr->size(); jj > 0; jj-- )
    {
      curr_point = loop_1_ptr->get_and_step();
      if ( curr_point->get_point_type() == CREATE_NEW_VERTEX )
        partition_points.append(curr_point);
    }
  }
    //Okay, we should be set now!
  return CUBIT_SUCCESS;
}
//--------------------------------------------------------------
// Private Function: point_intersects_case
// Description:  Given the fact that curr_point's PointType shows
//               that it is a vertex on both boundaries or a new vertex,
//               Resolve how the partition curves are extraced from the
//               boundary loop 2 on surface 1.
//--------------------------------------------------------------
CubitStatus ImprintBoundaryTool::point_intersects_case( ImprintPointData *curr_point,
                                                        ImprintPointData *next_point,
                                                        RefFace *ref_face,
                                                        PointLoopList &part_segs,
                                                        PointList *&new_part_line,
                                                        CubitBoolean &start_recording,
                                                        CubitBoolean surf_1 )
{
  PointType vertex_my_boundary = surf_1 ? VERTEX_ON_BOUNDARY_2 : VERTEX_ON_BOUNDARY_1;
  PointType edge_my_boundary = surf_1 ? ON_BOUNDARY_2 : ON_BOUNDARY_1;
  CubitBoolean do_nothing = CUBIT_FALSE;
  ImprintPointData *match_point_cur, *match_point_nex;
  if ( !start_recording )
  {
      //Okay.  This is an intersection.  The possible outcomes are:
      //1) This is the end of a paritition :-> do nothing.
      //2) This is the begining of a short partition (one segment)
      //3) This is the begining of a parition :-> start recording it.
      //4) This is a vertex imprint from our surface. :-> do nothing.
    
      //For case 1 resolution, we need to test to see if the next node is either:
      // a) on the boundary of both surf 1 and surf 2.
      // b) outside surf 1 and not intersecting surf_1's boundary.
      // If case a, then we need to find out if its matching point
      // and the matching point of the next point are on the same
      // segment.
    
      //For case 2, we need to test the underlying refentities to
      //see if we really are splitting a surface or on a the boundary.
    
      //For case 3, we need to test the next node to see if it is
      //only on loop 2's boundary and INTERIOR to the ref_face.
    
      //For case 4, that will be the counter examples of case 3 or
      //rather when their interior or refentity traversals fail.
      //Also if the next node is on both boundaries but not as
      //a vertex.  This is a simple overlap.
    
      //Both points are vertices.  Need to find out if this is a
      //single segment split or on the boundary.
    
    if ( next_point->get_point_type() == CREATE_NEW_VERTEX ||
         next_point->get_point_type() == VERTEX_ON_BOTH_BOUNDARIES )
    {
        //To decifer, determine this from the topology of the points
        //these match with on boundary 1.
      match_point_cur = curr_point->get_matching_point();
      match_point_nex = next_point->get_matching_point();
      int list_id_1 = match_point_cur->get_loop_pos();
      int list_id_2 = match_point_nex->get_loop_pos();
      if ( list_id_1 == 0 )
      {
        if (list_id_2 == match_point_nex->get_loop_size() - 1 )
        {
            //okay force this to work, we just went around the loop.
          list_id_2 = 1;
        }
      }
      if ( list_id_2 == 0 )
      {
        if (list_id_1 == match_point_nex->get_loop_size() - 1 )
        {
            //okay force this to work, we just went around the loop.
          list_id_1 = 1;
        }
      }
      int diff;
      if ( list_id_1 > list_id_2 )
        diff = list_id_1 - list_id_2;
      else
        diff = list_id_2 - list_id_1;
      if ( (match_point_cur->get_list_loop_pos() ==
            match_point_nex->get_list_loop_pos()) &&
           diff == 1 )
      {
          //CASE 1
          //This is the case that curr_point and next_point are
          //both create vertices or matched with vertices.
          //If the points that they match with, are also right
          //next to each other, then don't partition them since
          //they all match prefectly.
        do_nothing = CUBIT_TRUE;
      }
      else
      {
          //Also, this may be a case where the one boundary is
          //curved and the segment croses the curve portion.
          //We'll have to do a surface check mid-way through
          //the segment.
        CubitVector p0 = curr_point->coordinates();
        CubitVector p1 = next_point->coordinates();
        CubitVector point_q = p0 + (p1-p0)*.5;
        if ( on_surface(point_q, ref_face) )
        {
            //CASE 2
            //So they aren't connected.  This means that
            //we need to have a single segment where we
            //partition with just these two points.
          new_part_line = new PointList;
          allocatedPointLoops->append(new_part_line);
          new_part_line->append(curr_point);
          new_part_line->append(next_point);
          curr_point->set_start_partition();
          next_point->set_end_partition();
          part_segs.append(new_part_line);
          new_part_line = NULL;
          do_nothing = CUBIT_TRUE;
        }
        else
          do_nothing = CUBIT_TRUE;
      }
    }
    else if ( next_point->get_point_type() == ON_BOTH_BOUNDARIES )
    {
        //CASE 4.
        //hmmm.  This has got to be a simple imprint so do nothing here.
      do_nothing = CUBIT_TRUE;
    }
    else if ( next_point->get_point_type() == edge_my_boundary ||
              next_point->get_point_type() == vertex_my_boundary )
    {
        //We need to test to see if this point is on the ref_face or not.
        //If it is, then we need to start recording.  If it
        //isn't then we can do nothing.
      if ( on_surface(next_point, ref_face) )
          //CASE 3
        start_recording = CUBIT_TRUE;
      else
          //CASE 4
        do_nothing = CUBIT_TRUE;
    }
    else
    {
      PRINT_ERROR("Bad Logic ImprintBoundaryTool::find_graph_for_surf.\n");
      assert(0);
      return CUBIT_FAILURE;
    }
    if ( do_nothing )
      return CUBIT_SUCCESS;

    assert(start_recording);
      //Okay start recording.  Create a new list.  Put the curr_point
      //into that list and continue.
    curr_point->set_start_partition();
    new_part_line = new PointList;
    allocatedPointLoops->append(new_part_line);
    new_part_line->append(curr_point);
    part_segs.append(new_part_line);
  }
  else {
      //Now we need to terminate the recording.
    curr_point->set_end_partition();
    if ( !curr_point->get_start_partition() )
    {
      new_part_line->append(curr_point);
    }
    else{
        //case where circle loop imprints an edge at just one point.
        //do nothing.
    }
    new_part_line = NULL;
    start_recording = CUBIT_FALSE;
  }
  return CUBIT_SUCCESS;
}
//-------------------------------------------------------
// Private Function: on_surface
// Description:  Determines through the find_closest_point_trimmed
//               function if a point is inside or outside a surface.
//               This is not a very efficient method but since the
//               surface can be non-planar, we can't do a simple test.
//-------------------------------------------------------
CubitBoolean ImprintBoundaryTool::on_surface( ImprintPointData *point,
                                              RefFace *ref_face )
{
    //Determine if the vertex is on the curve within the feature size...
  CubitVector vert_point = point->coordinates();
  CubitVector point_on = vert_point;
    //Move the point to the curve and check the distance.
  ref_face->find_closest_point_trimmed(vert_point, point_on);
  if ( vert_point.within_tolerance(point_on, myTolerance) )
    return CUBIT_TRUE;
  else
    return CUBIT_FALSE;
}
CubitBoolean ImprintBoundaryTool::on_surface( CubitVector &vert_point,
                                              RefFace *ref_face )
{
    //Determine if the vertex is on the curve within the feature size...
  CubitVector point_on = vert_point;
    //Move the point to the curve and check the distance.
  ref_face->find_closest_point_trimmed(vert_point, point_on);
  if ( vert_point.within_tolerance(point_on, myTolerance) )
    return CUBIT_TRUE;
  else
    return CUBIT_FALSE;
}
CubitBoolean ImprintBoundaryTool::on_curve( ImprintPointData *point,
                                            RefEdge *ref_edge )
{
  CubitBoolean on_boundary = CUBIT_FALSE;
    //Determine if the vertex is on the curve within the feature size...
  CubitVector vert_point = point->coordinates();
  CubitVector point_on = point->coordinates();
    //Move the point to the curve and check the distance.
  ref_edge->closest_point_trimmed(vert_point, point_on);
  CubitVector distance = vert_point - point_on;
  double leng = distance.length();
  if ( leng < myTolerance )
  {
    on_boundary = CUBIT_TRUE;
  }
  return on_boundary;
}
//-------------------------------------------------------
// Private Function:are_connected
// Description:  Determines through topology traversals
//               if the two refentities are directly
//               connected.  The function assumes
//               that the RefEntities are either RefVerticies or
//               RefEdges.  The entites must be connected within
//               one entity for vertices or the same for
//               curves...
//-------------------------------------------------------
CubitBoolean ImprintBoundaryTool::are_connected( RefEntity *ent_1,
                                                 RefEntity *ent_2,
                                                 RefFace *ref_face )
{
  RefEdge *curve_cur, *curve_nex;
  RefVertex *vert_cur, *vert_nex;
  curve_cur = CAST_TO(ent_1, RefEdge);
  curve_nex = CAST_TO(ent_2, RefEdge);
  vert_cur = CAST_TO(ent_1, RefVertex);
  vert_nex = CAST_TO(ent_2, RefVertex);
    //just do some sanity checking first.  Really shouldn't happen...
  if ( !curve_cur && !vert_cur )
  {
    assert(curve_cur || vert_cur);
    PRINT_ERROR("messed up logic!!!! in ImprintBoundaryTool::are_connected\n");
    return CUBIT_FALSE;
  }
  if ( !curve_nex && !vert_nex )
  {
    assert(curve_nex || vert_nex);
    PRINT_ERROR("messed up logic!!!! in ImprintBoundaryTool::are_connected\n");
    return CUBIT_FALSE;
  }
  if ( curve_cur )
  {
    if ( curve_cur == curve_nex )
      return CUBIT_TRUE;
    if ( !vert_nex && curve_nex )
      return CUBIT_FALSE;
    if ( vert_nex && !curve_nex )
    {
      if ( curve_cur->start_vertex() == vert_nex ||
           curve_cur->end_vertex() == vert_nex )
        return CUBIT_TRUE;
      else
        return CUBIT_FALSE;
    }
  }
  else if ( vert_cur )
  {
    if ( vert_nex == vert_cur )
      return CUBIT_TRUE;
    if ( !curve_nex && vert_nex )
    {
      if ( vert_nex == vert_cur )
      {
        assert(0);
        PRINT_ERROR("Next and cur have same vertex in\n"
                    "ImprintBoundaryTool::are_connected.\n");
        return CUBIT_TRUE;
      }
      else if ( vert_cur->common_ref_edge(vert_nex, ref_face) )
        return CUBIT_TRUE;
      else
        return CUBIT_FALSE;
    }
    if ( !vert_nex && curve_nex )
    {
      if ( curve_nex->start_vertex() == vert_cur ||
           curve_nex->end_vertex() == vert_cur )
        return CUBIT_TRUE;
      else
        return CUBIT_FALSE;
    }
  }
  assert(0);
  PRINT_ERROR("ImprintBoundaryTool::are_connected failed\n");
  return CUBIT_FALSE;
}
CubitStatus ImprintBoundaryTool::imprint_boundary_vertices( PointList &partition_points,
                                                            CubitBoolean &modified_boundary )
{
  ImprintPointData *curr_point;
  RefEntity *ref_ent;
  int ii;
  RefEdge *ref_edge, *edge1, *edge2;  
  modified_boundary = CUBIT_FALSE;
  while( partition_points.size() > 0 )
  {
    curr_point = partition_points.pop();
    ref_ent = curr_point->owner();
    ref_edge = CAST_TO(ref_ent, RefEdge);
    if ( ref_edge == NULL )
    {
      PRINT_ERROR("Bad logic in point owners for ImprintBoundaryTool.\n");
      assert(ref_edge != NULL);
      return CUBIT_FAILURE;
    }
    CubitVector part_vec = curr_point->coordinates();
    CubitVector closest_point;
    ref_edge->closest_point(part_vec, closest_point);
    double dist = (part_vec - closest_point).length_squared();
    if ( dist > myTolerance*myTolerance )
    {
      PRINT_ERROR("Curve %d may have bad geometry.\n"
                  "While try to partition the curve the move_to operation\n"
                  "returned a move distance larger than the tolerance.\n"
                  "Please try healing or using mesh based geometry for curve %d.\n",
                  ref_edge->id(), ref_edge->id());
      return CUBIT_FAILURE;
    }
    PRINT_DEBUG_129("Partitioning curve %d\n",
                    ref_edge->id());
    RefVertex *vtx_ptr= PartitionTool::instance()->partition(ref_edge,
                                                             part_vec,
                                                             edge1, edge2);
    if ( vtx_ptr == NULL )
    {
      PRINT_ERROR("Partitioning curve failed.\n");
      return CUBIT_FAILURE;
    }
    if ( edge1->measure() < myTolerance ||
         edge2->measure() < myTolerance )
    {
      PRINT_WARNING("Virtual imprinting created an edge smaller than\n"
                    " the tolerance value.  This could be an error.\n"
                    " Check curve %d or %d, imprinting surfs %d and %d.\n",
                    edge1->id(), edge2->id(), refFacePtr1->id(),
                    refFacePtr2->id());
    }
    modified_boundary = CUBIT_TRUE;
    curr_point->owner(vtx_ptr);
      //Now readjust owners cause they can change if we partition more
      //than one curve.
    for ( ii = partition_points.size(); ii > 0; ii-- )
    {
      curr_point = partition_points.get_and_step();
      if( curr_point->owner() == ref_ent )
      {
        if ( on_curve(curr_point, edge1) )
          curr_point->owner(dynamic_cast<RefEntity*>(edge1));
        else
          curr_point->owner(dynamic_cast<RefEntity*>(edge2));
      }
    }    
  }
  return CUBIT_SUCCESS;
}
CubitStatus ImprintBoundaryTool::imprint_surface(RefFace *ref_face,
                                                 PointLoopList &part_segs,
                                                 DLIList <RefFace*> &results)
{
    //Now go through and create refedges for partitioning.
  int ii, jj, counter;
  DLIList<RefFace*> tmp_results;
  DLIList<RefEdge*> ref_edges;
  ImprintPointData *start_imp, *end_imp, *curr_imp, *matching_start, *matching_end;
  PointList *loop_ptr;
  DLIList<CubitVector*> curve_vectors;
  CubitBoolean first_on_loop;
  RefVertex *prev_end;
  RefVertex *first_v_on_loop = NULL;
  
  for ( ii = part_segs.size(); ii > 0; ii-- )
  {
    loop_ptr = part_segs.get_and_step();
    loop_ptr->reset();
    first_on_loop = CUBIT_TRUE;  
    prev_end = NULL;
    first_v_on_loop = NULL;
    for ( jj = loop_ptr->size(); jj > 0; jj-- )
    {
      start_imp = loop_ptr->get_and_step();
      if ( !start_imp->is_owner_vertex() )
      {
        PRINT_ERROR("Bad logic in imprint_surface for ImprintBoundaryTool\n");
        assert( start_imp->is_owner_vertex() );
        return CUBIT_FAILURE;
      }
      end_imp = NULL;
      curve_vectors.clean_out();
      counter = 0;
        //For the last point, only create another edge if the
        //the start is also an end...
      if ( jj == 1 && loop_ptr->get()->get_end_partition() != CUBIT_TRUE )
        continue;
        //Loop until we get to the end of the partition or
        //till we have to create a new ref_edge.
      while( loop_ptr->get()->get_end_partition() != CUBIT_TRUE )
      {
        counter++;
        curr_imp = loop_ptr->get_and_step();
        if (curr_imp->is_owner_vertex())
        {
          end_imp = curr_imp;
          counter--;
          loop_ptr->back();
          break;
        }
        else
        {
          CubitVector *temp = new CubitVector(curr_imp->coordinates());
          curve_vectors.append(temp);
        }
      }
      jj -= counter;
      CubitBoolean verts_same = CUBIT_FALSE;
      if ( end_imp == NULL )
      {
        if ( loop_ptr->get()->get_end_partition() == CUBIT_TRUE )
        {
            //Do this last segment.  This means that we have a complete loop
            //for partitioning.
          end_imp = loop_ptr->get();
          if ( curve_vectors.size() == loop_ptr->size() - 1)
            verts_same = CUBIT_TRUE;
        }
        else {
            //Either we have a screw-up or we have
            //a periodic boundary.
            //Lets say end_imp is the last in the loop.
          loop_ptr->reset();
          end_imp = loop_ptr->prev();
            //Now make sure that end_imp and start_imp are connected.
          RefEntity *ent_1 = start_imp->owner();
          RefEntity *ent_2 = end_imp->owner();
          if ( are_connected( ent_1, ent_2, ref_face ) )
          {
            verts_same = CUBIT_TRUE;
            end_imp = start_imp;
          }
          else {
            PRINT_ERROR("Problems extracting parition edge in ImprintBoundaryTool.\n");
            return CUBIT_FAILURE;
          }
        }
      }
        //Create vertices from the first and last vectors.
      CubitVector tmp_vec = start_imp->coordinates();
      matching_start = start_imp->get_matching_point();
      if ( matching_start != NULL )
      {
          //use this location instead.  This should be on 
          //the other boundary and will insure that the point
          // intersects it...
        tmp_vec = matching_start->coordinates();
      }

      RefVertex *first_ref = create_virtual_vertex(tmp_vec);
        //It is tempting to update the start_imp's owner at this point,
        //but DON'T.  These vertices belong to this surfaces boundary.  The
        //start_imp really belongs to the other boundary that is intersecting
        //with this boundary.  If you change the owner, you'll end up merging
        //vertices accross boundaries.  Don't do that...
      RefVertex *last_ref=NULL;
      if ( verts_same )
      {
        last_ref = first_ref;
        prev_end = last_ref;
      }
      else
      {
        tmp_vec = end_imp->coordinates();
        matching_end = end_imp->get_matching_point();
        if ( matching_end != NULL )
        {
            //use this location instead.  This should be on 
            //the other boundary and will insure that the point
            // intersects it...
          tmp_vec = matching_end->coordinates();
        }
        last_ref = create_virtual_vertex(tmp_vec);
          //again don't update the end_imp owner either., see above comment.
      }
      assert( first_ref && last_ref );
      RefEdge *ref_edge = NULL;
        //if curve_vectors is empty, this is just a straight line
        //curve.
      ref_edge = create_virtual_edge(first_ref, last_ref, curve_vectors);
      allocatedRefEdge->append(ref_edge);
     
      int tt;
      for ( tt = curve_vectors.size(); tt > 0; tt-- )
        delete curve_vectors.pop();
      if ( ref_edge == NULL )
      {
        PRINT_ERROR("Problem creating virtual edge in imprinting.\n");
        return CUBIT_FAILURE;
      }
      ref_edges.append(ref_edge);
      if ( first_on_loop )
      {
        first_on_loop = CUBIT_FALSE;
          //merge the start_imp's vertex with it's matching vertex.
        matching_start = start_imp->get_matching_point();
        if ( matching_start != NULL )
        {
          RefEntity *other_ent = matching_start->owner();
          RefVertex *matching_vert = CAST_TO(other_ent, RefVertex);
          if ( matching_vert == NULL )
          {
            PRINT_ERROR("Problems mergeing partition to boundary on imprint.\n");
            assert(matching_vert != NULL);
            return CUBIT_FAILURE;
          }
          CubitBoolean kept_1 = CUBIT_TRUE;
          CubitStatus temp_stat = merge_vertices(matching_vert, first_ref, kept_1);
          if ( temp_stat != CUBIT_SUCCESS )
            return CUBIT_FAILURE;
          RefVertex *owner_v = kept_1 ? matching_vert : first_ref;
            //Update the matching_start cause the matching imprint point data IS on
            //this boundary. (which is why we had to merge with it!).
          matching_start->owner(owner_v);
          if ( start_imp->get_start_partition() &&
               start_imp->get_end_partition() )
            first_v_on_loop = owner_v;
        }
        else
          first_v_on_loop = first_ref;
      }
      else if ( prev_end && last_ref != first_ref)
      {
          //merge the prev_end vert and first_ref verts.
        CubitBoolean kept_1 = CUBIT_TRUE;
        CubitStatus temp_stat = merge_vertices(prev_end, first_ref, kept_1);
        if ( temp_stat != CUBIT_SUCCESS )
          return CUBIT_FAILURE;
        prev_end = kept_1 ? prev_end : first_ref;
          //Don't update any imprint point data owners, nothing changed for this
          //boundary.
      }

      if ( last_ref != first_ref )
        prev_end = last_ref;
    }
      //Merge the end vertex to the boundary.
    loop_ptr->reset();
    end_imp = loop_ptr->prev();
    if ( end_imp->is_owner_vertex() &&
         (end_imp->get_point_type() == CREATE_NEW_VERTEX ||
          end_imp->get_point_type() == VERTEX_ON_BOTH_BOUNDARIES )
         && end_imp->get_matching_point() )
    {
      matching_end = end_imp->get_matching_point();
      RefEntity *other_ent = matching_end->owner();
      RefVertex *matching_vert = CAST_TO(other_ent, RefVertex);
        //use prev_end, hope its still valid!!!
      assert(prev_end != NULL);
      RefVertex *end_vertex = prev_end;
      if ( matching_vert == NULL || end_vertex == NULL )
      {
        PRINT_ERROR("Problems mergeing partition to boundary on imprint.\n");
        assert(matching_vert != NULL);
        assert(end_vertex != NULL);
        return CUBIT_FAILURE;
      }
      CubitBoolean kept_1 = CUBIT_TRUE;
      CubitStatus temp_stat = merge_vertices(matching_vert, end_vertex, kept_1);
      if ( temp_stat != CUBIT_SUCCESS )
        return CUBIT_FAILURE;
      RefVertex *owner_v = kept_1 ? matching_vert : end_vertex;
        //Again just update the matching_end since it is on this boundary.
      matching_end->owner(owner_v);
    }
    else if ( loop_ptr->get()->get_start_partition() &&
              loop_ptr->get()->get_end_partition() )
    {
      assert(first_v_on_loop != NULL);
        //Merge the prev_end and this vertex.
      CubitBoolean kept_1;
      CubitStatus temp_stat = merge_vertices(first_v_on_loop, prev_end, kept_1);
      if ( temp_stat != CUBIT_SUCCESS )
        return CUBIT_FAILURE;
    }
  }
  CubitStatus stat;
    //Okay, these should all be merged properly and ready to partition!
  if ( ref_edges.size() > 0 )
  {
      //Make sure that either the ref_edge_list touches the boundary at least
      //twice or is a complete loop...
    if ( valid_partition( ref_edges, ref_face ) )
    {
      PRINT_DEBUG_129("Partitioning surface %d\n",
                      ref_face->id());
      int tmp_debug = 0;
      if (tmp_debug)
      {
        GfxDebug::clear();
        GfxDebug::draw_ref_face(ref_face);
        int oo;
        for (oo=0; oo < ref_edges.size(); oo++ )
          GfxDebug::draw_ref_edge(ref_edges.next(oo), CUBIT_RED);
        GfxDebug::flush();
        GfxDebug::mouse_xforms();
      }
      stat = PartitionTool::instance()->partition(ref_face, ref_edges,
                                                  tmp_results, CUBIT_FALSE, NULL, CUBIT_TRUE);
      if ( stat != CUBIT_SUCCESS )
      {
        PRINT_ERROR("Partitioning surface %d failed\n", ref_face->id());
        return stat;
      }
      else 
      {
        results += tmp_results;
      }
    }
  }
  return CUBIT_SUCCESS;
}
CubitBoolean ImprintBoundaryTool::valid_partition(DLIList <RefEdge*> &ref_edges,
                                                  RefFace *ref_face )
{
    //Test to make sure the vertices of the ref_edges hit at least twice
    //on the boundary of the ref_face_ptr.
  DLIList <RefVertex*> boundary_verts;
  RefEdge *curr_edge;
  RefVertex *ref_vert1, *ref_vert2;
  int ii;
  //int vert_boundary_count = 0;
  
  for ( ii = ref_edges.size(); ii > 0; ii-- )
  {
    curr_edge = ref_edges.get_and_step();
    ref_vert1 = curr_edge->start_vertex();
    ref_vert2 = curr_edge->end_vertex();
    if ( ref_vert1 == ref_vert2 )
      continue;
    if ( ref_vert1->num_ref_edges() < 2 )
      return CUBIT_FALSE;
    if ( ref_vert1->num_ref_faces() > 0)
      boundary_verts.append(ref_vert1);
    if ( ref_vert2->num_ref_edges() < 2 )
      return CUBIT_FALSE;
    if ( ref_vert2->num_ref_faces() > 0 )
      boundary_verts.append(ref_vert2);
  }
  if ( ref_face->num_loops() > 1 )
  {
      //Walk on the surface.  Determine if this partition is
      //creating a sipeish partition.  We can't handle that right
      //now.
      //Start with a boundary_vert.  Determine the loop that it is
      //on.  Go on it accross the imprint.  Find its imprint.
    DLIList <Loop*> loop_list;
    DLIList <RefEdge*> ref_edge_list;
    Loop *start_loop = NULL;
    RefVertex *start_vert = NULL;
    while (boundary_verts.size())
    {
      ref_vert1 = boundary_verts.pop();
      loop_list.clean_out();
      ref_vert1->loops(loop_list);
      for ( ii = loop_list.size(); ii > 0; ii-- )
      {
        if ( loop_list.get()->get_ref_face_ptr() == ref_face )
          break;
        else
          loop_list.step();
      }
      if ( loop_list.get()->get_ref_face_ptr() != ref_face )
      {
        PRINT_ERROR("Problems in valid partition logic...");
        return CUBIT_FALSE;
      }
      start_loop = loop_list.get();
      start_vert = ref_vert1;
      do {
        loop_list.clean_out();
        ref_edge_list.clean_out();
        ref_vert1->ref_edges(ref_edge_list);
          //Find the ref_edge that has no surfaces.
        for ( ii = ref_edge_list.size(); ii > 0; ii--)
        {
          if ( ref_edge_list.get()->num_ref_faces() == 0 )
            break;
          else
            ref_edge_list.step();
        }
        if ( ref_edge_list.get()->num_ref_faces() != 0 )
        {
          PRINT_ERROR("Problems in valid partition logic...");
          return CUBIT_FALSE;
        }
        curr_edge = ref_edge_list.get();
        if ( curr_edge->start_vertex() ==
             curr_edge->end_vertex() )
          continue;
        if ( curr_edge->start_vertex() == ref_vert1 )
          ref_vert2 = curr_edge->end_vertex();
        else
          ref_vert2 = curr_edge->start_vertex();
        if ( ref_vert2->num_ref_faces() > 0 )
        {
            //This is a boundary vert.  Remove it...
          if ( !boundary_verts.move_to(ref_vert2) )
          {
              //If this has already been removed,
              //then we have a bad partition.
            return CUBIT_FALSE;
          }
          boundary_verts.remove(ref_vert2);
          loop_list.clean_out();
          ref_vert2->loops(loop_list);
          for ( ii = loop_list.size(); ii > 0; ii-- )
          {
            if ( loop_list.get()->get_ref_face_ptr() == ref_face )
              break;
            else
              loop_list.step();
          }
          if ( loop_list.get()->get_ref_face_ptr() != ref_face )
          {
            PRINT_ERROR("Problems in valid partition logic...");
            return CUBIT_FALSE;
          }
          if ( loop_list.get() == start_loop )
          {
              //This isn't a bad partition.  Break out to the outer loop.
              //The partition splits its own loop which is fine...
            break;
          }
          else if ( loop_list.get() != start_loop &&
                    boundary_verts.size() == 0 )
          {
            return CUBIT_FALSE;
          }
            //Okay, find a vertex on this loop that is a partition vertex.
          DLIList <RefVertex*> ref_vertices;
          loop_list.get()->ref_vertices(ref_vertices);
          for ( ii = ref_vertices.size(); ii > 0; ii-- )
          {
            ref_vert1 = ref_vertices.get_and_step();
            if ( ref_vert1 != ref_vert2 &&
                 boundary_verts.move_to(ref_vert1) )
            {
              break;
            }
          }
          if ( ref_vert1 == ref_vert2 ||
               !boundary_verts.move_to(ref_vert1) )
          {
              //if the vert had been removed or we have a loop
              //this is a bad partition.
            return CUBIT_FALSE;
          }
            //ref_vert1 is now set to continue on...
          boundary_verts.remove(ref_vert1);
        }
        else
          ref_vert1 = ref_vert2;
      } while (ref_vert1 != NULL );
      
    }
  }
  return CUBIT_TRUE;
} 

CubitStatus ImprintBoundaryTool::merge_vertices(DLIList <RefVertex*> &ref_verts)
{
    //This is simply a function to reduce code bloat...
  double tol_factor = myTolerance / GEOMETRY_RESABS;
  double old_geometry_factor = GeometryQueryTool::get_geometry_factor();
  GeometryQueryTool::set_geometry_factor( tol_factor );
    //merge the vertices to get a connected intersection graph.
  CubitStatus stat = MergeTool::instance()->merge_refvertices( ref_verts );
  GeometryQueryTool::set_geometry_factor( old_geometry_factor );
  if ( stat != CUBIT_SUCCESS )
    return stat;
  return CUBIT_SUCCESS;
}
CubitStatus ImprintBoundaryTool::merge_vertices(RefVertex *ref_vertex1,
                                                RefVertex *ref_vertex2,
                                                CubitBoolean &kept_1)
{

    //force merge the vertices
  RefVertex* result =  MergeTool::instance()->
                          force_merge( ref_vertex1, ref_vertex2 );
  if (result == ref_vertex1)
    kept_1 = CUBIT_TRUE;
  else if(result == ref_vertex2)
    kept_1 = CUBIT_FALSE;
  else
    return CUBIT_FAILURE;

  return CUBIT_SUCCESS;
}
void ImprintBoundaryTool::draw_seg(ImprintLineSegment *seg, int color)
{
  CubitVector start_point = seg->get_start()->coordinates();
  CubitVector end_point = seg->get_end()->coordinates();
  if ( color == -1 )
    color = CUBIT_GREEN;
  GfxDebug::draw_line( start_point.x(),start_point.y(),start_point.z(),
                          end_point.x(),end_point.y(),end_point.z(),
                          color);
  GfxDebug::flush();
}
void ImprintBoundaryTool::draw_loops(PointLoopList &boundary_loops)
{
  int ii, jj;
  PointList *point_list;
  ImprintPointData *imp_point;
  for (ii = 0; ii < boundary_loops.size(); ii++ )
  {
    point_list = boundary_loops.get_and_step();
    for ( jj = 0; jj < point_list->size(); jj++ )
    {
      imp_point = point_list->get_and_step();
      draw_point(imp_point);
    }
  }
}

void ImprintBoundaryTool::draw_point(ImprintPointData *imp_point, int color)
{
  CubitPoint *point = (CubitPoint*) imp_point;
  PointType p_type = imp_point->get_point_type();
  if (color == -1 )
  {
    switch ( p_type )
    {
      case VERTEX_ON_BOUNDARY_1:
        color = CUBIT_RED;
        break;
      case VERTEX_ON_BOUNDARY_2:
        color = CUBIT_ORANGE;
        break;
      case VERTEX_ON_BOTH_BOUNDARIES:
        color = CUBIT_MAGENTA;
        break;
      case CREATE_NEW_VERTEX:
        color = CUBIT_CYAN;
        break;
      case ON_BOUNDARY_1:
        color = CUBIT_GREEN;
        break;
      case ON_BOUNDARY_2:
        color = CUBIT_YELLOW;
        break;
      case ON_BOTH_BOUNDARIES:
        color = CUBIT_BLUE;
        break;
      default:
        color = CUBIT_BROWN;
    }
  }
//  point->draw(color);
  int id_val = point->id();
  GfxDebug::draw_label(id_val, point->x(), point->y(), point->z(), color);
  GfxDebug::flush();
}
void ImprintBoundaryTool::draw_end(ImprintLineSegment *seg)
{
  ImprintPointData *end = seg->get_end();
  CubitPoint *point = (CubitPoint*) end;
  PointType p_type = end->get_point_type();
  int color;
  switch ( p_type )
  {
    case VERTEX_ON_BOUNDARY_1:
      color = CUBIT_RED;
      break;
    case VERTEX_ON_BOUNDARY_2:
      color = CUBIT_ORANGE;
      break;
    case VERTEX_ON_BOTH_BOUNDARIES:
      color = CUBIT_MAGENTA;
      break;
    case CREATE_NEW_VERTEX:
      color = CUBIT_CYAN;
      break;
    case ON_BOUNDARY_1:
      color = CUBIT_GREEN;
      break;
    case ON_BOUNDARY_2:
      color = CUBIT_YELLOW;
      break;
    case ON_BOTH_BOUNDARIES:
      color = CUBIT_BLUE;
      break;
    default:
      color = CUBIT_BROWN;
  }
  int id_val = point->id();
  GfxDebug::draw_label(id_val, point->x(), point->y(), point->z(), color);
  GfxDebug::flush();
}
CubitStatus ImprintBoundaryTool::resolve_on_boundaries( PointLoopList &boundary_loops_1,
                                                        PointLoopList &boundary_loops_2)
{
  int ii, jj;
  PointList *loop_ptr;
  ImprintPointData *curr_point, *next_point, *prev_point, *matching_point;
  
    //loop around the boundaries, if the point goes from on both boundaries to just on boundary,
    //then we know we need to create a vertex there.
  for ( ii = boundary_loops_1.size(); ii > 0; ii-- )
  {
    loop_ptr = boundary_loops_1.get_and_step();
    for ( jj = loop_ptr->size(); jj > 0; jj-- )
    {
      curr_point = loop_ptr->get_and_step();
      if ( curr_point->get_point_type() == ON_BOTH_BOUNDARIES )
      {
        next_point = loop_ptr->get();
        prev_point = loop_ptr->prev(2);

        if ( (next_point->get_point_type() == ON_BOUNDARY_1 ||
              next_point->get_point_type() == VERTEX_ON_BOUNDARY_1) ||
             (prev_point->get_point_type() == ON_BOUNDARY_1 ||
              prev_point->get_point_type() == VERTEX_ON_BOUNDARY_1 ) )
          
        {
          curr_point->set_point_type(CREATE_NEW_VERTEX);
            //If curr_point is on both boundaries, then the
            //matching point cant be a vertex either so
            //create a new vertex there too.
          matching_point =curr_point->get_matching_point();
          matching_point->set_point_type(CREATE_NEW_VERTEX);
        }
      }
    }
  }
  for ( ii = boundary_loops_2.size(); ii > 0; ii-- )
  {
    loop_ptr = boundary_loops_2.get_and_step();
    for ( jj = loop_ptr->size(); jj > 0; jj-- )
    {
      curr_point = loop_ptr->get_and_step();
      if ( curr_point->get_point_type() == ON_BOTH_BOUNDARIES )
      {
        next_point = loop_ptr->get();
        prev_point = loop_ptr->prev(2);
        if ( (next_point->get_point_type() == ON_BOUNDARY_2 ||
              next_point->get_point_type() == VERTEX_ON_BOUNDARY_2) ||
             (prev_point->get_point_type() == ON_BOUNDARY_2 ||
              prev_point->get_point_type() == VERTEX_ON_BOUNDARY_2 ) )
        {
          curr_point->set_point_type(CREATE_NEW_VERTEX);
            //If curr_point is on both boundaries, then the
            //matching point cant be a vertex either so
            //create a new vertex there too.
          matching_point =curr_point->get_matching_point();
          matching_point->set_point_type(CREATE_NEW_VERTEX);
        }
      }
    }
  }
    
  return CUBIT_SUCCESS;
}
CubitBoolean ImprintBoundaryTool::resolve_match_conflict_other(ImprintPointData *this_point,
                                                               ImprintPointData *other,
                                                               CubitBoolean this_is_1 )
{
  CubitVector other_v = other->coordinates();
  CubitVector coords = this_point->coordinates();
  ImprintPointData *conflict_match = other->get_matching_point();
  CubitVector c_match_vec = conflict_match->coordinates();
  double dist_1 = (other_v - c_match_vec).length_squared();
  double dist_2 = (other_v - coords).length_squared();
  if ( dist_1 < dist_2 ) //this says leave it alone, and don't match these two.
  {
      //The conflict match and other should be left matching.  Don't
      //match this_point and other...
    return CUBIT_FALSE;
  }
  else
  {
      //Okay we need to clean this up.
      //It used to be that conflict_match and other were paired
      //up.  But now we need to undo that and match other and this_point.
      //To do that find out what the correct PointType setting should be
      //for conflict_match.  Try setting it back to its original.
    RefEntity* owner_conflict = conflict_match->owner();
      //Use the this_is_1 flag.  Assume conflict_match and
      //this_point are on the same boundary.
    if( CAST_TO(owner_conflict, RefEdge) )
    {
      if ( this_is_1 )
        conflict_match->set_point_type(ON_BOUNDARY_1);
      else
        conflict_match->set_point_type(ON_BOUNDARY_2);
    }
    else if ( CAST_TO(owner_conflict, RefVertex) )
    {
      if ( this_is_1 )
        conflict_match->set_point_type(VERTEX_ON_BOUNDARY_1);
      else
        conflict_match->set_point_type(VERTEX_ON_BOUNDARY_2);
    }
    else
    {
      PRINT_ERROR("Problems resolving matching point type.\n");
      return CUBIT_FALSE;
    }
      //Now do the same thing for other.  Clean it back to its
      //preset types. (except use the opposite of this_is_1 flag
      //since it is on the other boundary of this_point.)
    RefEntity* owner_other = other->owner();
    if( CAST_TO(owner_other,RefEdge) )
    {
      if ( this_is_1 )
        other->set_point_type(ON_BOUNDARY_2);
      else
        other->set_point_type(ON_BOUNDARY_1);
    }
    else if ( CAST_TO(owner_other,RefVertex))
    {
      if ( this_is_1 )
        other->set_point_type(VERTEX_ON_BOUNDARY_2);
      else
        other->set_point_type(VERTEX_ON_BOUNDARY_1);
    }
    else
    {
      PRINT_ERROR("Problems resolving matching point type.\n");
      return CUBIT_FALSE;
    }
    return CUBIT_TRUE;
  }
}
  
CubitBoolean ImprintBoundaryTool::resolve_match_conflict_this(ImprintPointData *this_point,
                                                              ImprintPointData *other,
                                                              CubitBoolean this_is_1 )
{
  ImprintPointData *conflict_match = this_point->get_matching_point();
  CubitVector other_v = other->coordinates();
  CubitVector coords = this_point->coordinates();
  CubitVector c_match_vec = conflict_match->coordinates();
  double dist_1 = (coords - c_match_vec).length_squared();
  double dist_2 = (coords - other_v).length_squared();
  if ( dist_1 < dist_2 ) //this says leave it alone, and don't match these two.
  {
      //The conflict match and this_point should be left matching.  Don't
      //match this_point and other...
    return CUBIT_FALSE;
  }
  else
  {
      //Okay we need to clean this up.
      //It used to be that conflict_match and this_point were paired
      //up.  But now we need to undo that and match other and this_point.
      //To do that find out what the correct PointType setting should be
      //for conflict_match.  Try setting it back to its original.
    RefEntity* owner_conflict = conflict_match->owner();
      //Use the this_is_1 flag.  Assume conflict_match and
      //other are on the same boundary.
    if( CAST_TO(owner_conflict,RefEdge) )
    {
      if ( this_is_1 )
        conflict_match->set_point_type(ON_BOUNDARY_2);
      else
        conflict_match->set_point_type(ON_BOUNDARY_1);
    }
    else if ( CAST_TO(owner_conflict,RefVertex) )
    {
      if ( this_is_1 )
        conflict_match->set_point_type(VERTEX_ON_BOUNDARY_2);
      else
        conflict_match->set_point_type(VERTEX_ON_BOUNDARY_1);
    }
    else
    {
      PRINT_ERROR("Problems resolving matching point type.\n");
      return CUBIT_FALSE;
    }
      //Now do the same thing for this_point.  Clean it back to its
      //preset types. (except use the opposite of this_is_1 flag
      //since it is on the other boundary of this_point.)
    RefEntity* owner_this = this_point->owner();
    if( CAST_TO(owner_this,RefEdge))
    {
      if ( this_is_1 )
        this_point->set_point_type(ON_BOUNDARY_1);
      else
        this_point->set_point_type(ON_BOUNDARY_2);
    }
    else if (CAST_TO(owner_this,RefVertex) )
    {
      if ( this_is_1 )
        this_point->set_point_type(VERTEX_ON_BOUNDARY_1);
      else
        this_point->set_point_type(VERTEX_ON_BOUNDARY_2);
    }
    else
    {
      PRINT_ERROR("Problems resolving matching point type.\n");
      return CUBIT_FALSE;
    }
    return CUBIT_TRUE;
  }
}
CubitStatus ImprintBoundaryTool::match_points( ImprintLineSegment *seg_1,
                                               ImprintLineSegment *seg_2,
                                               MatchType &type_0,
                                               MatchType &type_1,
                                               MatchType &type_2,
                                               MatchType &type_3)
{
  ImprintPointData *imp_point_0, *imp_point_1, *imp_point_2, *imp_point_3;
  ImprintPointData *imp_point_4, *imp_point_5, *imp_point_6, *imp_point_7;
  ImprintLineSegment *prev_seg_1 = seg_1->get_prev();
  ImprintLineSegment *next_seg_1 = seg_1->get_next();
  ImprintLineSegment *prev_seg_2 = seg_2->get_prev();
  ImprintLineSegment *next_seg_2 = seg_2->get_next();
  imp_point_0 = seg_1->get_start();
  imp_point_1 = seg_1->get_end();
  imp_point_4 = prev_seg_1->get_start();
  imp_point_5 = next_seg_1->get_end();
  imp_point_2 = seg_2->get_start();
  imp_point_3 = seg_2->get_end();
  imp_point_6 = prev_seg_2->get_start();
  imp_point_7 = next_seg_2->get_end();
  
  CubitVector point_0, point_1, point_2, point_3;
  CubitVector point_4, point_5, point_6, point_7;
  point_0 = imp_point_0->coordinates();
  point_1 = imp_point_1->coordinates();
  point_2 = imp_point_2->coordinates();
  point_3 = imp_point_3->coordinates();
  point_4 = imp_point_4->coordinates();
  point_5 = imp_point_5->coordinates();
  point_6 = imp_point_6->coordinates();
  point_7 = imp_point_7->coordinates();

    //There are four possible matches for each of the points (0,1,2,and3).
    //Find out which ones match.
  CubitBoolean match_02=CUBIT_FALSE, match_03=CUBIT_FALSE,
    match_06=CUBIT_FALSE, match_07=CUBIT_FALSE;
  CubitBoolean match_12=CUBIT_FALSE, match_13=CUBIT_FALSE,
    match_16=CUBIT_FALSE, match_17=CUBIT_FALSE;
  CubitBoolean match_24=CUBIT_FALSE, match_25=CUBIT_FALSE;
  CubitBoolean match_34=CUBIT_FALSE, match_35=CUBIT_FALSE;
    //compute the matching.
  if (point_0.within_tolerance(point_2, myTolerance) )
    match_02=CUBIT_TRUE;
  if (point_0.within_tolerance(point_3, myTolerance) )
    match_03=CUBIT_TRUE;
  if (point_0.within_tolerance(point_6, myTolerance) )
    match_06=CUBIT_TRUE;
  if (point_0.within_tolerance(point_7, myTolerance) )
    match_07=CUBIT_TRUE;

  if (point_1.within_tolerance(point_2, myTolerance) )
    match_12=CUBIT_TRUE;
  if (point_1.within_tolerance(point_3, myTolerance) )
    match_13=CUBIT_TRUE;
  if (point_1.within_tolerance(point_6, myTolerance) )
    match_16=CUBIT_TRUE;
  if (point_1.within_tolerance(point_7, myTolerance) )
    match_17=CUBIT_TRUE;

  if (point_2.within_tolerance(point_4, myTolerance) )
    match_24=CUBIT_TRUE;
  if (point_2.within_tolerance(point_5, myTolerance) )
    match_25=CUBIT_TRUE;
  if (point_3.within_tolerance(point_4, myTolerance) )
    match_34=CUBIT_TRUE;
  if (point_3.within_tolerance(point_5, myTolerance) )
    match_35=CUBIT_TRUE;
    //Now for those that match, compute the distances.
  double dist_02 = 0.0, dist_03 = 0.0, dist_06 = 0.0, dist_07 = 0.0;
  double dist_12 = 0.0, dist_13 = 0.0, dist_16 = 0.0, dist_17 = 0.0;
  double dist_24 = 0.0, dist_25 = 0.0;
  double dist_34 = 0.0, dist_35 = 0.0;
  if ( match_02 )
    dist_02 = (point_0-point_2).length_squared();
  if ( match_03 )
    dist_03 = (point_0-point_3).length_squared();
  if ( match_06 )
    dist_06 = (point_0-point_6).length_squared();
  if ( match_07 )
    dist_07 = (point_0-point_7).length_squared();

  if ( match_12 )
    dist_12 = (point_1-point_2).length_squared();
  if ( match_13 )
    dist_13 = (point_1-point_3).length_squared();
  if ( match_16 )
    dist_16 = (point_1-point_6).length_squared();
  if ( match_17 )
    dist_17 = (point_1-point_7).length_squared();

  if ( match_24 )
    dist_24 = (point_2-point_4).length_squared();
  if ( match_25 )
    dist_25 = (point_2-point_5).length_squared();

  if ( match_34 )
    dist_34 = (point_3-point_4).length_squared();
  if ( match_35 )
    dist_35 = (point_3-point_5).length_squared();

    //Okay, determine the solution by finding the closest
    //points, and then determining the final four results;

    //point 0
  double min_dist = CUBIT_DBL_MAX;
  type_0 = NO_MATCH;
  type_1 = NO_MATCH;
  type_2 = NO_MATCH;
  type_3 = NO_MATCH;
  if ( match_02 && dist_02 < min_dist )
  {
    min_dist = dist_02;
    type_0 = MATCH_0_2;
  }
  if ( match_03 && dist_03 < min_dist )
  {
    min_dist = dist_03;
    type_0 = MATCH_0_3;
  }
  if ( match_06 && dist_06 < min_dist )
  {
    min_dist = dist_06;
    type_0 = MATCH_0_6;
  }
  if ( match_07 && dist_07 < min_dist )
  {
    min_dist = dist_07;
    type_0 = MATCH_0_7;
  }
    //point 1
  min_dist = CUBIT_DBL_MAX;
  type_1 = NO_MATCH;
  if ( match_12 && dist_12 < min_dist )
  {
    min_dist = dist_12;
    type_1 = MATCH_1_2;
  }
  if ( match_13 && dist_13 < min_dist )
  {
    min_dist = dist_13;
    type_1 = MATCH_1_3;
  }
  if ( match_16 && dist_16 < min_dist )
  {
    min_dist = dist_16;
    type_1 = MATCH_1_6;
  }
  if ( match_17 && dist_17 < min_dist )
  {
    min_dist = dist_17;
    type_1 = MATCH_1_7;
  }
    //point 2
  min_dist = CUBIT_DBL_MAX;
  type_2 = NO_MATCH;
  if ( match_02 && dist_02 < min_dist )
  {
    min_dist = dist_02;
    type_2 = MATCH_0_2;
  }
  if ( match_12 && dist_12 < min_dist )
  {
    min_dist = dist_12;
    type_2 = MATCH_1_2;
  }
  if ( match_24 && dist_24 < min_dist )
  {
    min_dist = dist_24;
    type_2 = MATCH_2_4;
  }
  if ( match_25 && dist_25 < min_dist )
  {
    min_dist = dist_25;
    type_2 = MATCH_2_5;
  }
    //point 3
  min_dist = CUBIT_DBL_MAX;
  type_3 = NO_MATCH;
  if ( match_03 && dist_03 < min_dist )
  {
    min_dist = dist_03;
    type_3 = MATCH_0_3;
  }
  if ( match_13 && dist_13 < min_dist )
  {
    min_dist = dist_13;
    type_3 = MATCH_1_3;
  }
  if ( match_34 && dist_34 < min_dist )
  {
    min_dist = dist_34;
    type_3 = MATCH_3_4;
  }
  if ( match_35 && dist_35 < min_dist )
  {
    min_dist = dist_35;
    type_3 = MATCH_3_5;
  }
  return CUBIT_SUCCESS;
}

  
  
  
int ImprintBoundaryTool::num_coedges_on_face(RefEdge *edge_ptr,
                                             RefFace *ref_face)
{
  DLIList <CoEdge*> coedges;
  edge_ptr->co_edges(coedges);
  int ii;
  int count = 0;
  for ( ii = coedges.size(); ii > 0; ii-- )
  {
    if ( coedges.get_and_step()->get_ref_face() == ref_face )
      count++;
  }
  return count;
}

CubitStatus ImprintBoundaryTool::imprint_segments( SegLoopList &boundary_line_loops_1,
                                                   SegLoopList &boundary_line_loops_2,
                                                   PointLoopList &boundary_loops_1,
                                                   PointLoopList &boundary_loops_2)
{
  int ii, jj;
  ImprintLineSegment **surf_1_loop_heads, **surf_2_loop_heads;
  surf_1_loop_heads = new ImprintLineSegment*[boundary_line_loops_1.size()];
  surf_2_loop_heads = new ImprintLineSegment*[boundary_line_loops_2.size()];
  SegList *line_loop_1, *line_loop_2;
  ImprintLineSegment *curr_seg;
  
  for ( ii = 0; ii < boundary_line_loops_1.size(); ii++ )
  {
    line_loop_1 = boundary_line_loops_1.get_and_step();
    surf_1_loop_heads[ii] = line_loop_1->get();
  }
  for ( ii = 0; ii < boundary_line_loops_2.size(); ii++ )
  {
    line_loop_2 = boundary_line_loops_2.get_and_step();
    surf_2_loop_heads[ii] = line_loop_2->get();
  }
  KDDTree <ImprintLineSegment*> atree_1(myTolerance, true);
  KDDTree <ImprintLineSegment*> atree_2(myTolerance, true);
    //Insert all of the segments on for the boundary line loops into
    //their respective R-Tree.
  for ( ii = 0; ii < boundary_line_loops_1.size(); ii++ )
  {
    line_loop_1 = boundary_line_loops_1.get_and_step();
    for ( jj = 0; jj < line_loop_1->size(); jj++ )
    {
      curr_seg = line_loop_1->get_and_step();
      atree_1.add(curr_seg);
    }
  }
  atree_1.balance();
  for ( ii = 0; ii < boundary_line_loops_2.size(); ii++ )
  {
    line_loop_2 = boundary_line_loops_2.get_and_step();
    for ( jj = 0; jj < line_loop_2->size(); jj++ )
    {
      curr_seg = line_loop_2->get_and_step();
      atree_2.add(curr_seg);
    }
  }
  atree_2.balance();
    //Now for each segment on boundary 1, find the closest segments on
    //boundary 2.  Determine for each point on boundary 1 the closest
    //point on boundary 2 if it is within tolerance.  Set this
    //data in the pointmatched data in ImprintPointData.
  CubitStatus stat = find_matching_points(boundary_line_loops_1,
                                          atree_2 );
  if ( stat != CUBIT_SUCCESS )
    return stat;

    //Now for each segment on boundary 2, find the closest segments on
    //boundary 1.  Determine for each point on boundary 2 the closest
    //point on boundary 1 if it is within tolerance.  Set this
    //data in the pointmatched data in ImprintPointData.
  stat = find_matching_points(boundary_line_loops_2,
                              atree_1 );
  if ( stat != CUBIT_SUCCESS )
    return stat;

    //Now we need to go over each node on the boundaries and resolve
    //descripensies.  Find nodes that match each other and set them
    //to match each other. Also at this time.  Nodes that have
    //closest points on the interior of a segment will get the
    //segment split.  Also as segments are modified, update the
    //atree.
  int num_loops_1 = boundary_line_loops_1.size();
  int num_loops_2 = boundary_line_loops_2.size();
  stat = final_match( surf_1_loop_heads, num_loops_1,
                      surf_2_loop_heads, num_loops_2,
                      atree_1, atree_2);
  if ( stat != CUBIT_SUCCESS )
    return stat;

    //Unfortunately, we are not done.  We now need to find
    //cross intersections.  These are where two segments of different
    //surfaces intersect.  They don't intersect at the points themselves,
    //so they weren't captured previously.  If I could do this in any other
    //part of the previous code, I would.  I just can't figure out
    //how to do it without sacrificing robustness.
    //Note that only atree_2 is passed in.  After this point,
    //atree_1 will be out of date.
  stat = find_crossings(surf_1_loop_heads, num_loops_1,
                        surf_2_loop_heads, num_loops_2,
                        atree_2);
  if ( stat != CUBIT_SUCCESS )
    return stat;

    //Recaputure the linked list into the DLIList structure.
  update_boundary_loops( boundary_loops_1, surf_1_loop_heads);
  update_boundary_loops( boundary_loops_2, surf_2_loop_heads);

  delete [] surf_1_loop_heads;
  delete [] surf_2_loop_heads;

  stat = resolve_on_boundaries( boundary_loops_1,
                                boundary_loops_2);
  if ( stat != CUBIT_SUCCESS )
    return stat;
  if ( DEBUG_FLAG(129) )
  {
    GfxDebug::clear();
    draw_loops(boundary_loops_1);
    draw_loops(boundary_loops_2);
    GfxDebug::mouse_xforms();
  }
  
  return CUBIT_SUCCESS;
}
CubitStatus ImprintBoundaryTool::find_matching_points( SegLoopList &seg_loops,
                                                       AbstractTree <ImprintLineSegment*> &atree )
{
    //Loop over seg_loops.  For each segment loop, find the closest matching point
    //from the segments in the atree.  
  int ii, jj, kk;
  SegList *loop_ptr;
  ImprintLineSegment *curr_seg, *test_seg;
  ImprintPointData *imp_point_0, *imp_point_1;
  SegList close_segments;
  CubitBox curr_box;
  CubitStatus stat;
  for ( ii = 0; ii < seg_loops.size(); ii++ )
  {
    loop_ptr = seg_loops.get_and_step();
    for ( jj = 0; jj < loop_ptr->size(); jj++ )
    {
      curr_seg = loop_ptr->get_and_step();
      close_segments.clean_out();
      curr_box = curr_seg->bounding_box();
      atree.find(curr_box, close_segments);
	  if ( close_segments.size() == 0 )
		  continue;

      for ( kk = 0; kk < close_segments.size(); kk++ )
      {
        test_seg = close_segments.get_and_step();
        stat = find_closest_points(curr_seg, test_seg);
        if (stat != CUBIT_SUCCESS )
          return stat;
      }
        //Now given those closest points,
        //find the match for point 0.
      imp_point_0 = curr_seg->get_start();
      if (!imp_point_0->is_matched() )
        set_closest_point(imp_point_0);
        //Now given those closest points,
        //find the match for point 1.
      imp_point_1 = curr_seg->get_end();
      if (!imp_point_1->is_matched() )
        set_closest_point(imp_point_1);
    }
  }
  return CUBIT_SUCCESS;
}
void ImprintBoundaryTool::set_closest_point(ImprintPointData *imp_point_0)
{
  DLIList <ImprintMatchData*> closest_points;
  int kk;
  imp_point_0->get_match_list(closest_points);
  double min_dist = CUBIT_DBL_MAX;
  ImprintMatchData *closest_point = NULL, *close_point;
  for ( kk = 0; kk < closest_points.size(); kk++ )
  {
    close_point = closest_points.get_and_step();
    if ( close_point->get_closest_dist() < min_dist )
    {
      min_dist = close_point->get_closest_dist();
      closest_point = close_point;
    }
  }
  if ( closest_point )
  {
    imp_point_0->set_matched(closest_point);
  }
}

CubitStatus ImprintBoundaryTool::find_closest_points(ImprintLineSegment *curr_seg,
                                                     ImprintLineSegment *test_seg )
{
    //Given the two
  ImprintPointData *curr_data_a, *curr_data_b;
  curr_data_a = curr_seg->get_start();
  curr_data_b = curr_seg->get_end();
  ImprintMatchData *point_a_match = NULL;
  ImprintMatchData *point_b_match = NULL;
  
  CubitStatus stat1 = match_seg_points(curr_seg, test_seg, point_a_match,
                                       point_b_match);
  if ( stat1 != CUBIT_SUCCESS )
    return stat1;
  if ( !curr_data_a->is_matched() && 
	   point_a_match != NULL )
    curr_data_a->add_match_data(point_a_match);
  if ( !curr_data_b->is_matched() && 
	   point_b_match != NULL )
    curr_data_b->add_match_data(point_b_match);
  return CUBIT_SUCCESS;
}
CubitStatus ImprintBoundaryTool::match_seg_points( ImprintLineSegment *seg_1,
                                                   ImprintLineSegment *seg_2,
                                                   ImprintMatchData *&point_0_match,
                                                   ImprintMatchData *&point_1_match )
{
  ImprintPointData *imp_point_0, *imp_point_1, *imp_point_2, *imp_point_3;
  imp_point_0 = seg_1->get_start();
  imp_point_1 = seg_1->get_end();
  imp_point_2 = seg_2->get_start();
  imp_point_3 = seg_2->get_end();

  
  CubitVector point_0, point_1, point_2, point_3;
  point_0 = imp_point_0->coordinates();
  point_1 = imp_point_1->coordinates();
  point_2 = imp_point_2->coordinates();
  point_3 = imp_point_3->coordinates();
  CubitBoolean dont_test_0 = imp_point_0->is_matched();
  CubitBoolean dont_test_1 = imp_point_1->is_matched();
  point_0_match = NULL;
  point_1_match = NULL;
  if ( dont_test_0 && dont_test_1 )
    return CUBIT_SUCCESS;
  
  int debug = 0;
  if ( debug )
  {
    draw_seg(seg_1, CUBIT_GREEN);
    draw_seg(seg_2, CUBIT_YELLOW);
    draw_point(imp_point_0, CUBIT_GREEN);
    draw_point(imp_point_1, CUBIT_GREEN);
    draw_point(imp_point_2, CUBIT_YELLOW);
    draw_point(imp_point_3, CUBIT_YELLOW);
    GfxDebug::mouse_xforms();
  }
    //There are four possible matches for each of the points (0,1,2,and3).
    //Find out which ones match.
  CubitBoolean match_02=CUBIT_FALSE, match_03=CUBIT_FALSE;
  CubitBoolean match_12=CUBIT_FALSE, match_13=CUBIT_FALSE;
    //compute the matching.
  if ( !dont_test_0 )
  {
    if (point_0.within_tolerance(point_2, myTolerance) )
      match_02=CUBIT_TRUE;
    if (point_0.within_tolerance(point_3, myTolerance) )
      match_03=CUBIT_TRUE;
  }
  if ( !dont_test_1 )
  {
    if (point_1.within_tolerance(point_2, myTolerance) )
      match_12=CUBIT_TRUE;
    if (point_1.within_tolerance(point_3, myTolerance) )
      match_13=CUBIT_TRUE;
  }
    //Now for those that match, compute the distances.
  double dist_02 = 0.0, dist_03 = 0.0;
  double dist_12 = 0.0, dist_13 = 0.0;
  if ( match_02 )
    dist_02 = (point_0-point_2).length_squared();
  if ( match_03 )
    dist_03 = (point_0-point_3).length_squared();
  if ( match_12 )
    dist_12 = (point_1-point_2).length_squared();
  if ( match_13 )
    dist_13 = (point_1-point_3).length_squared();

  if ( match_02 && match_03 )
  {
    if ( dist_02 < dist_03 )
      match_03 = CUBIT_FALSE;
    else
      match_02 = CUBIT_FALSE;
  }
  if ( match_12 && match_13 )
  {
    if ( dist_12 < dist_13 )
      match_13 = CUBIT_FALSE;
    else
      match_12 = CUBIT_FALSE;
  }

    //Now find the closest points on the segment.
  CubitBoolean on_interior_0 = CUBIT_FALSE;
  CubitBoolean on_interior_1 = CUBIT_FALSE;
  CubitVector closest_point_0, closest_point_1;
  if ( !dont_test_0 && closest_point_interior_seg(imp_point_0,
                                                  seg_2,
                                                  closest_point_0 ) )
  {
    if ( point_0.within_tolerance(closest_point_0,
                                  myTolerance) )
    {
        //There is a close point on the interior.
      on_interior_0 = CUBIT_TRUE;
    }
  }
  if ( !dont_test_1 && closest_point_interior_seg(imp_point_1,
                                                  seg_2,
                                                  closest_point_1 ) )
  {
      //This means the closest_point is perpendicular
      //to the segment.  Now test if the point is within
      //tolerance.
    if ( point_1.within_tolerance(closest_point_1,
                                  myTolerance) )
    {
        //There is a close point on the interior.
      on_interior_1 = CUBIT_TRUE;
    }
  }
    //point 0
  //double min_dist = CUBIT_DBL_MAX;
  point_0_match = NULL;
  point_1_match = NULL;
  if ( match_02 )
  {
    point_0_match = new ImprintMatchData; 
    if ( on_interior_0 )
    {
        //set the distance to be equal to the
        //distance to the closest_point_0 and point_0 rather
        //than 2 and 0.  0 and 2 should snap together but
        //for finding the min distance, we want to make sure
        //we get the closest segment.
      dist_02 = (point_0 - closest_point_0).length_squared();
    }

    allocatedMatchData->append(point_0_match);
    point_0_match->set_closest_point(imp_point_2);
    point_0_match->set_closest_dist(dist_02);
    point_0_match->set_closest_seg(seg_2);
    seg_2->add_match_data(point_0_match);
  }
  else if ( match_03 )
  {
    point_0_match = new ImprintMatchData;
    if ( on_interior_0 )
    {
        //set the distance to be equal to the
        //distance to the closest_point_0 and point_0 rather
        //than 3 and 0.  0 and 3 should snap together but
        //for finding the min distance, we want to make sure
        //we get the closest segment.
      dist_03 = (point_0 - closest_point_0).length_squared();
    }

    allocatedMatchData->append(point_0_match);
    point_0_match->set_closest_point(imp_point_3);
    point_0_match->set_closest_dist(dist_03);
    point_0_match->set_closest_seg(seg_2);
    seg_2->add_match_data(point_0_match);
  }
  else if ( on_interior_0 && !match_02 && !match_03 )
  {
    double dist = (point_0 - closest_point_0).length_squared();
    point_0_match = new ImprintMatchData;
    allocatedMatchData->append(point_0_match);
    point_0_match->set_closest_seg(seg_2);
    seg_2->add_match_data(point_0_match);
    point_0_match->set_point_on(closest_point_0);
    point_0_match->set_closest_dist(dist);
  }
  if ( match_12 )
  {
    point_1_match = new ImprintMatchData;  
	if ( on_interior_1 )
    {
        //set the distance to be equal to the
        //distance to the closest_point_1 and point_1 rather
        //than 2 and 1.  1 and 2 should snap together but
        //for finding the min distance, we want to make sure
        //we get the closest segment.
      dist_12 = (point_1 - closest_point_1).length_squared();
    }

    allocatedMatchData->append(point_1_match);
    point_1_match->set_closest_point(imp_point_2);
    point_1_match->set_closest_dist(dist_12);
    point_1_match->set_closest_seg(seg_2);
    seg_2->add_match_data(point_1_match);
  }
  else if ( match_13 )
  {    
    point_1_match = new ImprintMatchData;
    if ( on_interior_1 )
    {
        //set the distance to be equal to the
        //distance to the closest_point_1 and point_1 rather
        //than 3 and 1.  1 and 3 should snap together but
        //for finding the min distance, we want to make sure
        //we get the closest segment.
      dist_13 = (point_1 - closest_point_1).length_squared();
    }

    allocatedMatchData->append(point_1_match);
    point_1_match->set_closest_point(imp_point_3);
    point_1_match->set_closest_dist(dist_13);
    point_1_match->set_closest_seg(seg_2);    
    seg_2->add_match_data(point_1_match);
  }
  else if ( on_interior_1 && !match_12 && !match_13 )
  {
    double dist = (point_1 - closest_point_1).length_squared();
    //double debug_dist = (point_0 - point_2).length();
    point_1_match = new ImprintMatchData;
    allocatedMatchData->append(point_1_match);
    point_1_match->set_closest_seg(seg_2);
    point_1_match->set_point_on(closest_point_1);
    point_1_match->set_closest_dist(dist);
    seg_2->add_match_data(point_1_match);
  }
  return CUBIT_SUCCESS;
}
CubitStatus ImprintBoundaryTool::final_match( ImprintLineSegment **surf_1_loop_heads,
                                              int num_loops_1,
                                              ImprintLineSegment **surf_2_loop_heads,
                                              int num_loops_2,
                                              AbstractTree <ImprintLineSegment*>& atree_1,
                                              AbstractTree <ImprintLineSegment*>& atree_2)
{
  CubitStatus stat;
    //Match the points on surface 1 with the points on surface 2.
    //Resolve cases where more than 1 point match to another point.
    //Also resolve where a point on surface 1 matches with a segment
    //on surface 2.  This will possibly alter the loop_heads of surface 2,
    //and the atree of surface2.
  stat = final_match_loop( surf_1_loop_heads, num_loops_1,
                           surf_2_loop_heads, num_loops_2,
                           atree_1, atree_2);
  if (stat != CUBIT_SUCCESS )
    return stat;
    //Now match the points on surface 2 with surface 1.
    //This should go faster since hopefully they were done
    //at once.  But there will be points that match with segments
    //that are left...
    //This will possibly alter the loop_heads of surface 1 and
    //the atree of surface1.
  stat = final_match_loop( surf_2_loop_heads, num_loops_2,
                           surf_1_loop_heads, num_loops_1,
                           atree_2, atree_1);
  if (stat != CUBIT_SUCCESS )
    return stat;
  
  return CUBIT_SUCCESS;
}
CubitStatus ImprintBoundaryTool::final_match_loop( ImprintLineSegment **surf_1_loop_heads,
                                                   int num_loops_1,
                                                   ImprintLineSegment **surf_2_loop_heads,
                                                   int num_loops_2,
                                                   AbstractTree <ImprintLineSegment*>& atree_1,
                                                   AbstractTree <ImprintLineSegment*>& atree_2)

{
  int ii;
  ImprintLineSegment *prev_seg, *curr_seg;
  CubitBoolean first = CUBIT_TRUE;
  ImprintPointData *start_point;
  DLIList<ImprintMatchData*> match_data_list;
  ImprintMatchData *start_match;
  CubitStatus stat;
  
  for(ii = 0; ii < num_loops_1; ii++ )
  {
    prev_seg = surf_1_loop_heads[ii]->get_prev();
    first = CUBIT_TRUE;
    while(prev_seg != surf_1_loop_heads[ii]->get_prev() || first )
    {
      first = CUBIT_FALSE;
      curr_seg = prev_seg->get_next();
      if ( DEBUG_FLAG(129) )
      {
        draw_seg(curr_seg, CUBIT_GREEN);
      }
        //Setup the next itteration.
      prev_seg = curr_seg;
        //Now on to the algorithm
        //Just check the start point for each segment (so that
        //we don't do each point more than once!
      start_point = curr_seg->get_start();
      if ( start_point->get_matching_point() != NULL )
        continue;
        //Get the match data.
      match_data_list.clean_out();
      start_point->get_match_list(match_data_list);
      if ( match_data_list.size() == 0 )
        start_match = NULL;
      else if ( match_data_list.size() == 1 )
        start_match = match_data_list.get();
      else
      {
        PRINT_ERROR("Unresolved match data in final match section.  This\n"
                    "is a bug!. Please report it.\n");
        return CUBIT_FAILURE;
      }
      if ( start_match )
      {
        CubitBoolean ignore = CUBIT_FALSE;
        if ( ignore_match(start_point, start_match, atree_1, ignore)
             != CUBIT_SUCCESS )
          return CUBIT_FAILURE;
        if ( ignore )
          continue;
        if ( start_match->get_closest_point() != NULL )
        {
          stat = match_on_point(start_point, start_match,
                                surf_2_loop_heads, num_loops_2,
                                atree_2);
          if ( stat != CUBIT_SUCCESS )
            return CUBIT_FAILURE;
        }
        else if ( start_match->get_point_on() != NULL )
        {
            //split the segment with start point and update
            //the loops heads for surf_2 if it is changed.
          stat = match_on_segment(start_point, start_match,
                                  surf_2_loop_heads, num_loops_2,
                                  atree_2);
          if ( stat != CUBIT_SUCCESS )
            return stat;
        }
        else
        {
            //There really can't be anything else at this point.  Either
            //it matches with a node, or it matches on a segment.
          PRINT_ERROR("Problems with point matching logic in imprinting.\n");
          PRINT_ERROR("MatchData found without closest point or segment.\n");
          PRINT_ERROR("Internal Imprint Problem, Please Report.\n");
          return CUBIT_FAILURE;
        }
      }
    }
  }
  return CUBIT_SUCCESS;
}
CubitStatus ImprintBoundaryTool::match_on_point( ImprintPointData *start_point,
                                                 ImprintMatchData *start_match,
                                                 ImprintLineSegment **loop_heads,
                                                 int num_loops,
                                                 AbstractTree <ImprintLineSegment*>& atree_2)
{
    //This point matches with a point on the other boundary.  Find
    //the other point, and make sure that the other point, also
    //matches with start_point.  If it doesn't, it is that the tolerance
    //is big for the features and that we have points that match up
    //improperly.
  ImprintMatchData *temp_match;
  ImprintPointData *closest_point;
  DLIList <ImprintMatchData*> match_data_list;
  closest_point = start_match->get_closest_point();
  closest_point->get_match_list(match_data_list);
  if ( match_data_list.size() == 1 )
  {
    temp_match = match_data_list.get();
    if ( temp_match->get_closest_point() == start_point )
    {
        //Great.  We have a solid match.  Setup the data
        //structure.
      just_match_two_points(start_point, closest_point);
    }
    else
    {
        // start_point is point D.
        //Basically there can be one of two conditions.
        //1) A--------------------------------------B
        // C----------D------------------------------E
        //This condition, D is within tolerance of A, but
        //A's closest point is C.  To resolve this, we need
        //to match D with the segment A-B.
        //
        //2)                     
        //                          D
        //                         /
        //                        /
        //                      E-----------C-------
        //                    B----------A--------------F
        //For this condition, D matches with A, but A matches
        //with the  C.  This case we want to ignore.

      ImprintLineSegment *line_seg = temp_match->get_closest_seg();
      if ( line_seg == NULL )
      {
        PRINT_ERROR("Bad logic in ImprintBoundaryTool::match_point, bug...\n");
        assert(line_seg != NULL);
        return CUBIT_FAILURE;
      }
      ImprintPointData *point_D = start_point;
      //ImprintPointData *point_A = closest_point;
      ImprintPointData *point_C = temp_match->get_closest_point();
        //Find if D and C are on the same segment.  If they
        //are assume we have condition 1, else assume condition 2.
      ImprintLineSegment *prev_seg = line_seg->get_prev();
      ImprintLineSegment *next_seg = line_seg->get_next();
      CubitBoolean ignore = CUBIT_TRUE;
      if ( (point_D == prev_seg->get_start() ||
            point_D == prev_seg->get_end() ) &&
           (point_C == prev_seg->get_start() ||
            point_C == prev_seg->get_end() ) )
      {
        ignore = CUBIT_FALSE;
      }
      else if ((point_D == line_seg->get_start() ||
                point_D == line_seg->get_end() ) &&
               (point_C == line_seg->get_start() ||
                point_C == line_seg->get_end() ) )
      {
        ignore = CUBIT_FALSE;
      }
      else if ((point_D == next_seg->get_start() ||
                point_D == next_seg->get_end() ) &&
               (point_C == next_seg->get_start() ||
                point_C == next_seg->get_end() ) )
      {
        ignore = CUBIT_FALSE;
      }
//      if ( ignore )
//        start_point->set_unmatched();
      if ( !ignore )
      {
        line_seg = start_match->get_closest_seg();
        ImprintLineSegment *other_seg;
        if ( closest_point == line_seg->get_start() )
          other_seg = line_seg->get_prev();
        else
          other_seg = line_seg->get_next();
          //Find out which segment start_point is closest to.
        CubitVector close_1, close_2, start_v = start_point->coordinates();
        CubitBoolean on_line_seg = CUBIT_FALSE;
        CubitBoolean on_other_seg = CUBIT_FALSE;
        if ( closest_point_interior_seg(start_point,
                                        line_seg, close_1 ) )
        {
          if ( start_v.within_tolerance(close_1, myTolerance) )
            on_line_seg = CUBIT_TRUE;
        }
        if ( closest_point_interior_seg(start_point,
                                        other_seg, close_2 ) )
        {
          if ( start_v.within_tolerance(close_2, myTolerance) )
            on_other_seg = CUBIT_TRUE;
        }
        ImprintLineSegment *the_seg = NULL;
        CubitVector the_point;
        if ( on_line_seg && !on_other_seg )
        {
          the_seg = line_seg;
          the_point = close_1;
        }
        else if ( !on_line_seg && on_other_seg )
        {
          the_seg = other_seg;
          the_point = close_2;
        }
        else if ( on_line_seg && on_other_seg )
        {
          double dist_1 = (start_v-close_1).length_squared();
          double dist_2 = (start_v-close_2).length_squared();
          if ( dist_1 < dist_2 )
          {
            the_seg = line_seg;
            the_point = close_1;
          }
          else
          {
            the_seg = other_seg;
            the_point = close_2;
          }
        }
        if ( the_point.within_tolerance(closest_point->coordinates(),
                                        GEOMETRY_RESABS*500 ) )
        {
          PRINT_ERROR("The tolerance is too large for virtual imprinting\n"
                      "Surface %d and %d\n", refFacePtr1->id(),
                      refFacePtr2->id());
          return CUBIT_FAILURE;
        }
        start_match->set_closest_point(NULL);
        start_match->set_closest_seg(the_seg);
        start_match->set_point_on(the_point);
        return match_on_segment(start_point, start_match,
                                loop_heads, num_loops,
                                atree_2);
      }
    }
  }
  else
  {
    if ( match_data_list.size() == 0 )
      return CUBIT_SUCCESS;
    else
    {
      PRINT_ERROR("More than one match in final match(%d).\n",match_data_list.size() );
      PRINT_ERROR("This is a bug in virtual imprinting, try reducing\n"
                  "the tolerance.\n");
      return CUBIT_FAILURE;
    }
  }
  return CUBIT_SUCCESS;
}

CubitStatus ImprintBoundaryTool::match_on_segment( ImprintPointData *start_point,
                                                   ImprintMatchData *start_match,
                                                   ImprintLineSegment **loop_heads,
                                                   int num_loops,
                                                   AbstractTree <ImprintLineSegment*>& atree_2)
{
    //Closest point is on a segment on the other boundary.
    //Split that segment at the closest point.
    //Adjust the data...
  ImprintLineSegment *new_seg_1,*new_seg_2;
  
  ImprintLineSegment* closest_seg = start_match->get_closest_seg();
  CubitVector closest_vec(*(start_match->get_point_on()));
  ImprintPointData *new_point = new ImprintPointData(start_point,
                                                     closest_vec);
  allocatedPointData->append(new_point);
  new_point->owner(closest_seg->owner());
  if ( start_point->is_owner_vertex() )
  {
    new_point->set_point_type(CREATE_NEW_VERTEX);
    if ( start_point->get_point_type() != CREATE_NEW_VERTEX )
      start_point->set_point_type(VERTEX_ON_BOTH_BOUNDARIES);
  }
  else
  {
    new_point->set_point_type(ON_BOTH_BOUNDARIES);
    start_point->set_point_type(ON_BOTH_BOUNDARIES);
  }

  new_seg_1 = new ImprintLineSegment( closest_seg->get_start(),
                                      new_point,
                                      closest_seg->owner());
  new_seg_2 = new ImprintLineSegment( new_point,
                                      closest_seg->get_end(),
                                      closest_seg->owner());
  allocatedLineData->append(new_seg_1);
  allocatedLineData->append(new_seg_2);
    //update the linked list structure to leave out the closest_seg
    //but include the new segs instead.
  update_linked_list(closest_seg, new_seg_1, new_seg_2);
    //update the ImprintMatchData that reference the closest_seg
    //to reference new_seg_1 or new_seg_2.
  CubitStatus stat = update_seg_matches(closest_seg, new_seg_1, new_seg_2, start_match);
  if ( stat != CUBIT_SUCCESS )
    return stat;
  closest_seg->set_inactive(CUBIT_TRUE);
  int jj;
  for (jj = 0; jj < num_loops; jj++ )
  {
    if ( loop_heads[jj] == closest_seg )
      loop_heads[jj] = new_seg_1;
  }
    //update the atree.
  atree_2.remove(closest_seg);
  atree_2.add(new_seg_1);
  atree_2.add(new_seg_2);
  return CUBIT_SUCCESS;

}

CubitStatus ImprintBoundaryTool::update_seg_matches(ImprintLineSegment *old_seg,
                                                    ImprintLineSegment *new_seg_1,
                                                    ImprintLineSegment *new_seg_2,
                                                    ImprintMatchData *curr_match_data)
{
  DLIList <ImprintMatchData*> seg_match_data;
  old_seg->get_match_data(seg_match_data);
  if ( seg_match_data.size() == 0 )
    return CUBIT_SUCCESS;
    //get the equation of the line.
  CubitVector split_coord = new_seg_1->get_end()->coordinates();
  CubitVector start_v = old_seg->get_start()->coordinates();
  CubitVector end_v = old_seg->get_end()->coordinates();
  CubitVector diff_end_start = end_v - start_v;
  double dot1 = diff_end_start%diff_end_start;
  if ( dot1 > -CUBIT_RESABS && dot1 < CUBIT_RESABS )
  {
      //This is bad.  It means that the start and end
      //vectors are really close together.
    PRINT_ERROR("Segment used for imprinting is corrupted.\n");
    PRINT_ERROR("Internal Imprint Problem, Please Report.\n");
    return CUBIT_FAILURE;
  }
  CubitVector diff_split_start = split_coord - start_v;
  double split_t = (diff_split_start%diff_end_start)/dot1;
  
  int ii;
  ImprintMatchData *match_data;
  CubitVector *closest_point, diff_curr_start, temp_v;
  double curr_t;
  for ( ii = seg_match_data.size(); ii > 0; ii-- )
  {
    match_data = seg_match_data.get_and_step();
    if ( match_data == curr_match_data )
      continue;
    closest_point = match_data->get_point_on();
    if (closest_point == NULL )
    {
	  temp_v = match_data->get_closest_point()->coordinates();
      closest_point = &temp_v;
    }
      //Find the parametric position on the old_seg for this
      //point.
    diff_curr_start = (*closest_point) - start_v;
    curr_t = (diff_curr_start%diff_end_start)/dot1;
    if ( curr_t < split_t )
    {
      match_data->set_closest_seg(new_seg_1);
      new_seg_1->add_match_data(match_data);
    }
    else if (curr_t > split_t )
    {
      match_data->set_closest_seg(new_seg_2);
      new_seg_2->add_match_data(match_data);
    }
    else
    {
        //hmm two match datas at the same point.  That is
        //an error.
      PRINT_ERROR("Bad logic in match data.\n");
      PRINT_ERROR("Internal Imprint Problem, Please Report.\n");
      return CUBIT_FAILURE;
    }
  }
  return CUBIT_SUCCESS;
}
  
void ImprintBoundaryTool::update_linked_list(ImprintLineSegment *old_seg,
                                             ImprintLineSegment *new_seg_1,
                                             ImprintLineSegment *new_seg_2 )
{
  new_seg_1->set_prev(old_seg->get_prev());
  old_seg->get_prev()->set_next(new_seg_1);
  new_seg_1->set_next(new_seg_2);
  new_seg_2->set_prev(new_seg_1);
  new_seg_2->set_next(old_seg->get_next());
  old_seg->get_next()->set_prev(new_seg_2);
}

void ImprintBoundaryTool::just_match_two_points(ImprintPointData *point_1,
                                                ImprintPointData *point_2 )
{
    //Okay, we need to match up point_2 and
    //point_1.          
    //Set the actual matching point for the imprinting
    //algorithm
  point_1->set_matching_point(point_2);
  point_2->set_matching_point(point_1);
    //Clear out the possible match lists.
  if ( point_1->id() == 681 ||
       point_1->id() == 150 )
  {
    int temp = 0;
    temp++;
  }
  if ( point_1->is_owner_vertex() &&
       !point_2->is_owner_vertex() )
  {
    point_1->set_point_type(VERTEX_ON_BOTH_BOUNDARIES);
    point_2->set_point_type(CREATE_NEW_VERTEX);
  }
  else if(!point_1->is_owner_vertex() &&
          point_2->is_owner_vertex() )
  {
    point_1->set_point_type(CREATE_NEW_VERTEX);
    point_2->set_point_type(VERTEX_ON_BOTH_BOUNDARIES);
  }
  else if ( point_1->is_owner_vertex() &&
            point_2->is_owner_vertex() )
  {
    point_1->set_point_type(VERTEX_ON_BOTH_BOUNDARIES);
    point_2->set_point_type(VERTEX_ON_BOTH_BOUNDARIES);
  }
  else // Both aren't vertices...
  {
    point_1->set_point_type(ON_BOTH_BOUNDARIES);
    point_2->set_point_type(ON_BOTH_BOUNDARIES);
  }
}
CubitStatus ImprintBoundaryTool::find_crossings( ImprintLineSegment **surf_1_loop_heads,
                                                 int num_loops_1,
                                                 ImprintLineSegment **surf_2_loop_heads,
                                                 int num_loops_2,
                                                 AbstractTree <ImprintLineSegment*>& atree )
{
  int ii, jj;
  ImprintLineSegment *curr_seg, *prev_seg, *test_seg;
  DLIList <ImprintLineSegment*> close_segments;
  CubitBoolean first = CUBIT_TRUE;
  CubitBox curr_box;
  ImprintPointData *imp_point_0, *imp_point_1;
  ImprintPointData *imp_point_2, *imp_point_3;
  
  for(ii = 0; ii < num_loops_1; ii++ )
  {
    prev_seg = surf_1_loop_heads[ii]->get_prev();
    first = CUBIT_TRUE;
    while(prev_seg != surf_1_loop_heads[ii]->get_prev() || first )
    {
      first = CUBIT_FALSE;
      curr_seg = prev_seg->get_next();
      if ( DEBUG_FLAG(129) )
      {
        draw_seg(curr_seg, CUBIT_GREEN);
      }
        //Setup the next itteration.
      prev_seg = curr_seg;
      curr_box = curr_seg->bounding_box();
      close_segments.clean_out();
      atree.find(curr_box, close_segments);
      imp_point_0 = curr_seg->get_start();
      imp_point_1 = curr_seg->get_end();
      for ( jj = 0; jj < close_segments.size(); jj++ )
      {
        test_seg = close_segments.pop();
        imp_point_2 = test_seg->get_start();
        imp_point_3 = test_seg->get_end();
          //don't test if these already match.
        if ( imp_point_0->get_matching_point() == imp_point_2 ||
             imp_point_0->get_matching_point() == imp_point_3 ||
             imp_point_1->get_matching_point() == imp_point_2 ||
             imp_point_1->get_matching_point() == imp_point_3 )
          continue;
        CubitBoolean modified = CUBIT_FALSE;
          //Test and resolve.  Note that only one of the atrees
          //is getting updated (this atree is for surface 2.  )
          //The loop heads are updated, and curr_seg will
          //get modified if the a crossing intersect occurs.
        CubitStatus stat = do_crossing( imp_point_0, imp_point_1,
                                        imp_point_2, imp_point_3,
                                        curr_seg, test_seg,
                                        surf_1_loop_heads, num_loops_1,
                                        surf_2_loop_heads, num_loops_2,
                                        atree, modified );
        if ( stat != CUBIT_SUCCESS )
          return CUBIT_FAILURE;
        if ( modified )
        {
            //curr_seg was modified in the do_crossing (the original
            //was split.).  Make sure we start again to check for
            //crossings with this item.  To do that, set prev seg
            //to it's previous segment.  The result is that the next
            //time through it will be curr_seg.
          prev_seg = curr_seg->get_prev();
        }
      }
    }
  }
  return CUBIT_SUCCESS;
}
CubitStatus ImprintBoundaryTool::do_crossing( ImprintPointData *imp_point_0,
                                              ImprintPointData *imp_point_1,
                                              ImprintPointData *imp_point_2,
                                              ImprintPointData *imp_point_3,
                                              ImprintLineSegment *&seg_1,
                                              ImprintLineSegment *seg_2,
                                              ImprintLineSegment **surf_1_loop_heads,
                                              int num_loops_1,
                                              ImprintLineSegment **surf_2_loop_heads,
                                              int num_loops_2,
                                              AbstractTree <ImprintLineSegment*>& atree,
                                              CubitBoolean &modified)
{
                                              
    //Okay, Test for the cross intersection.  If this isnt'
    //within tolerance, then these two segments couldn't
    //posibly intersect!
  IntersectionTool int_tool(GEOMETRY_RESABS); //use GEOMETRY_RESABS
    //rather than myTolerance for calculations...
  CubitVector point_0 = imp_point_0->coordinates();
  CubitVector point_1 = imp_point_1->coordinates();
  CubitVector point_2 = imp_point_2->coordinates();
  CubitVector point_3 = imp_point_3->coordinates();
  CubitVector closest_point_seg_1, closest_point_seg_2;
  double sc, tc;
  CubitStatus stat = int_tool.closest_points_on_segments(point_0, point_1,
                                                         point_2, point_3,
                                                         closest_point_seg_1,
                                                         closest_point_seg_2,
                                                         sc, tc);
  if (stat != CUBIT_SUCCESS )
  {
    PRINT_ERROR("Problems calculation closest points on "
                "segments for boundary imprinting.\n");
    return CUBIT_FAILURE;
  }
    //Make sure the closest points aren't the end points.  If they are
    //and we are within tolerance, it may be that the tolerance is too big
    //cause we shouldn't be at a cross if the closest point is an end point...
  if ( sc > 0. && sc < 1. && tc > 0. && tc < 1. &&
       closest_point_seg_1.within_tolerance(closest_point_seg_2,myTolerance) )
  {
    ImprintLineSegment *new_segments[4];
      //okay, these guys actually do intersect!
    PRINT_DEBUG_129("Found CROSS_INTERSECT\n");
    int debug = 0;
    if ( debug )
    {
      GfxDebug::clear();
      draw_point(imp_point_0);
      draw_point(imp_point_1);
      draw_point(imp_point_2);
      draw_point(imp_point_3);
      draw_seg(seg_1, CUBIT_RED);
      draw_seg(seg_2, CUBIT_YELLOW);
      GfxDebug::mouse_xforms();
    }
    ImprintPointData *new_point_seg_1, *new_point_seg_2;
    new_point_seg_1 = new ImprintPointData(closest_point_seg_1.x(),
                                           closest_point_seg_1.y(),
                                           closest_point_seg_1.z());
    allocatedPointData->append(new_point_seg_1);
    new_point_seg_1->owner(seg_1->owner());
    new_point_seg_1->set_point_type(CREATE_NEW_VERTEX);
    new_point_seg_2 = new ImprintPointData(closest_point_seg_2.x(),
                                           closest_point_seg_2.y(),
                                           closest_point_seg_2.z());
    allocatedPointData->append(new_point_seg_2);
    new_point_seg_2->owner(seg_2->owner());
    new_point_seg_2->set_point_type(CREATE_NEW_VERTEX);

    new_point_seg_2->set_matching_point(new_point_seg_1);
    new_point_seg_1->set_matching_point(new_point_seg_2);
    
    new_segments[0] = new ImprintLineSegment(imp_point_0, new_point_seg_1,
                                             seg_1->owner());
    new_segments[1] = new ImprintLineSegment(new_point_seg_1, imp_point_1,
                                             seg_1->owner());
    new_segments[2] = new ImprintLineSegment(imp_point_2, new_point_seg_2,
                                             seg_2->owner());
    new_segments[3] = new ImprintLineSegment(new_point_seg_2, imp_point_3,
                                             seg_2->owner());
    allocatedLineData->append(new_segments[0]);
    allocatedLineData->append(new_segments[1]);
    allocatedLineData->append(new_segments[2]);
    allocatedLineData->append(new_segments[3]);
      //now update the lists.
      //update list 1.
    update_list(new_segments, seg_1, seg_2, CUBIT_TRUE);
      //update list 2.
    update_list(new_segments, seg_1, seg_2, CUBIT_FALSE);
      //seg_1 and seg_2 will no longer be used
    seg_1->set_inactive(CUBIT_TRUE);
    seg_2->set_inactive(CUBIT_TRUE);
      //update the atree for the second surface.
    atree.remove(seg_2);
    atree.add(new_segments[2]);
    atree.add(new_segments[3]);
      //update the loop heads
    int jj;
    for (jj = 0; jj < num_loops_1;jj++ )
    {
      if ( surf_1_loop_heads[jj] == seg_1 )
        surf_1_loop_heads[jj] = new_segments[0];
    }
    for (jj = 0; jj < num_loops_2;jj++ )
    {
      if ( surf_2_loop_heads[jj] == seg_2 )
        surf_2_loop_heads[jj] = new_segments[2];
    }
      //now set seg_1 equal to new_segments[0]
      //so that we continue looping properly.
    seg_1 = new_segments[0];
    modified = CUBIT_TRUE;
  }
  
  return CUBIT_SUCCESS;
}
void ImprintBoundaryTool::update_boundary_loops( PointLoopList &boundary_loops,
                                                 ImprintLineSegment **surf_loop_heads)
{
  int ii;
  PointList *point_loop;
  ImprintLineSegment *prev_seg, *seg;
  
  CubitBoolean first = CUBIT_TRUE;
  boundary_loops.reset();
  for ( ii = 0; ii < boundary_loops.size(); ii++ )
  {
    point_loop = boundary_loops.get_and_step();
    point_loop->clean_out();
      //now from the linked list get the line loop.
    prev_seg = surf_loop_heads[ii]->get_prev();
    first = CUBIT_TRUE;
    while( prev_seg != surf_loop_heads[ii]->get_prev() || first )
    {
      first = CUBIT_FALSE;
      seg = prev_seg->get_next();
      point_loop->append(seg->get_end());
      if ( DEBUG_FLAG(129) )
        draw_end(seg);
      prev_seg = seg;
    }
  }
}

CubitStatus ImprintBoundaryTool::ignore_match( ImprintPointData *start_point,
                                               ImprintMatchData *start_match,
                                               AbstractTree <ImprintLineSegment*>& atree_1,
                                               CubitBoolean &ignore)
{
  ignore = CUBIT_TRUE;
    //assume that start_point is in atree_1 (here 1 and 2 don't neccesaryily mean
    //loopFacePtr1 and loopFacePtr2 because this is called for both faces.

    //First determine what type of match, point to point or point to segment.
  if ( start_match->get_closest_point() != NULL )
  {
     ImprintMatchData *temp_match;
     ImprintPointData *closest_point;
     DLIList <ImprintMatchData*> match_data_list;
     closest_point = start_match->get_closest_point();
     closest_point->get_match_list(match_data_list);
     if ( match_data_list.size() == 1 )
     {
       temp_match = match_data_list.get();
         //if they match don't ignore this.
       if ( temp_match->get_closest_point() == start_point )
       {
         ignore = CUBIT_FALSE;
         return CUBIT_SUCCESS;
       }
       else
       {
           //Now we have a conflict.  Start matches with closest_point,
           //but closest_point doesn't match with start.

           //Find out if there is a segment in atree_1 that is closer
           //to closest_point.  If there isn't, then we need to
           //match start with one of the two linesegments attached to closest
           //point.
           //If there is a closer line segment to closest_point, then don't
           //match.

           //First, define the search area as the area around the start_point.
         ImprintLineSegment *start_prev_seg = start_point->get_prev_seg();
         ImprintLineSegment *start_next_seg = start_point->get_next_seg();
         CubitBox curr_box = start_prev_seg->bounding_box();
         curr_box |= start_next_seg->bounding_box();
         DLIList<ImprintLineSegment*> closest_segments;
         atree_1.find(curr_box, closest_segments);
           //Find the segment that is closest to closest_point.
         ImprintLineSegment *closest_seg = NULL, *close_seg;
         double dist;
         double min_dist = CUBIT_DBL_MAX;
         int ii;
         CubitVector closest_vec, closest_point_vec = closest_point->coordinates();
         for ( ii = closest_segments.size(); ii > 0; ii-- )
         {
           close_seg = closest_segments.pop();
           if ( !closest_point_seg(closest_point, close_seg,
                                   closest_vec) )
             return CUBIT_FAILURE;
           dist = (closest_vec - closest_point_vec).length_squared();
           if ( dist < min_dist )
           {
             min_dist = dist;
             closest_seg = close_seg;
           }
         }
         if ( closest_seg == NULL )
         {
           PRINT_ERROR("Couldn't find a closest segment.(bug)\n");
           assert(closest_seg != NULL);
           return CUBIT_FAILURE;
         }
         if ( closest_seg != start_prev_seg ||
              closest_seg != start_next_seg )
         {
             //okay this we need to ignore.  Just say start point doesn't match with anything.
           start_point->set_unmatched();
           start_point->set_matching_point(NULL);
           ignore = CUBIT_TRUE;
           return CUBIT_SUCCESS;
         }
         else
         {
             //okay, now we actually want to match start point with one of the segments
             //attached to closest_point.
           ImprintLineSegment *next_close_seg, *prev_close_seg;
           next_close_seg = closest_point->get_next_seg();
           prev_close_seg = closest_point->get_prev_seg();
             //determine which is closer...
           CubitVector closest_next, closest_prev, start_vec = start_point->coordinates();
           if ( !closest_point_seg(start_point, next_close_seg,
                                   closest_next) )
             return CUBIT_FAILURE;
           if ( !closest_point_seg(start_point, prev_close_seg,
                                   closest_prev) )
             return CUBIT_FAILURE;
           double dist_prev, dist_next;
           CubitVector point_on;
           dist_next = (start_vec - closest_next).length_squared();
           dist_prev = (start_vec - closest_prev).length_squared();
           closest_seg = NULL;
           if ( dist_next < dist_prev )
           {
             closest_seg = next_close_seg;
             min_dist = dist_next;
             point_on = closest_next;
           }
           else
           {
             closest_seg = prev_close_seg;
             min_dist = dist_prev;
             point_on = closest_prev;
           }

            //make sure the point_on isn't the same location as the closest_point.
           CubitVector closest_point_vec = closest_point->coordinates();
           if ( closest_point_vec.within_tolerance(point_on, myTolerance/10.0) )
           {
             ignore = CUBIT_TRUE;
             start_point->set_unmatched();
             start_point->set_matching_point(NULL);
             return CUBIT_SUCCESS;
           }
             //make sure point_on and start point are within tolerances.
           dist = (point_on - start_vec).length_squared();
           if ( dist > myTolerance*myTolerance )
           {
             PRINT_ERROR("Bad logic in ignore_match in ImprintImprintTool.\n");
             assert(dist <= myTolerance*myTolerance);
             return CUBIT_FAILURE;
           }
           start_match->set_closest_point(NULL);
           start_point->set_matching_point(NULL);
           start_match->set_closest_dist(min_dist);
           start_match->set_point_on(point_on);
           start_match->set_closest_seg(closest_seg);
           ignore = CUBIT_FALSE;
           return CUBIT_SUCCESS;
         }
       }
     }
     else
     {
       if ( match_data_list.size() == 0 )
       {
         ignore = CUBIT_TRUE;
         return CUBIT_SUCCESS;
       }
       else
       {
         PRINT_ERROR("More than one match in final match(%d).\n",match_data_list.size() );
         PRINT_ERROR("This is a bug in loop imprinting, try reducing\n"
                     "the size for decomposition.\n");
         return CUBIT_FAILURE;
       }
     }
  }
  else if ( start_match->get_point_on() != NULL )
  {
    if ( start_point->get_matching_point() )
    {
      ImprintPointData *closest_point = start_point->get_matching_point();
      if ( closest_point->get_matching_point() == start_point )
      {
        ignore = CUBIT_TRUE;
        return CUBIT_SUCCESS;
      }
    }
      //this must mean we match on a segment.

      //Determine if there is a segment on 1, that is closer to this
      //point than start_point is.  If there is, ignore this match.
    CubitVector point_on;
    point_on.set(*(start_match->get_point_on()));
    ImprintLineSegment *start_prev_seg = start_point->get_prev_seg();
    ImprintLineSegment *start_next_seg = start_point->get_next_seg();
    CubitBox curr_box = start_prev_seg->bounding_box();
    curr_box |= start_next_seg->bounding_box();
    DLIList<ImprintLineSegment*> closest_segments;
    atree_1.find(curr_box, closest_segments);
      //Find the segment that is closest to closest_point.
    ImprintLineSegment *closest_seg = NULL, *close_seg;
    double dist;
    double min_dist = CUBIT_DBL_MAX;
    int ii;
    //dflush();
    CubitVector closest_vec;
    //CubitVector start_vec = start_point->coordinates();
   // dclear();
   // dseg(start_prev_seg);
   // dseg(start_next_seg);
    //dvec(point_on);
    //dmouse();
    for ( ii = closest_segments.size(); ii > 0; ii-- )
    {
      close_seg = closest_segments.pop();
      if ( !closest_point_seg(point_on, close_seg,
                              closest_vec) )
        return CUBIT_FAILURE;
      dist = (closest_vec - point_on).length_squared();
      if ( dist < min_dist )
      {
        min_dist = dist;
        closest_seg = close_seg;
      }
    }
    if ( closest_seg == NULL )
    {
      PRINT_ERROR("Couldn't find a closest segment.(bug)\n");
      assert(closest_seg != NULL);
      return CUBIT_FAILURE;
    }
    if ( closest_seg == start_prev_seg ||
         closest_seg == start_next_seg )
    {
      ignore = CUBIT_FALSE;
      return CUBIT_SUCCESS;
    }
    else
    {
      ignore = CUBIT_TRUE;
      start_point->set_unmatched();
      start_point->set_matching_point(NULL);
      return CUBIT_SUCCESS;
    }
    
  }
  
  return CUBIT_FAILURE;
}

//---------------------------------------------------------------------
// Private Function: create_virtual_vertex
// Creates a facet based vertex.
//---------------------------------------------------------------------
RefVertex* ImprintBoundaryTool::create_virtual_vertex(CubitVector &pos)
{
    //Create a facet point first.
  TBPoint *new_point;
  if ( FacetModifyEngine::instance()->make_facet_point(pos, new_point)
       != CUBIT_SUCCESS )
  {
    PRINT_ERROR("Vertex creation failed.\n");
    return (RefVertex*)NULL;
  }
  
  return GeometryQueryTool::instance()->make_free_RefVertex(new_point);
}
//---------------------------------------------------------------------
// Private Function: create_virtual_edge
// Creates a facet based vertex.
//---------------------------------------------------------------------
RefEdge* ImprintBoundaryTool::create_virtual_edge(RefVertex *start,
                                                  RefVertex *end,
                                                  DLIList<CubitVector*> &interior_points)
{
    //Create CubitPoints for all of these locations.
  DLIList<CubitPoint*> point_list;
  CubitPoint *start_p = (CubitPoint*) new CubitPointData(start->coordinates());
  CubitPoint *end_p = (CubitPoint*) new CubitPointData(end->coordinates());
  CubitPoint *new_point, *prev_point;
  int ii;
  CubitVector *curr_pos;
  interior_points.reset();
  CubitFacetEdge *new_edge;
  prev_point = start_p;
  DLIList<CubitFacetEdge*> edge_list;
  point_list.append(start_p);
  for ( ii = 0; ii < interior_points.size(); ii++ )
  {
    curr_pos = interior_points.get_and_step();
    new_point = (CubitPoint*) new CubitPointData(*curr_pos);
    point_list.append(new_point);
    new_edge = (CubitFacetEdge*) new CubitFacetEdgeData(prev_point,
                                                        new_point);
    edge_list.append(new_edge);
    prev_point = new_point;
  }
    //create the last edge.
  new_edge = (CubitFacetEdge*) new CubitFacetEdgeData(prev_point,
                                                      end_p);
  point_list.append(end_p);
  edge_list.append(new_edge);
  CurveFacetEvalTool *curve_facet_tool = new CurveFacetEvalTool;
  CubitStatus status = curve_facet_tool->initialize( edge_list, point_list );
  if ( status != CUBIT_SUCCESS )
  {
      PRINT_ERROR("Error in creating new facet curve.\n");
      return NULL;
  }

  bool start_vertex_is_free = true;
  if( start->num_parent_ref_entities() )
    start_vertex_is_free = false;

  bool end_vertex_is_free = true;
  if( end->num_parent_ref_entities() )
    end_vertex_is_free = false;

  DLIList<CoEdgeSM*> coedgesms;
  GeometryEntity *ge = start->get_geometry_entity_ptr();
  TBPoint *start_psm = CAST_TO(ge, TBPoint);
  ge = end->get_geometry_entity_ptr();
  TBPoint *end_psm = CAST_TO(ge, TBPoint);
  assert(start_psm && end_psm);
  FacetCurve *new_facet_curve = new FacetCurve(curve_facet_tool,
                                               start_psm,
                                               end_psm, coedgesms);
  if ( new_facet_curve == NULL )
  {
    PRINT_ERROR("Error in creating new facet curve.\n");
    assert(new_facet_curve != NULL);
    return (RefEdge*)NULL;
  }
  
  //if vertices are consumed....notify that they are gone.
  if( start_vertex_is_free )
    start->notify_all_observers( TOP_LEVEL_ENTITY_DESTRUCTED );
  if( end_vertex_is_free )
    end->notify_all_observers( TOP_LEVEL_ENTITY_DESTRUCTED );
  
  Curve *new_curve_ptr = (Curve*) new_facet_curve;
  return GeometryQueryTool::instance()->make_RefEdge(new_curve_ptr);
}


  
  
  

  
  
                                                       
