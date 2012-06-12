//----------------------------------------------------------
//Class: GeomMeasureTool
//Description: Does Geometric measureing and statistical
//             gathering on various topologies.
//Creator: David R. White
//Date: 7/9/2002
//----------------------------------------------------------


#include "GeomMeasureTool.hpp"
#include "RefEntity.hpp"
#include "RefVertex.hpp"
#include "RefEdge.hpp"
#include "RefFace.hpp"
#include "RefVolume.hpp"
#include "Body.hpp"
#include "RefGroup.hpp"
#include "RefEntityFactory.hpp"
#include "CoEdge.hpp"
#include "Loop.hpp"
#include "Shell.hpp"
#include "Chain.hpp"
#include "CastTo.hpp"
#include "GMem.hpp"
#include "Curve.hpp"
#include "Surface.hpp"
#include "GeometryQueryEngine.hpp"
#include "GeomSeg.hpp"
#include "GeomPoint.hpp"
#include "RTree.hpp"
#include "AbstractTree.hpp"
#include "IntersectionTool.hpp"
#include "CpuTimer.hpp"
#include "GeometryQueryTool.hpp"
#include "GeometryModifyTool.hpp"
#include "MergeTool.hpp"
#include "GeomDataObserver.hpp"
#include "Lump.hpp"
#include "ModelQueryEngine.hpp"
#include "ProgressTool.hpp"
#include "AppUtil.hpp"

#include <set>
#include <map>
#include <stack>
#include <algorithm>

void GeomMeasureTool::get_edges_from_list(DLIList <RefEntity*> &entity_list,
                                          DLIList <RefEdge*> &ref_edges )
{
    //assume the mark is zero.
  int ii,jj;
  RefEntity *entity_ptr;
  RefEdge *ref_edge;
  
  TopologyEntity *topo_ent;
  DLIList <RefEdge*> tmp_edges;
  for ( ii = entity_list.size(); ii > 0; ii-- )
  {
    tmp_edges.clean_out();
    entity_ptr = entity_list.get_and_step();
    topo_ent = CAST_TO(entity_ptr, TopologyEntity);
    if ( topo_ent == NULL )
      continue;
    topo_ent->ref_edges(tmp_edges);
    for ( jj = tmp_edges.size(); jj > 0; jj-- )
    {
      ref_edge = tmp_edges.pop();
      if ( !ref_edge->marked() )
      {
        ref_edge->marked(1);
        ref_edges.append(ref_edge);
      }
    }
  }
  for ( ii = ref_edges.size(); ii > 0; ii-- )
    ref_edges.get_and_step()->marked(0);
}

void GeomMeasureTool::get_bodies_from_list(DLIList <RefVolume*> &entity_list,
                                           DLIList <Body*> &ref_bodies )
{
    //assume the mark is zero.
  int ii, jj;
  RefVolume *entity_ptr;
  Body *ref_body;
  
  TopologyEntity *topo_ent;
  DLIList <Body*> tmp_bodies;
  for ( ii = entity_list.size(); ii > 0; ii-- )
  {
    tmp_bodies.clean_out();
    entity_ptr = entity_list.get_and_step();
    topo_ent = CAST_TO(entity_ptr, TopologyEntity);
    if ( topo_ent == NULL )
      continue;
    topo_ent->bodies(tmp_bodies);
    for ( jj = tmp_bodies.size(); jj > 0; jj-- )
    {
      ref_body = tmp_bodies.pop();
      if ( !ref_body->marked() )
      {
        ref_body->marked(1);
        ref_bodies.append(ref_body);
      }
    }
  }
  for ( ii = ref_bodies.size(); ii > 0; ii-- )
    ref_bodies.get_and_step()->marked(0);
}

void GeomMeasureTool::get_faces_from_list(DLIList <RefEntity*> &entity_list,
                                          DLIList <RefFace*> &ref_faces )
{
    //assume the mark is zero.
  int ii,jj;
  RefEntity *entity_ptr;
  RefFace *ref_face;
  
  TopologyEntity *topo_ent;
  DLIList <RefFace*> tmp_faces;
  for ( ii = entity_list.size(); ii > 0; ii-- )
  {
    tmp_faces.clean_out();
    entity_ptr = entity_list.get_and_step();
    topo_ent = CAST_TO(entity_ptr, TopologyEntity);
    topo_ent->ref_faces(tmp_faces);
    for ( jj = tmp_faces.size(); jj > 0; jj-- )
    {
      ref_face = tmp_faces.pop();
      if ( !ref_face->marked() )
      {
        ref_face->marked(1);
        ref_faces.append(ref_face);
      }
    }
  }
  for ( ii = ref_faces.size(); ii > 0; ii-- )
    ref_faces.get_and_step()->marked(0);
}
void GeomMeasureTool::get_volumes_from_list(DLIList <RefEntity*> &entity_list,
                                            DLIList <RefVolume*> &ref_volumes )
{
    //assume the mark is zero.
  int ii,jj;
  RefEntity *entity_ptr;
  RefVolume *ref_volume;
  
  TopologyEntity *topo_ent;
  DLIList <RefVolume*> tmp_volumes;
  for ( ii = entity_list.size(); ii > 0; ii-- )
  {
    tmp_volumes.clean_out();
    entity_ptr = entity_list.get_and_step();
    topo_ent = CAST_TO(entity_ptr, TopologyEntity);
    topo_ent->ref_volumes(tmp_volumes);
    for ( jj = tmp_volumes.size(); jj > 0; jj-- )
    {
      ref_volume = tmp_volumes.pop();
      if ( !ref_volume->marked() )
      {
        ref_volume->marked(1);
        ref_volumes.append(ref_volume);
      }
    }
  }
  for ( ii = ref_volumes.size(); ii > 0; ii-- )
    ref_volumes.get_and_step()->marked(0);
}
void GeomMeasureTool::measure_curve_length(DLIList <RefEntity*> &entity_list,
                                           double &smallest,
                                           RefEdge *&smallest_edge,
                                           double &largest,
                                           RefEdge *&largest_edge,
                                           double &average,
                                           double &sum)
{
  DLIList<RefEdge*> ref_edges;
  get_edges_from_list(entity_list, ref_edges);
  measure_curve_length(ref_edges, smallest, smallest_edge,
                       largest, largest_edge, average, sum);
  
}

void GeomMeasureTool::measure_curve_length(DLIList <RefEdge*> &ref_edges,
                                           double &smallest,
                                           RefEdge *&smallest_edge,
                                           double &largest,
                                           RefEdge *&largest_edge,
                                           double &average,
                                           double &sum)
{
  int ii;
  double curr_length, total_length = 0.0;
  smallest = CUBIT_DBL_MAX;
  largest = 0.0;
  RefEdge *curr_edge;
  for ( ii = ref_edges.size(); ii > 0; ii-- )
  {
    curr_edge = ref_edges.get_and_step();
    if( curr_edge->geometry_type() == POINT_CURVE_TYPE )
      continue;
    curr_length = curr_edge->measure();
    if ( curr_length < smallest )
    {
      smallest = curr_length;
      smallest_edge = curr_edge;
    }
    if ( curr_length > largest )
    {
      largest = curr_length;
      largest_edge = curr_edge;
    }
    total_length += curr_length;
  }
  if ( ref_edges.size() != 0 )
    average = total_length/ref_edges.size();
  sum = total_length;
  
}
void GeomMeasureTool::measure_face_curves( RefFace *ref_face,
                                           double &min_curve_length,
                                           double &max_curve_ratio,
                                           RefEdge *&min_ref_edge)
{
  DLIList<DLIList<CoEdge*>*> co_edge_loops;
  DLIList <CoEdge*> *co_edge_list;
  ref_face->co_edge_loops(co_edge_loops);
  CoEdge *curr_co_edge, *next_co_edge;
  min_curve_length = CUBIT_DBL_MAX;
  double curve_length, next_curve_length, ratio;
  max_curve_ratio = -CUBIT_DBL_MAX;
  CubitBoolean ratio_set = CUBIT_FALSE;
  int ii, jj;
  for ( ii = co_edge_loops.size(); ii > 0; ii-- )
  {
    co_edge_list = co_edge_loops.get_and_step();
    curve_length = co_edge_list->get()->get_ref_edge_ptr()->measure();
    for ( jj = co_edge_list->size(); jj > 0; jj-- )
    {
      curr_co_edge = co_edge_list->get_and_step();
      next_co_edge = co_edge_list->get();
      next_curve_length = next_co_edge->get_ref_edge_ptr()->measure();
        //find the smallest curve
      if ( curve_length < min_curve_length )
      {
        min_curve_length = curve_length;
        min_ref_edge = curr_co_edge->get_ref_edge_ptr();
      }
      if ( co_edge_list->size() > 1 )
      {
          //find the biggest change in curve lengths.
        if ( curve_length > next_curve_length )
        {
          if ( next_curve_length < CUBIT_RESABS )
            ratio = CUBIT_DBL_MAX;
          else
            ratio = curve_length / next_curve_length;
        }
        else
        {
          if ( curve_length < CUBIT_RESABS )
            ratio = CUBIT_DBL_MAX;
          else
            ratio = next_curve_length /curve_length;
        }
        curve_length = next_curve_length;
        if ( ratio > max_curve_ratio )
        {
          ratio_set = CUBIT_TRUE;
          max_curve_ratio = ratio;
        }
      }
    }
  }
  if ( !ratio_set )
    max_curve_ratio = -1.0;
}

        
  
  
CubitStatus GeomMeasureTool::measure_face_loops( RefFace *face,
                                                 double &min_distance_between_loops,
                                                 double &min_distance_in_one_loop,
                                                 double &min_angle,
                                                 double &max_angle,
                                                 double tolerance)
{
    //This is the most tricky of these functions.
    //To do this I'm going to facet the edges, then use an AbstractTree to
    //find the minimum distance between the loops and between a single
    //loop.
  PointLoopList boundary_point_loops;
  CubitStatus stat;
  stat = get_boundary_points(face, boundary_point_loops, tolerance);
  if ( stat != CUBIT_SUCCESS )
    return CUBIT_FAILURE;
  
  SegLoopList boundary_seg_loops;
  stat = convert_to_lines( boundary_point_loops,
                           boundary_seg_loops,
                           face);
  if ( stat != CUBIT_SUCCESS )
    return stat;

    //Now add the points to the AbstractTree.
  DLIList<AbstractTree<GeomSeg*>*> atree_list;
  AbstractTree<GeomSeg*> *curr_tree;
  SegList *seg_list;
  GeomSeg  *curr_seg, *next_seg;
  int ii,jj, kk;
  double angle;
  min_angle = CUBIT_DBL_MAX;
  max_angle = 0.0;
  double creation_time=0.0;
  CpuTimer creation_timer;
  const double angle_convert = 180.0/CUBIT_PI;
  boundary_seg_loops.reset();
  int num_segs = 0;
  for (ii = 0; ii < boundary_seg_loops.size(); ii++ )
  {
    seg_list = boundary_seg_loops.get_and_step();
    creation_timer.cpu_secs();
    curr_tree = new RTree<GeomSeg*>(GEOMETRY_RESABS);
    creation_time += creation_timer.cpu_secs();
    for ( jj = seg_list->size(); jj > 0; jj-- )
    {
        //build the r-tree.
      num_segs++;
      creation_timer.cpu_secs();
      curr_seg = seg_list->get_and_step();
      curr_tree->add(curr_seg);
      creation_time += creation_timer.cpu_secs();
        //calculate the interior angles.
      next_seg = seg_list->get();

      RefEntity* temp_ent = curr_seg->get_end()->owner();
        //only measure the angle at vertices.  This
        //should make things slightly faster.
      if ( CAST_TO(temp_ent, RefVertex) )
      {
        stat = interior_angle(curr_seg, next_seg,
                              face, angle);
        if ( stat != CUBIT_SUCCESS )
        {
          min_angle = 0.0;
          max_angle = 360.0;
        }
        else if ( angle < min_angle )
        {
          min_angle = angle;
        }
        else if ( angle > max_angle )
        {
          max_angle = angle;
        }
      }
    }
    atree_list.append(curr_tree);
  }
  min_angle *= angle_convert;
  max_angle *= angle_convert;

  if ( atree_list.size() <= 1 )
    min_distance_between_loops = -1.0;
  else
  {
      //determine the minimum distance between the loops.
    CpuTimer atree_time;
    atree_time.cpu_secs();
    PointList *curr_points;
    GeomPoint *curr_point;
    DLIList <GeomSeg*> nearest_neighbors;
    
    CubitVector curr_loc;
    double closest_dist;
    min_distance_between_loops = CUBIT_DBL_MAX;
    atree_list.reset();
    for ( ii = 0; ii < atree_list.size(); ii++ )
    {
      curr_points = boundary_point_loops.get_and_step();
      for ( jj = ii+1; jj < atree_list.size(); jj++ )
      {
        curr_tree = atree_list.next(jj);
          //now for every point in curr_points, find it's closest_point
          //in the curr_tree.
        for ( kk = 0; kk < curr_points->size(); kk++ )
        {
          curr_point = curr_points->get_and_step();
          curr_loc.set(curr_point->coordinates());
          nearest_neighbors.clean_out();
          curr_tree->k_nearest_neighbor(curr_loc,
                                        1, closest_dist, nearest_neighbors,
                                        dist_sq_point_data);
          if ( nearest_neighbors.size() == 0)
          {
              //hmm, not sure what is wrong here.
            PRINT_ERROR("Can't find closest point between loops.\n");
            return CUBIT_FAILURE;
          }
          if (closest_dist < min_distance_between_loops)
            min_distance_between_loops = closest_dist;
        }
      }
    }
    double time =  atree_time.cpu_secs();
    PRINT_DEBUG_139("Time to create atree: %f\n", creation_time);
    PRINT_DEBUG_139("Time for atree nn search: %f\n", time);
  }
  if ( min_distance_between_loops != -1.0 )
  min_distance_between_loops = sqrt(min_distance_between_loops);
    //Now find the min distance inside a loop.  This is tricky
    //in that we want to find distances that are across some area
    //and not just between adjacent curves.
  atree_list.reset();
  GeomSeg *tmp_seg;
  double dist;
  min_distance_in_one_loop = CUBIT_DBL_MAX;
  int k=5;
  CubitVector curr_loc;
  DLIList <GeomSeg*> nearest_neighbors;
  RefEntity *ent_1, *ent_2;
  TopologyEntity *top_1, *top_2;
  
  
  for ( ii = atree_list.size(); ii > 0; ii-- )
  {
    curr_tree = atree_list.get_and_step();
    seg_list = boundary_seg_loops.get_and_step();
    for ( jj = seg_list->size(); jj > 0; jj-- )
    {
      nearest_neighbors.clean_out();
      curr_seg = seg_list->get_and_step();
        //get the k closest segments.
      curr_loc.set(curr_seg->get_start()->coordinates());
      curr_tree->k_nearest_neighbor(curr_loc,
                                    k, dist, nearest_neighbors,
                                    dist_sq_point_data);
        //go through these nearest_neighbors and through out
        //the segments that are attached to curr_seg.
      for ( kk = 0; kk < nearest_neighbors.size(); kk++ )
      {
        tmp_seg = nearest_neighbors.get();
        if ( tmp_seg == curr_seg ||
             tmp_seg == curr_seg->get_prev() ||
             tmp_seg == curr_seg->get_next() )
          nearest_neighbors.remove();
      }
      if ( nearest_neighbors.size() == 0 )
        continue;
        //The min distance is the distance between the curr_loc
        //and the top nearest_neighbor.
      nearest_neighbors.reset();
      tmp_seg = nearest_neighbors.get();
      ent_1 = tmp_seg->owner();
      ent_2 = curr_seg->get_start()->owner();
      top_1 = CAST_TO(ent_1, TopologyEntity);
      top_2 = CAST_TO(ent_2, TopologyEntity);
      
      if ( top_1->is_directly_related(top_2) )
        continue;
      
        //we'll need to recompute it...
      dist = dist_sq_point_data(curr_loc, tmp_seg);
      if ( dist < min_distance_in_one_loop )
        min_distance_in_one_loop = dist;
    }
  }
  if ( min_distance_in_one_loop == CUBIT_DBL_MAX )
  {
    ent_1 = (RefEntity*)face;
    DLIList<RefEntity*> entities;
    entities.append(ent_1);
    double smallest, largest, ave, sum;
    RefEdge *small_e, *large_e;
    measure_curve_length(entities, smallest, small_e,
                         largest, large_e, ave, sum);
    min_distance_in_one_loop = smallest;
  }
  else
    min_distance_in_one_loop = sqrt(min_distance_in_one_loop);
    //clean up the memory.
  int ll,mm;
  atree_list.reset();
  boundary_seg_loops.reset();
  for ( ll = boundary_seg_loops.size(); ll > 0; ll-- )
  {
    delete atree_list.pop();
    seg_list = boundary_seg_loops.pop();
    for ( mm = seg_list->size(); mm > 0; mm-- )
      delete seg_list->pop();
    delete seg_list;
  }
		
  PointList *point_list;
  for ( ll = boundary_point_loops.size(); ll > 0; ll-- )
  {
    point_list = boundary_point_loops.pop();
    for ( mm = point_list->size(); mm > 0; mm-- )
    {
      GeomPoint *tmp_point = point_list->pop();
        //delete the CubitPointData
      delete tmp_point;
    }
    delete point_list;
  }
  return CUBIT_SUCCESS;
}
//-------------------------------------------------------------------
//Function: dist_point_data
//Computes the distance between the line segment and the vector position.
//This function is used to determine the nearest neigbhor search and
//gets passed to the AbstractTree.
//-------------------------------------------------------------------
double GeomMeasureTool::dist_sq_point_data( CubitVector &curr_point,
                                            GeomSeg *&curr_seg )
{
  double node[3];
  double point_0[3], point_1[3];
  node[0] = curr_point.x();
  node[1] = curr_point.y();
  node[2] = curr_point.z();
  CubitVector start_v = curr_seg->get_start()->coordinates();
  CubitVector end_v = curr_seg->get_end()->coordinates();
  point_0[0] = start_v.x();
  point_0[1] = start_v.y();
  point_0[2] = start_v.z();
  point_1[0] = end_v.x();
  point_1[1] = end_v.y();
  point_1[2] = end_v.z();

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
      dist = (curr_point -start_v).length_squared();
    }
    else
    {
      dist = (curr_point - end_v).length_squared();
    }
  }
  else
  {
      //I hate to do this but everything else is in
      //terms of distance squared.
    dist = dist*dist;
  }
  return dist;
}

CubitStatus GeomMeasureTool::convert_to_lines( PointLoopList &boundary_point_loops,
                                               SegLoopList &boundary_line_loops,
                                               RefFace *ref_face )
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
  GeomSeg *new_line, *prev = NULL, *next = NULL;
  SegList allocated_segs;
  RefEntity *entity_start, *entity_end, *seg_owner;
  boundary_point_loops.reset();
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
            PRINT_ERROR("GeomMeasureTool::convert_to_lines having problem finding owner \n" 
                        "       of boundary.  Could be indicative of corrupt geometry.\n");
            int ll;
            for ( ll = boundary_line_loops.size(); ll > 0; ll-- )
              delete boundary_line_loops.pop();
            for ( ll = allocated_segs.size(); ll > 0; ll-- )
              delete allocated_segs.pop();
            delete segment_list;
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
              int ll;
              for ( ll = boundary_line_loops.size(); ll > 0; ll-- )
                delete boundary_line_loops.pop();
              for ( ll = allocated_segs.size(); ll > 0; ll-- )
                delete allocated_segs.pop();
              delete segment_list;
              return CUBIT_FAILURE;
            }
            seg_owner = closest_edge;
          }
        }
      }
      new_line = new GeomSeg( point_list->get(), point_list->next(), seg_owner);
      allocated_segs.append(new_line);
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
    }
    else
      delete segment_list;
  }
  return CUBIT_SUCCESS;
}

CubitStatus GeomMeasureTool::get_boundary_points( RefFace *ref_face,
                                                  PointLoopList &boundary_point_loops,
                                                  double seg_length_tol)
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
  
  for ( ii = co_edge_loops.size(); ii > 0; ii-- )
  {
    co_edge_list_ptr = co_edge_loops.get_and_step();
    new_point_loop_ptr = new PointList;
    for ( jj = co_edge_list_ptr->size(); jj > 0; jj-- )
    {
      co_edge_ptr = co_edge_list_ptr->get_and_step();
      ref_edge_ptr = co_edge_ptr->get_ref_edge_ptr();
      tmp_point_list.clean_out();
      stat = get_curve_facets( ref_edge_ptr, tmp_point_list,
                               seg_length_tol );
      if ( stat != CUBIT_SUCCESS )
        return CUBIT_FAILURE;
      if ( tmp_point_list.size() == 0 )
        continue;
      tmp_point_list.reset();
        //the points are in order from start vertex to end vertex.
        //append them now according to the loop.
      sense = co_edge_ptr->get_sense();
      if ( CUBIT_FORWARD != sense )
        tmp_point_list.reverse();
        //Now take off the last point as it is a duplicate with the
        //other list...
      tmp_point_list.reset();
      delete tmp_point_list.pop();
      (*new_point_loop_ptr) += tmp_point_list;
    }
    boundary_point_loops.append(new_point_loop_ptr);
  }
    //clean up the list memory.
  for(ii = co_edge_loops.size(); ii>0; ii-- )
    delete co_edge_loops.pop();
  co_edge_loops.clean_out();
      
  return CUBIT_SUCCESS;
}

CubitStatus GeomMeasureTool::get_curve_facets( RefEdge* curve,
                                               PointList &segments,
                                               double seg_length_tol ) 
{
//  const double COS_ANGLE_TOL =  0.965925826289068312213715; // cos(15)
//  const double COS_ANGLE_TOL =  0.984807753012208020315654; // cos(10)
  //const double COS_ANGLE_TOL =  0.996194698091745545198705; // cos(5)
  GMem curve_graphics;
    //make this tol bigger than seg_length_tol, to
    //make sure the segments are larger than the tolerance.
  const double dist_tol = seg_length_tol;// + .05*seg_length_tol;
  //const double dist_tol_sqr = dist_tol*dist_tol;
  if ( curve->get_curve_ptr() == NULL )
  {
    PRINT_ERROR("Curve %d had no geometry, this is a bug!\n", curve->id());
    return CUBIT_FAILURE;
  }
  else if ( curve->get_curve_ptr()->get_geometry_query_engine() == NULL )
  {
    PRINT_ERROR("Curve %d had no geometry engine, this is a bug!\n",
                curve->id());
    return CUBIT_FAILURE;
  }
  curve->get_graphics( curve_graphics );
  
  GPoint* gp = curve_graphics.point_list();
  if ( gp == NULL )
  {
    if ( curve->measure() < GEOMETRY_RESABS ) 
    {
      if( curve->geometry_type() != POINT_CURVE_TYPE )
      {
        PRINT_INFO("Curve %d has a zero length, and won't be used for loop measurements.\n",
                    curve->id());
      }
    }
    else
      PRINT_ERROR("Curve %d had not faceting.  Ignoring this curve.\n",
                  curve->id());
    return CUBIT_SUCCESS;
  }
  GeomPoint *last = new GeomPoint( gp[0].x, gp[0].y, gp[0].z, curve );
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
    last = new GeomPoint(pos,curve);
    segments.append( last );
    lastv = last->coordinates();
  }
    // Now check if the segment list is reversed wrt the curve direction.
  segments.reset();
  if ( remove_second_to_end )
  {
    if ( segments.size() == 2 )
    {
        //just leave it in, let it be small...
    }
    else
    {
        //Remove the second to last one.  To do
        //this efficiently (don't do remove), pop
        //the last one, then save that and
        //re-add it after poping the second one.
      GeomPoint *temp = segments.pop();
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
  }
  if ( segments.size() < 2 )
  {
    PRINT_ERROR("Problem faceting boundary.\n");
    return CUBIT_FAILURE;
  }
  segments.get()->owner(curve->start_vertex());
  segments.prev()->owner(curve->end_vertex());
  return CUBIT_SUCCESS;
}
CubitStatus GeomMeasureTool::interior_angle(GeomSeg *first_seg,
                                            GeomSeg *next_seg,
                                            RefFace *face,
                                            double &angle )
{
  if ( first_seg->get_end() != next_seg->get_start() )
  {
    PRINT_ERROR("GeomMeasureTool::interior_angle expects segments\n"
                "to be connected (Users: this is a bug, please submit.\n");
    return CUBIT_FAILURE;
  }
  CubitVector point = first_seg->get_end()->coordinates();
  CubitVector to_prev = first_seg->get_start()->coordinates();
  to_prev -= point;
  CubitVector to_next = next_seg->get_end()->coordinates();
  to_next -= point;
  CubitVector surf_norm = face->normal_at(point);
  if ( surf_norm.length_squared() < CUBIT_RESABS )
  {
      //Try getting it at one of the other nodes...
    surf_norm = face->normal_at(first_seg->get_start()->coordinates() );
    if (surf_norm.length_squared() < CUBIT_RESABS )
    {
        //Try getting it at one of the other nodes...
      surf_norm = face->normal_at(next_seg->get_end()->coordinates() );
      if (surf_norm.length_squared() < CUBIT_RESABS )
      {
        PRINT_ERROR("Trying to get normal at surf %d.\n"
                    "       Normal length being returned equals zero.\n",
                    face->id() );
        angle = CUBIT_DBL_MAX;
        return CUBIT_FAILURE;
      }
    }
  }
  angle = surf_norm.vector_angle ( to_next, to_prev );
  return CUBIT_SUCCESS;
}

void GeomMeasureTool::measure_face_area (DLIList <RefEntity*> &entity_list,
                                         double &smallest,
                                         RefFace *&smallest_face,
                                         double &largest,
                                         RefFace *&largest_face,
                                         double &average,
                                         double &sum)
{
  DLIList<RefFace*> ref_faces;
  get_faces_from_list(entity_list, ref_faces);
  measure_face_area(ref_faces, smallest, smallest_face,
                    largest, largest_face, average, sum);
  
}

void GeomMeasureTool::measure_face_area (DLIList <RefFace*> &ref_faces,
                                         double &smallest,
                                         RefFace *&smallest_face,
                                         double &largest,
                                         RefFace *&largest_face,
                                         double &average,
                                         double &sum)
{
  int ii;
  double curr_area, total_area = 0.0;
  smallest = CUBIT_DBL_MAX;
  largest = 0.0;
  RefFace *curr_face;
  for( ii = ref_faces.size(); ii > 0; ii-- )
  {
    curr_face = ref_faces.get_and_step();
    curr_area = measure_area(curr_face);
    if ( curr_area < smallest )
    {
      smallest = curr_area;
      smallest_face = curr_face;
    }
    if( curr_area > largest )
    {
      largest = curr_area;
      largest_face = curr_face;
    }
    total_area += curr_area;
  }
  if (ref_faces.size() != 0 )
    average = total_area/ref_faces.size();
  sum = total_area;
  
}

void GeomMeasureTool::measure_volume_volume (DLIList <RefEntity*> &entity_list,
                                             double &smallest,
                                             RefVolume *&smallest_volume,
                                             double &largest,
                                             RefVolume *&largest_volume,
                                             double &average,
                                             double &sum)
{
  DLIList<RefVolume*> ref_volumes;
  get_volumes_from_list(entity_list, ref_volumes);
  measure_volume_volume(ref_volumes, smallest, smallest_volume,
                        largest, largest_volume, average, sum);
  
}

void GeomMeasureTool::measure_volume_volume (DLIList <RefVolume*> &ref_volumes,
                                             double &smallest,
                                             RefVolume *&smallest_volume,
                                             double &largest,
                                             RefVolume *&largest_volume,
                                             double &average,
                                             double &sum)
{
  int ii;
  double curr_volume, total_volume = 0.0;
  smallest = CUBIT_DBL_MAX;
  largest = 0.0;
  RefVolume *curr_solid;
  for ( ii = ref_volumes.size(); ii > 0; ii-- )
  {
    curr_solid = ref_volumes.get_and_step();
    curr_volume = curr_solid->measure();

    if ( curr_volume < smallest )
    {
      smallest = curr_volume;
      smallest_volume = curr_solid;
    }
    if ( curr_volume > largest )
    {
      largest = curr_volume;
      largest_volume = curr_solid;
    }
    total_volume += curr_volume;
  }
  if (ref_volumes.size() != 0 )
    average = total_volume/ref_volumes.size();
  sum = total_volume;
  
}

void GeomMeasureTool::measure_face_area_and_hydraulic_radius(RefFace *curr_face,
                                                             double &face_area,
                                                             double &face_hydraulic_radius)
{
  int ii;
  double total_edge_length = 0.0, curr_edge_length;
  RefEdge *curr_edge;
  DLIList <CoEdge*> my_co_edges;

  face_area = measure_area(curr_face);
  curr_face->co_edges(my_co_edges);
  for( ii = my_co_edges.size(); ii > 0; ii-- )
  {
    curr_edge = my_co_edges.get_and_step()->get_ref_edge_ptr();
    curr_edge_length = curr_edge->measure();
    total_edge_length += curr_edge_length;
  }
  if (total_edge_length != 0)
    face_hydraulic_radius = 4.0*(face_area / total_edge_length);
}

void GeomMeasureTool::measure_volume_volume_and_hydraulic_radius(RefVolume *curr_volume,
                                                                 double &volume_volume,
                                                                 double &volume_hydraulic_radius)
{
  int ii;
  double curr_face_area, total_surface_area = 0.0;
  DLIList <RefFace*> ref_face_list;
  RefFace *curr_face;

  volume_volume = curr_volume->measure();
  curr_volume->ref_faces(ref_face_list);
  
  for ( ii = ref_face_list.size(); ii > 0; ii-- )
  {
    curr_face = ref_face_list.get_and_step();
    curr_face_area = measure_area(curr_face);
    total_surface_area += curr_face_area;
  }
  if (total_surface_area != 0)
    volume_hydraulic_radius = 6.0*(volume_volume / total_surface_area);
 
}

void GeomMeasureTool::angles_between_volume_surfaces(RefVolume *curr_volume,
                                                     double &min_angle,
                                                     double &max_angle,
                                                     RefFace *&face_min_1,
                                                     RefFace *&face_min_2,
                                                     RefFace *&face_max_1,
                                                     RefFace *&face_max_2)
{
  int ii, jj, shared_angles = 0;
  double common_angle, sum_of_angles = 0.0;
  double smallest = CUBIT_DBL_MAX, largest = 0.0;
  DLIList <RefFace*> face_list_1, face_list_2;
  RefFace *curr_face_1, *curr_face_2;
#ifdef BOYD17
  DLIList <RefEdge*> my_ref_edges;
#endif
  RefEdge *common_edge;

  curr_volume->ref_faces(face_list_1);
  face_list_2 = face_list_1;

  for (ii = face_list_1.size(); ii > 0; ii--)
  {
    curr_face_1 = face_list_1.get_and_step();
    for(jj = face_list_2.size(); jj > 0; jj--)
    {
      curr_face_2 = face_list_2.get_and_step();

      if(curr_face_1 == curr_face_2)
        continue;
      else
      {
        common_edge = curr_face_1->common_ref_edge(curr_face_2);
        if( common_edge == NULL )
          continue;
        else
        {
          shared_angles++;
          
          if ( common_edge->get_curve_ptr()->geometry_type() != STRAIGHT_CURVE_TYPE )
          {
            int kk;
            double frac_pos = 0.00;
            for ( kk = 0; kk < 4; kk++ )
            {
              common_angle =  GeometryQueryTool::instance()->surface_angle(curr_face_1, curr_face_2,
                                                                           common_edge, curr_volume,
                                                                           frac_pos);
              frac_pos += .25;
              if( common_angle < smallest )
              {
                min_angle = common_angle * (180 / CUBIT_PI);
                smallest = common_angle;
                face_min_1 = curr_face_1;
                face_min_2 = curr_face_2;
              }
              if( common_angle > largest )
              {
                max_angle = common_angle * (180 / CUBIT_PI);
                largest = common_angle;
                face_max_1 = curr_face_1;
                face_max_2 = curr_face_2;
              }
            }
          }
          else
          {
            
            common_angle =  GeometryQueryTool::instance()->surface_angle(curr_face_1, curr_face_2,
                                                                         common_edge, curr_volume );
            if( common_angle < smallest )
            {
              min_angle = common_angle * (180 / CUBIT_PI);
              smallest = common_angle;
              face_min_1 = curr_face_1;
              face_min_2 = curr_face_2;
            }
            if( common_angle > largest )
            {
              max_angle = common_angle * (180 / CUBIT_PI);
              largest = common_angle;
              face_max_1 = curr_face_1;
              face_max_2 = curr_face_2;
            }
          }
          
          sum_of_angles += common_angle;
        }
      }
    }
  }
}

void GeomMeasureTool::merged_unmerged_surface_ratio(DLIList <RefVolume*> &ref_volumes,
                                                    int &merged, int &unmerged,
                                                    double &ratio)
  
{
  int ii, jj;
  RefVolume *curr_volume_1;
  DLIList <RefFace*> surface_list_1, surface_list_2;
  RefFace *curr_face_1;
  int merged_surfaces = 0, unmerged_surfaces = 0;
  
  for( ii = ref_volumes.size(); ii > 0; ii-- )
  {
    curr_volume_1 = ref_volumes.get_and_step();
    surface_list_1.clean_out();
    
    curr_volume_1->ref_faces(surface_list_1);
    for( jj = surface_list_1.size(); jj > 0; jj-- )
    {
      curr_face_1 = surface_list_1.get_and_step();
      if ( curr_face_1->marked() )
        continue;
      else
      {
        curr_face_1->marked(1);
        surface_list_2.append(curr_face_1);
      }
      if ( curr_face_1->num_ref_volumes() == 2 )
        merged_surfaces++;
      else
        unmerged_surfaces++;
    }
  }
  for ( ii = surface_list_2.size(); ii > 0; ii-- )
  {
    surface_list_2.get_and_step()->marked(0);
  }
  merged = merged_surfaces;
  unmerged = unmerged_surfaces;
  if ( unmerged_surfaces == 0 )
  {
    ratio = 100.0;
  }
  else
    ratio = (double) ((double)merged_surfaces / (double)unmerged_surfaces);
}

void GeomMeasureTool::report_intersected_volumes(DLIList <RefVolume*> &ref_vols,
                                                 DLIList <RefVolume*> &intersection_list)
{
  DLIList <RefVolume*> results;
  RefVolume *curr_vol;
  int i, j;
  ProgressTool *progress_ptr = NULL;
  int total_volumes = ref_vols.size();
  if (total_volumes > 5)
  {
    progress_ptr = AppUtil::instance()->progress_tool();
    assert(progress_ptr != NULL);
    progress_ptr->start(0, 100, "Overlapping Volumes Progress", 
      NULL, CUBIT_TRUE, CUBIT_TRUE);
  }
  double curr_percent = 0.0;

  for( i = 0; i < ref_vols.size(); i++ )
  {
    curr_vol = ref_vols.next(i);

    //is this body a multi-volume-body?
    Body *single_volume_body = curr_vol->body();
    DLIList<RefVolume*> body_vols;
    single_volume_body->ref_volumes( body_vols );
    if( body_vols.size() > 1 )
      single_volume_body = NULL;
      
    for( j = (i + 1); j < ref_vols.size(); j++ )
    {
      RefVolume *curr_vol2 = ref_vols.next(j);
      if ( CubitMessage::instance()->Interrupt() )
      {
          //interrpt.  We need to exit.
        if ( progress_ptr != NULL )
        {
          progress_ptr->end();
        }
          //just leave what has been calculated...
        return;
      }

      Body *single_volume_body2 = curr_vol2->body();
      DLIList<RefVolume*> body_vols2;
      single_volume_body2->ref_volumes( body_vols2 );
      if( body_vols2.size() > 1 )
        single_volume_body2 = NULL;

      //update the progress..
      if ( progress_ptr != NULL )
      {
        curr_percent = ((double)(i+1))/((double)(total_volumes));
        progress_ptr->percent(curr_percent);
      }
    
      //if both are single-volume-bodies
      if( single_volume_body && single_volume_body2 )
      {
        if( GeometryQueryTool::instance()->bodies_overlap( single_volume_body,
                                                       single_volume_body2 ) )
        {
          intersection_list.append( curr_vol );
          intersection_list.append( curr_vol2 );
        }
      }
      else if( GeometryQueryTool::instance()->volumes_overlap( curr_vol, curr_vol2 ) )
      {
       intersection_list.append( curr_vol );
       intersection_list.append( curr_vol2 );
      }
    }
  }
  
  if ( progress_ptr != NULL )
  {
    progress_ptr->end();
  }
}

void GeomMeasureTool::report_intersected_bodies(DLIList <Body*> &ref_bodies,
                                                DLIList <Body*> &intersection_list)
{
  Body *curr_body, *curr_body_2;
  int ii, jj;
  ProgressTool *progress_ptr = NULL;
  int total_bodies = ref_bodies.size();
  if (total_bodies > 5)
  {
    progress_ptr = AppUtil::instance()->progress_tool();
    assert(progress_ptr != NULL);
    progress_ptr->start(0, 100, "Overlapping Volumes Progress", 
      NULL, CUBIT_TRUE, CUBIT_TRUE);
  }
  double curr_percent = 0.0;
 
  
  for( ii = 0; ii < ref_bodies.size(); ii++ )
  {
    curr_body = ref_bodies.next(ii);

    for( jj = (ii + 1); jj < ref_bodies.size(); jj++ )
    {
      curr_body_2 = ref_bodies.next(jj);
      if ( CubitMessage::instance()->Interrupt() )
      {
          //interrpt.  We need to exit.
        if ( progress_ptr != NULL )
        {
          progress_ptr->end();
        }
          //just leave what has been calculated...
        return;
      }

      if (GeometryQueryTool::instance()->bodies_overlap(curr_body,
                                                        curr_body_2) )
      {
        intersection_list.append(curr_body);
        intersection_list.append(curr_body_2);
      }
    }
       //update the progress..
    if ( progress_ptr != NULL )
    {
      curr_percent = ((double)(ii+1))/((double)(total_bodies));
      progress_ptr->percent(curr_percent);
    }
  }
  if ( progress_ptr != NULL )
  {
    progress_ptr->end();
  }
}

void GeomMeasureTool::face_list_from_volume_list( DLIList <RefVolume*> &input_vols,
                                                  DLIList <RefFace*> &all_faces )
{
  DLIList <RefFace*> curr_faces;
  RefVolume *curr_vol;
  RefFace *curr_face;
  int ii, jj;

    //get the all RefFaces from input_vols.
    //put them into all_faces list.
  
  for ( ii = 0; ii < input_vols.size(); ii++ )
  {
    curr_vol = input_vols.get_and_step();
    curr_faces.clean_out();
    curr_vol->ref_faces(curr_faces);
    
    for( jj = 0; jj < curr_faces.size(); jj++ )
    {
      curr_face = curr_faces.get_and_step();
      if( curr_face->marked() )
        continue;
      else
      {
        curr_face->marked(1);
        all_faces.append(curr_face);
      }
    }
  }
  for ( ii = all_faces.size();ii >0; ii-- )
    all_faces.get_and_step()->marked(0);
}

RefFace* GeomMeasureTool::valid_start( DLIList <RefFace*> &all_faces )
{
  int ii;
  RefFace *start_face;

  for( ii = 0; ii < all_faces.size(); ii++ )
  {
    start_face = all_faces.get_and_step();

    if( !start_face->marked() && start_face->num_ref_volumes()== 1 )
      return start_face;
  }

  return NULL;
  
}

void GeomMeasureTool::find_shells( DLIList <RefVolume*> &input_vols,
                                   RefGroup *&owner_groups,
                                   int &number_of_shells)
{
  DLIList <DLIList<RefEntity*>*> list_of_lists;
  DLIList <RefEntity*> *new_list;
  DLIList <RefFace*> all_faces, stack, related_faces;
  DLIList <RefEdge*> curr_edges;
  RefFace *start_face, *curr_face, *related_face;
  RefEdge *curr_edge;
  int ii, jj;
  RefEntity *tmp_ent;
  number_of_shells = 0;
  DLIList <RefFace*> marked_faces;
  
  face_list_from_volume_list( input_vols, all_faces );
  
    // all_faces all marked - need unmarked
  for (ii = 0; ii < all_faces.size(); ii++ )
  {
    curr_face = all_faces.get_and_step();
    curr_face->marked(0);
  }
  while ( (start_face = valid_start(all_faces) ) != NULL )
  {
      //create a list.
    stack.clean_out();
    new_list = new DLIList<RefEntity*>;
    stack.append(start_face);
    start_face->marked(1);
    marked_faces.append(start_face);

    while (stack.size() > 0 )
    {
      curr_edges.clean_out();

      curr_face = stack.pop(); // get last face added
      tmp_ent = CAST_TO(curr_face, RefEntity);
      new_list->append(tmp_ent); // append that face to the newest list
      curr_face->ref_edges(curr_edges); // gets list of face's edges

      for( ii = 0; ii < curr_edges.size(); ii++ )  // loops edges
      {
        related_faces.clean_out();

        curr_edge = curr_edges.get_and_step(); // chooses next edge
        curr_edge->ref_faces(related_faces); // finds edge's related faces
        for( jj = 0; jj < related_faces.size(); jj++ ) // loops faces
        {
          related_face = related_faces.get_and_step(); // gets face
          if(!related_face->marked() && related_face->num_ref_volumes() == 1)
          {
            stack.append(related_face);
            related_face->marked(1);
            marked_faces.append(related_face);
          }
        }
      }
    }
    list_of_lists.append(new_list);
    number_of_shells++;
  }

  for (ii = 0; ii < marked_faces.size(); ii++ )
  {
    curr_face = marked_faces.get_and_step();
    curr_face->marked(0);
  }
  
  RefGroup *curr_group;
  owner_groups = RefEntityFactory::instance()->construct_RefGroup("volume_shells");
  for( ii = list_of_lists.size(); ii > 0; ii-- )
  {
    new_list = list_of_lists.pop();
    curr_group = RefEntityFactory::instance()->construct_RefGroup(*new_list);
    tmp_ent = CAST_TO(curr_group, RefEntity);
    owner_groups->add_ref_entity(tmp_ent);
    delete new_list;
  } 
}

void GeomMeasureTool::print_surface_measure_summary( DLIList <RefFace*> &ref_faces )
                                                  
{
  
  double min_dist_loops=CUBIT_DBL_MAX, min_dist_loop=CUBIT_DBL_MAX;
  double curr_min_dist_loops=CUBIT_DBL_MAX, curr_min_dist_loop=CUBIT_DBL_MAX;
  double min_angle=CUBIT_DBL_MAX, max_angle=-CUBIT_DBL_MAX;
  double curr_min_angle=CUBIT_DBL_MAX, curr_max_angle=-CUBIT_DBL_MAX;
  double curr_face_area, curr_face_hydraulic_radius;
  double min_face_area=CUBIT_DBL_MAX, min_face_hydraulic_radius=CUBIT_DBL_MAX;
  double max_face_area=-CUBIT_DBL_MAX;
  int counter = 0;
  double total_area= 0.0;
  double min_curve_length = CUBIT_DBL_MAX;
  double max_curve_ratio = -CUBIT_DBL_MAX;
  double curve_length, curve_ratio;
  int ii;
       
  RefFace *curr_face;
  RefFace *face_min_dist_loops=NULL, *face_min_dist_loop=NULL;
  RefFace *face_min_angle=NULL, *face_max_angle=NULL;
  RefFace *face_min_area=NULL, *face_min_hydraulic_radius=NULL;
  RefFace *face_max_area=NULL;
  RefFace *face_min_curve_length=NULL, *face_max_curve_ratio=NULL;
  RefEdge *smallest_edge=NULL, *curr_small_edge=NULL;
       
       
  for( ii = 0; ii < ref_faces.size(); ii++ )
  {
    curr_face = ref_faces.get_and_step();
    if ( curr_face->num_ref_edges() > 0 )
    {
      GeomMeasureTool::measure_face_loops(curr_face,
                                          curr_min_dist_loops,
                                          curr_min_dist_loop,
                                          curr_min_angle,
                                          curr_max_angle,
                                          GEOMETRY_RESABS*500);
      GeomMeasureTool::measure_face_area_and_hydraulic_radius(curr_face,
                                                              curr_face_area,
                                                              curr_face_hydraulic_radius);
      GeomMeasureTool::measure_face_curves(curr_face, curve_length,
                                           curve_ratio, curr_small_edge);
    }
    else
      curr_face_area = curr_face->measure();
         
    total_area += curr_face_area;
    counter++;
    if ( curr_face->num_ref_edges() > 0 )
    {
      if ( curve_length < min_curve_length )
      {
        min_curve_length = curve_length;
        smallest_edge = curr_small_edge;
        face_min_curve_length = curr_face;
      }
      if ( curve_ratio != -1.0 && curve_ratio > max_curve_ratio )
      {
        max_curve_ratio = curve_ratio;
        face_max_curve_ratio = curr_face;
      }
      if ( curr_min_dist_loops != -1.0 && curr_min_dist_loops < min_dist_loops )
      {
        face_min_dist_loops = curr_face;
        min_dist_loops = curr_min_dist_loops;
      }
      if ( curr_min_dist_loop < min_dist_loop )
      {
        face_min_dist_loop = curr_face;
        min_dist_loop = curr_min_dist_loop;
      }
      if ( curr_min_angle < min_angle )
      {
        face_min_angle = curr_face;
        min_angle = curr_min_angle;
      }
      if ( curr_max_angle > max_angle )
      {
        face_max_angle = curr_face;
        max_angle = curr_max_angle;
      }

      if ( curr_face_hydraulic_radius < min_face_hydraulic_radius )
      {
        min_face_hydraulic_radius = curr_face_hydraulic_radius;
        face_min_hydraulic_radius = curr_face;
      }

    }
    if ( curr_face_area < min_face_area )
    {
      min_face_area = curr_face_area;
      face_min_area = curr_face;
    } 
    if ( curr_face_area > max_face_area )
    {
      max_face_area = curr_face_area;
      face_max_area = curr_face;
    }
  }
  if ( counter > 0 )
  {
    PRINT_INFO("\n******Summary of Surface Measure Info*******\n\n");
    PRINT_INFO("Measured %d surfaces: Average Area: %f \n\n", counter, total_area/counter );
    PRINT_INFO("Minimum Surface Area is: %f on surface %d\n\n",
               min_face_area, face_min_area->id());
    PRINT_INFO("Maximum Surface Area is: %f on surface %d\n\n",
               max_face_area, face_max_area->id());
    if ( face_min_hydraulic_radius )
      PRINT_INFO("Minimum Hydraulic Radius (Area/Perimeter) is: %f on surface %d\n\n",
                 min_face_hydraulic_radius, face_min_hydraulic_radius->id());
    if ( smallest_edge )
      PRINT_INFO("Minimum Curve Length is: %f for curve %d on surface %d\n\n",
                 min_curve_length, smallest_edge->id(), face_min_curve_length->id());
    if ( face_max_curve_ratio )
      PRINT_INFO("Maximum Ratio of Adjacent Curve Lengths is: %f on surface %d\n\n",
                 max_curve_ratio, face_max_curve_ratio->id());
    if ( face_min_dist_loops && face_min_dist_loops->num_loops() > 1 )
    {
      PRINT_INFO("Minimum Distance Between Loops is: %f on surface %d\n\n",
                 min_dist_loops, face_min_dist_loops->id());
    }
    if ( face_min_dist_loop )
      PRINT_INFO("Minimum Distance Between a Single Loop is: %f on surface %d\n\n",
                 min_dist_loop, face_min_dist_loop->id());
    if ( face_min_angle )
      PRINT_INFO("Minimum Interior Angle Between Curves is: %f Degrees on surface %d\n\n",
                 min_angle, face_min_angle->id());
    if ( face_max_angle )
      PRINT_INFO("Maximum Interior Angle Between Curves is: %f Degrees on surface %d\n\n",
                 max_angle, face_max_angle->id());
	center(ref_faces);
    PRINT_INFO("************End of Summary***********\n\n");
  }

}
void GeomMeasureTool::find_adjacent_face_ratios(RefVolume *curr_volume, double &max_face_ratio,
                                                RefFace *&big_face,
                                                RefFace *&small_face)
{
    //Loop over all the faces in the model.  Find the maximum face area ratio.
    //(larger face area / smaller face area)
  
  DLIList <Shell*> shells;
  curr_volume->shells(shells);
  DLIList <RefFace*> ref_face_stack, adj_faces;
  Shell* curr_shell;
  RefFace *curr_ref_face, *adj_face;
    //AreaHashTuple is defined in the .hpp file.
#ifdef BOYD17
  DLIList <AreaHashTuple*> *hash_table, hash_index;
#endif
  DLIList <AreaHashTuple*> *hash_table;
  max_face_ratio = -CUBIT_DBL_MAX;
  
  
  int ii, jj, table_size;
  int hash_val;
  double tmp_area;
  double tmp_ratio;
  
  AreaHashTuple* curr_tuple, *adj_tuple;
  
  for ( ii = shells.size(); ii > 0; ii-- )
  {
    curr_shell = shells.get_and_step();
    ref_face_stack.clean_out();
    curr_shell->ref_faces(ref_face_stack);
    table_size = ref_face_stack.size();
      //build a hash table with the most simple form
      //of hashing based on the ref face id.
      //To get this more reliable, make table size a prime
      //number...
      //Do this so that we don't have to calculate the
      //areas more than once...
    hash_table = new DLIList<AreaHashTuple*> [table_size];
    for ( jj = 0; jj < ref_face_stack.size(); jj++ )
    {
      curr_ref_face = ref_face_stack.get_and_step();
      tmp_area = curr_ref_face->measure();
      curr_tuple = new AreaHashTuple();
      curr_tuple->myFace = curr_ref_face;
      curr_tuple->myArea = tmp_area;
      hash_val = curr_ref_face->id() % table_size;
      hash_table[hash_val].append(curr_tuple);
    }
    
    while ( ref_face_stack.size() )
    {
      curr_ref_face = ref_face_stack.pop();
        //Now get the adjacent faces.
      get_adjacent_faces(curr_ref_face, adj_faces);
      find_index(hash_table, table_size, curr_ref_face, curr_tuple);
      assert(curr_tuple != NULL);
      if ( curr_tuple == NULL )
        return;
      CubitBoolean curr_is_big = CUBIT_TRUE;
      for ( jj = 0; jj < adj_faces.size(); jj++ )
      {
        adj_face = adj_faces.get_and_step();
        find_index(hash_table, table_size, adj_face, adj_tuple);
        if ( adj_tuple == NULL )
          continue;
        curr_is_big = CUBIT_TRUE;
        if ( curr_tuple->myArea < adj_tuple->myArea )
        {
          tmp_ratio = adj_tuple->myArea / curr_tuple->myArea;
          curr_is_big = CUBIT_FALSE;
        }
        else
          tmp_ratio = curr_tuple->myArea / adj_tuple->myArea;
        if ( tmp_ratio > max_face_ratio )
        {
          max_face_ratio = tmp_ratio;
          if ( curr_is_big )
          {
            big_face = curr_ref_face;
            small_face = adj_face;
          }
          else
          {
            big_face = adj_face;
            small_face = curr_ref_face;
          }
        }
      }
    }
    for ( jj = 0; jj < table_size; jj++ )
    {
        //delete the area tuples.
      while ( hash_table[jj].size() > 0 )
        delete hash_table[jj].pop();
    }
      //delete the hash_table.
    delete [] hash_table;
  }
}
void GeomMeasureTool::find_index(DLIList<AreaHashTuple*> *hash_table,
                                 int table_size,
                                 RefFace *ref_face,
                                 AreaHashTuple *&curr_tuple )
{
  int ii;
  int hash_val = ref_face->id() % table_size;
  if ( hash_table[hash_val].size() == 0 )
  {
      //this ref face isn't store here.
    curr_tuple = NULL;
    return;
  }
  int curr_size = hash_table[hash_val].size();
    //resolve collisions by separate chaining method...
  for ( ii = 0;ii < curr_size; ii++ )
  {
    if ( hash_table[hash_val].next(ii)->myFace == ref_face )
    {
      curr_tuple = hash_table[hash_val].next(ii);
      return;
    }
  }
  curr_tuple = NULL;
  return;
}
  
void GeomMeasureTool::get_adjacent_faces(RefFace* curr_face,
                                         DLIList<RefFace*> &adjacent_faces)
{
  DLIList <RefEdge*> ref_edges;
  DLIList <RefFace*> tmp_ref_faces;
  curr_face->ref_edges(ref_edges);
  RefEdge *ref_edge;
  RefFace *tmp_ref_face;
  int ii, jj;
  for ( ii = 0; ii < ref_edges.size(); ii++ )
  {
    ref_edge = ref_edges.get_and_step();
    tmp_ref_faces.clean_out();
    ref_edge->ref_faces(tmp_ref_faces);
    for ( jj = 0; jj < tmp_ref_faces.size(); jj++ )
    {
      tmp_ref_face = tmp_ref_faces.get_and_step();
      if ( !tmp_ref_face->marked() && tmp_ref_face != curr_face )
      {
        tmp_ref_face->marked(1);
        adjacent_faces.append(tmp_ref_face);
      }
      else
        continue;
    }
  }
  for ( ii = adjacent_faces.size(); ii > 0; ii-- )
    adjacent_faces.get_and_step()->marked(0);
}
void GeomMeasureTool::print_volume_measure_summary(DLIList <RefVolume*> &ref_volumes)
{
  int ii;
  RefVolume *curr_volume;
  RefFace *curr_face;
  DLIList <RefFace*> ref_faces, tmp_faces;
    //First get a list of all the ref_faces.
    //Get them uniquely.
  for ( ii = 0; ii< ref_volumes.size(); ii++ )
  {
    curr_volume = ref_volumes.get_and_step();
    tmp_faces.clean_out();
    curr_volume->ref_faces(tmp_faces);
    int jj;
    for ( jj = 0; jj < tmp_faces.size(); jj++ )
    {
      curr_face = tmp_faces.get_and_step();
      if ( curr_face->marked() )
        continue;
      else
      {
        curr_face->marked(1);
        ref_faces.append(curr_face);
      }
    }
  }
    //reset the mark.
  for ( ii = 0; ii< ref_faces.size(); ii++ )
    ref_faces.get_and_step()->marked(0);

  double volume_volume, volume_hydraulic_radius;
  double min_angle=CUBIT_DBL_MAX, max_angle=-CUBIT_DBL_MAX;
  double total_volume = 0.0, min_volume_hydraulic_radius = CUBIT_DBL_MAX;
  double min_volume=CUBIT_DBL_MAX; //, max_volume=-CUBIT_DBL_MAX;
  double curr_min_angle = 361.0, curr_max_angle =361.;
  double max_face_ratio = -CUBIT_DBL_MAX, curr_face_ratio=CUBIT_DBL_MAX;
         

  RefVolume *volume_min_volume = NULL, *volume_min_hydraulic_radius=NULL;
  RefVolume *volume_min_angle=NULL, *volume_max_angle=NULL;
  RefVolume *volume_face_ratio= NULL;
  RefFace *max_small_face=NULL, *max_big_face=NULL;
  RefFace *small_face=NULL, *big_face=NULL;
  RefFace *face_min_1=NULL, *face_min_2=NULL;
  RefFace *face_max_1=NULL, *face_max_2=NULL;
  RefFace *face_min_angle_1=NULL, *face_min_angle_2=NULL;
  RefFace *face_max_angle_1=NULL, *face_max_angle_2=NULL;
  
  
  int counter = 0;
       
  for( ii = 0; ii< ref_volumes.size(); ii++ )
  {
    curr_volume = ref_volumes.get_and_step();
    GeomMeasureTool::measure_volume_volume_and_hydraulic_radius(curr_volume,
                                                                volume_volume,
                                                                volume_hydraulic_radius);
    if ( curr_volume->num_ref_faces() > 1 )
    {
      GeomMeasureTool::angles_between_volume_surfaces(curr_volume,
                                                      curr_min_angle,
                                                      curr_max_angle,
                                                      face_min_1, face_min_2,
                                                      face_max_1, face_max_2);
      GeomMeasureTool::find_adjacent_face_ratios(curr_volume,
                                                 curr_face_ratio,
                                                 big_face, small_face);
    }
    total_volume += volume_volume;
    counter++;
    if ( volume_volume < min_volume )
    {
      volume_min_volume = curr_volume;
      min_volume = volume_volume;
    }
    if ( volume_hydraulic_radius < min_volume_hydraulic_radius )
    {
      volume_min_hydraulic_radius = curr_volume;
      min_volume_hydraulic_radius = volume_hydraulic_radius;
    }
    if ( curr_volume->num_ref_faces() > 1 &&  curr_min_angle < min_angle )
    {
      volume_min_angle = curr_volume;
      face_min_angle_1 = face_min_1;
      face_min_angle_2 = face_min_2;
      min_angle = curr_min_angle;
    }
    if ( curr_volume->num_ref_faces() > 1 && curr_max_angle > max_angle )
    {
      volume_max_angle = curr_volume;
      max_angle = curr_max_angle;
      face_max_angle_1 = face_max_1;
      face_max_angle_2 = face_max_2;
    }
    if ( curr_volume->num_ref_faces() > 1 && curr_face_ratio > max_face_ratio )
    {
      volume_face_ratio = curr_volume;
      max_face_ratio = curr_face_ratio;
      max_small_face = small_face;
      max_big_face = big_face;
    }
  }
  if ( counter > 0 )
  {
    PRINT_INFO("*************Volume Summary************\n\n");
    PRINT_INFO("Measured %d Volumes: Average Volume was %f.\n\n", 
               counter, total_volume/counter);
    PRINT_INFO("Minimum Volume was %f in volume %d\n\n",
               min_volume, volume_min_volume->id());
    PRINT_INFO("Minimum Volume Hydraulic Radius (volume/surface area) was %f in volume %d\n\n",
               min_volume_hydraulic_radius, volume_min_hydraulic_radius->id());
    if ( volume_min_angle )
      PRINT_INFO("Minimum Interior Angle Between Two Surfaces was %f,\n"
                 "\tsurfaces %d and %d in volume %d\n\n",
                 min_angle, face_min_angle_1->id(), face_min_angle_2->id(), volume_min_angle->id());
    if ( volume_max_angle )
      PRINT_INFO("Maximum Interior Angle Between Two Surfaces was %f,\n"
                 "\tsurfaces %d and %d in volume %d\n\n",
                 max_angle, face_max_angle_1->id(), face_max_angle_2->id(), volume_max_angle->id());
    if ( volume_face_ratio )
      PRINT_INFO("Maximum Area Ratios Between Adjacent Faces was %f,\n"
                 "\tsurfaces %d and %d in volume %d.\n\n",
                 max_face_ratio, max_big_face->id(), max_small_face->id(), volume_face_ratio->id());
    int merged, unmerged;
    double ratio;     
    GeomMeasureTool::merged_unmerged_surface_ratio(ref_volumes, merged, unmerged, ratio);
    PRINT_INFO("Total number of merged surfaces:\t%d\n"
               "Total number of unmerged surfaces:\t%d\n"
               "Ratio of merged to unmerged surfaces:\t%f\n\n", merged, unmerged, ratio);
    RefGroup *group_of_shells;
    int number_shells;
    GeomMeasureTool::find_shells(ref_volumes, group_of_shells, number_shells);
    if ( number_shells > 0 )
    {
      if ( number_shells == 1 )
      {
        PRINT_INFO("Found %d shell of connected surfaces.\n", number_shells);
        PRINT_INFO("These shell is contained in the group %s\n",
                   group_of_shells->entity_name().c_str());
      }
      else
      {
        PRINT_INFO("Found %d shells of connected surfaces.\n", number_shells);
        PRINT_INFO("These shells are contained in the group %s\n",
                   group_of_shells->entity_name().c_str());
      }
    }
    PRINT_INFO("\n\n");
    PRINT_INFO("Now Printing Summaries of Surfaces contained by these volumes.\n");
    GeomMeasureTool::print_surface_measure_summary(ref_faces);
  }
  return;

}

// Find all of the surfaces in the given volumes that have narrow regions.
void GeomMeasureTool::find_surfs_with_narrow_regions(DLIList <RefVolume*> &ref_vols,
                                          double tol,
                                          DLIList <RefFace*> &surfs_with_narrow_regions)
{
  int j;

  int ii, jj;
  DLIList <RefFace*> ref_faces, temp_faces;
  RefVolume *ref_vol;
  RefFace *curr_face;
  for ( ii = 0; ii < ref_vols.size(); ii++ )
  {
    DLIList<RefFace*> faces;
    ref_vol = ref_vols.get_and_step();
    ref_vol->ref_faces(faces);
    for ( jj = faces.size(); jj > 0; jj-- )
    {
      curr_face = faces.get_and_step();
      curr_face->marked(0);
      temp_faces.append(curr_face);
    }
  }

  //uniquely add the faces.
  for ( jj = temp_faces.size(); jj > 0; jj-- )
  {
    curr_face = temp_faces.get_and_step();
    if ( curr_face->marked()== 0 )
    {
      curr_face->marked(1);
      ref_faces.append(curr_face);
    }
  }

  int num_faces = ref_faces.size();

  ProgressTool *progress_ptr = NULL;
  if (num_faces > 20)
  {
    progress_ptr = AppUtil::instance()->progress_tool();
    assert(progress_ptr != NULL);
    progress_ptr->start(0, 100, "Finding Surfaces with Narrow Regions", 
      NULL, CUBIT_TRUE, CUBIT_TRUE);
  }

  int total_faces = 0;
  double curr_percent = 0.0;

  for(j=0; j<num_faces; j++)
  {
    DLIList<CubitVector> split_pos1_list;
    DLIList<CubitVector> split_pos2_list;
    RefFace *cur_face = ref_faces.get_and_step();
    total_faces++;
    if ( progress_ptr != NULL )
    {
      curr_percent = ((double)(total_faces))/((double)(num_faces));
      progress_ptr->percent(curr_percent);
    }

    if ( CubitMessage::instance()->Interrupt() )
    {
        //interrpt.  We need to exit.
      if ( progress_ptr != NULL )
        progress_ptr->end();
        //just leave what has been calculated...
      return;
    }
    find_split_points_for_narrow_regions(cur_face,
      tol, split_pos1_list, split_pos2_list); 
    if(split_pos1_list.size() > 0)
      surfs_with_narrow_regions.append_unique(cur_face);
  }

  if ( progress_ptr != NULL )
    progress_ptr->end();
}

bool GeomMeasureTool::is_surface_narrow(RefFace *face, double small_curve_size)
{
  bool ret = true;
  DLIList<RefEdge*> edges;
  face->ref_edges(edges);
  RefVolume *vol = face->ref_volume();
  int i, j;
  double dist_sq = small_curve_size*small_curve_size;
  double proj_dist = 1.2 * small_curve_size;
  for(i=edges.size(); i>0 && ret == true; i--)
  {
    RefEdge *cur_edge = edges.get_and_step();
    double edge_length = cur_edge->measure();
    if(edge_length > small_curve_size)
    {
      int num_incs = (int)(edge_length/small_curve_size) + 1;
      double start, end;
      cur_edge->get_param_range(start, end);
      double t = start;
      bool one_bad = false;
      for(j=0; j<num_incs && ret == true; j++)
      {
        CubitVector pos1, tangent;
        cur_edge->position_from_u(t, pos1);
        cur_edge->tangent(pos1, tangent, face);
        CubitVector norm = face->normal_at(pos1, vol);
        CubitVector indir = norm * tangent;
        indir.normalize();
        CubitVector new_pos = pos1 + proj_dist * indir;
        CubitVector pt_on_surf;
        face->get_surface_ptr()->closest_point_trimmed(new_pos, pt_on_surf);
        if((pt_on_surf-pos1).length_squared() < dist_sq)
        {
          one_bad = false;
        }
        else // we found one out of small_curve range
        {
          if(one_bad)  // if we had already found one out of range this makes two in a row
            ret = false;
          else  // this is the first one out of range
            one_bad = true;
        }
      }
    }
  }

  return ret;
}

void GeomMeasureTool::find_split_points_for_narrow_regions(RefFace *face,
                                                          double size, 
                                                          DLIList<CubitVector> &split_pos1_list,
                                                          DLIList<CubitVector> &split_pos2_list)
{
  int k, i, j;
  double size_sq = size*size;
  DLIList<RefEdge*> edges;
  face->ref_edges(edges);
  while(edges.size() > 1)
  {
    RefEdge *cur_edge = edges.extract();
    for(k=edges.size(); k--;)
    {
      RefEdge *other_edge = edges.get_and_step();
      DLIList<CubitVector*> e1_pos_list, e2_pos_list;
      DLIList<RefVertex*> e1_vert_list, e2_vert_list;
      if(narrow_region_exists(cur_edge, other_edge, face, size,
        e1_pos_list, e2_pos_list, e1_vert_list, e2_vert_list))
      {
        e1_pos_list.reset();
        e2_pos_list.reset();
        e1_vert_list.reset();
        e2_vert_list.reset();

        // Loop through each pair of positions defining a split.
        for(i=e1_pos_list.size(); i--;)
        {
          int do_the_split = 1;
          RefVertex *e1_vert = e1_vert_list.get_and_step();
          RefVertex *e2_vert = e2_vert_list.get_and_step();
          CubitVector *e1_pos = e1_pos_list.get_and_step();
          CubitVector *e2_pos = e2_pos_list.get_and_step();

          // Snap to existing vertices if we are not already at at vertex.
          if(!e1_vert)
          {
            if((cur_edge->start_vertex()->coordinates() - *e1_pos).length_squared() < size_sq)
              e1_vert = cur_edge->start_vertex();
            else if((cur_edge->end_vertex()->coordinates() - *e1_pos).length_squared() < size_sq)
              e1_vert = cur_edge->end_vertex();
          }
          if(!e2_vert)
          {
            if((other_edge->start_vertex()->coordinates() - *e2_pos).length_squared() < size_sq)
              e2_vert = other_edge->start_vertex();
            else if((other_edge->end_vertex()->coordinates() - *e2_pos).length_squared() < size_sq)
              e2_vert = other_edge->end_vertex();
          }

          // We may have multiple edges separating these two vertices but the accumulated
          // length of them may still be within our narrow size so check this.  If this
          // is the case we will not want to consider this as a place to split and 
          // will as a result discard these positions.
          if(e1_vert && e2_vert)
          {
            double dist = size*sqrt(2.0);
            RefVertex *cur_vert = e1_vert;
            RefEdge *edge = cur_edge;
            double length = 0.0;
            while(edge && cur_vert != e2_vert && length <= dist)
            {
              edge = edge->get_other_curve(cur_vert, face);
              if(edge)
              {
                length += edge->get_arc_length();
                cur_vert = edge->other_vertex(cur_vert);
              }
            }
            if(length <= dist)
              do_the_split = 0;
            else
            {
              // We want to keep this split.
              split_pos1_list.append(e1_vert->coordinates());
              split_pos2_list.append(e2_vert->coordinates());
            }
          }
          else
          {
            // We want to keep this split.
            split_pos1_list.append(*e1_pos);
            split_pos2_list.append(*e2_pos);
          }
        }
      }
      while(e1_pos_list.size())
        delete e1_pos_list.pop();
      while(e2_pos_list.size())
        delete e2_pos_list.pop();
    }
  }

  // Make splits unique
  DLIList<CubitVector> unique_list1, unique_list2;
  split_pos1_list.reset();
  split_pos2_list.reset();
  for(i=split_pos1_list.size(); i--;)
  {
    CubitVector p1 = split_pos1_list.get_and_step();
    CubitVector p2 = split_pos2_list.get_and_step();
    int unique = 1;
    for(j=unique_list1.size(); j>0 && unique; j--)
    {
      CubitVector u1 = unique_list1.get_and_step();
      CubitVector u2 = unique_list2.get_and_step();
      if((p1.about_equal(u1) && p2.about_equal(u2)) ||
        (p1.about_equal(u2) && p2.about_equal(u1)))
      {
          unique = 0;
      }
    }
    if(unique)
    {
      unique_list1.append(p1);
      unique_list2.append(p2);
    }
  }
  split_pos1_list = unique_list1;
  split_pos2_list = unique_list2;
}

// Checks to see if at the given position the two edges are close together and 
// have the same tangent.
int GeomMeasureTool::is_narrow_region_at_point(RefEdge *e1,
                                               RefFace *face,
                                               const CubitVector &pt_on_e1,
                                               RefEdge *e2,
                                               const double &tol_sq,
                                               CubitVector &closest)
{
  int ret = 0;

  CubitVector tan_1, tan_2;
  e2->closest_point_trimmed(pt_on_e1, closest);
  double dist = (pt_on_e1-closest).length_squared();
  if(dist < tol_sq)
  {
    DLIList<CoEdge*> coes;
    e1->tangent(pt_on_e1, tan_1);
    e2->tangent(closest, tan_2);
    e1->get_co_edges(coes, face);
    if(coes.size() == 1)
    {
      if(coes.get()->get_sense() == CUBIT_REVERSED)
        tan_1 *= -1.0;
      coes.clean_out();
      e2->get_co_edges(coes, face);
      if(coes.size() == 1)
      {
        if(coes.get()->get_sense() == CUBIT_REVERSED)
          tan_2 *= -1.0;
        tan_1.normalize();
        tan_2.normalize();
        if(tan_1 % tan_2 < -0.9)
          ret = 1;
      }
    }
  }
  return ret;
}

int GeomMeasureTool::narrow_region_exists(RefFace *face,
                                          const double &tol)
{
  int k, ret = 0;
  DLIList<RefEdge*> edges;
  face->ref_edges(edges);
  int num_curves = edges.size();
  int num_small_curves = 0;
  while(edges.size() > 1 && !ret)
  {
    // Remove the current edge each time so that we aren't
    // doing redundant comparisons.
    RefEdge *cur_edge = edges.extract();
    if(cur_edge->get_arc_length() < tol)
      num_small_curves++;

    // Compare this edge with the remaining edges on the face.
    for(k=edges.size(); k && !ret; k--)
    {
      RefEdge *other_edge = edges.get_and_step();

      DLIList<CubitVector*> e1_pos_list, e2_pos_list;
      DLIList<RefVertex*> e1_vert_list, e2_vert_list;
      ret = narrow_region_exists(cur_edge, other_edge, face, tol,
        e1_pos_list, e2_pos_list, e1_vert_list, e2_vert_list);
      while(e1_pos_list.size())
        delete e1_pos_list.pop();
      while(e2_pos_list.size())
        delete e2_pos_list.pop();
    }
  }
  if(!ret)
  {
    if(edges.size() == 1 && edges.get()->get_arc_length() < tol)
      num_small_curves++;
  }

  ret = ret || (num_small_curves == num_curves);

  return ret;
}

int GeomMeasureTool::narrow_region_exists(RefFace *face,
                                          RefEdge *edge,
                                          const double &tol)
{
  int k, ret = 0;
  DLIList<RefEdge*> edges;
  face->ref_edges(edges);
  if(edges.move_to(edge))
  {
    edges.extract();

    // Compare this edge with the remaining edges on the face.
    for(k=edges.size(); k && !ret; k--)
    {
      RefEdge *other_edge = edges.get_and_step();

      DLIList<CubitVector*> e1_pos_list, e2_pos_list;
      DLIList<RefVertex*> e1_vert_list, e2_vert_list;
      ret = narrow_region_exists(edge, other_edge, face, tol,
        e1_pos_list, e2_pos_list, e1_vert_list, e2_vert_list);
      while(e1_pos_list.size())
        delete e1_pos_list.pop();
      while(e2_pos_list.size())
        delete e2_pos_list.pop();
    }
  }
  return ret;
}

int GeomMeasureTool::narrow_region_exists(
                                            RefEdge *e1,
                                            RefEdge *e2,
                                            RefFace *face,
                                            const double &tol,
                                            DLIList<CubitVector*> &e1_pos_list,
                                            DLIList<CubitVector*> &e2_pos_list,
                                            DLIList<RefVertex*> &e1_vert_list,
                                            DLIList<RefVertex*> &e2_vert_list)
{
  int ret = 0;
  double tol_sq = tol*tol;
  double max_dist_sq = 0.0;
  RefVertex *e1_start_vert = e1->start_vertex();
  RefVertex *e1_end_vert = e1->end_vertex();

  CubitVector closest;
  DLIList<RefVertex*> e1_verts, e2_verts;

  e1->ref_vertices(e1_verts);
  e2->ref_vertices(e2_verts);
  e1_verts.intersect_unordered(e2_verts);
  int num_shared_verts = e1_verts.size();
  RefVertex *shared_vert = NULL;
  RefEdge *edge1 = NULL;
  RefEdge *edge2 = NULL;
  if(num_shared_verts == 1)
  {
    shared_vert = e1_verts.get();
    DLIList<CoEdge*> coes;
    e1->get_co_edges(coes, face);
    if(coes.size() == 1)
    {
      RefVolume *vol = face->ref_volume();
      CubitSense facevol_sense = face->sense(vol);
      if((coes.get()->start_vertex() == shared_vert) ==
        (facevol_sense == CUBIT_FORWARD))
      {
        edge1 = e1;
        edge2 = e2;
      }
      else
      {
        edge1 = e2;
        edge2 = e1;
      }
    }
  }

  // Project cur endpoints onto other.
  int do_narrow_region_check = 1;
  if(num_shared_verts == 1 && shared_vert == e1_start_vert)
  {
    // Edges are next to each other.  Check the angle between them
    // before doing anything else.
    double interior_angle = edge1->angle_between(edge2, face);
    if(interior_angle > CUBIT_PI/4.0)
      do_narrow_region_check = 0;
  }
  if(do_narrow_region_check &&
     is_narrow_region_at_point(e1, face, e1_start_vert->coordinates(), e2, tol_sq, closest))
  {
    max_dist_sq = (closest - e1_start_vert->coordinates()).length_squared();
    e1_pos_list.append(new CubitVector(e1_start_vert->coordinates()));
    e2_pos_list.append(new CubitVector(closest));
    e1_vert_list.append(e1_start_vert);
    e2_vert_list.append(NULL);
  }
  do_narrow_region_check = 1;
  if(num_shared_verts == 1 && shared_vert == e1_end_vert)
  {
    // Edges are next to each other.  Check the angle between them
    // before doing anything else.
    double interior_angle = edge1->angle_between(edge2, face);
    if(interior_angle > CUBIT_PI/4.0)
      do_narrow_region_check = 0;
  }
  if(do_narrow_region_check &&
     is_narrow_region_at_point(e1, face, e1_end_vert->coordinates(), e2, tol_sq, closest))
  {
    double cur_dist_sq = (closest - e1_end_vert->coordinates()).length_squared();
    max_dist_sq = max_dist_sq > cur_dist_sq ? max_dist_sq : cur_dist_sq;
    e1_pos_list.append(new CubitVector(e1_end_vert->coordinates()));
    e2_pos_list.append(new CubitVector(closest));
    e1_vert_list.append(e1_end_vert);
    e2_vert_list.append(NULL);
  }

  if(e1_pos_list.size() < 2)
  {
    RefVertex *e2_start_vert = e2->start_vertex();
    RefVertex *e2_end_vert = e2->end_vertex();
    do_narrow_region_check = 1;
    if(num_shared_verts == 1 && shared_vert == e2_start_vert)
    {
      // Edges are next to each other.  Check the angle between them
      // before doing anything else.
      double interior_angle = edge1->angle_between(edge2, face);
      if(interior_angle > CUBIT_PI/4.0)
        do_narrow_region_check = 0;
    }
    if(do_narrow_region_check &&
       is_narrow_region_at_point(e2, face, e2_start_vert->coordinates(), e1, tol_sq, closest))
    {
      if(e1_pos_list.size() == 0 || !closest.about_equal(*e1_pos_list.get()))
      {
        double cur_dist_sq = (closest - e2_start_vert->coordinates()).length_squared();
        max_dist_sq = max_dist_sq > cur_dist_sq ? max_dist_sq : cur_dist_sq;
        e2_pos_list.append(new CubitVector(e2_start_vert->coordinates()));
        e1_pos_list.append(new CubitVector(closest));
        e2_vert_list.append(e2_start_vert);
        e1_vert_list.append(NULL);
      }
    }
    if(e1_pos_list.size() < 2)
    {
      do_narrow_region_check = 1;
      if(num_shared_verts == 1 && shared_vert == e2_end_vert)
      {
        // Edges are next to each other.  Check the angle between them
        // before doing anything else.
        double interior_angle = edge1->angle_between(edge2, face);
        if(interior_angle > CUBIT_PI/4.0)
          do_narrow_region_check = 0;
      }
      if(do_narrow_region_check &&
         is_narrow_region_at_point(e2, face, e2_end_vert->coordinates(), e1, tol_sq, closest))
      {
        if(e1_pos_list.size() == 0 || !closest.about_equal(*e1_pos_list.get()))
        {
          double cur_dist_sq = (closest - e2_end_vert->coordinates()).length_squared();
          max_dist_sq = max_dist_sq > cur_dist_sq ? max_dist_sq : cur_dist_sq;
          e2_pos_list.append(new CubitVector(e2_end_vert->coordinates()));
          e1_pos_list.append(new CubitVector(closest));
          e2_vert_list.append(e2_end_vert);
          e1_vert_list.append(NULL);
        }
      }
    }
  }

  if(e1_pos_list.size() == 2)
  {
    int w;
    int all_good = 1;
 //   double dist_tol = sqrt(max_dist_sq)*5.0;
    e1_pos_list.reset();
    e2_pos_list.reset();
    CubitVector *cur1 = e1_pos_list.get_and_step();
    CubitVector *cur2 = e1_pos_list.get();
    CubitVector *other1 = e2_pos_list.get_and_step();
    CubitVector *other2 = e2_pos_list.get();

    double len1 = e1->get_arc_length(*cur1, *cur2);
    if(len1 > tol)
    {
      double len2 = e2->get_arc_length(*other1, *other2);
      if(len2 > tol)
      {
        double cur_param1 = e1->u_from_position(*cur1);
        double cur_param2 = e1->u_from_position(*cur2);
        int num_divisions = 2;
        CubitVector cur_pos;
        double param_step = (cur_param2-cur_param1)/num_divisions;
        double cur_param = cur_param1 + param_step;
        for(w=1; w<num_divisions && all_good; w++)
        {
          e1->position_from_u(cur_param, cur_pos);
          cur_param += param_step;
          if(is_narrow_region_at_point(e1, face, cur_pos, e2, tol_sq, closest))
          {
            // Sanity check to make sure we aren't splitting off negative space.
            CubitVector mid = (cur_pos + closest)/2.0;
            CubitVector tmp_pt;
            face->get_surface_ptr()->closest_point_trimmed(mid, tmp_pt);
            if(!mid.about_equal(tmp_pt))
            {
              CubitVector norm = face->normal_at(tmp_pt);
              CubitVector dir(tmp_pt - mid);
              dir.normalize();
              if(fabs(norm % dir) < .9)
                all_good = 0;
            }
          }
        }
      }
      else
        all_good = 0;
    }
    else 
      all_good = 0;

    if(all_good)
      ret = 1;
    else
    {
      e1_pos_list.remove(cur1);
      e1_pos_list.remove(cur2);
      e2_pos_list.remove(other1);
      e2_pos_list.remove(other2);
      delete cur1;
      delete cur2;
      delete other1;
      delete other2;
    }
  }

  if(!ret && e1_pos_list.size() > 0)
  {
    int i;
    e1_pos_list.reset();
    e2_pos_list.reset();
    e1_vert_list.reset();
    e2_vert_list.reset();
    for(i=e1_pos_list.size(); i--;)
    {
      RefVertex *e1_vert = e1_vert_list.get_and_step();
      RefVertex *e2_vert = e2_vert_list.get_and_step();

      RefVertex *cur_vert = NULL;
      RefEdge *cur_edge = NULL, *other_edge = NULL;
      if(e1_vert)
      {
        cur_edge = e1;
        other_edge = e2;
        cur_vert = e1_vert;
      }
      else if(e2_vert)
      {
        cur_edge = e2;
        other_edge = e1;
        cur_vert = e2_vert;
      }
      if(cur_vert)
      {
        CubitVector next_pos;
        CubitVector prev_pos = cur_vert->coordinates();
        int num_incs = 20;
        double step = cur_edge->get_arc_length()/(double)num_incs;
        int still_good = 1;
        int cntr = 0;
        double direction = (cur_vert == cur_edge->start_vertex() ? 1.0 : -1.0);
        // Do coarse traversal along curve to see where we start deviating from
        // narrow.
        while(still_good)
        {
          cur_edge->point_from_arc_length(prev_pos, step*direction, next_pos);
          if(is_narrow_region_at_point(cur_edge, face, next_pos, other_edge, tol_sq, closest) &&
                      cntr < num_incs)
          {
            prev_pos = next_pos;
            cntr++;
          }
          else
            still_good = 0;
        }
        if(cntr < num_incs)
        {
          cntr = 0;
          double cur_arc_length = cur_edge->get_arc_length(prev_pos, next_pos);
          // Do bisection on remaining interval to zero in on point where
          // we go from narrow to non-narrow.
          CubitVector mid_pos;
          double close_enough = tol/20.0;
          while(cur_arc_length > close_enough && cntr < 100)
          {
            cntr++;  // prevent infinite looping
            cur_edge->point_from_arc_length(prev_pos, cur_arc_length*direction/2.0, mid_pos);
            if(is_narrow_region_at_point(cur_edge, face, mid_pos, other_edge, tol_sq, closest))
              prev_pos = mid_pos;
            else
              next_pos = mid_pos;
            cur_arc_length = cur_edge->get_arc_length(prev_pos, next_pos);
          }
          if(cur_edge->get_arc_length(cur_vert->coordinates(), prev_pos) > tol)
          {
            // end up with the position that guarantees a new split curve that
            // is smaller than the small curve size.
            other_edge->closest_point_trimmed(prev_pos, closest);
            if(cur_edge == e1)
            {
              e1_pos_list.append(new CubitVector(prev_pos));
              e2_pos_list.append(new CubitVector(closest));
            }
            else
            {
              e2_pos_list.append(new CubitVector(prev_pos));
              e1_pos_list.append(new CubitVector(closest));
            }
            e1_vert_list.append(NULL);
            e2_vert_list.append(NULL);
            ret = 1;
          }
        }
      }
    }
  }

  if(ret)
  {
    int i, j;
    e1_pos_list.reset();
    e2_pos_list.reset();
    e1_vert_list.reset();
    e2_vert_list.reset();
    for(i=e1_pos_list.size(); i--;)
    {
      CubitVector *e1_pos = e1_pos_list.get();
      CubitVector *e2_pos = e2_pos_list.get();
      int num_divisions = 6;
      CubitVector step = (*e2_pos - *e1_pos)/num_divisions;
      CubitVector cur_pos = *e1_pos + step;
      int removed = 0;
      for(j=1; j<num_divisions; j++)
      {
        CubitVector tmp_pt;
        face->get_surface_ptr()->closest_point_trimmed(cur_pos, tmp_pt);
        if(!cur_pos.about_equal(tmp_pt))
        {
          CubitVector norm = face->normal_at(tmp_pt);
          CubitVector dir(tmp_pt - cur_pos);
          dir.normalize();
          if(fabs(norm % dir) < .9)
          {
            removed = 1;
            j=num_divisions;
            e1_pos_list.remove();
            e2_pos_list.remove();
            e1_vert_list.remove();
            e2_vert_list.remove();
          }
        }
        else
          cur_pos += step;
      }
      if(!removed)
      {
        e1_pos_list.step();
        e2_pos_list.step();
        e1_vert_list.step();
        e2_vert_list.step();
      }
    }

    if(e1_pos_list.size() == 0)
      ret = 0;
  }

  return ret;
}

void GeomMeasureTool::find_small_curves( DLIList <RefVolume*> &ref_vols,
                                         double tol,
                                         DLIList <RefEdge*> &small_curves,
                                         DLIList <double> &small_lengths)
{
  int ii, jj;
  DLIList <RefEdge*> ref_edges, temp_edges;
  RefVolume *ref_vol;
  RefEdge *curr_edge;
  for ( ii = 0; ii < ref_vols.size(); ii++ )
  {
    ref_vol = ref_vols.get_and_step();
    ref_vol->ref_edges(temp_edges);
      //uniquely add the edges.
    for ( jj = temp_edges.size(); jj > 0; jj-- )
    {
      curr_edge = temp_edges.pop();
      if ( curr_edge->marked()== 0 )
      {
        curr_edge->marked(1);
        ref_edges.append(curr_edge);
      }
    }
  }

  int num_curves = ref_edges.size();
  ProgressTool *progress_ptr = NULL;
  if (num_curves> 20)
  {
    progress_ptr = AppUtil::instance()->progress_tool();
    assert(progress_ptr != NULL);
    progress_ptr->start(0, 100, "Finding Small Curves", 
      NULL, CUBIT_TRUE, CUBIT_TRUE);
  }

  int total_curves = 0;
  double curr_percent = 0.0;
  double length;
    //find the small curves and reset the marked flag.
  for ( ii = ref_edges.size(); ii > 0; ii-- )
  {
    total_curves++;
    if ( progress_ptr != NULL )
    {
      curr_percent = ((double)(total_curves))/((double)(num_curves));
      progress_ptr->percent(curr_percent);
    }
    if ( CubitMessage::instance()->Interrupt() )
    {
        //interrpt.  We need to exit.
      if ( progress_ptr != NULL )
        progress_ptr->end();
        //just leave what has been calculated...
      return;
    }

    curr_edge = ref_edges.get_and_step();

    //skip point curves
    if( curr_edge->geometry_type() == POINT_CURVE_TYPE )
      continue;

    curr_edge->marked(0);
    length = curr_edge->measure();
    if ( length <= tol )
    {
      small_curves.append(curr_edge);
      small_lengths.append(length);
    }
  }

  if ( progress_ptr != NULL )
    progress_ptr->end();

  return;
}

RefEdge* GeomMeasureTool::find_first_small_curve(RefVolume* vol,
                                         double tol)
{
  RefEdge *ret = NULL;
  int j;
  DLIList <RefEdge*> ref_edges;
  vol->ref_edges(ref_edges);
  for(j=ref_edges.size(); j > 0 && !ret; j--)
  {
    RefEdge *curr_edge = ref_edges.get_and_step();
    if(curr_edge->measure() <= tol)
      ret = curr_edge;
  }
  return ret;
}

void GeomMeasureTool::find_narrow_faces(DLIList<RefVolume*> &ref_vols,
                                        double small_curve_size,
                                        DLIList<RefFace*> &narrow_faces,
                                        DLIList<RefFace*> &surfs_to_ignore)
{
  int ii, jj;
  DLIList <RefFace*> ref_faces, temp_faces;
  RefVolume *ref_vol;
  RefFace *curr_face;

  for ( ii = 0; ii < ref_vols.size(); ii++ )
  {
    DLIList<RefFace*> faces;
    ref_vol = ref_vols.get_and_step();
    ref_vol->ref_faces(faces);
    for ( jj = faces.size(); jj > 0; jj-- )
    {
      curr_face = faces.get_and_step();
      curr_face->marked(0);
      temp_faces.append(curr_face);
    }
  }

  //uniquely add the faces.
  for ( jj = temp_faces.size(); jj > 0; jj-- )
  {
    curr_face = temp_faces.get_and_step();
    if ( curr_face->marked()== 0 )
    {
      curr_face->marked(1);
      ref_faces.append(curr_face);
    }
  }

  int num_faces = ref_faces.size();
  ProgressTool *progress_ptr = NULL;
  if (num_faces > 20)
  {
    progress_ptr = AppUtil::instance()->progress_tool();
    assert(progress_ptr != NULL);
    progress_ptr->start(0, 100, "Finding Narrow Surfaces", 
      NULL, CUBIT_TRUE, CUBIT_TRUE);
  }

  int total_faces = 0;
  double curr_percent = 0.0;

  for ( ii = ref_faces.size(); ii > 0; ii-- )
  {
    total_faces++;
    if ( progress_ptr != NULL )
    {
      curr_percent = ((double)(total_faces))/((double)(num_faces));
      progress_ptr->percent(curr_percent);
    }

    if ( CubitMessage::instance()->Interrupt() )
    {
        //interrpt.  We need to exit.
      if ( progress_ptr != NULL )
        progress_ptr->end();
        //just leave what has been calculated...
      return;
    }

    curr_face = ref_faces.get_and_step();
    if(!surfs_to_ignore.is_in_list(curr_face))
    {
      if(narrow_region_exists(curr_face, small_curve_size))
      {
        DLIList<CubitVector> split_pos1_list;
        DLIList<CubitVector> split_pos2_list;
        find_split_points_for_narrow_regions(curr_face,
          small_curve_size, split_pos1_list, split_pos2_list); 
        if(split_pos1_list.size() == 0)
          narrow_faces.append_unique(curr_face);
      }
    }
  }

  if ( progress_ptr != NULL )
    progress_ptr->end();

  return;
}

void GeomMeasureTool::find_small_faces( DLIList <RefVolume*> &ref_vols,
                                        double tol,
                                        DLIList <RefFace*> &small_faces)
{
  int ii, jj;
  DLIList <RefFace*> ref_faces, temp_faces;
  RefVolume *ref_vol;
  RefFace *curr_face;
  for ( ii = 0; ii < ref_vols.size(); ii++ )
  {
    DLIList<RefFace*> faces;
    ref_vol = ref_vols.get_and_step();
    ref_vol->ref_faces(faces);
    for ( jj = faces.size(); jj > 0; jj-- )
    {
      curr_face = faces.get_and_step();
      curr_face->marked(0);
      temp_faces.append(curr_face);
    }
  }

  //uniquely add the faces.
  for ( jj = temp_faces.size(); jj > 0; jj-- )
  {
    curr_face = temp_faces.get_and_step();
    if ( curr_face->marked()== 0 )
    {
      curr_face->marked(1);
      ref_faces.append(curr_face);
    }
  }

  int num_faces = ref_faces.size();
  ProgressTool *progress_ptr = NULL;
  if (num_faces > 20)
  {
    progress_ptr = AppUtil::instance()->progress_tool();
    assert(progress_ptr != NULL);
    progress_ptr->start(0, 100, "Finding Small Surfaces", 
      NULL, CUBIT_TRUE, CUBIT_TRUE);
  }

  int total_faces = 0;
  double curr_percent = 0.0;
  double area;
    //find the small curves and reset the marked flag.
  for ( ii = ref_faces.size(); ii > 0; ii-- )
  {
    total_faces++;
    if ( progress_ptr != NULL )
    {
      curr_percent = ((double)(total_faces))/((double)(num_faces));
      progress_ptr->percent(curr_percent);
    }

    if ( CubitMessage::instance()->Interrupt() )
    {
        //interrpt.  We need to exit.
      if ( progress_ptr != NULL )
        progress_ptr->end();
        //just leave what has been calculated...
      return;
    }

    curr_face = ref_faces.get_and_step();
    area = measure_area(curr_face);
    if ( area <= tol )
      small_faces.append(curr_face);
  }

  if ( progress_ptr != NULL )
    progress_ptr->end();

  return;
}

void GeomMeasureTool::find_small_faces_hydraulic_radius( DLIList <RefVolume*> &ref_vols,
                                                         double tol,
                                                         DLIList <RefFace*> &small_faces,
                                                         DLIList <double> &small_hyd_rad,
                                                         double &percent_planar,
                                                         double &percent_pl_co)
{
  int ii, jj;
  DLIList <RefFace*> ref_faces, temp_faces;
  RefVolume *ref_vol;
  RefFace *curr_face;
  for ( ii = 0; ii < ref_vols.size(); ii++ )
  {
    ref_vol = ref_vols.get_and_step();
    ref_vol->ref_faces(temp_faces);
      //uniquely add the faces.
    for ( jj = temp_faces.size(); jj > 0; jj-- )
    {
      curr_face = temp_faces.pop();
      if ( curr_face->marked()== 0 )
      {
        curr_face->marked(1);
        ref_faces.append(curr_face);
      }
    }
  }
  double area;
  double length;
  int total_faces = 0;
  int total_plane = 0;
  int total_cone = 0;
  double area_plane = 0.0;
  double total_area = 0.0;
    //find the small curves and reset the marked flag.
  DLIList <CoEdge*> co_edges;
  int num_faces = ref_faces.size();
  ProgressTool *progress_ptr = NULL;
  if (num_faces > 20)
  {
    progress_ptr = AppUtil::instance()->progress_tool();
    assert(progress_ptr != NULL);
    progress_ptr->start(0, 100, "Small Surface Progress", 
      NULL, CUBIT_TRUE, CUBIT_TRUE);
  }
  double curr_percent = 0.0;
  for ( ii = ref_faces.size(); ii > 0; ii-- )
  {
    curr_face = ref_faces.get_and_step();
    curr_face->marked(0);
    area = measure_area(curr_face);
    total_area += area;
    total_faces++;
      //update the progress..
    if ( progress_ptr != NULL )
    {
      curr_percent = ((double)(total_faces))/((double)(num_faces));
      progress_ptr->percent(curr_percent);
    }
    if ( CubitMessage::instance()->Interrupt() )
    {
        //interrpt.  We need to exit.
      if ( progress_ptr != NULL )
        progress_ptr->end();
        //just leave what has been calculated...
      return;
    }
    if ( curr_face->geometry_type() == PLANE_SURFACE_TYPE )
    {
      total_plane++;
      area_plane += area;
    }
    else if ( curr_face->geometry_type() == CONE_SURFACE_TYPE )
      total_cone++;
    co_edges.clean_out();
    curr_face->co_edges(co_edges);
    length = 0.0;
    for ( jj = co_edges.size(); jj > 0; jj-- )
    {
      length += co_edges.get_and_step()->get_ref_edge_ptr()->measure();
    }
      //reset the mark.
    curr_face->marked(0);
      //compute the hydraulic radius 4*(A/P).
    if ( length <= CUBIT_RESABS )
    {
      PRINT_INFO("Total Perimeter Length of Surface %d is less than tolerance.\n",
                 curr_face->id());
      
      continue;
    }
    area = 4.0*(area/length);
    if ( area <= tol )
    {
      small_faces.append(curr_face);
      small_hyd_rad.append(area);
    }
  }
  if ( total_faces > 0 )
  {
    percent_planar = 100.0*area_plane/total_area;
    percent_pl_co = 100.0*((double)total_plane+(double)total_cone)/(double)total_faces;
  }
  else
  {
    percent_planar = 0.;
    percent_pl_co = 0.;
  }
  if ( progress_ptr != NULL )
  {
    progress_ptr->end();
  }
  
  return;
}
void GeomMeasureTool::find_small_volumes( DLIList <RefVolume*> &ref_vols,
                                          double tol,
                                          DLIList <RefVolume*> &small_volumes)
{
  int ii;
  RefVolume *ref_vol;
  double volume;
  
  for ( ii = 0; ii < ref_vols.size(); ii++ )
  {
    ref_vol = ref_vols.get_and_step();
    volume = ref_vol->measure();
    if ( volume <= tol )
      small_volumes.append(ref_vol);
  }
  return;
}
void GeomMeasureTool::find_small_volumes_hydraulic_radius( DLIList <RefVolume*> &ref_vols,
                                                           double tol,
                                                           DLIList <RefVolume*> &small_volumes,
                                                           DLIList <double> &small_hyd_rad)
{
  int ii, jj;
  DLIList <RefFace*> ref_faces, temp_faces;
  RefVolume *ref_vol;
  RefFace *curr_face;
  double area;
  double hyd_volume;
  
  for ( ii = 0; ii < ref_vols.size(); ii++ )
  {
    ref_vol = ref_vols.get_and_step();
    ref_vol->ref_faces(temp_faces);
    area = 0.0;
    for ( jj = temp_faces.size(); jj > 0; jj-- )
    {
      curr_face = temp_faces.pop();
      area += measure_area(curr_face);
    }
    if ( area <= CUBIT_RESABS )
      continue;
    hyd_volume = 6.0*(ref_vol->measure() / area);
    if ( hyd_volume <= tol )
    {
      small_hyd_rad.append(hyd_volume);
      small_volumes.append(ref_vol);
    }
  }
  return;
}
///
/// Find the tangential meetings in the volume.
/// This specifically looks for surfaces that meet tangentially
/// that would be a problem.  Usually these are surfaces that
/// come into a side 180 degrees and on top there is a
/// sharpe angle.  Usually if there isn't a sharpe curve
/// angle these meetings are not bad.
/// Note that this function assumes that you are passing
/// in sets of curve edges that have small interior angles between them.
/// It also assumes that the edge pairs are ordered as they would
/// be found in the surfaces (first edge then next edge).
/// Basically, call the funciton, find_interior_curve_angles
/// before calling this function, and pass this function those results.
///
void GeomMeasureTool::find_sharp_tangential_meets( RefVolume *ref_volume,
                                                   DLIList <RefEdge*> &small_angle_edge_pairs,
                                                   DLIList <RefFace*> &tangential_surface_pairs )
{
    //Okay, given the small angle edge pairs.  See if there are surfaces that meet tangentially.
    //I'm really looking for something spefic.
  RefEdge *ref_edge_1, *ref_edge_2;
  int ii, jj;
  DLIList <Loop*> loops_1, loops_2;
  Loop *common_loop;
  const double rad_to_deg = 180.0/CUBIT_PI;
  RefFace *ref_face;
  
    //loop over for each edge pair.
  for ( ii = 0; ii < small_angle_edge_pairs.size(); ii += 2 )
  {
    ref_edge_1 = small_angle_edge_pairs.get_and_step();
    ref_edge_2 = small_angle_edge_pairs.get_and_step();
      //find the loop that these edges are both on.
    loops_1.clean_out();
    loops_2.clean_out();
    ref_edge_1->loops(loops_1);
    ref_edge_2->loops(loops_2);
    common_loop = NULL;
    for ( jj = loops_1.size(); jj > 0; jj-- )
    {
      if ( !loops_2.move_to(loops_1.get()) )
        loops_1.remove();
      else
        loops_1.step();
    }
    if ( loops_1.size() == 1 )
      common_loop = loops_1.get();
    else
    {
        //find the loop that has a co_edge from ref_edge_1
        //then a co_edge from ref_edge_2.
      DLIList <CoEdge*> co_edges;
      Loop *tmp_loop;
      int kk;
      CoEdge *co_edge;
      for ( jj = loops_1.size(); jj > 0; jj-- )
      {
        tmp_loop = loops_1.get_and_step();
        co_edges.clean_out();
        tmp_loop->ordered_co_edges(co_edges);
        for ( kk = co_edges.size(); kk > 0; kk-- )
        {
          co_edge = co_edges.get_and_step();
          if ( co_edge->get_ref_edge_ptr() == ref_edge_1 )
          {
            if ( co_edges.get()->get_ref_edge_ptr() == ref_edge_2 )
            {
              common_loop = tmp_loop;
              break;
            }
            else
                //not this loop..
              break;
          }
        }
        if ( common_loop != NULL )
          break;
      }
    }
    assert ( common_loop != NULL );
      //now find surface.
      //Then see if there are surfaces that meet at 180 degrees that
      //are 90* rotated from here.
    ref_face = common_loop->get_ref_face_ptr();
      //Okay, refface is the surface on which we have the bad angle.
      //One of these curves will have a concave angle and the other a
      //convex angle.  The surface on the other side of the convex angle
      //we need to check to see if it comes tangentially into another
      //surface (dihedral angle of 180.)
    RefFace *ref_face_1 = ref_edge_1->other_face(ref_face, ref_volume);
    if( ref_face_1 == NULL )
      continue;

    double surf_angle_1, surf_angle_2;
    RefFace *tangential_face = NULL;
    RefEdge *edge_next_to = NULL;
    surf_angle_1 =  GeometryQueryTool::instance()->surface_angle(ref_face,
                                                                 ref_face_1,
                                                                 ref_edge_1,
                                                                 ref_volume);
    surf_angle_1 *= rad_to_deg;
    RefFace *ref_face_2;
    ref_face_2 = ref_edge_2->other_face(ref_face, ref_volume);
    if( ref_face_2 == NULL )
      continue;
    surf_angle_2 = GeometryQueryTool::instance()->surface_angle(ref_face,
                                                                ref_face_2,
                                                                ref_edge_2,
                                                                ref_volume);
    surf_angle_2 *= rad_to_deg;
    if ( surf_angle_1 < 180.0 && surf_angle_2 > 180.0 )
    {
      tangential_face = ref_face_1;
      edge_next_to = ref_edge_1;
    }
    else if ( surf_angle_2 < 180.0 && surf_angle_1 > 180.0 )
    {
      tangential_face = ref_face_2;
      edge_next_to = ref_edge_2;
    }
    if ( tangential_face == NULL ) // This isn't the config we seek.
      continue;

      //Now, on this face, there is an edge, next to edge_next_to,
      //that will have an interior angle of 180.0.  If that is the
      //case then this is what we seek.
    RefVertex* ref_vert = ref_edge_1->common_ref_vertex(ref_edge_2,
                                                        ref_face);
    DLIList <RefEdge*> tmp_ref_edges;
    ref_vert->ref_edges(tmp_ref_edges);
      //Find the edge that is on tangential face and is not edge_next_to.
    RefEdge *tangential_edge=NULL, *tmp_edge;
    for ( jj = tmp_ref_edges.size(); jj > 0; jj-- )
    {
      tmp_edge = tmp_ref_edges.get_and_step();
      if ( tmp_edge != edge_next_to &&
           tmp_edge->is_directly_related(tangential_face) )
      {
        tangential_edge = tmp_edge;
        break;
      }
    }
    if ( tangential_edge == NULL )
        //not what we are looking for.
      continue;
      //Now get the other face and measure the dihedral angle.
    RefFace *other_face = tangential_edge->other_face(tangential_face,
                                                      ref_volume);
    double tan_angle;
    tan_angle = GeometryQueryTool::instance()->surface_angle(tangential_face,
                                                             other_face,
                                                             tangential_edge,
                                                             ref_volume);
    tan_angle *= rad_to_deg;
      //If they are within tolerance, then we have it...
    if ( tan_angle > 165.0 && tan_angle < 195.0 )
    {
      tangential_surface_pairs.append(tangential_face);
      tangential_surface_pairs.append(other_face);
    }
  }
}
    
      
    
void GeomMeasureTool::find_interior_curve_angles( RefVolume *ref_volume,
                                                  double upper_bound,
                                                  double lower_bound,
                                                  DLIList <RefEdge*> &large_edge_angles,
                                                  DLIList <RefEdge*> &small_edge_angles,
                                                  DLIList <double> &large_angles,
                                                  DLIList <double> &small_angles,
                                                  int &total_interior,
                                                  int &total_fuzzy)
{
  int ii, jj, kk;
#ifdef BOYD17
  DLIList <RefFace*> ref_faces, temp_faces;
#endif
  DLIList <RefFace*> ref_faces;
  RefFace *curr_face;
  total_interior = 0;
  total_fuzzy = 0;
  ref_volume->ref_faces(ref_faces);
  
    //Okay now loop over the edge loops of the reffaces
    //and find the interior angles.
  DLIList<Loop*> loops;
  Loop *curr_loop;
  DLIList <CoEdge*> co_edges;
  CoEdge *curr_edge, *next_edge;
  double angle = 0.0;
  const double RAD_TO_DEG = 180.0/CUBIT_PI;
  for ( ii = ref_faces.size(); ii > 0; ii-- )
  {
    curr_face = ref_faces.get_and_step(); 
    curr_face->loops(loops);
      //Now loop over these loops.
    for ( jj = loops.size(); jj > 0; jj-- )
    {
      curr_loop = loops.pop();
        //now loop over the edges in this
      co_edges.clean_out();
      curr_loop->ordered_co_edges(co_edges);
      for ( kk = co_edges.size(); kk > 0; kk-- )
      {
        curr_edge = co_edges.get_and_step();
        next_edge = co_edges.get();
	   
        if ( curr_edge == next_edge )
        {
            //This is a hardpoint
          if ( curr_edge->get_ref_edge_ptr()->geometry_type() == POINT_CURVE_TYPE )
          {
            angle = 360.0;
              //for now ignore this..
            continue;
          }
          else
          {
            RefEdge *the_edge = curr_edge->get_ref_edge_ptr();
              //Get points 1% from the ends in both directions.
            CubitVector point_start, point_end;
            CubitStatus stat = the_edge->position_from_fraction(0.01,
                                                                point_start);
            if ( stat != CUBIT_SUCCESS )
              angle = 180.0;
            else
            {
              stat = the_edge->position_from_fraction(0.99,
                                                      point_end);
              if ( stat != CUBIT_SUCCESS )
                angle = 180.0;
              else
              {
                CubitVector normal = curr_face->normal_at(point_start, NULL);
                CubitVector tangent_1, tangent_2;
                stat = the_edge->tangent( point_start, tangent_1, curr_face);
                CubitStatus stat2 = the_edge->tangent( point_end, tangent_2, curr_face);
                if ( stat != CUBIT_SUCCESS || stat2 != CUBIT_SUCCESS )
                  angle = 180.0;
                else
                {
                    //one of the tangents will be pointing the wrong
                    //direction, so reverse it based on the sense.
                  if ( curr_edge->get_sense() == CUBIT_REVERSED )
                    tangent_1 = -tangent_1;
                  else
                    tangent_2 = -tangent_2;
                  angle = normal.vector_angle(tangent_2, tangent_1);
                  angle *= RAD_TO_DEG;
                }
              }
            }
          }
        }
        else if ( curr_edge->get_ref_edge_ptr() == next_edge->get_ref_edge_ptr() )
        {
            //This is the end of a hard line. Set this angle because
            //it is tricky to measure it.
          angle = 360.0;
        }
        else if ( curr_edge != next_edge )
        {
          angle =
            GeometryQueryTool::instance()->geometric_angle(curr_edge,
                                                           next_edge);
          angle *= RAD_TO_DEG;
        }
        total_interior++;
          //count the number of angles that are fuzzy.
        if ( angle <= GEOM_END_LOWER )
          total_fuzzy++;
        else if ( angle >= GEOM_END_UPPER &&
                  angle <= GEOM_SIDE_LOWER )
          total_fuzzy++;
        else if ( angle >= GEOM_SIDE_UPPER &&
                  angle <= GEOM_CORNER_LOWER )
          total_fuzzy++;
        else if ( angle >= GEOM_CORNER_UPPER )
          total_fuzzy++;
        if ( angle <= lower_bound )
        {
          small_edge_angles.append(curr_edge->get_ref_edge_ptr());
          small_edge_angles.append(next_edge->get_ref_edge_ptr());
          small_angles.append(angle);
        }
        if ( angle >= upper_bound )
        {
          large_edge_angles.append(curr_edge->get_ref_edge_ptr());
          large_edge_angles.append(next_edge->get_ref_edge_ptr());
          large_angles.append(angle);
        }
      }
    }
  }
  return;
}

void GeomMeasureTool::find_dihedral_angles( DLIList <RefVolume*> &ref_vols,
                                            double upper_bound,
                                            double lower_bound,
                                            DLIList <RefFace*> &large_face_angles,
                                            DLIList <RefFace*> &small_face_angles,
                                            DLIList <double> &large_angles,
                                            DLIList <double> &small_angles,
                                            int &total_interior,
                                            int &total_fuzzy,
                                            int &total_not_flat)
{
  int ii;
#ifdef BOYD17
  DLIList <RefFace*> ref_faces, temp_faces;
#endif
  RefVolume *ref_vol;
  total_interior = 0;
  total_fuzzy = 0;
  total_not_flat = 0;
  for ( ii = 0; ii < ref_vols.size(); ii++ )
  {
    ref_vol = ref_vols.get_and_step();
    find_dihedral_angles(ref_vol, upper_bound,
                         lower_bound,
                         large_face_angles,
                         small_face_angles,
                         large_angles, small_angles,
                         total_interior, total_fuzzy, total_not_flat);
  }
  return;
}
void GeomMeasureTool::find_dihedral_angles(RefVolume *curr_volume,
                                           double upper_bound,
                                           double lower_bound,
                                           DLIList <RefFace*> &large_face_angles,
                                           DLIList <RefFace*> &small_face_angles,
                                           DLIList <double> &large_angles,
                                           DLIList <double> &small_angles,
                                           int &total_interior,
                                           int &total_fuzzy, int &total_not_flat)
{
  int ii, jj, shared_angles = 0;
  double common_angle;
  DLIList <RefFace*> face_list;
  RefFace *curr_face_1, *curr_face_2;
  RefFace *tmp_face;
  RefEdge *common_edge;
  DLIList <RefVolume*> ref_vols;
  DLIList<RefEdge*> vol_edges;
  curr_volume->ref_edges(vol_edges);

  //maybe the body is a sheet-body...don't print errors if so
  bool is_sheet_body = false;
  if( curr_volume->measure() == 0 )
    is_sheet_body = true;

  for ( ii = vol_edges.size(); ii > 0; ii-- )
  {
    common_edge = vol_edges.get_and_step();

    //skip point curves
    if( common_edge->geometry_type() == POINT_CURVE_TYPE )
      continue;

    curr_face_1 = NULL;
    curr_face_2 = NULL;
    face_list.clean_out();
    common_edge->ref_faces(face_list);
      //find the two faces that are part of this volume.
    for ( jj = face_list.size(); jj > 0; jj-- )
    {
      tmp_face = face_list.get_and_step();
      ref_vols.clean_out();
      tmp_face->ref_volumes(ref_vols);
      if ( ref_vols.move_to(curr_volume) )
      {
        if ( curr_face_1 == NULL )
          curr_face_1 = tmp_face;
        else if ( curr_face_2 == NULL )
        {
          curr_face_2 = tmp_face;
          break;
        }
      }
    }
    if ( curr_face_1 == NULL ||
         curr_face_2 == NULL )
    {
      if( is_sheet_body == false )
        PRINT_ERROR("Problems finding connected surfaces to a curve\n"
                   "\tfor measuring dihedral angles.\n");
      continue;
    }
    shared_angles++;
    common_angle =  GeometryQueryTool::instance()->surface_angle(curr_face_1,
                                                                 curr_face_2,
                                                                 common_edge,
                                                                 curr_volume);
    common_angle *= 180.0/CUBIT_PI;
    total_interior++;
      //use hard coded angles because we want these a little looser.
    if ( common_angle >= 200.0 ||
         common_angle <= 160.0 )
      total_not_flat++;
    
    if ( curr_face_1->geometry_type() != PLANE_SURFACE_TYPE ||
         curr_face_2->geometry_type() != PLANE_SURFACE_TYPE )
      total_fuzzy++;
    else if ( common_angle <= GEOM_END_LOWER )
      total_fuzzy++;
    else if ( common_angle >= GEOM_END_UPPER &&
              common_angle <= GEOM_SIDE_LOWER )
      total_fuzzy++;
    else if ( common_angle >= GEOM_SIDE_UPPER &&
              common_angle <= GEOM_CORNER_LOWER )
      total_fuzzy++;
    else if ( common_angle >= GEOM_CORNER_UPPER )
      total_fuzzy++;
    
    if( common_angle <= lower_bound )
    {
      small_face_angles.append(curr_face_1);
      small_face_angles.append(curr_face_2);
      small_angles.append(common_angle);
    }
    if( common_angle >= upper_bound )
    {
      large_face_angles.append(curr_face_1);
      large_face_angles.append(curr_face_2);
      large_angles.append(common_angle);
    }
  }
    //clean up the edge marks.
  for ( ii = vol_edges.size(); ii > 0; ii-- )
    vol_edges.get_and_step()->marked(0);
}
void GeomMeasureTool::find_close_loops(DLIList <RefVolume*> &ref_vols,
                                       DLIList <RefEdge*> &close_edges,
                                       DLIList <RefFace*> &close_loop_faces,
                                       DLIList <double> &small_lengths,
                                       double tol)
{
  int ii, jj;
  DLIList <RefFace*> ref_faces, temp_faces;
  RefVolume *ref_vol;
  RefFace *curr_face;
  for ( ii = 0; ii < ref_vols.size(); ii++ )
  {
    ref_vol = ref_vols.get_and_step();
    ref_vol->ref_faces(temp_faces);
      //uniquely add the faces.
    for ( jj = temp_faces.size(); jj > 0; jj-- )
    {
      curr_face = temp_faces.pop();
      if ( curr_face->marked()== 0 )
      {
        curr_face->marked(1);
        ref_faces.append(curr_face);
      }
    }
  }
  DLIList <RefEdge*> tmp_close_edges;
  DLIList <double> tmp_small_lengths;
  
  for ( ii = ref_faces.size(); ii > 0; ii-- )
  {
    curr_face = ref_faces.get_and_step();
    curr_face->marked(0);
    if ( curr_face->num_loops() < 2 )
      continue;
    tmp_close_edges.clean_out();
    tmp_small_lengths.clean_out();
    find_close_loops( curr_face, tmp_close_edges,
                      tmp_small_lengths, tol);
    if ( tmp_close_edges.size() > 0 )
    {
      assert(tmp_close_edges.size()%2 == 0 );
      tmp_close_edges.reset();
      tmp_small_lengths.reset();
      for ( jj = tmp_close_edges.size()/2; jj > 0; jj-- )
      {
        close_edges.append(tmp_close_edges.get_and_step());
        close_edges.append(tmp_close_edges.get_and_step());
        small_lengths.append(tmp_small_lengths.get_and_step());
        close_loop_faces.append(curr_face);
      }
    }
  }
  return;
}

  
void GeomMeasureTool::find_close_loops(RefFace *face,
                                       DLIList <RefEdge*> &close_edges,
                                       DLIList <double> &small_lengths,
                                       double tol)
{
    //only do faces with multiple loops since we are measuring the
    //distances between these loops.
  if ( face->num_loops() < 2 )
    return;
    //This is the most tricky of these functions.
    //To do this I'm going to facet the edges, then use an AbstractTree to
    //find the minimum distance between the loops and between a single
    //loop.
  PointLoopList boundary_point_loops;
  CubitStatus stat;
  stat = get_boundary_points(face, boundary_point_loops, GEOMETRY_RESABS*100);
  if ( stat != CUBIT_SUCCESS )
    return;
  
  SegLoopList boundary_seg_loops;
  stat = convert_to_lines( boundary_point_loops,
                           boundary_seg_loops,
                           face);
  if ( stat != CUBIT_SUCCESS )
    return;

    //Now add the points to the AbstractTree.
  DLIList<AbstractTree<GeomSeg*>*> atree_list;
  AbstractTree<GeomSeg*> *curr_tree;
  SegList *seg_list;
  GeomSeg  *curr_seg;
  int ii,jj, kk;
  //const double angle_convert = 180.0/CUBIT_PI;
  boundary_seg_loops.reset();
  int num_segs = 0;
  for (ii = 0; ii < boundary_seg_loops.size(); ii++ )
  {
    seg_list = boundary_seg_loops.get_and_step();
    curr_tree = new RTree<GeomSeg*>(GEOMETRY_RESABS);
    for ( jj = seg_list->size(); jj > 0; jj-- )
    {
        //build the r-tree.
      num_segs++;
      curr_seg = seg_list->get_and_step();
      curr_tree->add(curr_seg);
        //calculate the interior angles.
    }
    atree_list.append(curr_tree);
  }
    //determine the minimum distance between the loops.
  PointList *curr_points;
  GeomPoint *curr_point;
  DLIList <GeomSeg*> nearest_neighbors;
  RefEdge *closest_edge_1, *closest_edge_2;
  CubitVector curr_loc;
  double closest_dist;
  double min_for_loop_squared;
  atree_list.reset();
  for ( ii = 0; ii < atree_list.size(); ii++ )
  {
    curr_points = boundary_point_loops.get_and_step();
    for ( jj = ii+1; jj < atree_list.size(); jj++ )
    {
      min_for_loop_squared = CUBIT_DBL_MAX;
      closest_edge_1 = NULL;
      closest_edge_2 = NULL;
      curr_tree = atree_list.next(jj);
        //now for every point in curr_points, find it's closest_point
        //in the curr_tree.
      for ( kk = 0; kk < curr_points->size(); kk++ )
      {
        curr_point = curr_points->get_and_step();
        curr_loc.set(curr_point->coordinates());
        nearest_neighbors.clean_out();
        curr_tree->k_nearest_neighbor(curr_loc,
                                      1, closest_dist, nearest_neighbors,
                                      dist_sq_point_data);
        if ( nearest_neighbors.size() == 0)
        {
            //hmm, not sure what is wrong here.
          PRINT_ERROR("Can't find closest point between loops.\n");
        }
        if (closest_dist <= tol*tol && closest_dist < min_for_loop_squared )
        {
            //just store the smaller ones.  Don't store
          RefEntity *ref_ent = curr_point->owner();
          GeomSeg *near_seg = nearest_neighbors.get();
          RefEdge *ref_edge_1 = CAST_TO(ref_ent, RefEdge);
          if ( ref_edge_1 == NULL )
          {
              //assume this is a vertex.
            RefVertex *ref_vert = CAST_TO(ref_ent, RefVertex);
            if ( ref_vert == NULL )
            {
              PRINT_ERROR("problem with point ownership in closest loops.\n");
              continue;
            }
            DLIList <RefEdge*> ref_edges;
            ref_vert->ref_edges(ref_edges);
            int ll;
            RefEdge *the_edge = NULL;
            for ( ll = ref_edges.size(); ll > 0; ll-- )
            {
              if ( ref_edges.get()->is_directly_related(face) )
              {
                the_edge = ref_edges.get();
                break;
              }
              ref_edges.step();
            }
            if ( the_edge == NULL )
            {
              PRINT_ERROR("problem with point ownership in closest loops.\n");
              continue;
            }
            ref_edge_1 = the_edge;
          }
          ref_ent = near_seg->owner();
          RefEdge *ref_edge_2 = CAST_TO(ref_ent, RefEdge);
          if ( ref_edge_2 == NULL )
          {
            PRINT_ERROR("problem with point ownership in closest loops.\n");
            continue;
          }
          
          min_for_loop_squared = closest_dist;
          closest_edge_1 = ref_edge_1;
          closest_edge_2 = ref_edge_2;
        }
      }
      if ( closest_edge_1 != NULL )
      {
        close_edges.append(closest_edge_1);
        close_edges.append(closest_edge_2);
        small_lengths.append(sqrt(min_for_loop_squared));
      }
    }
  }
    //clean up the memory.
  int ll,mm;
  atree_list.reset();
  boundary_seg_loops.reset();
  for ( ll = boundary_seg_loops.size(); ll > 0; ll-- )
  {
    delete atree_list.pop();
    seg_list = boundary_seg_loops.pop();
    for ( mm = seg_list->size(); mm > 0; mm-- )
      delete seg_list->pop();
    delete seg_list;
  }
		
  PointList *point_list;
  for ( ll = boundary_point_loops.size(); ll > 0; ll-- )
  {
    point_list = boundary_point_loops.pop();
    for ( mm = point_list->size(); mm > 0; mm-- )
    {
      GeomPoint *tmp_point = point_list->pop();
        //delete the CubitPointData
      delete tmp_point;
    }
    delete point_list;
  }
  return;
}
double GeomMeasureTool::measure_area(RefFace* curr_face)
{
  double area;
  GeomDataObserver *g_d_o = GeomDataObserver::get(curr_face);
  if ( g_d_o == NULL )
    g_d_o = GeomDataObserver::create(curr_face);

  if ( g_d_o == NULL )
  {
    assert(g_d_o != NULL );
    return curr_face->measure();
  }
  if ( !g_d_o->is_measure_set() )
  {
    area  = (curr_face->get_surface_ptr())->measure();
    g_d_o->set_measure(area);
  }
  area = g_d_o->get_measure();
  return area;
}

///
///  This function simply gets the bad entities of the volume.  This
///  assumes the volume is from some solid modelar where this function
///  is defined. If not, it will just be empty...
///
void GeomMeasureTool::find_bad_geometry(RefVolume *volume,
                                        DLIList <RefEntity*> &bad_ents)
{
  Lump *lump_ptr = volume->get_lump_ptr();
  DLIList <TopologyEntity*> bad_top_ents;
  lump_ptr->validate(volume->entity_name(), bad_top_ents);
    //Now go through and get a RefEnttiy for each entity.
  RefEntity *ref_ent;
  TopologyEntity *topo_ent;
  GroupingEntity *gr_ent;
  SenseEntity *se_ent;
  
  int ii;
  for (ii = bad_top_ents.size(); ii > 0; ii-- )
  {
    topo_ent = bad_top_ents.pop();
    ref_ent = CAST_TO(topo_ent, RefEntity);
    if ( ref_ent != NULL )
    {
      bad_ents.append(ref_ent);
      continue;
    }
    gr_ent = CAST_TO(topo_ent, GroupingEntity);
    if ( gr_ent != NULL )
    {
      Chain *ch = CAST_TO(gr_ent, Chain);
      if ( ch != NULL )
      {
        RefVertex *start_v = ch->start_vertex();
        RefVertex *end_v = ch->end_vertex();
        RefEdge *ref_edge = start_v->common_ref_edge(end_v);
        bad_ents.append(ref_edge);
        continue;
      }
      Loop *lo = CAST_TO(gr_ent, Loop);
      if ( lo != NULL )
      {
        RefFace *rf = lo->get_ref_face_ptr();
        bad_ents.append(rf);
        continue;
      }
      Shell *sh = CAST_TO(gr_ent, Shell);
      if (sh != NULL )
      {
        RefVolume *rv = sh->get_ref_volume_ptr();
        bad_ents.append(rv);
        continue;
      }
      else
      {
          //this is incase there is some other entity we don't
          //know about, just get the body...
        Body *body = topo_ent->body();
        bad_ents.append(body);
        continue;
      }
    }
    se_ent = CAST_TO(topo_ent, SenseEntity);
    if ( se_ent != NULL )
    {
      BasicTopologyEntity *be = se_ent->get_basic_topology_entity_ptr();
      bad_ents.append(be);
      continue;
    }
      //If we are here we didn't get the right entity.  Just get
      //the body associated with this topology entity.
    Body *b = topo_ent->body();
    bad_ents.append(b);
  }
  return;
}

///
///  Find the irregular vertex valences.
///  Find things like vertices with valences greater than 4.
///  Assume for now that the volumes are not merged..
///
void GeomMeasureTool::find_irregular_valence( DLIList <RefVolume*> &ref_volumes,
                                              DLIList <RefVertex*> &irregular_vertices)
{
  RefVolume *ref_vol;
  int ii;
  for ( ii = 0; ii < ref_volumes.size(); ii++ )
  {
    ref_vol = ref_volumes.get_and_step();
    find_irregular_valence(ref_vol, irregular_vertices);
  }
}
///
///  Find the irregular vertex valences.
///  Find things like vertices with valences greater than 4.
///  Assume for now that the volumes are not merged..
///
void GeomMeasureTool::find_irregular_valence( RefVolume *ref_volume,
                                              DLIList <RefVertex*> &irregular_vertices)
{
  RefVertex *ref_vert;
  DLIList <RefVertex*> ref_vertices;
  ref_volume->ref_vertices(ref_vertices);
  int ii;
  for ( ii = ref_vertices.size(); ii > 0; ii-- )
  {
    ref_vert = ref_vertices.get_and_step();
    if ( ref_vert->num_ref_edges() > 4 )
      irregular_vertices.append(ref_vert);
  }
}

///
///  Find fillets and rounds.
///
#define FACE_BLEND      1
#define VERTEX_BLEND    2
#define BLEND_PROCESSED 3

void GeomMeasureTool::find_blends( RefVolume *ref_volume,
                                   DLIList <RefFace*> &blend_faces,
                                   DLIList <DLIList<RefFace*>*> &blend_groups )
{
    //Assume for now that blends have cartesian angles between
    //all the attached faces and itself.
    //Also a blend cannot be planar.  And directly accross from
    //one blend edge, there is usually another blended edge but
    //the angle between the two normals must be orthogonal. In other
    //words their must be some sort of transition.
  
    //mark all the faces as 0.  The mark is used here to find the
    // blend groups.  The function is_vertex_blend() marks faces
    // as being a vertex blend.  The group of faces is then traversed
    // until the adjacent surface (across the cross curve) is either
    // a vertex blend or not a blend.
  DLIList <RefFace*> ref_faces;
  DLIList <RefFace*> vertex_blends, face_blends;
  ref_volume->ref_faces(ref_faces);
  int ii, jj;
  for ( ii =0; ii < ref_faces.size(); ii++ )
    ref_faces.get_and_step()->marked(0);

    //First go through each face and each edge on each face.
  DLIList <RefEdge*> ref_edges;
  RefFace *ref_face;
  RefEdge *ref_edge, *other_edge;
  int num_vertex_blends = 0;
  int num_face_blends = 0;
  for ( ii = ref_faces.size(); ii > 0; ii-- )
  {
    ref_face = ref_faces.get_and_step();
      //Test the face to see if it is a face or vertex blend.
    other_edge = NULL;
    ref_edge = NULL;
    if ( is_face_blend(ref_face, ref_volume,
                       ref_edge, other_edge ) )
    {
      face_blends.append(ref_face);
      ref_face->marked(1);
      num_face_blends++;
    }
    else if ( is_vertex_blend(ref_face, ref_volume) )
    {
      vertex_blends.append(ref_face);
      ref_face->marked(2);
      num_vertex_blends++;  // keep track of the number of vertex blends
    }
  }
  if ( face_blends.size() + vertex_blends.size() == 0 )
  {
    return;
  }
    //Find out how many different groups of surfaces there
    //are that share curves.
  DLIList <RefFace*> *blend_group = NULL;
  DLIList <RefFace*> stack;
  RefFace *start_face = NULL, *other_face;
  
  // continue while there are blends to process
  while ( vertex_blends.size() + face_blends.size() > 0 || 
          start_face != NULL)
  {
    // if we don't have a start_face, get a new one.
    if (!start_face) 
    {
      // this is the start of a new group
      blend_group = new DLIList <RefFace*>;
      blend_groups.append(blend_group);

      // always prefer to start at a vertex blend if one exists
      if (vertex_blends.size() > 0 )
      {
        start_face = vertex_blends.pop();
      }
      else 
      {
        start_face = face_blends.pop();
      }
    }

    blend_group->append(start_face);  // add ref_face to group
    blend_faces.append(start_face);   // add the ref_face to the returned list
    start_face->marked(BLEND_PROCESSED);

    ref_edges.clean_out();
    start_face->ref_edges(ref_edges);
    for ( jj = 0; jj < ref_edges.size(); jj++ )
    {
      ref_edge = ref_edges.get_and_step();
      other_face = ref_edge->other_face(start_face, ref_volume);
      if (other_face == NULL)
      {
        start_face = NULL;
        break;
      }

      // if this blend has been processed and isn't in the current group
      if (other_face->marked() == BLEND_PROCESSED && 
          !blend_group->is_in_list(other_face) )
      {
        start_face = NULL;
        break;  // reached an end of this loop
      }
      else if (other_face->marked() == VERTEX_BLEND)
      {
        blend_group->append(other_face);  // add the ref_face to the group
        blend_faces.append(other_face);   // add the ref_face to the returned list
        vertex_blends.remove(other_face);
        other_face->marked(BLEND_PROCESSED);
        start_face = NULL;
        break;  // a vertex blend is also end of the line
      }
      else if (other_face->marked() == FACE_BLEND)
      {
        start_face = other_face;
        face_blends.remove(other_face);
        break;  // continue using this face as the new start face
      }
    }
    // if we traversed through all the edges of this blend without
    // finding another blend attached, this is the end of the chain 
    // and we need to find another starting place.
    if ( jj >= ref_edges.size() ) 
    {
      start_face = NULL;
    }
    
  }

    //cleanup marks.
  for ( ii =0; ii < ref_faces.size(); ii++ )
    ref_faces.get_and_step()->marked(0);
  for ( ii =0; ii < blend_faces.size(); ii++ )
    blend_faces.get_and_step()->marked(0);
    //Okay we should have it now.

    //REMEMBER TO DELETE the dllists in blend_faces.
  return;
}

// struct type object private to the next function
class NextBlend
{
  public:
  RefEdge* ref_edge;
  RefFace* ref_face;
};

void GeomMeasureTool::find_blends_from_edge( RefVolume* ref_volume, 
                                                    RefFace *start_face, 
                                                    RefEdge* start_edge,
                                                    std::vector <RefFace*> &blend_faces)
{
  // initialize the mark on all surfaces in the volume
  DLIList <RefFace*> ref_faces;
  RefEdge *next_edge, *other_edge;
  ref_volume->ref_faces(ref_faces);
  int ii;
  for ( ii =0; ii < ref_faces.size(); ii++ )
    ref_faces.get_and_step()->marked(0);

  std::stack<NextBlend> blend_stack;

  NextBlend blend;
  blend.ref_edge = start_edge;
  blend.ref_face = start_face;
  blend_faces.push_back(start_face);

  blend_stack.push(blend);

  while (blend_stack.size() > 0)
  {
    // this is an oddity with std::stack.  Get the top of the stack
    // and then you have to remove it from the stack with pop
    NextBlend next_blend = blend_stack.top();
    blend_stack.pop();
    RefEdge* ref_edge = next_blend.ref_edge;
    RefFace* ref_face = next_blend.ref_face;

    RefFace* next_face;
    next_face = ref_edge->other_face(ref_face, ref_volume);

    // the the next face exists and it's a face blend, save the current blend,
    // create a new blend object and add it to the stack
    if ( next_face && is_face_blend(next_face, ref_volume, next_edge, other_edge ) )
    {
      // if the next face is already processed, the chain loops back on itself.
      if (next_face->marked() == BLEND_PROCESSED)
      {
        for ( ii =0; ii < ref_faces.size(); ii++ )
          ref_faces.get_and_step()->marked(0);
        return;
      }

      blend_faces.push_back(next_face);
      next_face->marked(BLEND_PROCESSED);

      // the is_face_blend function assumes (poorly) rectangular surfaces
      // without any holes.  It returns the spring curves and we want
      // the cross curves. So get all four edges and remove the two
      // spring curves returned by is_face_blend
      DLIList <RefEdge*> ref_edges;
      next_face->ref_edges(ref_edges);
      ref_edges.remove(next_edge);
      ref_edges.remove(other_edge);

      if (ref_edges.get() != ref_edge)
      {
        next_blend.ref_edge = ref_edges.get();
      }
      else
      {
        next_blend.ref_edge = ref_edges.next(1);
      }
      next_blend.ref_face = next_face;
      blend_stack.push(next_blend);
    }
    else if ( next_face && is_vertex_blend(next_face, ref_volume) )
    {
      // we will stop the chain at a vertex blend
      blend_faces.push_back(next_face);
      next_face->marked(BLEND_PROCESSED);
    }
  }

  // clean up the marks
  std::vector<RefFace*>::iterator iter;
  for ( iter = blend_faces.begin(); iter != blend_faces.end(); iter++)
    (*iter)->marked(0);
}

// given a starting blend surface find a chain of blends from
// that surface.  
//
// Note that this function intentionally does _not_
// clear the blend_face list so that additional chains can be added.
CubitStatus GeomMeasureTool::find_blend_chains( RefFace *start_face,
                                std::vector<std::vector< RefFace*> > &blend_chains)
{

  if (start_face == NULL)
  {
    return CUBIT_FAILURE;
  }

  std::vector <RefFace*> blend_faces;

  // get the owning volume of this blend
  DLIList<RefEntity*> entity_list;
  RefVolume* ref_volume;
  int ii;

  start_face->get_parent_ref_entities(entity_list);

  // this indicates merged enitites and potential problems
  if (entity_list.size() > 1)
  {
    return CUBIT_FAILURE;
  }

  // make sure we're at the beginning of the list and get the first
  // and only item on the list and cast it to a RefVolume
  entity_list.reset();
  ref_volume = CAST_TO(entity_list.get(), RefVolume);
 
  if (!ref_volume)
  {
    return CUBIT_FAILURE;
  }

  // initialize the mark on all surfaces in the volume
  DLIList <RefFace*> ref_faces;
  ref_volume->ref_faces(ref_faces);

  RefEdge *spring_curve1, *spring_curve2;
  if ( is_face_blend(start_face, ref_volume, spring_curve1, spring_curve2 ) )
  {
    // the is_face_blend function assumes (poorly) rectangular surfaces
    // without any holes.  It returns the spring curves and we want
    // the cross curves. So get all four edges and remove the two
    // spring curves returned by is_face_blend
    DLIList <RefEdge*> ref_edges;
    start_face->ref_edges(ref_edges);
    ref_edges.remove(spring_curve1);
    ref_edges.remove(spring_curve2);

    // there is a special case where the blend is a periodic surface
    // meaning that there _are_ no cross curves
    if ( ref_edges.size() == 0 )
    {
      blend_faces.push_back(start_face);
    }
    else
    {
      // now find additional blends extending from either side of the blend
      blend_faces.clear();
      find_blends_from_edge( ref_volume, start_face, ref_edges.get(), blend_faces);
      find_blends_from_edge( ref_volume, start_face, ref_edges.next(1), blend_faces);

      // make sure that we have a unique list (the start surface is probably here twice)
      std::vector<RefFace*>::iterator new_end;
      std::sort( blend_faces.begin(), blend_faces.end() );
      new_end = std::unique( blend_faces.begin(), blend_faces.end() );
      blend_faces.erase(new_end, blend_faces.end());
    }

    blend_chains.push_back(blend_faces);
  }
  else if ( is_vertex_blend(start_face, ref_volume) )
  {
    DLIList<RefEdge*> ref_edges;
    start_face->ref_edges(ref_edges);

    for (ii = 0; ii < ref_edges.size(); ii++)
    {
      RefEdge* start_edge = ref_edges.get_and_step();
      blend_faces.clear();
      find_blends_from_edge( ref_volume, start_face, start_edge, blend_faces);

      blend_chains.push_back(blend_faces);
    }
  }
  else
  {
    // the given face is not a blend
    return CUBIT_FAILURE;
  }

  return CUBIT_SUCCESS;
}

//--------------------------------------------------------------------
//Function: Public, is_face_blend
//Description: Determines if a face is a blend surface, returns the
//   ref_edge on one side of the blend and other_edge on the opposite side.
//   For this type of blend, only ref_edge must be tangentially meeting
//   with another surface.  Other_edge must be oriented orthogonally to
//   it and may or may not blend with another surface.  This assumes
//   a rectangular blend surface, without holes.
//---------------------------------------------------------------------
CubitBoolean GeomMeasureTool::is_face_blend(RefFace *ref_face,
                                            RefVolume *ref_volume,
                                            RefEdge *&ref_edge,
                                            RefEdge *&other_edge)
{
    //first we know that blend surfaces are not planar.
    //so remove these first.
    //Also, don't look at faces with more than 2 loops.
  if ( ref_face->geometry_type() == PLANE_SURFACE_TYPE ||
       ref_face->num_loops() > 2 )
    return CUBIT_FALSE;
  
  CubitBoolean is_cartesian;
  DLIList<RefEdge*> ref_edges;
  ref_edge = NULL;
  other_edge = NULL;
  RefFace *other_face;
  int jj;
  double angle;
  ref_face->ref_edges(ref_edges);
  for ( jj = ref_edges.size(); jj > 0; jj-- )
  {
    ref_edge = ref_edges.get_and_step();

    //Weed-out case where one edge is shared between more 
    //than 2 surfaces of the same volume
    DLIList<RefFace*> tmp_faces;
    ref_edge->ref_faces( tmp_faces );

    if( tmp_faces.size() > 2 )
    {
      int kk;
      for(kk=tmp_faces.size(); kk--;)
      {
        if( !tmp_faces.get()->is_child( ref_volume ) )
          tmp_faces.change_to(NULL);
        tmp_faces.step();
      }
      tmp_faces.remove_all_with_value( NULL );
      if( tmp_faces.size() > 2 )
        //this isn't the type of surface we are looking for...
        continue;
    }

    other_face = ref_edge->other_face(ref_face, ref_volume);
    if ( other_face == NULL )
    {
        //this isn't the type of surface we are looking for...
      break;
    }
    angle = GeometryQueryTool::instance()->surface_angle(ref_face,
                                                         other_face,
                                                         ref_edge,
                                                         ref_volume);
    angle *= 180.0/CUBIT_PI;
    is_cartesian = CUBIT_TRUE;
    if ( angle <= GEOM_SIDE_LOWER ||
         angle >= GEOM_SIDE_UPPER )
      is_cartesian = CUBIT_FALSE;
      //Okay, we have one major criteria achieved, this edge is a cartesian meet.

    if ( !is_cartesian )
      continue;
      //Now we need to check the radius of curvatures between these
      //two surfaces. I'm not totally sure about this but I think we
      //don't want the same radius of curvature.
    double k1_s1, k2_s1, k1_s2, k2_s2;
    CubitVector mid_point;
    ref_edge->mid_point(mid_point);
    ref_face->get_principal_curvatures( mid_point, k1_s1, k2_s1, ref_volume);
    other_face->get_principal_curvatures( mid_point, k1_s2, k2_s2, ref_volume);
    if (( is_equal(k1_s1, k1_s2) || is_equal(k1_s1, k2_s2) ) &&
        ( is_equal(k2_s1, k1_s2) || is_equal(k2_s1, k2_s2) ) )
        //try a different edge.
      continue;

      //Okay, this looks pretty good.
      //Now try and find an edge that is on the "opposite" side of
      // this surface. Sort of assume this is like a mapable surface.
      //check the oposite side.
    other_edge = NULL;
    CubitVector opp_point;
    if ( !find_opposite_edge( ref_edge, ref_face, other_edge, opp_point) ||
         other_edge == NULL )
        //if we can't find an opposite edge, we need to not
        //continue checking this surface.
      return CUBIT_FALSE;
    
      //Okay, we have the point opposite.  If this is a blend surface
      //the normal of this point should be pointing 90 degrees from
      //us.
    CubitVector normal_1 = ref_face->normal_at(mid_point, ref_volume);
    CubitVector normal_2 = ref_face->normal_at(opp_point, ref_volume);
    normal_1.normalize();
    normal_2.normalize();
    double dot_product = normal_1 % normal_2;
    if ( dot_product >= .5 || dot_product <= -.5  )
        //try a different edge.
      continue;
    CubitVector tangent_1, tangent_2;
    ref_edge->tangent(mid_point, tangent_1, ref_face);
    other_edge->tangent(opp_point, tangent_2, ref_face);
    tangent_1.normalize();      
    tangent_2.normalize();      
    dot_product = tangent_1 % tangent_2;
      //basically these edges should be somewhat parallel at
      //this point.
    if ( dot_product < .5 && dot_product > -.5 )
      continue;
//This was called but never checked, not sure why.
//      //check the curvature opposite.
//    double k1_s3, k2_s3;
//    ref_face->get_principal_curvatures( opp_point, k1_s3, k2_s3, ref_volume);
      //we are done with this face. It is a blend face.
    return CUBIT_TRUE;
  }
  return CUBIT_FALSE;
}

//--------------------------------------------------------------------
//Function: Public, is_vertex_blend
//Description: Determines if a face is a vertex blend surface.
//   For this type of blend, all ref_edges must be meet tangentially
//   with another surface.  This assumes blend surface with no holes. 
//---------------------------------------------------------------------
CubitBoolean GeomMeasureTool::is_vertex_blend(RefFace *ref_face,
                                              RefVolume* ref_volume)
{
    //first we know that blend surfaces are not planar.
    //so remove these first.
    //Also, don't look at faces with more than 2 loops.
  if ( ref_face->geometry_type() == PLANE_SURFACE_TYPE ||
       ref_face->num_loops() > 2 )
    return CUBIT_FALSE;
  
  CubitBoolean is_cartesian;
  DLIList<RefEdge*> ref_edges;
  RefFace *other_face;
  RefEdge *ref_edge = NULL;
  int jj;
  double angle;
  ref_face->ref_edges(ref_edges);
  for ( jj = ref_edges.size(); jj > 0; jj-- )
  {
    ref_edge = ref_edges.get_and_step();

    //Weed-out case where one edge is shared between more 
    //than 2 surfaces of the same volume
    DLIList<RefFace*> tmp_faces;
    ref_edge->ref_faces( tmp_faces );

    if( tmp_faces.size() > 2 )
    {
      int kk;
      for(kk=tmp_faces.size(); kk--;)
      {
        if( !tmp_faces.get()->is_child( ref_volume ) )
          tmp_faces.change_to(NULL);
        tmp_faces.step();
      }
      tmp_faces.remove_all_with_value( NULL );
      if( tmp_faces.size() > 2 )
        //this isn't the type of surface we are looking for...
        continue;
    }

    other_face = ref_edge->other_face(ref_face, ref_volume);
    if ( other_face == NULL )
    {
        //this isn't the type of surface we are looking for...
      break;
    }
    angle = GeometryQueryTool::instance()->surface_angle(ref_face,
                                                         other_face,
                                                         ref_edge,
                                                         ref_volume);
    angle *= 180.0/CUBIT_PI;
    is_cartesian = CUBIT_TRUE;
    if ( angle <= GEOM_SIDE_LOWER ||
         angle >= GEOM_SIDE_UPPER )
      is_cartesian = CUBIT_FALSE;
      //Okay, we have one major criteria achieved, this edge is a cartesian meet.

    if ( !is_cartesian )
      return CUBIT_FALSE;
      //Now we need to check the radius of curvatures between these
      //two surfaces. I'm not totally sure about this but I think we
      // want the same radius of curvature.
    double k1_s1, k2_s1, k1_s2, k2_s2;
    CubitVector mid_point;
    ref_edge->mid_point(mid_point);
    ref_face->get_principal_curvatures( mid_point, k1_s1, k2_s1, ref_volume);
    other_face->get_principal_curvatures( mid_point, k1_s2, k2_s2, ref_volume);
    if (( is_equal(k1_s1, k1_s2) || is_equal(k1_s1, k2_s2) ) &&
        ( is_equal(k2_s1, k1_s2) || is_equal(k2_s1, k2_s2) ) )
        //try a different edge.
        continue;
    else
      return CUBIT_FALSE;
  }

  // if all edges are tangent and share curvatures then it must be a 
  // vertex blend.
  return CUBIT_TRUE;
}

CubitBoolean GeomMeasureTool::find_opposite_edge( RefEdge* ref_edge,
                                                  RefFace* ref_face,
                                                  RefEdge *&other_edge,
                                                  CubitVector &closest_point)
{
    //Here we want to find the edge that is "opposite" to ref_edge.
    //Opposite is defined by assuming the surface is a square or cylinder and
    //finding the edge on the other side, so that there is a point closest
    //to the mid point on the ref_edge.
  //double need_to_turn = CUBIT_PI;
  int ii;
  DLIList <Loop*> loops;
  DLIList <CoEdge*> co_edges;
  Loop *curr_loop;
  CoEdge *co_edge, *this_co_edge;
  ref_edge->get_co_edges(co_edges, ref_face);
  if ( co_edges.size() != 1 )
    return CUBIT_FALSE;
  this_co_edge = co_edges.get();
  co_edges.clean_out();
  ref_face->loops(loops);
  if ( loops.size() > 2 )
    return CUBIT_FALSE;
  CubitBoolean found_loop = CUBIT_FALSE;
    //find the loop where this_co_edge is on.
  for(ii = loops.size(); ii > 0; ii-- )
  {
    curr_loop = loops.get_and_step();
    co_edges.clean_out();
    curr_loop->ordered_co_edges(co_edges);
    if ( co_edges.move_to(this_co_edge) )
    {
      found_loop = CUBIT_TRUE;
      break;
    }
  }
  if ( !found_loop )
    return CUBIT_FALSE;
  CoEdge *opp_co_edge = NULL;
  if ( loops.size() == 1 )
  {
      //This is the normal case, I hope.
      //this co_edge, go around the interior until
      //we turn twice.  If there are angles greater than CUBIT_PI, this isn't
      //a traditional blend.  Everything should be convex.
    if ( co_edges.size() == 4 )
    {
        //Just get the one opposite...
      opp_co_edge = co_edges.next(2);
    }
    else
    {
      //double total_angle = 0.0;
      int turn_counter = 0;
      CoEdge *next_co_edge = NULL;
      for ( ii = co_edges.size(); ii > 0; ii-- )
      {
        co_edge = co_edges.get_and_step();
        next_co_edge = co_edges.get();
        double angle = GeometryQueryTool::instance()->geometric_angle(co_edge, next_co_edge);
        angle *= 180.0/CUBIT_PI;
          //we don't deal with small or concave angles.
        if ( angle < 10.0 )
          return CUBIT_FALSE;
        else if ( angle > 185.0 )
          return CUBIT_FALSE;
        else if ( angle < 140.0 )
          turn_counter++;
        if ( turn_counter == 2 )
        {
          opp_co_edge = next_co_edge;
          break;
        }
      }
    }
  }
  else
  {
      //Just get the other loop and get one coedge in it.
    co_edges.clean_out();
    curr_loop = loops.get();
    curr_loop->ordered_co_edges(co_edges);
    opp_co_edge = co_edges.get();
  }
  if ( opp_co_edge == NULL )
    return CUBIT_FALSE;
    //This is just here to test correctness to this point, it can be taken out
    //if need be.
  if ( !co_edges.move_to(opp_co_edge) )
  {
    PRINT_ERROR("logic messed-up in find_opposite edge\n");
    assert(0);
    return CUBIT_FALSE;
  }
    //okay, we now have the co_edges in the correct position and
    //we have a canidate co_edge.
    //Search from the co_edge with the point that has closest point
    //to the mid point of this_co_edge.
  other_edge = NULL;
  CubitVector mid_point;
  ref_edge->mid_point(mid_point);
  CubitVector tmp_closest_point;
  double min_dist_sq = CUBIT_DBL_MAX;
  double dist_sq;
  RefEdge *test_ref_edge;
  CoEdge *prev_co_edge;
  CubitBoolean first = CUBIT_TRUE;
  for ( ii = co_edges.size(); ii > 0; ii-- )
  {
    co_edge = co_edges.get_and_step();
    prev_co_edge = co_edges.prev(2);
    test_ref_edge = co_edge->get_ref_edge_ptr();
      //We don't want edges that are connected to ref edge.
      //If we have them, we have gone to far around the surface.
    if ( test_ref_edge->common_ref_vertex(ref_edge) )
      break;
    if ( !first )
    {
        //break out of the for loop if we turn again. We just want the closest
        //point opposite.
      double angle = GeometryQueryTool::instance()->
        geometric_angle(prev_co_edge, co_edge);
      angle *= 180.0/CUBIT_PI;
      if ( angle > 195.0 || angle < 145.0 )
        break;
    }
    else
      first = CUBIT_FALSE;
      //now find the closest point to the mid point of ref_edge.
    test_ref_edge->closest_point_trimmed( mid_point,
                                          tmp_closest_point);
    dist_sq = (mid_point - closest_point).length_squared();
      //find and keep the closest one.
    if ( dist_sq < min_dist_sq )
    {
      other_edge = test_ref_edge;
      min_dist_sq = dist_sq;
      closest_point = tmp_closest_point;
    }
  }
  if ( other_edge == NULL )
    return CUBIT_FALSE;

  return CUBIT_TRUE;
}
CubitBoolean GeomMeasureTool::is_equal(double v1, double v2)
{
  double smallv = v1-v2;
  if ( smallv < 0.0 )
  {
    if ( smallv >= -CUBIT_RESABS )
      return CUBIT_TRUE;
    else
      return CUBIT_FALSE;
  }
  else
  {
    if ( smallv <= CUBIT_RESABS )
      return CUBIT_TRUE;
    else
      return CUBIT_FALSE;
  }
}
CubitStatus GeomMeasureTool::get_centroid( RefFace *ref_face, CubitVector &centroid, double &tot_area )
{
  GMem g_mem;
  unsigned short norm_tol = 5;
  double dist_tol = -1.0;

  ref_face->get_geometry_query_engine()->
      get_graphics(ref_face->get_surface_ptr(), &g_mem, norm_tol, dist_tol );

  if(g_mem.fListCount < 1)
  {
      // Decrease tolerance and try again (we can get this for small features)
      norm_tol /= 2;
      ref_face->get_geometry_query_engine()->
          get_graphics(ref_face->get_surface_ptr(), &g_mem, norm_tol, dist_tol );
  }

  if(g_mem.fListCount < 1)
  {
      // Lets give up 
      PRINT_ERROR( "Unable to find the center of a surface\n" );
      return CUBIT_FAILURE;
  }

  // Initialize
  tot_area = 0.0;
  centroid.set( 0.0, 0.0, 0.0 );

  // Loop through the triangles
  double tri_area;
  GPoint p[3];
  GPoint* plist = g_mem.point_list();
  int* facet_list = g_mem.facet_list();
  int c = 0;
  for( ; c < g_mem.fListCount; )
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

      //Calculate area
      CubitVector vec1 = (p2 - p1);
	  CubitVector vec2 = (p3 - p1);

	  CubitVector cross = vec1 * vec2;
		
	  tri_area = 0.5* sqrt(cross.x() * cross.x() + cross.y() * cross.y() + cross.z() * cross.z());
		
		
	  centroid += (tri_area * center);

	  tot_area += tri_area;
        
    }
  if( tot_area == 0 )
    return CUBIT_FAILURE;

  centroid /= tot_area;
  return CUBIT_SUCCESS;
}    

CubitStatus
GeomMeasureTool::center( DLIList<RefFace*> ref_faces )
{
  int ii,id(0);
  double surf_area;
  double tot_area = 0.0;
  CubitVector surf_centroid;
  CubitVector centroid;
  centroid.set( 0.0, 0.0, 0.0 );

  for ( ii = ref_faces.size(); ii > 0; ii-- )
  {
	RefFace *ref_face = ref_faces.get_and_step();
	if( GeomMeasureTool::get_centroid( ref_face, surf_centroid, surf_area) != CUBIT_SUCCESS )
	  return CUBIT_FAILURE;
    centroid += (surf_area * surf_centroid);
	tot_area += surf_area;
	id = ref_face->id();
  }
  
  centroid /= tot_area;

  if( ref_faces.size() > 1 )
	PRINT_INFO("Centroid of Composite surface located at: ( %f , %f , %f )\n\n",
                 centroid.x(), centroid.y(), centroid.z());
  else
	PRINT_INFO("Centroid of Surface %i located at: ( %f , %f , %f )\n\n",
                  id, centroid.x(), centroid.y(), centroid.z());
  return CUBIT_SUCCESS;
}

CubitStatus GeomMeasureTool::find_near_coincident_vertices( 
                            DLIList<RefVolume*> &ref_volumes,
                            DLIList<RefVertex*> &ref_vertices_out,
                            DLIList<double> &distances,
                            double low_tol,
                            double high_tol,
                            bool filter_same_volume_cases)
{
  DLIList<RefVertex*> tmp_vert_list;
  DLIList<RefVertex*> ref_verts; 
  int i,j;
  for( i=ref_volumes.size(); i--; )
  {
    RefVolume *tmp_vol = ref_volumes.get_and_step();
    tmp_vert_list.clean_out();
    tmp_vol->ref_vertices( tmp_vert_list );
    ref_verts += tmp_vert_list;
  }

  //put all the vertices in a tree 
  AbstractTree <RefVertex*> *a_tree = new RTree<RefVertex*>( high_tol ); 
  for (i=ref_verts.size(); i--;)
    a_tree->add(ref_verts.get_and_step());

  std::multimap<double, dist_vert_struct> distance_vertex_map; 

  //for each vertex
  for (i=ref_verts.size(); i--;)
  {
    RefVertex *tmp_vert = ref_verts.get_and_step();
    RefVolume *v1 = tmp_vert->ref_volume();
    CubitVector vert_xyz = tmp_vert->coordinates();

    //get all close vertices
    DLIList<RefVertex*> close_verts;
    a_tree->find(tmp_vert->bounding_box(), close_verts);

    //if any vertex is between low_tol and high_tol
    //add it to the list
    DLIList<RefVertex*> near_coincident_verts;
    for( j=close_verts.size(); j--; )
    {
      RefVertex *close_vert = close_verts.get_and_step();
      if( close_vert == tmp_vert ) 
        continue;

      RefVolume *v2 = close_vert->ref_volume();
      bool check_distance = true;
      if(filter_same_volume_cases && v1 && v2 && v1 == v2)
        check_distance = false;
      if(check_distance)
      {
        double distance = vert_xyz.distance_between( close_vert->coordinates() );
        if( distance >= low_tol && distance <= high_tol )
        {
          dist_vert_struct tmp_struct;
          tmp_struct.dist = distance;
          tmp_struct.v1 = tmp_vert;
          tmp_struct.v2 = close_vert;
          distance_vertex_map.insert( std::multimap<double, dist_vert_struct>::
                                      value_type( distance, tmp_struct ));
        }
      }
    }

    a_tree->remove( tmp_vert );
  }

  std::multimap<double, dist_vert_struct>::reverse_iterator iter;
  
  iter = distance_vertex_map.rbegin();
  for(; iter!=distance_vertex_map.rend(); iter++ )
  {
    distances.append( (*iter).second.dist );
    ref_vertices_out.append( (*iter).second.v1 );
    ref_vertices_out.append( (*iter).second.v2 );
  }

  delete a_tree;

  return CUBIT_SUCCESS;
}

// This function is similar to find_near_coincident_vertices except for the
// fact that it will only find the closest vertex in a given volume for
// a vertex in another volume to be close to.  This tries to exclude the case where
// you would attempt to merge one vertex from one volume to two different
// vertices in another volume.
CubitStatus GeomMeasureTool::find_near_coincident_vertices_unique( 
                            DLIList<RefVolume*> &ref_volumes,
                            double high_tol,
                            std::map <RefVertex*, DLIList<dist_vert_struct*>*> &vert_dist_map)
{
  DLIList<RefVertex*> tmp_vert_list;
  DLIList<RefVertex*> ref_verts; 
  int i,j;
  for( i=ref_volumes.size(); i--; )
  {
    RefVolume *tmp_vol = ref_volumes.get_and_step();
    tmp_vert_list.clean_out();
    tmp_vol->ref_vertices( tmp_vert_list );
    ref_verts += tmp_vert_list;
  }

  //put all the vertices in a tree 
  AbstractTree <RefVertex*> *a_tree = new RTree<RefVertex*>( high_tol ); 
  for (i=ref_verts.size(); i--;)
    a_tree->add(ref_verts.get_and_step());

  //for each vertex
  for (i=ref_verts.size(); i--;)
  {
    RefVertex *tmp_vert = ref_verts.get_and_step();
    RefVolume *vol1 = tmp_vert->ref_volume();
    CubitVector vert_xyz = tmp_vert->coordinates();

    //get all close vertices
    DLIList<RefVertex*> close_verts;
    a_tree->find(tmp_vert->bounding_box(), close_verts);

    //if any vertex is between low_tol and high_tol
    //add it to the list
    DLIList<dist_vert_struct*> *near_coincident_verts = NULL;
    for( j=close_verts.size(); j--; )
    {
      RefVertex *close_vert = close_verts.get_and_step();
      if( close_vert == tmp_vert ) 
        continue;

      RefVolume *vol2 = close_vert->ref_volume();
      if(vol1 && vol2 && vol1 != vol2)
      {
        if(!near_coincident_verts)
        {
          near_coincident_verts = new DLIList<dist_vert_struct*>;
          vert_dist_map[tmp_vert] = near_coincident_verts;
        }
        double distance = vert_xyz.distance_between( close_vert->coordinates() );
        int h;
        bool found_entry_with_same_vol = false;
        for(h=near_coincident_verts->size(); h>0 && !found_entry_with_same_vol; h--)
        {
          dist_vert_struct* vds = near_coincident_verts->get_and_step();
          if(vds->vol2 == vol2)
          {
            found_entry_with_same_vol = true;
            if(distance < vds->dist)
            {
              vds->dist = distance;
              vds->vol2 = vol2;
              vds->v2 = close_vert;
            }
          }
        }
        if(!found_entry_with_same_vol)
        {
          dist_vert_struct *new_vds = new dist_vert_struct;
          new_vds->dist = distance;
          new_vds->v2 = close_vert;
          new_vds->vol2 = vol2;
          near_coincident_verts->append(new_vds);
        }
      }
    }
    a_tree->remove( tmp_vert );
  }

  delete a_tree;

  return CUBIT_SUCCESS;
}

struct dist_vert_vert_struct
{
  double dist;
  RefVertex *vert1;
  RefVertex *vert2;
};

struct dist_vert_curve_struct
{
  double dist;
  RefVertex *vert;
  RefEdge *edge;
//  bool operator<( const dist_vert_curve_struct& b ) const
//  {
//    return this->dist < b.dist;
 // }
};

struct vert_curve_dist_sort
{
  bool operator()( const dist_vert_curve_struct& a, const dist_vert_curve_struct& b ) const
  {
    return a.dist < b.dist;
  }
};

struct vert_curve_dist_sort_ptr
{
  bool operator()( dist_vert_curve_struct *a, dist_vert_curve_struct *b ) const
  {
    if( a->dist < b->dist )
      return true;
    else if( a->dist > b->dist )
      return false;
    else 
      return true;
  }
};

struct vert_vert_dist_sort_ptr
{
  bool operator()( dist_vert_vert_struct *a, dist_vert_vert_struct *b ) const
  {
    if( a->dist < b->dist )
      return true;
    else if( a->dist > b->dist )
      return false;
    else 
      return true;
  }
};

CubitStatus GeomMeasureTool::find_closest_vertex_curve_pairs(
                                  DLIList<RefVolume*> &vol_list,
                                  int &num_to_return,
                                  DLIList<RefVertex*> &vert_list,
                                  DLIList<RefEdge*> &curve_list,
                                  DLIList<double> &distances)
{
  DLIList<RefFace*> surfs;

  int i, total_num_entries = 0;
  for( i=vol_list.size(); i--; )
  {
    RefVolume *tmp_vol = vol_list.get_and_step();
    tmp_vol->ref_faces( surfs );
  }

  std::set<dist_vert_curve_struct*,vert_curve_dist_sort_ptr> distance_vertex_curve_set; 

  for(i=surfs.size(); i>0; i--)
  {
    RefFace *surf = surfs.get_and_step();
    DLIList<RefVertex*> surf_verts;
    surf->ref_vertices(surf_verts);
    DLIList<RefEdge*> surf_curves;
    surf->ref_edges(surf_curves);

    int j;
    for(j=surf_verts.size(); j>0; j--)
    {
      RefVertex *tmp_vert = surf_verts.get_and_step();
      CubitVector vert_xyz = tmp_vert->coordinates();
      CubitVector closest_pt;
      int k;
      for(k=surf_curves.size(); k>0; k--)
      {
        RefEdge *cur_edge = surf_curves.get_and_step();
        if(cur_edge->start_vertex() != tmp_vert &&
          cur_edge->end_vertex() != tmp_vert)
        {
          cur_edge->closest_point_trimmed(vert_xyz, closest_pt);
          if(!closest_pt.about_equal(cur_edge->start_coordinates()) &&
             !closest_pt.about_equal(cur_edge->end_coordinates()))
          {
            double dist_sq = vert_xyz.distance_between_squared(closest_pt);
            dist_vert_curve_struct *tmp_struct = new dist_vert_curve_struct;
            tmp_struct->dist = dist_sq;
            tmp_struct->vert = tmp_vert;
            tmp_struct->edge = cur_edge; 
            distance_vertex_curve_set.insert( tmp_struct );
            total_num_entries++;
          }
        }
      }
    }
  }

  std::set<dist_vert_curve_struct*, vert_curve_dist_sort_ptr>::iterator iter, upper_iter; 
  
  int local_num_to_return = num_to_return;
  if(num_to_return == -1)
  {
    local_num_to_return = total_num_entries;
  }
  int cntr = 0;
  iter = distance_vertex_curve_set.begin();
  for(; iter!=distance_vertex_curve_set.end() && cntr < local_num_to_return; iter++ )
  {
    distances.append( sqrt((*iter)->dist) );
    vert_list.append( (*iter)->vert );
    curve_list.append( (*iter)->edge );
    cntr++;
  }

  iter = distance_vertex_curve_set.begin();
  for(; iter!=distance_vertex_curve_set.end(); iter++ )
  {
    delete *iter;
  }

  return CUBIT_SUCCESS;
}

CubitStatus GeomMeasureTool::find_closest_vertex_vertex_pairs(
                                  DLIList<RefVolume*> &vol_list,
                                  int &num_to_return,
                                  DLIList<RefVertex*> &vert_list1,
                                  DLIList<RefVertex*> &vert_list2,
                                  DLIList<double> &distances)
{
  std::set<dist_vert_vert_struct*,vert_vert_dist_sort_ptr> distance_vertex_vertex_set; 

  int i, total_num_entries = 0;
  for( i=vol_list.size(); i--; )
  {
    RefVolume *tmp_vol = vol_list.get_and_step();
    DLIList<RefVertex*> vol_verts;
    tmp_vol->ref_vertices(vol_verts);
    while(vol_verts.size() > 1)
    {
      RefVertex *vert1 = vol_verts.pop();
      CubitVector vert1_xyz = vert1->coordinates();
      int j;
      for(j=vol_verts.size(); j>0; j--)
      {
        RefVertex *vert2 = vol_verts.get_and_step();
        double dist_sq = vert2->coordinates().distance_between_squared(vert1_xyz);
        dist_vert_vert_struct *tmp_struct = new dist_vert_vert_struct;
        tmp_struct->dist = dist_sq;
        tmp_struct->vert1 = vert1;
        tmp_struct->vert2 = vert2; 
        distance_vertex_vertex_set.insert( tmp_struct );
        total_num_entries++;
      }
    }
  }

  std::set<dist_vert_vert_struct*, vert_vert_dist_sort_ptr>::iterator iter, upper_iter; 
  
  int local_num_to_return = num_to_return;
  if(num_to_return == -1)
    local_num_to_return = total_num_entries;
  int cntr = 0;
  iter = distance_vertex_vertex_set.begin();
  for(; iter!=distance_vertex_vertex_set.end() && cntr < local_num_to_return; iter++ )
  {
    distances.append( sqrt((*iter)->dist) );
    vert_list1.append( (*iter)->vert1 );
    vert_list2.append( (*iter)->vert2 );
    cntr++;
  }

  iter = distance_vertex_vertex_set.begin();
  for(; iter!=distance_vertex_vertex_set.end(); iter++ )
  {
    delete *iter;
  }

  return CUBIT_SUCCESS;
}

CubitStatus GeomMeasureTool::find_near_coincident_vertex_curve_pairs( 
                                DLIList<RefVolume*> &ref_vols,
                                DLIList<RefEdge*> &ref_edges,
                                DLIList<RefVertex*> &ref_verts,
                                DLIList<double> &distances,
                                double low_tol,
                                double high_tol,
                                bool filter_same_volume_cases)
{
  //get all the curves and vertices of volumes in list
  DLIList<RefVertex*> verts;
  DLIList<RefEdge*> curves;

  RTree<RefEdge*> a_tree(high_tol);

  int i,j;
  for( i=ref_vols.size(); i--; )
  {
    RefVolume *tmp_vol = ref_vols.get_and_step();
    tmp_vol->ref_vertices( verts );
    
    curves.clean_out();
    tmp_vol->ref_edges( curves );
    for( j=curves.size(); j--; )
    {
      RefEdge *tmp_edge = curves.get_and_step();
      a_tree.add( tmp_edge );
    }
  }

  ProgressTool *progress_ptr = NULL;
  int total_verts = verts.size();
  if (total_verts > 5)
  {
    progress_ptr = AppUtil::instance()->progress_tool();
    assert(progress_ptr != NULL);
    progress_ptr->start(0, 100, "Finding Near Coincident Vertex-Curve Pairs", 
      NULL, CUBIT_TRUE, CUBIT_TRUE);
  }
    
  double curr_percent = 0.0;
  int processed_verts = 0;
  std::multimap<double, dist_vert_curve_struct> distance_vertex_curve_map; 

  //for each vertex
  for( i=verts.size(); i--; )
  {
    processed_verts++;
    if ( progress_ptr != NULL )
    {
      curr_percent = ((double)(processed_verts))/((double)(total_verts));
      progress_ptr->percent(curr_percent);
    }

    RefVertex *tmp_vert = verts.get_and_step();
    RefVolume *v1 = tmp_vert->ref_volume();
    CubitBox vertex_box ( tmp_vert->coordinates(), tmp_vert->coordinates() );
    DLIList<RefEdge*> close_curves;
    a_tree.find( vertex_box, close_curves );

    CubitVector vertex_xyz = tmp_vert->coordinates(); 

    for( j=close_curves.size(); j--; )
    {
      RefEdge *tmp_edge = close_curves.get_and_step();
      RefVolume *v2 = tmp_edge->ref_volume();

      bool check_distance = true;
      if(filter_same_volume_cases && v1 && v2 && v1 == v2)
        check_distance = false;
      if(check_distance)
      {
        CubitVector closest_location;
        tmp_edge->closest_point_trimmed( vertex_xyz, closest_location );
        double distance = closest_location.distance_between( vertex_xyz );

        if( distance >= low_tol && distance <= high_tol )
        {
          dist_vert_curve_struct tmp_struct;
          tmp_struct.dist = distance;
          tmp_struct.vert = tmp_vert;
          tmp_struct.edge = tmp_edge; 

          distance_vertex_curve_map.insert( std::multimap<double, dist_vert_curve_struct>::
                                            value_type( distance, tmp_struct ));
        }
      }
    }
  }

  if ( progress_ptr != NULL )
    progress_ptr->end();

  //std::set<dist_vert_curve_struct, vert_curve_dist_sort>::iterator iter, upper_iter; 
  std::multimap<double, dist_vert_curve_struct>::reverse_iterator iter;
  
  iter = distance_vertex_curve_map.rbegin();
  for(; iter!=distance_vertex_curve_map.rend(); iter++ )
  {
    distances.append( (*iter).second.dist );
    ref_verts.append( (*iter).second.vert );
    ref_edges.append( (*iter).second.edge );
  }

  return CUBIT_SUCCESS;
}


struct dist_vert_surf_struct
{
  double dist;
  RefVertex *vert;
  RefFace *face;
};

struct vert_surf_dist_sort
{
  bool operator()( dist_vert_surf_struct a, dist_vert_surf_struct b ) const
  {
    if( a.dist < b.dist )
      return true;
    else if( a.dist > b.dist )
      return false;
    else 
    {
      if( a.vert < b.vert )
        return true;
      else if( a.vert > b.vert )
        return false;
      else if( a.face < b.face )
        return true;
      else if( a.face > b.face )
        return false;
    }
    return false;
  }
};


CubitStatus GeomMeasureTool::find_near_coincident_vertex_surface_pairs( 
                                DLIList<RefVolume*> &ref_vols,
                                DLIList<RefFace*> &ref_faces,
                                DLIList<RefVertex*> &ref_verts,
                                DLIList<double> &distances,
                                double low_tol,
                                double high_tol,
                                bool filter_same_volume_cases)
{
  //get all the curves and vertices of volumes in list
  DLIList<RefVertex*> verts;
  DLIList<RefFace*> faces;

  AbstractTree<RefFace*> *a_tree = new RTree<RefFace*>( high_tol );

  int i,j;
  for( i=ref_vols.size(); i--; )
  {
    RefVolume *tmp_vol = ref_vols.get_and_step();
    tmp_vol->ref_vertices( verts );
    
    faces.clean_out();
    tmp_vol->ref_faces( faces );
    // Populate the Surface AbstractTree
    for( j=faces.size(); j--; )
    {
      RefFace *tmp_face = faces.get_and_step();
      a_tree->add( tmp_face );
    }
  }

  ProgressTool *progress_ptr = NULL;
  int total_verts = verts.size();
  if (total_verts > 50)
  {
    progress_ptr = AppUtil::instance()->progress_tool();
    assert(progress_ptr != NULL);
    progress_ptr->start(0, 100, "Finding Near Coincident Vertex-Surface Pairs", 
      NULL, CUBIT_TRUE, CUBIT_TRUE);
  }
  double curr_percent = 0.0;
  int processed_verts = 0;

  std::multimap<double, dist_vert_surf_struct> distance_vertex_surface_map; 

  //for each vertex
  for( i=verts.size(); i--; )
  {
    processed_verts++;
    if ( progress_ptr != NULL )
    {
      curr_percent = ((double)(processed_verts))/((double)(total_verts));
      progress_ptr->percent(curr_percent);
    }

    RefVertex *tmp_vert = verts.get_and_step();
    RefVolume *v1 = tmp_vert->ref_volume();
    CubitBox vertex_box ( tmp_vert->coordinates(), tmp_vert->coordinates() );
    DLIList<RefFace*> close_faces;
    a_tree->find( vertex_box, close_faces);

    CubitVector vertex_xyz = tmp_vert->coordinates(); 

    for( j=close_faces.size(); j--; )
    {
      RefFace *tmp_face = close_faces.get_and_step();
      RefVolume *v2 = tmp_face->ref_volume();

      bool check = true;
      if(filter_same_volume_cases && v1 && v2 && v1 == v2)
        check = false;

      if(check)
      {
        DLIList<RefVertex*> tmp_verts;
        tmp_face->ref_vertices( tmp_verts );
        if( tmp_verts.is_in_list( tmp_vert ) ) 
          continue;

        CubitVector closest_location;
        tmp_face->find_closest_point_trimmed( vertex_xyz, closest_location );
        double distance = closest_location.distance_between( vertex_xyz );

        if( distance > low_tol && distance < high_tol )
        {
          dist_vert_surf_struct tmp_struct;
          tmp_struct.dist = distance;
          tmp_struct.vert = tmp_vert;
          tmp_struct.face = tmp_face; 
          distance_vertex_surface_map.insert( std::multimap<double, dist_vert_surf_struct>::
                                      value_type( distance, tmp_struct ));
        }
      }
    }
  }

  if ( progress_ptr != NULL )
    progress_ptr->end();

  std::multimap<double, dist_vert_surf_struct>::reverse_iterator iter;

  iter = distance_vertex_surface_map.rbegin();
  for(; iter!=distance_vertex_surface_map.rend(); iter++ )
  {
    distances.append( (*iter).second.dist );
    ref_verts.append( (*iter).second.vert );
    ref_faces.append( (*iter).second.face);
  }

  delete a_tree;

  return CUBIT_SUCCESS;
}





