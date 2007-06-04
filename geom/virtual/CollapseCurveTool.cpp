//-------------------------------------------------------------------------
//- Filename:       CollapseCurveTool
//- Purpose:  To collapse small curves for preparing for mesh
//-      "collapse curve <id> vertex <id>\n",
//-
//- Creator:       Brett Clark
//- Creation date: 01/10/2006
//-------------------------------------------------------------------------

// ********** BEGIN STANDARD INCLUDES         **********
// ********** END STANDARD INCLUDES           **********
                                                                                
// ********** BEGIN MOTIF INCLUDES            **********
// ********** END MOTIF INCLUDES              **********
                                                                                
// ********** BEGIN ACIS INCLUDES             **********
// ********** END ACIS INCLUDES               **********
                                                                                
// ********** BEGIN CUBIT INCLUDES            **********
#include "GMem.hpp"                                                           
#include "RefVertex.hpp"
#include "RefEdge.hpp"
#include "DLIList.hpp"
#include "CoEdge.hpp"
#include "CoFace.hpp"
#include "Shell.hpp"
#include "Point.hpp"
#include "Loop.hpp"
#include "RefFace.hpp"
#include "PartitionTool.hpp"
#include "CompositeTool.hpp"
#include "CollapseCurveTool.hpp"

/*
*/
// ********** END CUBIT INCLUDES              **********

CollapseCurveTool *CollapseCurveTool::instance_ = NULL;

// ********** BEGIN PUBLIC FUNCTIONS          **********
CollapseCurveTool *CollapseCurveTool::instance()
{
  if (instance_ == NULL) 
    instance_ = new CollapseCurveTool();
                                                                                
  return instance_;
}

CubitStatus CollapseCurveTool::collapse_curve(DLIList <RefEdge*> ref_edge_list, 
                                              DLIList<RefVertex*> ref_vertex_list,
                                              int ignore_surfaces)
{
  CubitStatus result = CUBIT_SUCCESS;

  RefEdge *edge_to_collapse = ref_edge_list.get();
  RefVertex *keep_vertex = NULL;
  RefVertex *discard_vertex = NULL;
  CoEdge *exit_coedge = NULL;
  CoEdge *enter_coedge = NULL;
  DLIList<RefEdge*> curves_to_partition;
  DLIList<RefFace*> surfaces_to_partition;
  DLIList<CubitVector> partition_positions;
  DLIList<CubitVector> keep_vecs;
  DLIList<CubitVector> partition_curve_vecs;
  DLIList<RefFace*> adjacent_surfaces;
  DLIList<double> angles;
  DLIList<double> arc_lengths;
  DLIList<RefEdge*> ref_edges_on_discard_vertex;
  DLIList<RefEdge*> ref_edges_on_keep_vertex;
  DLIList<RefVertex*> new_partition_vertices;
  DLIList<RefEdge*> new_partition_edges;
  RefEdge *close_edge = NULL;
  RefEdge *far_edge = NULL;
  DLIList<RefEdge*> edges_to_composite;
  int keep_vertex_specified = 0;
  int discard_valency = 0;
  int keep_valency = 0;
  int finished = 0;
  double arc_length = 0;
  double u = 0, v = 0;
  CubitVector left_vec(0,0,0), right_vec(0,0,0);
  CubitVector position_on_left_curve(0,0,0), position_on_right_curve(0,0,0);
  CubitVector discard_pos(0,0,0), keep_pos(0,0,0);

  if(edge_to_collapse == NULL)
  {
    PRINT_ERROR("No curve was found to collapse.\n");
    finished = 1;
    result = CUBIT_FAILURE;
  }

  if(!finished)
  {
    // Make sure we have a vertex to keep.
    if(ref_vertex_list.size() > 0)
    {
      keep_vertex_specified = 1;
      keep_vertex = ref_vertex_list.get();
    }
    // We need to choose a vertex to keep.
    else
    {
      // Arbitrarily choose start vertex.
      keep_vertex = edge_to_collapse->start_vertex();
    }

    // Get the vertex that will go away.
    if(keep_vertex == edge_to_collapse->start_vertex())
    {
      discard_vertex = edge_to_collapse->end_vertex();
    }
    else if(keep_vertex == edge_to_collapse->end_vertex())
    {
      discard_vertex = edge_to_collapse->start_vertex();
    }
    else
    {
      PRINT_ERROR("Vertex to keep was not found on the curve being collapsed.\n");
      result = CUBIT_FAILURE;
      finished = 1;
    }
  }

  if(!finished)
  {
    // Get the valency of the keep and discard vertices.
    discard_vertex->ref_edges(ref_edges_on_discard_vertex);
    keep_vertex->ref_edges(ref_edges_on_keep_vertex);
    discard_valency = ref_edges_on_discard_vertex.size();
    keep_valency = ref_edges_on_keep_vertex.size();

    // Now do some error checking and also some logic to try to 
    // pick the best vertex to keep/discard.

    // If one of the valencies is 2 just composite the vertex out.
    if(discard_valency == 2)
    {
      CompositeTool::instance()->composite(ref_edges_on_discard_vertex);
      finished = 1;
    }
    else if (keep_valency == 2)
    {
      CompositeTool::instance()->composite(ref_edges_on_keep_vertex);
      finished = 1;
    }
    else
    {
      // Make sure that at least one of the vertices is qualified to
      // be the discard vertex (valency of 3 or 4).  The keep vertex
      // really doesn't have any restrictions.
      if((discard_valency > 4 || discard_valency < 3) &&
        (keep_valency > 4 || keep_valency < 3))
      {
        PRINT_ERROR("Cannot currently collapse curves where one of the vertices does not have a valency equal to 3 or 4.\n");
        result = CUBIT_FAILURE;
        finished = 1;
      }
      else
      {
        if(keep_vertex_specified)
        {
          if(discard_valency < 3 || discard_valency > 4)
          {
            PRINT_ERROR("Cannot currently collapse curves where the discard vertex valency is not 3 or 4.\n"
                  "Try specifying the keep vertex so that the discard vertex is 3 or 4.\n");
            result = CUBIT_FAILURE;
            finished = 1;
          }
        }
        else
        {
          // The user didn't specify a keep vertex so we can try to choose the best one.
          
          int swap_vertices = 0;
          if(discard_valency < 5 && discard_valency > 2 &&
            keep_valency < 5 && keep_valency > 2)
          {
            // Either vertex can be the discard/keep vertex so choose the one with the
            // lower valency as the discard vertex because it will require fewer operations.
            if(discard_valency > keep_valency)
            {
              swap_vertices = 1;
            }
          }
          else
          {
            // Only one of the vertices can be the discard vertex so pick it.
            if(discard_valency > 4 || discard_valency < 3)
            {
              swap_vertices = 1;
            }
          }

          if(swap_vertices)
          {
            // Swap the vertices.
            RefVertex *tmp = discard_vertex;
            discard_vertex = keep_vertex;
            keep_vertex = tmp;

            // Make sure to refresh the discard vertex edge list because
            // it is used below.
            ref_edges_on_discard_vertex.clean_out();
            discard_vertex->ref_edges(ref_edges_on_discard_vertex);
          }
        }
      }
    }
  }

  // Look at all of the coedges on the curve being collapsed and
  // throw out any that are on nonmanifold surfaces.  The collapse
  // curve may be on a boundary of a merged surface.  This is ok
  // but we want to make sure we don't involve this surface
  // in the partitions and composites so we will remove from the
  // list any coedges that are hooked to merged surfaces.
  if(!finished)
  {
    DLIList<CoEdge*> collapse_curve_coedges;
    edge_to_collapse->get_co_edges(collapse_curve_coedges);

    int initial_size = collapse_curve_coedges.size();
    while(collapse_curve_coedges.size() > 2 && initial_size > 0)
    {
      for(int g=collapse_curve_coedges.size(); g--;)
      {
        CoEdge *cur_coedge = collapse_curve_coedges.get_and_step();
        Loop *loop_ptr = cur_coedge->get_loop_ptr();
        if(loop_ptr)
        {
          RefFace *ref_face_ptr = loop_ptr->get_ref_face_ptr();
          if(ref_face_ptr)
          {
            DLIList<CoFace*> coface_list;
            ref_face_ptr->co_faces(coface_list);
            if(coface_list.size() > 1)
            {
              collapse_curve_coedges.remove(cur_coedge);
              g = 0;
            }
          }
        }
      }

      // Keep from looping infinitely.
      initial_size--;
    }

    // If we get to this point and have more than 2 coedges left
    // in the list we are not sure what to do.
    if(collapse_curve_coedges.size() != 2)
    {
      PRINT_ERROR("Currently can only collapse curves with manifold topology.\n");
      result = CUBIT_FAILURE;
      finished = 1;
    }
    else
    {
      // Get the collapse curve coedges entering and leaving the discard vertex.
      exit_coedge = collapse_curve_coedges.get_and_step();
      enter_coedge = collapse_curve_coedges.get();
      if(enter_coedge->end_vertex() != discard_vertex)
      {
        CoEdge *tmp = exit_coedge;
        exit_coedge = enter_coedge;
        enter_coedge = tmp;
      }
    }
  }

  // Next we need to explore the topology around the discard vertex
  // so that we can get the edges and surfaces that will be involved
  // with the partitioning and compositing.  We will identify a "left"
  // and "right" edge coming out of the discard vertex and adjacent
  // surfacs to these edges.
  if(!finished)
  {
    // Get these values for use later on.
    discard_pos = discard_vertex->get_point_ptr()->coordinates();
    keep_pos = keep_vertex->get_point_ptr()->coordinates();
    CubitVector discard_to_keep = keep_pos - discard_pos;
    arc_length = edge_to_collapse->get_arc_length();

    // Depending on whether the discard vertex has valency of 3 or 4 we will
    // either partition 1 or 2 of the curves coming into the discard vertex
    // respectively.  Set up lists so that below we can just loop through the
    // partitioning and compositing for either the 3 or 4 valency case.

    // "Left" and "right" will be defined as if you were standing on the
    // keep vertex looking at the discard vertex.
    Loop *left_loop = enter_coedge->get_loop_ptr();
    Loop *right_loop = exit_coedge->get_loop_ptr();

    // We need to get these two coedges because the chains defined by the "next" and 
    // "previous" pointers in the CoEdge objects are not circular (you can have a NULL
    // "previous" or "next" pointer).  We will manage the circularity manually.
    CoEdge *left_start = dynamic_cast<CoEdge*>(left_loop->get_first_sense_entity_ptr()); 
    CoEdge *right_end = dynamic_cast<CoEdge*>(right_loop->get_last_sense_entity_ptr()); 
    
    CoEdge *left_coedge = dynamic_cast<CoEdge*>(enter_coedge->next());
    if(left_coedge == NULL)
    {
      left_coedge = left_start;
    }
    CoEdge *right_coedge = dynamic_cast<CoEdge*>(exit_coedge->previous());
    if(right_coedge == NULL)
    {
      right_coedge = right_end;
    }

    RefEdge *left_edge = left_coedge->get_ref_edge_ptr();
    RefEdge *right_edge = right_coedge->get_ref_edge_ptr();
    RefFace *left_common_face = left_edge->common_ref_face(edge_to_collapse);
    RefFace *right_common_face = right_edge->common_ref_face(edge_to_collapse);

    double left_arc_length = left_edge->get_arc_length();
    double right_arc_length = right_edge->get_arc_length();
    int left_ok = 1;
    int right_ok = 1;

    DLIList<RefEdge*> left_face_edges;
    DLIList<RefEdge*> right_face_edges;
    left_common_face->ref_edges(left_face_edges);
    right_common_face->ref_edges(right_face_edges);
    if(left_face_edges.size() < 3 || right_face_edges.size() < 3)
    {
      PRINT_ERROR("Cannot collapse a curve that bounds a face with only two edges.\n");
      finished = 1;
    }
    else
    {
      // We use the length of the curve being collapsed as the distance
      // to partition the adjacent curves connected to it.  If the
      // adjacent curves are shorter than the curve being collapsed we
      // will bail because the user should probably be getting rid
      // of the adjacent edges instead or at least first.

      // Get the length of the adjacent curves.
      
      // If the adjacent curve is within resabs of the length of 
      // the curve being collapsed we will just snap to the end
      // of the adjacent curve for our partition position.
      if(arc_length > left_arc_length + GEOMETRY_RESABS)
      {
        left_ok = 0;
      }
      if(arc_length > right_arc_length + GEOMETRY_RESABS)
      {
        right_ok = 0;
      }

      // If neither curve is ok we need to bail.
      if(!left_ok && !right_ok)
      {
        PRINT_ERROR("Curve being collapsed is too long compared to adjacent curves.\n");
        finished = 1;
      }
    }

    // If it looks like the lengths of the adjacent curves are
    // ok we will go ahead and try to find a partition position
    // on them.
    if(!finished)
    {
      if(left_ok)
      {
         // First see if we can just use the end point of the
         // adjacent curve as the partition position.
         if(fabs(left_arc_length-arc_length) < GEOMETRY_RESABS)
         {
           if(left_edge->start_vertex() == discard_vertex)
             position_on_left_curve = left_edge->end_vertex()->coordinates();
           else
             position_on_left_curve = left_edge->start_vertex()->coordinates();
         }
         else
         {
           result = this->position_from_length(left_edge, discard_vertex,
                arc_length, position_on_left_curve);
         }
      }
      if(result == CUBIT_SUCCESS && right_ok)
      {
         // First see if we can just use the end point of the
         // adjacent curve as the partition position.
         if(fabs(right_arc_length-arc_length) < GEOMETRY_RESABS)
         {
           if(right_edge->start_vertex() == discard_vertex)
             position_on_right_curve = right_edge->end_vertex()->coordinates();
           else
             position_on_right_curve = right_edge->start_vertex()->coordinates();
         }
         else
         {
            result = this->position_from_length(right_edge, discard_vertex,
                  arc_length, position_on_right_curve);
         }
      }
      if(result == CUBIT_FAILURE)
      {
        PRINT_ERROR("Was unable to locate appropriate partition points on curves adjacent to curve being collapased.\n");
        finished = 1;
      }
    }
      
    if(!finished)
    {
      // Get the vectors from the discard vertex to the potential partition locations.
      CubitVector left_vec = position_on_left_curve - discard_pos;
      CubitVector right_vec = position_on_right_curve - discard_pos;

      // Calculate the angles between the left/right edge and the
      // edge being collapsed.  I am doing it this way rather than just
      // calling RefEdge::angle_between() because I want the angle
      // calculation done out at the potential partition location.
      // This will step over any small kinks in the curve near the
      // discard vertex that might otherwise give misleading angles
      // that don't represent what is happening out by the potential
      // partition position.
      left_common_face->u_v_from_position(discard_pos, u, v);
      CubitVector left_normal = left_common_face->normal_at(discard_pos, NULL, &u, &v);
      double left_angle = left_normal.vector_angle(left_vec, discard_to_keep);

      right_common_face->u_v_from_position(discard_pos, u, v);
      CubitVector right_normal = right_common_face->normal_at(discard_pos, NULL, &u, &v);
      double right_angle = right_normal.vector_angle(discard_to_keep, right_vec);

      // 3 valency case: 
      // We only need to partition one of the curves (left or right) and
      // the corresponding surface.
      if(ref_edges_on_discard_vertex.size() == 3)
      {
        int use_left = 0;
        // We can use either adjacent curve.
        if(left_ok && right_ok)
        {
          // Choose the side with the smaller angle as this will in general
          // give better angles for doing the partitioning on the surface.
          if(left_angle < right_angle)
          {
            use_left = 1;
          }
        }
        else if(left_ok)
        {
          use_left = 1;
        }
        if(use_left)
        {
          curves_to_partition.append(left_edge);
          surfaces_to_partition.append(left_common_face);
          angles.append(left_angle);
          partition_positions.append(position_on_left_curve);
          arc_lengths.append(arc_length);
 
          // These vectors will be used in calculating a bisector direction
          // below if necessary.
          keep_vecs.append(-discard_to_keep);
          partition_curve_vecs.append(left_vec);
        }
        else
        {
          curves_to_partition.append(right_edge);
          surfaces_to_partition.append(right_common_face);
          angles.append(right_angle);
          partition_positions.append(position_on_right_curve);
          arc_lengths.append(arc_length);
 
          // These vectors will be used in calculating a bisector direction
          // below if necessary.
          keep_vecs.append(discard_to_keep);
          partition_curve_vecs.append(-right_vec);
        }
      }
      else if(ref_edges_on_discard_vertex.size() == 4)
      {
        // We have to partition both left and right curves so make
        // sure we can.
        if(!left_ok || !right_ok)
        {
          PRINT_ERROR("One of the curves adjacent to the collapse curve is not long enough for the collapse operation.\n");
          finished = 1;
        }
        else
        {
          // Both curves (and surfaces) adjacent to the collapse curve (left and right)
          // will need to be partitioned so add them to the lists.

          curves_to_partition.append(left_edge);
          curves_to_partition.append(right_edge);
          surfaces_to_partition.append(left_common_face);
          surfaces_to_partition.append(right_common_face);
          angles.append(left_angle);
          angles.append(right_angle);
          keep_vecs.append(-discard_to_keep);
          keep_vecs.append(discard_to_keep);
          partition_curve_vecs.append(left_vec);
          partition_curve_vecs.append(-right_vec);
          partition_positions.append(position_on_left_curve);
          partition_positions.append(position_on_right_curve);
          arc_lengths.append(arc_length);
          arc_lengths.append(arc_length);
        }
      }
      else
      {
        PRINT_ERROR("Currently can only collapse curves with 3 or 4 valency vertices.\n");
        result = CUBIT_FAILURE;
        finished = 1;
      }
    }
  }

  if(!finished)
  {
    int num_reps = curves_to_partition.size();
    curves_to_partition.reset();
    surfaces_to_partition.reset();
    for(int n=0; n<num_reps && !finished; ++n)
    {
      RefEdge *side_edge = curves_to_partition.get_and_step();
      RefFace *side_face = surfaces_to_partition.get_and_step();

      // Now we need to find the face on the other side of the edge.
      // Because there may be edges on the boundaries of merged surfaces
      // we want to find the
      // face on the left/right edge that isn't the face shared by
      // the collapse edge and isn't nonmanifold.
      DLIList<CoEdge*> side_edge_coedges;
      side_edge->co_edges(side_edge_coedges);
      RefFace *adjacent_face = NULL;
      DLIList<RefFace*> possible_adj_faces;

      // First loop through and get all of the potential faces.
      for(int r=side_edge_coedges.size(); r--;)
      {
        CoEdge *cur_coedge = side_edge_coedges.get_and_step();
        if(cur_coedge &&
          cur_coedge->get_loop_ptr() &&
          cur_coedge->get_loop_ptr()->get_ref_face_ptr())
        {
          RefFace *cur_ref_face = cur_coedge->get_loop_ptr()->get_ref_face_ptr();
          if(cur_ref_face != side_face)
          {
            DLIList<CoFace*> coface_list;
            cur_ref_face->co_faces(coface_list);

            // Along with checking for whether the face is manifold we need
            // to check if it belongs to the same volume as the side face.
            // We have to check this because we can't composite faces
            // from different volumes and these faces will be involved in
            // a composite below.
            if(coface_list.size() == 1 &&
              side_face->ref_volume() == cur_ref_face->ref_volume())
            {
              possible_adj_faces.append(cur_ref_face);
            }
          }
        }
      }

      // If we ended up with more than one face in the list it isn't clear
      // what we should do so bail out.
      if(possible_adj_faces.size() != 1)
      {
        PRINT_ERROR("Couldn't figure out how to perform the collapse curve with the current topology.\n");
        result = CUBIT_FAILURE;
        finished = 1;
      }
      else
      {
        adjacent_surfaces.append(possible_adj_faces.get());
      }
    }
  }

  if(!finished)
  {
    // At this point we should know which curves and surfaces we need
    // to partition.

    int num_reps = curves_to_partition.size();
    curves_to_partition.reset();
    surfaces_to_partition.reset();
    adjacent_surfaces.reset();
    angles.reset();
    keep_vecs.reset();
    partition_curve_vecs.reset();
    partition_positions.reset();
    arc_lengths.reset();

    for(int i=0; i<num_reps && !finished; ++i)
    {
      RefEdge *curve_to_partition = curves_to_partition.get_and_step();

      // As we process each curve remove it from the list so that
      // at the end we can take the last one in the list and composite
      // it with the curve being collapsed.
      ref_edges_on_discard_vertex.remove(curve_to_partition);

      CubitVector position_on_curve = partition_positions.get_and_step();

      // Partition the curve adjacent to the collapse curve.
      RefEdge *edge1, *edge2;
      RefVertex *new_vertex = PartitionTool::instance()->
        partition( curve_to_partition, position_on_curve, edge1, edge2 );

      // Keep track of these in case we need to undo the partitioning.
      if(new_vertex)
      {
        new_partition_vertices.append(new_vertex);
      }

      // Get the two new curves and classify them as close and far.
      if(edge1 && edge2)
      {
        close_edge = edge1;
        far_edge = edge2;
        if(close_edge->start_vertex() != discard_vertex &&
          close_edge->end_vertex() != discard_vertex)
        {
          RefEdge *tmp = close_edge;
          close_edge = far_edge;
          far_edge = tmp;
        }
      }
      else
      {
        // If the partition didn't create a new vertex and two new curves
        // (the case when the partition position lands on an existing
        // vertex) just classify the curve as "close" and set the 
        // "far" curve to NULL. 
        close_edge = curve_to_partition;
        far_edge = NULL;
      }

      RefFace *surface_to_partition = surfaces_to_partition.get_and_step();

      DLIList<CubitVector*> positions;

      // Add the point on the curve we just partitioned.  The other
      // point we normally need to add is the keep_vertex.  However,
      // if the angle between the collapse curve and the curve we 
      // just partitioned dictates that we need to introduce more points 
      // (this would be cases where the angle is large and just partitioning the
      // surface with one curve--two points--would result in a very skinny
      // surface) we need to add them before adding the keep_vertex.  Therefore, check
      // that now and add the extra points if needed.
      positions.append(new CubitVector(position_on_curve));

      double cur_angle = angles.get_and_step();
      arc_length = arc_lengths.get_and_step();

      // These vectors are only used in the block below if we need to
      // add extra points but we need to always retrieve them from the lists
      // so that the pointer in the list stays in sync with the "for" loop.
      CubitVector keep_vec = keep_vecs.get_and_step();
      CubitVector partition_vec = partition_curve_vecs.get_and_step();

      // Greater than 3*PI/2--add 3 interior points.  Two of these
      // points will be generated by projecting from the discard vertex
      // into the surface normal to the vector from the discard vertex
      // to the keep point and new partition point respectively.  The third
      // point will be obtained by projecting from the discard point in the
      // direction of the bisector angle of the two previous projections.
      if(cur_angle > 4.71)
      {
        // Get the u,v position of the discard vertex on this surface.
        surface_to_partition->u_v_from_position(discard_pos, u, v);
        
        // Get the normal at the discard vertex.
        CubitVector normal = surface_to_partition->normal_at(discard_pos, NULL, &u, &v);
        normal.normalize();
        
        // We need to calculate a bisector direction between the
        // collape edge and the edge we just partitioned.  We will do 
        // this using the face normal at the discard vertex and the vectors 
        // we calculated previously with the partition points.  Cross
        // the face normal into the the direction vectors at the 
        // discard and partition points to get two vectors pointing
        // into the face and then average those two to get the
        // bisector direction.  This is all approximate but should 
        // be sufficient for locating another point for partitioning 
        // the surface.

        // I am not normalizing the result here because they should
        // be roughly the same length and it shouldn't affect 
        // the average too much.
        CubitVector vec1 = normal * partition_vec;
        CubitVector vec2 = normal * keep_vec;

        // Get the bisector direction.
        CubitVector bisector_dir = vec1 + vec2;
        bisector_dir.normalize();

        // Now normalise these because they will be used to
        // project two of the new interior points.
        vec1.normalize();
        vec2.normalize();

        CubitVector new_pos1 = discard_pos + (arc_length*vec1);
        CubitVector mid_pos = discard_pos + (arc_length * bisector_dir);
        CubitVector new_pos2 = discard_pos + (arc_length*vec2);

        // Use the u,v from the discard vertex because it is fairly
        // close to the new position and will at least provide a
        // meaningful starting point.
        double save_u = u, save_v = v;
        surface_to_partition->move_to_surface(new_pos1, &u, &v);
        u = save_u; v = save_v;
        surface_to_partition->move_to_surface(mid_pos, &u, &v);
        u = save_u; v = save_v;
        surface_to_partition->move_to_surface(new_pos2, &u, &v);

        // Add the new position to the list of partition points.
        positions.append(new CubitVector(new_pos1));
        positions.append(new CubitVector(mid_pos));
        positions.append(new CubitVector(new_pos2));
      }
      // Greater than 3*PI/4 and less than 3*PI/2--add one interior point
      else if(cur_angle > 2.4)
      {
        CubitVector third_pt;
        
        // Get the u,v position of the discard vertex on this surface.
        surface_to_partition->u_v_from_position(discard_pos, u, v);
        
        // Get the normal at the discard vertex.
        CubitVector normal = surface_to_partition->normal_at(discard_pos, NULL, &u, &v);
        normal.normalize();
        
        // We need to calculate a bisector direction between the
        // collape edge and the edge we just partitioned.  We will do 
        // this using the face normal at the discard vertex and the vectors 
        // we calculated previously with the partition points.  Cross
        // the face normal into the the direction vectors at the 
        // discard and partition points to get two vectors pointing
        // into the face and then average those two to get the
        // bisector direction.  This is all approximate but should 
        // be sufficient for locating another point for partitioning 
        // the surface.

        // I am not normalizing the result here because they should
        // be roughly the same length and it shouldn't affect 
        // the average too much.
        CubitVector vec1 = normal * keep_vec;
        CubitVector vec2 = normal * partition_vec;

        // Get the bisector direction.
        CubitVector bisector_dir = vec1 + vec2;
        bisector_dir.normalize();

        // Project from the discard vertex in the direction of the
        // bisector direction to get a new point for partitioning
        // the surface.
        CubitVector new_pos = discard_pos + (arc_length * bisector_dir);

        // Use the u,v from the discard vertex because it is fairly
        // close to the new position and will at least provide a
        // meaningful starting point.
        surface_to_partition->move_to_surface(new_pos, &u, &v);

        // Add the new position to the list of partition points.
        positions.append(new CubitVector(new_pos));
      }

      // Finally, add the keep_vertex to the list.
      positions.append(new CubitVector(keep_vertex->get_point_ptr()->coordinates()));
                                                                                    
      DLIList<RefEdge*> new_edges;
      PartitionTool::instance()->
          insert_edge( surface_to_partition, positions, CUBIT_FALSE, new_edges,
                        0, &arc_length);

      // Keep for later in case we need to clean up.
      if(new_edges.size() > 0)
      {
        new_partition_edges += new_edges;
      }

      delete positions.get_and_step();
      delete positions.get_and_step();

      if(cur_angle > 4.71)
      {
        delete positions.get_and_step();
        delete positions.get_and_step();
        delete positions.get_and_step();
      }
      else if(cur_angle > 2.4)
      {
        delete positions.get();
      }

      RefFace *new_small_face = edge_to_collapse->common_ref_face(close_edge);

      if(!new_small_face)
      {
        PRINT_ERROR("Failed to do the partition surface operation of the collapse curve.\n");
        result = CUBIT_FAILURE;
        finished = 1;
      }
      else
      {
        DLIList<RefFace*> result_faces;
        DLIList<RefFace*> faces_to_composite;

        // Used below if "ignore" keyword was specified.
        int new_small_face_id = new_small_face->id();

        faces_to_composite.append(new_small_face);
        faces_to_composite.append(adjacent_surfaces.get_and_step());

        RefFace *new_comp_face = CompositeTool::instance()->composite(faces_to_composite);

        CompositeSurface* csurf = NULL;
        
        if(new_comp_face)
        {
          csurf = dynamic_cast<CompositeSurface*>(new_comp_face->get_surface_ptr());
        }

        if(!new_comp_face || !csurf)
        {
          PRINT_ERROR("Failed to do the composite surface operation of the collapse curve.\n");
          result = CUBIT_FAILURE;
          finished = 1;
        }
        else
        {
          if(ignore_surfaces)
          {
            csurf->ignore_surface(new_small_face_id);
          }

          if(far_edge)
          {
            for(int k=new_edges.size(); k--;)
            {
              edges_to_composite.append(new_edges.get_and_step());
            }
            edges_to_composite.append(far_edge);

            if(edges_to_composite.size() > 1)
            {
              CompositeTool::instance()->composite(edges_to_composite);
            }
          }

          edges_to_composite.clean_out();
        }
      }
    }

    if(!finished)
    {
      ref_edges_on_discard_vertex.remove(edge_to_collapse);

      // Now there should only be one edge in the ref_edges_on_discard_vertex
      // list.  It should be the edge that hasn't had any modifications
      // done to it.  We finally want to composite it with the edge
      // being collapsed.
      if(ref_edges_on_discard_vertex.size() != 1)
      {
        PRINT_ERROR("Wasn't able to complete collapse operation.\n");
        result = CUBIT_FAILURE;
        finished = 1;
      }
      else
      {
        edges_to_composite.append(ref_edges_on_discard_vertex.get());
        edges_to_composite.append(edge_to_collapse);
        CompositeTool::instance()->composite(edges_to_composite);
      }
    }
  }

  if(result == CUBIT_FAILURE)
  {
    int i;
    for(i=new_partition_edges.size(); i--;)
    {
      CompositeTool::instance()->remove_edge(new_partition_edges.get_and_step(), true);
    }
    for(i=new_partition_vertices.size(); i--;)
    {
      CompositeTool::instance()->remove_vertex(new_partition_vertices.get_and_step(), true);
    }
  }

  return result;
}

CubitStatus CollapseCurveTool::position_from_length(RefEdge *edge,
                                                  RefVertex *root_vertex,
                                                  double  arc_length,
                                                  CubitVector& v_new)
{
  CubitStatus result;
  double sense = 1.0;
  if (root_vertex != edge->start_vertex())
    sense = -1.0;
  CubitVector v_root = root_vertex->get_point_ptr()->coordinates();
  result = edge->point_from_arc_length(v_root, arc_length*sense, v_new);
  return result;
}


// ********** END PUBLIC FUNCTIONS            **********
                                                                                
// ********** BEGIN PROTECTED FUNCTIONS       **********
// ********** END PROTECTED FUNCTIONS         **********
                                                                                
// ********** BEGIN PRIVATE FUNCTIONS         **********
CollapseCurveTool::CollapseCurveTool()
{
}

CollapseCurveTool::~CollapseCurveTool()
{
}

