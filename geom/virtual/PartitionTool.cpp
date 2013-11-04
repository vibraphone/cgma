//-------------------------------------------------------------------------
// Filename      : PartitionTool.cpp
//
// Purpose       : Functions for splitting geometry using VG
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/13/98
//-------------------------------------------------------------------------

const int PTT_EDGE_POINT_COUNT = 2;
//When determining if a RefEdge lies within a Loop, this
//is the number of points on that RefEdge that are checked.  
  
#include "PartitionTool.hpp"
//#include "VirtualQueryEngine.hpp"
#include "DLIList.hpp"
#include "GfxDebug.hpp"
#include "GMem.hpp"
#include "CastTo.hpp"

#include "Body.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include "Loop.hpp"
#include "Chain.hpp"
#include "Shell.hpp"

#include "CoEdge.hpp"
#include "CoVertex.hpp"
#include "CoFace.hpp"
#include "CoVolume.hpp"

#include "Point.hpp"
#include "Curve.hpp"
#include "IntersectionTool.hpp"

#include "ModelQueryEngine.hpp"

#include "PartitionLump.hpp"
#include "PartitionSurface.hpp"
#include "PartitionCurve.hpp"
#include "PartPTCurve.hpp"
#include "PartitionPoint.hpp"
//#include "PartitionLump.hpp"
#include "PartitionLoop.hpp"
#include "PartitionCoEdge.hpp"
#include "SegmentedCurve.hpp"

#include "CompositeTool.hpp"
#include "CompositeEngine.hpp"
#include "CompositePoint.hpp"
#include "GeometryUtil.hpp"
#include "RefEntityName.hpp"

#include "PartitionEngine.hpp"

#include "BodySM.hpp"
#include "SubCurve.hpp"

#include "RefEntityFactory.hpp"
#include "GeometryQueryTool.hpp"
#include "MergeTool.hpp"
#include "CubitUndo.hpp"

#include "TDUPtr.hpp"
#include <vector>

#include "CompositeCombineEvent.hpp"
#include "AppUtil.hpp"

PartitionTool* PartitionTool::instance_ = NULL;

//Debugging output
void pt_print_bte_list( int debug_flag, DLIList<BasicTopologyEntity*>& edges, 
                        const char* trailing_string = "\n"  );
void pt_print_edge_list( int debug_flag, DLIList<RefEdge*>& edges,
                         const char* trailing_string = "\n" );
void pt_print_face_list( int debug_flag, DLIList<RefFace*>& faces,
                         const char* trailing_string = "\n" );
void pt_print_loop( int debug_flag, Loop* loop_ptr, 
                    const char* trailing_string = "\n" );
void pt_print_shell( int debug_flag, Shell* shell_ptr, 
                     const char* trailing_string = "\n" );

//Helper class for surface partitioning.
struct LoopIntersection { 
  Loop* loop_ptr; 
  RefVertex* vtx_ptr; 
    //So we can put it in a DLIList:
  LoopIntersection( int i = 0 );


  bool operator<( const LoopIntersection& ) const { assert(0); return false; }
  bool operator>( const LoopIntersection& ) const { assert(0); return false; }
  bool operator<=( const LoopIntersection& ) const { assert(0); return false; }
  bool operator>=( const LoopIntersection& ) const { assert(0); return false; }
  bool operator==( const LoopIntersection& ) const { assert(0); return true; }
  bool operator!=( const LoopIntersection& ) const { assert(0); return true; }
};
  //So we can put it in a DLIList:
LoopIntersection::LoopIntersection( int ) : loop_ptr(0), vtx_ptr(0) {}

//-------------------------------------------------------------------------
// Purpose       : constructor/destructor
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 12/13/98
//-------------------------------------------------------------------------
PartitionTool::PartitionTool(){}
PartitionTool::~PartitionTool()
{
  instance_ = NULL;
}

//-------------------------------------------------------------------------
// Purpose       : Split a curve
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/15/99
//-------------------------------------------------------------------------
RefVertex* PartitionTool::partition( RefEdge* edge_ptr,
                                      const CubitVector& split_point,
                                      RefEdge*& first_new_edge,
                                      RefEdge*& second_new_edge,
                                      CubitBoolean  )
{
  int i, j;
  first_new_edge = second_new_edge = NULL;

  CubitBoolean is_free_curve = 
      edge_ptr->num_parent_ref_entities() == 0 ? CUBIT_TRUE : CUBIT_FALSE;

  DLIList<RefFace*> faces;
  edge_ptr->ref_faces( faces );
  RefVertex* input_start_vtx = edge_ptr->start_vertex();
  RefVertex* input_end_vtx = edge_ptr->end_vertex();

    // partition each merged curve in the refedge
  DLIList<TopologyBridge*> bridge_list, curve_bridges, coedge_bridges, point_bridges;
  edge_ptr->bridge_manager()->get_bridge_list( bridge_list );
  RefVertex* new_vertex = 0;
  for( i = bridge_list.size(); i--; )
  {
    Curve* curve_ptr = dynamic_cast<Curve*>(bridge_list.get_and_step());
    double u = curve_ptr->u_from_position( split_point );
    TBPoint* pt = PartitionEngine::instance().insert_point( curve_ptr, u );
    if( !pt )
    {
        // try to back out any partitions we've already made.
      bridge_list.clean_out();
      if( new_vertex )
      {
        new_vertex->bridge_manager()->get_bridge_list(bridge_list);
        new_vertex->deactivated(CUBIT_TRUE);
      }
      
      while( bridge_list.size() )
        PartitionEngine::instance().remove_point( 
          dynamic_cast<PartitionPoint*>(bridge_list.pop()) );
      
      break;
    }
    
    if( new_vertex ) 
      new_vertex->bridge_manager()->add_bridge( pt );
    else
      new_vertex = GeometryQueryTool::instance()->make_RefVertex(pt);
  }  
  
  if (!new_vertex) // failed -- no curves split
    return 0;
    
  
    // Get two "new" RefEdges to work with
  RefEdge* new_edges[2];  // two refedges split from input refedge
  RefVertex* vertices[2]; // the to-be other RefVertex on each edge (not new_vtx)
  bool start_vtx[2];      // for the first curve, is new_vtx the start of that curve
  CubitSense senses[2];   // bridge sense for first curve in each refedge
    // For each of two curves adjacent to the first split point
  TBPoint* point = new_vertex->get_point_ptr();
  point->get_parents(curve_bridges);
  assert(curve_bridges.size() == 2);
  curve_bridges.reset();
  for (i = 0; i < 2; i++)
  {
    Curve* curve = dynamic_cast<Curve*>(curve_bridges.get_and_step());
    new_edges[i] = dynamic_cast<RefEdge*>(curve->topology_entity());
    
      // need to create new RefEdge
    if (!new_edges[i])
    {
      if (curve->bridge_sense() == CUBIT_REVERSED)
        curve->reverse_bridge_sense();  // clear bridge sense
      if(is_free_curve)
        new_edges[i] = GeometryQueryTool::instance()->make_free_RefEdge(curve);
      else
        new_edges[i] = GeometryQueryTool::instance()->make_RefEdge(curve);
    }
    
      // get other vertex on refedge, save to use in determining 
      // which curve pairs should merge.
    bridge_list.clean_out();
    curve->get_children(bridge_list);
    bridge_list.move_to(point);
    bridge_list.reset();
    if (bridge_list.get() == point)
      start_vtx[i] = true;
    else if(bridge_list.step_and_get() == point)
      start_vtx[i] = false;
    else
      assert(0);
      
      // populate arrays with necessary info to determine relative 
      // sense of other curves to be merged with this one.
    vertices[i] = dynamic_cast<RefVertex*>(bridge_list.next()->topology_entity());
    senses[i] = curve->bridge_sense();
  }
  
    // Merge remaining curve pairs into the two refedges
  point_bridges.clean_out();
  new_vertex->bridge_manager()->get_bridge_list(point_bridges);
  point_bridges.move_to(point);
  assert(point_bridges.get() == point);
  for (i = point_bridges.size() - 1; i--; )  // for each pair of curves
  {
    TopologyBridge* point_bridge = point_bridges.step_and_get();
    curve_bridges.clean_out();
    point_bridge->get_parents(curve_bridges);
    assert(curve_bridges.size() == 2);
    
    while (curve_bridges.size())  // for each curve in the pair
    {
      Curve* curve = dynamic_cast<Curve*>(curve_bridges.pop());
      bridge_list.clean_out();
      curve->get_children(bridge_list);
      bridge_list.reset();
      bool is_start;
      if (bridge_list.get() == point_bridge)
        is_start = true;
      else if(bridge_list.step_and_get() == point_bridge)
        is_start = false;
      else
        assert(0);  // bad bridge connectivity
      
      RefVertex* vtx = dynamic_cast<RefVertex*>(bridge_list.next()->topology_entity());
      
      for (j = 0; j < 2 && vertices[j] != vtx; j++);
      if (j == 2)
        continue; // can't merge -- something went wrong
      
      if (vertices[0] == vertices[1]) // closed curve, need to do geometric check
      {
        CubitVector center = curve->center_point();
        double d1 = (center - new_edges[0]->center_point()).length_squared();
        double d2 = (center - new_edges[1]->center_point()).length_squared();
        j = d1 < d2 ? 0 : 1;
      }

      RefEdge* edge = new_edges[j];
      if (curve->topology_entity() != edge) // if not already merged ...
      {
        if (curve->bridge_manager())  // if in a different edge, unmerge
          curve->bridge_manager()->remove_bridge(curve);

        curve->get_parents(coedge_bridges);  // need to unmerge coedges too
        while (coedge_bridges.size())
        {
          TopologyBridge* coe_bridge = coedge_bridges.pop();
          if (coe_bridge->owner())  
            coe_bridge->bridge_manager()->remove_bridge(coe_bridge);
        }

        edge->bridge_manager()->add_bridge(curve); // merge
      }
      
        // set relative sense
      bool edge_reversed = (is_start != start_vtx[j]);
      bool curv_reversed = curve->bridge_sense() != senses[j];
      if (edge_reversed != curv_reversed)
        curve->reverse_bridge_sense();
    }
  }
          
    // Rebuild ModelEntity topology
  for( i = faces.size(); i--; )
  { 
    RefFace* face = faces.get_and_step();
    bridge_list.clean_out();
    face->bridge_manager()->get_bridge_list(bridge_list);
    assert(bridge_list.size());
    for (int j = bridge_list.size(); j--; )
    {
      TopologyBridge* bridge = bridge_list.get_and_step();
      Surface* surf = dynamic_cast<Surface*>(bridge);
      assert(!!surf);
      GeometryQueryTool::instance()->make_RefFace(surf);
    }
  }
  
  DLIList<RefEntity*> vtx_list;
  vtx_list.append(new_vertex);
  notify_partition( vtx_list, first_new_edge, second_new_edge, edge_ptr );

    // tweak some stuff so the resulting RefEdges have the same
    // sense as the input RefEdge.

    // Get the edges in the order we want.
    // If curve was not closed, put whichever contains the original start vertex first.
  if (input_start_vtx != input_end_vtx)
    i = new_edges[0]->other_vertex(input_start_vtx) ? 0 : 1;
    //If closed, sense of original edge didn't change so use that to determine the order
  else if (new_edges[0] == edge_ptr)
    i = new_edges[0]->start_vertex() == input_start_vtx ? 0 : 1;
  else
    i = new_edges[1]->start_vertex() == input_start_vtx ? 1 : 0;
  first_new_edge = new_edges[i];
  second_new_edge = new_edges[1-i];
    // Adjust sense of new RefEdge
    // The sense of the input edge should not have changed, so don't mess with
    // it.  Flip the new edge if it doesn't end with the correct vertex.
  if (first_new_edge != edge_ptr && first_new_edge->start_vertex() != input_start_vtx)
  {
    first_new_edge->reverse_tangent();
  }
  if (second_new_edge != edge_ptr && second_new_edge->end_vertex() != input_end_vtx)
  {
    second_new_edge->reverse_tangent();
  }

  PRINT_DEBUG_76("First new curve is %d, second new curve is %d.\n",
                 first_new_edge ? first_new_edge->id() : 0,
                 second_new_edge ? second_new_edge->id() : 0 );
  return new_vertex;
}

RefVertex* PartitionTool::partition( RefEdge* edge_ptr,
                                      const CubitVector& position,
                                      DLIList<RefEdge*>& result_list,
                                      CubitBoolean b )
{
  RefEdge *part1 = 0, *part2 = 0;
  RefVertex* result = partition( edge_ptr, position, part1, part2, b );
  if( result && part1 )
    result_list.append(part1);
  if( result && part2 )
    result_list.append(part2);
  return result;
}
  
//-------------------------------------------------------------------------
// Purpose       : Split a surface
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/20/99
//-------------------------------------------------------------------------
CubitStatus PartitionTool::partition( RefFace* face_ptr,
                                      RefEdge* edge_ptr, 
                                      RefFace*& new_face1,
                                      RefFace*& new_face2,
                                      CubitBoolean )
{
  PRINT_DEBUG_86( "PT:  PartitionTool::partition( RefFace %d, RefEdge %d )\n",
      	      	      	      	  face_ptr->id(), edge_ptr->id() );
  new_face1 = new_face2 = 0;

  DLIList<CubitVector*> segments;
  DLIList<RefEdge*> new_edges;
  GMem edge_facets;
  
  Curve* curve = edge_ptr->get_curve_ptr();
  if( ! curve->get_geometry_query_engine()->get_graphics( curve, &edge_facets ) )
  {
    return CUBIT_FAILURE;
  }

  int j;
    //Use the exact end points to avoid precision problems with faceting.
  CubitVector start_coord = edge_ptr->start_vertex()->coordinates();
  CubitVector end_coord = edge_ptr->end_vertex()->coordinates();      
  segments.append( new CubitVector ( start_coord ) );
  for( j = 1; j < (edge_facets.pointListCount-1); j++ )
  {
    segments.append( new CubitVector( edge_facets.point_list()[j].x,
                                      edge_facets.point_list()[j].y,
                                      edge_facets.point_list()[j].z ) );
  }
  segments.append( new CubitVector ( end_coord ) );

  new_face2 = insert_edge( face_ptr, segments, CUBIT_TRUE, new_edges );
  new_face1 = face_ptr;

  while( segments.size() )
    delete segments.pop();

  return new_face2 ? CUBIT_SUCCESS : CUBIT_FAILURE;
}

RefFace* PartitionTool::insert_edge( RefFace* face_ptr,
                                     DLIList<CubitVector*>& segments,
                                     CubitBoolean       is_meshed,
                                     DLIList<RefEdge*>& new_edges,
                                     int level_of_recursion,
                                     const double *tolerance_length)
{
  int i,j,k;
  const double TOL_SQR = GEOMETRY_RESABS*GEOMETRY_RESABS; 
  DLIList<Body*> bodies;
  face_ptr->bodies(bodies);
  CubitStatus status = CUBIT_SUCCESS;
  
  //Partition edges first
  if (CUBIT_FALSE == is_meshed)
  {
    DLIList<RefEdge*> edge_list;
    face_ptr->ref_edges(edge_list);
    CubitVector pos = *segments.get();
    CubitVector closest;
    RefEdge * the_edge = NULL;
    DLIList <CubitVector *> vect_list;
    DLIList<RefEdge*> temp_list;
    int count = 0;
    for (i = 0 ; count < 2 && i < segments.size(); i++)
    {
      for (j = 0; j < edge_list.size(); j++) 
      { 
        the_edge = edge_list.get_and_step();
        the_edge->closest_point_trimmed(pos, closest);

        if ( (pos - closest).length_squared() > TOL_SQR )
          continue;

        //closest point need to be on the edge
/*
        the_edge->get_graphics(polyline, 0.0);
        int size = polyline.point_list_size();
        GPoint* pts = polyline.point_list();
        CubitVector pt1(pts[0]);
        CubitVector pt2(pts[1]);
        CubitBox edge_box (pt1, pt2);
        if (size > 1)
        {
          for (k = 2; k < size; k ++)
          {
            CubitVector pt3(pts[k]);
            CubitVector pt4(pts[k-1]);
            CubitBox box2(pt4, pt3);
            edge_box |= box2;
          }
        }

        CubitBox position_box(closest);
        if (size == 1 || (size > 1 && !position_box.overlap(0, edge_box)))      
            continue;
*/  
        CubitBox position_box(closest);
        CubitBox edge_box = the_edge->bounding_box(); 
        if (!position_box.overlap(GEOMETRY_RESABS, edge_box))
//        if (!position_box.overlap(0, edge_box))
          continue;

        delete segments.change_to(new CubitVector (closest));
        if((the_edge->start_coordinates() - closest).length_squared() < TOL_SQR 
          ||(the_edge->end_coordinates() - closest).length_squared() < TOL_SQR)
          break;

        vect_list.clean_out();
        vect_list.append(&closest);
        temp_list.clean_out();
        status = partition(the_edge, vect_list, temp_list);

        // explicitly test for status here
        if (!status == CUBIT_SUCCESS)
        {
          PRINT_ERROR("PartitionTool::insert_edge failed at partition...\n");
        }

        // I've seen cases where the resulting partition position is not
        // what we passed in so look at the resulting curves and get
        // the end point that is closest to the point we used for 
        // partitioning the curve.
        double smallest_dist_sq = CUBIT_DBL_MAX;
        RefVertex *closest_vert = NULL;
        for(k=temp_list.size(); k--;)
        {
          RefEdge *current_edge = temp_list.get_and_step();
          RefVertex *start_vert = current_edge->start_vertex();
          RefVertex *end_vert = current_edge->end_vertex();
          double len_sq = (start_vert->coordinates()-closest).length_squared();
          if(len_sq < smallest_dist_sq)
          {
            smallest_dist_sq = len_sq;
            closest_vert = start_vert;
          }
          len_sq = (end_vert->coordinates()-closest).length_squared();
          if(len_sq < smallest_dist_sq)
          {
            smallest_dist_sq = len_sq;
            closest_vert = end_vert;
          }
        }

        if(closest_vert)
          delete segments.change_to(new CubitVector(closest_vert->coordinates()));

        temp_list.remove(the_edge);
        DLIList<RefEntity*> ref_list;
        CAST_LIST_TO_PARENT(temp_list, ref_list);
        notify_partition( ref_list, the_edge, temp_list.get(), the_edge );

        //need to debug later if this happens.
        if (!status == CUBIT_SUCCESS)
        {
          PRINT_ERROR("PartitionTool::insert_edge failed at assert...\n");
        }
        assert (CUBIT_SUCCESS == status);
        edge_list.remove(the_edge);
        count ++;
        break;
      }
    pos = *segments.step_and_get();
    }
  }
  
  if( CubitUndo::get_undo_enabled() )
  {
    DLIList<RefFace*> tmp_face_list(1);
    tmp_face_list.append( face_ptr );
    CubitUndo::save_state_with_cubit_file( tmp_face_list );
  }

  DLIList<TopologyBridge*> bridge_list;
  face_ptr->bridge_manager()->get_bridge_list( bridge_list );
  bridge_list.reset();
  bridge_list.step();
  for( i = bridge_list.size(); i > 1; i-- )
    face_ptr->bridge_manager()->remove_bridge( bridge_list.get_and_step() );
  
  bridge_list.reset();
  DLIList<Surface*> result_surfaces;
  DLIList<Curve*> curves(2), new_curves;;
  DLIList<Surface*> curve_surfs(2);
  int removed_composite_curve = 0;
  for( i = bridge_list.size(); i>0; i--)
  {
    Surface* old_surf = dynamic_cast<Surface*>(bridge_list.get_and_step());

    Surface* new_surf = PartitionEngine::instance().
      insert_curve( old_surf, segments, curves, tolerance_length );
    
    // do it twice in case the segments changed during the first insert_curve
    // because the segments intersect with the old_surf boundary. This is to 
    // ensure that new vertices are generated at those intersecting points. 
    if (!new_surf && level_of_recursion == 0)
    {
      level_of_recursion++;
      RefFace *return_face = insert_edge( face_ptr, segments, CUBIT_FALSE, new_edges, level_of_recursion);
      
      if( CubitUndo::get_undo_enabled() && return_face == NULL ) 
        CubitUndo::remove_last_undo();

      return return_face; 
    }

    if(!new_surf)
    {
      CompositeSurface* cs = dynamic_cast<CompositeSurface*>(old_surf);
      if(cs)
      {
        Surface *tmp_srf = cs->get_surface(0);
        GeometryQueryEngine *gqe = tmp_srf->get_geometry_query_engine();
        double tmp_tol = gqe->get_sme_resabs_tolerance();
        DLIList<Curve*> hidden_curves;
        cs->get_hidden_curves(hidden_curves);
        int k;
        for(k=hidden_curves.size(); k--;)
        {
          Curve *cur_curve = hidden_curves.get_and_step();
          CubitVector start_vec, end_vec;
          double min, max;
          IntersectionTool int_tool;
          int curve_on_polyline = 1;
          cur_curve->get_param_range(min, max);
          cur_curve->position_from_u(min, start_vec);
          if(!int_tool.point_on_polyline(start_vec, segments, &tmp_tol))
            curve_on_polyline = 0;
          if(curve_on_polyline)
          {
            cur_curve->position_from_u(max, end_vec);
            if(!int_tool.point_on_polyline(end_vec, segments, &tmp_tol))
              curve_on_polyline = 0;
          }
          if(curve_on_polyline)
          {
            int num_mid_pts = 2;
            double dt = ((double)(max-min))/((double)(num_mid_pts+1));
            double cur_t = min + dt;
            int n;
            CubitVector cur_vec;
            for(n=0; n<num_mid_pts && curve_on_polyline; ++n)
            {
              cur_curve->position_from_u(cur_t, cur_vec);
              if(!int_tool.point_on_polyline(cur_vec, segments, &tmp_tol))
                curve_on_polyline = 0;
              cur_t += dt;
            }
          }
          if(curve_on_polyline)
          {
            DLIList<Curve*> cur_curve_crvs;
            CompositeCurve *ccurve = dynamic_cast<CompositeCurve*>(cur_curve);
            if(ccurve)
            {
              int f, num_crvs = ccurve->num_curves();
              for(f=0; f<num_crvs; ++f)
                cur_curve_crvs.append(ccurve->get_curve(f));
            }
            else
            {
              cur_curve_crvs.append(cur_curve);
            }
            CompositeEngine::instance().restore_curve(cur_curve);
            ccurve = dynamic_cast<CompositeCurve*>(cur_curve_crvs.get()->owner());
            if(ccurve)
            {
              cur_curve_crvs.clean_out();
              cur_curve_crvs.append(ccurve);
            }
            removed_composite_curve = 1;
            int y, z;
            for(z=cur_curve_crvs.size(); z--;)
            {
              cur_curve = cur_curve_crvs.get_and_step();
              DLIList<TopologyBridge*> pts;
              cur_curve->get_children_virt(pts);
              for(y=pts.size(); y--;)
              {
                DLIList<TopologyBridge*> crvs;
                pts.get()->get_parents_virt(crvs);
                if(crvs.size() == 2)
                {
                  TopologyBridge *other_curve;
                  if(crvs.get() == cur_curve)
                    other_curve = crvs.step_and_get();
                  else
                    other_curve = crvs.get();
                  int p;
                  for(p=curves.size(); p--;)
                  {
                    if(curves.get() == other_curve)
                    {
                      CompositeEngine::instance().remove_point(dynamic_cast<TBPoint*>(pts.get()));
                      curves.remove(curves.get());
                      p = 0;
                    }
                    else
                      curves.step();
                  }
                }
                pts.step();
              }
            }
          }
        }
      }
    }
    if ( !curves.size() && !removed_composite_curve )
    {
      status = CUBIT_FAILURE;
      break;
    }
    for ( int j = curves.size(); j--; )
    {
      curve_surfs.clean_out();
      curves.step_and_get()->surfaces(curve_surfs);
      result_surfaces.merge_unique(curve_surfs);
    }
    new_curves += curves;
    curves.clean_out();

  }
  
  for( i = bodies.size(); i--; )
  {
    Body* body = bodies.get_and_step();
    GeometryQueryTool::instance()->make_Body(body->get_body_sm_ptr());
  }
  
  DLIList<RefFace*> result_faces;
  for( i = result_surfaces.size(); i--; )
    result_faces.append( 
      GeometryQueryTool::instance()->make_RefFace(
        result_surfaces.get_and_step() ) );
  
  MergeTool::instance()->merge_reffaces( result_faces );
  result_faces.clean_out();
  for( i = result_surfaces.size(); i--; )
    result_faces.append_unique( dynamic_cast<RefFace*>(result_surfaces.get_and_step()->topology_entity()) );
  
  result_faces.move_to( face_ptr );
  RefFace* result = result_faces.size() ? result_faces.next() : 0;
  if(!result && removed_composite_curve)
    result = face_ptr;
  for( i = new_curves.size(); i--; )
  {
    Curve* curve = new_curves.get_and_step();
    RefEdge* new_edge = dynamic_cast<RefEdge*>(curve->topology_entity());
    new_edges.append_unique(new_edge);
  }
  
  if ( !removed_composite_curve && !new_edges.size() )
  {
    if( CubitUndo::get_undo_enabled() ) 
      CubitUndo::remove_last_undo();

    return 0;
  }
    
  DLIList<RefEntity*> ref_list;
  CAST_LIST_TO_PARENT(new_edges, ref_list);
  notify_partition( ref_list, face_ptr, result, face_ptr );
  
  if( CubitUndo::get_undo_enabled() && status == CUBIT_FAILURE ) 
    CubitUndo::remove_last_undo();

  return status ? result : 0;
}

DLIList<Surface*>* PartitionTool::group_merged_surfaces(
                                    DLIList<RefFace*>& face_list,
                                    int& num_result_sets )
{
  // make sure all surfaces are merged the same number of times
  int i;
  num_result_sets = 0;
  int bridge_count = face_list.get()->bridge_manager()->number_of_bridges();
  for ( i = face_list.size(); i--; )
    if( face_list.step_and_get()->bridge_manager()->number_of_bridges() != bridge_count )
      return 0;
 
  // get the bridge list (merged Surface*) for the first RefFace
  DLIList<Surface*>* results = new DLIList<Surface*>[bridge_count];
  DLIList<Lump*> surf_owners(bridge_count);
  face_list.reset();
  RefFace* face = face_list.get_and_step();
  DLIList<TopologyBridge*> bridge_list(bridge_count);
  face->bridge_manager()->get_bridge_list(bridge_list);
  
  // get the list of lumps involved in the merge from the first RefFace
  bridge_list.reset();
  DLIList<Lump*> lump_list;
  for ( i = 0; i < bridge_list.size(); i++ )
  {
    Surface* surface = dynamic_cast<Surface*>(bridge_list.get_and_step());
    assert(!!surface);
    lump_list.clean_out();
    surface->lumps(lump_list);
    Lump* lump = 0;
    if ( lump_list.size() > 1 )
    {
      delete [] results;
      return 0;
    }
    else if (lump_list.size())
      lump = lump_list.get();
    
    if( surf_owners.is_in_list( lump ) )
    {
      delete [] results;
      return 0;
    }
    
    surf_owners.append(lump);
    results[i].append(surface);
  }
  
  // add the rest of the Surfaces into separate lists according to
  // the Lump they belong to
  for ( i = 1; i < face_list.size(); i++ )
  {
    face = face_list.get_and_step();
  
    bridge_list.clean_out();
    face->bridge_manager()->get_bridge_list(bridge_list);
    bridge_list.reset();
    
    for ( int j = bridge_list.size(); j--; )
    {
      Surface* surface = dynamic_cast<Surface*>(bridge_list.get_and_step());
      assert(!!surface);
      lump_list.clean_out();
      surface->lumps(lump_list);
      Lump* lump = 0;
      if ( lump_list.size() > 1 )
      {
        delete [] results;
        return 0;
      }
      else if (lump_list.size())
        lump = lump_list.get();

      if ( !surf_owners.move_to(lump) ||
           results[surf_owners.get_index()].size() != i )
      {
        delete [] results;
        return 0;
      }
      
      results[surf_owners.get_index()].append(surface);
    }
  }
  
  // return the number of sets - one set for each Lump that was
  // involved in the merge
  num_result_sets = bridge_count;
  return results;
}


CubitStatus PartitionTool::insert_edge( 
                                     DLIList<RefFace*>& input_faces,
                                     DLIList<CubitVector*>& segments,
                                     DLIList<RefFace*>& output_faces,
                                     DLIList<RefEdge*>& new_edges,
                                     CubitBoolean do_split_curves)
{
  int i, j, w, set_count;
  const double TOL_SQR = GEOMETRY_RESABS*GEOMETRY_RESABS; 
  
  // group surfaces by Lump
  DLIList<Surface*> *surface_sets = group_merged_surfaces( input_faces, set_count );
  if( !surface_sets )
  {
    PRINT_ERROR("Failed to group merged Surfaces in PartitionTool::insert_edge()\n");
    return CUBIT_FAILURE;
  }
  
  input_faces.reset();
  if(do_split_curves)
  {
    for ( w = input_faces.size(); w--; )
    {
      RefFace* face_ptr = input_faces.get_and_step();
      DLIList<RefEdge*> edge_list;
      face_ptr->ref_edges(edge_list);
      CubitVector pos = *segments.get();
      CubitVector closest;
      RefEdge * the_edge = NULL;
      DLIList <CubitVector *> vect_list;
      DLIList<RefEdge*> temp_list;
      int count = 0;
      for (i = 0 ; count < 2 && i < segments.size(); i++)
      {
        for (j = 0; j < edge_list.size(); j++) 
        { 
          the_edge = edge_list.get_and_step();
          the_edge->closest_point(pos, closest);

          if ( (pos - closest).length_squared() > TOL_SQR )
            continue;

        //closest point need to be on the edge
/*
        the_edge->get_graphics(polyline, 0.0);
        int size = polyline.point_list_size();
        GPoint* pts = polyline.point_list();
        CubitVector pt1(pts[0]);
        CubitVector pt2(pts[1]);
        CubitBox edge_box (pt1, pt2);
        if (size > 1)
        {
          for (k = 2; k < size; k ++)
          {
            CubitVector pt3(pts[k]);
            CubitVector pt4(pts[k-1]);
            CubitBox box2(pt4, pt3);
            edge_box |= box2;
          }
        }

        CubitBox position_box(closest);
        if (size == 1 || (size > 1 && !position_box.overlap(0, edge_box)))      
            continue;
*/  
          CubitBox position_box(closest);
          CubitBox edge_box = the_edge->bounding_box(); 
          if (!position_box.overlap(0, edge_box))
            continue;

          delete segments.change_to(new CubitVector (closest));
          if((the_edge->start_coordinates() - closest).length_squared() < TOL_SQR 
            ||(the_edge->end_coordinates() - closest).length_squared() < TOL_SQR)
            break;

          vect_list.clean_out();
          vect_list.append(&closest);
          CubitStatus status = partition(the_edge, vect_list, temp_list);
          //need to debug later if this happens.
          //assert (CUBIT_SUCCESS == status);
          if( CUBIT_FAILURE == status )
             return status;
          edge_list.remove(the_edge);
          count ++;
          break;
        }
        pos = *segments.step_and_get();
      }
    }
  }

  input_faces.reset();

  // unmerge all the surfaces by removing all but the first Lump's
  // surfaces from the bridge manager
  DLIList<TopologyBridge*> bridge_list;
  DLIList<Surface*> surf_list;
  for ( i = input_faces.size(); i--; )
  {
     RefFace* face_ptr = input_faces.get_and_step();
     bridge_list.clean_out();
     face_ptr->bridge_manager()->get_bridge_list( bridge_list );
     CAST_LIST( bridge_list, surf_list, Surface );
     surf_list -= surface_sets[0];
     for ( int j = surf_list.size(); j--; )
      face_ptr->bridge_manager()->remove_bridge( surf_list.get_and_step() );
  }
  
  // perform the partition operation - insert curve - on each Lump's set of surfaces
  DLIList<Surface*> new_surfs, all_new_surfs;
  DLIList<Curve*> new_curves, all_new_curves;
  CubitStatus result = CUBIT_SUCCESS;
  for ( i = 0; i < set_count; i++ )
  {
    new_surfs.clean_out();
    new_curves.clean_out();
    if ( !PartitionEngine::instance().
      insert_curve( surface_sets[i], segments, new_surfs, new_curves ) )
      result = CUBIT_FAILURE;
    all_new_surfs += new_surfs;
    all_new_curves += new_curves;
  }
  delete [] surface_sets;

  // Do not create a new body in the failure case
  if (CUBIT_FAILURE == result)
    return result;
  
  // get the list of all bodies involved in the partition operation
  DLIList<Body*> bodies, all_bodies;
  for ( i = input_faces.size(); i--; )
  {
    bodies.clean_out();
    input_faces.step_and_get()->bodies( bodies );
    all_bodies.merge_unique(bodies);
  }
  
  
  for( i = all_bodies.size(); i--; )
  {
    Body* body = all_bodies.get_and_step();
    BodySM* bodysm = dynamic_cast<BodySM*>(body->bridge_manager()->topology_bridge());
    if( bodysm )
      GeometryQueryTool::instance()->make_Body( bodysm );
  }
  
  // re-merge if surfaces were merged before the partition operation
  if( set_count > 1 )
  {
    DLIList<RefFace*> result_faces;
    for( i = all_new_surfs.size(); i--; )
      result_faces.append( 
        GeometryQueryTool::instance()->make_RefFace(
          all_new_surfs.get_and_step() ) );
    MergeTool::instance()->merge_reffaces( result_faces );
  }
  
  output_faces.clean_out();
  for( i = all_new_surfs.size(); i--; )
    output_faces.append_unique( dynamic_cast<RefFace*>(all_new_surfs.get_and_step()->topology_entity()) );
  output_faces.merge_unique(input_faces);
  
  for( i = all_new_curves.size(); i--; )
  {
    Curve* curve = all_new_curves.get_and_step();
    new_edges.append( dynamic_cast<RefEdge*>(curve->topology_entity()) );
  }
  return result;
}
  
//-------------------------------------------------------------------------
// Purpose       : Create point curves on RefFace
//
// Special Notes : 
//
// Creator       : 
//
// Creation Date :  
//-------------------------------------------------------------------------
CubitStatus PartitionTool::make_point_curves( RefFace* face_ptr,
                                       DLIList<CubitVector> &positions,
                                       DLIList<RefVertex*> &new_vertices)
{
  int i, j;
  DLIList<TopologyBridge*> bridge_list;
  DLIList<Surface*> surface_list;
  DLIList<TBPoint*> new_points;
  CubitStatus rval = CUBIT_SUCCESS;
  
    // Get list of surfaces from RefFace
  face_ptr->bridge_manager()->get_bridge_list( bridge_list );
  CAST_LIST(bridge_list, surface_list, Surface);
  assert(bridge_list.size() == surface_list.size());
  surface_list.reset();
  
  if( CubitUndo::get_undo_enabled() )
  {
    DLIList<RefFace*> tmp_faces(1);
    tmp_faces.append( face_ptr );
    CubitUndo::save_state_with_cubit_file( tmp_faces );
  }

    // For each position in the list
  positions.reset();
  for (i = positions.size(); i--; )
  {
    CubitVector pos = positions.get_and_step();
    RefVertex* result_vtx = 0;
    RefEdge* point_curve = 0;
    CubitSense bridge_sense = CUBIT_UNKNOWN;
    
      // Create point-curve on each surface
    surface_list.last();
    for (j = surface_list.size(); j--; )
    {
      Surface *new_surf, *old_surf = surface_list.step_and_get();
      TBPoint* new_pt = PartitionEngine::instance().
                        insert_point_curve( old_surf, pos, new_surf );
        
      if (!new_pt)
      {
        PRINT_ERROR("Error imprint surface at position (%f,%f,%f)\n",
          pos.x(), pos.y(), pos.z());
        rval = CUBIT_FAILURE;
        continue;
      }
      
        // If surface was replaced with a virtual surface, update
        // the list
      new_points.append(new_pt);
      surface_list.change_to(new_surf);
      
        // Merge each result point into a single RefVertex
      if (!result_vtx)
        result_vtx = GeometryQueryTool::instance()->make_RefVertex(new_pt);
      else if (!new_pt->owner())
        result_vtx->bridge_manager()->add_bridge(new_pt);
        
        // Merge each result point-curve into a single RefEdge.
        // If the result was not a point-curve (input position was
        // on some existing point or curve), just let the merge
        // code clean up, merging/unmerging depending on tolerance.
      bridge_list.clean_out();
      new_pt->get_parents(bridge_list);
      if (bridge_list.size() != 1)
        continue;
      Curve* curve = dynamic_cast<Curve*>(bridge_list.get());
      if (curve->geometry_type() != POINT_CURVE_TYPE)
        continue;
      if (!point_curve)
      {
        bridge_sense = new_surf->bridge_sense();
        if (curve->bridge_sense() != bridge_sense)
          curve->reverse_bridge_sense();
		RefEdge::suppress_edge_length_warning(true);
        point_curve = GeometryQueryTool::instance()->make_RefEdge(curve);
		RefEdge::suppress_edge_length_warning(false);
      }
      else if(!curve->owner())
      {
        point_curve->bridge_manager()->add_bridge(curve);
        CubitSense merge_sense = new_surf->bridge_sense() == bridge_sense ?
                                 CUBIT_FORWARD : CUBIT_REVERSED ;
        if (curve->bridge_sense() != merge_sense)
          curve->reverse_bridge_sense();
      }
    }
  }
  
    // Update DAG topoloy for TB changes
  DLIList<RefFace*> result_faces;
  surface_list.reset();
  for (i = surface_list.size(); i--; )
  {
    RefFace* face = GeometryQueryTool::instance()
      ->make_RefFace(surface_list.get_and_step());
    result_faces.append_unique(face);
  }
  
    // Attempt to re-merge anything that was unmerged
    // during DAG update
  if( result_faces.size() > 1 )
  {
    MergeTool::instance()->merge_reffaces( result_faces );
  }
  
    // Get list of RefVertices that were actually created after
    // all the merging, unmerging, re-merging, etc.
  new_points.reset();
  for (i = new_points.size(); i--; )
    new_vertices.append( 
      dynamic_cast<RefVertex*>(new_points.get_and_step()->topology_entity()));
  
  if (new_vertices.size() > 1)
    new_vertices.uniquify_ordered();  


  if( CubitUndo::get_undo_enabled() && new_points.size() == 0 )
    CubitUndo::remove_last_undo();

  return rval; 
}


//-------------------------------------------------------------------------
// Purpose       : Imprint a point on a surface
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/28/03
//-------------------------------------------------------------------------
RefVertex* PartitionTool::make_point_curve( RefFace* face_ptr,
                                     const CubitVector& position )
{
  DLIList<CubitVector> positions(1);
  DLIList<RefVertex*> results(1);
  positions.append(position);
  
  if( CubitUndo::get_undo_enabled() )
    CubitUndo::save_state_with_cubit_file( face_ptr );

  if (make_point_curves(face_ptr, positions, results) && results.size())
    return results.get();
  else
  {
    if( CubitUndo::get_undo_enabled() )
      CubitUndo::remove_last_undo();
    return 0;
  }
}


//-------------------------------------------------------------------------
// Purpose       : Partition an Edge
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/24/99
//-------------------------------------------------------------------------
/*
CubitStatus PartitionTool::partition( RefEdge* edge_ptr,
                                      DLIList<RefVertex*>& split_points,
                                      DLIList<RefEdge*>& new_edges,
                                      CubitBoolean ignore_mesh )
{
  DLIList<CubitVector*> vect_list;
  split_points.reset();
  for( int i = split_points.size(); i--; )
    vect_list.append( new CubitVector( split_points.get_and_step()->coordinates() ) );
  CubitStatus result = partition( edge_ptr, vect_list, new_edges, ignore_mesh );
  while( vect_list.size() )
    delete vect_list.pop();
  return result;
}
*/
CubitStatus PartitionTool::partition( RefEdge* edge_ptr,
                                      DLIList<CubitVector*>& split_points,
                                      DLIList<RefEdge*>& new_edges,
                                      CubitBoolean ignore_mesh )
{
  DLIList<RefVertex*> new_vertices;
  return partition( edge_ptr, split_points, new_vertices, new_edges, ignore_mesh );
}                                      

CubitStatus PartitionTool::partition( RefEdge* edge_ptr,
                                      DLIList<CubitVector*>& split_points,
                                      DLIList<RefVertex*>&  new_vertices,
                                      DLIList<RefEdge*>& new_edges,
                                      CubitBoolean ignore_mesh )
{
  if( split_points.size() < 1 ) return CUBIT_FAILURE;
  CubitStatus result = CUBIT_SUCCESS;
  RefVertex* new_vtx = 0;

  if( split_points.size() == 1 ) {
    RefEdge *edge1_ptr = 0, *edge2_ptr = 0;  
    new_vtx = partition( edge_ptr, *(split_points.get()), edge1_ptr, edge2_ptr, ignore_mesh );
    if( !new_vtx )
      return CUBIT_FAILURE;
    
    new_vertices.append(new_vtx);  
    if( edge1_ptr ) new_edges.append( edge1_ptr );
    if( edge2_ptr ) new_edges.append( edge2_ptr );
    return CUBIT_SUCCESS;
  }

  int i, j, count = split_points.size();

  CubitVector** vtx_array = new CubitVector*[count];
  double* param_array = new double[count];
  for( i = 0; i < count; i++ )
  {
    vtx_array[i] = split_points.get_and_step();
    param_array[i] = edge_ptr->u_from_position( *(vtx_array[i]) );
  }

  for( i = 0; i < count - 1; i++ )
  {
    int index = i;
    for( j = i + 1; j < count; j++ )
      if( param_array[j] < param_array[index] ) index = j;
    double temp_dbl = param_array[i];
    CubitVector* temp_vtx = vtx_array[i];
    param_array[i] = param_array[index];
    vtx_array[i] = vtx_array[index];
    param_array[index] = temp_dbl;
    vtx_array[index] = temp_vtx;
  }
  delete [] param_array;
  
  RefVertex* end_vtx = edge_ptr->end_vertex();
  for( i = 0; i < count; i++ )
  {
    RefEdge *first_new_edge, *secnd_new_edge;
    new_vtx = partition( edge_ptr, *(vtx_array[i]), first_new_edge, 
                               secnd_new_edge, ignore_mesh );
    if( ! new_vtx ) 
    {
      result = CUBIT_FAILURE;
      continue;
    }
    
    new_vertices.append(new_vtx);
    
    if( ! secnd_new_edge ) secnd_new_edge = edge_ptr;
    if( NULL == secnd_new_edge->other_vertex(end_vtx) )
    {
      edge_ptr = first_new_edge;
      new_edges.append( secnd_new_edge );
    }
    else
    {
      edge_ptr = secnd_new_edge;
      new_edges.append( first_new_edge );
    }
  }
  new_edges.append( edge_ptr );
  delete [] vtx_array;
  
  return result;
}


    

//-------------------------------------------------------------------------
// Purpose       : Combine all possible PartitionEdges
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 08/23/99
//-------------------------------------------------------------------------
CubitStatus PartitionTool::unpartitionAll( DLIList<RefEdge*>& partition_edges,
                                           DLIList<RefEdge*>& restored_edges)
{
  int i;
  DLIList<RefVertex*> vertex_list, tmp_list;
  for ( i = partition_edges.size(); i--; )
  {
    tmp_list.clean_out();
    partition_edges.step_and_get()->ref_vertices( tmp_list );
    vertex_list += tmp_list;
  }
  
  for ( i = vertex_list.size(); i--; )
    vertex_list.step_and_get()->marked(0);
  for ( i = vertex_list.size(); i--; )
  {
    RefVertex* vtx = vertex_list.step_and_get();
    vtx->marked( vtx->marked() + 1 );
  }
  for ( i = vertex_list.size(); i--; )
  {
    RefVertex* vtx = vertex_list.step_and_get();
    
    TBPoint* pt = vtx->get_point_ptr();
    CompositePoint* comp = dynamic_cast<CompositePoint*>(pt);
    if ( comp )
      pt = comp->get_point();
    PartitionPoint* part = dynamic_cast<PartitionPoint*>(pt);
      
    if( vtx->marked() != 2 || !can_remove(vtx) ||
        !part || part->real_point() )
      vertex_list.change_to(0);
    vtx->marked(0);
  }
  vertex_list.remove_all_with_value(0);
  
  if (vertex_list.size() == 0)
    return CUBIT_SUCCESS;
 
  RefEdge* edge;
  std::vector< TDUPtr<RefEdge> > result_list;
  for ( i = vertex_list.size(); i--; )
  {
    edge = CompositeTool::instance()
      ->remove_vertex( vertex_list.get_and_step(), true );
    if ( edge )
    {
      result_list.push_back( edge );
    }
  }
 
  unsigned int j;
  for ( j = 0; j < result_list.size(); j++ )
    if( result_list[j] )
      restored_edges.append_unique( result_list[j].ptr() );
  
  return restored_edges.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
}


//-------------------------------------------------------------------------
// Purpose       : Unpartition faces and cleanup.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/26/99
//-------------------------------------------------------------------------
CubitStatus PartitionTool::unpartitionAll( 
  DLIList<RefFace*>& partition_faces,
  DLIList<RefFace*>& restored_faces )
{
  DLIList<ModelEntity*> query_input, query_output;
  CAST_LIST_TO_PARENT( partition_faces, query_input );
  ModelQueryEngine::instance()->query_model(
    query_input, DagType::co_edge_type(), query_output );
  
  DLIList<CoEdge*> coedge_list;
  DLIList<RefEdge*> refedge_list;
  CAST_LIST( query_output, coedge_list, CoEdge );
  assert( query_output.size() == coedge_list.size() );
  int i;
  for ( i = coedge_list.size(); i--; )
    refedge_list.append( coedge_list.step_and_get()->get_ref_edge_ptr() );
  
  for ( i = refedge_list.size(); i--; )
    refedge_list.step_and_get()->marked(0);
  for ( i = refedge_list.size(); i--; )
  {
    RefEdge* edge = refedge_list.step_and_get();
    edge->marked( edge->marked() + 1 );
  }
  for ( i = refedge_list.size(); i--; )
  {
    RefEdge* edge = refedge_list.step_and_get();
    PartPTCurve* point_curve = dynamic_cast<PartPTCurve*>(edge->get_curve_ptr());
    SegmentedCurve* seg_curve = dynamic_cast<SegmentedCurve*>(edge->get_curve_ptr());
    if((point_curve || (seg_curve && edge->marked() == 2)) &&
               can_remove(edge))
    {
      // Try to remove this partition.
    }
    else
    {
      // Don't try to remove this partition.
      refedge_list.change_to(0);
    }
    edge->marked(0);
  }
  refedge_list.remove_all_with_value(0);
  
  if (refedge_list.size() == 0)
    return CUBIT_SUCCESS;
 
  RefFace* face;
  std::vector< TDUPtr<RefFace> > result_list;
  DLIList<TopologyBridge*> bridge_list;
  for ( i = refedge_list.size(); i--; )
  {
    bool is_partition = false;
    RefEdge* edge = refedge_list.get_and_step();
    edge->bridge_manager()->get_bridge_list( bridge_list );
    while ( bridge_list.size() )
    {
      TopologyBridge* bridge = bridge_list.pop();
      if ( dynamic_cast<PartitionEntity*>(bridge) ) {
        is_partition = true;
        break;
      }
      CompositeCurve* comp = dynamic_cast<CompositeCurve*>(bridge);
      if ( comp ) {
        for ( int j = 0; j < comp->num_curves(); j++ )
          bridge_list.append(comp->get_curve(j));
      }
    }
    if ( !is_partition )
      continue;
    
    face = CompositeTool::instance()->remove_edge( edge, true );
    if ( face )
    {
      result_list.push_back( face );
    }
  }
 
  unsigned int j;
  for ( j = 0; j < result_list.size(); j++ )
    if( result_list[j] )
      restored_faces.append_unique( result_list[j].ptr() );
  
  return restored_faces.size() ? CUBIT_SUCCESS : CUBIT_FAILURE;
}


          

CubitStatus PartitionTool::partition_face_by_curves( RefFace* face_to_split, 
                                      const DLIList<Curve*>& split_curves,
                                      DLIList<RefFace*>& result_set,
                                      CubitBoolean do_split_curves,
                                      DLIList<RefEdge*>* new_edges,
                                      CubitBoolean )
{
  int i, j, numpts;
  CubitStatus status = CUBIT_SUCCESS;
  DLIList<RefEdge*> new_edge_list;
  DLIList<RefFace*> faces[2];
  DLIList<CubitVector*> segments;
  GMem gmem;
  
  if( CubitUndo::get_undo_enabled() )
  {
    DLIList<RefFace*> tmp_ref_face_list(1);
    tmp_ref_face_list.append( face_to_split );
    CubitUndo::save_state_with_cubit_file( tmp_ref_face_list );
  }

  faces[0].append( face_to_split );
  for (i = 0; i < split_curves.size(); i++ )
  {
    Curve* curve = split_curves.next(i);
    GeometryQueryEngine* engine = curve->get_geometry_query_engine();

    if (!engine->get_graphics( curve, &gmem ))
    {
      status = CUBIT_FAILURE;
      continue;
    }

    numpts = gmem.pointListCount;
/*
This code causes a failure in virtual imprint because it compares
against the number of graphics curve facets.  I think this is wrong but
for now I will comment out this code.  In the end, I think that
partitioning should not be done on the faceted model. 
If this hasn't been addressed in 6 months just kill this code. KGM 1/5/06

    // arbitrary # of subdivisions = 20 --KGM
    const double NUM_SEGS = 20.0;

    // if the graphics facets give too many points
    // (which happens if the model is large right now)
    // trim it to some more reasonable number
    if ( numpts > static_cast<int>(NUM_SEGS) )
    {
    
      // this independent of the curve sense
      double t0 = curve->start_param();
      double t1 = curve->end_param();

      double inc = (t1 - t0)/NUM_SEGS;
      double t;
      CubitVector* p1;
      for (t = t0; t < t1;  t+= inc) 
      {
        p1 = new CubitVector;
        CubitStatus status = curve->position_from_u( t, *p1 );
        if (status == CUBIT_SUCCESS)
          segments.append( p1 );
      }

      // make sure the end point is EXACTLY represented
      p1 = new CubitVector;
      CubitStatus status = curve->position_from_u( t1, *p1 );
      if (status == CUBIT_SUCCESS)
      {
        segments.append( p1 );
      }
      else
      {
        PRINT_ERROR("Bad end point evaluation on curve\n");
      }
      
      bool tmp_debug = false;
      if ( tmp_debug )  // KGM
      {
        for ( int i = 0; i < segments.size(); i++ )
        {
          CubitVector* point = segments.get_and_step();
          GfxDebug::draw_point( *point, CUBIT_RED );
        }
        GfxDebug::flush();
      }
    }
    else 
      */
    {
      if (curve->bridge_sense() == CUBIT_REVERSED)
      {
        CubitVector* p1 = new CubitVector;
        curve->position_from_u( curve->end_param(), *p1 );
        segments.append( p1 );
      }
      else
      {
        CubitVector* p1 = new CubitVector;
        curve->position_from_u( curve->start_param(), *p1 );
        segments.append( p1 );
      }
    
      --numpts;  
      for (j = 1; j < numpts; j++)
      {
        const GPoint& p = gmem.point_list()[j];
        segments.append( new CubitVector( p.x, p.y, p.z ) );
      }
      
      if (curve->bridge_sense() == CUBIT_REVERSED)
      {
        CubitVector* p1 = new CubitVector;
        curve->position_from_u( curve->start_param(), *p1 );
        segments.append( p1 );
      }
      else
      {
        CubitVector* p1 = new CubitVector;
        curve->position_from_u( curve->end_param(), *p1 );
        segments.append( p1 );
      }
    }
    
    int index = i % 2;
    status = insert_edge( faces[index], segments, faces[1-index], new_edge_list, do_split_curves );
    if (status)
    {
      if (new_edges)
        *new_edges += new_edge_list;
    }
    else
    {
      faces[1-index] = faces[index];
    }
    
    while (segments.size())
      delete segments.pop();
  }

  result_set = faces[split_curves.size()%2];

  //surface might not have gotton split
  if( CubitUndo::get_undo_enabled() && result_set.size() == 1)
    CubitUndo::remove_last_undo();
  
  return status;
}

//-------------------------------------------------------------------------
// Purpose       : Do multiple partitioning of a RefFace
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 11/23/99
//-------------------------------------------------------------------------
CubitStatus PartitionTool::partition( RefFace* face_to_split, 
                                      const DLIList<RefEdge*>& split_edges,
                                      DLIList<RefFace*>& result_set,
                                      CubitBoolean do_split_curves,
                                      DLIList<RefEdge*>* new_edges,
                                      CubitBoolean )
{
  int i;
  CubitStatus status = CUBIT_SUCCESS;
  
  // get the curves out of the edge list
  DLIList<Curve*> split_curves;
  for (i = 0; i < split_edges.size(); i++ )
  {
    RefEdge* edge = split_edges.next(i);
    split_curves.append( edge->get_curve_ptr() );
  }

  // partition based on the Curves
  status = partition_face_by_curves( face_to_split, split_curves, 
                                     result_set, do_split_curves, new_edges );

  return status;
}

CubitStatus PartitionTool::partition( RefFace* face, 
                                      DLIList<RefEdge*>& edges,
				      RefFace*& new1,
				      RefFace*& new2,
				      CubitBoolean im )
{
  DLIList<RefFace*> result_faces;
  CubitStatus result = partition( face, edges, result_faces, im );
  result_faces.reset();
  new1 = result_faces.get_and_step();
  new2 = result_faces.get_and_step();
  return result;
}



//-------------------------------------------------------------------------
// Purpose       : Partition a volume.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 11/24/99
//-------------------------------------------------------------------------
CubitStatus PartitionTool::partition( RefVolume* volume_ptr,
                                      DLIList<RefFace*>& split_faces,
                                      RefVolume*& first_new_volume,
                                      RefVolume*& second_new_volume,
                                      CubitBoolean  )
{
  int i;
  
  Lump* lump = volume_ptr->get_lump_ptr();
  DLIList<Surface*> surface_list;
  for( i = split_faces.size(); i--; )
    surface_list.append(split_faces.step_and_get()->get_surface_ptr());
  
  DLIList<Lump*> result_list;
  for( i = surface_list.size(); result_list.size() < 2 && i--; )
  {
    Lump* result = PartitionEngine::instance().
      insert_surface( surface_list.step_and_get(), lump );
      
    if(!result)
    {
      RefFace* face = dynamic_cast<RefFace*>(surface_list.get()->topology_entity());
      PRINT_ERROR("Insertion of surface %d into volume %d failed.\n",
        face ? face->id() : 0, volume_ptr->id() );
      continue;
    }
    
    lump = result;      
    result_list.append_unique(lump);
  }
  
  DLIList<Body*> body_list;
  volume_ptr->bodies(body_list);
  assert(body_list.size() == 1);
  Body* body_ptr = body_list.get();
  BodySM* bodysm = body_ptr->get_body_sm_ptr();
  GeometryQueryTool::instance()->make_Body(bodysm);
  
  result_list.reset();
  if( ! result_list.size() )
    return CUBIT_FAILURE;

  first_new_volume = volume_ptr;
  second_new_volume = dynamic_cast<RefVolume*>(result_list.get()->topology_entity());

  return CUBIT_SUCCESS;
}
CubitStatus PartitionTool::partition( RefVolume* volume_ptr,
                                      DLIList<CubitFacet*>& split_faces,
                                      RefVolume*& first_new_volume,
                                      RefVolume*& second_new_volume,
                                      DLIList<RefFace*>& new_surfaces,
                                      CubitBoolean  )
{
  Lump* lump = volume_ptr->get_lump_ptr();
  
  Surface* result = PartitionEngine::instance().
      insert_surface( split_faces, lump );
      
  DLIList<Body*> body_list;
  volume_ptr->bodies(body_list);
  assert(body_list.size() == 1);
  Body* body_ptr = body_list.get();
  BodySM* bodysm = body_ptr->get_body_sm_ptr();
  GeometryQueryTool::instance()->make_Body(bodysm);
  
  if( ! result )
  {
    PRINT_ERROR("Failed to insert faceted surface into volume %d\n",
      volume_ptr->id());
    return CUBIT_FAILURE;
  }

  RefFace* result_face = dynamic_cast<RefFace*>(result->topology_entity());
  assert(!!result_face);
  new_surfaces.append(result_face);
  
  DLIList<RefVolume*> face_vols;
  result_face->ref_volumes(face_vols);
  face_vols.reset();
  first_new_volume = face_vols.get();
  second_new_volume = face_vols.next();
  
if ( DEBUG_FLAG(88) ) {
  GfxDebug::display_all();
  if( PartitionSurface* ps = dynamic_cast<PartitionSurface*>(result) )
    ps->draw_facets(CUBIT_GREEN);
  DLIList<RefEdge*> edges;
  DLIList<RefFace*> faces;
  result_face->ref_edges(edges);
  for ( int i = edges.size(); i--; )
  {
    faces.clean_out();
    RefEdge* edge = edges.get_and_step();
    edge->ref_faces(faces);
    if ( faces.size() > 1 )
      GfxDebug::highlight_ref_edge( edge );
  }
  GfxDebug::flush();
}

  return CUBIT_SUCCESS;
}


//-------------------------------------------------------------------------
// Purpose       : Delete a volume
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/11/03
//-------------------------------------------------------------------------
CubitStatus PartitionTool::destroy_volume_partition( RefVolume* volume )
{
  Body* body = volume->get_body_ptr();
  PartitionLump* lump = dynamic_cast<PartitionLump*>(volume->get_lump_ptr());
  if( !lump ) return CUBIT_FAILURE;
  CubitStatus result = PartitionEngine::instance().destroy_lump(lump);
  TopologyBridge* tb = body->bridge_manager()->topology_bridge();
  GeometryQueryTool::instance()->make_Body(dynamic_cast<BodySM*>(tb));
  return result;
}


CubitStatus PartitionTool::unpartitionAll( DLIList<RefVolume*>& passed_list,
                                           DLIList<RefVolume*>& restored_list )
{
  int i;
  
  if (!passed_list.size())
    return CUBIT_FAILURE;
  
    // Get owning body
  DLIList<Body*> body_list;
  Body* body_ptr = 0;
  for( i = passed_list.size(); i--; )
  {
    RefVolume* vol_ptr = passed_list.step_and_get();
    body_list.clean_out();
    vol_ptr->bodies(body_list);
    assert(body_list.size() == 1);
    if( !body_ptr )
      body_ptr = body_list.get();
    else if( body_ptr != body_list.get() )
    {
      PRINT_ERROR("Invalid volume list passed to PartitionTool::unpartition.\n"
                  "All volumes must belong to the same body.\n");
      return CUBIT_FAILURE;
    }
  }
  
    // Get lump list and real, partitioned lump
  DLIList<PartitionLump*> lump_list;
  TopologyBridge* real_lump = 0;
  for( i = passed_list.size(); i--; )
  {
    Lump* lump = passed_list.step_and_get()->get_lump_ptr();
    PartitionLump* partlump = dynamic_cast<PartitionLump*>(lump);
    if( !partlump )
    {
      PRINT_ERROR("Invalid volume list passed to PartitionTool::unpartition.\n"
                  "Volume %d is not a partition.\n", 
                  passed_list.get()->id());
      return CUBIT_FAILURE;
    }
    if( !real_lump )
      real_lump = partlump->partitioned_entity();
    else if( real_lump != partlump->partitioned_entity() )
    {
      PRINT_ERROR("Invalid volume list passed to PartitionTool::unpartition.\n"
                  "All volumes must be partitions of the same real volume.\n");
      return CUBIT_FAILURE;
    }
    
    lump_list.append(partlump);
  }
  
    // Find all faces volumes have in common
  DLIList<RefFace*> all_faces, vol_faces;
  passed_list.last();
  for( i = passed_list.size(); i--; )
  {
    vol_faces.clean_out();
    passed_list.step_and_get()->ref_faces( vol_faces );
    all_faces += vol_faces;
  }
  for( i = all_faces.size(); i--; )
    all_faces.step_and_get()->marked(-1);
  for( i = all_faces.size(); i--; )
  {
    RefFace* face = all_faces.step_and_get();
    face->marked( face->marked() + 1 );
  }
  for( i = all_faces.size(); i--; )
    if( all_faces.step_and_get()->marked() )
      all_faces.get()->marked(0);
    else
      all_faces.change_to(0);
  all_faces.remove_all_with_value(0);
  
  for( i = all_faces.size(); i--; )
    if ( !can_remove( all_faces.step_and_get() ) )
      all_faces.change_to(0);
  all_faces.remove_all_with_value(0);
  
  
    // Find all partition surfaces
  DLIList<PartitionSurface*> psurf_list;
  DLIList<TopologyBridge*> bridge_list;
  for( i = all_faces.size(); i--; )
  {
    bridge_list.clean_out();
    all_faces.step_and_get()->bridge_manager()->get_bridge_list(bridge_list);
    for( int j = bridge_list.size(); j--; )
    {
      PartitionSurface* psurf = dynamic_cast<PartitionSurface*>(bridge_list.step_and_get());
      if(psurf && psurf->partitioned_entity() == real_lump)
      {
        psurf_list.append(psurf);
      }
    }
  }
    
    // Remove partitioning surfaces
  for( i = psurf_list.size(); i--; )
  {
    PartitionSurface* surf = psurf_list.step_and_get();

    DLIList<Lump*> lumps;
    surf->lumps(lumps);
    assert(lumps.size() > 0  && lumps.size() <= 2);  // should be one or two surfs

    RefVolume* volume1 = NULL;
    RefVolume* volume2 = NULL;

    lumps.reset();
    if(lumps.size() == 2)
    {
      TopologyEntity *topo = lumps.get_and_step()->topology_entity();
      volume1 = CAST_TO(topo, RefVolume);
      assert(volume1 != NULL);

      topo = lumps.get_and_step()->topology_entity();
      volume2 = CAST_TO(topo, RefVolume);
      assert(volume2 != NULL);
    }
    
    //RefFace* face = dynamic_cast<RefFace*>(surf->owner());
    Lump* lump = PartitionEngine::instance().remove_surface(surf);
    if( !lump )
    {
      RefFace* face = dynamic_cast<RefFace*>(surf->topology_entity());
      PRINT_ERROR("Failed to remove surface %d\n", face ? face->id() : 0);
    }
    else
    {
      RefVolume *survivor = dynamic_cast<RefVolume*>(lump->topology_entity());
      RefVolume *dead = lumps.size() < 2 ? 0 : survivor == volume1 ? volume2 : volume1;
      CompositeTool::instance()->update_combined_vols( survivor, dead );
    }
  }
  
    // Find a surviving volume to return
  for ( i = passed_list.size(); i--; )
  {
    if( passed_list.step_and_get()->get_lump_ptr() )
    {
      restored_list.append( passed_list.get() );
    }
  }
    
    // update DAG and return
  BodySM* bodysm = body_ptr->get_body_sm_ptr();
  GeometryQueryTool::instance()->make_Body(bodysm);
  GeometryQueryTool::instance()->cleanout_deactivated_geometry();
  
  return CUBIT_SUCCESS;;
}
  


/************** Methods for debugging output ******************/

void pt_print_bte_list( int debug_flag, DLIList<BasicTopologyEntity*>& edges, 
                        const char* trailing_string  )
{
  if( DEBUG_FLAG( debug_flag ) )
  {
    PRINT_DEBUG(debug_flag)("{ ");
    edges.reset();
    for( int i = edges.size(); i > 1; i-- )
      PRINT_DEBUG(debug_flag)("%d, ",edges.get_and_step()->id() );
    if( edges.size() )
      PRINT_DEBUG(debug_flag)("%d",edges.get_and_step()->id() );
    PRINT_DEBUG(debug_flag)(" }%s",trailing_string ? trailing_string : "");
  }
}
void pt_print_edge_list( int debug_flag, DLIList<RefEdge*>& edges,
                         const char* trailing_string )
{
  if( DEBUG_FLAG( debug_flag ) )
  {
    DLIList<BasicTopologyEntity*> bte_list;
    bte_list.reset();
    CAST_LIST_TO_PARENT( edges, bte_list);
    pt_print_bte_list( debug_flag, bte_list, trailing_string );
  }
}
void pt_print_face_list( int debug_flag, DLIList<RefFace*>& faces,
                         const char* trailing_string )
{
  if( DEBUG_FLAG( debug_flag ) )
  {
    DLIList<BasicTopologyEntity*> bte_list;
    bte_list.reset();
    CAST_LIST_TO_PARENT( faces, bte_list);
    pt_print_bte_list( debug_flag, bte_list, trailing_string );
  }
}
void pt_print_loop( int debug_flag, Loop* loop_ptr, 
                    const char* trailing_string )
{
  if( DEBUG_FLAG( debug_flag ) )
  {
    DLIList<RefEdge*> edges;
    loop_ptr->ordered_ref_edges( edges );
    pt_print_edge_list( debug_flag, edges, trailing_string );
  }
}
void pt_print_shell( int debug_flag, Shell* shell_ptr, 
                     const char* trailing_string )
{
  if( DEBUG_FLAG( debug_flag ) )
  {
    DLIList<RefFace*> faces;
    shell_ptr->ref_faces( faces );
    pt_print_face_list( debug_flag, faces, trailing_string );
  }
}

void PartitionTool::notify_partition( 
  DLIList<RefEntity*> & /*partitioning_entities*/,
  BasicTopologyEntity *first_partitioned_entity,
  BasicTopologyEntity *second_partitioned_entity,
  BasicTopologyEntity *old_entity )
{
  int i;

    // are the first and second entities unique?
  int first_unique = 
    first_partitioned_entity != old_entity &&
    first_partitioned_entity != NULL;
  
  int second_unique = 
    second_partitioned_entity != first_partitioned_entity &&
    second_partitioned_entity != old_entity &&
    second_partitioned_entity != NULL;
  
      // add unique entities to graphics
      // attach names to entities.
  if ( first_unique )
  {
    RefEntityName::instance()
      ->copy_refentity_names( old_entity, first_partitioned_entity );
  }
  if ( second_unique )
  {
    RefEntityName::instance()
      ->copy_refentity_names( old_entity, second_partitioned_entity );
  }
  
    // notify observers of this 
  AppUtil::instance()->send_event(old_entity, GEOMETRY_TOPOLOGY_MODIFIED );
  
    // notify containing entities
  BasicTopologyEntity* entity_ptr = old_entity;
  if( first_unique ) 
    entity_ptr = first_partitioned_entity;
  else if( second_unique )
    entity_ptr = second_partitioned_entity;
    
  DLIList<RefFace*> face_list;
  RefFace *ref_face;
  DLIList<RefVolume*> vol_list;
  RefVolume *ref_volume;
  DLIList<Body*> body_list;
  Body *body;
  int has_parents = CUBIT_TRUE;
  int dim = entity_ptr->dimension();
  switch ( dim )
  {
    case 0:
        // ?
        // no break
    case 1:
        // does old_entity still have pointers to the ref_faces?
      entity_ptr->ref_faces( face_list );
      if ( face_list.size() == 0 && dim == 1 )
        has_parents = CUBIT_FALSE;
      for( i = face_list.size(); i > 0; i-- )
      {
        ref_face =  face_list.get_and_step();
        AppUtil::instance()->send_event( ref_face, TOPOLOGY_MODIFIED );
      }

        // no break      
    case 2:
      entity_ptr->ref_volumes( vol_list );
      if ( vol_list.size() == 0 && dim == 2 )
        has_parents = CUBIT_FALSE;
      for( i = vol_list.size(); i > 0; i-- )
      {
        ref_volume = vol_list.get_and_step();
        AppUtil::instance()->send_event( ref_volume, TOPOLOGY_MODIFIED );
      }
      
        // no break      
    case 3:
      entity_ptr->bodies( body_list );
      if ( body_list.size() == 0 && dim == 3 )
        has_parents = CUBIT_FALSE;
      for( i = body_list.size(); i > 0; i-- )
      {
        body = body_list.get_and_step();
        AppUtil::instance()->send_event( body, TOPOLOGY_MODIFIED );
      }
      break;
    default:
      break;
  }
  
  if( !has_parents )
  {
    // TODO -- should call GMT::make_free_RefEdge where appropriate or
    // GQE could be changed to recoginize it is making a free RefEdge.
    // When that is done, we no longer need these calls to the static observers.
    if( second_unique )
    {
      AppUtil::instance()->send_event( second_partitioned_entity, CubitEvent(FREE_REF_ENTITY_GENERATED));
    }
    if( first_unique )
    {
      AppUtil::instance()->send_event( first_partitioned_entity, CubitEvent(FREE_REF_ENTITY_GENERATED));
    }
  }
  
}

//-------------------------------------------------------------------------
// Purpose       : These are overloaded by PartitionToolMesh
//                 to prevent un-partitioning where the mesh cannot
//                 be updated.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/12/03
//-------------------------------------------------------------------------
CubitStatus PartitionTool::can_remove( RefVertex* ) { return CUBIT_SUCCESS; }
CubitStatus PartitionTool::can_remove( RefEdge* ) { return CUBIT_SUCCESS; }
CubitStatus PartitionTool::can_remove( RefFace* ) { return CUBIT_SUCCESS; }
