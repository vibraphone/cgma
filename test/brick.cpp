#include "GeometryModifyTool.hpp"
#include "GeometryQueryTool.hpp"
#include "ModelQueryEngine.hpp"
#include "RefVertex.hpp"
#include "RefEdge.hpp"
#include "RefFace.hpp"
#include "Body.hpp"
#include "AppUtil.hpp"
#include "Loop.hpp"
#include "CoEdge.hpp"
#include "InitCGMA.hpp"

#include <assert.h>
#include <algorithm>


/* Check that CGM's presentation of the modeler topology
   is correct by creating a simple cube and checking that
   everything is correctly connected.
 */


void check_valid_edge( RefEdge* edge );
void check_valid_face( RefFace* face );
void check_valid_loop( Loop* loop );

int main( int argc, char* argv[] )
{
    // Start up CGM
  CubitStatus result = InitCGMA::initialize_cgma();
  if (CUBIT_SUCCESS != result) return 1;

    // Create a brick
  Body* brick = GeometryModifyTool::instance()->brick( 2, 2, 2 );
  assert(brick != 0);
  
    // Get all child RefEntities
  DLIList<ModelEntity*> query_results;
  DLIList<RefFace*> surfaces;
  DLIList<RefEdge*> curves;
  DLIList<RefVertex*> points;
  
  ModelQueryEngine::instance()->query_model( *brick, DagType::ref_face_type(), query_results );
  CAST_LIST( query_results, surfaces, RefFace );
  query_results.clean_out();
  ModelQueryEngine::instance()->query_model( *brick, DagType::ref_edge_type(), query_results );
  CAST_LIST( query_results, curves, RefEdge );
  query_results.clean_out();
  ModelQueryEngine::instance()->query_model( *brick, DagType::ref_vertex_type(), query_results );
  CAST_LIST( query_results, points, RefVertex );
  query_results.clean_out();
  
    // Check expected number of child RefEntities
  assert( surfaces.size() == 6 );
  assert( curves.size() == 12 );
  assert( points.size() == 8 );
  
    // Cube should be 2x2x2 and centered at origin, so each vertex
    // should have the coords (+/-1, +/-1, +/-1).  Check that this
    // is the case and put vertex pointers in an array in the canonical
    // order for hexahedra.
  RefVertex* corners[8] = { 0, 0, 0, 0, 0, 0, 0, 0};
  for (int i = 0; i < points.size(); ++i) {
    RefVertex* vtx = points.get_and_step();
    CubitVector pos = vtx->coordinates();
    // all coordinate values are plus or minus one
    assert( fabs(fabs(pos.x()) - 1) < 1e-12 );
    assert( fabs(fabs(pos.y()) - 1) < 1e-12 );
    assert( fabs(fabs(pos.z()) - 1) < 1e-12 );
    int idx;
    if (pos.x() < 0) 
      idx = pos.y() < 0 ? 0 : 3;
    else
      idx = pos.y() < 0 ? 1 : 2;
    if (pos.z() > 0)
      idx += 4;
    // No duplicate vertices allowed
    assert( corners[idx] == 0 );
    corners[idx] = vtx;
  }
  
    // Now figure out which of the 12 edges each RefEdge is by
    // checking the position of each attached vertex in the
    // above 'corners' array,  Put edge pointers in an array in
    // the canonical edge order for a heaxahedron.
  RefEdge* edges[12] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  for (int i = 0; i < curves.size(); ++i) {
    RefEdge* edge = curves.get_and_step();
    int start_i = std::find( corners, corners + 8, edge->start_vertex() ) - corners;
    int end_i = std::find( corners, corners + 8, edge->end_vertex() ) - corners;
    assert( start_i < 8 );
    assert( end_i < 8 );
    int idx;
      // If edge is lateral edge from bottom face to top
    if (start_i < 4 && end_i >= 4) {
      idx = start_i + 4;
      assert( end_i == idx );
    }
      // If edge is lateral edge from top face to bottom
    else if (start_i >= 4 && end_i < 4) {
      idx = end_i + 4;
      assert( start_i == idx );
    }
      // If edge is part of bottom face
    else if (start_i < 4) {
      if ((start_i + 1)%4 == end_i) 
        idx = start_i;
      else if ((end_i + 1)%4 == start_i)
        idx = end_i;
      else
        assert( false );
    }
      // If edge is part of top face
    else {
      if ((start_i - 3) % 4 + 4 == end_i)
        idx = start_i + 4;
      else if ((end_i - 3) % 4 + 4 == start_i)
        idx = end_i + 4;
      else 
        assert( false );
    }
    // Check no duplicate edges
    assert( 0 == edges[idx] );
    edges[idx] = edge;    
  }
  
  // Face connectivity, specified by edge numbers
  int std_faces[6][4] = { { 0, 5, 8, 4 },
                          { 1, 6, 9, 5 },
                          { 2, 7, 10, 6 },
                          { 3, 4, 11, 7 },
                          { 3, 2, 1, 0 },
                          { 8, 9, 10,11 } };
  RefFace* faces[6] = { 0, 0, 0, 0, 0, 0 };
  DLIList< DLIList<RefEdge*>* > loops;
  for (int i = 0; i < surfaces.size(); ++i) {
    // get ref edges
    RefFace* face = surfaces.get_and_step();
    face->ref_edge_loops( loops );
    assert( loops.size() == 1 );
    DLIList<RefEdge*> loop(*loops.get());
    delete loops.get();
    loops.clean_out();
    
    // match loop of edges to one of the expected faces
    assert( loop.size() == 4 );
    int k;
    for (k = 0; k < 6; ++k) {
      int j;
      for (j = 0; j < 4; ++j) 
        if (!loop.is_in_list(edges[std_faces[k][j]]))
          break;
      if (j == 4)
        break;
    }
    assert( k < 6 ); // fail of didn't match any expected face
    assert(0 == faces[k]); // no duplicates
    faces[k] = face;
  }
  
  for (int i = 0; i < 12; ++i)
    check_valid_edge( edges[i] );
  for (int i = 0; i < 6; ++i) 
    check_valid_face( faces[i] );
  for (int i = 0; i < 6; ++i) {
    assert(faces[i]->number_of_Loops() == 1);
    check_valid_loop( dynamic_cast<Loop*>(faces[i]->get_first_grouping_entity_ptr()) );
  }

  return 0;
}

void check_valid_edge( RefEdge* edge )
{
    // Check that start and end vertices are in the correct order
    // by comparing vertex coordinates to edge position from start
    // and end parameters.
  double start_u, end_u;
  edge->get_param_range( start_u, end_u );
  CubitStatus s;
  CubitVector start_p, end_p;
  s = edge->position_from_u( start_u, start_p );
  assert(CUBIT_SUCCESS == s);
  assert( (start_p - edge->start_vertex()->coordinates()).length() < 1e-6 );
  s = edge->position_from_u( end_u, end_p );
  assert(CUBIT_SUCCESS == s);
  assert( (end_p - edge->end_vertex()->coordinates()).length() < 1e-6 );
  
    // Check reverse query from vertex to edge brings us back to the
    // input edge.
  DLIList<RefEntity*> ents;
  edge->start_vertex()->get_parent_ref_entities( ents );
  assert( ents.is_in_list( edge ) );
  ents.clean_out();
  edge->end_vertex()->get_parent_ref_entities( ents );
  assert( ents.is_in_list( edge ) );
}

void check_valid_face( RefFace* face )
{
  DLIList<RefEntity*> edges, parents;
  face->get_child_ref_entities( edges );
  
  double u_min, u_max, v_min, v_max, u, v;
  face->get_param_range_U( u_min, u_max );
  face->get_param_range_V( v_min, v_max );
  
    // Check that each edge is on the face
  for (int i = 0; i < edges.size(); ++i) {
    RefEdge* edge = dynamic_cast<RefEdge*>(edges.get_and_step());
    assert(!!edge);
    
      // Get three points to sample at
    CubitVector p[3];
    p[0] = edge->start_vertex()->coordinates();
    p[1] = edge->center_point();
    p[2] = edge->end_vertex()->coordinates();
    
      // For each sample point
    for (int j = 0; j < 3; ++j) {
        // Check that point is on geometric surface
      CubitVector close(p[j]);
      face->move_to_surface( close );
      assert( (close - p[j]).length() < 1e-6 );
        // Check that point is within UV bounds
      CubitStatus s = face->u_v_from_position( p[j], u, v );
      assert( CUBIT_SUCCESS == s );
      assert( (u - u_min) > -1e-6 );
      assert( (u_max - u) > -1e-6 );
      assert( (v - v_min) > -1e-6 );
      assert( (v_max - v) > -1e-6 );
    }
  }
  
    // Check reverse query from edge to face brings us back to the
    // input face.
  for (int i = 0; i < edges.size(); ++i) {
    parents.clean_out();
    RefEntity* edge = edges.get_and_step();
    edge->get_parent_ref_entities( parents );
    assert( parents.is_in_list( face ) );
  }
}

// Check that coedges in loop chare expected vertices.
void check_valid_loop( Loop* loop )
{
  DLIList<CoEdge*> coedges;
  loop->ordered_co_edges( coedges );
  CoEdge* prev = coedges.get_and_step();
  for (int i = 0; i < coedges.size(); ++i) {
    CoEdge* curr = coedges.get_and_step();
    RefVertex* common_1 = prev->get_sense() == CUBIT_REVERSED ?
                          prev->get_ref_edge_ptr()->start_vertex() :
                          prev->get_ref_edge_ptr()->end_vertex() ;
    RefVertex* common_2 = curr->get_sense() == CUBIT_REVERSED ?
                          curr->get_ref_edge_ptr()->end_vertex() :
                          curr->get_ref_edge_ptr()->start_vertex();
    assert( common_1 == common_2 );
    prev = curr;
  }
}


    
    
    
