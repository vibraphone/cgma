//-------------------------------------------------------------------------
// Filename      : ChordalAxis.cpp
//
// Purpose       :
//
// Special Notes : Delaunay Triangulation Using a Uniform Grid
//
// Creator       : Jit Ken Tan
//
// Creation Date : 01/20/03
//-------------------------------------------------------------------------

#include "ChordalAxis.hpp"
#include "DLIList.hpp"
#include "CubitPoint.hpp"
#include "CubitPointData.hpp"
#include "CubitFacet.hpp"
#include "CubitFacetData.hpp"
#include "CubitFacetEdge.hpp"
#include "CubitFacetEdgeData.hpp"
#include "CubitMessage.hpp"
#include "CpuTimer.hpp"
#include "CubitVector.hpp"
#include "assert.h"
#include "FacetorTool.hpp"
#include "TDChordal.hpp"
#include "ParamTool.hpp"
#include "PlanarParamTool.hpp"

//-------------------------------------------------------------------------
// Purpose       : Constructor
//
// Special Notes : 
//
// Creator       : Jit Ken Tan
//
// Creation Date : 01/20/03
//-------------------------------------------------------------------------
ChordalAxis::ChordalAxis()
{
  startNodes = NULL;
  endNodes = NULL;
}

//-------------------------------------------------------------------------
// Purpose       : Destructor
//
// Special Notes : 
//
// Creator       : Jit Ken Tan
//
// Creation Date : 01/20/03
//-------------------------------------------------------------------------
ChordalAxis::~ChordalAxis()
{
  if(startNodes)
    delete [] startNodes;
  if(endNodes)
    delete [] endNodes;

}

//-------------------------------------------------------------------------
// Purpose       : Create a delaunay constrained triangulation
//
// Special Notes : 
//
// Creator       : Jit Ken Tan
//
// Creation Date : 01/20/03
//-------------------------------------------------------------------------

CubitStatus ChordalAxis::create_triangulation(DLIList <CubitPoint*> boundary_loop_list){ 
		     
  int ii;
  //int num_boundary_edges = 0;
  int dummy_variable = -1;
  int surf_id = 1;
  CubitStatus status = CUBIT_SUCCESS;
  DLIList<CubitPoint*> bounding_nodes;
  ParamTool *p_tool = NULL;
 
  //computing number of boundary edges
  numBoundaryEdges = boundary_loop_list.size();

  //finding the start and end nodes of edges
  startNodes = new CubitPoint* [numBoundaryEdges];
  endNodes = new CubitPoint* [numBoundaryEdges];

  for(ii = boundary_loop_list.size(); ii > 0; ii--){
    startNodes[ii-1] = boundary_loop_list.get_and_step();
    endNodes[ii-1] = boundary_loop_list.get();
  }

   //obtaining a list of boundary_nodes
  for(ii = 0; ii < numBoundaryEdges; ii++){
    bounding_nodes.append(startNodes[ii]);
  }

  FacetorTool<int,CubitFacet, CubitFacetEdge, CubitPoint, CubitFacetData, CubitPointData, int> facetor = FacetorTool<int,CubitFacet, CubitFacetEdge, CubitPoint, CubitFacetData, CubitPointData, int> (&surf_id, bounding_nodes, startNodes, endNodes, numBoundaryEdges, &dummy_variable, p_tool);
  
  status = facetor.mesh_surfwoIP(triangulation);

  return status;
}

//-------------------------------------------------------------------------
// Purpose       : Obtain the TDChordal Tool Data
//
// Special Notes : 
//
// Creator       : Jit Ken Tan
//
// Creation Date : 01/20/03
//-------------------------------------------------------------------------
TDChordal * ChordalAxis::get_tool_data(CubitFacet *curr_facet){

  ToolData *td = curr_facet-> get_TD( TDChordal::is_chordal );
  TDChordal *td_chordal = dynamic_cast<TDChordal *> (td);
  return td_chordal;
}

//-------------------------------------------------------------------------
// Purpose       : Perform the entire chordal axis extraction
//
// Special Notes : 
//
// Creator       : Jit Ken Tan
//
// Creation Date : 01/20/03
//-------------------------------------------------------------------------
CubitStatus ChordalAxis::chordal_axis_transform(DLIList <CubitPoint*> boundary_loop_list,DLIList <CubitFacet *> &facet_list, DLIList <CubitFacetEdge *> &axis){
  CubitStatus status = CUBIT_SUCCESS;

  // initialize alpha
  pruningAlpha = 0;
  
  status = create_triangulation(boundary_loop_list);
  if(status != CUBIT_SUCCESS){
    PRINT_ERROR("Delaunay Triangulation failed\n");
    return status;
  }

  facet_list = triangulation;

   status = flagging_boundary();
  if(status != CUBIT_SUCCESS){
    PRINT_ERROR("Flagging of boundary failed\n");
    return status;
  }

  status = process_triangles();
  if(status != CUBIT_SUCCESS){
    PRINT_ERROR("Processing of triangles failed\n");
    return status;
  }

  if(pruningAlpha > 0){
    status = pruning();
    if(status != CUBIT_SUCCESS){
      PRINT_ERROR("Pruning failed\n");
      return status;
    }
  }

  status =  axis_extraction(axis);
  if(status != CUBIT_SUCCESS){
    PRINT_ERROR("Axis extraction failed\n");
    return status;
  }

  axis.reset();

  
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Consolidate the junction triangles in our mesh 
//
// Special Notes : If pruning is implemented, the pruned triangles are ignored //                 in determining junction tris
//
// Creator       : Jit Ken Tan
//
// Creation Date : 01/20/03
//-------------------------------------------------------------------------
CubitStatus ChordalAxis::junction_tris_collection(){
  junctionTris.clean_out();
  int ii;
  for (ii = triangulation.size(); ii > 0; ii--){
    CubitFacet *current_facet = triangulation.get_and_step();
    TDChordal *td_chordal = get_tool_data(current_facet);
    TriType tri_class = td_chordal->get_tritype();

    if(tri_class == JUNCTION)
      junctionTris.append(current_facet);
  }
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Perform pruning of the chordal axis
//
// Special Notes : 
//
// Creator       : Jit Ken Tan
//
// Creation Date : 01/20/03
//-------------------------------------------------------------------------
CubitStatus ChordalAxis::pruning(){
  int ii;
  CubitStatus status = CUBIT_SUCCESS;
  bool pruning_performed = FALSE;

  junction_tris_collection();
  
  for(ii = 0; ii < junctionTris.size(); ii++){
    CubitFacet *curr_facet = junctionTris.get_and_step();
    TDChordal *td_chordal = get_tool_data(curr_facet);
    
    //error, return failure
    if(td_chordal == NULL){
      PRINT_ERROR("Tool Data of facet doesn't exist\n");
      return CUBIT_FAILURE;
    }

    //triangle is already pruned
    if(td_chordal->get_tritype() == DISCARDED)
      continue;
    
    int jj;
    for(jj = 0; jj < 3; jj++){
      mark_visited(curr_facet);
      status = junction_triangle_pruning(curr_facet, jj, pruning_performed);
      if(status == CUBIT_FAILURE){
	PRINT_ERROR("Junction triangle pruning function failed\n");
	return status;
      }
    }

    unmark_visited(curr_facet);

  }

  if(pruning_performed == TRUE){
    status = pruning();
    if(status == CUBIT_FAILURE)
      return status;
  }
  
  return CUBIT_SUCCESS;
  
}

//-------------------------------------------------------------------------
// Purpose       : Pruning commences with a junction triangle
//
// Special Notes : 
//
// Creator       : Jit Ken Tan
//
// Creation Date : 01/20/03
//-------------------------------------------------------------------------
CubitStatus ChordalAxis::junction_triangle_pruning(CubitFacet *curr_facet, int edge_index, bool &pruning_performed){
  
  CubitStatus status = CUBIT_SUCCESS;
  CubitFacetEdge *curr_edge = curr_facet->edge(edge_index);
  CubitPoint *start_node = curr_edge->start_node();
  CubitPoint *end_node = curr_edge->end_node();
  CubitVector start_vector = start_node->coordinates();
  CubitVector end_vector = end_node->coordinates();
  CubitVector edge_vector = end_vector - start_vector;
  double edge_length = edge_vector.length();
  CubitVector unit_edge_vector = edge_vector/edge_length;

  CubitFacet *adj_triangle = curr_facet->shared_facet(start_node, end_node);

  if(adj_triangle == NULL){
    PRINT_INFO("Error in Collection of Junction Triangles\n");
    return CUBIT_FAILURE;
  }
  
  DLIList <CubitFacet *> pruned_triangles;
  double pruning_distance = pruningAlpha*edge_length;
  CubitStatus prune = CUBIT_SUCCESS;

  status = search_triangles_to_prune(adj_triangle,curr_edge, start_node, unit_edge_vector, pruning_distance, pruned_triangles, prune );
  if(status == CUBIT_FAILURE){
    PRINT_ERROR("The search for triangles to Prune has failed\n");
    return status;
  }  

  // We pruned the triangles
  if(prune == CUBIT_SUCCESS){
    pruning_performed = TRUE;
    flag_boundary(curr_facet, start_node,end_node);
    determine_tritype(curr_facet);
    int kk;
    for(kk = pruned_triangles.size(); kk > 0; kk--){
      CubitFacet *pruning_triangle = pruned_triangles.get_and_step();
      status = mark_tri_discarded(pruning_triangle);
      if(status == CUBIT_FAILURE){
	PRINT_ERROR("Marking tris as discarded has failed\n");
	return CUBIT_FAILURE;
      }
    }
    
  }

  // unmarked the visited triangles
  int kk;
  for(kk = pruned_triangles.size(); kk > 0; kk--){
    CubitFacet *visited_triangle = pruned_triangles.get_and_step();
    status = unmark_visited(visited_triangle);
    if(status == CUBIT_FAILURE){
      PRINT_ERROR("Marking tris as unvisited has failed\n");
      return CUBIT_FAILURE;
    }
  }

  return CUBIT_SUCCESS;

}

//-------------------------------------------------------------------------
// Purpose       : Marked a triangle as discarded
//
// Special Notes : 
//
// Creator       : Jit Ken Tan
//
// Creation Date : 01/20/03
//-------------------------------------------------------------------------
CubitStatus ChordalAxis::mark_tri_discarded(CubitFacet *curr_facet){

  TDChordal *td_chordal = get_tool_data(curr_facet);
  
  //error, return failure
  if(td_chordal == NULL)
    return CUBIT_FAILURE;
  
  td_chordal->set_tritype(DISCARDED);
 
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Searching for triangles to prune
//
// Special Notes : 
//
// Creator       : Jit Ken Tan
//
// Creation Date : 01/20/03
//-------------------------------------------------------------------------
CubitStatus ChordalAxis::search_triangles_to_prune(CubitFacet *curr_facet, CubitFacetEdge *curr_edge, CubitPoint *start_node, CubitVector unit_edge_vector, double pruning_distance, DLIList <CubitFacet *> &pruned_triangles, CubitStatus &prune ){
  
  CubitStatus status = CUBIT_SUCCESS;
  
  int ii, other_index;
  other_index = curr_facet->other_index(curr_edge->start_node(),curr_edge->end_node());


  if(other_index == -1){
    PRINT_ERROR("Can't find the other_index to the current Facet\n");
    return CUBIT_FAILURE;
  }

  //Computing perpendicular distance of other point from edge
  
  CubitVector other_point_vec = curr_facet->point(other_index)->coordinates();
  CubitVector edge_to_point_vec = other_point_vec - start_node->coordinates();
  double edge_to_point_dist = (edge_to_point_vec * unit_edge_vector).length();
  
  //make edge_to_point distance positive
  if(edge_to_point_dist < 0){
    edge_to_point_dist = -edge_to_point_dist;
  }
 
  if(edge_to_point_dist > pruning_distance){
    prune = CUBIT_FAILURE;
    return CUBIT_SUCCESS;
  }

  status = mark_visited(curr_facet);
  if(status == CUBIT_FAILURE){
    PRINT_ERROR("Marking facet as unvisited has failed\n");
    return status;
  }

  pruned_triangles.append(curr_facet);

  for(ii = 0; ii < 3 && prune == CUBIT_SUCCESS; ii++){
    CubitFacetEdge *curr_edge = curr_facet->edge(ii);
    CubitPoint *curr_edge_start = curr_edge->start_node();
    CubitPoint *curr_edge_end = curr_edge->end_node();
    CubitFacet *adj_facet = curr_facet->shared_facet(curr_edge_start, curr_edge_end);
    
    // Check if there is an adjacent triangle
    if(adj_facet == NULL )
      continue;
    
    // Check if adjacent triangle is visited  
    bool visited = TRUE;
    status = get_visited(adj_facet, visited);
    if(status == CUBIT_FAILURE){
      PRINT_ERROR("getting visited status of triangle has failed\n");
      return status;
    }
    if(visited == TRUE)
      continue;

    status = search_triangles_to_prune(adj_facet,curr_edge, start_node, unit_edge_vector, pruning_distance, pruned_triangles, prune);
    if(status == CUBIT_FAILURE){
      PRINT_ERROR("Search for triangles has failed\n");
      return CUBIT_FAILURE;
    }
  }

  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Flagging all the triangles with boundary edges
//
// Special Notes : 
//
// Creator       : Jit Ken Tan
//
// Creation Date : 01/20/03
//-------------------------------------------------------------------------
CubitStatus ChordalAxis::flagging_boundary(){
  int ii;
  for (ii = 0; ii < numBoundaryEdges; ii++){
    CubitFacet *f1, *f2;
    startNodes[ii] -> shared_facets(endNodes[ii], f1, f2);
    if(!f1 || f2)
      return CUBIT_FAILURE;
    flag_boundary(f1, startNodes[ii], endNodes[ii]);
  }
  
  return CUBIT_SUCCESS;

}

//-------------------------------------------------------------------------
// Purpose       : Mark Triangle as visited
//
// Special Notes : 
//
// Creator       : Jit Ken Tan
//
// Creation Date : 01/20/03
//-------------------------------------------------------------------------
CubitStatus ChordalAxis::mark_visited(CubitFacet *curr_facet){
  ToolData *td = curr_facet->get_TD( TDChordal::is_chordal );
  TDChordal *td_chordal = dynamic_cast<TDChordal *> (td);

  
  if (td_chordal == NULL) {
    return CUBIT_FAILURE;
  }

  td_chordal->mark_visited();
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Mark Triangle as un-visited
//
// Special Notes : 
//
// Creator       : Jit Ken Tan
//
// Creation Date : 01/20/03
//-------------------------------------------------------------------------
CubitStatus ChordalAxis::unmark_visited(CubitFacet *curr_facet){
  ToolData *td = curr_facet->get_TD( TDChordal::is_chordal );
  TDChordal *td_chordal = dynamic_cast<TDChordal *> (td);

  
  if (td_chordal == NULL) {
    return CUBIT_FAILURE;
  }
  
  td_chordal->unmark_visited();
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Check if the triangle has been visited
//
// Special Notes : 
//
// Creator       : Jit Ken Tan
//
// Creation Date : 01/20/03
//-------------------------------------------------------------------------
CubitStatus ChordalAxis::get_visited(CubitFacet *curr_facet, bool &visited){
  ToolData *td = curr_facet->get_TD( TDChordal::is_chordal );
  TDChordal *td_chordal = dynamic_cast<TDChordal *> (td);

  
  if (td_chordal == NULL) {
    return CUBIT_FAILURE;
  }
  
  visited = td_chordal->get_visited();
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Mark Triangle with boundary edge
//
// Special Notes : 
//
// Creator       : Jit Ken Tan
//
// Creation Date : 01/20/03
//-------------------------------------------------------------------------
CubitStatus ChordalAxis::flag_boundary(CubitFacet *curr_facet, CubitPoint *start_node, CubitPoint *end_node){
  ToolData *td = curr_facet->get_TD( TDChordal::is_chordal );
  TDChordal *td_chordal = dynamic_cast<TDChordal *> (td);
  int index;

  if (td_chordal == NULL) {
    td_chordal = new TDChordal();
    curr_facet -> add_TD( td_chordal);
  }
  
  int dummy;
  index =  curr_facet->edge_index(start_node, end_node, dummy);
  return  (td_chordal->flag_boundary_edge(index));

}

//-------------------------------------------------------------------------
// Purpose       : Classifying the triangles based on number of boundary edges
//
// Special Notes : 
//
// Creator       : Jit Ken Tan
//
// Creation Date : 01/20/03
//-------------------------------------------------------------------------
CubitStatus ChordalAxis::process_triangles(){
  int ii;
  CubitStatus status = CUBIT_SUCCESS;

  for( ii = triangulation.size(); ii > 0; ii--){
    CubitFacet *curr_facet = triangulation.get_and_step();
    ToolData *td = curr_facet-> get_TD( TDChordal::is_chordal );
    TDChordal *td_chordal = dynamic_cast<TDChordal *> (td);
   
    // it is a junction triangle if there is no tool attached
    if (td_chordal == NULL){
      td_chordal = new TDChordal();
      curr_facet-> add_TD(td_chordal);
    }          
    
    status = td_chordal -> determine_tritype();
    if(status == CUBIT_FAILURE)
      return status;
 
  }
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Determine triangle type based on number of boundary edges
//
// Special Notes : 2 boundary edges -> termination triangle
//                 1 boundary edges -> sleeve triangle
//                 0 boundary edges -> junction triangle
//
// Creator       : Jit Ken Tan
//
// Creation Date : 01/20/03
//-------------------------------------------------------------------------
CubitStatus ChordalAxis::determine_tritype(CubitFacet *curr_facet){
  CubitStatus status = CUBIT_SUCCESS;
  ToolData *td = curr_facet-> get_TD( TDChordal::is_chordal );
  TDChordal *td_chordal = dynamic_cast<TDChordal *> (td);
    
  if (td_chordal == NULL){
    return CUBIT_FAILURE;
  }          
    
  status = td_chordal->determine_tritype();

  return status;
}

//-------------------------------------------------------------------------
// Purpose       : Perform the chordal axis extraction
//
// Special Notes : 
//
// Creator       : Jit Ken Tan
//
// Creation Date : 01/20/03
//-------------------------------------------------------------------------
CubitStatus ChordalAxis::axis_extraction(DLIList <CubitFacetEdge *> &axis){
 
  triangulation.reset();
  int ii;
  TriType tri_class;
  CubitStatus status;

  for (ii = triangulation.size(); ii > 0; ii--){
    CubitFacet *current_facet = triangulation.get_and_step();
    TDChordal *td_chordal = get_tool_data(current_facet);
    tri_class = td_chordal->get_tritype();
    if(tri_class == UNDEFINED){
      PRINT_ERROR("Tri_class is undefined\n");
      return CUBIT_FAILURE;
    } else if(tri_class == SLEEVE) {

      status = sleeve_triangle(current_facet, td_chordal, axis);
      if(status == CUBIT_FAILURE){
	PRINT_ERROR("Error in Sleeve Triangles\n");
	return status;
      }

    } else if(tri_class == JUNCTION){
      junction_triangle(current_facet, td_chordal,axis);
    }
  }
  
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Perform chordal axis extraction given a sleeve triangle
//
// Special Notes : 
//
// Creator       : Jit Ken Tan
//
// Creation Date : 01/20/03
//-------------------------------------------------------------------------
CubitStatus ChordalAxis::sleeve_triangle(CubitFacet* curr_facet,TDChordal* td_chordal, DLIList <CubitFacetEdge *> &axis){
  
  int index1, index2;
  DLIList<int> edge_list;
  CubitPoint *point11, *point12, *point21, *point22;
  CubitVector vector11, vector12, vector21, vector22;
  CubitPoint *mid_point1, *mid_point2;
  CubitVector mid1, mid2;
  CubitFacetEdge *segment;

  //get non_boundary edges
  td_chordal->get_non_boundary_edges(edge_list);

  if(edge_list.size() != 2){
    PRINT_ERROR("Error in Sleeve Triangles: Boundary Edges Wrong. \n");
  }
  
  index1 = edge_list.get_and_step();
  index2 = edge_list.get_and_step();

  curr_facet->get_edge_pts(index1, point11, point12);
  curr_facet->get_edge_pts(index2, point21, point22);

  vector11 = point11->coordinates();
  vector12 = point12->coordinates();
  vector21 = point21->coordinates();
  vector22 = point22->coordinates();

  mid1 = (vector11 + vector12)/2;
  mid2 = (vector21 + vector22)/2;

  mid_point1 = (CubitPoint *) new CubitPointData(mid1);
  mid_point2 = (CubitPoint *) new CubitPointData(mid2);    
 
  segment = (CubitFacetEdge *) new CubitFacetEdgeData(mid_point1, mid_point2);
  
  axis.append(segment);
  
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Perform chordal axis extraction given a junction triangle
//
// Special Notes : 
//
// Creator       : Jit Ken Tan
//
// Creation Date : 01/20/03
//-------------------------------------------------------------------------
CubitStatus ChordalAxis::junction_triangle(CubitFacet *current_facet, TDChordal* td_chordal, DLIList <CubitFacetEdge *> &axis){
 
  CubitPoint *point11, *point12, *point21, *point22, *point31, *point32;
  CubitVector vector11, vector12, vector21, vector22, vector31, vector32;
  double mag_edge_one, mag_edge_two, mag_edge_three;
  CubitPoint *mid_point1, *mid_point2, *mid_point3;
  CubitVector mid1, mid2, mid3;
  CubitVector center;
  CubitPoint *center_point;
  CubitFacetEdge *segment1, *segment2, *segment3;
  CubitVector facet_vec1, facet_vec2, facet_vec3;
  bool obtuse = false;
  int largest_edge_index = -1;

  // obtain the points of the edges 
  current_facet->get_edge_pts(0, point11, point12);
  current_facet->get_edge_pts(1, point21, point22);
  current_facet->get_edge_pts(2, point31, point32);
  
  //get vector equivalents of the points
  vector11 = point11->coordinates();
  vector12 = point12->coordinates();
  vector21 = point21->coordinates();
  vector22 = point22->coordinates();
  vector31 = point31->coordinates();
  vector32 = point32->coordinates();

  
  //calculate the length of each edge
  mag_edge_one = (vector11 - vector12).length();
  mag_edge_two = (vector21 - vector22).length();
  mag_edge_three = (vector31 - vector32).length();
    
  //Checking whether triangle is acute or obtuse
  if(pow(mag_edge_one,2) >= pow(mag_edge_two,2) + pow(mag_edge_three,2)){
    obtuse = true;
    largest_edge_index = 1;
  } else if(pow(mag_edge_two,2) >= pow(mag_edge_one,2) + pow(mag_edge_three,2)){
    obtuse = true;
    largest_edge_index = 2;
  }else if (pow(mag_edge_three,2) >= pow(mag_edge_one,2) +  pow(mag_edge_two,2)){
    obtuse = true;
    largest_edge_index = 3;
  } else {
    obtuse = false;
  }

  mid1 =(vector11+vector12)/2;
  mid2 =(vector21+vector22)/2;
  mid3 =(vector31+vector32)/2;
  
  mid_point1 = (CubitPoint *) new CubitPointData(mid1);
  mid_point2 = (CubitPoint *) new CubitPointData(mid2);
  mid_point3 = (CubitPoint *) new CubitPointData(mid3);

  // Extracting the chordal axis.
  if(obtuse == false){
    //triangle is acute.
    //Connect the circumcenter and the midpoints of the edges
    facet_vec1 = (current_facet->point(0)) -> coordinates();
    facet_vec2 = (current_facet->point(1)) -> coordinates();
    facet_vec3 = (current_facet->point(2)) -> coordinates();
    
    center = get_center(facet_vec1, facet_vec2, facet_vec3);

    center_point = (CubitPoint *) new CubitPointData(center);

    segment1 = (CubitFacetEdge *) new CubitFacetEdgeData(mid_point1, center_point);
    segment2 = (CubitFacetEdge *) new CubitFacetEdgeData(mid_point2, center_point);
    segment3 = (CubitFacetEdge *) new CubitFacetEdgeData(mid_point3, center_point);

    axis.append(segment1);
    axis.append(segment2);
    axis.append(segment3);

  
  } else {
    //triangle is obtuse

    switch(largest_edge_index){
    case 1:
      segment1 = (CubitFacetEdge *) new CubitFacetEdgeData(mid_point1, mid_point2);
      segment2 = (CubitFacetEdge *) new CubitFacetEdgeData(mid_point1, mid_point3);
      break;
    case 2:
      segment1 = (CubitFacetEdge *) new CubitFacetEdgeData(mid_point2, mid_point1);
      segment2 = (CubitFacetEdge *) new CubitFacetEdgeData(mid_point2, mid_point3);
      break;
    case 3:
      segment1 = (CubitFacetEdge *) new CubitFacetEdgeData(mid_point3, mid_point1);
      segment2 = (CubitFacetEdge *) new CubitFacetEdgeData(mid_point3, mid_point2);
      break;
    default:
      //Not a possible scenario
      assert(0); 
      break;
    }
    
    axis.append(segment1);
    axis.append(segment2);
  }

  return CUBIT_SUCCESS;

}

//-------------------------------------------------------------------------
// Purpose       : get circumcenter given 3 points of a circle
//
// Special Notes : 
//
// Creator       : Jit Ken Tan
//
// Creation Date : 01/20/03
//-------------------------------------------------------------------------
CubitVector ChordalAxis::get_center(CubitVector v1, CubitVector v2, CubitVector v3){
  double x1,x2,x3,y1,y2,y3;
  double ma, mb, cx, cy;

  x1 = v1.x();
  y1 = v1.y();
  x2 = v2.x();
  y2 = v2.y();
  x3 = v3.x();
  y3 = v3.y();
  
  //for debugging
  bool center_found = false;

  if(x2 != x1 && x3 != x2 ){
    ma = (y2-y1)/(x2-x1);
    mb = (y3-y2)/(x3-x2);
    cx = (ma*mb*(y1-y3) + mb*(x1+x2) - ma*(x2+x3))/(2*(mb-ma));
    
    if(ma != 0)
      cy = (-1/ma)*(cx - (x1+x2)/2)+ (y1+y2)/2;
    else
      cy = (-1/mb)*(cx - (x3+x2)/2)+ (y3+y2)/2;

    center_found = true;

  } else if (x3 != x1 && x2 != x3) {
    ma = (y3-y1)/(x3-x1);
    mb = (y2-y3)/(x2-x3);
    cx = (ma*mb*(y1-y2) + mb*(x1+x3) - ma*(x2+x3))/(2*(mb-ma));
    if(ma != 0)
      cy = (-1/ma)*(cx - (x1+x3)/2)+ (y1+y3)/2;
    else 
      cy = (-1/mb)*(cx - (x3+x2)/2)+ (y3+y2)/2;

    center_found = true;

  } else if ( x1 != x2 && x1 != x3 ){
    ma = (y1-y2)/(x1-x2);
    mb = (y3-y1)/(x3-x1);    
    cx = (ma*mb*(y2-y3) + mb*(x1+x2) - ma*(x1+x3))/(2*(mb-ma));
    if(ma != 0)
      cy = (-1/ma)*(cx - (x1+x2)/2)+ (y1+y2)/2;
    else
      cy = (-1/mb)*(cx - (x1+x3)/2)+ (y1+y3)/2;

    center_found = true;
  } else {
    //scenario not possible
    assert(0);
  }

  assert(center_found);

  return CubitVector(cx,cy,0);
}



