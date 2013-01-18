///-----------------------------------------------------------------------------------
// Class: PointLoopFacetor
// Description:
// Creates a set of CubitFacets which is a Delauney triagulation of the boundary
// points.  No interior points are added.  The points form a closed loop, and may
// consist of several loops.  The first loop must be the exterior loop, followed by
// any internal loops (or holes).
// Creator: David R. White
// Owner: David R. White
// Creation Date: 3/1/2003
//------------------------------------------------------------------------------------
#include "PointLoopFacetor.hpp"
#include "DLIList.hpp"
#include "CubitPoint.hpp"
#include "CubitFacet.hpp"
#include "CubitFacetData.hpp"
#include "CubitFacetEdge.hpp"
#include "CubitPointData.hpp"
#include "LoopParamTool.hpp"
#include "TDVector.hpp"
#include "FacetorTool.hpp"
#include <stdio.h>
CubitStatus PointLoopFacetor::generate_facets( PointLoopList &boundary_loops,
                                               DLIList<CubitFacet*> &resulting_triangles)
{
    //First transform the boundary points to the XY plane.
  LoopParamTool *l_tool = new LoopParamTool;
  if ( l_tool->new_space_LoopParam(boundary_loops) != CUBIT_SUCCESS )
  {
    delete l_tool;
    
    PRINT_ERROR("Faild to generate loop u-v space for generating facets of boundary points.\n");
    return CUBIT_FAILURE;
  }
    //Now transform the poins.
  int ii, jj;
  if ( l_tool->transform_loopspoints_to_uv(boundary_loops) != CUBIT_SUCCESS )
  {
    write_xy(boundary_loops);
      //Now transform the points back to their original location and remomve the TDVector.
    CubitPoint *curr_point;
    ToolData *td;
    TDVector *td_vec;
    CubitVector orig_vec;
    PointList bounding_points, *curr_list;
    for ( ii = 0; ii < boundary_loops.size(); ii++)
    {
      curr_list = boundary_loops.get_and_step();
      bounding_points += (*curr_list);
    }
    for ( ii = 0; ii < bounding_points.size(); ii++ )
    {
      curr_point = bounding_points.get_and_step();
      td = curr_point->get_TD(&TDVector::is_td_vector);
      if ( td == NULL )
      {
        PRINT_ERROR("No TDVector on CubitPoint, can't transform back to x-y-z.\n");
        assert(td !=NULL);
      }
      td_vec = CAST_TO(td, TDVector);
      orig_vec = td_vec->get_vector();
      curr_point->set(orig_vec);
      curr_point->delete_TD(&TDVector::is_td_vector);
    }
    write_xyz(boundary_loops);
    
    delete l_tool;
    PRINT_ERROR("Faild to transform loop to u-v space for generating facets of boundary points.\n");
    return CUBIT_FAILURE;
  }
  delete l_tool;
    //now generate a triangulation for these nodes.
  int num_boundary_edges = 0;
  PointList bounding_points, *curr_list;
  for ( ii = 0; ii < boundary_loops.size(); ii++)
  {
    curr_list = boundary_loops.get_and_step();
    num_boundary_edges += curr_list->size();
    bounding_points += (*curr_list);
  }
  CubitPoint **start_nodes = new CubitPoint* [num_boundary_edges];
  CubitPoint **end_nodes = new CubitPoint* [num_boundary_edges];
  int array_index = 0;
  for ( ii = 0; ii < boundary_loops.size(); ii++ )
  {
    curr_list = boundary_loops.get_and_step();
    for ( jj = 0; jj < curr_list->size(); jj++ )
    {
      assert(array_index < num_boundary_edges);
      start_nodes[array_index] = curr_list->get_and_step();
      end_nodes[array_index] = curr_list->get();
      array_index++;
    }
  }
  if ( array_index != bounding_points.size() )
  {
    PRINT_ERROR("Problem building edge arrays for faceting.\n");
    assert(array_index == bounding_points.size());
    return CUBIT_FAILURE;
  }
  int dummy_variable = -1;
  int surf_id = 1;
  FacetorTool<int,CubitFacet, CubitFacetEdge,
              CubitPoint, CubitFacetData, CubitPointData, int>
    facetor = FacetorTool<int,CubitFacet, CubitFacetEdge, CubitPoint,
              CubitFacetData, CubitPointData, int> (&surf_id, bounding_points,
                                                    start_nodes, end_nodes, num_boundary_edges,
                                                    &dummy_variable, l_tool);
  if (facetor.mesh_surfwoIP(resulting_triangles) != CUBIT_SUCCESS )
  {
    PRINT_ERROR("The facetor tool failed.\n");
    return CUBIT_FAILURE;
  }
    //Now transform the points back to their original location and remomve the TDVector.
  CubitPoint *curr_point;
  ToolData *td;
  TDVector *td_vec;
  CubitVector orig_vec;
  for ( ii = 0; ii < bounding_points.size(); ii++ )
  {
    curr_point = bounding_points.get_and_step();
    td = curr_point->get_TD(&TDVector::is_td_vector);
    if ( td == NULL )
    {
      PRINT_ERROR("No TDVector on CubitPoint, can't transform back to x-y-z.\n");
      assert(td !=NULL);
    }
    td_vec = CAST_TO(td, TDVector);
    orig_vec = td_vec->get_vector();
    curr_point->set(orig_vec);
    curr_point->delete_TD(&TDVector::is_td_vector);
  }
    //Now update the facet planes...
  CubitFacet *facet;
  for ( ii = 0; ii < resulting_triangles.size(); ii++ )
  {
    facet = resulting_triangles.get_and_step();
    facet->reset_bounding_box();
    facet->update_plane();
  }
  return CUBIT_SUCCESS;
}
void PointLoopFacetor::write_xyz(PointLoopList &boundary_loops)
{
  
  FILE *fp = fopen("Loop.txt", "w");
  if (!fp)
  {
    PRINT_ERROR("Couldn't open temp file for writing.\n");
    return;
  }
  PointList *node_loop;
  CubitPoint *point_ptr;
  CubitVector coords;
  fprintf(fp, "Real locations\n");
  int ii, jj;
  for ( ii = boundary_loops.size(); ii > 0; ii-- )
  {
    node_loop = boundary_loops.get_and_step();
    fprintf(fp, "Loop %d, contains %d nodes\n", ii, node_loop->size());
    node_loop->reset();
    for(jj=0; jj<node_loop->size(); jj++)
    {
      point_ptr = node_loop->get_and_step();
      coords = point_ptr->coordinates();
      fprintf(fp,"%d\t%f\t%f\t%f\n", jj, coords.x(), coords.y(), coords.z());
    }
  }
  fclose(fp);
}
void PointLoopFacetor::write_xy(PointLoopList &boundary_loops)
{
  
  FILE *fp = fopen("Loop_transformed.txt", "w");
  if (!fp)
  {
    PRINT_ERROR("Couldn't open temp file for writing.\n");
    return;
  }
  PointList *node_loop;
  CubitPoint *point_ptr;
  CubitVector coords;
  fprintf(fp, "Transformed locations\n");
  int ii, jj;
  for ( ii = boundary_loops.size(); ii > 0; ii-- )
  {
    node_loop = boundary_loops.get_and_step();
    fprintf(fp, "Loop %d, contains %d nodes\n", ii, node_loop->size());
    node_loop->reset();
    for(jj=0; jj<node_loop->size(); jj++)
    {
      point_ptr = node_loop->get_and_step();
      coords = point_ptr->coordinates();
      fprintf(fp,"%d\t%f\t%f\n", jj, coords.x(), coords.y());
      
    }
  }
  fclose(fp);
}
