//-------------------------------------------------------------------------
// Filename      : ChordalAxis.hpp
//
// Purpose       :
//
// Special Notes : 
//
// Creator       : Jit Ken Tan
//
// Creation Date : 09/02/02
//-------------------------------------------------------------------------

#ifndef CHORDALAXIS_HPP
#define CHORDALAXIS_HPP

#include "CubitDefines.h"
#include "DLIList.hpp"
#include "CubitFacet.hpp"
#include "CubitFacetEdge.hpp"
#include "CubitPoint.hpp"
#include "GfxDebug.hpp"

class TDChordal;
class CubitVector;
class ChordalAxis
{
 
  public:

  ChordalAxis();

  ~ChordalAxis();
				 
  //- Perform the entire chordal axis extraction
  CubitStatus chordal_axis_transform(DLIList <CubitPoint*> boundary_loop_list,DLIList <CubitFacet *> &facet_list, DLIList <CubitFacetEdge *> &axis);
   
  private:
  
  //- create a delaunay constrained triangulation
  CubitStatus create_triangulation(DLIList <CubitPoint*> boundary_loop_list);

  //- extract the chordal axis for a junction triangle
  CubitStatus junction_triangle(CubitFacet *current_facet, TDChordal * td_chordal, DLIList <CubitFacetEdge *> &axis);
   
  //-Get circumcenter given 3 points of the circle
  CubitVector get_center(CubitVector v1, CubitVector v2, CubitVector v3);

  //- Flag the boundary edges
  CubitStatus flagging_boundary();

  //- Flag a boundary edge
  CubitStatus flag_boundary( CubitFacet *curr_facet, CubitPoint *start_node, CubitPoint *end_node );

  //- Classified each triangle based on its boundary edges
  CubitStatus process_triangles();

  //- perform extraction
  CubitStatus axis_extraction( DLIList <CubitFacetEdge *> &axis );
  
  //- Obtain TDChordal Tool data
  TDChordal *get_tool_data(CubitFacet* curr_facet);
  
  //- Perform chordal axis extraction given a sleeve triangle
  CubitStatus sleeve_triangle(CubitFacet* curr_facet,TDChordal* td_chordal, DLIList <CubitFacetEdge *> &axis);
  
  //- Perform pruning of the chordal axis
  CubitStatus pruning();
  
  //- Searching for triangles to prune
  CubitStatus search_triangles_to_prune(CubitFacet *adj_triangle, CubitFacetEdge *curr_edge, CubitPoint *start_node, CubitVector unit_edge_vector, double pruning_distance, DLIList <CubitFacet *> &pruned_triangles, CubitStatus &prune );
  
  //- Pruning commences with the junction triangle
  CubitStatus junction_triangle_pruning(CubitFacet *curr_facet, int edge_index, bool &pruning_performed);
  
  //- Mark Triangle as discarded
  CubitStatus mark_tri_discarded(CubitFacet *curr_facet);

  //- Mark Triangle as visited
  CubitStatus mark_visited(CubitFacet *curr_facet);
  
  //- Mark Triangle as unvisited
  CubitStatus unmark_visited(CubitFacet *curr_facet);
  
  //- Check if the triangle is visited
  CubitStatus get_visited(CubitFacet *curr_facet, bool &visited);

  //- Determine the triangle type based on the number of boundary edges
  CubitStatus determine_tritype(CubitFacet *curr_facet);

  //- Consolidate the junction triangles in our mesh
  CubitStatus junction_tris_collection();

  //- for debugging purposes
  void debug_draw_facet(CubitFacet *curr_facet, int color){
    DLIList<CubitFacetEdge *> edge_list;
    curr_facet->edges(edge_list);
    int ii;
    
    for(ii =0; ii<edge_list.size(); ii++){
      GfxDebug::draw_facet_edge(edge_list.get_and_step(), color);
    }

  };

  //- for debugging purposes
  void print_point(CubitPoint *curr_point){
    PRINT_INFO("Point  = (%f, %f) \n", curr_point->x(), curr_point->y());
  }

  int numBoundaryEdges;
  CubitPoint **startNodes;
  CubitPoint **endNodes;
  DLIList <CubitFacet*> triangulation;
  DLIList <CubitFacet*> junctionTris;
  double pruningAlpha;
};

#endif

