//- Class: FacetorUtil
//- Description: Utility functions used by FactorTool and available to other classes.
//               FacetorUtil uses template functions that accepts both CubitNode, CubitTri, CubitEdge classes 
//-              and CubitPoint, CubitFacet, CubitFacetEdge classes.  CubitPoint, CubitFacet, and 
//-              CubitFacetEdge classes have been modified to account for the difference between the mdb
//=              entities and geometric entities.
//-------------------------------------------------------------------------
// Filename      : FacetorTool.hpp
//
// Purpose       : 2D Delaunay Mesher
//
// Creator       : Christopher Hynes
//
// Creation Date : 5/31/2002
//
// Owner         : Steve Owen
//-------------------------------------------------------------------------

#ifndef FACETOR_UTIL_HPP
#define FACETOR_UTIL_HPP


#include "CubitDefines.h"
class ParamTool;

template<class SURF, class TRI, class EDGE, class NODE, class TRICHILD, class NODECHILD, class SIZEFUNC>
class FacetorUtil
{
public:
  FacetorUtil(){}
  ~FacetorUtil(){}

    //- Insert one node into a Delaunay mesh
  static CubitStatus insert_node(NODE *node_ptr,
                                 DLIList<TRI*>& facet_list,
                                 SURF* surface_ptr,
                                 int& curr_visit_flag,
                                 TRI*& last_tri);

    // Insert a node at the circumcenter of a tri
  static CubitStatus insert_at_circumcenter(TRI *tri_ptr,
                                            DLIList<TRI*>& facet_list,
                                            TRI*& start_tri,
                                            int& curr_visit_flag,
                                            SURF* surface_ptr,
                                            DLIList<TRI*>* sorted_lists,
                                            const int num_lists,
                                            const double interval,
                                            const double quality_angle,
                                            SIZEFUNC* sizing_function,
                                            ParamTool* p_tool);

    // Compute the angles at the triangle vertices and classify
    // by its worst triangle
  static CubitStatus classify_tri_by_angle(TRI* tri_ptr,
                                           DLIList<TRI*>* sorted_lists,
                                           const int num_lists,
                                           const double interval,
                                           const double quality_angle);

private:

    //Get the tri_visited flag.
  static CubitBoolean tri_visited(TRI *tri_ptr,
                                  int curr_visit_flag);

    //Set the tri_visited flag.
  static void tri_visited(TRI *tri_ptr,
                          CubitBoolean visited,
                          int curr_visit_flag);
  
    //Search every TRI until the one that the_point is in
    //is found.
  static CubitStatus exhaustive_locate_point(CubitVector& the_point,
                                             DLIList<TRI*>& facet_list,
                                             TRI*& tri_ptr);

  static CubitStatus point_in_circumcircle(CubitVector& the_point,
                                           TRI* tri_ptr,
                                           DLIList<TRI*>& tri_list,
                                           SURF* surface_ptr,
                                           int curr_visit_flag);

    //Find the TRI in facet_list that the_point is in.
    //If starting_tri is not NULL the search will begin
    //with it.
  static CubitStatus locate_point(CubitVector& the_point,
                                  DLIList<TRI*>& facet_list,
                                  TRI* starting_tri,
                                  SURF* owning_surface,
                                  TRI*& tri_ptr);

    //check the void created in Bowyer-Watson to ensure the boundary
    //seems valid.  
  static CubitStatus valid_void( NODE * /*point_ptr*/,
                                 DLIList<TRI *> &tri_list,
                                 SURF* surface_ptr,
                                 int curr_visit_flag);

  static CubitStatus natural_neighbor_tris(CubitVector& the_point,
                                           DLIList<TRI*>& facet_list,
                                           TRI*& start_tri,
                                           SURF* surface_ptr,
                                           int& curr_visit_flag,
                                           DLIList<TRI*>& tri_list);

   //- Bowyer-Watson insertion into an existing Delaunay Mesh
  static CubitStatus bowyer_watson_insert(NODE* point_ptr,
                                          DLIList<TRI*>& tri_list,
                                          DLIList<TRI*>& facet_list,
                                          int& curr_visit_flag,
                                          SURF* surf_ptr,
                                          TRI*& last_tri);

    //Return true if all three points in tri_ptr are colinear.
    //Return false otherwise.
  static CubitBoolean are_nodes_colinear(TRI* tri_ptr);

  static CubitVector& circumcenter(TRI* tri_ptr);

  static double get_size(CubitVector &cc, TRI *tri_ptr);

    // get and set the tri sort array index that this tri is in
  static void tri_sort_list( TRI *tri_ptr, int sort_list_index );
  static int tri_sort_list( TRI *tri_ptr );

  static double radius(TRI* tri_ptr);  
};

#include "FacetorUtil.cpp"
#endif
