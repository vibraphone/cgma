//- Class:       FacetDataUtil
//- Description: static library of general functions for querying and/or modifying
//-              facet entities
//- Owner:       Steve Owen
//- Checked by:
//- Version:

#ifndef __FACETDATAUTIL__
#define __FACETDATAUTIL__

#include "DLIList.hpp"
#include "RTree.hpp"

class CubitPoint;
class CubitFacet;
class CubitFacetEdge;
class CubitQuadFacet;
class FacetEntity;
class CubitVector;
class CubitBox;

class FacetDataUtil
{
public:
  // get edges from a facet list that are common to "count" facets
  // For example, passing in:
  //   count = 1 returns the boundary edges of a set of facets
  //   count = 2 returns the interior, manifold edges of a set of facets
  //   count = n>2 returns non-manifold edges shared by n facets
  //
  static void edges_by_count(DLIList<CubitFacet*> &facets,
 										   unsigned int count,
 										   DLIList<CubitFacetEdge*> &edges);

  // get an ordered list of points and edges around the boundary of a set of facets
  // The output "chain" consists of point - edge - point - edge - ... around the
  // boundary of the input facets.  Having the ordered points as well as the edges
  // helps provides sense information for the bounding edges.
  static void ordered_point_edge_bdry(DLIList<CubitFacet*> &facets,
                                      DLIList<FacetEntity*> &point_edge_chain);

  // given an ordered point - edge boundary and a start and end point on the boundary,
  // return the point - edge chain between the two points, inclusive.
  static CubitStatus partial_chain(DLIList<FacetEntity*> &point_edge_chain,
                                   FacetEntity* point1,
                                   FacetEntity* point2,
                                   DLIList<FacetEntity*> &chain_between);

  // given a set of facets, get the unique list of points used to define the facets.
  static void get_facet_points(DLIList<CubitFacet*> &cubit_facets,
                               DLIList<CubitPoint*> &facet_points);

  // given a set of facets, return an unordered list of points at its boundary
  static void get_boundary_points(DLIList<CubitFacet*> &facet_list,
                                  DLIList<CubitPoint*> &point_list);

  // return an unordered list of edges at the boundary of a set of facets
  static void get_boundary_edges(DLIList<CubitFacet*> &facet_list,
                                 DLIList<CubitFacetEdge*> &edge_list);

  // get the unique set of points from a facet list
  static void get_points( DLIList<CubitFacet*> &facet_list,
                          DLIList<CubitPoint*> &point_list);

  // make a complete copy of facets and points
  static void copy_facets(DLIList<CubitFacet*> &old_facet_list,
                          DLIList<CubitFacet*> &new_facet_list,
                          DLIList<CubitPoint*> &new_point_list,
                          DLIList<CubitFacetEdge*> &new_edge_list );

  // populates the edge list from the list of facets - creates the
  // edges if necessary
  static void get_edges( DLIList<CubitFacet *> &facet_list,
                         DLIList<CubitFacetEdge *> &edge_list );

  // compute the quality of a triangle defined by three vertices
  static double quality(CubitVector &c1, CubitVector &c2, CubitVector &c3,
                              CubitVector &surf_normal);

  // collapse edges that are shorter than 10% of the average edge length
  static CubitStatus collapse_short_edges( DLIList<CubitFacet*> &facet_list,
                                           CubitBoolean non_manifold_only );

  static void check(DLIList<CubitFacetEdge *> &edge_list,
                    DLIList<CubitFacet *> &facet_list);

  static void draw_edge( CubitFacetEdge *edge );

  static int edge_compare(CubitFacetEdge *&ea, CubitFacetEdge *&eb);

  // for triangulated, manifold, non-self-intersecting polyhedra
  static CubitStatus is_point_in_polyhedron(DLIList<CubitFacet *> &tfacet_list,
                                     CubitVector &point_coords,
                                     CubitPointContainment &is_point_in);

  //  can return CUBIT_FAILURE if triangles are very thin
  //  used by is_point_in_polyhedron()
  static CubitStatus pt_in_tri_2d(double xpt, double ypt,
                       double x0, double y0,
		       double x1, double y1,
		       double x2, double y2,
		       CubitPointContainment &is_point_in);


  //  Returns a normalized ray always pointing in the (+++) direction.
  //  used by is_point_in_polyhedron()
  static void random_positive_ray(CubitVector& ray);

  //  used by is_point_in_polyhedron()
  static bool ray_intersects_boundingbox(CubitVector& point, CubitVector& ray, const CubitBox& bbox);

  // write a facet list to a facet file
  static CubitStatus write_facets( const char *file_name,
                                   DLIList<CubitFacet *> &facet_list);

  // group facets into continuous lists
  static CubitStatus split_into_shells( DLIList<CubitFacet *> &tfacet_list,
                                 DLIList<CubitQuadFacet *> &qfacet_list,
                                 DLIList<DLIList<CubitFacet *> *> &shell_list,
                                 CubitBoolean &is_water_tight);

  // attempt to stitch up the model to make water-tight
  static CubitStatus stitch_facets( DLIList<DLIList<CubitFacet *> *> &shell_list,
                             double tol, CubitBoolean &is_water_tight,
                             CubitBoolean write_results = CUBIT_TRUE);

  // delete a list of facets and their child facet entities
  static void delete_facets(DLIList<DLIList<CubitFacet*>*> &shell_list);
  static void delete_facets(DLIList<CubitFacet*> &facet_list);
  static void delete_facet(CubitFacet *facet_ptr);
  static void destruct_facet_no_delete(CubitFacet *facet_ptr);
  

  // determine intersection of a segment with a facet
  // returns CUBIT_PNT_UNKNOWN: if segment is in plane of facet
  //         CUBIT_PNT_OUTSIDE: if segment does not intersect
  //         CUBIT_PNT_INSIDE: if segment intersects inside facet
  //         CUBIT_PNT_BOUNDARY: if segment intersects a vertex or edge
  static CubitPointContainment intersect_facet(CubitVector &start, CubitVector &end,
                                               CubitFacet *facet_ptr,
                                               CubitVector &qq,
                                               CubitVector &ac,
                                               CubitBoolean bound = CUBIT_TRUE);

   //  get axis-aligned bounding box of list of points
  static CubitStatus get_bbox_of_points(DLIList<CubitPoint*>& point_list, CubitBox& bbox);

  //  get square of distance from point to closest point on a line segment;
  //  returns closest point.
  static CubitVector squared_distance_to_segment(CubitVector &p, 
                                                  CubitVector &p0, 
                                                  CubitVector &p1,
                                                  double &distance2);

  //  Get the intersection of the line defined by point1 and point2 with
  //  bbox.  Returns 0,1 or 2 for the number of intersections.  A line
  //  in one of the planes of the box will return 0.
  static int get_bbox_intersections(CubitVector& point1,
                                    CubitVector& point2,
                                    const CubitBox& bbox,
                                    CubitVector& intersection_1,
                                    CubitVector& intersection_2);
  
  // mark facets and their children
  static void mark_facets( DLIList<FacetEntity *> &facet_list, int mark_value );
  
  // used with stitch_facets to merge points that are coincident
  static CubitStatus merge_coincident_vertices( DLIList<DLIList<CubitFacet *> *> &shell_list,
                                         double tol,
                                         int &mpmerge,
                                         int &nemerge,
                                         DLIList<CubitPoint *> &unmerged_points);  

  static CubitStatus merge_coincident_vertices( DLIList<CubitPoint*> &points );
  
  // used with stitch_facets to merge two points together that have already been
  // determined to be coincident.  Also handles their adjacent edges
  static CubitStatus merge_points( CubitPoint *pt0, CubitPoint *pt1, int &nemerge, 
                                    RTree<CubitFacetEdge*> *r_tree=NULL);

  static CubitStatus find_facet(  DLIList<CubitFacet *> temp_split_facets, 
                                  CubitPoint *pnt0, 
                                  CubitPoint *pnt1, 
                                  DLIList<CubitFacet *> &select_facets );


private:

 // used with split_into_shells to group facets
  static CubitStatus get_adj_facets_on_shell( CubitFacet *start_facet_ptr,
                                       DLIList<CubitFacet *> *shell_ptr,
                                       CubitBoolean &is_water_tight,
                                       int mydebug );

  // used with stitch_facets to merge points that are colinear with edges
  static CubitStatus merge_colinear_vertices( DLIList<DLIList<CubitFacet *> *> &shell_list,
                                         double tol,
                                         DLIList<CubitPoint *> &merge_points,
                                         int &mpmerge,
                                         int &nemerge,
                                         int &nnomerge);


};

#endif

