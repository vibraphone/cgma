//-------------------------------------------------------------------------
// Filename      : PartSurfFacetTool.hpp
//
// Purpose       : Code for initializing PartitionSurface facetting and
//                 other misc. facet-related code moved from PartitionSurf.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/11/04
//-------------------------------------------------------------------------


#ifndef PART_SURF_FACET_TOOL_HPP
#define PART_SURF_FACET_TOOL_HPP

#include <set>
#include <DLIList.hpp>

class CubitVector;
class Curve;

class PartitionSurface;
class PartitionCurve;
class PartitionPoint;

class CubitFacet;
class CubitFacetEdge;
class CubitPoint;

class CubitFacetData;
class CubitFacetEdgeData;
class CubitPointData;

class PartSurfFacetTool
{
  
public:
  
  /** Interface for initializing PartitionSurface facets */

  PartSurfFacetTool(PartitionSurface* surface);

  CubitStatus init_facet_data( DLIList<CubitFacetData*>& facets );
  
  /** Debugging functions */
  static void validate_facets(PartitionSurface* surf);

  /** Remove small edges */
  
  /** Misc. static utility functions */

  static CubitStatus split_edge( CubitFacetEdge*    edge_to_split,
                                 const CubitVector& split_location,
                                 CubitFacet*        edge_facet,
                                 CubitPoint*&       new_split_point,
                                 CubitFacetEdge*&   new_edge,
                                 CubitFacet*&       new_facet );
  //- Split a facet edge and any adjacent facets.
  //I edge_to_split
  //I- The edge to split
  //I split_location
  //I- The location at which to split the edge
  //I edge_facet 
  //I- Optional argument - facet adjacent to edge
  //O new_split_point
  //O- The new point splitting the edge
  //O new_Edge
  //O- The new edge (from new_split_point to the original end point)
  //O new_facet
  //O- If edge_facet is not NULL, the new facet resulting from
  //O- splitting edge_facet.
               
                                 
  static CubitStatus collapse_edge( CubitPoint* keep_point,
                                    CubitPoint* dead_point,
                                    DLIList<CubitFacetData*>* unowned = 0 );
  //- Combine the two passed points, removing the edge connecting
  //- them and collapsing any facets adjacent to that edge.

  static CubitFacet* closest_facet( const CubitVector& input_position,
                                    const DLIList<CubitFacetData*>& facets,
                                    CubitVector& result_position );
  //- Find the facet and point on that facet closest to the passed position.
  //R CubitFacet*
  //R- The closest facet
  //I input_position
  //I- The location to test
  //O facets
  //O- The set of facets to test
  //O result_position
  //O- The closest position on the returned facet.
  

  static CubitStatus seam_curve( DLIList<CubitFacetEdgeData*>& new_edge_list,
                                 PartitionCurve *curve,
                                 DLIList <CubitFacetData*> &adjacent_facets,
                                 DLIList <CubitFacetEdgeData*>* dead_edge_ptrs = 0 );
  //- Given two lists of facet edges, the passed one and the one associated
  //- with the passed curve, split and merge edges to seam the two lists
  //- together.
  //I new_edge_list
  //I- The list of edges to seam
  //I curve
  //I- The curve who's edges are to be seamed with the input list.
  //I adjacent_facets
  //I- The list of un-owned facets which may need to be modified during
  //I- this operation.  For any owned facets, the owning PartitionSurface
  //I- will be updated appropriately.  For any un-owned facets, this list
  //I- will be updated. 

  static CubitPoint* split_edge_closest( CubitFacetEdgeData* edge,
                                         const CubitVector& pos,
                                         double tolerance,
                                         CubitFacetEdgeData*& new_edge,
                                         DLIList<CubitFacetData*>& new_facets );
  //- If the input position is sufficiently close to the end of the passed
  //- edge, return the corresponding end point of the edge.  Otherwise
  //- split the edge and any adjacent facets.
  //-
  //- Effects:  For each adjacent facet that is split, if the facet has
  //-           an owning PartitionSurface, the owning surface will be
  //-           updated for the facet split.  OTHERWISE the new facet
  //-           will be appended to new_facets.
  //-
  //- Notes:    The edge is always split such that the original edge
  //-           is between the original start point and the new point
  //-           and the new edge contains the new point and the 
  //-           original end point.
  //-
  //R CubitPoint*
  //R- The new point used to split the edge, or the closest edge end
  //R- point if the edge was not split.
  //I edge
  //I- The edge to split
  //I pos
  //I- The location at which to split the edge
  //O new_edge
  //O- The new edge, or NULL if the edge was not split.
  //O new_facets
  //O- New facets created by splitting un-owned facets adjacent to 
  //O- the input edge.
  
  static CubitStatus fix_move_point( CubitPoint* point,
                                     const CubitVector& new_pos,
                                     const DLIList<CubitFacetData*>& facets,
                                     DLIList<CubitFacetData*>& old_facets,
                                     DLIList<CubitFacetData*>& new_facets,
                                     PartitionSurface* owning_surf = 0 );
  //- Given a patch of facets, a point on the perimeter of that patch
  //- of facets and a new position for that point, check if any re-facetting
  //- is required keep from inverting facets, and if so refacet.
  //I point
  //I- The point to move
  //I new_pos
  //I- The desired new position for the point
  //I facets
  //I- The patch of facets
  //O old_facets
  //O- Facets removed during refacetting.
  //O new_facets
  //O- New facets created during refacetting.
  //I owning_surf
  //I- Optional argument.  If specified, function assumes that all
  //I- facets in passed list are owned by this surface.  Uses facet
  //I- owner rather than list search to see if facet is in list.
  

  
  static CubitStatus get_facet_points_and_edges( 
                         const DLIList<CubitFacetData*>& facets,
                         DLIList<CubitPoint*>& boundary_points,
                         DLIList<CubitPoint*>& interior_points,
                         DLIList<CubitFacetEdge*>& boundary_edges,
                         DLIList<CubitFacetEdge*>& interior_edges );
  //- Given a patch of facets, get all points and edges classified
  //- as either boundary or interior.


private:
 
  static void edge_facets( CubitFacetEdge* edge,
                           const DLIList<CubitFacet*>& input_facets,
                           DLIList<CubitFacet*>& output_list );
  
  static void edge_facets( PartitionSurface* surface,
                           CubitFacetEdge* edge, 
                           DLIList<CubitFacet*>& facets );

  /** Helper functions for init_facet_data */
 
  static void closest_pt_on_facet( CubitFacet* facet,
                                   const CubitVector& input_position,
                                   CubitVector& result_position );

  
  CubitStatus seam_curves( DLIList<PartitionCurve*>& curve_list,
                           DLIList<CubitFacetEdge*>& edge_list,
                           DLIList<CubitFacetData*>& adjacent_facets );
  //- Do seam_curve for the set of all partitions of a real curve.
  
  static CubitStatus associate_points( DLIList<CubitPoint*>& facet_points,
                                       DLIList<PartitionPoint*>& geom_points );
  //- For each pointin geom_points, find closest point in facet_points
  //- and attach that facet point to the geometric point.  If the geometric
  //- point already has an attached facet point, the facet points are merged
  //- such that the point from facet_points is destroyed during the merge.
  
  static Curve* get_real_curve( DLIList<PartitionCurve*>& curve_list,
                                PartitionPoint*& start_point,
                                PartitionPoint*& end_point );
  //- Verify that the passed set of curves are exactly the set of
  //- partitions of a single real curve.  Return the real curve and
  //- the PartitionPoints corresponding to the start and end of the
  //- real curve.
  
  CubitStatus get_boundary_chain( CubitPoint* start_point,
                                         CubitFacetEdge* start_edge,
                                         CubitPoint* end_point,
                                         DLIList<CubitFacetEdge*>& result );
  //- Helper function for init_facet_data().
  //- NOTE: expects marks on facet edges to be set by init_facet_data.().
  //-
  //- Given an edge and a point on that edge, find the chain of
  //- boundary edges beginning with the edge and ending at the passesd
  //- end point
  //- Assumes:  all boundary edges marked with a '1'.
  //-           chain of edges may not cross a vertex with more than
  //-            two adjacent boundaryt edges.
  //-           chain of edges may not corss a vertex owned by a 
  //-            partition point.
  
  CubitStatus seam_nonmanifold_curves( DLIList<PartitionCurve*>& partitions,
                                       DLIList<CubitFacetData*>& facet_list );
  //- Helper function for init_facet_data()
  //- Associate facet edges in the interior of the passed facet patch
  //- with the passed list of partitions of a real, non-manifold curve.
  //- NOTE: Assumes edges internal to the facet list have been marked
  //-       with a '2' by the caller.


  std::set<CubitFacetEdge*> boundary_set, interior_set;
  PartitionSurface *const mySurface;
};

inline PartSurfFacetTool::PartSurfFacetTool(PartitionSurface* surface)
  : mySurface(surface) {}

#endif
