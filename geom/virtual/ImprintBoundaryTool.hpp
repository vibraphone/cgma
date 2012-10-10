//---------------------------------------------------------------------------
//- Class:          ImprintBoundaryTool
//- Description:    Given two surfaces, imprints the boundary loops on each
//-                 other.  Returns the intersection graph.
//-
//- Owner:          David R. White
//- Created:        4/15/2002
//---------------------------------------------------------------------------
#ifndef IMPRINT_BOUNDARY_TOOL_HPP
#define IMPRINT_BOUNDARY_TOOL_HPP

#include "GeometryDefines.h"
#include "AbstractTree.hpp"
#include "CubitDefines.h"
#include "DLIList.hpp"
class CubitVector;
class RefEntity;
class RefFace;
class RefEdge;
class RefVertex;
class ImprintPointData;
class ImprintLineSegment;
class ImprintMatchData;
template <class Y> class KDDTree;

template <class X> class DLIList;
typedef DLIList <ImprintPointData*> PointList;
typedef DLIList <PointList*> PointLoopList;
typedef DLIList <ImprintLineSegment*> SegList;
typedef DLIList <SegList*> SegLoopList;

enum IntersectResult { UNSET_INTERSECT = -1,
                       NO_INTERSECT = 0, 
                       CROSS_INTERSECT,//Segments intersect like a cross or plus sign, not near end points.
                       T_INTERSECT_0, //tells which node intersects of segments {0,1} and {2,3}
                       T_INTERSECT_1,
                       T_INTERSECT_2,
                       T_INTERSECT_3,
                       L_INTERSECT_0_2, // Segments intersect like an L or just at one end.
                       L_INTERSECT_0_3,// Numbers denote which points of segments intersect of segments {0,1} and {2,3}
                       L_INTERSECT_1_2,
                       L_INTERSECT_1_3,
                       OVERLAP_ALL_0_2_3_1, //One line fits entirely inside the other line.  The number denotes
                       OVERLAP_ALL_2_0_1_3, //which line is the super segment of segments {0,1} or {2,3}
                       OVERLAP_ALL_2_1_0_3, 
                       OVERLAP_ALL_0_3_2_1, 
                       OVERLAP_PART_0_3, //Lines overlap, number denote which points are exterior of overlap.
                       OVERLAP_PART_0_2,
                       OVERLAP_PART_1_2,
                       OVERLAP_PART_1_3,
                       OVERLAP_JOIN_02_1_3,//Lines overlap, but are joined at one end.
                       OVERLAP_JOIN_02_3_1,//First two digits denote nodes that are joined (or ontop of each other)
                       OVERLAP_JOIN_03_1_2,//Next digit denotes node on segment.
                       OVERLAP_JOIN_03_2_1,//Last digit denotes node at end that is not on anything.
                       OVERLAP_JOIN_13_0_2,
                       OVERLAP_JOIN_13_2_0,
                       OVERLAP_JOIN_12_0_3,
                       OVERLAP_JOIN_12_3_0,
                       SEGS_EQUAL_0_2, //Linesegments are equal and match according to number.
                       SEGS_EQUAL_0_3};

//  MatchType enum nomen-clature. Match points on segments 1 and 2.
//  But use the prev and next of each one also.  This is
//  for resolve small segments relative to myTolerance.
//     4---prev_seg_1--|0|------seg_1----|1|----next_seg_1---5
//
//     7---next_seg_2--|3|------seg_2----|2|----prev_seg_2---6 
//
// OR
//
//     4---prev_seg_1--|0|------seg_1----|1|----next_seg_1---5
//
//     6---prev_seg_2--|2|------seg_2----|3|----next_seg_2---7
enum MatchType { NO_MATCH,
                  MATCH_0_2,
                  MATCH_0_3,
                  MATCH_0_6,
                  MATCH_0_7,
                  MATCH_1_2,
                  MATCH_1_3,
                  MATCH_1_6,
                  MATCH_1_7,
                  MATCH_2_4,
                  MATCH_2_5,
                  MATCH_3_4,
                  MATCH_3_5,
                  ON_SEG_1,
                  ON_SEG_2
};
// Enum for returning during a loop.  Continue on, return or get out sucessful or
// report and error.
enum LoopEnum {RETURN_FAILURE, RETURN_SUCCESS, CONTINUE
};

class ImprintBoundaryTool
{
private:
  double myTolerance;
  RefFace *refFacePtr1;
  RefFace *refFacePtr2;
  CubitBoolean modBound1;
  CubitBoolean modBound2;
  
    //- Tolerance used for imprinting and snapping imprints together.
  PointList *allocatedPointData;
  SegList *allocatedLineData;
  PointLoopList *allocatedPointLoops;
  SegLoopList *allocatedLineLoops;
  DLIList <ImprintMatchData*> *allocatedMatchData;
  DLIList <RefEdge*> *allocatedRefEdge;
  

  CubitStatus imprint_boundaries(PointLoopList &boundary_loops_1,
                                 PointLoopList &boundary_loops_2,
                                 RefFace *ref_face_1,
                                 RefFace *ref_face_2,
                                 DLIList <RefFace*> &results);
    //- Intersects the two boundary point loops.  May delete point data or add point data to
    //- the loops.

  CubitStatus get_boundary_points( RefFace *ref_face,
                                   PointLoopList &boundary_point_loops );
    //- Gets the boundary points for the given surface.

  CubitStatus get_curve_facets( RefEdge* curve, PointList& segments );
    //- Gets the points that facet this curve.
  CubitStatus convert_to_lines( PointLoopList &boundary_point_loops,
                                SegLoopList &boundary_line_loops,
                                RefFace *ref_face,
                                const CubitBoolean boundary_1);
  CubitStatus imprint_segments( SegLoopList &boundary_line_loops_1,
                                SegLoopList &boundary_line_loops_2,
                                PointLoopList &boundary_loops_1,
                                PointLoopList &boundary_loops_2);
    //- Does the imprinting of the segments through finding the closest
    //- segments and closest points to match.  Utilizes the AbstractTree for
    //- efficiency. 

  CubitStatus find_matching_points( SegLoopList &seg_loops,
                                    AbstractTree <ImprintLineSegment*>& atree );
    //- Finds the matching points for a segments with those segments stored
    //- in the AbstractTree.  The seg_loops contain the segments for one surface
    //- and the AbstractTree contains them for the other surface.
  
  void set_closest_point(ImprintPointData *imp_point_0);
    //- Given the point, and the match data that has been added to it, find
    //- the point that is closest, and set that for the point.  This is
    //- in preparation to actually setting the matching point...

  CubitStatus find_closest_points(ImprintLineSegment *curr_seg,
                                  ImprintLineSegment *test_seg );
    //- For the two segments.  Finds the closest points between the
    //- segments for curr_seg.  (Doesn't do anything for test_seg.)

  CubitStatus match_seg_points( ImprintLineSegment *seg_1,
                                ImprintLineSegment *seg_2,
                                ImprintMatchData *&point_0_match,
                                ImprintMatchData *&point_1_match );
    //- Finds the match information between seg_1 and seg_2.  The
    //- ImprintMatchData returned stores how the points for seg_1
    //- line up with seg_2 (do they match with points on seg2, do
    //- they match along seg_2, do they not match at all...).

  CubitStatus final_match( ImprintLineSegment **surf_1_loop_heads,
                           int num_loops_1,
                           ImprintLineSegment **surf_2_loop_heads,
                           int num_loops_2,
                           AbstractTree <ImprintLineSegment*>& atree_1,
                           AbstractTree <ImprintLineSegment*>& atree_2);
    //- Finally sets the matching point information for both of the loops.
    //- Makes sure that if a point matches one boundary with a point
    //- on the other boundary, that that point also matches with it.
    //- Does this for both surfaces.

  CubitStatus final_match_loop(ImprintLineSegment **surf_1_loop_heads,
                               int num_loops_1,
                               ImprintLineSegment **surf_2_loop_heads,
                               int num_loops_2,
                               AbstractTree <ImprintLineSegment*>& atree_1,
                               AbstractTree <ImprintLineSegment*>& atree_2);
    //- Resolve the matching points for loop_1, with the points in surface 2.
    //- If the match is determined to be on a segment, update the AbstractTree by
    //- removing split segments and adding the new segments.

  CubitStatus match_on_point(ImprintPointData *start_point,
                             ImprintMatchData *start_match,
                             ImprintLineSegment **loop_heads,
                             int num_loops,
                             AbstractTree <ImprintLineSegment*>& atree_2);
    //- Match the start_point with the point stored in start match.  Set
    //- the matching point data.  This function may call match_on_segment if
    //- a point is matched up incorrectly.

  CubitStatus match_on_segment(ImprintPointData *start_point,
                               ImprintMatchData *start_match,
                               ImprintLineSegment **loop_heads,
                               int num_loops,
                               AbstractTree <ImprintLineSegment*>& atree_2);
    //- Match the start point with the segment stored in startmatch.  Split
    //- the segment and update the loops heads and AbstractTree accordingly.

  CubitStatus update_seg_matches(ImprintLineSegment *old_seg,
                                 ImprintLineSegment *new_seg_1,
                                 ImprintLineSegment *new_seg_2,
                                 ImprintMatchData *curr_match_data);
    //- As the segments get split, match data points to it.  Update
    //- the match data that may point to the old segment to point to
    //- the new_seg_1 or new_seg_2.  curr_match_data is not updated
    //- since it was already used (passed in to ignore).

  void update_linked_list(ImprintLineSegment *old_seg,
                          ImprintLineSegment *new_seg_1,
                          ImprintLineSegment *new_seg_2 );
    //- With new_seg_1 and new_seg_2 replacing old_seg, add them
    //- in the proper places to the segment linked list.

  void just_match_two_points(ImprintPointData *point_1,
                             ImprintPointData *point_2 );
    //- Does the low level data assignments for matching point_1 and point_2.

  CubitStatus find_crossings(ImprintLineSegment **surf_1_loop_heads,
                             int num_loops_1,
                             ImprintLineSegment **surf_2_loop_heads,
                             int num_loops_2,
                             AbstractTree <ImprintLineSegment*>& atree );
    //- Finds the crossings (segments that don't having any matching
    //- end points, but cross each other.).  Splits the segments that
    //- cross and updates the data lists etc.  Uses the atree of surface2.

  CubitStatus do_crossing(  ImprintPointData *imp_point_0,
                            ImprintPointData *imp_point_1,
                            ImprintPointData *imp_point_2,
                            ImprintPointData *imp_point_3,
                            ImprintLineSegment *&seg_1,
                            ImprintLineSegment *seg_2,
                            ImprintLineSegment **surf_1_loop_heads,
                            int num_loops_1,
                            ImprintLineSegment **surf_2_loop_heads,
                            int num_loops_2,
                            AbstractTree <ImprintLineSegment*>& atree,
                            CubitBoolean &modified);

    //- Computes the intersection of the two segments and does tha actual
    //- Crossing data updates.

  void update_boundary_loops( PointLoopList &boundary_loops,
                              ImprintLineSegment **surf_loop_heads);
    //- Updates the boundary_loops or rather replaces them with the
    //- data stored in the line segment linked lists.  Remember all
    //- the allocated memory is tracked by the class and deleted for
    //- the class so no memory is cleaned at this time.
  
    CubitStatus intersect_segments( ImprintLineSegment *seg_1,
                                  ImprintLineSegment *seg_2,
                                  IntersectResult &int_result,
                                  ImprintLineSegment **new_segments);
    //- Finds the intersection between seg_1 and seg_2.  The possible intersections
    //- are enumerated in the IntersectResult enum and the function assigns
    //- the appropriate result to int_result.  Also, if the intersection
    //- status requires it, the line segments are broken up into the new_segments
    //- array.  This array is assumed to be previously allocated and to
    //- be of size 4.  The first two positioins in the array are reserved for
    //- breaking up segment 1, the last two are for segment two.  They are
    //- assumed to be properly ordered (the first positions line segment should
    //- have the same start point as seg_1).  There is one exception to this.  In
    //- the case of OVERLAP_ALL_*, one line segment fits entirely inside another
    //- line segment.  The result is that one line segment gets split into *3* segments
    //- while the other is not split at all.  In that case, the new_segments array,
    //- will have three entries.  If the first segment is split into three, the new
    //- segments will occupy positions 0, 1, and 2. If it is the second, the new segments
    //- will occupy positions, 1, 2, and 3.  Again these are also properly ordered.
  
  void update_list(ImprintLineSegment **new_segments,
                   ImprintLineSegment *seg_1,
                   ImprintLineSegment *seg_2,
                   const CubitBoolean list_1 );
    //- utility function to update the linked list structure of the
    //- new segments.  List 1 flag tells which list to update.

  CubitStatus determine_on_boundary( ImprintLineSegment *seg_1,
                                     ImprintLineSegment *seg_2,
                                     MatchType &type_0, MatchType &type_1,
                                     MatchType &type_2, MatchType &type_3 );
    //- Checks the nodes to find if they are on seg_1 or 2.  They also
    //- determine based on the MatchType if the nodes should be matched
    //- with the segment or not.
  
  CubitStatus case_0_2_equal( ImprintLineSegment *seg_1,
                              ImprintLineSegment *seg_2,
                              ImprintPointData* imp_point_0,
                              ImprintPointData* imp_point_1,
                              ImprintPointData* imp_point_2,
                              ImprintPointData* imp_point_3,
                              MatchType type_1, MatchType type_3,
                              IntersectResult &int_result,
                              ImprintLineSegment **new_segments);
    //- This function resolves the intersections cases where points 0 and
    //- 2 are known to be equal and the other two nodes are not within tolerance.
    //- Basically two types of intersections can be occuring.  Either one of two
    //- OVERLAP_JOIN cases or an L_INTERSECT case.
  
  CubitStatus case_0_3_equal( ImprintLineSegment *seg_1,
                              ImprintLineSegment *seg_2,
                              ImprintPointData* imp_point_0,
                              ImprintPointData* imp_point_1,
                              ImprintPointData* imp_point_2,
                              ImprintPointData* imp_point_3,
                              MatchType type_1, MatchType type_2,
                              IntersectResult &int_result,
                              ImprintLineSegment **new_segments);
    //- This function resolves the intersections cases where points 0 and
    //- 2 are known to be equal and the other two nodes are not within tolerance.
    //- Basically two types of intersections can be occuring.  Either one of two
    //- OVERLAP_JOIN cases or an L_INTERSECT case.
  
  CubitStatus case_1_2_equal( ImprintLineSegment *seg_1,
                              ImprintLineSegment *seg_2,
                              ImprintPointData* imp_point_0,
                              ImprintPointData* imp_point_1,
                              ImprintPointData* imp_point_2,
                              ImprintPointData* imp_point_3,
                              MatchType type_0, MatchType type_3,
                              IntersectResult &int_result,
                              ImprintLineSegment **new_segments);
    //- This function resolves the intersections cases where points 0 and
    //- 2 are known to be equal and the other two nodes are not within tolerance.
    //- Basically two types of intersections can be occuring.  Either one of two
    //- OVERLAP_JOIN cases or an L_INTERSECT case.

  CubitStatus case_1_3_equal( ImprintLineSegment *seg_1,
                              ImprintLineSegment *seg_2,
                              ImprintPointData* imp_point_0,
                              ImprintPointData* imp_point_1,
                              ImprintPointData* imp_point_2,
                              ImprintPointData* imp_point_3,
                              MatchType type_0, MatchType type_2,
                              IntersectResult &int_result,
                              ImprintLineSegment **new_segments);
    //- This function resolves the intersections cases where points 0 and
    //- 2 are known to be equal and the other two nodes are not within tolerance.
    //- Basically two types of intersections can be occuring.  Either one of two
    //- OVERLAP_JOIN cases or an L_INTERSECT case.

  CubitStatus set_type_for_equal(ImprintPointData *pair_1_1,
                                 ImprintPointData *pair_1_2,
                                 ImprintPointData *pair_2_1,
                                 ImprintPointData *pair_2_2);
    //- Given the two pairs of points that are equal on the
    //- respective boundaries (_1_1 is pair 1, point on boundary 1,
    //- _1_2 is pair 1, point on boundary 2, _2_1 is pair 2, point
    //- on boundary 2, etc...), set up the correct PointTypes for
    //- the nodes for the EQUAL case.

  CubitStatus set_type_for_L(ImprintPointData *imp_point_1,
                             ImprintPointData *imp_point_2);
    //- Given point_1 on boundary 1 and point_2 on boundary 2, set
    //- the correct point type for each point given the two points
    //- meet at an L condition.
  
  CubitBoolean on_interior_segment(ImprintPointData *point,
                          ImprintLineSegment *line );
    //- Determines if the point is within myTolerance of the line segment.

  CubitStatus closest_point_seg( ImprintPointData *point,
                                 ImprintLineSegment *line,
                                 CubitVector &closest_point );
    //- Finds the closest point to the line segment.

  CubitStatus closest_point_seg( CubitVector &point_v,
                                 ImprintLineSegment *line,
                                 CubitVector &closest_point );
    //- Finds the closest point to the line segment.
  CubitBoolean closest_point_interior_seg( ImprintPointData *point,
                                           ImprintLineSegment *line,
                                           CubitVector &closest_point );
    //- Finds the closest point to the line segment.

  CubitStatus find_graph_for_surf(PointLoopList &boundary_loops_1,
                                  PointLoopList &boundary_loops_2,
                                  RefFace *ref_face,
                                  PointLoopList &part_segs,
                                  PointList &partition_points,
                                  CubitBoolean surf_1);
    //-Given the boundary_line_loops for a surface
    //-that have already been intersected against another
    //-surface, and the boundary_line_loops for that surface,
    //-determine the line segments that are needed
    //-to perform the actual imprint.  Also find points
    //-That are need to partition or imprint the existing boundaries.
    //-These points will also include the points of the segments
    //-that touch against the boundary, if necessary.  The
    //-RefFace is pretty important as it is used to determine
    //-if the points are inside or outside the surface.
  
  CubitStatus point_intersects_case( ImprintPointData *curr_point,
                                     ImprintPointData *next_point,
                                     RefFace *ref_face,
                                     PointLoopList &part_segs,
                                     PointList *&new_part_line,
                                     CubitBoolean &start_recording,
                                     CubitBoolean surf_1 );
    //-Given the fact that curr_point's PointType shows
    //-that it is a vertex on both boundaries or a new vertex,
    //-Resolve how the partition curves are extraced from the
    //-boundary loop 2 on surface 1.

  CubitBoolean on_surface( ImprintPointData *point,
                           RefFace *ref_face );
    //-Determines through the find_closest_point_trimmed
    //-function if a point is inside or outside a surface.
    //-This is not a very efficient method but since the
    //-surface can be non-planar, we can't do a simple test.

  CubitBoolean on_surface( CubitVector &vert_point,
                           RefFace *ref_face );
    //-Determines through the find_closest_point_trimmed
    //-function if a point is inside or outside a surface.
    //-This is not a very efficient method but since the
    //-surface can be non-planar, we can't do a simple test.

  CubitBoolean on_curve( ImprintPointData *point,
                         RefEdge *ref_edge );

    //- Determine if the point is within tolerance of the ref_edge.
  
  CubitBoolean are_connected( RefEntity *ent_1,
                              RefEntity *ent_2,
                              RefFace *ref_face );
    //-Determines through topology traversals
    //-if the two refentities are directly
    //-connected.  The function assumes
    //-that the RefEntities are either RefVerticies or
    //-RefEdges.  The entites must be connected within
    //-one entity for vertices or the same for
    //-curves...

  CubitStatus imprint_boundary_vertices( PointList &partition_points,
                                         CubitBoolean &modified_boundary );
    //-Given the data in the partition_points, partition the curves that the data points
    //-point to as their owner's.  Reasign owners as RefVertices.
  
  CubitStatus imprint_surface(RefFace *ref_face,
                              PointLoopList &part_segs,
                              DLIList <RefFace*> &results);
    //- Given the sorted partition segments, create refedges, then partition the
    //- surface.  Also merge appropriate vertices.  Relies on owner data, and matching
    //- data in the ImprintPointData class.

  CubitBoolean valid_partition(DLIList <RefEdge*> &ref_edges, RefFace *ref_face );
    //- Tests to see if the ref_edges are a valid parition, are there hard lines or unmerged
    //- lines in the lists of ref_edges.

  CubitStatus merge_vertices(DLIList <RefVertex*> &ref_verts);
    //- Do the normal merging but use myTolerance as the merge tolerance.
  
  CubitStatus merge_vertices(RefVertex *ref_vertex1,
                             RefVertex *ref_vertex2,
                             CubitBoolean &kept_1);
    //- Merge the two vertices.  Do a force for speed.
    //- Assigns kept_1 if ref_vertex1 is the vertex that survives the
    //- merge.

  CubitStatus resolve_on_boundaries( PointLoopList &boundary_loops_1,
                                     PointLoopList &boundary_loops_2);
  
  CubitBoolean resolve_match_conflict_this(ImprintPointData *this_point,
                                           ImprintPointData *other,
                                           CubitBoolean this_is_1 );
  CubitBoolean resolve_match_conflict_other(ImprintPointData *this_point,
                                           ImprintPointData *other,
                                           CubitBoolean this_is_1 );
  CubitStatus match_points( ImprintLineSegment *seg_1,
                            ImprintLineSegment *seg_2,
                            MatchType &type_0,
                            MatchType &type_1,
                            MatchType &type_2,
                            MatchType &type_3);
    //Computes how the end points between the segments match.  If
    //The segments are small, points may match with previous or next nodes
    //along the boundary and are indicated by the MatchType.
  
  void draw_seg(ImprintLineSegment *seg, int color = -1);
    //- draws the end node colored according to it's point type.
  void draw_point(ImprintPointData *imp_point, int color = -1);
    //- draws the end node colored according to it's point type.
  void draw_end(ImprintLineSegment *seg);
    //- draws the end node colored according to it's point type.
  void draw_loops(PointLoopList &boundary_loop);
    //- draws the loops of nodes.

  int num_coedges_on_face(RefEdge *edge_ptr,
                          RefFace *ref_face);
    //- Counts number of coedges this edge_ptr has that are also part of
    //- just this refface.
  

  CubitStatus ignore_match( ImprintPointData *start_point,
                            ImprintMatchData *start_match,
                            AbstractTree <ImprintLineSegment*>& atree_1,
                            CubitBoolean &ignore);
    ///
    /// Determine if the match should be ingored or not (if creating
    /// the match would make the loop self intersect.)
    ///

  RefEdge* create_virtual_edge(RefVertex *start,
                               RefVertex *end,
                               DLIList<CubitVector*> &interior_points);
  RefVertex* create_virtual_vertex(CubitVector &pos);
  


  
public:
  ImprintBoundaryTool(RefFace *ref_face_1,
                      RefFace *ref_face_2,
                      double tol = GEOMETRY_RESABS);
  ~ImprintBoundaryTool();
  CubitBoolean modified_bound_1()
    {return modBound1;}
  CubitBoolean modified_bound_2()
    {return modBound2;}
      
  CubitStatus imprint(DLIList <RefFace*> &results,
                      CubitBoolean merge = CUBIT_FALSE);
    //- Main interface function for imprinting the boundary loops
    //- of two refFaces.  It is assumed that if the reffaces intersect
    //- at any where else other than the boundary loops, this will
    //- not be found.
  
  
};

#endif
