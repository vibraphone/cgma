//-------------------------------------------------------------------------
// Filename      : FacetProjectTool.hpp
//
// Purpose       : Project curve segments onto a facetted surface.
//
// Special Notes : Used for surface partitioning.
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/08/02
//-------------------------------------------------------------------------

#ifndef FACET_PROJECT_TOOL_HPP
#define FACET_PROJECT_TOOL_HPP

#include <vector>

#include "CubitDefines.h"
#include "DLIList.hpp"

class PST_Edge;
class PST_CoEdge;
class PST_Face;
class PST_Point;
class PST_Entity;

class CubitVector;

/**
  * Project a polyline onto a manifold patch of (triangular) facets.
  *
  * @see FacetProjectTool::project
  */
class FacetProjectTool
{
  public:

  FacetProjectTool();
  ~FacetProjectTool();
  
  /**
   * Project a polyline onto a manifold patch of (triangular) facets.
   *
   * @param segments    Input list of ordered polyline vertices.
   * @param coordinates Input list of facet vertex coordinate triples 
   *                     ordered as { {x,y,z}, {x,y,z}, ...}
   * @param connections Input list of facets specified as { {p1,p2,p3}, ... } 
   *                     where each value p in the list corresponds to the 
   *                     vertex at coordinates[3*p]
   * @param duddedFacets Output list of facets that were modified during the operation
   *                     where each value i in the list corresponds to the
   *                     facet at connections[3*i]
   * @param newFacets   Output list of new facets replacing those specified in
   *                     duddedFacets where each value p corresponds to
   *                     the vertex at coordinates[3*p] or the vertex at
   *                     newPoints[3*p - coordinates.size()].
   * @param newFacetsIndex This list is populated with indices specifying
   *                     which facets in newFacets replace each facet in 
   *                     duddedFacets.  The facets in newFacets are 
   *                     grouped by the dudded facet they replace.  For 
   *                     each dudded facet specified at duddedFacets[i],
   *                     the new facets replacing that facet begin at
   *                     newFacets[a] and end at newFacets[b-1] where
   *                     a = newFacetsIndex[i-1] and b = newFacetsIndex[i].
   *                     Thus this list will always be the same size as
   *                     duddedFacets and the last value will always be
   *                     newFacets.size().
   * @param newPoints   Output list of new facet vertices created during 
   *                     the operation ordered as { {x,y,z}, {x,y,z}, ...}
   * @param edges       Output list of facet vertices specifying the projection
   *                     of the passed polyline onto the facets. The list may 
   *                     contain multiple polylines if the projection was 
   *                     self-intersecting or intersected the boundary of the 
   *                     facet patch in more than two places.  The format is:
   *                     {n1, {p1.1, p1.2, .., p1.n1}, n2, {p2.1, p2.2, .., p2.n2), ... 0 }
   *                     The list contains a count 'n' followed by n vertex 
   *                     indices.  The last value is always zero.  Each vertex
   *                     index p corresponds the vertex at coordinates[3*p] or 
   *                     the vertex at newPoints[3*p - coordinates.size()].
   * @param segmentPoints The list of facet points corresponding to input
   *                     segment points.  Point indices are the specified in
   *                     the same fashion as those in 'edges'.  A value of
   *                     -1 indicates that the segment point did not project
   *                     into the interior of the bounded surface.
   * @return            A return value of CUBIT_FAILURE indicates an internal
   *                     error in FacetProjectTool.  If the passed polyline
   *                     does not project onto any part of the facet patch,
   *                     the return value will be CUBIT_SUCCESS and the
   *                     output lists will be empty.
   * @see FSInt
   */
  CubitStatus project(
    DLIList<CubitVector*>& segments, 	        // in
    const std::vector<double>& coordinates,	// in
    const std::vector<int>& connections,	  // in
    std::vector<int>& duddedFacets,		      // out
    std::vector<int>& newFacets, 		        // out
    std::vector<int>& newFacetsIndex,		    // out
    std::vector<double>& newPoints,		      // out
    std::vector<int>& edges, 		            // out
    std::vector<int>& segmentPoints,        // out
    const double *tolerance_length = NULL
  );
  
  void debug( DLIList<CubitVector*>& segments,
              std::vector<double>& coordinates,
              std::vector<int>& connections );
    
  protected:
  
  void cleanup();
  
  PST_Point* project( const CubitVector& point,
                      const double* tolerance_length = NULL );
    //- Project a single point into the surface facetting.
    //- Returns resulting facet point.  
    //- Used for creating point-curves on partition surfaces.
      
    //- Use get_result_set(..) to get the results of the projection.
  
  CubitStatus get_result_set( DLIList<PST_CoEdge*>& edges );
    //- Get one of multiple possible sections of projected curve
    //- from the above project() method.  Returns CUBIT_FAILURE
    //- when there are no more results available.
    
  CubitStatus do_projection( const DLIList<CubitVector*>& segments ,
                             CubitBoolean& point_changed,
                             const double *tolerance_length = NULL);
    //- Project a segmented curve onto the facetted surface.
    //- project(..) will return CUBIT_FAILURE if none of the
    //- segmentes projected onto the surface.
    //- Use get_result_set(..) to get the results of the projection.

  CubitStatus project( PST_Point* &start, PST_Point* &end );
    //- Given two points within the facet patch, project the
    //- connecting segment onto the facets.  Assumes the
    //- projected segment will not intersect the boundary of
    //- the patch anywhere except the passed points.
  
  PST_Face* closest_face( const CubitVector& position );
    //- Find the facet for which the projection of the
    //- passed position onto the facet lies within the
    //- boundary of the facet and is closest to the passed
    //- position.  May return NULL if the position projects
    //- into no facets.
  
  PST_Edge* closest_edge( const CubitVector& position );
    //- Find the edge for which the closest point on the line
    //- of the edge to the passed position lies within the
    //- bounded edge and is closest to the passed position.
    //- May return NULL if the projection of the position 
    //- lies outside of all bounded edges.
    
  PST_Point* closest_point( const CubitVector& position );
    //- Find the point closest to the passed position.  Does
    //- not check points on the boundary of the facet patch.
    //- Will return null if all points are on the boundary.
  
  static bool project_to_face( PST_Face* face, 
                               const CubitVector& position,
                               CubitVector& result ) ;
    //- Project the pasesd position onto the passed face.  Returns
    //- false and does not change result if the projection of the
    //- point onto the plane of the face lies outside the face.
    
  static bool inside_face( PST_Face* face, const CubitVector& position );
  
  PST_Point* insert_in_face( PST_Face* face, 
                             const CubitVector& position );
    //- Insert a point in the interior of a triangle, splitting 
    //- the triangle into three triangles.
                             
  PST_Point* split_edge( PST_Edge* edge, double t );
    //- Split an edge, creating a new point and edge.  Adjacent
    //- faces are split with an edge from the new point to
    //- the point opposite the original edge on the face.
                             
  static
  void closest_on_lines( const CubitVector& base1,
                         const CubitVector& direction1,
			 const CubitVector& base2,
			 const CubitVector& direction2,
			 double& t1, double& t2 ); 
    //- Find the closest point on each of two skew lines to
    //- the other.  Result is returned as fraction of the 
    //- direction vector (from the base point) for each line.
    
  static
  double closest_on_line( const CubitVector& base,
                          const CubitVector& direction,
			  const CubitVector& from_position );
    //- Find the closest point on a line from the passed position.
    //- Returns result location as a multiple of the direction vector
    //- from the base point.
  
  CubitStatus populate_data( const DLIList<PST_Face*>& facets );
    //- Used to be the class constructor.  Fills in the private data lists 
    //- and sets marks on facet entities
  
  static
  bool within_tolerance ( PST_Edge* edge, const CubitVector& point, double tolerance );
  
  static
  CubitStatus next_around_point( PST_Point* point,
                                 const CubitVector& segment_end,
                                 PST_Face*& closest_face,
                                 PST_Edge*& closest_edge,
                                 CubitBoolean & is_boundary_edge,
                                 PST_Edge *last_closest_edge);
  
  private:
  
  DLIList<PST_Face*> facetList;
  DLIList<PST_Edge*> edgeList;
  DLIList<PST_Point*> pointList;
    //- Lists of all entities in the facet patch
  
  DLIList<PST_Point*> segPoints;
    //- List of facet points corresponding to input segment points.
    //- Some entries may be NULL.
    
  DLIList<PST_Edge*> boundaryList;
    //- List of edges on the boundary of the facet patch
    
  DLIList<PST_CoEdge*> imprintResult;
    //- List of resulting imprinted edges.
    
  void finalize_results();
};

#endif

