#include "FacetProjectTool.hpp"
#include "PST_Data.hpp"
#include "GeometryDefines.h"

#include "GfxDebug.hpp" /* for debugging output */

const int VG_FACET_DEBUG = 145;
const int VG_FACET_BOUNDARY_COLOR = CUBIT_RED;
const int VG_FACET_RAY_COLOR = CUBIT_BLUE;
const int VG_FACET_FACET_COLOR = CUBIT_GREEN;
const int VG_FACET_IMPRINT_COLOR = CUBIT_WHITE;
const int VG_FACET_SEGMENT_COLOR = CUBIT_CYAN;
#define VG_FACET_PRINT PRINT_DEBUG(VG_FACET_DEBUG)


//-------------------------------------------------------------------------
// Purpose       : Constructor
//
// Special Notes :  
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/11/02
//-------------------------------------------------------------------------
FacetProjectTool::FacetProjectTool( ) 
{ }


static int parent_compare_facets( PST_Face*& f1, PST_Face*& f2 )
{
  return (f1->parent < f2->parent) ? -1 :
         (f1->parent > f2->parent) ?  1 : 0;
}
static int sequence_compare_pts( PST_Point*& p1, PST_Point*& p2 )
{
  return (p1->sequence < p2->sequence) ? -1 :
         (p1->sequence > p2->sequence) ?  1 : 0;
}

//-------------------------------------------------------------------------
// Purpose       : takes a list of edge segments and a surface and returns
//                 the new facets, dudded facets, new points, and new edges
// Special Notes :  
//                 
//
// Creator       : John Fowler
//
// Creation Date : 12/03/02
//-------------------------------------------------------------------------
CubitStatus FacetProjectTool::project(
    DLIList<CubitVector*>& segments,            	// in
    const std::vector<double>& coordinates,		// in
    const std::vector<int>& connections,		// in
    std::vector<int>& duddedFacets, 			// out
    std::vector<int>& newFacets,			// out
    std::vector<int>& newFacetsIndex,		// out
    std::vector<double>& newPoints, 			// out
    std::vector<int>& edges,				// out
    std::vector<int>& segmentPoints,
    const double *tolerance_length
 )
{
  int i;
//   const double TOL_SQR = GEOMETRY_RESABS*GEOMETRY_RESABS; 
                     
  assert(connections.size()%3 == 0);  //  For triangles these must be triples
  assert(coordinates.size()%3 == 0);  //  Coordinates must come in triples   
  DLIList<PST_Edge*> facet_edges;

    // make the PST edges and faces from the coordinates and connections.
  PST_Edge::make_facets( coordinates, connections, GEOMETRY_RESABS, facet_edges );
  DLIList<PST_Face*> faces;
  PST_Edge::faces( facet_edges, faces );
  populate_data(faces);

    // Order the orignal points by sequence number.  New points
    // will be appended in order.
  pointList.sort(&sequence_compare_pts);
      
    // Do the work
  CubitBoolean point_changed;
  do_projection(segments, point_changed, tolerance_length);
  
    // fill in segmentPoints
  assert (segPoints.size() == segments.size());
  segmentPoints.resize( segPoints.size() );
  segPoints.reset();
  CubitStatus success = CUBIT_SUCCESS;
  for ( i = 0; i < segPoints.size(); i++ )
  {
    PST_Point* pt = segPoints.get_and_step();
    segmentPoints[i] = pt ? pt->sequence : -1;
    if (point_changed)
      success = CUBIT_FAILURE;
  }
  
  if (!success)
  {
     cleanup();
     return success;
  }
                                                                                
    //  Now put the points with sequence > coordinates.size()/3 in newPoints.
  int orig_num_points = coordinates.size()/3;
  int num_new_points = pointList.size() - orig_num_points;
  pointList.reset();
  pointList.step(orig_num_points);
  newPoints.resize(num_new_points*3);
  std::vector<double>::iterator ditor = newPoints.begin();
  while( num_new_points-- ) {
    PST_Point* pt = pointList.get_and_step();
    *ditor++ = pt->x();
    *ditor++ = pt->y();
    *ditor++ = pt->z();
  }

  if ( facetList.size() > 1 ) {  //  If only one facet, skip finding dudded facets
    // Sort (group) facets by parent.
    facetList.sort(&parent_compare_facets);

    // Fill in the duddedFacets, newFacets, and newFacetsIndex vectors.  

    int orig_num_facets = connections.size()/3;
    int new_facet_index = 0;
    DLIList<PST_Point*> facet_pts(3);

    facetList.reset();
    duddedFacets.clear();
    newFacetsIndex.clear();
    newFacets.clear();

    for ( i = facetList.size(); i > 0; ) {
      PST_Face* dudded_facet = 0;
      if( facetList.get()->parent != facetList.next()->parent ) {

          // unmodified original facet - skip it
        PST_Face* facet = facetList.get();
        assert(facet->sequence < orig_num_facets &&
               facet->sequence == facet->parent);
        facetList.step();
        i--;
      } 
      else {
          // new facets
        
          // put all new facets split from the
          // same original facet (having same parent),
          // including the orignal facet which was 
          // recycled, into newFacets
        PST_Face* facet = 0;
        do {
          facet = facetList.get_and_step();
          i--;
        
          if ( facet->parent == facet->sequence ) {
            assert(!dudded_facet);
            dudded_facet = facet;
          }
      
          facet_pts.clean_out();
          facet->append_points(facet_pts);
          facet_pts.reset();
          newFacets.push_back(facet_pts.get_and_step()->sequence);
          newFacets.push_back(facet_pts.get_and_step()->sequence);
          newFacets.push_back(facet_pts.get_and_step()->sequence);
          new_facet_index += 3;
        } while ( (facet->parent == facetList.get()->parent) && (i > 0) );
      
          // add replaced facet to duddedFacets
        assert(dudded_facet && dudded_facet->sequence < orig_num_facets);
        duddedFacets.push_back(dudded_facet->sequence);
          // add end position in new_facets for the
          // set of replacement facets to newFacetsIndex
        newFacetsIndex.push_back(new_facet_index);
      }
    }  
  }
  
  DLIList<PST_CoEdge*> coedge_list;
  edges.clear();
  finalize_results();
  while( get_result_set(coedge_list) ) {

      // add count and first point
    coedge_list.reset();
    edges.push_back(coedge_list.size() + 1);
    PST_Point* pt = coedge_list.get()->start_point();
    edges.push_back(pt->sequence);
    
      // add all other points
    for( i = coedge_list.size(); i--; ) {
      pt = coedge_list.get_and_step()->end_point();
      edges.push_back(pt->sequence);
    }
  
      // clean up for next iteration
    coedge_list.clean_out();
  }
  edges.push_back(0);  //  terminating zero for next segment length      
  
  cleanup();
  return CUBIT_SUCCESS;
}

 
//-------------------------------------------------------------------------
// Purpose       : fills in the data
//
// Special Notes : populate private data lists and set marks on facet 
//                 entities
//
// Creator       : Jason Kraftcheck, modified by John Fowler --
//		   used to be the class constructor
//
// Creation Date : 05/11/02, 12/03/02
//-------------------------------------------------------------------------
CubitStatus FacetProjectTool::populate_data( const DLIList<PST_Face*>& facets )
{
  facetList = facets;
  int i;
  
    // Mark all connected edges with a 3, and initialize 
    // parent number to sequence number.
  facetList.last();
  for( i = facetList.size(); i--; )
  {
    PST_Face* face = facetList.step_and_get();
    face->parent = face->sequence;
    PST_CoEdge* coe = face->first_coedge();
    do
    {
      coe->edge()->mark = 3;
      coe = coe->next();
    } while( coe != face->first_coedge() );
  }
  
    // Decrement the edge mark by one for each
    // face that contains the edge.
  for( i = facetList.size(); i--; )
  {
    PST_Face* face = facetList.step_and_get();
    PST_CoEdge* coe = face->first_coedge();
    do
    {
      coe->edge()->mark--;
      coe = coe->next();
    } while( coe != face->first_coedge() );
  }
  
    // Populate edgeList and boundaryList with edges
  for( i = facetList.size(); i--; )
  {
    PST_Face* face = facetList.step_and_get();
    PST_CoEdge* coe = face->first_coedge();
    do
    {
      PST_Edge* edge = coe->edge();
      	// If mark is 2, edge was only in one face in the
        // passed list, and is thus on the boundary of the
        // passed patch of faces.
      if( edge->mark == 2 )
      {
        boundaryList.append( edge );
        edgeList.append( edge );
        edge->mark = 0;
      }
      	// If mark is 1, the edge is an interior edge in the
        // patch of faces
      else if( edge->mark == 1 )
      {
        edgeList.append( edge );
        edge->mark = 0;
      }
      	// If the mark on the edge is zero, we have already
        // added it to the list(s).
      else
      {
        assert( edge->mark == 0 );  
      }
      
      coe = coe->next();
    } while( coe != face->first_coedge() );
  }
  
    // Get all the points
  PST_Edge::points( edgeList, pointList );
  
    // Clear marks on all faces connected to our edge
    // list (the passed faces and any faces connected
    // at the boundary edges.)
  for( i = edgeList.size(); i--; )
  {
    PST_Edge* edge = edgeList.step_and_get();
    edge->mark = 0;
    if( edge->forward()->face() )
      edge->forward()->face()->mark = 0;
    if( edge->reverse()->face() )
      edge->reverse()->face()->mark = 0;
  }

    // clear marks on faces
  for( i = facetList.size(); i--; )
    facetList.step_and_get()->mark = 0;
    // Set mark to 1 on boundary edges
  for( i = boundaryList.size(); i--; )
    boundaryList.step_and_get()->mark = 1;
    
  if( DEBUG_FLAG( VG_FACET_DEBUG ) )
  {
//    PST_Face::draw_faces( facetList, VG_FACET_FACET_COLOR, false );
    PST_Edge::debug_draw_edges( boundaryList, VG_FACET_BOUNDARY_COLOR );
  }

  return CUBIT_SUCCESS;
} 
 
//-------------------------------------------------------------------------
// Purpose       : Find the face for which the projection of the passed
//                 position onto the face is closest to the passed
//                 position. 
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/11/02
//-------------------------------------------------------------------------
PST_Face* FacetProjectTool::closest_face( const CubitVector& position )
{
  PST_Face* result = 0;
  double dist_sqr = 0;
  CubitVector projected;
  
  for( int i = facetList.size(); i--; )
  {
    PST_Face* facet = facetList.step_and_get();
    if( project_to_face( facet, position, projected ) )
    {
      double d = (projected - position).length_squared();
      if( !result || d < dist_sqr )
      {
        result = facet;
        dist_sqr = d;
      }
    }
  }
  
  return result;
}


//-------------------------------------------------------------------------
// Purpose       : Find the edge for which the projection of the passed
//                 position onto the edge is closest to the passed
//                 position.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/11/02
//-------------------------------------------------------------------------
PST_Edge* FacetProjectTool::closest_edge( const CubitVector& position )
{
  PST_Edge* result = 0;
  double dist_sqr = 0;
  
  for( int i = edgeList.size(); i--; )
  {
    PST_Edge* edge = edgeList.step_and_get();
    
    double t = edge->closest_on_line( position );
    if( (t > -CUBIT_RESABS) && (t < (1.0 + CUBIT_RESABS)) )
    {
      double d = (edge->position(t) - position).length_squared();
      if( !result || d < dist_sqr )
      {
        result = edge;
        dist_sqr = d;
      }
    }
  }

  return result;
}

//-------------------------------------------------------------------------
// Purpose       : Find the point closest to the passed position.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/11/02
//-------------------------------------------------------------------------
PST_Point* FacetProjectTool::closest_point( const CubitVector& position )
{
  PST_Point* result = 0;
  double dist_sqr = 0;
  
  for( int i = pointList.size(); i--; )
  {
    PST_Point* pt = pointList.step_and_get();
    double d = (position - *pt).length_squared();
    if( !result || d < dist_sqr )
    {
      dist_sqr = d;
      result = pt;
    }
  }
  
  return result;
}


//-------------------------------------------------------------------------
// Purpose       : Project the passed position onto the plane of the 
//                 passed facet.  Returns false if projection is outside
//                 the facet.
//
// Special Notes : "result" is unchanged if false is returned.
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/11/02
//-------------------------------------------------------------------------
bool FacetProjectTool::project_to_face( PST_Face* triangle,
                                        const CubitVector& position,
                                        CubitVector& result )
{
    // First check if projection of the point will be within
    // the bounds of the triangle.  
  if( ! inside_face( triangle, position ) )
    return false;
  
    // Calculate the actual projection
  result = triangle->plane().project( position );
  return true;
}

//-------------------------------------------------------------------------
// Purpose       : Destructor.  Clear marks on all facet entities.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/11/02
//-------------------------------------------------------------------------
FacetProjectTool::~FacetProjectTool()
{
  cleanup();
}


//-------------------------------------------------------------------------
// Purpose       : clean out internal data
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/13/03
//-------------------------------------------------------------------------
void FacetProjectTool::cleanup()
{
  while( pointList.size() )
    delete pointList.pop();
  facetList.clean_out();
  edgeList.clean_out();
  boundaryList.clean_out();
  imprintResult.clean_out();
  segPoints.clean_out();
}


//-------------------------------------------------------------------------
// Purpose       : Inset a point in the interior of a facet.
//
// Special Notes : Assumes input facet is a triangle, and splits facet
//                 as necessary to ensure that resulting facets are tris.
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/11/02
//-------------------------------------------------------------------------
PST_Point* FacetProjectTool::insert_in_face( PST_Face* face, 
                                             const CubitVector& position )
{
  PST_Point* new_point = new PST_Point( position );
  PST_CoEdge* ce1 = face->first_coedge();
  PST_CoEdge* ce2 = ce1->next();
  PST_CoEdge* ce3 = ce2->next();

    // insert non-mainifold edge in face after ce1
  PST_Edge* edge1 = PST_Edge::insert_in_face( new_point, ce1 );
  if( ! edge1 )
  {
    delete new_point;
    return 0;
  }
  new_point->sequence = pointList.size();
  pointList.append(new_point);
  edgeList.append( edge1 );
  if( DEBUG_FLAG( VG_FACET_DEBUG ) )
    edge1->debug_draw( VG_FACET_FACET_COLOR );
  
    // connect non-manifold edge to opposite side of the face
    // to make manifold topology.
  PST_Edge* edge2 = PST_Edge::split_face( ce2, edge1->start_point() == new_point ?
                                               edge1->reverse() : edge1->forward() );   
  if( ! edge2 ) return new_point;
  edgeList.append( edge2 );
  PST_Face* new_face = edge2->other(face);
  new_face->sequence = facetList.size();
  new_face->parent = face->parent;
  facetList.append( new_face );
  if( face->owner() )
    face->owner()->notify_split( face, new_face );
    
  if( DEBUG_FLAG( VG_FACET_DEBUG ) )
    edge2->debug_draw( VG_FACET_FACET_COLOR );
  
    // Split face so that all faces are triangles.
  PST_Face* split_face = ce3->face();
  PST_Edge* edge3 = PST_Edge::split_face( new_point, ce3->end_point(), split_face );
  if( edge3 )
  {
    edgeList.append( edge3 );
    new_face = edge3->other(split_face);
    new_face->sequence = facetList.size();
    new_face->parent = face->parent;
    facetList.append( new_face );
    if( split_face->owner() )
      split_face->owner()->notify_split( split_face, new_face );
  }
  
  return new_point;
}


//-------------------------------------------------------------------------
// Purpose       : Add a point where the passed position is projected
//                 into the facet patch.
//
// Special Notes : May return NULL if projection does not lie within 
//                 facet patch.
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/11/02
//-------------------------------------------------------------------------
PST_Point* FacetProjectTool::project( const CubitVector& pos,
                                      const double* tolerance_length )
{
  PST_Point* point = closest_point( pos );
  PST_Edge *  edge = closest_edge ( pos );
  PST_Face *  face = closest_face ( pos );
  
    // find closest facet
  double face_dist = CUBIT_DBL_MAX;
  CubitVector face_pos;
  if( face )
  {
    project_to_face( face, pos, face_pos );
    face_dist = (pos - face_pos).length_squared();
  }
  
    // find closest edge
  double edge_dist = CUBIT_DBL_MAX;
  double edge_pos  = CUBIT_DBL_MAX;
  if( edge )
  {
    edge_pos = edge->closest_on_edge( pos );
    edge_dist = (pos - edge->position(edge_pos)).length_squared();
  }
  
    // find closest point
  double point_dist = CUBIT_DBL_MAX;
  if( point )
  {
    point_dist = (pos - *point).length_squared();
  }
  
    // if the position is closer to a facet than a point or edge...
  if( face && (face_dist < point_dist) && (face_dist < edge_dist) )
  {
      // get a tolerance based on a sample face size. I attempted
      // to do this with a representative facet size but still 
      // ended up with problems.  While this is slower is more
      // robust.  KGM
    double facet_length;

    if(tolerance_length)
    {
      facet_length = *tolerance_length;
    }
    else
    {
      facet_length = face->bounding_length();
    }

    // use the square of the facet length times an arbitrary factor to
    // determine the tolerance
    double tol_sqr = facet_length*facet_length*.00001;

      // check if projected position is within tolerance of
      // any bounding edge of the face.
    PST_CoEdge* coe = face->first_coedge();
    bool insert = true;
    do {
      PST_Edge* e = coe->edge();
      CubitVector p = e->position( e->closest_on_line( face_pos ) );
      double d = (p - face_pos).length_squared();

      // use a relative tolerance here - KGM
      if( d < tol_sqr ) 
      {
        edge = e;
        edge_pos = edge->closest_on_edge( pos );
        face_dist = d;
        edge_dist = d;
        insert = false;
      }
      coe = coe->next();
    } while( coe != face->first_coedge() );
    
      // If projected position was sufficiently far from the edges
      // of the facet, insert an interior point in the facet.  Otherwise
      // fall through and split edge.
    if( insert )
    {
      return insert_in_face( face, face_pos );
    }
  }
  
    // If the point was closer to an edge than a point
  if( edge && (edge_dist <= point_dist) )
  {
      // If within tolerance of an end point of the edge, return
      // the end point.
    CubitVector pt = edge->position( edge_pos );
    if( (pt - *edge->start_point()).length_squared() < GEOMETRY_RESABS*GEOMETRY_RESABS )
      point = edge->start_point();
    else if( (pt - *edge->end_point()).length_squared() < GEOMETRY_RESABS*GEOMETRY_RESABS )
      point = edge->end_point();
      // Otherwise split the edge at the projected position
    else point = split_edge( edge, edge_pos );
  }
  
  if( DEBUG_FLAG( VG_FACET_DEBUG ) && point )
    point->debug_draw( VG_FACET_IMPRINT_COLOR );
  
    // If this was not set in the above if-block, then it is either
    // the closest point in the facet patch or NULL if the position
    // projected outside the facet patch.
  return point;
}    
  
//-------------------------------------------------------------------------
// Purpose       : Split an edge at the specified parameter value on the
//                 edge.  Splits connected facets to maintain triangular
//                 facets.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/11/02
//-------------------------------------------------------------------------
PST_Point* FacetProjectTool::split_edge( PST_Edge* edge, double t )
{
  assert( t > 0.0 && t < 1.0 );
  
    // create point
  PST_Point* point = new PST_Point( edge->position(t) );
  
    // get connected faces 
  PST_Face* face1 = edge->forward()->face();
  PST_Face* face2 = edge->reverse()->face();
  
    // point on each triangle opposite the edge to be split
  PST_Point* face1_pt = face1 ? edge->forward()->next()->end_point() : 0;
  PST_Point* face2_pt = face2 ? edge->reverse()->next()->end_point() : 0;

  PST_Edge* face1_edge = 0, *face2_edge = 0;
  
    // split edge edge
  PST_Edge* split_edge = edge->split( point );
  if( edge->owner() )
    edge->owner()->notify_split( edge, split_edge );
    // if we split a boundary edge, put the new edge in the
    // list of boundary edges too.
  if( edge->mark )
  {
    split_edge->mark = edge->mark;
    boundaryList.append( split_edge );
  }
    // add edge to our list of edges
  edgeList.append( split_edge );
    // If the first face was not a boundary face (different
    // than the boundary of our patch of facets), then split
    // it so that it remains a triangle.
  if( face1 ) 
  {
    face1_edge = PST_Edge::split_face( point, face1_pt, face1 );
    
    if( face1_edge )
    {
      if( DEBUG_FLAG( VG_FACET_DEBUG ) )
        face1_edge->debug_draw( VG_FACET_FACET_COLOR );

      edgeList.append( face1_edge );
	    PST_Face* new_face = face1_edge->other(face1);
      new_face->sequence = facetList.size();
      new_face->parent = face1->parent;
      facetList.append(new_face );
      
      if( face1->owner() )
        face1->owner()->notify_split( face1, face1_edge->other(face1) );
    }
  }
    // If the second face is not a boundary face, split it
    // so that we still have triangles.
  if( face2 )
  {
    face2_edge = PST_Edge::split_face( point, face2_pt, face2 );
    if( face2_edge )
    {
      if( DEBUG_FLAG( VG_FACET_DEBUG ) )
        face2_edge->debug_draw( VG_FACET_FACET_COLOR );

      edgeList.append( face2_edge );
	    PST_Face* new_face = face2_edge->other(face2);
      new_face->sequence = facetList.size();
      new_face->parent = face2->parent;
      facetList.append( new_face );
      
      if( face2->owner() )
        face2->owner()->notify_split( face2, face2_edge->other(face2) );
    }
  }
  point->sequence = pointList.size();
  pointList.append(point);
  
  return point;
}



//-------------------------------------------------------------------------
// Purpose       : This is the main or outer-loop function for the
//                 polyline projection.  It is called after the input data
//                 is converted to the working representation.  The passed
//                 list of CubitVectors is the polyline.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : ??
//-------------------------------------------------------------------------
CubitStatus 
FacetProjectTool::do_projection( const DLIList<CubitVector*>& passed_list, 
                                 CubitBoolean & point_changed,
                                 const double *tolerance_length)
{
  const double TOL_SQR = GEOMETRY_RESABS*GEOMETRY_RESABS;
  point_changed = CUBIT_FALSE;
  DLIList<CubitVector*> segments( passed_list );
  segments.reset();
  CubitVector end = *(segments.get_and_step());
  PST_Point* last_point = 0;
  
  for( int i = segments.size(); i > 1; i-- )
  {
    CubitVector start(end);
    end = *(segments.get_and_step());
    
    if( DEBUG_FLAG( VG_FACET_DEBUG ) )
    {
      GfxDebug::draw_line( 
      	(float)start.x(), (float)start.y(), (float)start.z(),
        (float)end.x()  , (float)end.y()  , (float)end.z(),
        VG_FACET_SEGMENT_COLOR );
      GfxDebug::flush();
      //GfxDebug::mouse_xforms();
    }

    
      // get projection of start point if we don't have it from
      // the previous segment.
    if( !last_point )
    {
      last_point = project( start, tolerance_length );
      segPoints.append(last_point); // last_point may be NULL.
    }
      
      // no projection for entire segment?
    if( !last_point  )
      continue;
      
      // create end point
    PST_Point* end_pt = project( end, tolerance_length );
    segPoints.append( end_pt );
    if ( !end_pt )
      continue;
    
    PST_Point* start_point = last_point;
    PST_Point* end_point = end_pt;
    project( last_point, end_pt );
    if ((*start_point - *last_point).length_squared() > TOL_SQR)
    {
       segPoints.move_to(start_point);
       segPoints.change_to(last_point);
       point_changed = CUBIT_TRUE;
    }
    if ((*end_pt - *end_point).length_squared() > TOL_SQR)
    {
       segPoints.move_to(end_point);
       segPoints.change_to(end_pt);
       point_changed = CUBIT_TRUE;
    }
    last_point = end_pt;
  
  } // end for( i = segments )
  
  return CUBIT_SUCCESS;
}


//-------------------------------------------------------------------------
// Purpose       : Given two existing points in a triangulation, add edges
//                 to the triangulation to connect those two points by the
//                 most direct path that does not leave the triangulation.
//
// Special Notes : Thisis the second step of the facet-polyline projection
//                 algorithm.  First each point of the polyline is 
//                 projected onto the facetting.  Second this function is
//                 called to find/create the set of facet edges connecting
//                 those points.
//
// Creator       : Jason Kraftcheck
//
// Creation Date : ??
//-------------------------------------------------------------------------
CubitStatus FacetProjectTool::project( PST_Point* &start, PST_Point* &end )
{
  PST_Edge* edge,* closest_edge;
  PST_Face* closest_face;
  PST_Point* point;
  PST_Point* start_point = start;
#ifdef BOYD17
  DLIList<PST_Face*> face_list;
  DLIList<PST_Edge*> boundary_edges;
#endif
  CubitStatus stat;
  
  const double TOL_SQR = GEOMETRY_RESABS*GEOMETRY_RESABS;
  
  while( start_point != end )
  {
      // If the points share an edge, that edge is the result.
      // Return it.
    CubitBoolean is_boundary_edge;
    if( (edge = start_point->common( end )) )
    {
      PST_CoEdge* coedge = edge->end_point() == end ? 
                           edge->forward() : edge->reverse();
      imprintResult.append( coedge );

      if( DEBUG_FLAG( VG_FACET_DEBUG ) )
      {
      	edge->debug_draw( VG_FACET_IMPRINT_COLOR );
      }
	
      return CUBIT_SUCCESS;
    }
    
      // Draw the current segment
    if( DEBUG_FLAG( VG_FACET_DEBUG ) )
    {
      GfxDebug::draw_line( 
      	(float)start_point->x(), (float)start_point->y(), 
        (float)start_point->z(),
        (float)end->x()  , (float)end->y()  , (float)end->z(),
        VG_FACET_SEGMENT_COLOR );
      GfxDebug::flush();
      GfxDebug::mouse_xforms();
    }
    
      // Get face or edge adjacent to start and closest to the polyline segment
    stat = next_around_point( start_point, *end, closest_face, closest_edge, 
                              is_boundary_edge );
    if (!stat || (closest_face && closest_edge))
    {
      assert(false);
      return stat;
    }
    
    const CubitVector seg_tan = *end - *start_point;
    double t_seg = 0.0;
    point = 0;
    
    if (closest_face)
    {
        // calculate intersection with triangle edge opposite
        // 'start_point'.
      double t_edge;
      edge = closest_face->opposite(start_point);
      PST_Point* edge_start = edge->start_point();
      PST_Point* edge_end   = edge->end_point();
      closest_on_lines( *edge_start, *edge_end - *edge_start,
                        *start_point,  seg_tan, t_edge, t_seg);
      
      const CubitVector  seg_pt = *start_point + t_seg  *  seg_tan;
      const CubitVector edge_pt = (1.0 - t_edge) * *edge_start
                                       + t_edge  * *edge_end;
      
        // if line segment intersects opposite edge of triangle...
      if ( (t_seg < 1.0) || ((seg_pt - edge_pt).length_squared() < TOL_SQR) )
      {
          // Check if intersection is near start of opposite edge
        if ((t_edge <= 0.0) || 
            ((*edge_start - edge_pt).length_squared() < TOL_SQR) ||
            ((*edge_start - seg_pt ).length_squared() < TOL_SQR) )
          point = edge_start;
          // Check if intersection is near end of opposite edge
        else if ((t_edge >= 1.0) || 
                 ((*edge_end - edge_pt).length_squared() < TOL_SQR) ||
                 ((*edge_end -  seg_pt).length_squared() < TOL_SQR))
          point = edge_end;
          // Otherwise intersection is in the interior of the edge
        else
          point = split_edge(edge, t_edge);
      }
        // otherwise segment end is inside triangle
      else
      {
        t_seg = 1.0;
        const CubitVector proj = closest_face->plane().project(*end);
        PST_Edge *edge1, *edge2;
        closest_face->two_edges (start_point, edge1, edge2);
          // Too close to start position - skip segment
        if ( (*start_point - proj).length_squared() < TOL_SQR )
          point = start_point;
          // If close to side of triangle, fall through to edge
          // handling code
        else if (within_tolerance(edge1, *end, GEOMETRY_RESABS) ||
                 within_tolerance(edge1, proj, GEOMETRY_RESABS))
          closest_edge = edge1;
          // If close to side of triangle, fall through to edge
          // handling code
        else if (within_tolerance(edge2, *end, GEOMETRY_RESABS) ||
                 within_tolerance(edge2, proj, GEOMETRY_RESABS))
          closest_edge = edge2;
          // Insert new point in triangle interior...
        else
          point = insert_in_face (closest_face, proj);
      }
    } // if(closest_face)
    
    if (closest_edge)
    {
      PST_Point *const edge_end = closest_edge->other(start_point);
        // If edge and segment ends are equal...
      if ( (*edge_end - *end).length_squared() < TOL_SQR )
      {
        point = edge_end;
        edge = closest_edge;
        if (is_boundary_edge)
          end = edge_end;
      }
        // Otherwise calculate closest point
      else
      {
          // edge and segment share a start point so just need
          // to calculate the projection of each on the other.
        const CubitVector edge_tan = *edge_end - *start_point;
        const double dot_prod = seg_tan % edge_tan;
        const double t_seg  = dot_prod /  seg_tan.length_squared();
        const double t_edge = dot_prod / edge_tan.length_squared();
        
        const CubitVector  seg_pt = *start_point + t_seg  *  seg_tan;
        const CubitVector edge_pt = *start_point + t_edge * edge_tan;
        
          // Too close to start point -- skip segment
        if (t_edge <= 0.0 || t_seg <= 0.0 ||
            (*start_point - edge_pt).length_squared() < TOL_SQR ||
            (*start_point -  seg_pt).length_squared() < TOL_SQR )
          point = start_point;
          // If segment end is near or passed edge end, 
          // then the edge end is the intersection with this triangle
        else if (t_edge > 1.0 || 
                 (*edge_end - edge_pt).length_squared() < TOL_SQR ||
                 (*edge_end -  seg_pt).length_squared() < TOL_SQR)
        {
          //make sure the facet edge is not a boundary edge.
          if (is_boundary_edge && start_point == start)
            start =  edge_end;
          point = edge_end;
        }
          // Otherwise closest point to segment end is in the interior
          // of the edge -- split the edge.
        else
          point = split_edge( closest_edge, t_edge );
      }
    } // if(closest_edge)

    if (point == start_point) // skip perpendicular or very short segment
      continue;

    edge = start_point->common( point );
    if (!edge)
    {
      assert(false);
      continue;
    }
    
    if ( !is_boundary_edge || start_point != start )
    { 
       PST_CoEdge* coedge = edge->end_point() == point ?
                            edge->forward() : edge->reverse();
       imprintResult.append(coedge);
    }

    if( DEBUG_FLAG( VG_FACET_DEBUG ) )
    {
      edge->debug_draw( VG_FACET_IMPRINT_COLOR );
      point->debug_draw( VG_FACET_IMPRINT_COLOR );
    }

    start_point = point;
  } // while(start_point != end)
  
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Check if a point is within the passed tolerance
//                 of a bounded edge.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/04/03
//-------------------------------------------------------------------------
bool FacetProjectTool::within_tolerance( PST_Edge* edge,
                                         const CubitVector& point,
                                         double tolerance )
{
  const double t_edge = edge->closest_on_edge( point );
  const CubitVector edge_pos = edge->position(t_edge);
  return (edge_pos - point).length_squared() < tolerance*tolerance;
}

//-------------------------------------------------------------------------
// Purpose       : Given a line segment beginning at a point in the
//                 triangulation, find the face or edge adjacent to the
//                 point and closest to the line segment.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/04/03
//-------------------------------------------------------------------------
CubitStatus FacetProjectTool::next_around_point ( PST_Point* start,
                                        const CubitVector& seg_end,
                                        PST_Face*& closest_face,
                                        PST_Edge*& closest_edge,
                                        CubitBoolean & is_boundary_edge )
{
  DLIList<PST_Face*> face_list;
  DLIList<PST_Edge*> boundary_edges;
  is_boundary_edge = CUBIT_FALSE; 
  double shortest_dist_sqr = CUBIT_DBL_MAX;
  closest_face = 0;
  closest_edge = 0;
  
  const CubitVector tangent = seg_end - *start;
   
    // To traverse the edges in clockwise order about the point,
    // it is necessary to start with the appropriate boundary edge
    // if there is one.  If this point is not on the boundary of 
    // the facet patch, any edge is acceptable.

    // First mark all faces adjacent to the vertex.  This is to make
    // sure disjoint sets of faces (sets of faces with no edges in
    // common) are not missed when traversing edge-face adjacencies.
  start->faces( &face_list );
  for (int i = face_list.size(); i--; )
    face_list.step_and_get()->mark = true;

    // Loop until all faces in face_list have been un-marked
  boundary_edges.clean_out();
  closest_edge = 0;
  closest_face = 0;
  while (face_list.size())
  {
    PST_Face* face = face_list.pop();
    if (!face->mark)
      continue;

    PST_CoEdge* coedge = face->first_coedge();
    if (coedge->edge()->other(start) == 0)
      coedge = coedge->next();

      // search for the 'first' coedge in the clockwise list.
    PST_CoEdge* first = coedge;
    while (coedge->face())
    {
      if (coedge->end_point() == start)
        coedge = coedge->other();
      else
        coedge = coedge->previous();

      if (coedge == first) // no start, point is not on boundary
        break;
    }

      // keep boundary edges encountered for later
    if (!coedge->face())
      boundary_edges.append(coedge->edge());

      // Of the two edges on the face that are adjacent to the
      // point, begin with the one that is first in clockwise
      // order around the point.
    if (coedge->start_point() == start)
      coedge = coedge->other();
    first = coedge;

      // Traverse adjacent facets in clockwise order
      // around point.
    bool clockwise_of_prev = false;
    while ( (face = coedge->face()) )
    {
      face->mark = false;

        // get vectors for edges of face pointing outwards
        // from vertex.
      assert(coedge->end_point() == start && coedge->next()->start_point() == start);
      CubitVector prev = coedge->other()->direction();
      CubitVector next = coedge->next()->direction();

        // calculate if segment is to the 'inside' of each
        // adjacent edge.
      CubitVector normal = next * prev;  // ccw 
      bool inside_prev = ((tangent * prev) % normal) >= 0.0;
      bool inside_next = ((next * tangent) % normal) >= 0.0;

        // If edge is inside face, check distance from face
      if (inside_prev && inside_next)
      {
          // calculate distance from end point to plane of face
        const double plane_coeff = normal % *start;
        const double end_coeff = normal % seg_end;
        const double numerator = end_coeff - plane_coeff;
        const double dist_sqr = (numerator * numerator) / normal.length_squared();

        if (dist_sqr < shortest_dist_sqr)
        {
          closest_face = face;
          closest_edge = 0;
          shortest_dist_sqr = dist_sqr;
        }

        clockwise_of_prev = false;
      }
        // If the edge is above a peak formed by the shared edge of two faces
      else if (inside_next && clockwise_of_prev)
      {
          // calculate distance from end of segment to line of facet edge
        const double dot_prod = tangent % prev;
        const double dist_sqr = 
          tangent.length_squared() - dot_prod * dot_prod / prev.length_squared();

        if (dist_sqr <= shortest_dist_sqr)
        {
          closest_face = 0;
          closest_edge = coedge->edge();
          shortest_dist_sqr = dist_sqr;
        }
        clockwise_of_prev = false;
      }
        // If clockwise of the face, note for next iteration
      else if (inside_prev)
      {
        clockwise_of_prev = true;
      }
      else
      {
        clockwise_of_prev = false;
      }  

      coedge = coedge->next()->other();
      if (coedge == first)
        break;
    } // while(coedge->face())
    
    if (!coedge->face())
      boundary_edges.append(coedge->edge());
  
  } // while(face_list.size())
  
    // If a closest entity was found, then done.
  if (closest_face || closest_edge)
    return CUBIT_SUCCESS;
    
    // Otherwise return closest boundary edge
  while (boundary_edges.size())
  {
    PST_Edge *const edge = boundary_edges.pop();
    const CubitVector edge_tan = *(edge->other(start)) - *start;
    const double dot_prod = tangent % edge_tan;
    if (dot_prod <= 0.0)
      continue;
    
    const double dist_sqr =
      tangent.length_squared() - dot_prod * dot_prod / edge_tan.length_squared();
    if (dist_sqr < shortest_dist_sqr)
    {
      closest_edge = edge;
      is_boundary_edge = CUBIT_TRUE;
      shortest_dist_sqr = dist_sqr;
    }
  }
  
  return closest_edge ? CUBIT_SUCCESS : CUBIT_FAILURE;
}

  
  
    
  

//-------------------------------------------------------------------------
// Purpose       : Remove duplciate edges from results if input polyline
//                 was self-intersecting.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : ??
//-------------------------------------------------------------------------
void FacetProjectTool::finalize_results()
{
  imprintResult.reset();
  for ( int i = imprintResult.size(); i--; )
  {
    PST_Edge* edge = imprintResult.step_and_get()->edge();
    if (edge->mark)
      imprintResult.change_to(0);
    else
      edge->mark = 1;
  }
  imprintResult.remove_all_with_value(0);
  imprintResult.reset();
}

//-------------------------------------------------------------------------
// Purpose       : Return non-intersecting chains of result edges, in the
//                 same sequence as the input polyline and with the same
//                 orientation.
//
// Special Notes : This is called after the projection is done to get the
//                 resulting facet edge chains.
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/04/03
//-------------------------------------------------------------------------
CubitStatus FacetProjectTool::get_result_set( DLIList<PST_CoEdge*>& coedges )
{
  if ( !imprintResult.size() )
    return CUBIT_FAILURE;
  
  imprintResult.reset();
  PST_CoEdge* prev = imprintResult.get();
  imprintResult.change_to(0);
  coedges.append(prev);
  
  for ( int i = imprintResult.size() - 1; i--; )
  {
    PST_CoEdge* coedge = imprintResult.step_and_get();
    PST_Point* pt = prev->end_point();
    if (pt != coedge->start_point())
      break;
    
    int count = 1;
    for (PST_Edge* pt_edge = pt->next(coedge->edge());
         pt_edge != coedge->edge();
         pt_edge = pt->next(pt_edge))
    {
      if (pt_edge->mark)
        count++;
    }
     
    if (count != 2)
      break;
    
    coedges.append(coedge);
    prev = coedge;
    imprintResult.change_to(0);
  }
  
  imprintResult.remove_all_with_value(0);
  return CUBIT_SUCCESS;
}


double FacetProjectTool::closest_on_line( const CubitVector& B,
                                      const CubitVector& M,
                                      const CubitVector& P ) 
{ 
  const double dist_tol_sqr = GEOMETRY_RESABS*GEOMETRY_RESABS;
  double D = M.length_squared();
  return (D < dist_tol_sqr ) ? 0.0 : (M % (P - B)) / D;
}


void FacetProjectTool::closest_on_lines( const CubitVector& B1,
                                     const CubitVector& M1,
                                     const CubitVector& B2,
                                     const CubitVector& M2,
                                     double& t1, double& t2 ) 
{
  CubitVector N = M2 * (M2 * M1);
  double d = N % M1;
  if( fabs(d) < CUBIT_RESABS ) //parallel lines
    t1 = 0.0;
  else
    t1 = ((N % B2) - (N % B1)) / d;

  t2 = closest_on_line( B2, M2, B1 + t1 * M1 );
}
  
      
bool FacetProjectTool::inside_face( PST_Face* face, const CubitVector& pos )
{      
    // This code checks if the position is in the interior of 
    // the prism formed by sweeping the facet infinetly in the
    // direction of, and opposite direction of the facet normal.
  PST_CoEdge* coe = face->first_coedge();
  do {
    CubitVector vect( pos );
    vect -= *coe->start_point();
    vect *= coe->direction();
    if( vect % face->normal() > -CUBIT_RESABS )
      return false;
    coe = coe->next();
  } while( coe != face->first_coedge() );
  return true;
}

/*
static void draw_pt( const std::vector<double>& list, int index, int color )
{
  GfxDebug::draw_point( 
    (float)list[index], (float)list[index+1], (float)list[index+2], color );
}

static void draw_edge( const std::vector<double>& list, int i1, int i2, int color )
{
  i1 *= 3; i2 *= 3;
  GfxDebug::draw_line(
    (float)list[i1], (float)list[i1+1], (float)list[i1+2],
    (float)list[i2], (float)list[i2+1], (float)list[i2+2],
    color );
}

static void draw_tri( const std::vector<double>& pts,
                      const std::vector<int>& tris, int index, 
                      const std::vector<int>& seq, int color )
{
  for ( int i = 0; i < 3; i++ )
  {
    int i1 = tris[index + i];
    int i2 = tris[index + (i + 1) % 3];
    if ( !seq[i1] || !seq[i2] || abs(seq[i1]-seq[i2]) != 1 )
      draw_edge( pts, i1, i2, color );
  }
}
*/

void FacetProjectTool::debug( DLIList<CubitVector*>& segments,
                              std::vector<double>& coordinates,
                              std::vector<int>& connections )
{
#if 1 /* test with old interface */  
  DLIList<PST_Edge*> edges;
  PST_Edge::make_facets(coordinates, connections, GEOMETRY_RESABS, edges);
  
  DLIList<PST_Face*> faces;
  PST_Edge::faces( edges, faces );
  populate_data(faces);
  
  CubitBoolean point_changed;
  do_projection( segments, point_changed );
  
  PST_Edge::debug_draw_edges( edgeList, CUBIT_BLUE );
  
  DLIList<PST_CoEdge*> coedges;
  finalize_results();
  while( get_result_set( coedges ) ) {
    edges.clean_out();
    coedges.reverse();
    while (coedges.size()) edges.append(coedges.pop()->edge());
    PST_Edge::debug_draw_edges( edges, CUBIT_WHITE );
    PST_Edge::debug_draw_points( edges, CUBIT_WHITE );
  }

#else /* test with new interface */

  int i, j;

    // do projection
  std::vector<double> new_pts;
  std::vector<int> dudded, new_facets, facet_index, polyline;
  if( !project( segments, coordinates, connections, dudded, 
                new_facets, facet_index, new_pts, polyline ) )
  {
    PRINT_ERROR("FacetProjectTool::project returned CUBIT_FAILURE\n");
    return;
  } 
  
  assert( connections.size() % 3 == 0 );
  assert( coordinates.size() % 3 == 0 );
  assert( new_facets.size() % 3 == 0 );
  assert( new_pts.size() % 3 == 0 );
  assert( dudded.size() == facet_index.size() );
  int num_old_tri = connections.size() / 3;
  int num_new_tri = new_facets.size() / 3;
  int num_old_pts = coordinates.size() / 3;
  int num_new_pts = new_pts.size() / 3;
  int num_tot_tri = num_new_tri + num_old_tri;
  int num_tot_pts = num_new_pts + num_old_pts;
  
  PRINT_INFO("%d input triangles using %d points.\n", num_old_tri, num_old_pts );
  PRINT_INFO("%d new triangles and %d new poitns.\n", num_new_tri, num_new_pts );
  PRINT_INFO("%d total triangles using %d total points.\n", num_tot_tri, num_tot_pts );
  PRINT_INFO("connections.size() = %d\n", (int)connections.size());
  PRINT_INFO("coordinates.size() = %d\n", (int)coordinates.size());
  PRINT_INFO("dudded.size() = %d\n", (int)dudded.size());
  PRINT_INFO("new_facets.size() = %d\n", (int)new_facets.size());
  PRINT_INFO("facet_index.size() = %d\n", (int)facet_index.size());
  PRINT_INFO("new_pts.size() = %d\n", (int)new_pts.size());
  PRINT_INFO("polyline.size() = %d\n", (int)polyline.size());
  
  for ( i = 0; i < (int)new_facets.size(); i++ )
    if ( new_facets[i] >= num_tot_pts || new_facets[i] < 0 )
      PRINT_ERROR("Invalid value %d at new_facets[%d]\n", new_facets[i], i);
  for ( i = 0; i < (int)dudded.size(); i++ )
    if ( dudded[i] >= num_old_tri || dudded[i] < 0 )
      PRINT_ERROR("Invalid value %d at dudded[%d]\n", dudded[i], i);
  for ( i = 0; i < (int)polyline.size()-1; )
  {
    int count = polyline[i];
    if ( count < 2 )
      PRINT_ERROR("Invalid polyline length %d at polyline[%d]\n", count, i );
    i++;
    for ( j = 0; j < count; j++, i++ )
    {
      if( i >= (int)polyline.size() )
        PRINT_ERROR("Ran out of indices in polyline.\n");
      if ( polyline[i] >= num_tot_pts || polyline[i] < 0 )
        PRINT_ERROR("Invalid value %d at dudded[%d]\n", dudded[i], i);
    }
  }
  
    // draw original points
  for ( i = 0; i < (int)coordinates.size(); i += 3 )
    draw_pt( coordinates, i, CUBIT_BLUE );
/*  
    // draw new points
  for ( i = 0; i < (int)new_pts.size(); i += 3 )
    draw_pt( new_pts, i, CUBIT_RED ); 
*/  
    // concatenate point lists
  for ( i = 0; i < (int)new_pts.size(); i++ )
    coordinates.push_back(new_pts[i]);
  
    // make reverse mapping for dudded facets
  std::vector<bool> dead(connections.size()/3);
  for ( i = 0; i < (int)dead.size(); i++ )
    dead[i] = false;
  for ( i = 0; i < (int)dudded.size(); i++ )
    dead[dudded[i]] = true;
    
    // make polyline sequence list (this won't work quite
    // right for self-intersecting polylines -- *shrug*)
  std::vector<int> sequence(coordinates.size() / 3);
  for ( i = 0; i < (int)sequence.size(); i++ )
    sequence[i] = 0;
  j = 0;
  for ( i = 0; i < (int)polyline.size(); i++ )
    sequence[polyline[i]] = ++j;
  
  
    // draw orginal facets
  for ( i = 0; i < (int)dead.size(); i++ )
    if( ! dead[i] )
      draw_tri( coordinates, connections, 3*i, sequence, CUBIT_BLUE );
  
    // draw new facets
  for ( i = 0; i < (int)new_facets.size(); i += 3 )
    draw_tri( coordinates, new_facets, i, sequence, CUBIT_RED );
  
    // draw polyline
  for ( i = 0; i < (int)polyline.size()-1; i++)
  {
    int count = polyline[i];
    if ( count < 2 )
      PRINT_ERROR("Invalid polyline length %d at polyline[%d]\n", count, i );
    int first = ++i;
    for ( j = 0; j < count-1; j++, i++ )
    {
      draw_edge( coordinates, polyline[i], polyline[i+1], CUBIT_WHITE );
    }
    if( polyline[first] == polyline[i] )
      PRINT_INFO("Polyline is closed.\n");
  }

  GfxDebug::flush();
#endif
  cleanup();  
}
