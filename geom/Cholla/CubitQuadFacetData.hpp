#ifndef CUBITQUADFACETDATA_HPP
#define CUBITQUADFACETDATA_HPP

#include "CubitQuadFacet.hpp"
class CubitFacet;
class CubitFacetEdge;
class CubitPoint;

class CubitQuadFacetData : public CubitQuadFacet
{
private:
  CubitFacet *myTriangleFacets[2];
  int triToQuadIndex[2][3];

public:
  CubitQuadFacetData( CubitFacet *tri_facets[2],
                  CubitPoint *points[4] );
  CubitQuadFacetData( CubitPoint *points[4] );
  CubitQuadFacetData( CubitPoint *point0, CubitPoint *point1, 
                      CubitPoint *point2, CubitPoint *point3 );
  CubitQuadFacetData( CubitFacetEdge *e0, CubitFacetEdge *e1,
                      CubitFacetEdge *e2, CubitFacetEdge *e3 );
  ~CubitQuadFacetData();

  int tri_to_quad_index( int tri_index, int pt_index ) 
    { return triToQuadIndex[ tri_index ][ pt_index ]; }
    // return the point index on the quad facet gicen one of its triangles
    // and the point index on the triangle (indecies start at 0)

  CubitFacet *get_tri_facet( int index )
    { return myTriangleFacets[ index ]; }
    // return one of the underlying triangles (index should be 0 or 1) 
  CubitFacet *get_tri_facet_at_point( CubitPoint *point_ptr );
    // return the underlying triangle of which the point_ptr is a vertex

  void remove_tri_facets(  );
    // used prior to delete if tri facets have already been deleted by
    // another mechanism.

  void points( CubitPoint *points[4] );
    // return all four point
  void points( DLIList<CubitPoint*> &point_list )
    {CubitQuadFacet::points(point_list);}
    //implement the other points function... just call parent.
  CubitPoint *point( int index );
    // return the specified point
  CubitFacetEdge *edge( int index );
    // return the specified edge

};

#endif


