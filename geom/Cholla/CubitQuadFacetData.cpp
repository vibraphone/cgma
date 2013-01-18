//
// File: CubitQuadFacetData.cpp
//
// Owner: sjowen
//

#include "CubitQuadFacetData.hpp"
#include "CubitPoint.hpp"
#include "CubitFacetEdge.hpp"
#include "CubitFacetEdgeData.hpp"
#include "CubitFacet.hpp"
#include "CubitFacetData.hpp"

//===========================================================================
//  Function: CubitQuadFacetData
//  Purpose:  constructor
//  Notes:    defines a quad facet given two existing triangle facets.
//            The points on the facet are also passed in so the quad
//            can be oriented with respect to its triangles
//  Date:     4/11/01
//  Author:   sjowen
//===========================================================================
CubitQuadFacetData::CubitQuadFacetData(
  CubitFacet *tri_facets[2],
  CubitPoint *points[4] )
{
  myTriangleFacets[0] = tri_facets[0];
  myTriangleFacets[1] = tri_facets[1];

  CubitBoolean found = CUBIT_FALSE;
  CubitPoint *pt;
  int ii, jj, kk;

  for (jj=0; jj<2; jj++)
  {
    for (kk=0; kk<3; kk++)
    {
      pt = myTriangleFacets[jj]->point(kk);
      found = CUBIT_FALSE;
      for (ii=0; ii<4 && !found; ii++)
      {
        if (points[ii] == pt)
        {
          triToQuadIndex[jj][kk] = ii;
          found = CUBIT_TRUE;
        }
      }

      // If assertion occurs, the that points passed to the constructor do 
      // not match the points on the triangle facets

      assert(found == CUBIT_TRUE);
    }
  }
}

//===========================================================================
//  Function: CubitQuadFacetData
//  Purpose:  constructor
//  Notes:    defines a quad facet given its ordered point array.
//            Two triangle facets are created
//  Date:     4/11/01
//  Author:   sjowen
//===========================================================================
CubitQuadFacetData::CubitQuadFacetData(
  CubitPoint *point0, CubitPoint *point1, 
  CubitPoint *point2, CubitPoint *point3 )
{
  
  // The orientation of the triangles in the quad is arbitrary.
  // We may want to split based on some other criteria (ie. angles)

  myTriangleFacets[0] = new CubitFacetData( point0, point1, point2 );
  myTriangleFacets[1] = new CubitFacetData( point2, point3, point0 );

  triToQuadIndex[0][0] = 0;
  triToQuadIndex[0][1] = 1;
  triToQuadIndex[0][2] = 2;
  triToQuadIndex[1][0] = 2;
  triToQuadIndex[1][1] = 3;
  triToQuadIndex[1][2] = 0;
}

//===========================================================================
//  Function: CubitQuadFacetData
//  Purpose:  constructor
//  Notes:    defines a quad facet given its ordered edge array.
//            Two triangle facets are created
//  Date:     4/11/01
//  Author:   sjowen
//===========================================================================
CubitQuadFacetData::CubitQuadFacetData(
  CubitFacetEdge *e0, CubitFacetEdge *e1,
  CubitFacetEdge *e2, CubitFacetEdge *e3 )
{
  // create the diagonal edge

  CubitPoint *p0 = e3->shared_point(e0);
  CubitPoint *p2 = e1->shared_point(e2);
  assert(p0 != NULL && p2 != NULL);

  CubitFacetEdge *e4 = (CubitFacetEdge *) new CubitFacetEdgeData( p0, p2 );

  // create the triangles

  myTriangleFacets[0] = new CubitFacetData( e1, e4, e0 );
  myTriangleFacets[1] = new CubitFacetData( e3, e4, e2 );

  triToQuadIndex[0][0] = 0;
  triToQuadIndex[0][1] = 1;
  triToQuadIndex[0][2] = 2;
  triToQuadIndex[1][0] = 2;
  triToQuadIndex[1][1] = 3;
  triToQuadIndex[1][2] = 0;
}

//===========================================================================
//  Function: CubitQuadFacetData
//  Purpose:  constructor
//  Notes:    defines a quad facet given its ordered point array.
//            Two triangle facets are created
//  Date:     4/11/01
//  Author:   sjowen
//===========================================================================
CubitQuadFacetData::CubitQuadFacetData(
  CubitPoint *points[4] )
{
  
  // The orientation of the triangles in the quad is arbitrary.
  // We may want to split based on some other criteria (ie. angles)

  myTriangleFacets[0] = new CubitFacetData( points[0], points[1], points[2] );
  myTriangleFacets[1] = new CubitFacetData( points[2], points[3], points[0] );

  triToQuadIndex[0][0] = 0;
  triToQuadIndex[0][1] = 1;
  triToQuadIndex[0][2] = 2;
  triToQuadIndex[1][0] = 2;
  triToQuadIndex[1][1] = 3;
  triToQuadIndex[1][2] = 0;
}

//===========================================================================
//  Function: ~CubitQuadFacetData
//  Purpose:  destructor
//  Date:     4/11/01
//  Author:   sjowen
//===========================================================================
CubitQuadFacetData::~CubitQuadFacetData()
{
  if (myTriangleFacets[0] == NULL && myTriangleFacets[1] == NULL)
    return;

  assert(myTriangleFacets[0] && myTriangleFacets[1]);

  // determine the dialogonal edge
  CubitPoint *p0 = this->point(0);
  CubitPoint *p2 = this->point(2);
  CubitFacetEdge *dedge = NULL;
  if (p0 != NULL && p2 != NULL)
  {
    dedge = p0->shared_edge(p2);
  }

  // delete both triangles
  delete myTriangleFacets[0];
  delete myTriangleFacets[1];

  // delete the dialgonal edge
  if (dedge != NULL)
  { 
    delete dedge;
  }
}

//===========================================================================
//  Function: ~CubitQuadFacetData
//  Purpose:  destructor
//  Date:     4/11/01
//  Author:   sjowen
//===========================================================================
void CubitQuadFacetData::remove_tri_facets(  )
{
  myTriangleFacets[0] = NULL;
  myTriangleFacets[1] = NULL;
}

//===========================================================================
//  Function: points
//  Purpose:  get the points from the facet
//  Date:     4/11/01
//  Author:   sjowen
//===========================================================================
void CubitQuadFacetData::points(
  CubitPoint *thepoints[4] )
{

  thepoints[triToQuadIndex[0][0]] = myTriangleFacets[0]->point(0);
  thepoints[triToQuadIndex[0][1]] = myTriangleFacets[0]->point(1);
  thepoints[triToQuadIndex[0][2]] = myTriangleFacets[0]->point(2);
  thepoints[triToQuadIndex[1][0]] = myTriangleFacets[1]->point(0);
  thepoints[triToQuadIndex[1][1]] = myTriangleFacets[1]->point(1);
  thepoints[triToQuadIndex[1][2]] = myTriangleFacets[1]->point(2);
}

//===========================================================================
//  Function: point
//  Purpose:  get the point from the facet
//  Date:     11/28/2002
//  Author:   sjowen
//===========================================================================
CubitPoint *CubitQuadFacetData::point( int index )
{
  int ii, jj;
  for (ii=0; ii<2; ii++)
  {
    for (jj=0; jj<3; jj++)
    {
      if (triToQuadIndex[ii][jj] == index)
      {
        return myTriangleFacets[ii]->point(jj);
      }
    }
  }
  assert(0);  // index is probably out of range
  return NULL;
}

//===========================================================================
//  Function: edge
//  Purpose:  get the edge from the facet
//  Date:     11/28/2002
//  Author:   sjowen
//===========================================================================
CubitFacetEdge *CubitQuadFacetData::edge( int index )
{
  int ii, jj;
  CubitPoint *p0 = NULL;
  CubitPoint *p1 = NULL;
  int index1 = (index + 1) % 4;
  for (ii=0; ii<2; ii++)
  {
    for (jj=0; jj<3; jj++)
    {
      if (triToQuadIndex[ii][jj] == index)
      {
        p0 = myTriangleFacets[ii]->point(jj);
      }
      else if ( triToQuadIndex[ii][jj] == index1 )
      {
        p1 = myTriangleFacets[ii]->point(jj);
      }
    }
  }
  assert(p0 != NULL && p1 != NULL);  // index is probably out of range
  
  return p0->shared_edge( p1 );
}


//===========================================================================
//  Function: get_tri_facet
//  Purpose:  return the underlying triangle of which the point_ptr is a vertex
//  Date:     11/28/2002
//  Author:   sjowen
//===========================================================================
CubitFacet *CubitQuadFacetData::get_tri_facet_at_point( CubitPoint *point_ptr )
{
  int ii, jj;
  for (ii=0; ii<2; ii++)
  {
    for (jj=0; jj<3; jj++)
    { 
      if (myTriangleFacets[ii]->point(jj) == point_ptr)
      {
        return myTriangleFacets[ii]; 
      }
    }
  }
  assert(0);  // point isn't on facet
  return NULL;
}

