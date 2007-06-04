#ifndef FACETERFACETDATA_HPP
#define FACETERFACETDATA_HPP

// Include for CubitBoolean
#include "CubitDefines.h"
#include "MemoryManager.hpp"
#include "CubitPoint.hpp"
#include "DLIList.hpp"
#include "CubitBox.hpp"
#include "CubitFacet.hpp"
class CubitVector;
class CubitBox;
class CubitPlane;
class CubitFacetEdge;
class RefEntity;

class FaceterFacetData : public CubitFacet
{
private:
  CubitPoint *pointArray[3];
  CubitFacetEdge *edgeArray[3];
  int edgeUse[3];
  static MemoryManager memoryManager;
    //- memory management object
  int entityId;
  RefEntity *myOwner;
  
    //- for debug tracking...

public:
  FaceterFacetData(CubitPoint *p1, CubitPoint *p2, CubitPoint *p3);
  ~FaceterFacetData();
  virtual int id(){return entityId;}
  virtual void set_id(int ii){entityId = ii;}

  void owner(RefEntity *owner)
    {myOwner = owner;}
  RefEntity* owner()
    {return myOwner;}
  SetDynamicMemoryAllocation(memoryManager)
    //- class specific new and delete operators
  CubitStatus closest_point( const CubitVector &point, 
                                     CubitVector &closest_point);
    //- Sets the closest point on the plane defined by 
    //- this facet to the point in space.
    //- If the normal length to the facet is 0, it will return CUBIT_FAILURE.

  CubitPoint* point( int index );
    //Get the point at the specified index.
    //asserts that the index is in range, for a
    //triangle, 0 <= index <= 2.

  
  void set_point( CubitPoint *the_point, int index ) { pointArray[index] = the_point; } 
    //- sets the point into the facet.

  CubitFacetEdge *edge( int index );
  void edge( CubitFacetEdge *the_edge, int index )
    { edgeArray[index] = the_edge; };
    //- get and set edge pointers

  void edge_use( int direction, int index ) { edgeUse[index] = direction; };
  int edge_use( int index );
    //- get and set the edge uses

  CubitPoint* split_edge( CubitPoint* edge_pt1, 
                          CubitPoint* edge_pt2,
                          const CubitVector& position );
    //R CubitPoint
    //R- The new CubitPoint created.
    //I edge_pt1, edge_pt2
    //I- The end points of an edge of this triangle.
    //I position
    //I- The position at which to split the edge.
    //O new_facet1
    //O- The new facet resulting from splitting the edge on
    //O- this triangle.
    //O new_facet2
    //O- The new facet resulting from splitting the same edge
    //O- on the ajacent triangle.  If there is no other 
    //O- triangle sharing the edge, NULL will be passed back.
    //- Split an edge on this triangle and the other triangle
    //- sharing the edge, if it exists.
    //-
    //- Note:  No check is done on the location of the split
    //-        position.
    
  CubitPoint* insert_point( const CubitVector& position,
                            CubitFacet*& new_tri1,
                            CubitFacet*& new_tri2 );
    //R CubitPoint
    //R- The new CubitPoint created.
    //I position
    //I- The position at which to insert a point in the triangle.
    //O new_facet1, new_facet2
    //O- The two new facets created.
    //- Insert a point in the interior of this triangle.
    //-
    //- Note:  No check is done on the location of the split
    //-        position.
  
  void flip();
    //- reorient the facet
};

inline CubitPoint* FaceterFacetData::point( int index )
{
  assert( (index >= 0) && (index < 3) );
  if (!is_backwards())
    return pointArray[index];
  else
  {
    switch(index)
    {
    case 0: return pointArray[0];
    case 1: return pointArray[2];
    case 2: return pointArray[1];
    }
  }
  return NULL;
}

inline CubitFacetEdge* FaceterFacetData::edge( int index )
{
  assert( (index >= 0) && (index < 3) );
  if (!is_backwards())
    return edgeArray[index];
  else
  {
    switch(index)
    {
    case 0: return edgeArray[0];
    case 1: return edgeArray[2];
    case 2: return edgeArray[1];
    }
  }
  return NULL;
}

inline int FaceterFacetData::edge_use( int index ) 
{ 
  assert( (index >= 0) && (index < 3) );
  if (!is_backwards())
    return edgeUse[index]; 
  else
  {
    switch(index)
    {
    case 0: return -edgeUse[0];
    case 1: return -edgeUse[2];
    case 2: return -edgeUse[1];
    }
  }
  return 0;
}


#endif
