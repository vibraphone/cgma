#ifndef CUBITFACETDATA_HPP 
#define CUBITFACETDATA_HPP 
 
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
 
class CubitFacetData : public CubitFacet 
{ 
private: 
  CubitPoint *pointArray[3]; 
  CubitFacetEdge *edgeArray[3]; 
  int edgeUse[3]; 
  static MemoryManager memoryManager; 
    //- memory management object 
  int entityId; 
    //- for debug tracking...
  void allocate_edge(CubitPoint *p1, CubitPoint *p2, int edge_index);
  void define_point(CubitFacetEdge *e1, CubitFacetEdge *e2, 
                    int point_index);
  void define_bounding_box();
 
public: 
  CubitFacetData( CubitPoint *p1, CubitPoint *p2, 
                  CubitPoint *p3);
  CubitFacetData( CubitPoint *p1, CubitPoint *p2,
                  CubitPoint *p3, int *tool_data);
  CubitFacetData( CubitFacetEdge *e1, CubitFacetEdge *e2, 
                  CubitFacetEdge *e3);
  ~CubitFacetData();
  
  void destruct_facet_internals();
  
  virtual int id(){return entityId;} 
  virtual void set_id( int ii ) { entityId = ii; }
       
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
 
   
  void set_point( CubitPoint *the_point, int index );  
    //- sets the point into the facet. 
 
  CubitFacetEdge *edge( int index ); 
  void edge( CubitFacetEdge *the_edge, int index ); 
    //- get and set edge pointers 
 
  void edge_use( int direction, int index ); 
  int edge_use( int index ); 
    //- get and set the edge uses 
  int sense(int index); 
 
  CubitPoint* split_edge( int edge_index, const CubitVector& position );
  CubitPoint* split_edge( CubitPoint* edge_pt1, CubitPoint* edge_pt2, 
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
  
  CubitStatus flip_edge( int edge_index );
  CubitStatus flip_edge( CubitFacetEdge *edge );
    //- Given an edge shared by exactly two facets, 
    //- flip the edge such that it connects the 
    //- "other" two facet points instead.
    //-
    //-  *----*        *----*
    //-  |\   |        |   /|
    //-  | \  |  ==>   |  / |  
    //-  |  \ |        | /  |
    //-  |   \|        |/   |
    //-  *----*        *----*
    //-
    //- NOTE:  If you are trying to reverse the sense of
    //-        an edge, this is NOT the function you want!
   
  void flip(); 
    //- reorient the facet 

}; 


inline CubitPoint* CubitFacetData::point( int index ) 
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
 
inline void CubitFacetData::set_point( CubitPoint* the_point, int index ) 
{ 
  assert( (index >= 0) && (index < 3) ); 
  if (!is_backwards()) 
    pointArray[index] = the_point; 
  else 
  { 
    switch(index) 
    { 
    case 0: pointArray[0] = the_point; break; 
    case 1: pointArray[2] = the_point; break; 
    case 2: pointArray[1] = the_point; break; 
    } 
  } 
} 
 
inline CubitFacetEdge* CubitFacetData::edge( int index ) 
{ 
  assert( (index >= 0) && (index < 3) ); 
  if (!is_backwards()) { 
    return edgeArray[index]; 
  } 
  else 
  { 
    switch(index) 
    { 
    case 0:return edgeArray[0]; 
    case 1:return edgeArray[2]; 
    case 2:return edgeArray[1]; 
    } 
  } 
  return NULL; 
} 
 
 
inline void CubitFacetData::edge( CubitFacetEdge *the_edge, int index ) 
{ 
  assert( (index >= 0) && (index < 3) ); 
  if (!is_backwards()) 
    edgeArray[index] = the_edge; 
  else 
  { 
    switch(index) 
    { 
    case 0: edgeArray[0] = the_edge; break; 
    case 1: edgeArray[2] = the_edge; break; 
    case 2: edgeArray[1] = the_edge; break; 
    } 
  } 
} 
 
inline int CubitFacetData::edge_use( int index )  
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

//- This function is defined so that Meshing and Geometry Entities 
//- return the same values as needed in the FacetorTool
inline int CubitFacetData::sense( int index ) { 
  if(edge_use(index) == 1)
    return CUBIT_FORWARD;
  else if(edge_use(index) == -1)
    return CUBIT_REVERSED;
  else 
    return CUBIT_UNKNOWN;
} 
 
inline void CubitFacetData::edge_use( int direction, int index )  
{  
  assert( (index >= 0) && (index < 3) ); 
  if (!is_backwards()) 
    edgeUse[index] = direction;  
  else 
  { 
    switch(index) 
    { 
    case 0: edgeUse[0] = -direction; break; 
    case 1: edgeUse[2] = -direction; break; 
    case 2: edgeUse[1] = -direction; break; 
    } 
  } 
} 
 
#endif 


