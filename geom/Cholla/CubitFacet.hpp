//----------------------------------------------------------------------------- 
// 
//  File: CubitFacet.hpp 
// 
//  Purpose:  Facet Class used for mesh based geometry and other tools. 
// 
//  Notes:    A facet is a triangle used for defining surfaces. The 
//            default construction is a linear representation, although 
//            b-spline control points can be defined to represent a 
//            higher order triangular b-spline patch. There is also  
//            support for quad facets by using the CubitQuadFacet 
// 
//            Note that the main private data components (ie. points, edges) 
//            are not contained in this class. 
//            This data should be defined within the child classes inherited 
//            from this class.  The current Cubit data class that inherits 
//            from CubitFacet is CubitFacetData.  This is done so that  
//            other applications using CubitFacets can use their own 
//            facet data, but take advantage of the CGM/Cubit functionality. 
//            Please do not add private data to this class unless it is only 
//            applicable to Cubit algorithms; instead add the 
//            data to the children and access through virtual functions. 
//  
//            Do not create a CubitFacet directly.  For example, don't do:  
//              CubitFacet *cf = new CubitFacet(...); 
//            You should instead create the appropriate child class, and  
//            cast it to a CubitFacet for use.  For example: 
//              CubitFacet *cf = (CubitFacet *) new CubitFacetData(...); 
//                         
//----------------------------------------------------------------------------- 
 
#ifndef CUBITFACET_HPP 
#define CUBITFACET_HPP 
 
#include "CubitDefines.h" 
#include "MemoryManager.hpp" 
#include "DLIList.hpp" 
#include "CubitBox.hpp" 
#include "FacetEntity.hpp" 
#include "CubitVector.hpp" 
#include "CubitPlane.hpp" 
 
class CubitFacetEdge; 
class CubitPoint; 
 
class CubitFacet : public FacetEntity 
{ 
public: 
  static const int point_edge_conn[30][2]; 
  static const int point_facet_conn[16][3]; 
  static const double my_points[15][3]; 
 
protected: 
 
  CubitPlane *cachedPlane; 
    //- cached plane for this facet 
 
  CubitBox bBox; 
    //- The bounding box for the facet 
 
  double facetWeight; 
    //- for interpolation 
 
  CubitVector *patchCtrlPts; 
    //- control points for triangular bezier batch 
 
  int markedFlag; 
    //- generic marker for list sorting.  Assume 0. Users must clean after use. 
 
  int isFlat; 
    //- facet is flat 
 
  int isBackwards; 
    //- facet is oriented backwards wrt. orientation of surface 
 
  int toolID; 
    //- id of facet eval tool 
 
public: 
  CubitFacet(); 
 
//------------------------------------------------------------------------ 
// The following functions are pure virtual and must be  
// implemented in the child class (ie. CubitFacetData or equivalent) 
 
  virtual ~CubitFacet(); 
 
  CubitVector normal();  
    //- Returns the normal to this facet. 

  CubitVector update_normal( void );
  //- update the plane and return the new normal


  const CubitBox& bounding_box() {return bBox;}; 
  void bounding_box( CubitBox &box ) 
  { bBox = box; } 
  void reset_bounding_box(); 
    //- returns/sets the bounding box of the facet 
 
  void weight(double facweight) { facetWeight = facweight; } 
  double weight() { return facetWeight; } 
    //- get/set weight 
 
  void get_control_points( CubitVector points[6] ); 
  void set_control_points( CubitVector points[6] ); 
  void set_control_points( const double *pt_array );
  CubitVector *control_points() {return patchCtrlPts;};
    //- get and set the bezier patch internal control points 
 
  int eval_order() 
    { return (patchCtrlPts == NULL) ? 0 : 4; } 
    //- return the evaluation order of the facet.  Currently only 
    //- order 0 (linear) and 4 (b-spline patch) are implemented 
     
  virtual CubitPoint* point( int index ) = 0; 
    //Get the point at the specified index. 
    //asserts that the index is in range, for a 
    //triangle, 0 <= index <= 2. 
 
  void marked(int mark) {markedFlag = mark;}; 
  int marked() {return markedFlag;}; 
    //- generic marking functions. 
 
  int is_flat(); 
  void is_flat(int flat) {isFlat = flat;}; 
    //- get/set isFlat 
     
  const CubitPlane& plane(); 
  void plane(CubitPlane &this_plane); 
  void update_plane();  
    //- Cache and return the plane containing this triangle.  
 
  int is_backwards() {return isBackwards;}; 
  void is_backwards( int flipped ) { isBackwards = flipped; } 
    //- return whether the facet was reoriented 
 
  virtual int id() = 0; 
  virtual void set_id( int ii ) = 0; 
  int tool_id(){return toolID;} 
  void set_tool_id(int tool_id){toolID = tool_id;} 
       
  virtual CubitFacetEdge *edge( int index ) = 0; 
  virtual void edge( CubitFacetEdge *the_edge, int index ) = 0; 
    //- get and set edge pointers 
 
  virtual void edge_use( int direction, int index ) = 0; 
  virtual int edge_use( int index ) = 0; 
    //- get and set the edge uses.  Direction is 1 or -1 based 
    //- upon the orientation of the edge with respect to a CCW orientation 
    //- of the facet 
 
  virtual void flip() = 0; 
    //- reorient the facet 
 
//--------------------------------------------------------------------- 
//  The following virtual functions are optional.  Implement them 
//  based upon the application's need.  If they are called without 
//  being implemented, an assertion will occur. 
 
  virtual int sense(int /*index*/){ assert(0); return -1;}
  //- Gives the direction of the edge given the edge index

  virtual CubitPoint* split_edge( CubitPoint* edge_pt1,  
                                  CubitPoint* edge_pt2, 
                                  const CubitVector& position ); 
    //R CubitPoint 
    //R- The new CubitPoint created. 
    //I edge_pt1, edge_pt2 
    //I- The end points of an edge of this triangle. 
    //I position 
    //I- The position at which to split the edge. 
    //- Split an edge on this triangle and the other triangle 
    //- sharing the edge, if it exists. 
    //- 
    //- Note:  No check is done on the location of the split 
    //-        position. 
     
  virtual CubitPoint* insert_point( const CubitVector& position, 
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
 
  virtual CubitStatus closest_point( const CubitVector &point,  
                                     CubitVector &closest_point); 
    //- Sets the closest point on the plane defined by  
    //- this facet to the point in space. 
    //- If the normal length to the facet is 0, it will return CUBIT_FAILURE. 
 
//---------------------------------------------------------------------- 
// The following functions access the data in a child class 
 
  CubitVector center(); 
    //- Returns the center of the facet. 
 
  void debug_draw(int color=-1, int flush_it = 1, int draw_uv=0); 
    //- draws the current facet. 
 
  CubitStatus closest_point_trimmed( const CubitVector &point,  
                                     CubitVector &closest_point, 
                                     CubitPoint *&next_edge_p1,  
                                     CubitPoint *&next_edge_p2 ); 
    //- returns the closest point on this facet to the point in space. 
    //- If the point is not in the facet, then the edge to which it was closest 
    //- will be set and returned.  Other wise they will be "NULLED". 
 
  void get_edge_1( CubitPoint *&p1, CubitPoint *&p2 ) 
    { p1 = point(0); p2 = point(1);} 
  void get_edge_2( CubitPoint *&p1, CubitPoint *&p2 ) 
    { p1 = point(1); p2 = point(2);} 
  void get_edge_3( CubitPoint *&p1, CubitPoint *&p2 ) 
    { p1 = point(2); p2 = point(0);} 
    //- Get the points that define the three edges of the facet. 
  void get_edge_pts( int index, CubitPoint *&p1, CubitPoint *&p2 ); 
    //- get the points that define an edge 
    //- The index of the edge corresponds to the point index on 
    //- the facet opposite the edge 
 
  CubitFacetEdge *edge_from_pts( CubitPoint *p1, CubitPoint *p2, int &index ); 
    //- return the edge pointer based on its end points 
  int edge_index( CubitPoint *p1, CubitPoint *p2, int &sense ); 
  int edge_index( CubitFacetEdge *edge ); 
    //- return the index of the edge on the facet 
  int point_index( CubitPoint *pt ); 
    //- return the index of the point on the facet 
     
  void opposite_edge( CubitPoint* point,  
                      CubitPoint*& p1, CubitPoint*& p2 ); 
    //- Get the edge on the triangle opposite of the passed point. 
    //- p1 and p2 will be passed back as NULL if point is not a 
    //- point on this triangle. 

  
  CubitPoint *opposite_point( CubitFacetEdge *edge );
  //- returns the opposite point of input edge.  Input edge should be incident on facet. 
  //- If input edge is not incident on facet then any of the one point can be returned
 
  
  void points( CubitPoint *&p0, CubitPoint *&p1, CubitPoint *&p2 ) 
    { p0 = point(0); p1 = point(1); p2 = point(2); } 
  void tri_nodes( CubitPoint *&p0, CubitPoint *&p1, CubitPoint *&p2 )  
    { points( p0, p1, p2); } 
  void points(DLIList<CubitPoint*> &point_list ) 
    { for ( int i = 0; i < 3; i++ ) 
      point_list.append(point(i)); 
    } 
  void facets(DLIList<CubitFacet*> &facet_list ) 
    { facet_list.append( this ); } 
  void edges(DLIList<CubitFacetEdge*> &edge_list ) 
    { for ( int i = 0; i < 3; i++ ) 
      edge_list.append(edge(i)); 
    } 
 
  CubitBoolean contains(CubitPoint *p1); 
    //- Returns TRUE/FALSE if the point defines this facet. 

  CubitFacetEdge* shared_edge( CubitFacet *cubit_facet );
 
  CubitFacet* shared_facet( CubitPoint *p1, CubitPoint *p2 ); 
    //- Returns the "other" facet attached to this edge. (if  
    //- there is one., other wise returns NULL). 
  void shared_facets( CubitPoint *p1, CubitPoint *p2, 
                      DLIList <CubitFacet *> &adj_facet_list); 
    //- Returns all adjacent facets to this edge (not including 
    //- this facet) 
  
  CubitFacet *adjacent( int &index, int *tool_data )  
    //- Returns all adjacent facets to this edge index(not including 
    //- this facet) 
  {
		//convert MESHING TRI edge index to GEOMETRY FACET edge index
		index = (index+2)%3;
		CubitPoint *p1, *p2; 
   
		get_edge_pts( index, p1, p2 ); 
		return shared_facet_on_surf( p1, p2, *tool_data ); 
  } 

  CubitPoint* next_node(CubitPoint *current_point)
  {
	int index = point_index(current_point);
	return point((index+1)%3);
  }
    //- return the adjacent triangle at the edge index 
  CubitFacet* shared_facet_on_surf( CubitPoint *p1, CubitPoint *p2, int tool_id ); 
    //- same as above except also matches tool_id with toolID 
 
  double min_diagonal(); 
    //- Returns the length of the minimum diagonal of the triangle. 
   
  double angle( CubitPoint *pt ); 
    //- return the angle (radians) at one of the points on the facet 
 
  void update_bezier_bounding_box( ); 
    //- Update the facet bounding box based on its control polygon 
 
  CubitFacetEdge *prev_edge( CubitFacetEdge *edge ); 
  CubitFacetEdge *next_edge( CubitFacetEdge *edge ); 
    //- return the prvious or next edge on the facet 
  
  int other_index( CubitPoint* pt1, CubitPoint* pt2 ); 
 
  CubitStatus evaluate_position( const CubitVector &start_position, 
                                 CubitVector *eval_point, 
                                 CubitVector *eval_normal = NULL );   
 
  CubitStatus evaluate( CubitVector &areacoord,  
                        CubitVector *eval_point, 
                        CubitVector *eval_normal = NULL );  
 
  void get_parents(DLIList<FacetEntity *> &){}; 
    // dummy for this class 
 
  CubitFacetEdge *next_edge_at_point( CubitFacetEdge *edge_ptr, 
                                      CubitPoint *point_ptr ); 
    // return the next edge on the triangle at the point 
 
  CubitStatus get_edge_control_points( CubitVector P[3][5] ); 
    // get the control points on the facet edges 
 
  double area(); 
    // return the area of the facet 

  double aspect_ratio();


  CubitStatus init_patch(  );
    // computes the interior control points for the facet and
    // stores them with the facet.  Assumes edge control points
    // have already been computed

  void add_edge(CubitFacetEdge *edge);
   //- add an existing edge to a facet

  void unlink_from_children( void );
 
}; 
 
inline void CubitFacet::plane(CubitPlane &this_plane)  
{ 
  this_plane = plane(); 
} 
 
#endif 


