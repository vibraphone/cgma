/** \class ChollaCurve.hpp
//- Class: ChollaCurve
//- Owner: Steven J. Owen
//- Description: Maintains a list of mesh.  This is used to store
//-        the exterior skin while creating a geometry from a mesh.
//- Created: 12/3/00
//- Checked By:
//- Version:
*/
#ifndef ChollaCurve_HPP
#define ChollaCurve_HPP

#define MYLENGTH_UNINITIALIZED -1.0

#include "DLIList.hpp"
#include "ChollaEntity.hpp"

class CubitVector;
class CubitFacetEdge;
class CubitPoint;
class ChollaPoint;
class ChollaSurface;
class ChollaVolume;
class FacetEntity;
class CurveFacetEvalTool;

class ChollaCurve : public ChollaEntity
{
private:

  DLIList<ChollaSurface*> surfaceList;
  DLIList<FacetEntity*> curveEdgeList;
  void *myCurve;
  CurveFacetEvalTool *myEvalTool;
  CubitPoint *startPoint;
  CubitPoint *endPoint;
  DLIList<ChollaPoint*> pointList;
  int flag;
  int blockID;
  int id;
  double myLength;
  ChollaCurve *myMergePartner;

  CubitStatus determine_ends();
    //- determine the end nodes for the curve (once the curvEdgeList is full)

  CubitFacetEdge *next_edge( CubitPoint *node_ptr, CubitFacetEdge *edge_ptr );
    //- return the next edge on the curve (NULL if at end)
  
  
public:


  ChollaCurve( int block_id );
    //- default constructor

  ~ChollaCurve();
    // destructor

  int get_block_id()
  { return blockID; }

  void remove_td_associativity( ChollaSurface *fsm_ptr );
    // delete the asociativity of this curve with all edge's tool datas

  void add_facet(FacetEntity *exterior_edge)
    {myLength = MYLENGTH_UNINITIALIZED;
    curveEdgeList.append(exterior_edge);}
    //- add an edge to the curve

  int add_facet_unique(FacetEntity *exterior_edge)
  {myLength = MYLENGTH_UNINITIALIZED;
  return curveEdgeList.append_unique(exterior_edge);}
    //- add an edge to this curve - check to see if it already there before adding

  void remove_facet( FacetEntity *facet_edge )
  {myLength = MYLENGTH_UNINITIALIZED;
  curveEdgeList.remove( facet_edge ); }
    //- remove a facet_edge from underlying backing

  CubitStatus replace_facet( FacetEntity *remove_edge, FacetEntity *replace_edge );
    //- replace a facet_edge from underlying backing

  
  CubitStatus insert_facet( FacetEntity *old_edge, FacetEntity *new_edge );
    //- inserts new_edge after old_edge
    
  CubitStatus is_contain( FacetEntity *edge );
  //- check to see if an edge is in the list

  DLIList<FacetEntity*> &get_facet_list()
    {return curveEdgeList;}
    //- get the list of edges that make up this curve
  DLIList<FacetEntity*> *get_facet_list_ptr()
    {return &curveEdgeList;}
  
  //- return the length of the curve
  double length();
  
  int num_edges() {return curveEdgeList.size();}
    //- return the number of edges in the curve

  void add_surface( ChollaSurface *fsm_ptr )
    {surfaceList.append_unique( fsm_ptr );}
    //- associate a surface with this curve

  inline void remove_surface( ChollaSurface *fsm_ptr)
    {surfaceList.remove(fsm_ptr);}
    //- remove a suface from the curve

  DLIList<ChollaSurface*> &get_surfaces()
    {return surfaceList;}
  void get_surfaces( DLIList<ChollaSurface *> &surf_list )
    { surf_list += surfaceList; }
    //- get the list of surfaces attached to this curve

  DLIList<ChollaSurface*> *get_surface_list_ptr()
    {return &surfaceList;} 

  void add_point( ChollaPoint *fpm_ptr )
    {pointList.append_unique( fpm_ptr );}
    //- associate a point with this curve
   
  inline void remove_point( ChollaPoint *fpm_ptr)
    {pointList.remove(fpm_ptr);}
    //- remove a point from the curve

  DLIList<ChollaPoint*> &get_points()
    {return pointList;}
    //- get the list of points attached to this curve
  
  void get_facet_points( DLIList<CubitPoint *> &point_list, CubitBoolean inclusive);
    //- return the facet points on this chollacurve
  
  CubitBoolean is_in_volume( ChollaVolume *chvol_ptr );
    //- return whether this curve is contained within the specified volume
  
  CubitBoolean is_in_surface( ChollaSurface *chsurf_ptr );
    //- return whether this surface is contained within the specified surface

  void assign_geometric_curve(void *curv)
    {myCurve = curv;}
    //- set the geometric curve associated with the ChollaCurve
    //- myCurve is void * so we can deal with any type of geometry engine

  void* get_geometric_curve()
    {return myCurve;}
    //- return the geometric curve associated with the ChollaCurve

  void assign_eval_tool(CurveFacetEvalTool *curv_eval_tool)
    {myEvalTool = curv_eval_tool;}
    //- set the curve facet eval tool associated with the ChollaCurve

  CurveFacetEvalTool* get_eval_tool()
    {return myEvalTool;}
    //- return the curve facet eval tool associated with the ChollaCurve

  CubitStatus get_ends( CubitVector &start, CubitVector &end );
    //- return the end node locations
  CubitStatus get_ends( CubitPoint *&start_ptr, CubitPoint *&end_ptr );
    //- return the end nodes

  void set_start( CubitPoint *start_pt )
    { startPoint = start_pt; }
    // set the start point on the curve

  void set_end( CubitPoint *end_pt )
    { endPoint = end_pt; }
    //- set the end point on the curve

  CubitStatus split_curve( DLIList<ChollaCurve*> &facet_curve_list );
    //- split this curve into multiple ChollaCurve where there are
    //- discontinuous strings of edges.  Define start end end points
    //- for each curve while we are at it.
  
  CubitStatus build_curve_from_edges(CubitPoint *start_point,
                                     int periodic, int max_edges,
                                     CubitFacetEdge *start_edge_ptr,
                                     ChollaCurve *parent_curve);

  CubitStatus order_edges();
    //- put the current edges in the curve's list in order from start to end

  CubitStatus feature_angle( double min_dot );
    //- compute angles at nodes on the curve to see if we need to split
    //- the curve.  Mark the node tooldata hitflag if the node will
    //- break the curve (this is refernced in next_edge) 

  //  disassociate from cholla points
  CubitStatus disassociate_from_points( void );
  
  // disassociate from cholla surface
  CubitStatus disassociate_from_surfaces( void);

  // build a new curve_facet_eval_tool and assign it to ChollaCurve
  CubitStatus build_curve_facet_eval_tool( void );

  int get_flag( ) { return flag; }
  void set_flag( int flg ) { flag = flg; }
  int get_id() {return id;}
  void debug_draw();
  void print();

  // find adjacent edges at a point that lie on the curve
  bool adj_facet_edges( CubitPoint *cubit_pnt, CubitFacetEdge *&adj_edge1, CubitFacetEdge *&adj_edge2 );
  
  // return whether the curve contains the point
  CubitBoolean has_point( ChollaPoint *pt );
  
  // verify that all points on this curve have this curve as an adjacency
  CubitStatus verify_points();
  
  // clear the edge list
  void clean_out_edges(){ curveEdgeList.clean_out(); }

  int num_volumes();

  ChollaCurve *merge_partner(){ return myMergePartner; }

  void set_merge_partner( ChollaCurve *merge_partner )
    { myMergePartner = merge_partner;}  

};

#endif

