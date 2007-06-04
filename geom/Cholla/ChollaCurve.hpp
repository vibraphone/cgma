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

#include "DLIList.hpp"
#include "ChollaEntity.hpp"

class CubitVector;
class CubitFacetEdge;
class CubitPoint;
class ChollaPoint;
class ChollaSurface;
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
    {curveEdgeList.append(exterior_edge);}
    //- add an edge to the curve

  int add_facet_unique(FacetEntity *exterior_edge)
    {return curveEdgeList.append_unique(exterior_edge);}
    //- add an edge to this curve - check to see if it already there before adding

  DLIList<FacetEntity*> &get_facet_list()
    {return curveEdgeList;}
    //- get the list of edges that make up this curve
  DLIList<FacetEntity*> *get_facet_list_ptr()
    {return &curveEdgeList;}

  void add_surface( ChollaSurface *fsm_ptr )
    {surfaceList.append_unique( fsm_ptr );}
    //- associate a surface with this curve

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

  DLIList<ChollaPoint*> &get_points()
    {return pointList;}
    //- get the list of points attached to this curve

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

  CubitStatus feature_angle( double min_dot );
    //- compute angles at nodes on the curve to see if we need to split
    //- the curve.  Mark the node tooldata hitflag if the node will
    //- break the curve (this is refernced in next_edge) 

  int get_flag( ) { return flag; }
  void set_flag( int flg ) { flag = flg; }
  int get_id() {return id;}
  void debug_draw();
  void print();
};

#endif

