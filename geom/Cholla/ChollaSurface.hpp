//- Class: ChollaSurface
//- Owner: Steven J. Owen
//- Description: Maintains a list of mesh.  This is used to store
//-        the exterior skin while creating a geometry from a mesh.
//- Created: 5/25/00
//- Checked By:
//- Version:

#ifndef ChollaSurface_HPP
#define ChollaSurface_HPP


#include "DLIList.hpp"
#include "ChollaEntity.hpp"

class ChollaCurve;
class FacetEntity;
class CubitFacet;
class CubitPoint;
class CubitFacetEdge;
class FacetEvalTool;

class ChollaSurface : public ChollaEntity
{
private:
   
  int id;
  int blockId;

  DLIList<FacetEntity*> surfaceElemList;
  DLIList<ChollaCurve*> curveList;
  void *mySurface;
  FacetEvalTool *myEvalTool;

  CubitStatus get_adj_facets( FacetEntity *start_face_ptr, 
                             DLIList<FacetEntity*> &face_list,
                             int mydebug = 0);
    // recursive function that creates a list of all faces connected
    // the passed in face that are part of this surface (or shell)

  void check_faceting();
  
public:
   
  ChollaSurface(int block_id);
  ~ChollaSurface();

  void get_facets(DLIList<FacetEntity*> &surface_facets)
    {surface_facets += surfaceElemList;}
  void get_points( DLIList<CubitPoint *> &point_list );

  void assign_geometric_surface(void *surf)
    {mySurface = surf;}
  void* get_geometric_surface()
    {return mySurface;}
  void assign_eval_tool(FacetEvalTool *eval_tool_ptr)
    {myEvalTool = eval_tool_ptr;}
  FacetEvalTool* get_eval_tool()
    {return myEvalTool;}

  void add_facet(FacetEntity *exterior_face)
    {surfaceElemList.append(exterior_face);}
  int add_mesh_unique(FacetEntity *exterior_face)
    {return surfaceElemList.append_unique(exterior_face);} 
  DLIList<FacetEntity*> &get_facet_list()
    {return surfaceElemList;}
  DLIList<FacetEntity*> *get_facet_list_ptr()
    {return &surfaceElemList;}
  void get_curves( DLIList<ChollaCurve*> &bcm_list )
    {bcm_list = curveList; }
  void add_curve( ChollaCurve *bcm_ptr )
    {curveList.append(bcm_ptr);}
  void add_curve_unique( ChollaCurve *bcm_ptr )
    {curveList.append_unique(bcm_ptr);}
  void remove_curve( ChollaCurve *bcm_ptr)
    {curveList.remove(bcm_ptr);}
  int get_block_id()
    {return blockId;}
  void set_block_id(int flag)
    { blockId = flag; }

  int get_id(){return id;}

  CubitStatus split_surface( DLIList<ChollaSurface*> &block_surface_list );
    // split this surface into multiple ChollaSurface where there are
    // discontinuous faces.

  CubitStatus feature_angle( double min_dot, 
                             DLIList<CubitFacetEdge *> &feature_edge_list);
    // mark all edges that exceed the specified feature angle
    // min_dot is the minimum dot product between adjacent face normals

  CubitStatus non_manifold_edges( DLIList<CubitFacetEdge *> &feature_edge_list );
    // mark all edges that are non-manifold (have more than 2 adj facets

  CubitStatus clean_features( );
    // clean up edges that do not form complete loops as a result 
    // of feature angle
  
  void init_hit_flags();
    // initialize the hit flags to the block surface id - used for 
    // traversing the faces

  CubitStatus update_boundary_tool_data();
    //- update the surface IDs on the boundary facet tooldatas

  void reset_facet_flags();
    //- reset the marked flags on the facets

  CubitBoolean is_adjacent( ChollaSurface *other_surf );
    //- determine if other_surf is adjacent to this surface

  DLIList<DLIList<CubitFacetEdge *>*> *get_loop_edges( );
    // return the ordered list of edges on the boundary of this surface 

  void debug_draw();
};

#endif





