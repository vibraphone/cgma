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

class ChollaVolume;
class ChollaCurve;
class ChollaPoint;
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
  CubitBoolean myFlag;

  DLIList<FacetEntity*> surfaceElemList;
  DLIList<ChollaCurve*> curveList;
  DLIList<ChollaVolume*> volList;
  void *mySurface;
  FacetEvalTool *myEvalTool;
  ChollaSurface *myMergePartner;
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
  
  ChollaSurface *merge_parter(){ return myMergePartner; }
  void set_merge_partner( ChollaSurface *merge_partner )
    { myMergePartner = merge_partner;}
  
  CubitBoolean get_flag(){ return myFlag; }
  void set_flag( CubitBoolean stat ){ myFlag = stat; }                       

  void add_facet(FacetEntity *exterior_face)
    {surfaceElemList.append(exterior_face);}

  int add_mesh_unique(FacetEntity *exterior_face)
    {return surfaceElemList.append_unique(exterior_face);} 
  
  void remove_facet( FacetEntity *facet )
  {surfaceElemList.remove( facet );}
  
  bool is_contain( FacetEntity *facet );

  DLIList<FacetEntity*> &get_facet_list()
    {return surfaceElemList;}
  DLIList<FacetEntity*> *get_facet_list_ptr()
    {return &surfaceElemList;}
  void get_vertices( DLIList<ChollaPoint *> &chpt_list );
  void get_curves( DLIList<ChollaCurve*> &bcm_list )
    {bcm_list = curveList; }
  void add_curve( ChollaCurve *bcm_ptr )
    {curveList.append(bcm_ptr);}
  void add_curve_unique( ChollaCurve *bcm_ptr )
    {curveList.append_unique(bcm_ptr);}
  void remove_curve( ChollaCurve *bcm_ptr)
    {curveList.remove(bcm_ptr);}
  void get_volumes( DLIList<ChollaVolume*> &cholla_vol_list )
    {cholla_vol_list = volList; }
  void add_volume( ChollaVolume *cholla_vol_ptr )
    {volList.append(cholla_vol_ptr);}
  void add_volume_unique( ChollaVolume *cholla_vol_ptr )
    {volList.append_unique(cholla_vol_ptr);}
  void remove_volume( ChollaVolume *cholla_vol_ptr)
    {volList.remove(cholla_vol_ptr);}
  int num_volumes(){return volList.size();}
  CubitBoolean is_in_volume( ChollaVolume *chvol_ptr );
  
  int get_block_id()
    {return blockId;}
  void set_block_id(int flag)
    { blockId = flag; }

  int get_id(){return id;}

  CubitStatus get_adj_facets( FacetEntity *start_face_ptr, 
                             DLIList<FacetEntity*> &face_list,
                             int mydebug = 0,
                             bool bound_check = false,
                             bool feature_edge_check = true);
    // recursive function that creates a list of all faces connected
    // the passed in face that are part of this surface (or shell)


  CubitStatus split_surface( DLIList<ChollaSurface*> &block_surface_list );
    // split this surface into multiple ChollaSurface where there are
    // discontinuous faces.

  CubitStatus feature_angle( double min_dot, 
                             DLIList<CubitFacetEdge *> &feature_edge_list);
    // mark all edges that exceed the specified feature angle
    // min_dot is the minimum dot product between adjacent face normals
	
	CubitStatus add_preexisting_feature_edges( DLIList<CubitFacetEdge *> &feature_edge_list);
	  // edges that were marked previously in function ChollaEngine::mark_features
	  // are added to the feature edge list

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
  
  void flip_facets();
    // invert all the facets on this surface
  
  void debug_draw();
};

#endif





