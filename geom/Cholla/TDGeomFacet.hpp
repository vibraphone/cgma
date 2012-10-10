//- Class: TDGeomFacet
//- Owner: Steve Owen
//- Description: Data for reading in a mesh and creating a geometry.
//-              This is assumed to be stored on a facet entity.
//- Checked By: 
//- Version:

#ifndef TD_GEOM_MESH_HPP


#define TD_GEOM_MESH_HPP

#include "CubitDefines.h"
#include "DLIList.hpp"
#include "ToolData.hpp"
#include "CubitPoint.hpp"
class FacetEntity;
class ChollaEntity;
class ChollaVolume;
class ChollaSurface;
class ChollaCurve;
class ChollaPoint;
class CubitFacet;
class CubitFacetEdge;


class TDGeomFacet : public ToolData
{
private:

  static MemoryManager memoryManager;
    //- memory management object

  int blockId;

  int hitFlag;
    //- Used for looping.

  DLIList<ChollaSurface*> ChollaSurfaceList;
    //- list of surfaces this element is on.

  DLIList<ChollaCurve*> ChollaCurveList;
    //- list of curves this element (edge) is on
  
  DLIList<ChollaPoint*> ChollaPointList;
    //- list chollapoints this element (point) is on 
  
  DLIList<CubitPoint*> myPoints;
    //- list of associated points for this entity (could be one).

  CubitVector normal;
    //- normal or tangent vector of a facet or facet edge

  DLIList<CubitFacetEdge *> *partnerEdgeList;
    //- other edges that share the same geometry

  DLIList<CubitPoint *> *partnerPointList;
    //- other points that share the same geometry

public:

  TDGeomFacet();
    //- constructor

  ~TDGeomFacet();

  static int is_geom_facet(const ToolData* td);
  
  int get_block_id()
    {return blockId;}

  void set_block_id(int new_id)
    {blockId = new_id;}

  int get_hit_flag()
    {return hitFlag;}

  void set_hit_flag(int flag)
    {hitFlag = flag;}
  
  void add_cholla_owner(ChollaEntity *cholla_entity);

  void add_cholla_surf(ChollaSurface *f_s_m);

  void get_cholla_surfs(DLIList<ChollaSurface*> &surf_list)
    {surf_list =  ChollaSurfaceList;}

  void remove_cholla_surfs()
    {ChollaSurfaceList.clean_out();}

  void remove_cholla_surf( ChollaSurface *chsurf_ptr )
    {ChollaSurfaceList.remove( chsurf_ptr );}

  void add_cholla_curve(ChollaCurve *chcurv_ptr)
    {ChollaCurveList.append_unique( chcurv_ptr );}

  void get_cholla_curves(DLIList<ChollaCurve*> &curv_list)
    {curv_list =  ChollaCurveList;}

  void remove_cholla_curves()
    {ChollaCurveList.clean_out();}

  void remove_cholla_curve( ChollaCurve *chcurv_ptr )
    {ChollaCurveList.remove( chcurv_ptr );}
  
  void add_cholla_point(ChollaPoint *chpt_ptr)
    {ChollaPointList.append_unique( chpt_ptr );}

  void get_cholla_points(DLIList<ChollaPoint*> &point_list)
    {point_list =  ChollaPointList;}
  
  void remove_cholla_points()
    {ChollaPointList.clean_out();}
  
  void remove_cholla_point( ChollaPoint *point_ptr )
    { ChollaPointList.remove( point_ptr ); }
  
  void add_point(CubitPoint *point)
    {myPoints.append(point);}

  void delete_point(CubitPoint *point)
    {myPoints.omit(point);}

  void get_points(DLIList<CubitPoint*> &point_list)
    {point_list += myPoints;}

  CubitPoint* get_first_point()
    {
      myPoints.reset();
      if ( myPoints.size() )
        return myPoints.get();
      else
        return (CubitPoint*) NULL;
    }

  void add_partner_edge( CubitFacetEdge *partner )
    { if (partnerEdgeList == NULL)
        partnerEdgeList = new DLIList<CubitFacetEdge *>;
      partnerEdgeList->append( partner ); };

  void get_partner_edges( DLIList <CubitFacetEdge *> &partner_list )
    { if (partnerEdgeList != NULL) partner_list += *partnerEdgeList; };

  int num_partner_edges( )
    { if (partnerEdgeList == NULL)
        return 0;
      else
        return partnerEdgeList->size(); }

  void add_partner_point( CubitPoint *partner )
    { if (partnerPointList == NULL)
        partnerPointList = new DLIList<CubitPoint *>;
      partnerPointList->append( partner ); };

  void get_partner_points( DLIList <CubitPoint *> &partner_list )
    { if (partnerPointList != NULL) partner_list += *partnerPointList; };

  int num_partner_points( )
    { if (partnerPointList == NULL)
        return 0;
      else
        return partnerPointList->size(); }
  
  CubitBoolean is_in_volume( ChollaVolume *chvol_ptr );

  SetDynamicMemoryAllocation(memoryManager)
    //- class specific new and delete operators
    
  static void set_memory_allocation_increment(int increment = 0)
    {memoryManager.set_memory_allocation_increment(increment);}
    //- set block memory size increment
  
  static void destroy_memory()
    {memoryManager.destroy_memory();}
    //- destroy all memory allocted to this object

  static CubitStatus add_geom_facet(FacetEntity *facet_entity, int block_id);
  static CubitStatus add_geom_facet(CubitFacet *facet_ptr, int block_id);
  static CubitStatus add_geom_facet(CubitFacetEdge *edge_ptr, int block_id);
  static CubitStatus add_geom_facet(CubitPoint *point_ptr, int block_id);
  
  static TDGeomFacet* get_geom_facet(FacetEntity *facet_entity);
  static TDGeomFacet* get_geom_facet(CubitFacet *facet_ptr);
  static TDGeomFacet* get_geom_facet(CubitFacetEdge *edge_ptr);
  static TDGeomFacet* get_geom_facet(CubitPoint *point_ptr);
  
  static int get_block_id(FacetEntity *mesh_entity);
  static int get_block_id(CubitFacet *facet_ptr);
  static int get_block_id(CubitFacetEdge *edge_ptr);
  static int get_block_id(CubitPoint *point_ptr);

  static int get_hit_flag(FacetEntity *mesh_entity);
  static void set_hit_flag(FacetEntity *mesh_entity, int new_val);
  CubitVector get_normal( ) { return normal; }
  void set_normal( CubitVector &norm ) { normal = norm; }

  CubitBoolean is_partner( CubitFacetEdge *edge_ptr );
  CubitBoolean is_partner( CubitPoint *point_ptr );

  void reset_TD_as_new();
  
  int geo_type();
  
};
    

#endif // TD_GEOM_MESH_HPP


