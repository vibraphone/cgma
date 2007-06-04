/** \class  TDFacetBoundaryPoint
//- Class:       TDFacetBoundaryPoint
//- Description: Tool data for storing additional information at 
//-              the boundary of a facet set
//- Owner:       Steve Owen
//- Checked by:
//- Version:
*/
#ifndef TD_FACET_BOUNDARY_POINT_HPP


#define TD_FACET_BOUNDARY_POINT_HPP

#include "CubitDefines.h"
#include "ToolData.hpp"
#include "MemoryManager.hpp"
#include "DLIList.hpp"
#include "CubitVector.hpp"
#include "CastTo.hpp"
class CubitPoint;
class CubitFacet;
class CubitQuadFacet;
class CubitFacetEdge;
class CubitTransformMatrix;
/** \struct BoundaryPointData
    No comments.
*/
typedef struct BoundaryPointData
{
  DLIList<CubitFacet *> surfFacetList;
  CubitVector normal;
  double uVal, vVal, sizeVal;
  int surfID;
} BoundaryPointData;

class TDFacetBoundaryPoint : public ToolData
{
private:

  static MemoryManager memoryManager;
    //- memory management object

  CubitPoint *pointPtr;
  DLIList <BoundaryPointData *> pointDataList;

  void init_normal( BoundaryPointData *bpd_ptr );
    // initialize normals at point
  BoundaryPointData *get_bpd( CubitFacet *facet );
    // return a boundary point data

public:

  TDFacetBoundaryPoint();
    //- constructor

  ~TDFacetBoundaryPoint();

  static int is_facet_boundary_point(const ToolData* td)
     {return (CAST_TO(const_cast<ToolData*>(td), TDFacetBoundaryPoint) != NULL);}
  
  void add_surf(int new_id);

  CubitPoint *get_point()
    { return pointPtr; }
 
  void set_point( CubitPoint *point_ptr )
    { pointPtr = point_ptr; }

  SetDynamicMemoryAllocation(memoryManager)
    //- class specific new and delete operators
    
  static void set_memory_allocation_increment(int increment = 0)
    {memoryManager.set_memory_allocation_increment(increment);}
    //- set block memory size increment
  
  static void destroy_memory()
    {memoryManager.destroy_memory();}
    //- destroy all memory allocted to this object

  static CubitStatus add_facet_boundary_point(CubitPoint *point_ptr);
    // basic initializer
  
  static CubitStatus add_facet_boundary_point( CubitPoint *point_ptr,
                                               CubitFacet *facet_ptr,
                                               CubitVector &normal );
  static CubitStatus add_facet_boundary_point( CubitPoint *point_ptr,
                                               CubitQuadFacet *qfacet_ptr,
                                               CubitVector &normal );
    // initializes with a single facet and normal

  static TDFacetBoundaryPoint* get_facet_boundary_point(CubitPoint *point_ptr);
    // get the fbp from the point

  CubitStatus get_normal( int surf_id, CubitVector &normal );
  CubitStatus get_normal( CubitFacet *adj_facet, CubitVector &normal );
  CubitStatus get_normal( CubitFacetEdge *edge_ptr, CubitVector &normal );
  CubitStatus set_normal( int surf_id, CubitVector &normal );
    //recalculate normals.  Call this if an attached facet changes in some way.
  CubitStatus reset_normals();
    //- get and set normals at a specific surface or facet

  CubitStatus set_uv( CubitFacet *adj_facet, double u, double v );
  CubitStatus set_uvs( CubitFacet *adj_facet, double u, double v, double s );
  double u( CubitFacet *adj_facet );
  double v( CubitFacet *adj_facet );
  double s( CubitFacet *adj_facet );
  CubitStatus get_uv( CubitFacet *adj_facet, double &u, double &v ); 
    //- get and set the u-v values for this surface
  CubitStatus get_uvs( CubitFacet *adj_facet, double &u, double &v, double &s ); 
    //- get and set the u-v values and size for this surface


  void add_surf_facets( DLIList<CubitFacet *> adj_facet_list );
    //- add a group of facets that are adjacent to this point
    //- that are on the same surface
  void set_surf_id( CubitFacet *facet_ptr, int surf_id );
    //- set id of surface associated with facet

  CubitStatus merge_normals( CubitFacet *facet0, CubitFacet *facet1);
    //- merge the normals from facets

  CubitStatus rotate_normal( CubitTransformMatrix &rotmat );
    //- rotate the normals at the point

  CubitStatus get_boundary_point_data_size( int &size_int_data,
                                            int &size_double_data );
  CubitStatus get_boundary_point_data(int *int_data, double *double_data,
                                      int &iidx, int &didx );
   //- retreive info for dumping to a CUB file

  static CubitStatus new_facet_boundary_point(CubitPoint **points,
                    CubitFacet **facets,int &iidx,int &didx,
                    int *int_data,double *double_data);
  void initialize(CubitFacet **facets, 
            int &iidx, int &didx, int *int_data,double *double_data);
    //- restore data from a CUB file
};
    

#endif // TD_FACET_BOUNDARY_EDGE_HPP


