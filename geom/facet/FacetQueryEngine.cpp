//-------------------------------------------------------------------------
// Filename      : FacetQueryEngine.cpp
//
// Purpose       : Implementation of the FacetQueryEngine class.
//                 This class provides facet-based implementations
//                 of various virtual functions in the GeometryQueryEngine
//                 hierarchy.
//
// Special Notes :
//
// Creator       : David R. White
//
// Creation Date : 7/17/00
//
//-------------------------------------------------------------------------
#include "FacetQueryEngine.hpp"
#include "FacetModifyEngine.hpp"
#include "FacetboolInterface.hpp"
#include "TopologyEntity.hpp"
#include "TopologyBridge.hpp"
#include "RefEntity.hpp"
#include "Body.hpp"
#include "Shell.hpp"
#include "Loop.hpp"
#include "Chain.hpp"
#include "CoEdge.hpp"
#include "CoFace.hpp"
#include "RefVolume.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include "GeometryEntity.hpp"
#include "DLIList.hpp"
#include "CubitBox.hpp"
#include "CubitString.hpp"
#include "FacetPoint.hpp"
#include "FacetCurve.hpp"
#include "FacetCoEdge.hpp"
#include "FacetLoop.hpp"
#include "FacetSurface.hpp"
#include "FacetShell.hpp"
#include "FacetLump.hpp"
#include "FacetBody.hpp"
#include "CubitFacetEdge.hpp"
#include "CubitFacetEdgeData.hpp"
#include "CubitFacet.hpp"
#include "CubitFacetData.hpp"
#include "CubitQuadFacet.hpp"
#include "CubitQuadFacetData.hpp"
#include "CubitPoint.hpp"
#include "GMem.hpp"
#include "FacetEvalTool.hpp"
#include "CurveFacetEvalTool.hpp"
#include "CubitPointData.hpp"
#include "GeometryQueryTool.hpp"
#include "debug.hpp"
#include "CubitObserver.hpp"
#include "ChollaEngine.hpp"
#include "ChollaSurface.hpp"
#include "ChollaCurve.hpp"
#include "ChollaPoint.hpp"
#include "Cholla.h"
#include "CubitFileIOWrapper.hpp"
#include "CubitFile.hpp"
#include "TDFacetBoundaryPoint.hpp"
#include "GfxDebug.hpp"
#include "KDDTree.hpp"
#include "RTree.hpp"
#include "FacetDataUtil.hpp"
#include "GridSearchTree.hpp"
#include <stdio.h>
#include <errno.h>

using namespace NCubitFile;

FacetQueryEngine* FacetQueryEngine::instance_ = NULL;
int FacetQueryEngine::hashPointSize = 0;
DLIList<CubitPoint*> *FacetQueryEngine::hashPointArray = NULL;

const int FacetQueryEngine::FQE_MAJOR_VERSION = 10;
const int FacetQueryEngine::FQE_MINOR_VERSION = 0;
const int FacetQueryEngine::FQE_SUBMINOR_VERSION = 0;

//================================================================================
// Description:
// Author     :
// Date       :
//================================================================================
FacetQueryEngine* FacetQueryEngine::instance()
{
  if (instance_ == NULL ) {
      instance_ = new FacetQueryEngine;
   }
  return instance_;
}

//================================================================================
// Description:  default constructor
// Author     :
// Date       :
//================================================================================
FacetQueryEngine::FacetQueryEngine()
{
  GeometryQueryTool::instance()->add_gqe( this );
}

//================================================================================
// Description:  destructor
// Author     :
// Date       :
//================================================================================
FacetQueryEngine::~FacetQueryEngine()
{
  instance_ = NULL;
}

//================================================================================
// Description:  can_delete_bodies
// Author     : sjowen
// Date       : 4/25/02
//================================================================================
CubitBoolean FacetQueryEngine::can_delete_bodies(DLIList<Body*>body_list)
{
  CubitBoolean delete_ok = CUBIT_TRUE;
  int ii;
  for (ii=0; ii<body_list.size() && delete_ok; ii++)
  {
    Body *body_ptr = body_list.get_and_step();
    // Extract the BODY from Body
    FacetBody* fbody_ptr = CAST_TO(body_ptr->get_body_sm_ptr(), FacetBody);
    if (fbody_ptr)
    {
      delete_ok = fbody_ptr->can_be_deleted(body_list);
    }
  }
  return delete_ok;
}

int FacetQueryEngine::get_major_version()
{
  return FQE_MAJOR_VERSION;
}

int FacetQueryEngine::get_minor_version()
{
  return FQE_MINOR_VERSION;
}

int FacetQueryEngine::get_subminor_version()
{
  return FQE_SUBMINOR_VERSION;
}

CubitString FacetQueryEngine::get_engine_version_string()
{
  CubitString version_string = "Facet Geometry Engine version ";
  version_string += CubitString(get_major_version());
  version_string += CubitString(".");
  version_string += CubitString(get_minor_version());
  version_string += CubitString(".");
  version_string += CubitString(get_subminor_version());
  
  return version_string;
}

//================================================================================
// Description:
// Author     :
// Date       :
//================================================================================
const type_info& FacetQueryEngine::entity_type_info() const
{
   return typeid(*this);
}

//================================================================================
// Description:
// Author     :
// Date       :
//================================================================================
Point* FacetQueryEngine::make_Point( GeometryType ,
                                        CubitVector const& ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return (Point*) NULL;
}

//================================================================================
// Description:
// Author     :
// Date       :
//================================================================================
Curve* FacetQueryEngine::make_Curve(Curve *) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return (Curve*) NULL;
}

//================================================================================
// Description:
// Author     :
// Date       :
//================================================================================
Curve* FacetQueryEngine::make_Curve( Point const* ,
                                        Point const* ,
                                        RefFace* ,
                                        CubitVector * ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return (Curve*) NULL;
}

//================================================================================
// Description:
// Author     :
// Date       :
//================================================================================
Curve* FacetQueryEngine::make_Curve( GeometryType ,
                                        Point const* ,
                                        Point const* ,
                                        DLIList<CubitVector*>& ,
                                        RefFace*  ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return (Curve*) NULL;
}

//================================================================================
// Description:
// Author     :
// Date       :
//================================================================================
Curve* FacetQueryEngine::make_Curve( GeometryType ,
                                        Point const* ,
                                        Point const* ,
                                        CubitVector const* ,
                                        CubitSense ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return (Curve*) NULL;
}

//================================================================================
// Description:
// Author     :
// Date       :
//================================================================================
Surface* FacetQueryEngine::make_Surface( Surface *,
                                            DLIList<Loop*> &,
                                            CubitBoolean  ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return (Surface*) NULL;
}

//================================================================================
// Description:
// Author     :
// Date       :
//================================================================================
Surface* FacetQueryEngine::make_Surface( GeometryType ,
                                            DLIList<Curve*>& ,
                                            DLIList<Loop*> &,
                                            Surface *) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return (Surface*) NULL;
}

//================================================================================
// Description:
// Author     :
// Date       :
//================================================================================
Lump* FacetQueryEngine::make_Lump( GeometryType ,
                                      DLIList<Surface*>&  ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return (Lump*) NULL;
}

//================================================================================
// Description:
// Author     :
// Date       :
//================================================================================
BodySM* FacetQueryEngine::make_BodySM( Surface * ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return (BodySM*) NULL;
}

//================================================================================
// Description:
// Author     :
// Date       :
//================================================================================
BodySM* FacetQueryEngine::make_BodySM( DLIList<Lump*>&  ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return (BodySM*) NULL;
}

//================================================================================
// Description: create a new body by copying
// Author     : sjowen
// Date       : 9/7/01
//================================================================================
Body* FacetQueryEngine::copy_body( Body *body_ptr )
{
  BodySM* bodysm_ptr = body_ptr->get_body_sm_ptr();
  FacetBody *facet_body_ptr = CAST_TO(bodysm_ptr, FacetBody);
  if (!facet_body_ptr)
  {
    PRINT_ERROR("Attempt to copy mesh-based geometry Body.  This body is not MBG.");
    return (Body*)NULL;
  }
  BodySM* osme_body_ptr = (BodySM*)facet_body_ptr->copy();
  if (!osme_body_ptr)
  {
    PRINT_ERROR("Failed to copy mesh-based geometry Body");
    return (Body*)NULL;
  }
  Body *new_body_ptr = GeometryQueryTool::instance()->make_Body( osme_body_ptr );
  return new_body_ptr;
}

//================================================================================
// Description: reflect body about an axis
// Author     : sjowen
// Date       : 9/7/01
//================================================================================
CubitStatus FacetQueryEngine::reflect( BodySM *bodysm,
                                       const CubitVector& axis)
{
  FacetBody *fbody = CAST_TO(bodysm, FacetBody);
  if (!fbody)
  {
    PRINT_ERROR("Attempt to reflect mesh-based geometry Body.  This body is not MBG.");
    return CUBIT_FAILURE;
  }

  fbody->reflect( axis.x(), axis.y(), axis.z() );

  return CUBIT_SUCCESS;
}

CubitStatus FacetQueryEngine::get_connected_patch(
                      DLIList<FacetSurface*>& remaining_surfs,
                      DLIList<FacetSurface*>& output_patch )
{
  DLIList<FacetSurface*> stack, curve_surfs;
  std::set<FacetSurface*> marked;
  DLIList<FacetCurve*> curves;
  int debug_this = 0;
  if(debug_this){
    GfxDebug::clear();
  }
  int i;
  for (i = remaining_surfs.size(); i--; )
    marked.insert( remaining_surfs.get_and_step() );

  // Choose one surface to start with
  FacetSurface* surf = remaining_surfs.pop();
  marked.erase( surf );
  stack.append( surf );

  // Get all connected surfaces
  while( stack.size() )
  {
    surf = stack.pop();
    if(debug_this){
      surf->get_eval_tool()->debug_draw_facets(CUBIT_BLUE);
      GfxDebug::mouse_xforms();
    }
    output_patch.append( surf );

    surf->get_curves( curves );
    while (curves.size())
    {
      FacetCurve* curve = curves.pop();
      curve->get_surfaces( curve_surfs );
      if(debug_this && curve_surfs.size() < 2){
        PRINT_INFO("Curve is not connected to 2 surfaces.\n");
        curve->get_eval_tool()->debug_draw_facet_edges(CUBIT_MAGENTA);
      }
      
        
      while (curve_surfs.size())
      {
        FacetSurface* curve_surf = curve_surfs.pop();
        if (marked.erase(curve_surf))
          stack.append( curve_surf );
      }
    }
  }

  // Remove output surfaces from input list
  // At this point, "marked" contains all the surfaces not in "ouput_patch"
  remaining_surfs.last();
  for (i = remaining_surfs.size(); i--; )
    if (marked.find(remaining_surfs.step_and_get()) == marked.end())
      remaining_surfs.change_to( 0 );
    remaining_surfs.remove_all_with_value( 0 );

  return CUBIT_SUCCESS;

}


//================================================================================
// Description:
// Author     :
// Date       :
//================================================================================
CubitStatus FacetQueryEngine::get_graphics( Surface* surface_ptr,
                                                     int& number_triangles,
                                                     int& number_points,
                                                     int& number_facets,
                                                     GMem* gMem,
                                                     unsigned short ,
                                                     double, double ) const
{
    //  get the FacetSurface.
  int i;
  FacetSurface *facet_surf_ptr = CAST_TO(surface_ptr,FacetSurface);
  if( ! facet_surf_ptr )
    return CUBIT_FAILURE;

    // get the facets for the surface
  DLIList<CubitFacet*> surface_facets;
  DLIList<CubitPoint*> surface_points;
  facet_surf_ptr->get_my_facets(surface_facets, surface_points);

    // set return values, and if GMem is NULL return
    // (caller just wanted to know the counts.)
  number_facets = surface_facets.size() * 4;
  number_points = surface_points.size();
  number_triangles = surface_facets.size();
  if( !gMem )
    return CUBIT_SUCCESS;


    // Allocate space in GMem
  gMem->allocate_tri(surface_facets.size());
  gMem->fListCount = number_facets;
  gMem->pointListCount = number_points;

    // Copy points to GMem
  surface_points.reset();
  GPoint* pt_itor = gMem->point_list();
  for ( i = 0; i < number_points; i++ )
  {
    CubitPoint* point = surface_points.get_and_step();
    pt_itor->x = (float)point->x();
    pt_itor->y = (float)point->y();
    pt_itor->z = (float)point->z();
    point->marked(i);
    pt_itor++;
  }

    // Copy facets to Gmem
  surface_facets.reset();
  int* f_itor = gMem->facet_list();
  for ( i = 0; i < number_triangles; i++ )
  {
    *(f_itor++) = 3;
    CubitFacet* facet = surface_facets.get_and_step();
    for ( int j = 0; j < 3; j++ )
      *(f_itor++) = facet->point(j)->marked();
  }

    // cleanup
  while( surface_points.size() )
    surface_points.pop()->marked(0);

  return CUBIT_SUCCESS;
}

//================================================================================
// Description:
// Author     :
// Date       :
//================================================================================
CubitStatus FacetQueryEngine::get_graphics( Curve* curve_ptr,
                                            int& num_points,
                                            GMem* gMem,
                                            double /*tolerance*/ ) const
{
  //  get the FacetCurve.
  FacetCurve *facet_curv_ptr = CAST_TO(curve_ptr,FacetCurve);

  DLIList<CubitFacetEdge*> curve_facets;
  facet_curv_ptr->get_facets(curve_facets);
  int number_facets = curve_facets.size();
  DLIList<CubitPoint*> curve_points;
  facet_curv_ptr->get_points(curve_points);
  curve_points.reset();
  curve_facets.reset();
  num_points = curve_points.size();
  GPoint *new_point_list = new GPoint[num_points];
  int *new_facet_list = new int [number_facets*3];
  int ii;
  for ( ii = 0; ii < num_points; ii++ )
  {
    new_point_list[ii].x = curve_points.get()->x();
    new_point_list[ii].y = curve_points.get()->y();
    new_point_list[ii].z = curve_points.get()->z();
      //mark the point with the index into the point list.
      //This is very important to make sure we can index these
      //points when we do the facet array.
    curve_points.get_and_step()->marked(ii);
  }
  for ( ii = 0; ii < number_facets; ii++ )
  {
      //All our facets are segments.  So the first value is 2.
    int index_count = 3*ii;
    new_facet_list[index_count] = 2;
    int jj = index_count + 1;
    int ll, kk;
    CubitFacetEdge *facet_ptr = curve_facets.get_and_step();
    CubitPoint *temp_point;
    for ( kk = jj, ll = 0; kk < jj+2; kk++, ll++)
    {
      temp_point = facet_ptr->point(ll);
      int index_to_point = temp_point->marked();
      if ( index_to_point < 0 || index_to_point > num_points )
      {
        PRINT_ERROR("Bad logic in converting edge facets to drawing facets.\n");
        return CUBIT_FAILURE;
      }
      new_facet_list[kk] = index_to_point;
    }
  }
  gMem->replace_point_list(new_point_list, num_points, num_points);
  gMem->replace_facet_list(new_facet_list, number_facets*3, number_facets*3);
  gMem->points_consolidated(CUBIT_TRUE);

  return CUBIT_SUCCESS;
}

//================================================================================
// Description:
// Author     :
// Date       :
//================================================================================
CubitStatus FacetQueryEngine::get_isoparametric_points(Surface* ,
                                                          int &nu, int &nv,
                                                          GMem*&) const
{
  nu = nv = 0;
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

CubitStatus FacetQueryEngine::get_u_isoparametric_points(Surface* ,
                                                            double, int&,
                                                            GMem*&) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

CubitStatus FacetQueryEngine::get_v_isoparametric_points(Surface* ,
                                                            double, int&,
                                                            GMem*&) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//================================================================================
// Description:
// Author     :
// Date       :
//================================================================================
CubitStatus FacetQueryEngine::transform_vec_position( CubitVector const& ,
                                                         BodySM *,
                                                         CubitVector & ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//================================================================================
// Description:
// Author     :
// Date       :
//================================================================================
CubitStatus FacetQueryEngine::get_intersections( Curve* , CubitVector& ,
                                                 CubitVector&,
                                                    DLIList<CubitVector*>& ,
                                                    CubitBoolean,
                                                    CubitBoolean )
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

CubitStatus FacetQueryEngine::get_intersections( Curve* , Curve* ,
                                                    DLIList<CubitVector*>& ,
                                                    CubitBoolean,
                                                    CubitBoolean )
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}
CubitStatus FacetQueryEngine::get_intersections( Curve*,
                                                   Curve*,
                                                   DLIList<CubitVector*>&,
                                                   double, CubitBoolean )
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

CubitStatus
FacetQueryEngine::get_intersections( Curve* /*ref_edge*/, Surface* /*ref_face*/,
                                        DLIList<CubitVector*>& ,
                                        CubitBoolean )
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
   return CUBIT_FAILURE;
}

//================================================================================
// Description: Find extrema position on an entity list
// Author     :
// Date       :
//================================================================================
CubitStatus
FacetQueryEngine::entity_extrema( DLIList<GeometryEntity*> &,
                                  const CubitVector *,
                                  const CubitVector *,
                                  const CubitVector *,
                                  CubitVector &,
                                  GeometryEntity *& )
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//================================================================================
// Description: Find distance between two entities and closest positions.
// Author     :
// Date       :
//================================================================================
CubitStatus
FacetQueryEngine::entity_entity_distance( GeometryEntity *,
                                          GeometryEntity *,
                                          CubitVector &, CubitVector &,
                                          double & )
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===========================================================================
//Function Name: save_temp_geom_file
//Member Type:  PUBLIC
//Description:  function called for save/restore to save temporary FACET file
//              If file is created, CubitString 'created_file' and
//              'create_file_type' are set
//Author:       Corey Ernst
//Date:         1/18/2003
//===========================================================================

CubitStatus FacetQueryEngine::save_temp_geom_file( DLIList<TopologyBridge*>& ref_entity_list,
                                                   const char *file_name,
                                                   const CubitString &cubit_version,
                                                   CubitString &created_file,
                                                   CubitString &created_file_type)
{
  int size_before = ref_entity_list.size();
  CubitString temp_filename(file_name);
  temp_filename += ".mbg";

  if( export_solid_model( ref_entity_list, temp_filename.c_str(), "FACET",
                          cubit_version ) == CUBIT_FAILURE )
  {
    PRINT_ERROR( "Error occured while trying to save temporary MESH_BASED_GEOMETRY file\n");
    return CUBIT_FAILURE;
  }

  int size_after = ref_entity_list.size();

  if( size_before > size_after )
  {
    created_file +=  temp_filename;
    created_file_type += "FACET";
  }
  return CUBIT_SUCCESS;
}

//===========================================================================
//Function Name:export_solid_model
//Member Type:  PUBLIC
//Description:  function called for save/restore to save temporary FACET file
//Author:       Corey Ernst
//Date:         1/18/2003
//===========================================================================

CubitStatus FacetQueryEngine::export_solid_model( DLIList<TopologyBridge*>& ref_entity_list,
                                                     const char* file_name,
                                                     const char* file_type,
                                                     const CubitString &,
                                                     const char*)
{
  if( strcmp( file_type, "FACET" ) )
    return CUBIT_SUCCESS;

  DLIList<FacetBody*>    facet_bodies;
  DLIList<FacetLump*>    facet_lumps;
  DLIList<FacetShell*>   facet_shells;
  DLIList<FacetSurface*> facet_surfaces;
  DLIList<FacetLoop*>    facet_loops;
  DLIList<FacetCoEdge*>  facet_coedges;
  DLIList<FacetCurve*>   facet_curves;
  DLIList<FacetPoint*>   facet_points;

  DLIList<TopologyBridge*> ref_entities_handled;

  int i;
  //Collect all free entities (bodies, curves, vertices )
  ref_entity_list.reset();
  for(i=ref_entity_list.size(); i>0; i--)
  {
    TopologyBridge* ref_entity_ptr = ref_entity_list.get();
    CubitBoolean handled = CUBIT_TRUE;

    //if it is a Vertex
    if( FacetPoint* pt = CAST_TO( ref_entity_ptr, FacetPoint) )
      facet_points.append( pt );
    //if it is a Curve
    else if( FacetCurve* curve = CAST_TO( ref_entity_ptr, FacetCurve) )
      facet_curves.append( curve );
    /*
    //if it is a Surface -- I don't think you can ever have a free surface
    //without it being a Body
    else if( FacetSurface* surf = CAST_TO( ref_entity_ptr, FacetSurface) )
      facet_surfaces.append( surf );
   */
    //if it is a Body
    else if( FacetBody* body = CAST_TO( ref_entity_ptr, FacetBody ) )
      facet_bodies.append( body );
    else
      handled = CUBIT_FALSE;

    if( handled == CUBIT_TRUE )
    {
      ref_entities_handled.append( ref_entity_ptr );
      ref_entity_list.change_to(NULL);
    }

    ref_entity_list.step();
  }

  ref_entity_list.remove_all_with_value(NULL);

  int free_body_count = facet_bodies.size();
  int free_curve_count = facet_curves.size();
  int free_point_count = facet_points.size();

  //if nothing to write out...return
  if( free_body_count == 0 && free_curve_count == 0 && free_point_count == 0)
    return CUBIT_SUCCESS;

  //get file pointer
  FILE *file_ptr = fopen( file_name, "wb" );

  //get all child entities from the bodies, curves, & vertices
  gather_all_facet_entities( facet_bodies, facet_lumps,
                             facet_shells, facet_surfaces,
                             facet_loops, facet_coedges,
                             facet_curves, facet_points );

  //create a wrapper object for writing
  CIOWrapper file_writer( file_ptr );

  // write out file type "MESHED_BASED_GEOMETRY"
  file_writer.BeginWriteBlock(0);
  file_writer.Write( "MESH_BASED_GEOMETRY", 19 );

  // write out Endian value
  UnsignedInt32 endian_value = CCubitFile::mintNativeEndian;
  file_writer.Write( &endian_value, 1 );

  // write out version #
  UnsignedInt32 version = 1;
  file_writer.Write( &version, 1 );


  //save the facets (geometry info )
  CubitStatus status;
  status = save_facets( file_ptr, facet_surfaces, facet_curves, facet_points );
  if( status == CUBIT_FAILURE ) return CUBIT_FAILURE;

  //write out topology and attributes
  status = write_topology( file_ptr,
                           facet_bodies, facet_lumps,
                           facet_shells, facet_surfaces,
                           facet_loops, facet_coedges,
                           facet_curves, facet_points );
  if( status == CUBIT_FAILURE ) return CUBIT_FAILURE;

  if( free_body_count || free_curve_count || free_point_count )
      PRINT_INFO( "\nExported:" );

   int flg = 0;

   if( free_body_count )
   {
      if( flg )PRINT_INFO( "         " );else flg=1;
      if( DEBUG_FLAG( 153 ) )
      {
        if( free_body_count == 1 )
           PRINT_INFO( "%4d Facet Body\n", free_body_count );
        else
           PRINT_INFO( "%4d Facet Bodies\n", free_body_count );
      }
      
      if( facet_lumps.size() == 1 )
         PRINT_INFO( "%4d Facet Volume\n", facet_lumps.size() );
      else
         PRINT_INFO( "%4d Facet Volumes\n", facet_lumps.size() );
   }
   if( free_curve_count )
   {
      if( flg )PRINT_INFO( "         " );else flg=1;
      if( free_curve_count == 1 )
         PRINT_INFO( "%4d Facet Curve\n", free_curve_count );
      else
         PRINT_INFO( "%4d Facet Curves\n", free_curve_count );
   }
   if( free_point_count )
   {
      if( flg )PRINT_INFO( "         " );else flg=1;
      if( free_point_count == 1 )
         PRINT_INFO( "%4d Facet Point\n", free_point_count );
      else
         PRINT_INFO( "%4d Facet Points\n", free_point_count );
   }
   PRINT_INFO( "\n" );



  fclose( file_ptr );

  return CUBIT_SUCCESS;
}

CubitStatus FacetQueryEngine::export_solid_model( DLIList<TopologyBridge*>& bridge_list,
						  char*& p_buffer,
						  int& n_buffer_size,
						  bool b_export_buffer)
{
  return CUBIT_FAILURE;
}

CubitStatus
FacetQueryEngine::gather_all_facet_entities( DLIList<FacetBody*> &facet_bodies,
                                             DLIList<FacetLump*> &facet_lumps,
                                             DLIList<FacetShell*> &facet_shells,
                                             DLIList<FacetSurface*> &facet_surfaces,
                                             DLIList<FacetLoop*> &facet_loops,
                                             DLIList<FacetCoEdge*> &facet_coedges,
                                             DLIList<FacetCurve*> &facet_curves,
                                             DLIList<FacetPoint*> &facet_points )
{
  int i;


  //Collect FacetLumps from FacetBody
  for(i=0; i<facet_bodies.size(); i++)
  {
    FacetBody *facet_body = facet_bodies.get_and_step();
    facet_body->get_lumps( facet_lumps );
  }


  //Collect FacetShells from FacetLumps
  for(i=0; i<facet_lumps.size(); i++)
  {
    FacetLump *facet_lump = facet_lumps.get_and_step();
    facet_lump->get_shells( facet_shells );
  }


  //Collect FacetSurfaces from FacetShells
  for(i=0; i<facet_shells.size(); i++)
  {
    FacetShell *facet_shell = facet_shells.get_and_step();
    facet_shell->get_surfaces( facet_surfaces );
  }
  facet_surfaces.uniquify_unordered();


  //Collect FacetLoops from FacetSurfaces
  for(i=0; i<facet_surfaces.size(); i++)
  {
    FacetSurface *facet_surface = facet_surfaces.get_and_step();
    facet_surface->get_loops( facet_loops );
  }

  //Collect FacetCoEdges from FacetLoops
  for(i=0; i<facet_loops.size(); i++)
  {
    FacetLoop *facet_loop = facet_loops.get_and_step();
    facet_loop->get_coedges( facet_coedges );
  }


  //Collect FacetCurves from FacetCoEdges
  for( i=0; i<facet_coedges.size(); i++)
  {
    FacetCoEdge *facet_coedge = facet_coedges.get_and_step();
    facet_coedge->get_curves( facet_curves );
  }
  facet_curves.uniquify_unordered();


  //Collect FacetPoints from FacetCurves
  for( i=0; i<facet_curves.size(); i++)
  {
    FacetCurve *facet_curve = facet_curves.get_and_step();
    facet_curve->get_points( facet_points );
  }

  //uniquify lists
  facet_points.uniquify_unordered();


  return CUBIT_SUCCESS;

}


CubitStatus
FacetQueryEngine::write_topology( FILE *file_ptr,
                                  DLIList<FacetBody*> &facet_bodies,
                                  DLIList<FacetLump*> &facet_lumps,
                                  DLIList<FacetShell*> &facet_shells,
                                  DLIList<FacetSurface*> &facet_surfaces,
                                  DLIList<FacetLoop*> &facet_loops,
                                  DLIList<FacetCoEdge*> &facet_coedges,
                                  DLIList<FacetCurve*> &facet_curves,
                                  DLIList<FacetPoint*> &facet_points )
{

  int i;

  //create a wrapper object for writing
  CIOWrapper file_writer( file_ptr );

  //-----------------write FacetPoints--------------
  UnsignedInt32 size = facet_points.size();
  //write out number of FacetPoints
  file_writer.Write( &size, 1 );
  facet_points.reset();
  for( i=0; i<facet_points.size(); i++)
  {
    FacetPoint *curr_point = facet_points.get_and_step();
    int id = curr_point->get_cubit_point()->id();

    file_writer.Write( reinterpret_cast<UnsignedInt32*>(&id), 1 );
    if( curr_point->save_attribs(file_ptr) == CUBIT_FAILURE )
      return CUBIT_FAILURE;
  }

  //-----------------write FacetCurves--------------
  size = facet_curves.size();
  //write out number of FacetCurves
  file_writer.Write( &size, 1 );
  facet_curves.reset();
  for( i=0; i<facet_curves.size(); i++)
  {
    FacetCurve *curr_curve = facet_curves.get_and_step();
    Point *s_point, *e_point;
    s_point = curr_curve->start_point();
    e_point = curr_curve->end_point();

    int data_to_write[4];

    // get start&end points implicit ids
    FacetPoint *temp_point = NULL;
    temp_point = CAST_TO( s_point, FacetPoint );
    if( !temp_point ) assert(0);
    int found;
    found = facet_points.where_is_item( temp_point );
    if( found == -1)
      assert(0);
    data_to_write[0] = found;

    temp_point = CAST_TO( e_point, FacetPoint );
    if( !temp_point ) assert(0);
    found = facet_points.where_is_item( temp_point );
    if( found == -1)
      PRINT_ERROR("Problem saving Facet Curves\n");
    data_to_write[1] = found;

    //convert Sense info to integer
    if( curr_curve->get_sense() == CUBIT_UNKNOWN )
      data_to_write[2] = -1;
    else
      data_to_write[2] = (curr_curve->get_sense() == CUBIT_REVERSED) ? 1 : 0;

    data_to_write[3] = curr_curve->get_eval_tool()->get_output_id();

    //write the data
    file_writer.Write( reinterpret_cast<UnsignedInt32*>(data_to_write), 4 );

    if( curr_curve->save_attribs(file_ptr) == CUBIT_FAILURE )
      return CUBIT_FAILURE;
  }

  //-----------------write FacetCoedges--------------
  size = facet_coedges.size();
  // write out number of FacetCurves
  file_writer.Write( &size, 1 );
  facet_coedges.reset();
  for( i=0; i<facet_coedges.size(); i++)
  {
    FacetCoEdge *curr_coedge = facet_coedges.get_and_step();
    Curve *curve_sm;
    curve_sm = curr_coedge->curve();

    FacetCurve *temp_curve = NULL;
    temp_curve = CAST_TO( curve_sm, FacetCurve );

    int data_to_write[2];

    // get implicit id of this curve
    int found;
    found = facet_curves.where_is_item( temp_curve );
    if( found == -1)
      PRINT_ERROR("Problem saving Facet CoEdges\n");
    data_to_write[0] = found;

    // convert sense info to integer
    if( curr_coedge->get_sense() == CUBIT_UNKNOWN )
      data_to_write[1] = -1;
    else
      data_to_write[1] = (curr_coedge->get_sense() == CUBIT_REVERSED) ? 1 : 0;

    // write out the data
    file_writer.Write( reinterpret_cast<UnsignedInt32*>(data_to_write), 2 );

  }

  //-----------------write FacetLoops--------------
  size = facet_loops.size();
  // write out number of FacetLoops
  file_writer.Write( &size, 1 );
  facet_loops.reset();
  for( i=0; i<facet_loops.size(); i++)
  {
    FacetLoop *curr_loop = facet_loops.get_and_step();
    DLIList<FacetCoEdge*> coedge_list;
    curr_loop->get_coedges( coedge_list );

    // get number of coedges in this loop
    UnsignedInt32 *data_to_write;
    size = coedge_list.size();
    data_to_write = new UnsignedInt32[ size + 1 ];
    data_to_write[0] = size;

    UnsignedInt32 j;
    // get implicit ids of coedges
    coedge_list.reset();
    for( j=1; j<size+1; j++)
    {
      FacetCoEdge *temp_coedge = coedge_list.get_and_step();
      int found;
      found = facet_coedges.where_is_item( temp_coedge );
      if( found == -1)
        PRINT_ERROR("Problem saving Facet Loops\n");
      data_to_write[j] = found;
    }

    // write out the data
    file_writer.Write( data_to_write, size + 1);
    delete [] data_to_write;
  }

  //-----------------write FacetSurfaces--------------
  size = facet_surfaces.size();
  // write out number of FacetSurfaces
  file_writer.Write( &size, 1 );
  facet_surfaces.reset();
  for( i=0; i<facet_surfaces.size(); i++)
  {
    FacetSurface *curr_surface = facet_surfaces.get_and_step();

    DLIList<FacetLoop*> loop_list;
    curr_surface->get_loops( loop_list );

    int num_loops = loop_list.size();
    int data_to_write[6];

    // convert sense info to integer
    // if( curr_surface->get_relative_surface_sense() == CUBIT_UNKNOWN )
//       data_to_write[0] = -1;
//     else
//       data_to_write[0] = (curr_surface->get_relative_surface_sense() == CUBIT_REVERSED) ? 1 : 0;
    data_to_write[0]=0;
    // get "useFacets"
    data_to_write[1] = 1;

    // get output id of FacetEvalTool
    data_to_write[2] = curr_surface->get_eval_tool()->get_output_id();

    // get Shell Sense stuff
    CubitSense sense0;
    
    curr_surface->get_shell_sense( sense0 );
    if( sense0 == CUBIT_UNKNOWN )
      data_to_write[3] = -1;
    else
      data_to_write[3] = (sense0 == CUBIT_REVERSED) ? 1 : 0;

//    if( sense1 == CUBIT_UNKNOWN )
      data_to_write[4] = -1;
//    else
//      data_to_write[4] = (sense1 == CUBIT_REVERSED) ? 1 : 0;

    // get number of loops
    data_to_write[5] = num_loops;

    file_writer.Write( reinterpret_cast<UnsignedInt32*>(data_to_write), 6 );

    // get implicit ids of loops
    if( num_loops > 0 )
    {
      int *loop_ids = new int[num_loops];
      int j;
      loop_list.reset();
      for( j=0; j<num_loops; j++)
      {
       FacetLoop *temp_loop = loop_list.get_and_step();
       int found;
       found = facet_loops.where_is_item( temp_loop );
       if( found == -1 )
         PRINT_ERROR("Problem saving Facet Surfaces\n");
       loop_ids[j] = found;
      }

      // write out data
      file_writer.Write( reinterpret_cast<UnsignedInt32*>(loop_ids), num_loops );
      delete [] loop_ids;
    }

    if( curr_surface->save_attribs(file_ptr) == CUBIT_FAILURE )
      return CUBIT_FAILURE;
  }

  //-----------------write FacetShells--------------
  size = facet_shells.size();
  // write out number of FacetShells
  file_writer.Write( &size, 1 );
  facet_shells.reset();
  for( i=0; i<facet_shells.size(); i++)
  {
    FacetShell *curr_shell= facet_shells.get_and_step(); //number of surfaces
    DLIList<FacetSurface*> temp_facet_surf_list;
    curr_shell->get_surfaces( temp_facet_surf_list );

    // get number of surfaces in this shell
    UnsignedInt32 *data_to_write;
    int num_surfs = temp_facet_surf_list.size();
    data_to_write = new UnsignedInt32[ num_surfs + 1];
    data_to_write[0] = num_surfs;

    // get implicit ids of surfaces
    int j;
    temp_facet_surf_list.reset();
    for( j=1; j<num_surfs+1; j++)
    {
      FacetSurface *temp_facet_surface = temp_facet_surf_list.get_and_step();
      int found;
      found = facet_surfaces.where_is_item( temp_facet_surface );
      if( found == -1 )
        PRINT_ERROR("Problem saving Facet Shells\n");
      data_to_write[j] = found;
    }

    // write the data
    file_writer.Write( data_to_write, num_surfs + 1 );
    delete [] data_to_write;
  }

  //-----------------write FacetLumps--------------
  size = facet_lumps.size();
  // write out number of FacetLumps
  file_writer.Write( &size, 1 );
  facet_lumps.reset();
  for( i=0; i<facet_lumps.size(); i++)
  {
    FacetLump *curr_lump = facet_lumps.get_and_step();

    DLIList<FacetShell*> temp_facet_shell_list;
    curr_lump->get_shells( temp_facet_shell_list );

    // get number of shells in this lump
    UnsignedInt32 *data_to_write;
    int num_shells= temp_facet_shell_list.size();
    data_to_write = new UnsignedInt32[ num_shells+ 1];
    data_to_write[0] = num_shells;

    //get implicit ides of the lumps in this shell
    int j;
    temp_facet_shell_list.reset();
    for( j=1; j<num_shells+1; j++)
    {
      FacetShell *temp_facet_shell = temp_facet_shell_list.get_and_step();
      int found;
      found = facet_shells.where_is_item( temp_facet_shell );
      if( found == -1 )
        PRINT_ERROR("Problem saving Facet Lumps\n");
      data_to_write[j] = found;
    }

    //write the data
    file_writer.Write( data_to_write, num_shells + 1 );
    delete [] data_to_write;
    if( curr_lump->save_attribs(file_ptr) == CUBIT_FAILURE )
      return CUBIT_FAILURE;
  }

  //-----------------write FacetBodies--------------
  size = facet_bodies.size();
  // write out number of FacetBodies
  file_writer.Write( &size, 1 );
  facet_bodies.reset();
  for( i=0; i<facet_bodies.size(); i++)
  {
    FacetBody *curr_body = facet_bodies.get_and_step();

    DLIList<FacetLump*> temp_facet_lump_list;
    curr_body->get_lumps( temp_facet_lump_list );

    // get the number of lumps in this body
    UnsignedInt32 *data_to_write;
    int num_lumps = temp_facet_lump_list.size();
    data_to_write = new UnsignedInt32[ num_lumps + 1];
    data_to_write[0] = num_lumps;

    // get the implicit ids of the lumps in this body
    int j;
    temp_facet_lump_list.reset();
    for( j=1; j<num_lumps+1; j++)
    {
      FacetLump *temp_facet_lump = temp_facet_lump_list.get_and_step();
      int found;
      found = facet_lumps.where_is_item( temp_facet_lump );
      if( found == -1 )
        PRINT_ERROR("Problem saving Facet Bodies\n");
      data_to_write[j] = found;
    }

    // write the data
    file_writer.Write( data_to_write, num_lumps + 1 );
    delete [] data_to_write;

    // write the transformation matrix of this body
    CubitTransformMatrix trans_matrix;
    curr_body->get_transforms( trans_matrix );

    UnsignedInt32 num_rows  = trans_matrix.num_rows();
    UnsignedInt32 num_cols  = trans_matrix.num_cols();
    UnsignedInt32 rows_and_cols[2];
    rows_and_cols[0] = num_rows;
    rows_and_cols[1] = num_cols;

    file_writer.Write( rows_and_cols, 2 );

    double *trans_matrix_array;
    trans_matrix_array = new double[ num_rows*num_cols ];

    //fill up the array row-by-row
    unsigned u, k = 0;
    for(u=0; u<num_rows; u++)
    {
      for(k=0; k<num_cols; k++)
        trans_matrix_array[(u*num_cols)+k] = trans_matrix.get(u,k);
    }

    file_writer.Write( trans_matrix_array, u*k );
    delete [] trans_matrix_array;
    if( curr_body->save_attribs(file_ptr) == CUBIT_FAILURE )
      return CUBIT_FAILURE;
  }

  return CUBIT_SUCCESS;
}


CubitStatus
FacetQueryEngine::restore_topology( FILE *file_ptr,
                                    unsigned int endian,
                                    int num_points,
                                    CubitPoint **points_array,
                                    int num_cfet,
                                    CurveFacetEvalTool** cfev_array,
                                    int num_fet,
                                    FacetEvalTool** fev_array,
                                    DLIList<TopologyBridge*> &imported_entities )
{

  //get file pointer
  unsigned int size;
  unsigned int i, j;
  int id, k;

  FacetPoint **facet_points;
  FacetCurve **facet_curves;
  FacetCoEdge **facet_coedges;
  FacetLoop **facet_loops;
  FacetSurface **facet_surfaces;
  FacetShell **facet_shells;
  FacetLump **facet_lumps;
  FacetBody **facet_bodies;

  FacetPoint *tmp_point;
  FacetCurve *tmp_curve;
  FacetCoEdge *tmp_coedge;
  FacetLoop *tmp_loop;
  FacetSurface *tmp_surf;
  FacetShell *tmp_shell;
  FacetLump *tmp_lump;
  FacetBody*tmp_body;

  int num_facet_points = 0;
  int num_facet_curves = 0;
  int num_facet_coedges = 0;
  int num_facet_loops = 0;
  int num_facet_surfaces = 0;
  int num_facet_shells = 0;
  int num_facet_lumps = 0;
  int num_facet_bodies = 0;

  //create a wrapper object for writing
  CIOWrapper file_reader( endian, file_ptr );

  //-----------------Read FacetPoints------------------
  // Read number of FacetPoints
  file_reader.Read( &size, 1 );
  // Allocate memory for FacetPoints
  facet_points = new FacetPoint*[size];
  for(i=0; i<size; i++)
  {
    file_reader.Read( reinterpret_cast<UnsignedInt32*>(&id), 1 );
    if( id >= num_points || id < 0 )
    {
      delete [] facet_points;
      return CUBIT_FAILURE;
    }

    tmp_point = new FacetPoint( points_array[ id ] );
    if( tmp_point == NULL )
      return CUBIT_FAILURE;

    num_facet_points++;
    facet_points[i] = tmp_point;
    tmp_point->restore_attribs( file_ptr, endian );
  }

  //-----------------Read FacetCurves------------------
  //Read number of FacetCurves
  file_reader.Read( &size, 1 );
  // Allocate memory for FacetCurves
  facet_curves = new FacetCurve*[size];
  for(i=0; i<size; i++)
  {
    int data[4];
    file_reader.Read( reinterpret_cast<UnsignedInt32*>(data), 4 );
    if( data[0] >= num_facet_points || data[0] < 0 ||
        data[1] >= num_facet_points || data[1] < 0 ||
        data[3] >= num_cfet ||  data[3] < 0 )
    {
      delete [] facet_curves;
      return CUBIT_FAILURE;
    }

    CubitSense sense;
    if (data[2] == -1 )
      sense = CUBIT_UNKNOWN;
    else
      sense = data[2] ? CUBIT_REVERSED : CUBIT_FORWARD;

    tmp_curve = new FacetCurve( cfev_array[ data[3] ], facet_points[ data[0] ],
                                facet_points[ data[1] ], sense );
    if( tmp_curve == NULL )
      return CUBIT_FAILURE;

    // Add curve to FacetPoints
    facet_points[data[0]]->add_curve( tmp_curve );
    facet_points[data[1]]->add_curve( tmp_curve );

    num_facet_curves++;
    facet_curves[i] = tmp_curve;
    tmp_curve->restore_attribs( file_ptr, endian );
  }

  //-----------------Read FacetCoEdges------------------
  //Read number of FacetCoEdges
  file_reader.Read( &size, 1 );
  facet_coedges = new FacetCoEdge*[ size ];
  for(i=0; i<size; i++)
  {
    int data[2];
    file_reader.Read( reinterpret_cast<unsigned int*>(data), 2 );

    if( data[0] >= num_facet_curves || data[0] < 0 )
    {
      delete [] facet_coedges;
      return CUBIT_FAILURE;
    }

    CubitSense sense;
    if (data[1] == -1 )
      sense = CUBIT_UNKNOWN;
    else
      sense = data[1] ? CUBIT_REVERSED : CUBIT_FORWARD;

    tmp_coedge =  new FacetCoEdge( facet_curves[ data[0] ], sense );
    if( tmp_coedge == NULL )
      return CUBIT_FAILURE;

    facet_coedges[i] = tmp_coedge;
    num_facet_coedges++;

    // Add FacetCoEdge to FacetCurve
    facet_curves[ data[0] ]->add_coedge( tmp_coedge );
  }

  //-----------------Read FacetLoops------------------
  //Read number of FacetLoops
  file_reader.Read( &size, 1 );
  DLIList<CoEdgeSM*> temp_coedge_list;
  facet_loops = new FacetLoop*[ size ];
  for(i=0; i<size; i++)
  {
    temp_coedge_list.clean_out();

    unsigned int num_coedges = 0;
    file_reader.Read( &num_coedges, 1);

    int* coedge_ids = new int[ num_coedges];
    file_reader.Read( (unsigned*)coedge_ids, num_coedges );

    for( j=0; j<num_coedges; j++)
    {
      if( coedge_ids[j] >= num_facet_coedges || coedge_ids[j] < 0 )
      {
        delete [] facet_loops;
        return CUBIT_FAILURE;
      }
      temp_coedge_list.append( facet_coedges[ coedge_ids[j] ] );
    }

    tmp_loop = new FacetLoop( temp_coedge_list );
    if( tmp_loop == NULL)
      return CUBIT_FAILURE;

    num_facet_loops++;
    facet_loops[i] = tmp_loop;

    for( j=0; j<num_coedges; j++)
      facet_coedges[ coedge_ids[j] ]->add_loop( tmp_loop );

    delete [] coedge_ids;
  }

  //-----------------Read FacetSurfaces------------------
  //Read number of FacetSurfaces
  file_reader.Read( &size, 1 );
  facet_surfaces = new FacetSurface*[size];
  DLIList<LoopSM*> temp_loops;
  for(i=0; i<size; i++)
  {
    temp_loops.clean_out();
    int data[6];
    file_reader.Read( (unsigned int*)data, 6);

    CubitSense sense, sense0, sense1;
    if (data[0] == -1 )
      sense = CUBIT_UNKNOWN;
    else
      sense = data[0] ? CUBIT_REVERSED : CUBIT_FORWARD;

    CubitBoolean useFacets;
    useFacets = data[1] ? CUBIT_TRUE : CUBIT_FALSE;

    // make sure FacetEvalTool ID is in range
    if( data[2] >= num_fet || data[2] < 0 )
    {
      delete [] facet_surfaces;
      return CUBIT_FAILURE;
    }

    if (data[3] == -1 )
      sense0 = CUBIT_UNKNOWN;
    else
      sense0 = data[3] ? CUBIT_REVERSED : CUBIT_FORWARD;

    if (data[4] == -1 )
      sense1 = CUBIT_UNKNOWN;
    else
      sense1 = data[4] ? CUBIT_REVERSED : CUBIT_FORWARD;

    int num_loops = data[5];
    int* loop_ids = new int[ num_loops ];
    file_reader.Read( (unsigned*)loop_ids, num_loops );

    for( k=0; k<num_loops; k++)
    {
      if( loop_ids[k] >= num_facet_loops ||
          loop_ids[k] < 0 )
      {
        delete [] loop_ids;
        return CUBIT_FAILURE;
      }
      temp_loops.append( facet_loops[ loop_ids[k] ] );
    }
    tmp_surf = new FacetSurface( fev_array[ data[2] ],
                                 sense, sense0,
                                 useFacets, temp_loops );
    if( tmp_surf == NULL)
      return CUBIT_FAILURE;

    facet_surfaces[i] = tmp_surf;
    num_facet_surfaces++;

    // Add FacetSurface to FacetLoops
    for( k=0; k<num_loops; k++)
      facet_loops[ loop_ids[k] ]->add_surface( tmp_surf );

    delete [] loop_ids;
    tmp_surf->restore_attribs( file_ptr, endian );
  }

  //-----------------Read FacetShells------------------
  //Read number of FacetShells
  file_reader.Read( &size, 1 );
  facet_shells = new FacetShell*[size];
  DLIList<Surface*> temp_surfs;
  for(i=0; i<size; i++)
  {
    temp_surfs.clean_out();
    unsigned int num_surfs = 0;
    file_reader.Read( &num_surfs, 1 );

    int* surface_ids = new int[ num_surfs ];
    file_reader.Read( (unsigned*)surface_ids, num_surfs );

    for( j=0; j<num_surfs; j++)
    {
      if( surface_ids[j] >= num_facet_surfaces ||
          surface_ids[j] < 0 )
      {
        delete [] facet_shells;
        return CUBIT_FAILURE;
      }

      temp_surfs.append( facet_surfaces[ surface_ids[j] ] );
    }

    tmp_shell = new FacetShell( temp_surfs );
    facet_shells[i] = tmp_shell;
    num_facet_shells++;

    // Add this shell to surfaces
    for( j=0; j<num_surfs; j++)
      facet_surfaces[ surface_ids[j] ]->add_shell( tmp_shell );

    delete [] surface_ids;
  }

  //-----------------Read FacetLumps------------------
  //Read number of FacetLumps
  file_reader.Read( &size, 1 );
  facet_lumps = new FacetLump*[size];
  DLIList<ShellSM*> temp_shells;
  for(i=0; i<size; i++)
  {
    temp_shells.clean_out();

    unsigned int num_shells = 0;
    file_reader.Read( &num_shells, 1 );

    int* shell_ids = new int[ num_shells ];
    file_reader.Read( (unsigned*)shell_ids, num_shells );

    for( j=0; j<num_shells; j++)
    {
      if( shell_ids[j] >= num_facet_shells )
      {
        delete [] facet_lumps;
        return CUBIT_FAILURE;
      }
      temp_shells.append( facet_shells[ shell_ids[j] ] );
    }

    tmp_lump = new FacetLump( temp_shells );
    if( tmp_lump == NULL )
      return CUBIT_FAILURE;

    facet_lumps[i] = tmp_lump;
    num_facet_lumps++;

    for( j=0; j<num_shells; j++)
      facet_shells[ shell_ids[j] ]->add_lump( tmp_lump );

    delete [] shell_ids;
    tmp_lump->restore_attribs( file_ptr, endian );
  }

  //-----------------Read FacetBodies ------------------
  //Read number of FacetBodies
  file_reader.Read( &size, 1 );
  facet_bodies = new FacetBody*[size];
  DLIList<Lump*> temp_lumps;
  for(i=0; i<size; i++)
  {
    temp_lumps.clean_out();

    unsigned int num_lumps= 0;
    file_reader.Read( &num_lumps, 1 );

    int* lump_ids = new int[ num_lumps ];
    file_reader.Read( (unsigned*)lump_ids, num_lumps );

    for( j=0; j<num_lumps; j++)
    {
      if( lump_ids[j] >= num_facet_lumps )
      {
        delete [] facet_bodies;
        return CUBIT_FAILURE;
      }
      temp_lumps.append( facet_lumps[ lump_ids[j] ] );
    }

    tmp_body = new FacetBody( temp_lumps );
    if( tmp_body == NULL )
      return CUBIT_FAILURE;

    facet_bodies[i] = tmp_body;
    num_facet_bodies++;

    // Add this FacetBody to FacetLumps
    for( j=0; j<num_lumps; j++)
      facet_lumps[ lump_ids[j] ]->add_body( tmp_body );

    delete [] lump_ids;

    //read in trans matrix
    unsigned int rows_and_cols[2];
    file_reader.Read( rows_and_cols, 2 );

    unsigned int num_rows = rows_and_cols[0];
    unsigned int num_cols = rows_and_cols[1];

    if( num_rows || num_cols)
    {
      CubitTransformMatrix trans_matrix;

      double *trans_matrix_array;
      trans_matrix_array = new double[ num_rows*num_cols ];
      file_reader.Read( trans_matrix_array, num_rows*num_cols );

      unsigned c;
      for(j=0; j<num_rows; j++ )
      {
        for(c=0; c<num_cols; c++)
          trans_matrix.add(j,c, trans_matrix_array[(j*num_cols)+c] );
      }
      tmp_body->set_transforms( trans_matrix );
      delete [] trans_matrix_array;
    }
    tmp_body->restore_attribs( file_ptr, endian );
  }

  // Here is where we determine if the entites are free or not
  // bodies, all bodies are free
  for(k=0; k<num_facet_bodies; k++)
    imported_entities.append( facet_bodies[k] );

  // surfaces are free if they are not in a shell
  DLIList<FacetShell*> shell_list;
  for(k=0; k<num_facet_surfaces; k++)
  {
    shell_list.clean_out();
    facet_surfaces[k]->get_shells( shell_list );
    if( shell_list.size() == 0 )
      imported_entities.append( facet_surfaces[k] );
  }

  // curves are free if they are not associate with a coedge
  DLIList<FacetCoEdge*> coedge_list;
  for(k=0; k<num_facet_curves; k++)
  {
    coedge_list.clean_out();
    facet_curves[k]->get_coedges( coedge_list );
    if( coedge_list.size() == 0 )
      imported_entities.append( facet_curves[k] );
  }

  // points are free if they are not associate with a curve
  DLIList<FacetCurve*> curve_list;
  for(k=0; k<num_facet_points; k++)
  {
    curve_list.clean_out();
    facet_points[k]->get_curves( curve_list );
    if( curve_list.size() == 0 )
      imported_entities.append( facet_points[k] );
  }

  // clean up
  delete [] facet_points;
  delete [] facet_curves;
  delete [] facet_coedges;
  delete [] facet_loops;
  delete [] facet_surfaces;
  delete [] facet_shells;
  delete [] facet_lumps;
  delete [] facet_bodies;

  return CUBIT_SUCCESS;
}


CubitStatus
FacetQueryEngine::import_temp_geom_file(FILE* file_ptr,
                                        const char* /*file_name*/,
                                        const char* file_type,
                                        DLIList<TopologyBridge*> &bridge_list )
{
  //make sure that file_type == "FACET"
  if( !strcmp( file_type,"FACET") )
    return import_solid_model( file_ptr, file_type, bridge_list );
  else
    return CUBIT_FAILURE;
}

CubitStatus FacetQueryEngine::import_solid_model(
                                                 const char* file_name ,
                                                 const char* file_type,
                                                 DLIList<TopologyBridge*> &imported_entities,
                                                 CubitBoolean print_results,
                                                 const char* logfile_name,
                                                 CubitBoolean heal_step,
                                                 CubitBoolean import_bodies,
                                                 CubitBoolean import_surfaces,
                                                 CubitBoolean import_curves,
                                                 CubitBoolean import_vertices,
                                                 CubitBoolean free_surfaces)
{
  errno = 0;
  FILE *file_ptr = fopen(file_name, "rb");
  if (!file_ptr)
  {
    PRINT_ERROR("Cannot open file: %s (%s)\n", file_name, strerror(errno) );
    return CUBIT_FAILURE;
  }

  CubitStatus status = import_solid_model(file_ptr, file_type, imported_entities, print_results,
                            logfile_name, heal_step, import_bodies, import_surfaces,
                            import_curves, import_vertices, free_surfaces);
  fclose(file_ptr);
  return status;
}

CubitStatus FacetQueryEngine::import_solid_model(FILE *file_ptr,
                                                 const char* /*file_type*/,
                                                 DLIList<TopologyBridge*> &imported_entities,
                                                 CubitBoolean ,
                                                 const char* ,
                                                 CubitBoolean,
                                                 CubitBoolean,
                                                 CubitBoolean,
                                                 CubitBoolean,
                                                 CubitBoolean,
                                                 CubitBoolean )

{
  CubitPoint **points_array = NULL;
  CurveFacetEvalTool **cfet_array = NULL;
  FacetEvalTool **fet_array = NULL;

  int num_points, num_edges, num_facets;
  int num_cfet, num_fet;

  // read in the file type "MESHED_BASED_GEOMETRY"
  char fileType[19] = {0};

  if( fread( fileType, 1, 19, file_ptr) != 19 )
  {
    PRINT_ERROR("Trouble reading in file type for MBG\n");
    return CUBIT_FAILURE;
  }
  
  if( strncmp( fileType, "MESH_BASED_GEOMETRY", 19 ) )
  {
    PRINT_ERROR("Not MESH_BASED_GEOMETRY file type\n");
    return CUBIT_FAILURE;
  }
  
  // read in the endian value
  NCubitFile::CIOWrapper file_reader(file_ptr, 19, 0);

  // read in version #
  UnsignedInt32 version;
  file_reader.Read( &version, 1 );

  //Read in points/edges/facets
  CubitStatus status;
  status = restore_facets( file_ptr, file_reader.get_endian(),
                           num_points, num_edges,
                           num_facets, points_array, num_cfet,
                           num_fet, cfet_array, fet_array );
  if( status == CUBIT_FAILURE)
  {
    PRINT_ERROR("Problems restore facets\n");
    return CUBIT_FAILURE;
  }

  //Restore Topology
  status = restore_topology( file_ptr, file_reader.get_endian(),
                             num_points, points_array,
                             num_cfet, cfet_array, num_fet,
                             fet_array, imported_entities);
  if( status == CUBIT_FAILURE)
  {
    PRINT_ERROR("Problems restore MDB topology\n");
    return CUBIT_FAILURE;
  }


  if(cfet_array != NULL)
    delete [] cfet_array;
  if(fet_array != NULL)
    delete [] fet_array;
  if(points_array != NULL)
    delete [] points_array;

  return CUBIT_SUCCESS;
}

CubitStatus FacetQueryEngine::import_solid_model(DLIList<TopologyBridge*> &imported_entities,
						 const char* pBuffer,
						 const int n_buffer_size)
{
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : restore_facets
// Member Type: PUBLIC
// Description: restore facets and eval tools onto the list of entities
// Author     : sjowen
// Date       : 1/26/03
//===============================================================================
CubitStatus FacetQueryEngine::restore_facets(
  FILE *fp,  // CUB file we are currently reading
  unsigned int endian,
  int &num_points,
  int &num_edges,
  int &num_facets,
  CubitPoint **&points,
  int &num_cfet,
  int &num_fet,
  CurveFacetEvalTool **&cfet_array,
  FacetEvalTool **&fet_array)
{
  CubitStatus rv = CUBIT_SUCCESS;

  CubitFacet **facets = NULL;
  CubitFacetEdge **edges = NULL;

  // read facets from the file and build arrays of facet entities

  rv = read_facets( fp, endian,
                    num_facets, num_edges, num_points,
                    facets,     edges,   points );

  // create the CurveFacetEvalTools and FacetEval Tools

  if (rv != CUBIT_FAILURE )
  {
    rv = restore_eval_tools( fp, endian,
                             num_facets, num_edges, num_points,
                             facets,     edges,  points,
                             num_cfet, num_fet,
                             cfet_array, fet_array );
  }

  if (facets != NULL)
    delete [] facets;
  if(edges != NULL)
    delete [] edges;
  return rv;

}


//===============================================================================
// Function   : read_facets
// Member Type: PUBLIC
// Description: read facets from the file and create facet entities
// Author     : sjowen
// Date       : 1/26/03
//===============================================================================
CubitStatus FacetQueryEngine::read_facets(
  FILE *fp,
  unsigned int endian,
  int &num_facets,
  int &num_edges,
  int &num_points,
  CubitFacet **&facets,
  CubitFacetEdge **&edges,
  CubitPoint **&points )
{

  NCubitFile::CIOWrapper cio(endian, fp );
  int ii;

  // read points
  UnsignedInt32 npoints;
  double uu,vv,ss;
  cio.Read(&npoints, 1);
  num_points = (int)npoints;
  if (num_points > 0)
  {
    double* coord_array  = new double [num_points * 3];
    double* uvs_array = new double [num_points * 3];
    cio.Read(coord_array, npoints*3);
    cio.Read(uvs_array, npoints*3);

    // create CubitPoints
    CubitPoint *point_ptr;
    //CubitVector normal;
    points = new CubitPoint * [num_points];
    for(ii=0; ii<num_points; ii++)
    {
      point_ptr = (CubitPoint *) new CubitPointData( coord_array[ii*3],
                                                     coord_array[ii*3+1],
                                                     coord_array[ii*3+2] );
      points[ii] = point_ptr;
      uu = uvs_array[ii*3];
      vv = uvs_array[ii*3+1];
      ss = uvs_array[ii*3+2];
      point_ptr->set_uvs(uu, vv, ss);
    }

    //Clean up
    delete [] coord_array;
    delete [] uvs_array;
  }


  //Read in Normals
  int nnormals;
  cio.Read(reinterpret_cast<UnsignedInt32*>(&nnormals), 1);
  if( nnormals > 0 )
  {
    double* normal_array = new double [nnormals * 3];
    int* normal_ids = new int[nnormals * 3];
    cio.Read(normal_array, nnormals * 3);
    cio.Read(reinterpret_cast<UnsignedInt32*>(normal_ids), nnormals);
    CubitVector normal;
    for(ii=0; ii<nnormals; ii++)
    {
      normal.x( normal_array[ii*3] );
      normal.y( normal_array[ii*3+1] );
      normal.z( normal_array[ii*3+2] );
      points[ normal_ids[ii] ]->normal( normal );
    }

    //Clean up
    delete [] normal_array;
    delete [] normal_ids;
  }


  // read edges and edge control points

  UnsignedInt32 nedges;
  cio.Read(&nedges, 1);
  num_edges = (int)nedges;
  if (num_edges > 0)
  {
    UnsignedInt32 *edge_vert_array = new UnsignedInt32 [nedges*2];
    cio.Read(edge_vert_array, nedges*2);
    UnsignedInt32 nctrl_pts;
    cio.Read(&nctrl_pts, 1);
    double *control_points = NULL;
    if (nctrl_pts > 0)
    {
      control_points = new double [nctrl_pts*3];
      cio.Read(control_points, nctrl_pts*3);
    }

    unsigned id0, id1, ii;
    edges = new CubitFacetEdge * [num_edges];
    CubitFacetEdge *edge_ptr;
    for(ii=0; ii<nedges; ii++)
    {
      id0 = edge_vert_array[ii*2];
      id1 = edge_vert_array[ii*2+1];

      edge_ptr = (CubitFacetEdge *) new CubitFacetEdgeData( points[id0], points[id1] );
      edges[ii] = edge_ptr;
      if (nctrl_pts > 0)
      {
        edge_ptr->set_control_points(&control_points[ii*NUM_EDGE_CPTS*3]);
      }
    }

    //Clean up
    delete [] edge_vert_array;
    delete [] control_points;
    edge_vert_array = NULL;
    control_points = NULL;
  }

  // read the facets and the facet control points

  UnsignedInt32 nfacets;
  cio.Read(&nfacets, 1);
  num_facets = (int)nfacets;
  if(num_facets > 0)
  {
    UnsignedInt32 *facet_edge_array = new UnsignedInt32 [nfacets*3];
    cio.Read(facet_edge_array, nfacets*3);
    int *int_data = new int [nfacets*2];
    cio.Read(reinterpret_cast<UnsignedInt32*>(int_data), nfacets*2);

    UnsignedInt32 nctrl_pts;
    cio.Read(&nctrl_pts, 1);
    double *control_points = NULL;
    if (nctrl_pts > 0)
    {
      control_points = new double [nctrl_pts*3];
      cio.Read(control_points, nctrl_pts*3);
    }

    unsigned id0, id1, id2, ii;
    CubitFacet *facet_ptr;
    facets = new CubitFacet * [num_facets];
    for (ii=0; ii<nfacets; ii++)
    {
      id0 = facet_edge_array[ii*3];
      id1 = facet_edge_array[ii*3+1];
      id2 = facet_edge_array[ii*3+2];
      facet_ptr = (CubitFacet *) new CubitFacetData(edges[id0], edges[id1], edges[id2]);
      facets[ii] = facet_ptr;
      facet_ptr->is_flat( int_data[ii*2] );
      facet_ptr->is_backwards( int_data[ii*2+1] );

      if(nctrl_pts > 0)
      {
        facet_ptr->set_control_points(&control_points[ii*NUM_TRI_CPTS*3]);
      }
    }

    //Clean up
    delete [] facet_edge_array;
    delete [] control_points;
    delete [] int_data;
    facet_edge_array = NULL;
    control_points = NULL;
    int_data = NULL;
  }

  // read the extra info at the surface boundaries

  UnsignedInt32 num_c_zero_points = 0;
  cio.Read(&num_c_zero_points, 1);
  if (num_c_zero_points > 0)
  {
    UnsignedInt32 c_zero_int_data_size;

    cio.Read(&c_zero_int_data_size, 1);
    if (c_zero_int_data_size <= 0)
      return CUBIT_FAILURE;
    UnsignedInt32 *c_zero_int32_data = new UnsignedInt32 [c_zero_int_data_size];
    cio.Read(c_zero_int32_data, c_zero_int_data_size);
    int *c_zero_int_data = new int [c_zero_int_data_size];
    unsigned int jj;
    for (jj=0; jj<c_zero_int_data_size; jj++)
      c_zero_int_data[jj] = (int) c_zero_int32_data[jj];

    UnsignedInt32 c_zero_double_data_size;

    cio.Read(&c_zero_double_data_size, 1);
    if (c_zero_double_data_size <= 0)
      return CUBIT_FAILURE;
    double *c_zero_double_data = new double [c_zero_double_data_size];
    cio.Read(c_zero_double_data, c_zero_double_data_size);

    // create the facet boundary tool datas and assign to points

    int didx = 0;
    int iidx = 0;
    UnsignedInt32 zz;
    for (zz=0; zz<num_c_zero_points; zz++)
    {
      if (didx >= (int)c_zero_double_data_size ||
          iidx >= (int)c_zero_int_data_size)
        return CUBIT_FAILURE;
      TDFacetBoundaryPoint::new_facet_boundary_point( points, facets,
        iidx, didx, c_zero_int_data, c_zero_double_data );
    }

    //Clean up
    delete [] c_zero_int_data;
    delete [] c_zero_int32_data;
    delete [] c_zero_double_data;
  }

  return CUBIT_SUCCESS;
}




//-------------------------------------------------------------------------
// Purpose       : Deletes all solid model entities associated with the
//                 Bodies in the input list.
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 4/23/01
//-------------------------------------------------------------------------
void FacetQueryEngine::delete_solid_model_entities(DLIList<BodySM*>&BodyList) const
{
  BodySM* BodyPtr = NULL;
  for (int i = 0; i < BodyList.size(); i++ )
  {
    BodyPtr = BodyList.get_and_step();
    this->delete_solid_model_entities(BodyPtr);
  }

  return;
}

//-------------------------------------------------------------------------
// Purpose       : Delete a FacetBody and child entities.
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/29/03
//-------------------------------------------------------------------------
CubitStatus
FacetQueryEngine::delete_solid_model_entities( BodySM* bodysm ) const
{
  FacetBody* fbody = dynamic_cast<FacetBody*>(bodysm);
  if (!fbody)
    return CUBIT_FAILURE;

  DLIList<FacetLump*> lumps;
  DLIList<FacetShell*> shells;
  DLIList<FacetSurface*> surfaces;

  fbody->get_lumps(lumps);
  fbody->disconnect_all_lumps();
  delete fbody;

  for (int i = lumps.size(); i--; )
  {
    FacetLump* lump = lumps.get_and_step();

    shells.clean_out();
    lump->get_shells(shells);
    lump->disconnect_all_shells();
    delete lump;

    for (int j = shells.size(); j--; )
    {
      FacetShell* shell = shells.get_and_step();

      surfaces.clean_out();
      shell->get_surfaces(surfaces);
      shell->disconnect_all_surfaces();
      delete shell;

      for (int k = surfaces.size(); k--; )
      {
        FacetSurface* surface = surfaces.get_and_step();
        if (!surface->has_parent_shell())
          delete_solid_model_entities(surface);
      }
    }
  }

  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Delete a FacetSurface and child entities.
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/29/03
//-------------------------------------------------------------------------
CubitStatus
FacetQueryEngine::delete_solid_model_entities( Surface* surface ) const
{
  FacetSurface* fsurf = dynamic_cast<FacetSurface*>(surface);
  if (!fsurf || fsurf->has_parent_shell())
    return CUBIT_FAILURE;

  DLIList<FacetLoop*> loops;
  DLIList<FacetCoEdge*> coedges;

  fsurf->get_loops(loops);
  fsurf->disconnect_all_loops();
  delete fsurf;

  for (int i = loops.size(); i--; )
  {
    FacetLoop* loop = loops.get_and_step();

    coedges.clean_out();
    loop->get_coedges(coedges);
    loop->disconnect_all_coedges();
    delete loop;

    for (int j = coedges.size(); j--; )
    {
      FacetCoEdge* coedge = coedges.get_and_step();
      FacetCurve* curve = dynamic_cast<FacetCurve*>(coedge->curve());
      if (curve)
      {
        curve->disconnect_coedge(coedge);
        if (!curve->has_parent_coedge())
          delete_solid_model_entities(curve);
      }

      delete coedge;
    }
  }

  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Delete a FacetCurve and child entities.
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/29/03
//-------------------------------------------------------------------------
CubitStatus
FacetQueryEngine::delete_solid_model_entities( Curve* curve ) const
{
  FacetCurve* fcurve = dynamic_cast<FacetCurve*>(curve);
  if (!fcurve || fcurve->has_parent_coedge())
    return CUBIT_FAILURE;

  FacetPoint* start = dynamic_cast<FacetPoint*>(fcurve->start_point());
  FacetPoint*   end = dynamic_cast<FacetPoint*>(fcurve->end_point()  );

  if (start == end )
      end = NULL;

  if (start)
  {
    start->disconnect_curve(fcurve);
    if (!start->has_parent_curve())
      delete_solid_model_entities(start);
  }

  if (end)
  {
    end->disconnect_curve(fcurve);
    if (!end->has_parent_curve())
      delete_solid_model_entities(end);
  }

  delete curve;
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Delete a FacetPoint and child entities.
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/29/03
//-------------------------------------------------------------------------
CubitStatus
FacetQueryEngine::delete_solid_model_entities( Point* point ) const
{
  FacetPoint* fpoint = dynamic_cast<FacetPoint*>(point);
  if (!fpoint || fpoint->has_parent_curve())
    return CUBIT_FAILURE;

  delete point;
  return CUBIT_SUCCESS;
}

CubitStatus FacetQueryEngine::fire_ray(BodySM *,
                                          const CubitVector &,
                                          const CubitVector &,
                                          DLIList<double>&,
                                          DLIList<GeometryEntity*> *) const
{
  PRINT_ERROR("FacetQueryEngine::fire_ray not yet implemented.\n");
  return CUBIT_FAILURE;
}
  //- fire a ray at the specified body, returning the entities hit and
  //- the parameters along the ray; return CUBIT_FAILURE if error

double FacetQueryEngine::get_sme_resabs_tolerance() const
{
  PRINT_ERROR("FacetQueryEngine::get_sme_resabs_tolerance not yet implemented.\n");
  return CUBIT_FAILURE;
}
// Gets solid modeler's resolution absolute tolerance

double FacetQueryEngine::set_sme_resabs_tolerance( double )
{
  PRINT_ERROR("FacetQueryEngine::set_sme_resabs_tolerance not yet implemented.\n");
  return CUBIT_FAILURE;
}

CubitStatus FacetQueryEngine::set_int_option( const char* , int )
{
  PRINT_ERROR("FacetQueryEngine::set_int_option not yet implemented.\n");
  return CUBIT_FAILURE;
}

CubitStatus FacetQueryEngine::set_dbl_option( const char* , double )
{
  PRINT_ERROR("FacetQueryEngine::set_dbl_option not yet implemented.\n");
  return CUBIT_FAILURE;
}

CubitStatus FacetQueryEngine::set_str_option( const char* , const char* )
{
  PRINT_ERROR("FacetQueryEngine::set_str_option not yet implemented.\n");
  return CUBIT_FAILURE;
}
  //- Set solid modeler options

//===========================================================================
//Function Name: make_facets
//Member Type:  PUBLIC
//Description:  creates Cubit quad facet entities from the hash points and
//              connectivity
//===========================================================================
CubitStatus FacetQueryEngine::make_facets(
  int *conn,        // conectivity array (size = 4 * nfacets)
  int nfacets,      // total number of facets (tri+quad)
  DLIList<CubitQuadFacet *> &facet_list ) // return the facet list
{
  CubitQuadFacet *facet_ptr = NULL;
  CubitStatus rv = CUBIT_SUCCESS;
  CubitPoint *point0, *point1, *point2, *point3;

  // create the facet array

  for(int ii=0; ii<nfacets; ii++)
  {
    if (conn[ii*4+2] != conn[ii*4+3])
    {
      point0 = get_hash_point(conn[ii*4]);
      point1 = get_hash_point(conn[ii*4+1]);
      point2 = get_hash_point(conn[ii*4+2]);
      point3 = get_hash_point(conn[ii*4+3]);

      facet_ptr = new CubitQuadFacetData( point0, point1, point2, point3 );

      if (!facet_ptr)
      {
        rv = CUBIT_FAILURE;
        return rv;
      }
      facet_list.append( facet_ptr );
    }
  }
  return rv;
}

//===========================================================================
//Function Name: make_facets
//Member Type:  PUBLIC
//Description:  creates Cubit tri facet entities from the hash points and
//              connectivity
//===========================================================================
CubitStatus FacetQueryEngine::make_facets(
  int *conn,        // conectivity array (size = 4 * nfacets)
  int nfacets,      // total number of facets (tri+quad)
  DLIList<CubitFacet *> &facet_list ) // return the facet list
{
  CubitFacet *facet_ptr = NULL;
  CubitStatus rv = CUBIT_SUCCESS;
  CubitPoint *point0, *point1, *point2;

  // create the facet array

  for(int ii=0; ii<nfacets; ii++)
  {
    if (conn[ii*4+2] == conn[ii*4+3])
    {
      point0 = get_hash_point(conn[ii*4]);
      point1 = get_hash_point(conn[ii*4+1]);
      point2 = get_hash_point(conn[ii*4+2]);
      if( (point0 == point1) || (point0 == point2) || (point1 == point2) ){
          PRINT_ERROR("Point used more than once in a single facet.  This is not allowed.\n");
          return CUBIT_FAILURE;
      }
      if( !point0 || !point1 || !point2 ){
          PRINT_ERROR("Point in facet is undefined.  This is not allowed.\n");
          return CUBIT_FAILURE;
      }
      
      facet_ptr = new CubitFacetData( point0, point1, point2 );

      if (!facet_ptr)
      {
        rv = CUBIT_FAILURE;
        return rv;
      }
      facet_list.append( facet_ptr );
    }
  }

 return rv;
}

//===========================================================================
//Function Name: ensure_is_ascii_stl_file
//Member Type:
//Description: returns CUBIT_TRUE in is_ascii if fp is an ascii stl file
//Author: Plamen Stoyanov (USF)
//===========================================================================
CubitStatus FacetQueryEngine::ensure_is_ascii_stl_file(FILE * fp, CubitBoolean &is_ascii)
{

  char line[128]="";

  if (fgets(line, 128, fp)==NULL)
  {
    return CUBIT_FAILURE;
  }
  if (fgets(line, 128, fp)==NULL)
  {
    return CUBIT_FAILURE;
  }
  if (strlen(line)==127)
  {
    if (fgets(line, 128, fp)==NULL)
    {
      return CUBIT_FAILURE;
    }
  }


  unsigned int dummy_int=0;

  while (isspace(line[dummy_int])&& dummy_int<strlen(line)) dummy_int++;

  if (strlen(line)-dummy_int>5)
  {
    if (tolower(line[dummy_int++])=='f' &&
      tolower(line[dummy_int++])=='a' &&
      tolower(line[dummy_int++])=='c' &&
      tolower(line[dummy_int++])=='e' &&
      tolower(line[dummy_int])=='t')
    {
      if (fgets(line, 128, fp)==NULL)
      {
        return CUBIT_FAILURE;
      }
      dummy_int=0;
      while (isspace(line[dummy_int])&& dummy_int<strlen(line))
      {
        dummy_int++;
      }
      if (strlen(line)-dummy_int>5)
      {
        if (tolower(line[dummy_int++])=='o' &&
          tolower(line[dummy_int++])=='u' &&
          tolower(line[dummy_int++])=='t' &&
          tolower(line[dummy_int++])=='e' &&
          tolower(line[dummy_int])=='r')
        {
          if (fgets(line, 128, fp)==NULL)
          {
            return CUBIT_FAILURE;
          }
          dummy_int=0;
          while (isspace(line[dummy_int])&& dummy_int<strlen(line)) {
            dummy_int++;
          }
          if (strlen(line)-dummy_int>6)
          {
            if (tolower(line[dummy_int++])=='v' &&
              tolower(line[dummy_int++])=='e' &&
              tolower(line[dummy_int++])=='r' &&
              tolower(line[dummy_int++])=='t' &&
              tolower(line[dummy_int++])=='e'	&&
              tolower(line[dummy_int])=='x')
            {
              is_ascii=CUBIT_TRUE;
            }
          }
        }
      }
    }
  }
  return CUBIT_SUCCESS;
}

//===========================================================================
//Function Name: read_facets_stl
//Member Type:
//Description:  read facets from stl file combining vertices within tolerance
//distance
//Author: Plamen Stoyanov (USF)
//===========================================================================
CubitStatus FacetQueryEngine::read_facets_stl_tolerance(
                                              DLIList<CubitFacet *> &tfacet_list,
                                              DLIList<CubitPoint *> & /*point_list*/,
                                              const char * file_name,
                                              int &npoints,
                                              int &ntri,
                                              long& seek_address,
                                              double tolerance
                                              )
{
  
  FILE *fp = fopen(file_name, "r");
  if (fp == NULL)
  {
    PRINT_ERROR("Could not open file %s for reading\n", file_name);
    seek_address = 0;
    return CUBIT_FAILURE;
  }

  DLIList <CubitPoint *> file_points;
  CubitPoint *current_point;
  int points_in_file=0;
  int ii;

  CubitBoolean is_ascii=CUBIT_FALSE;
  if (!ensure_is_ascii_stl_file(fp, is_ascii))
  {
    seek_address = 0;
    fclose(fp);
    return CUBIT_FAILURE;
  }

  if (is_ascii==CUBIT_TRUE)
  {

    PRINT_INFO("Reading facets...\n");
    fclose(fp);
    fp = fopen(file_name, "r");
    fseek(fp,seek_address,SEEK_SET);
    
    char line[128], junk[30];
    int numverts;
    double xx[3], yy[3], zz[3];
    int linenumber, num;
    bool done, error_found, eof_found;

    linenumber = 0;

    strcpy(line,"");
    done = false;
    error_found = false;
    eof_found = false;
    while ( (strstr(line,"endsolid") == 0) && (error_found == false) ) {
      numverts = 0;    
      while ( numverts < 3 ) {
        if ( fgets(line,127,fp) == 0 ) {
          linenumber++;
          eof_found = true;
          break; // EOF
        }
        linenumber++;        
//      makelowercase(line);
        int len, ij;
        len = strlen(line);
        for ( ij = 0; ij < len; ij++ ) 
          line[ij] = tolower(line[ij]);

        if ( strstr(line,"endsolid") != 0 ) {
          done = true;
          break; // End of part definition
        }
        if ( strstr(line,"vertex") != 0 ) {
          num = sscanf(line,"%s %le %le %le",junk,&xx[numverts],&yy[numverts],&zz[numverts]);
          if ( num != 4 ) {
            error_found = true;
            break; // error in reading vertices
          }
          numverts += 1;
        }
      } // end of while ( numverts < 3 )  

      if ( (eof_found == true) || (done == true)  || 
           (error_found == true) || (numverts != 3) ) break;

      current_point = (CubitPoint *) new CubitPointData(xx[0],yy[0],zz[0]);
      current_point->set_id(points_in_file++);
      file_points.append(current_point);
      current_point = (CubitPoint *) new CubitPointData(xx[1],yy[1],zz[1]);
      current_point->set_id(points_in_file++);
      file_points.append(current_point);
      current_point = (CubitPoint *) new CubitPointData(xx[2],yy[2],zz[2]);
      current_point->set_id(points_in_file++);
      file_points.append(current_point);

    }
    seek_address = 0;
    if ( eof_found == true ) {
      PRINT_WARNING("Premature end-of-file on STL file. Body may be incomplete.\n");
      fclose(fp);
      goto end_read_file_points;
    }
    if ( error_found == true ) {
      PRINT_WARNING("Error found while reading line number %d of file %s. Body may be incomplete.\n",
              linenumber,file_name);
      fclose(fp);       
      goto end_read_file_points;
    }
    if ( done == false ) {
      PRINT_WARNING("Error found while reading line number %d of file %s. Body may be incomplete.\n",
              linenumber,file_name); 
      fclose(fp);         
      goto end_read_file_points;
    }  

    if ( (eof_found == false) && (error_found == false) && (done == true) ) {
      while ( fgets(line,127,fp) != 0 ) {
        if ( (strstr(line,"solid") != 0) && (strstr(line,"endsolid") == 0) ) {
          seek_address = ftell(fp);
          break;
        }
      }
    }
    fclose(fp);
    goto end_read_file_points;
  }
  else
  {
    fclose(fp);
    // file is closed so that it can be opened as binary
    fp = fopen(file_name, "rb");
    if (fp == NULL)
    {
      PRINT_ERROR("Could not open file %s for reading\n", file_name);
      return CUBIT_FAILURE;
    }

    char dummy;
        // iterates through the facets of the file
    float cur[12];		// an array to hold 48 bytes representing 1AVS facet

    fseek(fp, 80, SEEK_SET);
    fread(&ntri, 4, 1, fp);

    PRINT_INFO ("Reading facets...\n");
    for (ii=0; ii<ntri; ii++) {
      // read in 1 facet
      if (!fread(cur, 4, 12, fp))
      {
        PRINT_INFO ("Abnormal file termination %s \n", file_name); break;
      }
    	fread(&dummy, 1, 1, fp);
       fread(&dummy, 1, 1, fp);

        //make point
      current_point = (CubitPoint *) new CubitPointData(cur[3],cur[4],cur[5]);
      current_point->set_id(points_in_file++);
      file_points.append(current_point);
      current_point = (CubitPoint *) new CubitPointData(cur[6],cur[7],cur[8]);
      current_point->set_id(points_in_file++);
      file_points.append(current_point);
      current_point = (CubitPoint *) new CubitPointData(cur[9],cur[10],cur[11]);
      current_point->set_id(points_in_file++);
      file_points.append(current_point);
    }
  }
  fclose(fp);
  
// at this point all points from the file are in file_points
end_read_file_points:

  // grid search tree to hold points
  GridSearchTree * node_grid = new GridSearchTree (tolerance);

  CubitPoint
             *point0,
             *point1,
             *point2;
  CubitFacet
             *facet_ptr;
  ntri=0;
  npoints=0;

  if (file_points.size() % 3 != 0)
  {
    PRINT_INFO("File Error.");
    return CUBIT_FAILURE;
  }

  for (ii = file_points.size(); ii>0; ii-=3)
  {

    // get three points from the file_points list
    // and compare them against the data already in the
    // grid search tree; if it contains a point within
    // epsilon tolerance this point would replace the
    // current point
    point0 = node_grid->fix(file_points.get_and_step());
    point1 = node_grid->fix(file_points.get_and_step());
    point2 = node_grid->fix(file_points.get_and_step());


    if (point0 && point1 && point2 && point0 != point1 && point1!=point2 && point2!=point0 )
    {

      facet_ptr = new CubitFacetData(point0, point1, point2);

      if (facet_ptr)
      {
        tfacet_list.append(facet_ptr);
        ntri++;
      }

    }
  }

  //delete node_grid;
  return CUBIT_SUCCESS;
}


//===========================================================================
//Function Name: read_facets_stl
//Member Type:
//Description:  read facets from stl file
//Author: Plamen Stoyanov (USF)
//===========================================================================

CubitStatus FacetQueryEngine::read_facets_stl(
                                              DLIList<CubitFacet *> &tfacet_list,
                                              DLIList<CubitPoint *> &point_list,
                                              const char * file_name,
                                              int &npoints,
                                              int &ntri,
                                              long& seek_address
                                              )
{

  ntri = 0;
  npoints = 0;
  FILE *fp = fopen(file_name, "r");
  if (fp == NULL)
  {
    PRINT_ERROR("Could not open file %s for reading\n", file_name);
    seek_address = 0;
    return CUBIT_FAILURE;
  }

  typedef std::map< CubitPoint * , int, CubitPointComparator > vMap;

  vMap mm;			// binary search tree to hold the vertices for efficiency
  vMap::iterator pos;
  bool append_to_facet_list=CUBIT_TRUE;

  CubitBoolean is_ascii=CUBIT_FALSE;

  if (!ensure_is_ascii_stl_file(fp, is_ascii))
  {
    seek_address = 0;
    fclose(fp);
    return CUBIT_FAILURE;
  }

  if (is_ascii==CUBIT_TRUE)
  {
    CubitPoint *point0,*point1,*point2;
    CubitFacet *facet_ptr = NULL;

    PRINT_INFO("Reading facets...\n");
    fclose(fp);
    fp = fopen(file_name, "r");
    fseek(fp,seek_address,SEEK_SET);
    char line[128], junk[30];
    int numverts;
    double xx[3], yy[3], zz[3];
    int linenumber, num;
    bool done, error_found, eof_found;

    linenumber = 0;

    strcpy(line,"");
    done = false;
    error_found = false;
    eof_found = false;
    while ( (strstr(line,"endsolid") == 0) && (error_found == false) ) {
      numverts = 0;    
      while ( numverts < 3 ) {
        if ( fgets(line,127,fp) == 0 ) {
          linenumber++;
          eof_found = true;
          break; // EOF
        }
        linenumber++;        
//      makelowercase(line);
        int len, ij;
        len = strlen(line);
        for ( ij = 0; ij < len; ij++ ) 
          line[ij] = tolower(line[ij]);

        if ( strstr(line,"endsolid") != 0 ) {
          done = true;
          break; // End of part definition
        }
        if ( strstr(line,"vertex") != 0 ) {
          num = sscanf(line,"%s %le %le %le",junk,&xx[numverts],&yy[numverts],&zz[numverts]);
          if ( num != 4 ) {
            error_found = true;
            break; // error in reading vertices
          }
          numverts += 1;
        }
      } // end of while ( numverts < 3 )  

      if ( (eof_found == true) || (done == true)  || 
           (error_found == true) || (numverts != 3) ) break;
      point0 = (CubitPoint *) new CubitPointData(xx[0],yy[0],zz[0]);
      pos=mm.find(point0);

      if (pos==mm.end())
      {
        mm.insert ( vMap::value_type(point0, npoints));
        point0->set_id( npoints ++ );
        //point_list is output of functionAVS
        point_list.append(point0);
      }
      else
      {
        delete point0;
        point0=((*pos).first);
      }
      point1 = (CubitPoint *) new CubitPointData(xx[1],yy[1],zz[1]);
      pos=mm.find(point1);
      if (pos==mm.end())
      {
        mm.insert ( vMap::value_type(point1, npoints));
        point1->set_id( npoints ++ );
        //point_list is output of function
        point_list.append(point1);
      }
      else
      {
        delete point1;
        point1=((*pos).first);
      }
      point2 = (CubitPoint *) new CubitPointData(xx[2],yy[2],zz[2]);
      pos=mm.find(point2);
      if (pos==mm.end())
      {
        mm.insert ( vMap::value_type(point2, npoints));
        point2->set_id( npoints ++ );
        //point_list is output of function
        point_list.append(point2);
      }
      else
      {
        delete point2;
        point2=((*pos).first);
      }

      if (point0 && point1 && point2 && point0!=point1 && point1!=point2 && point2 !=point0)
      {
        facet_ptr = new CubitFacetData( point0, point1, point2 );
        append_to_facet_list=CUBIT_TRUE;
      }
      else
      {
        append_to_facet_list=CUBIT_FALSE;
      }
    
      if (!facet_ptr)
      {
        seek_address = 0;
        fclose(fp);
        return CUBIT_FAILURE;
      }

      /// APPEND facet to output facet_ptr list
      if (append_to_facet_list)
      {
        tfacet_list.append( facet_ptr );
        ntri++;
      }
      
    }
    seek_address = 0;
    if ( eof_found == true ) {
      PRINT_WARNING("Premature end-of-file on STL file. Body may be incomplete.\n");
      fclose(fp);
      return CUBIT_SUCCESS;
    }
    if ( error_found == true ) {
      PRINT_WARNING("Error found while reading line number %d of file %s. Body may be incomplete.\n",
              linenumber,file_name);
      fclose(fp);
      return CUBIT_SUCCESS;
    }
    if ( done == false ) {
      PRINT_WARNING("Error found while reading line number %d of file %s. Body may be incomplete.\n",
              linenumber,file_name); 
      fclose(fp); 
      return CUBIT_SUCCESS;
    }  

    if ( (eof_found == false) && (error_found == false) && (done == true) ) {
      while ( fgets(line,127,fp) != 0 ) {
        if ( (strstr(line,"solid") != 0) && (strstr(line,"endsolid") == 0) ) {
          seek_address = ftell(fp);
//          PRINT_INFO("This STL file is a multipart file.  Only the first part was read.\n");
          break;
        }
      }
    }
    fclose(fp);
      
    return CUBIT_SUCCESS;
  }
  else
  {
    fclose(fp);
    // file is closed so that it can be opened as binary
    fp = fopen(file_name, "rb");
    if (fp == NULL)
    {
      PRINT_ERROR("Could not open file %s for reading\n", file_name);
      return CUBIT_FAILURE;
    }


    CubitFacet *facet_ptr = NULL;
    CubitPoint *point0,
              *point1,
                *point2;

    char dummy;
    int		   ii=0;    // iterates through the facets of the file
    float cur[12];		// an array to hold 48 bytes representing 1 facet
    npoints=0;          // only tri facets in stl files

    fseek(fp, 80, SEEK_SET);
    fread(&ntri, 4, 1, fp);
    int distinct_ntri=ntri;

    PRINT_INFO ("Reading facets...\n");
    for (ii=0; ii<ntri; ii++) {
      // read in 1 facet
      if (!fread(cur, 4, 12, fp))
      {
        PRINT_INFO ("Abnormal file termination %s \n", file_name); break;
      }
    		fread(&dummy, 1, 1, fp);
        fread(&dummy, 1, 1, fp);

        //make point
        point0 = (CubitPoint *) new CubitPointData( cur[3], cur[4], cur[5] );

        pos=mm.find(point0);
        if (pos==mm.end())
        {
          mm.insert ( vMap::value_type(point0, npoints));
          point0->set_id( npoints ++ );
          //point_list is output of function
          point_list.append(point0);
        }
        else
        {
          delete point0;
          point0=((*pos).first);
        }

        //make point
        point1 = (CubitPoint *) new CubitPointData( cur[6], cur[7], cur[8] );

        pos=mm.find(point1);
        if (pos==mm.end())
        {
          mm.insert ( vMap::value_type(point1, npoints));
          point1->set_id( npoints ++ );
          //point_list is output of function
          point_list.append(point1);
        }
        else
        {
          delete point1;
          point1=((*pos).first);
        }

        //make point
        point2 = (CubitPoint *) new CubitPointData( cur[9], cur[10], cur[11] );

        pos=mm.find(point2);
        if (pos==mm.end())
        {
          mm.insert ( vMap::value_type(point2, npoints));
          point2->set_id( npoints ++ );
          //point_list is output of function
          point_list.append(point2);
        }
        else
        {
          delete point2;
          point2=((*pos).first);
        }

        // this is to avoid a facet with all points on same line, which crashes CUBIT
        // because of assertion
        if (point0 && point1 && point2 && point0!=point1 && point1!=point2 && point2 !=point0)
        {
          facet_ptr = new CubitFacetData( point0, point1, point2 );
          append_to_facet_list=CUBIT_TRUE;
        }
        else
        {
          append_to_facet_list=CUBIT_FALSE;
        }

        if (!facet_ptr)
        {
          fclose(fp);
          return CUBIT_FAILURE;
        }

        if (append_to_facet_list)
        {
          tfacet_list.append( facet_ptr );
        }
        else
        {
          distinct_ntri--;
        }
    }
    ntri=distinct_ntri;
    fclose(fp);
    return CUBIT_SUCCESS;
  }
}


//===========================================================================
//Function Name: import_facets
//Member Type:  PUBLIC
//Description:  import facets and create geometry
//===========================================================================
CubitStatus FacetQueryEngine::import_facets(
  const char *file_name,                // to be read
  CubitBoolean use_feature_angle,       // to define where surfaces are broken
  double feature_angle,        // angle where surfaces are broken (degrees)
  double tolerance,             // for stl files loading
  int interp_order,            // Facet representation 1= linear, 4= Bezier
  CubitBoolean smooth_non_manifold,  // continutiy accross non-manifold edges
  CubitBoolean split_surfaces,   // break facet rep at surface boundaries
  CubitBoolean stitch,             // attempt to stitch together facets at the boundary
  CubitBoolean improve,      // smooth, swap and collapse to improve quality
  DLIList<CubitQuadFacet *> &quad_facet_list,  // return quad facets
  DLIList<CubitFacet *> &tri_facet_list,  // return tri facets
  DLIList<Surface *> &surface_list,  // return list of surfaces
  FacetFileFormat file_format )
{
  DLIList <CubitFacet *>tfacet_list;
  DLIList <CubitQuadFacet *>qfacet_list;
  DLIList <CubitPoint *>point_list;

  ShellSM *shell_ptr;
  DLIList<ShellSM*> shell_list;
  Lump *lump_ptr;
  DLIList<Lump*> lump_list;
  BodySM *bodysm_ptr;
  Body *body_ptr;
  GeometryQueryTool *gti = GeometryQueryTool::instance();
  FacetShell* facet_shell;
  DLIList<DLIList<CubitFacet *> *> shell_facet_list;
  int ishell, ii;
  CubitBoolean is_water_tight = CUBIT_TRUE;

  // read the facets from a file

  int *conn = NULL;
  int npoints = 0;
  int nfacets = 0;
  int nquad = 0;
  int ntri = 0;

  int prev_vert;
  int prev_edge;
  int prev_face;
  int prev_vol;
  int prev_bod;

  long stl_seek_address=0;
  bool stl_multiple_parts = true;
  
  CubitStatus rv;
  //  Here we add the capability to read stl files with more than one part.
  while ( stl_multiple_parts == true ) {
    stl_multiple_parts = false;
    prev_vert = gti->num_ref_vertices();
    prev_edge = gti->num_ref_edges();
    prev_face = gti->num_ref_faces();
    prev_vol  = gti->num_ref_volumes();
    prev_bod  = gti->num_bodies();
    switch (file_format)
    {
    case STL_FILE:
      tfacet_list.clean_out(); // In case there are multiple parts in the stl file
      point_list.clean_out();
      shell_list.clean_out();
      lump_list.clean_out();
      shell_facet_list.clean_out();
      if (tolerance>0)
      {
        rv = read_facets_stl_tolerance(tfacet_list, point_list, file_name, npoints, ntri, stl_seek_address, tolerance);
      }
      else
      {
        rv = read_facets_stl(tfacet_list, point_list, file_name, npoints, ntri, stl_seek_address);
      }
      PRINT_INFO("  %d facets read.\n", ntri);
      if ( (rv == CUBIT_SUCCESS) && (stl_seek_address) > 0 ) stl_multiple_parts = true;
      nfacets = ntri;

      if (rv != CUBIT_SUCCESS)
      {
        goto end_import_facets;
      }
      break;

    case CUBIT_FACET_FILE:
    case AVS_FILE:
      rv = read_facets( file_name, conn, npoints, nquad, ntri, file_format );

      nfacets = ntri + nquad;

      if (rv != CUBIT_SUCCESS)
        goto end_import_facets;

      // make cubit facet entities from the points and connectivity

      if (nquad  > 0)
        rv = make_facets(conn, nfacets, qfacet_list);
      if (ntri > 0)
        rv = make_facets(conn, nfacets, tfacet_list);
      if (rv != CUBIT_SUCCESS)
        goto end_import_facets;
      get_all_hash_points(point_list);
      delete_hash_points();
      break;
    case CHOLLA_FILE:
      rv = read_cholla_file( file_name, feature_angle, point_list, tfacet_list );
      nquad = 0;
      ntri = tfacet_list.size();
      npoints = point_list.size();
      break;
    case FROM_FACET_LIST:
      tfacet_list = tri_facet_list;
      qfacet_list = quad_facet_list;
      if (tfacet_list.size() + qfacet_list.size() == 0)
      {
        PRINT_ERROR("No facets found to build geometry\n");
        rv = CUBIT_FAILURE;
        goto end_import_facets;
      }
      break;
    default:
      assert(0); // unrecognized file format
      break;
    }

    if (tfacet_list.size() + qfacet_list.size() == 0)
    {
      PRINT_ERROR("No facets read from file %s.\n", file_name);
      rv = CUBIT_FAILURE;
      goto end_import_facets;
    }
    else
    {
      PRINT_INFO("Building facet-based geometry from %d facets...\n",
                  tfacet_list.size() + qfacet_list.size() );
    }

//  if (fix)
//  {
//    rv = check_facets( point_list, tfacet_list );
//  }

    if (0)
    {
      // Call function to generate an x-y-z file for a CTH/SPH simulation
      // the following is an exampl call...
      //make_sph( point_list, tfacet_list, 100.0, "fem1-sph100.xyz" );
    }

    if (0)
    {
      // This is an example of using the export_facets function.  Writes
      // a facet file of all facets in the list

        //commented out because otherwise it's a compiler warning.

//     char filename[128];
//     strcpy(filename, "my_test.facets");
//     export_facets(tfacet_list,filename);
    }


    // split the facets into shells if needed

    rv = FacetDataUtil::split_into_shells(tfacet_list, qfacet_list,
                                          shell_facet_list, is_water_tight);
    if (rv != CUBIT_SUCCESS)
    {
      PRINT_ERROR("Error processing facets from %s.\n", file_name);
      goto end_import_facets;
    }

    // if the facets aren't watertight, see if they can be merged

    if (!is_water_tight && stitch)
    {
      rv = FacetDataUtil::stitch_facets(shell_facet_list,
                                        GEOMETRY_RESABS,
                                        is_water_tight);
      if (rv != CUBIT_SUCCESS)
      {
        PRINT_WARNING("Couldn't stitch facets.\n");
      }
    }

    DLIList <CubitFacet *> *facet_list_ptr;
    if (improve)
    {
      for (ishell = 0; ishell < shell_facet_list.size(); ishell++)
      {
        facet_list_ptr = shell_facet_list.get_and_step();
        rv = FacetDataUtil::collapse_short_edges( *facet_list_ptr,  CUBIT_TRUE );
        if (rv != CUBIT_SUCCESS)
        {
          PRINT_WARNING("Couldn't improve facets.\n");
        }
      }
    }

    // create the surface geometry

    if (!use_feature_angle)
      feature_angle = -1.0;

    for (ishell = 0; ishell < shell_facet_list.size(); ishell++)
    {
      DLIList <Surface *> shell_surfaces;
      DLIList <CubitPoint *> mypoint_list;
      facet_list_ptr = shell_facet_list.get_and_step();
      rv = FacetModifyEngine::instance()->build_facet_surface( NULL,
                               *facet_list_ptr, mypoint_list,
                               feature_angle, interp_order,
                               smooth_non_manifold,
                               split_surfaces, shell_surfaces);
      if (rv != CUBIT_SUCCESS || shell_surfaces.size() == 0)
      {
        PRINT_ERROR("Couldn't build facet based geometry from facets in %s\n", file_name);
        rv = CUBIT_FAILURE;
        goto end_import_facets;
      }

      // make a shell out of these surfaces

      rv = FacetModifyEngine::instance()->make_facet_shell(shell_surfaces, shell_ptr);
      if ( shell_ptr == NULL || rv != CUBIT_SUCCESS )
      {
        PRINT_ERROR("Problems building facet based shell entity.\n");
        rv = CUBIT_FAILURE;
        goto end_import_facets;
      }

      //Set the sense for the surfaces (will be cofaces) on this shell.
      //Assumption: The sense is always forward when creating geom from facets.
      // (This may not be correct -especially with multiple shells in a body)

      facet_shell = CAST_TO( shell_ptr, FacetShell );
      for( ii = shell_surfaces.size(); ii > 0; ii-- )
      {
        Surface* surf = shell_surfaces.get_and_step();
        FacetSurface* facet_surf = CAST_TO( surf, FacetSurface );
        facet_surf->set_shell_sense( facet_shell, CUBIT_FORWARD );
      }

      surface_list += shell_surfaces;
      shell_list.append(shell_ptr);
    }

    // make a body out of it

    rv = FacetModifyEngine::instance()->make_facet_lump(shell_list,lump_ptr);
    if ( lump_ptr == NULL || rv != CUBIT_SUCCESS )
    {
      PRINT_ERROR("Problems building facet based lump entity.\n");
      rv = CUBIT_FAILURE;
      goto end_import_facets;
    }
    lump_list.append(lump_ptr);
    rv = FacetModifyEngine::instance()->make_facet_body(lump_list,bodysm_ptr);
    body_ptr = GeometryQueryTool::instance()->make_Body(bodysm_ptr);

    if ( body_ptr == NULL || rv != CUBIT_SUCCESS )
    {
      PRINT_ERROR("Problems building facet based body entity.\n");
      rv = CUBIT_FAILURE;
      goto end_import_facets;
    }

    if(gti->num_bodies() - prev_bod > 1)
         PRINT_INFO("Bodies successfully created.\n");
    else
       PRINT_INFO("Body successfully created.\n");
    PRINT_INFO("  Number of new vertices = %d\n", gti->num_ref_vertices() - prev_vert);
    PRINT_INFO("  Number of new curves = %d\n", gti->num_ref_edges() - prev_edge);
    PRINT_INFO("  Number of new surfaces = %d\n", gti->num_ref_faces() - prev_face);
    PRINT_INFO("  Number of new shells = %d\n", shell_facet_list.size());
    PRINT_INFO("  Number of new volumes = %d\n", gti->num_ref_volumes() - prev_vol);
    PRINT_INFO("  Number of new bodies = %d\n", gti->num_bodies() - prev_bod);

    if (!is_water_tight)
    {
      PRINT_WARNING("Volume generated does not completely close. 3D meshing (ie. hex/tet) will not be permitted.\n");
      PRINT_INFO("Hint: In some cases the \"stitch\" option on the import command may correct the problem.\n");
    }
  }  // end while ( stl_multiple_parts == true )

end_import_facets:
  if (conn != NULL)
    delete [] conn;
  for (ii=0; ii<shell_list.size(); ii++)
    delete shell_facet_list.get_and_step();

  quad_facet_list = qfacet_list;
  tri_facet_list = tfacet_list;
  return rv;
}



//===========================================================================
//Function Name: read_facets
//Member Type:  PUBLIC
//Description:  read facets from facet file
//===========================================================================
CubitStatus FacetQueryEngine::read_facets( const char * file_name,
                                           int *&conn,
                                           int &npoints,
                                           int &nquad, int &ntri,
                                           FacetFileFormat file_format )
{
  // open the file
  FILE *fp = fopen(file_name, "r");
  if (fp == NULL)
  {
    PRINT_ERROR("Could not open file %s for reading\n", file_name);
    return CUBIT_FAILURE;
  }

  PRINT_INFO("Reading facets...\n");

  // read the number of nodes

  int id;
  int nfacets = 0;
  int iline = 1;
  char line[128];
  if (fgets(line, 128, fp) == NULL)
  {
    PRINT_ERROR("Format error in facet file %s on line %d\n", file_name, iline);
    fclose( fp );
    return CUBIT_FAILURE;
  }

  int n = sscanf(line, "%d %d", &npoints, &nfacets);
  if (n < 1 || n > 2 && file_format == CUBIT_FACET_FILE )
  {
    PRINT_ERROR("Format error in facet file %s on line %d\n", file_name, iline);
    fclose( fp );
    return CUBIT_FAILURE;
  }
  if (npoints <= 0)
  {
    PRINT_ERROR("Expecting number of nodes in facet file %s on line %d\n", file_name, iline);
    fclose( fp );
    return CUBIT_FAILURE;
  }
  if (n==1)
  {
    nfacets = 0;
  }
  else if (nfacets <= 0)
  {
    PRINT_ERROR("Format error in facet file %s on line %d reading number of facets\n", file_name, iline);
    fclose( fp );
    return CUBIT_FAILURE;
  }

  if (init_hash_points( npoints ) == CUBIT_FAILURE)
  {
    PRINT_ERROR("Can't allocate memory for points in facet file %s on line %d\n", file_name, iline);
    fclose( fp );
    return CUBIT_FAILURE;
  }

  // read the nodes

  int ii;
  double xx, yy, zz;
  CubitPoint *new_point;
  for (ii=0; ii<npoints; ii++)
  {
    iline++;
    if (fgets(line, 128, fp)== NULL)
    {
      PRINT_ERROR("Format error in facet file %s on line %d\n", file_name, iline);
      fclose( fp );
      return CUBIT_FAILURE;
    }

    n = sscanf(line, "%d %lf %lf %lf", &id, &xx, &yy, &zz );
    if (n != 4)
    {
      PRINT_ERROR("Format error in facet file %s on line %d\n", file_name, iline);
      PRINT_INFO("  Expecting 1 integer and 3 doubles, but instead read %d values\n", n);
      fclose( fp );
      return CUBIT_FAILURE;
    }
    new_point = (CubitPoint *) new CubitPointData( xx, yy, zz );
    new_point->set_id( id );
    add_hash_point( new_point );
  }

    // read the number of facets

  if (nfacets == 0)
  {
    iline++;
    if (fgets(line, 128, fp) == NULL)
    {
      PRINT_ERROR("Format error in facet file %s on line %d reading number of facets\n", file_name, iline);
      fclose( fp );
      return CUBIT_FAILURE;
    }

    n = sscanf(line, "%d", &nfacets);
    if (n != 1)
    {
      PRINT_ERROR("Format error in facet file %s on line %d reading number of facets\n", file_name, iline);
      PRINT_INFO("  Expecting <num facets>\n");
      fclose( fp );
      return CUBIT_FAILURE;
    }
    if (nfacets <= 0)
    {
      PRINT_ERROR("Format error in facet file %s on line %d reading number of facets\n", file_name, iline);
      fclose( fp );
      return CUBIT_FAILURE;
    }
  }

  conn = new int [4*nfacets];
  if (!conn)
  {
    PRINT_ERROR("Can't allocate memory for facets in facet file %s on line %d\n", file_name, iline);
    fclose( fp );
    return CUBIT_FAILURE;
  }
  ntri = nquad = 0;

  // read the facets

  for (ii=0; ii<nfacets; ii++)
  {
    iline++;
    if (fgets(line, 128, fp) == NULL)
    {
      PRINT_ERROR("Format error in facet file %s on line %d\n", file_name, iline);
      fclose( fp );
      return CUBIT_FAILURE;
    }

    if( AVS_FILE == file_format )
    {
      n = sscanf( line, "%d %*d %*s %d %d %d %d", &id, &conn[4*ii],
                  &conn[4*ii+1], &conn[4*ii+2], &conn[4*ii+3] );

      if( n < 4 || n > 5 )
      {
        PRINT_ERROR("Format error in avs file %s on line %d reading facets\n", file_name, iline);
        PRINT_INFO("Expecting 6 or 7 values, but instead read %d values\n", n);
        fclose( fp );
        return CUBIT_FAILURE;
      }

        //duplicate the point
      if( n == 4 )
      {
        conn[4*ii+3] = conn[4*ii+2];
        ntri++;
      }
      else
      {
        nquad++;
      }
    }
    else
    {
      n = sscanf(line, "%d %d %d %d %d", &id, &conn[4*ii], &conn[4*ii+1],
                 &conn[4*ii+2], &conn[4*ii+3] );

      if (n < 4 || n > 5)
      {
        PRINT_ERROR("Format error in facet file %s on line %d reading facets\n", file_name, iline);
        PRINT_INFO("  Expecting 4 or 5 integers, but instead read %d values\n", n);
        PRINT_INFO("  For example:  <id> <i0> <i1> <i2> [<i3>]\n");
        fclose( fp );
        return CUBIT_FAILURE;
      }

        // for triangles -- duplicate the point
      if (n==4)
      {
        conn[4*ii+3] = conn[4*ii+2];
        ntri++;
      }
      else
      {
        nquad++;
      }
    }
  }

  fclose( fp );
  return CUBIT_SUCCESS;
}

//===========================================================================
//  Function: read_cholla_file
//  Purpose:  import the face mesh from the cholla file
//  Date:     11/28/2002
//  Author:   sjowen
//===========================================================================
CubitStatus FacetQueryEngine::read_cholla_file( const char *file_name,
                              double &feature_angle,
                              DLIList<CubitPoint *> &point_list,
                              DLIList<CubitFacet *> &facet_list )

{

  double angle = 135.0;
  int num_tri = 0;
  int num_quad = 0;
  int num_edge = 0;
  int num_vert = 0;
  int* tri_edge = NULL;
  int* quad_edge = NULL;
  int* edge_vert = NULL;
  double* vert = NULL;
  double* edge_ctrl_pts = NULL;
  double* tri_ctrl_pts = NULL;
  double* quad_ctrl_pts = NULL;
  int results_included = 0;

  // importMesh and resolveFaceVectors are in Cholla.h

  importMesh(file_name, &results_included, &angle, &num_tri, &num_quad, &num_edge,
             &num_vert, &tri_edge, &quad_edge, &edge_vert, &vert,
             &edge_ctrl_pts, &tri_ctrl_pts, &quad_ctrl_pts);

  feature_angle = angle;

  // create the points

  CubitPoint **point_array = new CubitPoint * [num_vert];
  double x, y, z;
  int ii;

  for (ii=0; ii<num_vert; ii++)
  {
    x = vert[ii*3];
    y = vert[ii*3 + 1];
    z = vert[ii*3 + 2];
    point_array[ii] = (CubitPoint *) new CubitPointData( x, y, z );
    point_list.append( point_array[ii] );
  }

  // create the edges

  CubitFacetEdge **edge_array = new CubitFacetEdge * [num_edge];
  int ip, jp;
  for (ii=0; ii<num_edge; ii++)
  {
    ip = edge_vert[ii*2];
    jp = edge_vert[ii*2 + 1];
    assert(ip < num_vert && jp < num_vert && ip >= 0 && jp >= 0);
    edge_array[ii] = (CubitFacetEdge *) new CubitFacetEdgeData( point_array[ip],
                                                                point_array[jp] );
  }

  // create the tri facets

  int begin_face;
  int jj, iedge;
  CubitFacetEdge *edges[4];
  CubitFacet *facet_ptr = NULL;
  CubitQuadFacet *qfacet_ptr = NULL;
  for (ii=0; ii<num_tri; ii++)
  {
    begin_face = 3 * ii;
    for(jj=0; jj<3; jj++)
    {
      iedge = tri_edge[begin_face+jj];
      edges[jj] = edge_array[iedge];
    }
    facet_ptr = (CubitFacet *)
      new CubitFacetData(edges[0], edges[1], edges[2]);
    facet_list.append(facet_ptr);
  }

  // create the quad facets

  for (ii=0; ii<num_quad; ii++)
  {
    begin_face = 4 * ii;
    for(jj=0; jj<4; jj++)
    {
      iedge = quad_edge[begin_face+jj];
      edges[jj] = edge_array[iedge];
    }
    qfacet_ptr = (CubitQuadFacet *)
      new CubitQuadFacetData(edges[0], edges[1], edges[2], edges[3]);
    facet_list.append(qfacet_ptr->get_tri_facet(0));
    facet_list.append(qfacet_ptr->get_tri_facet(1));
  }

  // delete the temp arrays

  delete [] point_array;
  delete [] edge_array;
  delete [] tri_edge;
  delete [] quad_edge;
  delete [] edge_vert;
  delete [] vert;
  delete [] edge_ctrl_pts;
  delete [] tri_ctrl_pts;
  delete [] quad_ctrl_pts;

  return CUBIT_SUCCESS;
}


//===========================================================================
//Function Name: is_close
//
//Member Type:  PRIVATE
//Description:  determine if one of the facets in the list is within a
//              certain distance of the point.
//===========================================================================
CubitBoolean
FacetQueryEngine::is_close(CubitVector &this_point,
                              DLIList<CubitFacet *>&facet_list,
                              CubitFacet *&lastFacet,
                              double tol)
{
  CubitBoolean isclose = CUBIT_FALSE;
  CubitBox bbox;
  CubitVector close_point;
  CubitFacet *facet;
  int ii;

  // set the first facet to be checked as the last one located

  if (lastFacet) {
    if (!facet_list.move_to(lastFacet)) {
      facet_list.reset();
    }
  }
  else {
    facet_list.reset();
  }

  CubitBoolean done = CUBIT_FALSE;

  // define a bounding box around the point

  CubitVector ptmin( this_point.x() - tol,
		                 this_point.y() - tol,
					           this_point.z() - tol );
	CubitVector ptmax( this_point.x() + tol,
		                 this_point.y() + tol,
					           this_point.z() + tol );

  for ( ii = facet_list.size(); ii > 0 && !done; ii-- ) {
	  facet = facet_list.get_and_step();

    // Try to trivially reject this facet with a bounding box test

	  bbox = facet->bounding_box();
	  if (ptmax.x() < bbox.minimum().x() ||
		    ptmin.x() > bbox.maximum().x()) {
      continue;
    }
    if (ptmax.y() < bbox.minimum().y() ||
		    ptmin.y() > bbox.maximum().y()) {
      continue;
    }
    if (ptmax.z() < bbox.minimum().z() ||
		    ptmin.z() > bbox.maximum().z()) {
      continue;
    }

    // Only facets that pass the bounding box test will get past here!

    // Project point to plane of the facet and determine its area coordinates

    CubitVector pt_on_plane;
    double dist_to_plane;
    FacetEvalTool::project_to_facet_plane( facet, this_point, pt_on_plane, dist_to_plane );

    CubitVector areacoord;
    FacetEvalTool::facet_area_coordinate( facet, pt_on_plane, areacoord );

    // If sign of areacoords are all positive then its inside the triangle
    // and we are done - go interpolate the point. (use an absolute
    // tolerance since the coordinates arenormalized)

    if (areacoord.x() > -GEOMETRY_RESABS &&
        areacoord.y() > -GEOMETRY_RESABS &&
        areacoord.z() > -GEOMETRY_RESABS) {
      FacetEvalTool::eval_facet( facet, areacoord, &close_point, NULL );
    }

    // otherwise find the closest vertex or edge to the projected point

    else if (areacoord.x() < GEOMETRY_RESABS) {
      if (areacoord.y() < GEOMETRY_RESABS) {
        FacetEvalTool::eval_point( facet, 2, close_point );
      }
      else if(areacoord.z() < GEOMETRY_RESABS) {
        FacetEvalTool::eval_point( facet, 1, close_point );
      }
      else {
        FacetEvalTool::eval_edge( facet, 1, 2, pt_on_plane, close_point );
      }
    }
    else if (areacoord.y() < GEOMETRY_RESABS) {
      if (areacoord.z() < GEOMETRY_RESABS) {
        FacetEvalTool::eval_point( facet, 0, close_point );
      }
      else {
        FacetEvalTool::eval_edge( facet, 2, 0, pt_on_plane, close_point);
      }
    }
    else {
      FacetEvalTool::eval_edge( facet, 0, 1, pt_on_plane, close_point );
    }

    // keep track of the minimum distance

    double dist = sqrt(sqr(close_point.x() - this_point.x()) +
                  sqr(close_point.y() - this_point.y()) +
                  sqr(close_point.z() - this_point.z()));
    if (dist <= tol) {
      done = CUBIT_TRUE;
      isclose = CUBIT_TRUE;
      lastFacet = facet;
    }
	}
  return isclose;
}




//=============================================================================
// Function Name: export_facets
// Member Type:  PUBLIC
// Description:  export a list of facets to a facet file for debugging purposes
// Author: sjowen
// Date: 10/29/01
//=============================================================================
CubitStatus FacetQueryEngine::export_facets(DLIList<CubitFacet*> &facet_list,
                                               char *filename)
{

  // open the file for writing

  FILE *fp = fopen(filename, "w");
  if (!fp)
  {
    PRINT_ERROR("Couldn't open file %s for writing.\n", filename);
    return CUBIT_FAILURE;
  }

  // make list of all points

  CubitPoint *point_ptr;
  CubitFacet *facet_ptr;
  DLIList<CubitPoint*>temp_point_list;
  int ii, jj;
  for (ii=0; ii<facet_list.size(); ii++)
  {
    facet_ptr = facet_list.get_and_step();
    for(jj=0; jj<3; jj++)
    {
      point_ptr = facet_ptr->point(jj);
      point_ptr->marked(-1);
      temp_point_list.append(point_ptr);
    }
  }

  // assign ids to points and make list unique

  int id = 0;
  DLIList<CubitPoint*>point_list;
  for (ii=0; ii<temp_point_list.size(); ii++)
  {
    point_ptr = temp_point_list.get_and_step();
    if (point_ptr->marked() == -1)
    {
      point_ptr->marked(id++);
      point_list.append( point_ptr );
    }
  }

  // write the points

  CubitVector coords;
  fprintf(fp, "%d\n", point_list.size());
  point_list.reset();
  for(ii=0; ii<point_list.size(); ii++)
  {
    point_ptr = point_list.get_and_step();
    coords = point_ptr->coordinates();
    fprintf(fp,"%d %f %f %f\n", ii, coords.x(), coords.y(), coords.z());
  }

  // write the facets

  int indx[3];
  fprintf(fp, "%d\n", facet_list.size());
  for (ii=0; ii<facet_list.size(); ii++)
  {
    facet_ptr = facet_list.get_and_step();
    for (jj=0; jj<3; jj++)
    {
      indx[jj] = facet_ptr->point(jj)->marked();
    }
    fprintf(fp,"%d %d %d %d\n", ii, indx[0], indx[1], indx[2]);
  }

  //- unmark the points

  for (ii=0; ii<point_list.size(); ii++)
    point_list.get_and_step()->marked(0);

  fclose(fp);
  return CUBIT_SUCCESS;
}

//=============================================================================
//Function:  init_hash_point (PRIVATE)
//Description: create a hash array of all points.  They are hashed based on the
//             id of the point
//Author: sjowen
//Date: 4/3/01
//=============================================================================
CubitStatus FacetQueryEngine::init_hash_points( int num_points )
{
  /* === find the next highest prime number */

  hashPointSize = num_points / 10;
  int i;
  if (hashPointSize < 13) hashPointSize = 13;
  else
  {
    i=2;
    while (i<hashPointSize*0.5 + 1) {
      if (hashPointSize % i == 0) {
        i=2;
        hashPointSize++;
      }
      else {
        i++;
      }
    }
  }
  hashPointArray = new DLIList<CubitPoint*>[hashPointSize];

  return CUBIT_SUCCESS;
}

//=============================================================================
//Function:  add_hash_point (PRIVATE)
//Description: hash the point into the hash table
//Author: sjowen
//Date: 4/3/02
//=============================================================================
CubitStatus FacetQueryEngine::add_hash_point( CubitPoint *point_ptr )
{
  int id = point_ptr->id();
  int key = get_hash_key( id );
  hashPointArray[key].append( point_ptr );
  return CUBIT_SUCCESS;
}

//=============================================================================
//Function:  get_hash_point (PRIVATE)
//Description: retreive the CubitPoint pointer from the id
//Author: sjowen
//Date: 4/3/02
//=============================================================================
CubitPoint *FacetQueryEngine::get_hash_point( int id )
{
  int key = get_hash_key( id );
  int ii;
  int found = 0;
  CubitPoint *point_ptr;
  for (ii=0; ii<hashPointArray[key].size() && !found; ii++)
  {
    point_ptr = hashPointArray[key].get_and_step();
    if (point_ptr->id() == id)
    {
      return point_ptr;
    }
  }
  return (CubitPoint*)NULL;
}

//=============================================================================
//Function:  delete_hash_points (PRIVATE)
//Description: delete the hash points stuff
//Author: sjowen
//Date: 4/3/02
//=============================================================================
void FacetQueryEngine::delete_hash_points( )
{
  if (hashPointArray)
   delete [] hashPointArray;
  hashPointArray = NULL;
  hashPointSize = 0;
}

//=============================================================================
//Function:  get_hash_key (PRIVATE)
//Description:
//Author: sjowen
//Date: 4/3/02
//=============================================================================
int FacetQueryEngine::get_hash_key( int id )
{
  int key = id % hashPointSize;
  return key;
}

//===========================================================================
//Function Name: get_all_hash_points
//Member Type:  PUBLIC
//Description:  makes a single list from all points in the hash table
//===========================================================================
CubitStatus FacetQueryEngine::get_all_hash_points(
  DLIList<CubitPoint *> &point_list ) // return point list
{
  CubitStatus rv = CUBIT_SUCCESS;

  int ii, jj;
  CubitPoint *point_ptr;
  for(jj=0; jj<hashPointSize; jj++)
  {
    for (ii=0; ii<hashPointArray[jj].size(); ii++)
    {
      point_ptr = hashPointArray[jj].get_and_step();
      point_list.append(point_ptr);
    }
  }

  return rv;
}

//===========================================================================
//Function Name: save_facets
//Member Type:  PUBLIC
//Description:  save the facets from the entities to a Cubit file
//Author:       sjowen
//Date:         1/16/2003
//===========================================================================
CubitStatus FacetQueryEngine::save_facets(
  FILE *fp,  // file we are dumping to
  DLIList<FacetSurface*> facet_surfaces,
  DLIList<FacetCurve*> facet_curves,
  DLIList<FacetPoint*> facet_points ) // save facets from these entities
{
  CubitStatus rv = CUBIT_SUCCESS;

  DLIList<CubitFacet *> facet_list;
  DLIList<CubitFacetEdge *> edge_list;
  DLIList<CubitPoint *> point_list;

  // get a unique list of all facets, edges and points
  rv = gather_facets(facet_surfaces, facet_curves, facet_points,
                     facet_list, edge_list, point_list );

  if (rv != CUBIT_SUCCESS)
    return rv;

  // dump the facet entities the file

  rv = dump_facets( fp, facet_list, edge_list, point_list );
  if (rv != CUBIT_SUCCESS)
    return rv;

  // dump the CurveFacetEvalTools and FacetEval Tools
  rv = save_eval_tools( fp, facet_surfaces, facet_curves );

  return rv;
}

//===========================================================================
//Function Name: save_eval_tools
//Member Type:  PUBLIC
//Description:  go through the RefEntities and save the CurveFacetEvalTools
//              and FacetEvalTools.  These are essentially lists of facets
//              that are owned by the curve or surface
//Authors:      sjowen & cdernst
//Date:         1/21/2003
//===========================================================================
CubitStatus FacetQueryEngine::save_eval_tools(
  FILE *fp,        // cubit file we are dumping to
  DLIList<FacetSurface*> facet_surfaces,
  DLIList<FacetCurve*>  facet_curves )
{

  int ii;
  FacetSurface *fsurf_ptr;
  FacetCurve *fcurv_ptr;

  //separate out eval tools and set ids
  DLIList<CurveFacetEvalTool*> cfe_tools;
  DLIList<FacetEvalTool*> fe_tools;
  int ft_id = 0;
  int cft_id = 0;

  facet_surfaces.reset();
  for (ii=facet_surfaces.size(); ii--; )
  {
    fsurf_ptr = facet_surfaces.get_and_step();
    FacetEvalTool *feval_tool = fsurf_ptr->get_eval_tool();
    if( feval_tool )
    {
      feval_tool->set_output_id( ft_id++ );
      fe_tools.append( feval_tool );
    }
  }

  facet_curves.reset();
  for (ii=facet_curves.size(); ii--; )
  {
    fcurv_ptr = facet_curves.get_and_step();
    CurveFacetEvalTool *ceval_tool = fcurv_ptr->get_eval_tool();
    if( ceval_tool )
    {
      ceval_tool->set_output_id( cft_id++ );
      cfe_tools.append( ceval_tool );
    }
  }

  //Write out number of FacetEvalTools
  int count = fe_tools.size();
  NCubitFile::CIOWrapper cio(fp);
  cio.Write( reinterpret_cast<UnsignedInt32*>(&count), 1 );

  //Write out each FacetEvalTool
  fe_tools.reset();
  for( ii=0; ii<fe_tools.size(); ii++)
  {
    FacetEvalTool *feval_tool = fe_tools.get_and_step();
    feval_tool->save( fp );
  }

  //Write out number of CurveFacetEvalTools
  count = cfe_tools.size();
  cio.Write( reinterpret_cast<UnsignedInt32*>(&count), 1 );

  //Write out each CurveFacetEvalTool
  cfe_tools.reset();
  for( ii=0; ii<cfe_tools.size(); ii++)
  {
    CurveFacetEvalTool *ceval_tool = cfe_tools.get_and_step();
    ceval_tool->save( fp );
  }

  return CUBIT_SUCCESS;
}

//===============================================================================
// Function   : restore_eval_tools
// Member Type: PUBLIC
// Description: restore the curve and surface eval tools
// Author     : sjowen
// Date       : 1/26/03
//===============================================================================
CubitStatus FacetQueryEngine::restore_eval_tools(
  FILE *fp,  // CUB file we are currently reading
  unsigned int endian,
  int num_facets, // list of facet entities already read from the CUB file
  int num_edges,  // that we will assign to the eval tools.  Index into
  int num_points, // the array is the same as ID
  CubitFacet **facets,
  CubitFacetEdge **edges,
  CubitPoint **points,
  int &num_cfet,
  int &num_fet,
  CurveFacetEvalTool **&cfeval_tools,
  FacetEvalTool **&feval_tools)   // return list of curve and surface tools
{
  NCubitFile::CIOWrapper cio(endian, fp);

  int ii;
  FacetEvalTool *fsurf_ptr;
  CurveFacetEvalTool *fcurv_ptr;

  // read the number of facet eval tools
  // we are about to read from the CUB file
  cio.Read(reinterpret_cast<UnsignedInt32*>(&num_fet), 1);

  feval_tools = new FacetEvalTool *[num_fet];

  // read all the facet surface eval tools

  for (ii=0; ii<num_fet; ii++)
  {
    /*
    // read the topo ID from the CUB file.  This will be stored as the
    // FacetEvalTool toolID and serve as the means to associate this
    // eval tool with it's topology entity
    UnsignedInt32 surf_id = 0;
    cio.Read(&surf_id, 1);
    int topo_id = (int)surf_id; */

    // read the facet eval tool data and create the tool
    fsurf_ptr  = new FacetEvalTool();
    fsurf_ptr->restore(fp, endian, num_facets,
                       num_edges, num_points,
                       facets, edges, points);
    if (fsurf_ptr == NULL)
    {
      PRINT_ERROR("Error restoring mesh-based geometry\n");
      return CUBIT_FAILURE;
    }
    feval_tools[ii] = fsurf_ptr;
  }

  // read the number of curves
  cio.Read(reinterpret_cast<UnsignedInt32*>(&num_cfet), 1);
  cfeval_tools = new CurveFacetEvalTool *[num_cfet];

  // read all the curve facet eval tools.  Same thing as surfaces

  for (ii=0; ii<num_cfet; ii++)
  {
    // read the facet eval tool data and create the tool
    fcurv_ptr = new CurveFacetEvalTool();
    fcurv_ptr->restore(fp, endian, num_edges, num_points,
                       edges, points, num_fet, feval_tools );
    if (fcurv_ptr == NULL)
    {
      PRINT_ERROR("Error restoring mesh-based geometry\n");
      return CUBIT_FAILURE;
    }

    cfeval_tools[ii] = fcurv_ptr;
  }

  return CUBIT_SUCCESS;
}


//===========================================================================
//Function Name: dump_facets
//Member Type:  PUBLIC
//Description:  dump the facets and associated data to the cubit file
//Authors:      sjowen & cdernst
//Date:         1/16/2003
//===========================================================================
CubitStatus FacetQueryEngine::dump_facets(
  FILE *fp,                          // file we are dumping to
  DLIList<CubitFacet *> &facet_list,  // return unique list of facets
  DLIList<CubitFacetEdge *> &edge_list,  // return unique list of edges
  DLIList<CubitPoint *> &point_list )  // return unique list of points
{

  NCubitFile::CIOWrapper cio(fp);
  typedef NCubitFile::UnsignedInt32 UnsignedInt32;

  // Gather CubitPoint Data
  int npoints = point_list.size();
  int nnormals = 0;
  double* uvs_array = new double [npoints * 3];
  double* coord_array  = new double [npoints * 3];
  //double* du_array = new double [npoints * 3];
  //double* dv_array = new double [npoints * 3];
  int c_zero_points = 0;
  int c_zero_int_data_size = 0;
  int c_zero_double_data_size = 0;
  int ii;
  DLIList<CubitPoint*> has_td_boundary_point;
  DLIList<CubitPoint*> points_with_normal;

  TDFacetBoundaryPoint *td_fbp;
  CubitPoint *cp_ptr;
  CubitVector coord, *normal;
  point_list.reset();
  for (ii=0; ii<point_list.size(); ii++)
  {
    cp_ptr = point_list.get_and_step();

    // write out uVal, vVal, sizeVal
    uvs_array[ii*3] = cp_ptr->u();
    uvs_array[ii*3+1] = cp_ptr->v();
    uvs_array[ii*3+2] = cp_ptr->size();
    // coordinates
    coord = cp_ptr->coordinates();
    coord_array[ii*3]   = coord.x();
    coord_array[ii*3+1] = coord.y();
    coord_array[ii*3+2] = coord.z();
    // surfNormal
    if( cp_ptr->normal_ptr() )
      points_with_normal.append( cp_ptr );

    /*
    // surfU
    du = cp_ptr->du();
    du_array[ii*3]   = du.x();
    du_array[ii*3+1] = du.y();
    du_array[ii*3+2] = du.z();
    // surfV
    dv = cp_ptr->dv();
    dv_array[ii*3]   = dv.x();
    dv_array[ii*3+1] = dv.y();
    dv_array[ii*3+2] = dv.z(); */

    if ((td_fbp = TDFacetBoundaryPoint::get_facet_boundary_point( cp_ptr ))!= NULL)
    {
      has_td_boundary_point.append( cp_ptr );
      c_zero_points++;
      td_fbp->get_boundary_point_data_size( c_zero_int_data_size,
                                            c_zero_double_data_size );
    }
  }

  // Normals
  nnormals = points_with_normal.size();
  double* normal_array = NULL;
  int* normal_ids = NULL;
  points_with_normal.reset();
  if( nnormals > 0 )
  {
    normal_array = new double [nnormals* 3];
    normal_ids = new int[nnormals];
    points_with_normal.reset();
    for( ii=0; ii<nnormals; ii++)
    {
      cp_ptr = points_with_normal.get_and_step();
      normal_ids[ii] = cp_ptr->id();
      normal = cp_ptr->normal_ptr();
      normal_array[ii*3]   = normal->x();
      normal_array[ii*3+1] = normal->y();
      normal_array[ii*3+2] = normal->z();
    }
  }

  // write arrays  --CubitPoint
  cio.Write(reinterpret_cast<UnsignedInt32*>(&npoints), 1);
  if( npoints > 0 )
  {
    cio.Write(coord_array, npoints*3);
    cio.Write(uvs_array, npoints*3);
  }

  // clean up
  delete [] uvs_array;
  delete [] coord_array;
  uvs_array = NULL;
  coord_array = NULL;

  // write normals & ids of points to which normals belong
  cio.Write( reinterpret_cast<UnsignedInt32*>(&nnormals), 1 );
  if( nnormals > 0 )
  {
    cio.Write(normal_array, nnormals*3);
    cio.Write(reinterpret_cast<UnsignedInt32*>(normal_ids), nnormals);
    delete [] normal_array;
    delete [] normal_ids;
    normal_array = NULL;
    normal_ids = NULL;
  }

  // Gather CubitFacetEdge Data
  int jj, idx;
  int nedges = edge_list.size();
  int nctrl_pts = nedges * NUM_EDGE_CPTS;
  CubitFacetEdge *edge_ptr;
  CubitVector *edge_ctrl_pts;
  double* control_points = new double [nctrl_pts * 3];
  UnsignedInt32* edge_vert_array = new UnsignedInt32 [nedges * 2];
  for(ii=0; ii<nedges; ii++)
  {
    // pointArray[2]
    edge_ptr = edge_list.get_and_step();
    edge_vert_array[ii*2] = edge_ptr->point(0)->id();
    edge_vert_array[ii*2+1] = edge_ptr->point(1)->id();
    edge_ctrl_pts = edge_ptr->control_points();

    // controlPoints[3]
    if (nctrl_pts > 0)
    {
      if (!edge_ctrl_pts)
      {
        nctrl_pts = 0;
        continue;
      }
      for(jj=0; jj<NUM_EDGE_CPTS; jj++)
      {
        idx = (ii*NUM_EDGE_CPTS+jj)*3;
        control_points[idx]   = edge_ctrl_pts[jj].x();
        control_points[idx+1] = edge_ctrl_pts[jj].y();
        control_points[idx+2] = edge_ctrl_pts[jj].z();
      }
    }
  }

  // write arrays  --CubitFacetEdge
  cio.Write(reinterpret_cast<UnsignedInt32*>(&nedges), 1);
  if( nedges > 0 )
  {
    cio.Write(edge_vert_array, nedges*2 );
    cio.Write(reinterpret_cast<UnsignedInt32*>(&nctrl_pts), 1);
    if (nctrl_pts > 0)
      cio.Write(control_points, nctrl_pts*3);
  }

  // clean up
  delete [] edge_vert_array;
  delete [] control_points;
  edge_vert_array = NULL;
  control_points = NULL;

  // Gather CubitFacet Data
  int nfacets = facet_list.size();
  nctrl_pts = nfacets * NUM_TRI_CPTS;
  CubitFacet *facet_ptr;
  CubitVector *facet_ctrl_pts;
  control_points = new double [nctrl_pts * 3];
  UnsignedInt32* facet_edge_array = new UnsignedInt32 [nfacets * 3];
  int *int_data = new int[ nfacets * 2 ];
  facet_list.reset();
  for (ii=0; ii<nfacets; ii++)
  {
    facet_ptr = facet_list.get_and_step();

    // is Flat and isBackwards
    int_data[ii*2] = facet_ptr->is_flat();
    int_data[ii*2+1] = facet_ptr->is_backwards();

    // edgeArray[3]
    facet_edge_array[ii*3]   = facet_ptr->edge(0)->id();
    facet_edge_array[ii*3+1] = facet_ptr->edge(1)->id();
    facet_edge_array[ii*3+2] = facet_ptr->edge(2)->id();
    facet_ctrl_pts = facet_ptr->control_points();
    if(nctrl_pts > 0)
    {
      if (!facet_ctrl_pts)
      {
        nctrl_pts = 0;
        continue;
      }
      for(jj=0; jj<NUM_TRI_CPTS; jj++)
      {
        idx = (ii*NUM_TRI_CPTS+jj)*3;
        control_points[idx]   = facet_ctrl_pts[jj].x();
        control_points[idx+1] = facet_ctrl_pts[jj].y();
        control_points[idx+2] = facet_ctrl_pts[jj].z();
      }
    }
  }

  // write arrays  --CubitFacet
  cio.Write(reinterpret_cast<UnsignedInt32*>(&nfacets), 1);
  if( nfacets > 0 )
  {
    cio.Write(facet_edge_array, nfacets*3);
    cio.Write( reinterpret_cast<UnsignedInt32*>(int_data), nfacets*2 );
    cio.Write(reinterpret_cast<UnsignedInt32*>(&nctrl_pts), 1);
    if (nctrl_pts > 0)
      cio.Write(control_points, nctrl_pts*3);
  }

  // clean up
  delete [] facet_edge_array;
  delete [] control_points;
  delete [] int_data;
  facet_edge_array = NULL;
  control_points = NULL;
  int_data = NULL;


  // At points along the boundary (C0 continuity) there may be
  // multiple normals.  Write this data too.  First gether it
  // into two arrays and then dump.

  cio.Write(reinterpret_cast<UnsignedInt32*>(&c_zero_points), 1);
  if (c_zero_points > 0)
  {
    int *c_zero_int_data = new int [c_zero_int_data_size];
    double *c_zero_double_data = new double [c_zero_double_data_size];
    point_list.reset();
    int iidx = 0;
    int didx = 0;
    point_list.reset();
    for (ii=0; ii<point_list.size(); ii++)
    {
      cp_ptr = point_list.get_and_step();
      td_fbp = TDFacetBoundaryPoint::get_facet_boundary_point( cp_ptr );
      if (td_fbp != NULL)
      {
        td_fbp->get_boundary_point_data( c_zero_int_data,
                                         c_zero_double_data,
                                         iidx, didx );
      }
    }

    // dump the int data array

    cio.Write(reinterpret_cast<UnsignedInt32*>(&c_zero_int_data_size), 1);
    //convert to UnsignedInt32
    UnsignedInt32 *c_zero_int32_data = new UnsignedInt32 [c_zero_int_data_size];
    for (ii=0; ii<c_zero_int_data_size; ii++)
      c_zero_int32_data[ii] = (UnsignedInt32)c_zero_int_data[ii];
    cio.Write(c_zero_int32_data, c_zero_int_data_size);

    // dump the double array

    cio.Write(reinterpret_cast<UnsignedInt32*>(&c_zero_double_data_size), 1);
    cio.Write(c_zero_double_data, c_zero_double_data_size);

    // clean up
    delete [] c_zero_int_data;
    delete [] c_zero_int32_data;
    delete [] c_zero_double_data;
  }

  return CUBIT_SUCCESS;
}

//===========================================================================
//Function Name: gather_facets
//Member Type:  PUBLIC
//Description:  get a unique list of facets from the entities
//Author:       sjowen
//Date:         1/16/2003
//===========================================================================
CubitStatus FacetQueryEngine::gather_facets(
  DLIList<FacetSurface *> facet_surfaces,
  DLIList<FacetCurve *> facet_curves,
  DLIList<FacetPoint *> facet_points,
  DLIList<CubitFacet *> &facet_list,  // return unique list of facets
  DLIList<CubitFacetEdge *> &edge_list,  // return unique list of edges
  DLIList<CubitPoint *> &point_list)   // return unique list of points
{
  int ii, jj, kk;
  FacetSurface *fsurf_ptr;
  FacetCurve *fcurv_ptr;
  FacetPoint *fpoint_ptr;
  DLIList<CubitFacet *>temp_facets;
  DLIList<CubitPoint *>temp_points;
  DLIList<CubitFacetEdge *>temp_edges;

  CubitFacet *facet_ptr;
  CubitFacetEdge *edge_ptr;
  CubitPoint *cp_ptr;

  //unmark facets/edges/points on surface
  for(ii=facet_surfaces.size(); ii--;)
  {
    fsurf_ptr = facet_surfaces.get_and_step();
    temp_facets.clean_out();
    temp_points.clean_out();
    fsurf_ptr->get_my_facets(temp_facets, temp_points);
    for(jj=0; jj<temp_facets.size(); jj++)
    {
      facet_ptr = temp_facets.get_and_step();
      facet_ptr->marked( 0 );
      for (kk=0; kk<3; kk++)
      {
        edge_ptr = facet_ptr->edge( kk );
        edge_ptr->marked( 0 );
      }
    }
    for (jj=0; jj<temp_points.size(); jj++)
    {
      cp_ptr = temp_points.get_and_step();
      cp_ptr->marked( 0 );
    }
  }

  //unmark facet-edges/points on curves
  for(ii=facet_curves.size(); ii--; )
  {
    fcurv_ptr = facet_curves.get_and_step();
    temp_edges.clean_out();
    fcurv_ptr->get_facets( temp_edges );
    for (jj=0; jj<temp_edges.size(); jj++)
    {
      edge_ptr = temp_edges.get_and_step();
      edge_ptr->marked( 0 );
      edge_ptr->point( 0 )->marked( 0 );
      edge_ptr->point( 1 )->marked( 0 );
    }
  }

  //unmark facet-points on points
  for(ii=facet_points.size(); ii--; )
  {
    fpoint_ptr = facet_points.get_and_step();
    if (fpoint_ptr != NULL)
    {
      cp_ptr = fpoint_ptr->get_cubit_point();
      cp_ptr->marked( 0 );
    }
  }


  // make unique lists
  for(ii=facet_surfaces.size(); ii--;)
  {
    fsurf_ptr = facet_surfaces.get_and_step();
    temp_facets.clean_out();
    temp_points.clean_out();
    fsurf_ptr->get_my_facets(temp_facets, temp_points);
    for(jj=0; jj<temp_facets.size(); jj++)
    {
      facet_ptr = temp_facets.get_and_step();
      if (0 == facet_ptr->marked())
      {
        facet_ptr->marked( 1 );
        facet_list.append( facet_ptr );
      }
      for (kk=0; kk<3; kk++)
      {
        edge_ptr = facet_ptr->edge( kk );
        if (0 == edge_ptr->marked())
        {
          edge_ptr->marked( 1 );
          edge_list.append( edge_ptr );
        }
      }
    }
    for (jj=0; jj<temp_points.size(); jj++)
    {
      cp_ptr = temp_points.get_and_step();
      if (0 == cp_ptr->marked())
      {
        cp_ptr->marked(1);
        point_list.append( cp_ptr );
      }
    }
  }

  for(ii=facet_curves.size(); ii--;)
  {
    fcurv_ptr = facet_curves.get_and_step();
    temp_edges.clean_out();
    fcurv_ptr->get_facets( temp_edges );
    for (jj=0; jj<temp_edges.size(); jj++)
    {
      edge_ptr = temp_edges.get_and_step();
      if (0 == edge_ptr->marked())
      {
        edge_ptr->marked(1);
        edge_list.append( edge_ptr );
        for (kk=0; kk<2; kk++)
        {
          cp_ptr = edge_ptr->point( kk );
          if (0 == cp_ptr->marked())
          {
            cp_ptr->marked(1);
            point_list.append(cp_ptr);
          }
        }
      }
    }
  }

  for(ii=facet_points.size(); ii--;)
  {
    fpoint_ptr = facet_points.get_and_step();
    cp_ptr = fpoint_ptr->get_cubit_point();
    if (0 == cp_ptr->marked())
    {
      cp_ptr->marked( 1 );
      point_list.append( cp_ptr );
    }
  }

  // set the IDs of the facet entities and clear the marks
  facet_list.reset();
  for (ii=0; ii<facet_list.size(); ii++)
  {
    facet_ptr = facet_list.get_and_step();
    facet_ptr->set_id(ii);
    facet_ptr->marked(0);
  }
  edge_list.reset();
  for (ii=0; ii<edge_list.size(); ii++)
  {
    edge_ptr = edge_list.get_and_step();
    edge_ptr->set_id(ii);
    edge_ptr->marked(0);
  }
  for (ii=0; ii<point_list.size(); ii++)
  {
    cp_ptr = point_list.get_and_step();
    cp_ptr->set_id(ii);
    cp_ptr->marked(0);
  }

  return CUBIT_SUCCESS;
}

//=============================================================================
//Function:   create_super_facet_bounding_box(PUBLIC)
//Description: Find the bounding box of a list of BodySMs
//Author: jdfowle
//Date: 12/15/03
//=============================================================================
CubitStatus FacetQueryEngine::create_super_facet_bounding_box(
                                DLIList<BodySM*>& body_list,
                                CubitBox& super_box )
{
BodySM *bodySM;
int i;
CubitStatus status = CUBIT_SUCCESS;

  body_list.reset();
  for ( i = 0; i < body_list.size(); i++ ) {
    bodySM = body_list.get_and_step();  
    create_facet_bounding_box(bodySM,super_box);
  }

  return status;
}

//=============================================================================
//Function:   create_facet_bounding_box(PUBLIC)
//Description: Find the bounding box of a BodySM
//Author: jdfowle
//Date: 12/15/03
//=============================================================================
CubitStatus FacetQueryEngine::create_facet_bounding_box(
                                BodySM* bodySM,
                                CubitBox& bbox )
{
  FacetBody *fbody_ptr;
  int j;
  DLIList<FacetSurface*> facet_surf_list;
  DLIList<CubitFacet*> facet_list;
  DLIList<CubitPoint*> point_list;
  FacetSurface *facet_surface;
  CubitBox surf_bbox, total_box;
  CubitStatus status = CUBIT_FAILURE;

    fbody_ptr = dynamic_cast<FacetBody *>(bodySM);
    fbody_ptr->get_surfaces(facet_surf_list);
    for ( j = 0; j < facet_surf_list.size(); j++ ) {
      facet_surface = facet_surf_list.get_and_step();
      facet_list.clean_out();
      point_list.clean_out();
      facet_surface->get_my_facets(facet_list,point_list);
      status = FacetDataUtil::get_bbox_of_points(point_list,surf_bbox);
      if ( j == 0 ) total_box = surf_bbox;
      else 
        total_box |= surf_bbox;
    }
  bbox |= total_box;  
  return status;
}

const char* fqe_xform_err = "Transform not implemented for facet geometry.\n";
CubitStatus FacetQueryEngine::restore_transform( BodySM* body )
{
  FacetBody* facetbod = dynamic_cast<FacetBody*>(body);
  return facetbod ? facetbod->restore( ) : CUBIT_FAILURE;
}
CubitStatus FacetQueryEngine::translate( BodySM* body, const CubitVector& d )
{
  FacetBody* facetbod = dynamic_cast<FacetBody*>(body);
  return facetbod ? facetbod->move( d.x(), d.y(), d.z() ) : CUBIT_FAILURE;
}
CubitStatus FacetQueryEngine::rotate( BodySM* body, const CubitVector& v, double a )
{
  FacetBody* facetbod = dynamic_cast<FacetBody*>(body);
  return facetbod ? facetbod->rotate( v.x(), v.y(), v.z(), a ) : CUBIT_FAILURE;
}
CubitStatus FacetQueryEngine::scale( BodySM* body, double factor )
{
  FacetBody* facetbod = dynamic_cast<FacetBody*>(body);
  return facetbod ? facetbod->scale( factor ) : CUBIT_FAILURE;
}
CubitStatus FacetQueryEngine::scale( BodySM* body, const CubitVector& f )
{
  FacetBody* facetbod = dynamic_cast<FacetBody*>(body);
  return facetbod ? facetbod->scale( f.x(), f.y(), f.z() ) : CUBIT_FAILURE;
}

CubitStatus FacetQueryEngine::translate( GeometryEntity* , const CubitVector&  )
{
  PRINT_ERROR(fqe_xform_err);
  return CUBIT_FAILURE;
}
CubitStatus FacetQueryEngine::rotate( GeometryEntity* , const CubitVector& , double  )
{
  PRINT_ERROR(fqe_xform_err);
  return CUBIT_FAILURE;
}
CubitStatus FacetQueryEngine::scale( GeometryEntity* , double  )
{
  PRINT_ERROR(fqe_xform_err);
  return CUBIT_FAILURE;
}
CubitStatus FacetQueryEngine::scale( GeometryEntity* , const CubitVector&  )
{
  PRINT_ERROR(fqe_xform_err);
  return CUBIT_FAILURE;
}
CubitStatus FacetQueryEngine::reflect( GeometryEntity* , const CubitVector&  )
{
  PRINT_ERROR(fqe_xform_err);
  return CUBIT_FAILURE;
}




//===============================================================================
// Function   : bodies_overlap
// Member Type: PUBLIC
// Description: determine if facet-based bodies overlap
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitBoolean FacetQueryEngine::bodies_overlap (BodySM * /*body_ptr_1*/,
                                                BodySM * /*body_ptr_2*/ ) const
{
  PRINT_ERROR("Currently, Cubit is unable to determine overlap for Facet-based geometry.\n");
  return CUBIT_FALSE;
    //The following should work, but due to bug 5599, we are not
    // using it.  Having this code in the default would cause the
    // GeometryPowerTool to rely on Facetbool, which isn't robust
    // yet.
    /*
FacetboolInterface *fbint;
CubitFacetboolOp op = CUBIT_FB_INTERSECTION;
bool keep_old = true;
CubitStatus status;
BodySM* body_out = 0;
bool intersection_found;

  fbint = new FacetboolInterface;  
  status = fbint->dofacetboolean_2bodies(body_ptr_1,body_ptr_2,body_out,
                                       keep_old,intersection_found,op);
  if ( status == CUBIT_FAILURE ) return CUBIT_FALSE;
  else return CUBIT_TRUE;
    */
}

CubitBoolean FacetQueryEngine::volumes_overlap (Lump * /*lump_ptr_1*/,
                                                Lump * /*lump_ptr_2*/ ) const
{
  PRINT_ERROR("Currently, Cubit is unable to determine overlap for Facet-based geometry.\n");
  return CUBIT_FALSE;
}

//EOF
