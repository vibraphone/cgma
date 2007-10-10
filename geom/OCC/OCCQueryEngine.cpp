//-------------------------------------------------------------------------
// Filename      : OCCQueryEngine.cpp
//
// Purpose       : Implementation of the OCCQueryEngine class.
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
#include "config.h"
#include "BRep_Tool.hxx"
#include "gp_Pnt.hxx"
#include "OCCQueryEngine.hpp"
#include "OCCModifyEngine.hpp"
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
#include "OCCPoint.hpp"
#include "OCCCurve.hpp"
#include "OCCCoEdge.hpp"
#include "OCCLoop.hpp"
#include "OCCSurface.hpp"
#include "OCCShell.hpp"
#include "OCCLump.hpp"
#include "OCCBody.hpp"
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
#include <stdio.h>
#include <errno.h>

#include <BRep_Builder.hxx>
#include <BRepTools.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_CompSolid.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Vertex.hxx>
#include <Bnd_Box.hxx>
#include <BndLib_AddSurface.hxx>
#include <Precision.hxx>
#include <BRepAdaptor_Surface.hxx>
using namespace NCubitFile;

OCCQueryEngine* OCCQueryEngine::instance_ = NULL;
int OCCQueryEngine::hashPointSize = 0;
DLIList<CubitPoint*> *OCCQueryEngine::hashPointArray = NULL;

const int OCCQueryEngine::OCCQE_MAJOR_VERSION = 6;
const int OCCQueryEngine::OCCQE_MINOR_VERSION = 1;
const int OCCQueryEngine::OCCQE_SUBMINOR_VERSION = 0;

TopTools_DataMapOfShapeInteger *OCCQueryEngine::OCCMap = new TopTools_DataMapOfShapeInteger;
TopTools_DataMapOfShapeInteger *OCCQueryEngine::OCCMapr = new TopTools_DataMapOfShapeInteger;
DLIList<TopologyBridge*> *OCCQueryEngine::CGMList = new DLIList<TopologyBridge*>;

//================================================================================
// Description:
// Author     :
// Date       :
//================================================================================
OCCQueryEngine* OCCQueryEngine::instance()
{
  if (instance_ == NULL ) {
      instance_ = new OCCQueryEngine;
   }
  return instance_;
}

//================================================================================
// Description:  default constructor
// Author     :
// Date       :
//================================================================================
OCCQueryEngine::OCCQueryEngine()
{
  GeometryQueryTool::instance()->add_gqe( this );
}

//================================================================================
// Description:  destructor
// Author     :
// Date       :
//================================================================================
OCCQueryEngine::~OCCQueryEngine()
{
  instance_ = NULL;
}

//================================================================================
// Description:  can_delete_bodies
// Author     : sjowen
// Date       : 4/25/02
//================================================================================
CubitBoolean OCCQueryEngine::can_delete_bodies(DLIList<Body*>body_list)
{
  CubitBoolean delete_ok = CUBIT_TRUE;
  int ii;
  for (ii=0; ii<body_list.size() && delete_ok; ii++)
  {
    Body *body_ptr = body_list.get_and_step();
    // Extract the BODY from Body
    OCCBody* fbody_ptr = CAST_TO(body_ptr->get_body_sm_ptr(), OCCBody);
    if (fbody_ptr)
    {
      delete_ok = fbody_ptr->can_be_deleted(body_list);
    }
  }
  return delete_ok;
}

int OCCQueryEngine::get_major_version()
{
  return OCCQE_MAJOR_VERSION;
}

int OCCQueryEngine::get_minor_version()
{
  return OCCQE_MINOR_VERSION;
}

int OCCQueryEngine::get_subminor_version()
{
  return OCCQE_SUBMINOR_VERSION;
}

CubitString OCCQueryEngine::get_engine_version_string()
{
  CubitString version_string = "OCC Geometry Engine version ";
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
const type_info& OCCQueryEngine::entity_type_info() const
{
   return typeid(*this);
}

//================================================================================
// Description:
// Author     :
// Date       :
//================================================================================
Point* OCCQueryEngine::make_Point( GeometryType ,
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
Curve* OCCQueryEngine::make_Curve(Curve *) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return (Curve*) NULL;
}

//================================================================================
// Description:
// Author     :
// Date       :
//================================================================================
Curve* OCCQueryEngine::make_Curve( Point const* ,
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
Curve* OCCQueryEngine::make_Curve( GeometryType ,
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
Curve* OCCQueryEngine::make_Curve( GeometryType ,
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
Surface* OCCQueryEngine::make_Surface( Surface *,
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
Surface* OCCQueryEngine::make_Surface( GeometryType ,
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
Lump* OCCQueryEngine::make_Lump( GeometryType ,
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
BodySM* OCCQueryEngine::make_BodySM( Surface * ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return (BodySM*) NULL;
}

//================================================================================
// Description:
// Author     :
// Date       :
//================================================================================
BodySM* OCCQueryEngine::make_BodySM( DLIList<Lump*>&  ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return (BodySM*) NULL;
}

//================================================================================
// Description: create a new body by copying
// Author     : sjowen
// Date       : 9/7/01
//================================================================================
Body* OCCQueryEngine::copy_body( Body *body_ptr )
{
  BodySM* bodysm_ptr = body_ptr->get_body_sm_ptr();
  OCCBody *facet_body_ptr = CAST_TO(bodysm_ptr, OCCBody);
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
CubitStatus OCCQueryEngine::reflect( BodySM *bodysm,
                                       const CubitVector& axis)
{
  OCCBody *fbody = CAST_TO(bodysm, OCCBody);
  if (!fbody)
  {
    PRINT_ERROR("Attempt to reflect mesh-based geometry Body.  This body is not MBG.");
    return CUBIT_FAILURE;
  }

  fbody->reflect( axis.x(), axis.y(), axis.z() );

  return CUBIT_SUCCESS;
}


//================================================================================
// Description:
// Author     :
// Date       :
//================================================================================
CubitStatus OCCQueryEngine::get_graphics( Surface* surface_ptr,
                                                     int& number_triangles,
                                                     int& number_points,
                                                     int& number_facets,
                                                     GMem* gMem,
                                                     unsigned short ,
                                                     double,
                                                     double ) const
{
  return CUBIT_SUCCESS;
}

//================================================================================
// Description:
// Author     :
// Date       :
//================================================================================
CubitStatus OCCQueryEngine::get_graphics( Curve* curve_ptr,
                                            int& num_points,
                                            GMem* gMem,
                                            double /*tolerance*/ ) const
{
  //  get the OCCCurve.
  //OCCCurve *facet_curv_ptr = CAST_TO(curve_ptr,OCCCurve);

  return CUBIT_SUCCESS;
}

//================================================================================
// Description:
// Author     :
// Date       :
//================================================================================
CubitStatus OCCQueryEngine::get_isoparametric_points(Surface* ,
                                                          int &nu, int &nv,
                                                          GMem*&) const
{
  nu = nv = 0;
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

CubitStatus OCCQueryEngine::get_u_isoparametric_points(Surface* ,
                                                            double, int&,
                                                            GMem*&) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

CubitStatus OCCQueryEngine::get_v_isoparametric_points(Surface* ,
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
CubitStatus OCCQueryEngine::transform_vec_position( CubitVector const& ,
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
CubitStatus OCCQueryEngine::get_intersections( Curve* , CubitVector& ,
                                                 CubitVector&,
                                                    DLIList<CubitVector*>& ,
                                                    CubitBoolean,
                                                    CubitBoolean )
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

CubitStatus OCCQueryEngine::get_intersections( Curve* , Curve* ,
                                                    DLIList<CubitVector*>& ,
                                                    CubitBoolean,
                                                    CubitBoolean )
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}
CubitStatus OCCQueryEngine::get_intersections( Curve*,
                                                   Curve*,
                                                   DLIList<CubitVector*>&,
                                                   double, CubitBoolean )
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

CubitStatus
OCCQueryEngine::get_intersections( Curve* /*ref_edge*/, Surface* /*ref_face*/,
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
OCCQueryEngine::entity_extrema( DLIList<GeometryEntity*> &,
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
OCCQueryEngine::entity_entity_distance( GeometryEntity *,
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

CubitStatus OCCQueryEngine::save_temp_geom_file( DLIList<TopologyBridge*>& ref_entity_list,
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

CubitStatus OCCQueryEngine::export_solid_model( DLIList<TopologyBridge*>& ref_entity_list,
                                                     const char* file_name,
                                                     const char* file_type,
                                                     const CubitString &,
                                                     const char*)
{
  if( strcmp( file_type, "FACET" ) )
    return CUBIT_SUCCESS;

  DLIList<OCCBody*>    facet_bodies;
  DLIList<OCCLump*>    facet_lumps;
  DLIList<OCCShell*>   facet_shells;
  DLIList<OCCSurface*> facet_surfaces;
  DLIList<OCCLoop*>    facet_loops;
  DLIList<OCCCoEdge*>  facet_coedges;
  DLIList<OCCCurve*>   facet_curves;
  DLIList<OCCPoint*>   facet_points;

  DLIList<TopologyBridge*> ref_entities_handled;

  int i;
  //Collect all free entities (bodies, curves, vertices )
  ref_entity_list.reset();
  for(i=ref_entity_list.size(); i>0; i--)
  {
    TopologyBridge* ref_entity_ptr = ref_entity_list.get();
    CubitBoolean handled = CUBIT_TRUE;

    //if it is a Vertex
    if( OCCPoint* pt = CAST_TO( ref_entity_ptr, OCCPoint) )
      facet_points.append( pt );
    //if it is a Curve
    else if( OCCCurve* curve = CAST_TO( ref_entity_ptr, OCCCurve) )
      facet_curves.append( curve );
    /*
    //if it is a Surface -- I don't think you can ever have a free surface
    //without it being a Body
    else if( OCCSurface* surf = CAST_TO( ref_entity_ptr, OCCSurface) )
      facet_surfaces.append( surf );
   */
    //if it is a Body
    else if( OCCBody* body = CAST_TO( ref_entity_ptr, OCCBody ) )
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

CubitStatus
OCCQueryEngine::gather_all_facet_entities( DLIList<OCCBody*> &facet_bodies,
                                             DLIList<OCCLump*> &facet_lumps,
                                             DLIList<OCCShell*> &facet_shells,
                                             DLIList<OCCSurface*> &facet_surfaces,
                                             DLIList<OCCLoop*> &facet_loops,
                                             DLIList<OCCCoEdge*> &facet_coedges,
                                             DLIList<OCCCurve*> &facet_curves,
                                             DLIList<OCCPoint*> &facet_points )
{
  int i;


  //Collect OCCLumps from OCCBody
  for(i=0; i<facet_bodies.size(); i++)
  {
    OCCBody *facet_body = facet_bodies.get_and_step();
    facet_body->get_lumps( facet_lumps );
  }


  //Collect OCCShells from OCCLumps
  for(i=0; i<facet_lumps.size(); i++)
  {
    OCCLump *facet_lump = facet_lumps.get_and_step();
    facet_lump->get_shells( facet_shells );
  }


  //Collect OCCSurfaces from OCCShells
  for(i=0; i<facet_shells.size(); i++)
  {
    OCCShell *facet_shell = facet_shells.get_and_step();
    facet_shell->get_surfaces( facet_surfaces );
  }
  facet_surfaces.uniquify_unordered();


  //Collect OCCLoops from OCCSurfaces
  for(i=0; i<facet_surfaces.size(); i++)
  {
    OCCSurface *facet_surface = facet_surfaces.get_and_step();
    facet_surface->get_loops( facet_loops );
  }

  //Collect OCCCoEdges from OCCLoops
  for(i=0; i<facet_loops.size(); i++)
  {
    OCCLoop *facet_loop = facet_loops.get_and_step();
    facet_loop->get_coedges( facet_coedges );
  }


  //Collect FacetCurves from OCCCoEdges
  for( i=0; i<facet_coedges.size(); i++)
  {
    OCCCoEdge *facet_coedge = facet_coedges.get_and_step();
    facet_coedge->get_curves( facet_curves );
  }
  facet_curves.uniquify_unordered();


  //Collect OCCPoints from FacetCurves
  for( i=0; i<facet_curves.size(); i++)
  {
    OCCCurve *facet_curve = facet_curves.get_and_step();
    facet_curve->get_points( facet_points );
  }

  //uniquify lists
  facet_points.uniquify_unordered();


  return CUBIT_SUCCESS;

}


CubitStatus
OCCQueryEngine::write_topology( FILE *file_ptr,
                                  DLIList<OCCBody*> &facet_bodies,
                                  DLIList<OCCLump*> &facet_lumps,
                                  DLIList<OCCShell*> &facet_shells,
                                  DLIList<OCCSurface*> &facet_surfaces,
                                  DLIList<OCCLoop*> &facet_loops,
                                  DLIList<OCCCoEdge*> &facet_coedges,
                                  DLIList<OCCCurve*> &facet_curves,
                                  DLIList<OCCPoint*> &facet_points )
{

  int i;

  //create a wrapper object for writing
  CIOWrapper file_writer( file_ptr );

  //-----------------write OCCPoints--------------
  UnsignedInt32 size = facet_points.size();
  //write out number of OCCPoints
  file_writer.Write( &size, 1 );
  facet_points.reset();
  for( i=0; i<facet_points.size(); i++)
  {
    OCCPoint *curr_point = facet_points.get_and_step();

    //file_writer.Write( reinterpret_cast<UnsignedInt32*>(&id), 1 );
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
    OCCCurve *curr_curve = facet_curves.get_and_step();
    Point *s_point, *e_point;
    s_point = curr_curve->start_point();
    e_point = curr_curve->end_point();

    int data_to_write[4];

    // get start&end points implicit ids
    OCCPoint *temp_point = NULL;
    temp_point = CAST_TO( s_point, OCCPoint );
    if( !temp_point ) assert(0);
    int found;
    found = facet_points.where_is_item( temp_point );
    if( found == -1)
      assert(0);
    data_to_write[0] = found;

    temp_point = CAST_TO( e_point, OCCPoint );
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
    OCCCoEdge *curr_coedge = facet_coedges.get_and_step();
    Curve *curve_sm;
    curve_sm = curr_coedge->curve();

    OCCCurve *temp_curve = NULL;
    temp_curve = CAST_TO( curve_sm, OCCCurve );

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

  //-----------------write OCCLoops--------------
  size = facet_loops.size();
  // write out number of OCCLoops
  file_writer.Write( &size, 1 );
  facet_loops.reset();
  for( i=0; i<facet_loops.size(); i++)
  {
    OCCLoop *curr_loop = facet_loops.get_and_step();
    DLIList<OCCCoEdge*> coedge_list;
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
      OCCCoEdge *temp_coedge = coedge_list.get_and_step();
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

  //-----------------write OCCSurfaces--------------
  size = facet_surfaces.size();
  // write out number of OCCSurfaces
  file_writer.Write( &size, 1 );
  facet_surfaces.reset();
  for( i=0; i<facet_surfaces.size(); i++)
  {
    OCCSurface *curr_surface = facet_surfaces.get_and_step();

    DLIList<OCCLoop*> loop_list;
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
    //data_to_write[2] = curr_surface->get_eval_tool()->get_output_id();

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
       OCCLoop *temp_loop = loop_list.get_and_step();
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

  //-----------------write OCCShells--------------
  size = facet_shells.size();
  // write out number of OCCShells
  file_writer.Write( &size, 1 );
  facet_shells.reset();
  for( i=0; i<facet_shells.size(); i++)
  {
    OCCShell *curr_shell= facet_shells.get_and_step(); //number of surfaces
    DLIList<OCCSurface*> temp_facet_surf_list;
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
      OCCSurface *temp_facet_surface = temp_facet_surf_list.get_and_step();
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

  //-----------------write OCCLumps--------------
  size = facet_lumps.size();
  // write out number of OCCLumps
  file_writer.Write( &size, 1 );
  facet_lumps.reset();
  for( i=0; i<facet_lumps.size(); i++)
  {
    OCCLump *curr_lump = facet_lumps.get_and_step();

    DLIList<OCCShell*> temp_facet_shell_list;
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
      OCCShell *temp_facet_shell = temp_facet_shell_list.get_and_step();
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
    OCCBody *curr_body = facet_bodies.get_and_step();

    DLIList<OCCLump*> temp_facet_lump_list;
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
      OCCLump *temp_facet_lump = temp_facet_lump_list.get_and_step();
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
OCCQueryEngine::restore_topology( FILE *file_ptr,
                                    unsigned int endian,
                                    int num_points,
                                    CubitPoint **points_array,
                                    int num_cfet,
                                    CurveFacetEvalTool** cfev_array,
                                    int num_fet,
                                    FacetEvalTool** fev_array,
                                    DLIList<TopologyBridge*> &imported_entities )
{

  return CUBIT_SUCCESS;
}


CubitStatus
OCCQueryEngine::import_temp_geom_file(FILE* file_ptr,
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

CubitStatus OCCQueryEngine::import_solid_model(
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
  TopoDS_Shape *aShape = new TopoDS_Shape;
  BRep_Builder aBuilder;
  Standard_Boolean result = BRepTools::Read(*aShape, (char*) file_name, aBuilder);
  if (result==0) return CUBIT_FAILURE;
  OCCBody *body = new OCCBody(aShape);
//  CGMList->append(body);
  imported_entities.append(body);
//  OCCMap->Bind(*aShape, CGMList->where_is_item(body));
  populate_topology_bridge_solid(*aShape, imported_entities);
  return CUBIT_SUCCESS;
}

CubitStatus OCCQueryEngine::populate_topology_bridge_solid(TopoDS_Shape aShape, DLIList<TopologyBridge*> &imported_entities)
{
	TopExp_Explorer Ex;
	for (Ex.Init(aShape, TopAbs_SOLID); Ex.More(); Ex.Next()) {
		TopoDS_Solid *posolid =  new TopoDS_Solid;
 		*posolid = TopoDS::Solid(Ex.Current());
		OCCLump *lump;
		if (!OCCMap->IsBound(*posolid)) {
			printf("Adding solid\n");
			lump = new OCCLump(posolid);
			CGMList->append(lump);
//			imported_entities.append(lump);
			OCCMap->Bind(*posolid, CGMList->where_is_item(lump));
			populate_topology_bridge_shell(*posolid, imported_entities);
		} else {
			lump = (OCCLump*)((*CGMList)[OCCMap->Find(*posolid)]);
		}
	}
	return CUBIT_SUCCESS;
}

CubitStatus OCCQueryEngine::populate_topology_bridge_shell(TopoDS_Shape aShape, DLIList<TopologyBridge*> &imported_entities)
{
	TopExp_Explorer Ex;
	for (Ex.Init(aShape, TopAbs_SHELL); Ex.More(); Ex.Next()) {
		TopoDS_Shell *poshell = new TopoDS_Shell;
		*poshell = TopoDS::Shell(Ex.Current());
		OCCShell *shell;
		if (!OCCMap->IsBound(*poshell)) {
			printf("Adding shell\n");
			shell = new OCCShell(poshell);
			CGMList->append(shell);
//			imported_entities.append(shell);
			OCCMap->Bind(*poshell, CGMList->where_is_item(shell));
			populate_topology_bridge_face(*poshell, imported_entities);
		} else {
			shell = (OCCShell*)(*CGMList)[OCCMap->Find(*poshell)];
		}
		shell->add_lump((OCCLump*)(*CGMList)[OCCMap->Find(aShape)]);
	}
	return CUBIT_SUCCESS;
}

CubitStatus OCCQueryEngine::populate_topology_bridge_face(TopoDS_Shape aShape, DLIList<TopologyBridge*> &imported_entities)
{
	TopExp_Explorer Ex;
	for (Ex.Init(aShape, TopAbs_FACE); Ex.More(); Ex.Next()) {
		TopoDS_Face *poface = new TopoDS_Face;
		*poface = TopoDS::Face(Ex.Current());
		OCCSurface *surface;
		if (!OCCMap->IsBound(*poface)) {
			printf("Adding face\n");
/*			BRepAdaptor_Surface asurface(*poface);
			Bnd_Box aBox;
			double min[3], max[3];
			BndLib_AddSurface::Add(asurface, Precision::Approximation(), aBox);
			aBox.Get( min[0], min[1], min[2], max[0], max[1], max[2]);
			printf(".  .  .  box: %lf %lf %lf, %lf %lf %lf\n", min[0], min[1], min[2], max[0], max[1], max[2]);*/

			surface = new OCCSurface(poface);
			CGMList->append(surface);
//			imported_entities.append(surface);
			OCCMap->Bind(*poface, CGMList->where_is_item(surface));
			populate_topology_bridge_wire(*poface, imported_entities);
		} else {
			surface = (OCCSurface*)(*CGMList)[OCCMap->Find(*poface)];
		}
		surface->add_shell((OCCShell*)(*CGMList)[OCCMap->Find(aShape)]);
	}
	return CUBIT_SUCCESS;
}

CubitStatus OCCQueryEngine::populate_topology_bridge_wire(TopoDS_Shape aShape, DLIList<TopologyBridge*> &imported_entities)
{
	TopExp_Explorer Ex;
	for (Ex.Init(aShape, TopAbs_WIRE); Ex.More(); Ex.Next()) {
		TopoDS_Wire *powire = new TopoDS_Wire;
		*powire = TopoDS::Wire(Ex.Current());
		OCCLoop *loop;
		if (!OCCMap->IsBound(*powire)) {
			printf("Adding wire\n");
			loop = new OCCLoop(powire);
			CGMList->append(loop);
//			imported_entities.append(loop);
			OCCMap->Bind(*powire, CGMList->where_is_item(loop));
			populate_topology_bridge_edge(*powire, imported_entities);
		} else {
			loop = (OCCLoop*)(*CGMList)[OCCMap->Find(*powire)];
		}
		loop->add_surface((OCCSurface*)(*CGMList)[OCCMap->Find(aShape)]);
	}
	return CUBIT_SUCCESS;
}

CubitStatus OCCQueryEngine::populate_topology_bridge_edge(TopoDS_Shape aShape, DLIList<TopologyBridge*> &imported_entities)
{
	TopExp_Explorer Ex;
	TopTools_DataMapOfShapeInteger *Map, *Mapr;
	for (Ex.Init(aShape, TopAbs_EDGE); Ex.More(); Ex.Next()) {
		TopoDS_Edge *poedge = new TopoDS_Edge;
		*poedge = TopoDS::Edge(Ex.Current());
		Curve *curve;
		OCCCoEdge *coedge;
		TopoDS_Shape test = *poedge;
		if (poedge->Orientation() == TopAbs_FORWARD) {
			Map = OCCMap;
			Mapr = OCCMapr;
		} else if (poedge->Orientation() == TopAbs_REVERSED) {
			Map = OCCMapr;
			Mapr = OCCMap;
		} else {
			printf("Oooooop!\n");
		}
		if (!Mapr->IsBound(*poedge)) {
			printf("Adding edge\n");
			curve = new OCCCurve(poedge);
		} else {
			curve = ((OCCCoEdge*)(*CGMList)[Mapr->Find(*poedge)])->curve();
		}
		if (!Map->IsBound(*poedge)) {
			printf("Adding coedge\n");
			coedge = new OCCCoEdge(poedge, curve);
			CGMList->append(coedge);
//			CGMList->append(curve);
			Map->Bind(*poedge, CGMList->where_is_item(coedge));
			populate_topology_bridge_vertex(*poedge, imported_entities, curve);
		} else {
			coedge = (OCCCoEdge*)(*CGMList)[Map->Find(*poedge)];
		}
		coedge->add_loop((OCCLoop*)occ_to_cgm(aShape));
	}
	return CUBIT_SUCCESS;
}

CubitStatus OCCQueryEngine::populate_topology_bridge_vertex(TopoDS_Shape aShape, DLIList<TopologyBridge*> &imported_entities, Curve *curve)
{
	TopExp_Explorer Ex;
	for (Ex.Init(aShape, TopAbs_VERTEX); Ex.More(); Ex.Next()) {
		TopoDS_Vertex *povertex = new TopoDS_Vertex;
		*povertex = TopoDS::Vertex(Ex.Current());
		OCCPoint *point;
		if (!OCCMap->IsBound(*povertex)) {
			printf("Adding vertex\n");
                        gp_Pnt pt = BRep_Tool::Pnt(*povertex);
			point = new OCCPoint(pt);
			CGMList->append(point);
//			imported_entities.append(point);
			OCCMap->Bind(*povertex, CGMList->where_is_item(point));
		} else {
			point = (OCCPoint*)(*CGMList)[OCCMap->Find(*povertex)];
		}
	}
	return CUBIT_SUCCESS;
}

TopologyBridge* OCCQueryEngine::occ_to_cgm(TopoDS_Shape shape)
{
	
	if ((shape.ShapeType() != TopAbs_EDGE) || (shape.Orientation() == TopAbs_FORWARD)) return (*CGMList)[OCCMap->Find(shape)];
	else return (*CGMList)[OCCMapr->Find(shape)];
/*	if (OCCMap->IsBound(shape))*/ return (*CGMList)[OCCMap->Find(shape)];
//	else return (*CGMList)[(*OCCMapt)[shape]];
}

CubitStatus OCCQueryEngine::import_solid_model(FILE *file_ptr,
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


//===============================================================================
// Function   : restore_facets
// Member Type: PUBLIC
// Description: restore facets and eval tools onto the list of entities
// Author     : sjowen
// Date       : 1/26/03
//===============================================================================
CubitStatus OCCQueryEngine::restore_facets(
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

  return rv;

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
void OCCQueryEngine::delete_solid_model_entities(DLIList<BodySM*>&BodyList) const
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
// Purpose       : Delete a OCCBody and child entities.
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/29/03
//-------------------------------------------------------------------------
CubitStatus
OCCQueryEngine::delete_solid_model_entities( BodySM* bodysm ) const
{
  OCCBody* fbody = dynamic_cast<OCCBody*>(bodysm);
  if (!fbody)
    return CUBIT_FAILURE;

  DLIList<OCCLump*> lumps;
  DLIList<OCCShell*> shells;
  DLIList<OCCSurface*> surfaces;

  fbody->get_lumps(lumps);
  fbody->disconnect_all_lumps();
  delete fbody;

  for (int i = lumps.size(); i--; )
  {
    OCCLump* lump = lumps.get_and_step();

    shells.clean_out();
    lump->get_shells(shells);
    lump->disconnect_all_shells();
    delete lump;

    for (int j = shells.size(); j--; )
    {
      OCCShell* shell = shells.get_and_step();

      surfaces.clean_out();
      shell->get_surfaces(surfaces);
      shell->disconnect_all_surfaces();
      delete shell;

      for (int k = surfaces.size(); k--; )
      {
        OCCSurface* surface = surfaces.get_and_step();
        if (!surface->has_parent_shell())
          delete_solid_model_entities(surface);
      }
    }
  }

  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Delete a OCCSurface and child entities.
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/29/03
//-------------------------------------------------------------------------
CubitStatus
OCCQueryEngine::delete_solid_model_entities( Surface* surface ) const
{
  OCCSurface* fsurf = dynamic_cast<OCCSurface*>(surface);
  if (!fsurf || fsurf->has_parent_shell())
    return CUBIT_FAILURE;

  DLIList<OCCLoop*> loops;
  DLIList<OCCCoEdge*> coedges;

  fsurf->get_loops(loops);
  fsurf->disconnect_all_loops();
  delete fsurf;

  for (int i = loops.size(); i--; )
  {
    OCCLoop* loop = loops.get_and_step();

    coedges.clean_out();
    loop->get_coedges(coedges);
    loop->disconnect_all_coedges();
    delete loop;

    for (int j = coedges.size(); j--; )
    {
      OCCCoEdge* coedge = coedges.get_and_step();
      OCCCurve* curve = dynamic_cast<OCCCurve*>(coedge->curve());
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
// Purpose       : Delete a OCCCurve and child entities.
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/29/03
//-------------------------------------------------------------------------
CubitStatus
OCCQueryEngine::delete_solid_model_entities( Curve* curve ) const
{
  OCCCurve* fcurve = dynamic_cast<OCCCurve*>(curve);
  if (!fcurve || fcurve->has_parent_coedge())
    return CUBIT_FAILURE;

  OCCPoint* start = dynamic_cast<OCCPoint*>(fcurve->start_point());
  OCCPoint*   end = dynamic_cast<OCCPoint*>(fcurve->end_point()  );

  if (start == end )
      end = NULL;

  if (start)
  {
    //start->disconnect_curve(fcurve);
    //if (!start->has_parent_curve())
    //  delete_solid_model_entities(start);
  }

  if (end)
  {
    //end->disconnect_curve(fcurve);
    //if (!end->has_parent_curve())
    //  delete_solid_model_entities(end);
  }

  delete curve;
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Delete a OCCPoint and child entities.
//
// Special Notes :
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/29/03
//-------------------------------------------------------------------------
CubitStatus
OCCQueryEngine::delete_solid_model_entities( Point* point ) const
{
  //OCCPoint* fpoint = dynamic_cast<OCCPoint*>(point);

  delete point;
  return CUBIT_SUCCESS;
}

CubitStatus OCCQueryEngine::fire_ray(BodySM *,
                                          const CubitVector &,
                                          const CubitVector &,
                                          DLIList<double>&,
                                          DLIList<GeometryEntity*> *) const
{
  PRINT_ERROR("OCCQueryEngine::fire_ray not yet implemented.\n");
  return CUBIT_FAILURE;
}
  //- fire a ray at the specified body, returning the entities hit and
  //- the parameters along the ray; return CUBIT_FAILURE if error

double OCCQueryEngine::get_sme_resabs_tolerance() const
{
  PRINT_ERROR("OCCQueryEngine::get_sme_resabs_tolerance not yet implemented.\n");
  return CUBIT_FAILURE;
}
// Gets solid modeler's resolution absolute tolerance

double OCCQueryEngine::set_sme_resabs_tolerance( double )
{
  PRINT_ERROR("OCCQueryEngine::set_sme_resabs_tolerance not yet implemented.\n");
  return CUBIT_FAILURE;
}

CubitStatus OCCQueryEngine::set_int_option( const char* , int )
{
  PRINT_ERROR("OCCQueryEngine::set_int_option not yet implemented.\n");
  return CUBIT_FAILURE;
}

CubitStatus OCCQueryEngine::set_dbl_option( const char* , double )
{
  PRINT_ERROR("OCCQueryEngine::set_dbl_option not yet implemented.\n");
  return CUBIT_FAILURE;
}

CubitStatus OCCQueryEngine::set_str_option( const char* , const char* )
{
  PRINT_ERROR("OCCQueryEngine::set_str_option not yet implemented.\n");
  return CUBIT_FAILURE;
}
  //- Set solid modeler options


//===========================================================================
//Function Name: ensure_is_ascii_stl_file
//Member Type:
//Description: returns CUBIT_TRUE in is_ascii if fp is an ascii stl file
//Author: Plamen Stoyanov (USF)
//===========================================================================
CubitStatus OCCQueryEngine::ensure_is_ascii_stl_file(FILE * fp, CubitBoolean &is_ascii)
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
CubitStatus OCCQueryEngine::read_facets_stl_tolerance(
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

  fclose( fp );
  return CUBIT_SUCCESS;
}


//===========================================================================
//Function Name: read_facets_stl
//Member Type:
//Description:  read facets from stl file
//Author: Plamen Stoyanov (USF)
//===========================================================================

CubitStatus OCCQueryEngine::read_facets_stl(
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
  fclose( fp );
  return CUBIT_SUCCESS;
}



//===========================================================================
//Function Name: read_facets
//Member Type:  PUBLIC
//Description:  read facets from facet file
//===========================================================================
CubitStatus OCCQueryEngine::read_facets( const char * file_name,
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

  fclose( fp );
  return CUBIT_SUCCESS;
}

//===========================================================================
//  Function: read_cholla_file
//  Purpose:  import the face mesh from the cholla file
//  Date:     11/28/2002
//  Author:   sjowen
//===========================================================================
CubitStatus OCCQueryEngine::read_cholla_file( const char *file_name,
                              double &feature_angle,
                              DLIList<CubitPoint *> &point_list,
                              DLIList<CubitFacet *> &facet_list )

{
  return CUBIT_SUCCESS;
}


//=============================================================================
// Function Name: export_facets
// Member Type:  PUBLIC
// Description:  export a list of facets to a facet file for debugging purposes
// Author: sjowen
// Date: 10/29/01
//=============================================================================
CubitStatus OCCQueryEngine::export_facets(DLIList<CubitFacet*> &facet_list,
                                               char *filename)
{

  // open the file for writing

  FILE *fp = fopen(filename, "w");
  if (!fp)
  {
    PRINT_ERROR("Couldn't open file %s for writing.\n", filename);
    return CUBIT_FAILURE;
  }

  fclose(fp);
  return CUBIT_SUCCESS;
}


//===========================================================================
//Function Name: save_facets
//Member Type:  PUBLIC
//Description:  save the facets from the entities to a Cubit file
//Author:       sjowen
//Date:         1/16/2003
//===========================================================================
CubitStatus OCCQueryEngine::save_facets(
  FILE *fp,  // file we are dumping to
  DLIList<OCCSurface*> facet_surfaces,
  DLIList<OCCCurve*> facet_curves,
  DLIList<OCCPoint*> facet_points ) // save facets from these entities
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
CubitStatus OCCQueryEngine::save_eval_tools(
  FILE *fp,        // cubit file we are dumping to
  DLIList<OCCSurface*> facet_surfaces,
  DLIList<OCCCurve*>  facet_curves )
{

  int ii;
  OCCSurface *fsurf_ptr;
  OCCCurve *fcurv_ptr;

  //separate out eval tools and set ids
  DLIList<CurveFacetEvalTool*> cfe_tools;
  DLIList<FacetEvalTool*> fe_tools;
  int ft_id = 0;
  int cft_id = 0;

  facet_surfaces.reset();
  for (ii=facet_surfaces.size(); ii--; )
  {
    fsurf_ptr = facet_surfaces.get_and_step();
    FacetEvalTool *feval_tool; //= fsurf_ptr->get_eval_tool();
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
CubitStatus OCCQueryEngine::restore_eval_tools(
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
CubitStatus OCCQueryEngine::dump_facets(
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


//=============================================================================
//Function:   create_super_facet_bounding_box(PUBLIC)
//Description: Find the bounding box of a list of BodySMs
//Author: jdfowle
//Date: 12/15/03
//=============================================================================
CubitStatus OCCQueryEngine::create_super_facet_bounding_box(
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
CubitStatus OCCQueryEngine::create_facet_bounding_box(
                                BodySM* bodySM,
                                CubitBox& bbox )
{
  OCCBody *fbody_ptr;
  int j;
  DLIList<OCCSurface*> facet_surf_list;
  DLIList<CubitFacet*> facet_list;
  DLIList<CubitPoint*> point_list;
  OCCSurface *facet_surface;
  CubitBox surf_bbox, total_box;
  CubitStatus status = CUBIT_FAILURE;

    fbody_ptr = dynamic_cast<OCCBody *>(bodySM);
    fbody_ptr->get_surfaces(facet_surf_list);
    for ( j = 0; j < facet_surf_list.size(); j++ ) {
      facet_surface = facet_surf_list.get_and_step();
      facet_list.clean_out();
      point_list.clean_out();
      //facet_surface->get_my_facets(facet_list,point_list);
      status = FacetDataUtil::get_bbox_of_points(point_list,surf_bbox);
      if ( j == 0 ) total_box = surf_bbox;
      else 
        total_box |= surf_bbox;
    }
  bbox |= total_box;  
  return status;
}

const char* fqe_xform_err = "Transform not implemented for facet geometry.\n";
CubitStatus OCCQueryEngine::restore_transform( BodySM* body )
{
  OCCBody* facetbod = dynamic_cast<OCCBody*>(body);
  return facetbod ? facetbod->restore( ) : CUBIT_FAILURE;
}
CubitStatus OCCQueryEngine::translate( BodySM* body, const CubitVector& d )
{
  OCCBody* facetbod = dynamic_cast<OCCBody*>(body);
  return facetbod ? facetbod->move( d.x(), d.y(), d.z() ) : CUBIT_FAILURE;
}
CubitStatus OCCQueryEngine::rotate( BodySM* body, const CubitVector& v, double a )
{
  OCCBody* facetbod = dynamic_cast<OCCBody*>(body);
  return facetbod ? facetbod->rotate( v.x(), v.y(), v.z(), a ) : CUBIT_FAILURE;
}
CubitStatus OCCQueryEngine::scale( BodySM* body, double factor )
{
  OCCBody* facetbod = dynamic_cast<OCCBody*>(body);
  return facetbod ? facetbod->scale( factor ) : CUBIT_FAILURE;
}
CubitStatus OCCQueryEngine::scale( BodySM* body, const CubitVector& f )
{
  OCCBody* facetbod = dynamic_cast<OCCBody*>(body);
  return facetbod ? facetbod->scale( f.x(), f.y(), f.z() ) : CUBIT_FAILURE;
}

CubitStatus OCCQueryEngine::translate( GeometryEntity* , const CubitVector&  )
{
  PRINT_ERROR(fqe_xform_err);
  return CUBIT_FAILURE;
}
CubitStatus OCCQueryEngine::rotate( GeometryEntity* , const CubitVector& , double  )
{
  PRINT_ERROR(fqe_xform_err);
  return CUBIT_FAILURE;
}
CubitStatus OCCQueryEngine::scale( GeometryEntity* , double  )
{
  PRINT_ERROR(fqe_xform_err);
  return CUBIT_FAILURE;
}
CubitStatus OCCQueryEngine::scale( GeometryEntity* , const CubitVector&  )
{
  PRINT_ERROR(fqe_xform_err);
  return CUBIT_FAILURE;
}
CubitStatus OCCQueryEngine::reflect( GeometryEntity* , const CubitVector&  )
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
CubitBoolean OCCQueryEngine::bodies_overlap (BodySM * body_ptr_1,
                                                BodySM * body_ptr_2 ) const
{

   return CUBIT_TRUE;
}

CubitBoolean OCCQueryEngine::volumes_overlap (Lump *lump1, Lump *lump2 ) const
{
}

//EOF
