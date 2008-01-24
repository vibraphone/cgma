//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

// Filename      : OCCModifyEngine.cpp
//
// Purpose       : ModifyEngine for OCC geometry
//
// Special Notes : Modeled after GeometryModifyEngine and AcisModifyEngine.
//
// Author        : Jane Hu
//
// Creation Date : 1/08
//
//-------------------------------------------------------------------------
#include "config.h"
#include "gp_Pnt.hxx"
#include "TopoDS_Shape.hxx"
#include "TopoDS.hxx"
#include "TopologyBridge.hpp"
#include "OCCModifyEngine.hpp"
#include "OCCQueryEngine.hpp"
#include "CubitMessage.hpp"
#include "CubitDefines.h"

#include "CubitUtil.hpp"
#include "CubitPoint.hpp"
#include "CubitPointData.hpp"
#include "CubitFacet.hpp"
#include "CubitQuadFacet.hpp"
#include "CubitFacetData.hpp"
#include "CubitFacetEdge.hpp"
#include "CubitFacetEdgeData.hpp"
#include "GeometryQueryTool.hpp"
#include "GeometryModifyTool.hpp"
#include "ChollaSurface.hpp"
#include "ChollaCurve.hpp"
#include "ChollaPoint.hpp"
#include "OCCPoint.hpp"
#include "OCCCurve.hpp"
#include "CurveFacetEvalTool.hpp"
#include "FacetEvalTool.hpp"
#include "OCCCoEdge.hpp"
#include "OCCLoop.hpp"
#include "OCCSurface.hpp"
#include "OCCShell.hpp"
#include "OCCLump.hpp"
#include "OCCBody.hpp"
#include "ChollaEngine.hpp"
#include "TDGeomFacet.hpp"
#include "CubitFileIOWrapper.hpp"
#include "TDFacetBoundaryPoint.hpp"
#include "Cholla.h"
#include "Body.hpp"
#include "GfxDebug.hpp"
#include "RefFace.hpp"
#include "FacetDataUtil.hpp"
#include "FBDataUtil.hpp"
#include "FBIntersect.hpp"
#include "IntegerHash.hpp"
#include "CpuTimer.hpp"
#include "AppUtil.hpp"
#include "SphereEvaluator.hpp"
#include "CylinderEvaluator.hpp"

OCCModifyEngine* OCCModifyEngine::instance_ = 0;

//===============================================================================
// Function   : OCCModifyEngine
// Member Type: PUBLIC
// Description: constructor
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
OCCModifyEngine::OCCModifyEngine()
{
//  assert( !instance_ );

    // add this modify engine to geometrymodifytool
  GeometryModifyTool::instance()->add_gme(this);
}


//===============================================================================
// Function   : ~OCCModifyEngine
// Member Type: PUBLIC
// Description: destructor
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
OCCModifyEngine::~OCCModifyEngine() 
{
        instance_ = 0;
}

//===============================================================================
// Function   : make_Point
// Member Type: PUBLIC
// Description: make a geometric entity point
// Author     : Jane Hu 
// Date       : 10/07
//===============================================================================
Point* OCCModifyEngine::make_Point( CubitVector const& point) const
{
  gp_Pnt pt = gp_Pnt( point.x(), point.y(), point.z());
  TopoDS_Vertex theVertex = BRepBuilderAPI_MakeVertex(pt);

  // Create a new PointACIS object
  return OCCQueryEngine::instance()->populate_topology_bridge( theVertex );
}

//===============================================================================
// Function   : make_Curve
//              This function creates a curve given an existing curve, copy. 
// Member Type: PUBLIC
// Description: make a curve
// Author     : Jane Hu
// Date       : 01/08
//===============================================================================
Curve* OCCModifyEngine::make_Curve(Curve * curve_ptr) const
{
  OCCCurve* occ_curve = CAST_TO(curve_ptr, OCCCurve);
  if (!occ_curve)
  {
     PRINT_ERROR("Cannot create an OCC curve from the given curve.\n"
                 "Possible incompatible geometry engines.\n");
     return (Curve *)NULL;
  }
 
  TopoDS_Edge *theEdge = occ_curve->get_TopoDS_Edge();  
 
  TopoDS_Shape newShape = theEdge->EmptyCopied();
 
  TopoDS_Edge newEdge = TopoDS::Edge(newShape);

  return OCCQueryEngine::instance()->populate_topology_bridge(newEdge);
}

//===============================================================================
// Function   : make_Curve
// Member Type: PUBLIC
// Description: make a curve by projecting a straight line defined by 
//              point1_ptr, and point2_ptr onto face_ptr, third_point
//              is used for curves that could be periodic to dertermine
//              the correct direction.
// Author     : Jane Hu
// Date       : 01/08
//===============================================================================
Curve* OCCModifyEngine::make_Curve( Point const* point1_ptr,
                             Point const* point2_ptr,
                             Surface* face_ptr,
                             const CubitVector * third_point) const
{
  assert (point1_ptr != NULL && point2_ptr != NULL);
  GeometryType type = STRAIGHT_CURVE_TYPE;
  CubitBoolean closed = CUBIT_FALSE;
  DLIList<CubitVector*> mid_points;
  Curve* curve = NULL;
  if (point1_ptr != point2_ptr)
    curve = make_Curve(type, point1_ptr, point2_ptr, mid_points);
  else //could be a closed shape
  {
    if(third_point != NULL && face_ptr != NULL) 
    {
       closed = CUBIT_TRUE;
       //curve = make_Curve(type, point1_ptr, third_point, mid_points);
    }
  }

  Curve* new_curve = NULL;
  if(face_ptr == NULL)
    return curve;
 
  new_curve = 
	CAST_TO(curve, OCCCurve)->project_curve(face_ptr, closed, third_point); 

  delete curve;
  return new_curve;
}

//===============================================================================
// Function   : make_Curve
// Member Type: PUBLIC
// Description: make a curve
// Author     : Jane Hu
// Date       : 01/08
//===============================================================================
Curve* OCCModifyEngine::make_Curve( GeometryType curve_type,
                             Point const* point1_ptr,
                             Point const* point2_ptr,
                             DLIList<CubitVector*>& vector_list,
                             Surface* face_ptr) const
{
  assert(point1_ptr != NULL && point2_ptr != NULL);
 /* 
  OCCPoint const* occ_point1 = CAST_TO(point1_ptr, OCCPoint);
  OCCPoint const* occ_point2 = CAST_TO(point2_ptr, OCCPoint);

  if (occ_point1 == NULL || occ_point2 == NULL)
  {
     PRINT_ERROR("Cannot create an OCC curve from the given points.\n"
                 "Possible incompatible geometry engines.\n");
     return (Curve *)NULL;
  }
 */   
  //project all points on the surface
  if (face_ptr != NULL)
  {
  }    
     
  return (Curve*) NULL;
}

//===============================================================================
// Function   : make_Curve
// Member Type: PUBLIC
// Description: make a curve
// Author     : Jane Hu 
// Date       : 01/08
//===============================================================================
Curve* OCCModifyEngine::make_Curve( GeometryType curve_type,
                             Point const* point1_ptr,
                             Point const* point2_ptr,
                             CubitVector const* intermediate_point_ptr,
                             CubitSense sense) const
{
  assert (point1_ptr != NULL && point2_ptr != NULL);
  DLIList<CubitVector*> mid_points;
  CubitVector mid_point = *intermediate_point_ptr;
  if (intermediate_point_ptr != NULL )
    mid_points.append(&mid_point);

  Point const* tmp_point = point1_ptr;
  if (sense == CUBIT_REVERSED)
  {
     point1_ptr = point2_ptr;
     point2_ptr = tmp_point;
  }   
  return make_Curve(curve_type, point1_ptr, point2_ptr, mid_points);
}

//===============================================================================
// Function   : make_Surface
// Member Type: PUBLIC
// Description: make a surface
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
Surface* OCCModifyEngine::make_Surface( Surface * /*old_surface_ptr*/,
                                 CubitBoolean /*extended_from*/) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return (Surface*) NULL;
}

//===============================================================================
// Function   : make_Surface
// Member Type: PUBLIC
// Description: make a surface
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
Surface* OCCModifyEngine::make_Surface( GeometryType /*surface_type*/,
                                 DLIList<Curve*>& /*curve_list*/,
                                 Surface * /*old_surface_ptr*/,
                                 bool /*check_edges*/) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return (Surface*) NULL;
}

//===============================================================================
// Function   : make_Lump
// Member Type: PUBLIC
// Description: make a lump
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
Lump* OCCModifyEngine::make_Lump( DLIList<Surface*>& /*surface_list*/ ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return (Lump*) NULL;
}

//===============================================================================
// Function   : make_BodySM
// Member Type: PUBLIC
// Description: make a BodySM
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
BodySM* OCCModifyEngine::make_BodySM( Surface * ) const
    {return NULL ;}

//===============================================================================
// Function   : make_BodySM
// Member Type: PUBLIC
// Description: make a BodySM
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
BodySM* OCCModifyEngine::make_BodySM( DLIList<Lump*>& /*lump_list*/ ) const
    {return NULL ;}


//===============================================================================
// Function   : sphere
// Member Type: PUBLIC
// Description: build a sphere with facets
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
BodySM* OCCModifyEngine::sphere(double radius) const
{
  return (BodySM*) NULL;
}


//===============================================================================
// Function   : brick
// Member Type: PUBLIC
// Description: build a brick with facets
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
BodySM* OCCModifyEngine::brick( double wid, double dep, double hi ) const
{
  return (BodySM*)NULL;
}


//===============================================================================
// Function   : brick
// Member Type: PUBLIC
// Description: create a brick with facets given center axes and extension
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
BodySM* OCCModifyEngine::brick( const CubitVector &/*center*/, 
                                  const CubitVector* /*axes[3]*/,
                                  const CubitVector &/*extension*/) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return (BodySM*) NULL;
}

//===============================================================================
// Function   : prism
// Member Type: PUBLIC
// Description: create a prism with facets
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
BodySM* OCCModifyEngine::prism( double /*height*/, int /*sides*/, double /*major*/,
                               double /*minor*/) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return (BodySM*) NULL;
}

//===============================================================================
// Function   : pyramid
// Member Type: PUBLIC
// Description: create a pyramid with facets
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
BodySM* OCCModifyEngine::pyramid( double /*height*/, int /*sides*/, double /*major*/,
                                 double /*minor*/, double /*top*/) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return (BodySM*) NULL;
}

//===============================================================================
// Function   : cylinder
// Member Type: PUBLIC
// Description: create a cylinder with facets
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
BodySM* OCCModifyEngine::cylinder( double hi, double r1, double r2, double r3 ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return (BodySM*) NULL;
}

//===============================================================================
// Function   : torus
// Member Type: PUBLIC
// Description: create a torus with facets
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
BodySM* OCCModifyEngine::torus( double r1, double r2 ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return (BodySM*) NULL;
}

//===============================================================================
// Function   : planar_sheet
// Member Type: PUBLIC
// Description: create a planar_sheet with facets
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
BodySM* OCCModifyEngine::planar_sheet ( const CubitVector& /*p1*/,
                                       const CubitVector& /*p2*/,
                                       const CubitVector& /*p3*/,
                                       const CubitVector& /*p4*/ ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return (BodySM*) NULL;
}

//===============================================================================
// Function   : copy_body
// Member Type: PUBLIC
// Description: copy a facet-based body
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
BodySM* OCCModifyEngine::copy_body ( BodySM* /*bodyPtr*/) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return (BodySM*) NULL;
}

CubitStatus OCCModifyEngine::stitch_surfs(
                      DLIList<BodySM*>& surf_bodies,
                      BodySM*& stitched_body) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_SUCCESS;
}

//===============================================================================
// Function   : subtract
// Member Type: PUBLIC
// Description: subtract boolean operation on facet-based bodies
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus     OCCModifyEngine::subtract(DLIList<BodySM*> &tool_body_list,
                                            DLIList<BodySM*> &from_bodies,
                                            DLIList<BodySM*> &new_bodies,
                                            bool /*imprint*/,
                                            bool keep_old) const
{
  return CUBIT_SUCCESS; 
}

//===============================================================================
// Function   : imprint
// Member Type: PUBLIC
// Description: imprint boolean operation on facet-based bodies
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus     OCCModifyEngine::imprint(BodySM* /*BodyPtr1*/, BodySM* /*BodyPtr2*/,
                                           BodySM*& /*newBody1*/, BodySM*& /*newBody2*/,
                                           bool  /*keep_old*/) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : imprint
// Member Type: PUBLIC
// Description: imprint boolean operation on facet-based bodies
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus     OCCModifyEngine::imprint(DLIList<BodySM*> &from_body_list ,
                                           DLIList<BodySM*> &new_from_body_list,
                                           bool keep_old,
                                           DLIList<TopologyBridge*>* ,
                                           DLIList<TopologyBridge*>*) const
{
  CubitStatus success = CUBIT_SUCCESS;

  return success;
}

//===============================================================================
// Function   : imprint
// Member Type: PUBLIC
// Description: imprint boolean operation on facet-based bodies
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus     OCCModifyEngine::imprint( DLIList<BodySM*> &body_list,
                                           DLIList<Curve*> &ref_edge_list,
                                           DLIList<BodySM*>& new_body_list,
                                           bool keep_old,
                                           bool show_messages) const
{
  CubitStatus success = CUBIT_SUCCESS;
  return success;
}

//===============================================================================
// Function   : imprint
// Member Type: PUBLIC
// Description: imprint boolean operation on facet-based bodies
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus     OCCModifyEngine::imprint( DLIList<Surface*> &/*ref_face_list*/,
                                           DLIList<Curve*> &/*ref_edge_list*/,
                                           DLIList<BodySM*>& /*new_body_list*/,
                                           bool /*keep_old_body*/ ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : imprint
// Member Type: PUBLIC
// Description: imprint boolean operation on facet-based bodies
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus     OCCModifyEngine::imprint( DLIList<Surface*> &/*surface_list*/,
                                           DLIList<DLIList<Curve*>*> &/*curve_lists_list*/,
                                           BodySM*& /*new_body*/,
                                           bool /*keep_old_body*/ ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : imprint
// Member Type: PUBLIC
// Description: imprint boolean operation on facet-based bodies
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus     OCCModifyEngine::imprint( DLIList<BodySM*> &/*body_list*/,
                                           DLIList<CubitVector*> &/*vector_list*/,
                                           DLIList<BodySM*>& /*new_body_list*/,
                                           bool keep_old /*keep_old_body*/,
                                           DLIList<TopologyBridge*>*,
                                           DLIList<TopologyBridge*>* ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : imprint_projected_edges
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus     OCCModifyEngine::imprint_projected_edges( DLIList<Surface*> &/*ref_face_list*/,
                                                           DLIList<Curve*> &/*ref_edge_list*/,
                                                           DLIList<BodySM*>& /*new_body_list*/,
                                                           bool /*keep_old_body*/,
                                                           bool /*keep_free_edges*/) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : imprint_projected_edges
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus     OCCModifyEngine::imprint_projected_edges(DLIList<Surface*> &/*ref_face_list*/,
                                                           DLIList<BodySM*> &/*body_list*/,
                                                           DLIList<Curve*> &/*ref_edge_list*/,
                                                           DLIList<BodySM*>& /*new_body_list*/,
                                                           bool /*keep_old_body*/,
                                                           bool /*keep_free_edges*/) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : project_edges
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus     OCCModifyEngine::project_edges( DLIList<Surface*> &/*ref_face_list*/,
                                                 DLIList<Curve*> &/*ref_edge_list_in*/,
                                                 DLIList<Curve*> &/*ref_edge_list_new*/,
                                                 bool /*print_error*/ ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}


//===============================================================================
// Function   : intersect
// Member Type: PUBLIC
// Description: intersect boolean operation between facet-based bodies
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus     OCCModifyEngine::intersect(BodySM*  tool_body_ptr,
                                             DLIList<BodySM*>  &from_bodies,
                                             DLIList<BodySM*>  &new_bodies,
                                             bool  keep_old) const
{

  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : chop
// Member Type: PUBLIC
// Description: chop boolean operation between facet-based bodies
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus      OCCModifyEngine::chop(DLIList<BodySM*>& bodies, 
                                         DLIList<BodySM*> &intersectBodies, 
                                         DLIList<BodySM*> &outsideBodies,
                                         BodySM*& leftoversBody,
                                         bool keep_old ,
                                         bool nonreg) const
{
 PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE; 
}

//===============================================================================
// Function   : unite
// Member Type: PUBLIC
// Description: unite boolean operation between facet-based bodies
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus     OCCModifyEngine::unite(DLIList<BodySM*> &bodies, 
                                         DLIList<BodySM*> &newBodies,
                                         bool keep_old) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE; 
}


//===============================================================================
// Function   : thicken
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine::thicken(DLIList<BodySM*>& /*bodies*/, 
                                       DLIList<BodySM*>& /*new_bodies*/,
                                       double /*depth*/,
                                       bool /*both*/) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}


//===============================================================================
// Function   : flip_normals
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine :: flip_normals( DLIList<Surface*>& face_list ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}


//===============================================================================
// Function   : sweep_translational
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine:: sweep_translational(
  DLIList<GeometryEntity*>& /*ref_ent_list*/,
  DLIList<BodySM*>& /*result_body_list*/,
  const CubitVector& /*sweep_vector*/,
  double /*draft_angle*/,
  int /*draft_type*/,
  bool /*switchside*/,
  bool /*rigid*/,
  Surface* stop_surf,
  BodySM* to_body) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : sweep_perpendicular
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine:: sweep_perpendicular(
  DLIList<GeometryEntity*>& /*ref_ent_list*/,
  DLIList<BodySM*>& /*result_body_list*/,
  double /*distance*/,
  double /*draft_angle*/,
  int /*draft_type*/,
  bool /*switchside*/,
  bool /*rigid*/,
  Surface* stop_surf,
  BodySM* to_body) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : sweep_rotational
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine:: sweep_rotational(
  DLIList<GeometryEntity*>& /*ref_ent_list*/,
  DLIList<BodySM*>& /*result_body_list*/,
  const CubitVector& /*point*/,
  const CubitVector& /*direction*/,
  double /*angle*/,
  int /*steps*/,
  double /*draft_angle*/,
  int /*draft_type*/,
  bool /*switchside*/,
  bool /*make_solid*/,
  bool /*rigid*/,
  Surface* stop_surf,
  BodySM* to_body) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : sweep_along_curve
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine::sweep_along_curve(
  DLIList<GeometryEntity*>& /*ref_ent_list*/,
  DLIList<BodySM*>& /*result_body_list*/,
  DLIList<Curve*>& /*ref_edge_list*/,
  double /*draft_angle*/,
  int /*draft_type*/,
  bool /*rigid*/,
  Surface* stop_surf,
  BodySM* to_body) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//HEADER- Webcut-related functions

//===============================================================================
// Function   : webcut
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine::webcut(DLIList<BodySM*>& webcut_body_list,
                              const CubitVector &v1,
                              const CubitVector &v2,
                              const CubitVector &v3,
                              DLIList<BodySM*>& results_list,
                              bool imprint ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : webcut
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus    OCCModifyEngine::webcut(DLIList<BodySM*>& /*webcut_body_list*/,
                                 BodySM const* /*tool_body*/,
                                 DLIList<BodySM*>& /*results_list*/,
                                 bool /*imprint*/ ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : webcut_across_translate
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus    OCCModifyEngine::webcut_across_translate( DLIList<BodySM*>& /*body_list*/,
                                                          Surface* /*plane_surf1*/,
                                                          Surface* /*plane_surf2*/,
                                                          DLIList<BodySM*>& /*results_list*/,
                                                          bool /*imprint*/ ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : webcut_with_sheet
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine::webcut_with_sheet(DLIList<BodySM*> & /*webcut_body_list*/,
                                                 BodySM * /*sheet_body*/,
                                                 DLIList<BodySM*> & /*new_bodies*/,
                                                 bool /*imprint*/ )
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : webcut_with_extended_surf
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine::webcut_with_extended_surf(DLIList<BodySM*> & /*webcut_body_list*/,
                                                         Surface * /*extend_from*/,
                                                         DLIList<BodySM*> & /*new_bodies*/,
                                                         int & /*num_cut*/,
                                                         bool /*imprint*/ )
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : webcut_with_cylinder
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine::webcut_with_cylinder(DLIList<BodySM*> &webcut_body_list,
                                            double radius,
                                            const CubitVector &axis,
                                            const CubitVector &center,
                                            DLIList<BodySM*>& results_list,
                                            bool imprint )
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : webcut_with_brick
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine::webcut_with_brick( 
                                      DLIList<BodySM*>& /*webcut_body_list*/, 
                                      const CubitVector &/*center*/,
                                      const CubitVector* /*axes[3]*/, 
                                      const CubitVector &/*extension*/,
                                      DLIList<BodySM*> &/*results_list*/,
                                      bool /*imprint*/ )
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : webcut_with_planar_sheet
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine::webcut_with_planar_sheet( 
                                          DLIList<BodySM*>& /*webcut_body_list*/,
                                          const CubitVector &/*center*/,
                                          const CubitVector* /*axes[2]*/,
                                          double /*width*/, 
                                          double /*height*/,
                                          DLIList<BodySM*> &/*results_list*/,
                                          bool /*imprint*/ )
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : webcut_with_curve_loop
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine::webcut_with_curve_loop(
                                              DLIList<BodySM*> &/*webcut_body_list*/,
                                              DLIList<Curve*> &/*ref_edge_list*/,
                                              DLIList<BodySM*>& /*results_list*/,
                                              bool /*imprint*/)
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : section
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine::section( DLIList<BodySM*> &/*section_body_list*/,
                                        const CubitVector &/*point_1*/,
                                        const CubitVector &/*point_2*/,
                                        const CubitVector &/*point_3*/,
                                        DLIList<BodySM*>& /*new_body_list*/,
                                        bool /*keep_normal_side*/,
                                        bool /*keep_old*/,
                                        bool /*keep_both_sides*/)
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : split_body
// Member Type: PUBLIC
// Description: Splits multiple lumps in one body into separate bodies
// Author     : Corey Ernst 
// Date       : 08/04
//===============================================================================
CubitStatus OCCModifyEngine::split_body( BodySM *body_ptr,
                                           DLIList<BodySM*> &new_bodies )
{
  return CUBIT_SUCCESS;
}


//===============================================================================
// Function   : reverse_body
// Member Type: PUBLIC
// Description: Turn body inside-out
// Author     : Jason Kraftcheck
// Date       : 05/25/04
//===============================================================================
CubitStatus OCCModifyEngine::reverse_body( BodySM* body_ptr )
{
  return CUBIT_SUCCESS;
}
    


//===============================================================================
// Function   : split_periodic
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine::split_periodic( BodySM * /*body_ptr*/,
                                               BodySM *& /*new_body*/ )
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : regularize_body
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus    OCCModifyEngine::regularize_body( BodySM * /*body_ptr*/,
                                                   BodySM *& /*new_body_ptr*/ )
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : regularize_refentity
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus  OCCModifyEngine::regularize_entity( GeometryEntity * /*old_entity_ptr*/,  
                                                      BodySM *& /*new_body_ptr*/ )
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : offset_curves
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine::offset_curves( DLIList<Curve*>& /*ref_edge_list*/, 
                                              DLIList<Curve*>&,
                                              double /*offset_distance*/,
                                              const CubitVector& /*offset_direction*/, 
                                              int /*gap_type*/ )
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : trim_curve
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
Curve* OCCModifyEngine::trim_curve( Curve* /*trim_curve*/, 
                                      const CubitVector& /*trim_vector*/,
                                      const CubitVector& /*keep_vector*/,
                                      bool )
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return 0;
}

//===============================================================================
// Function   : create_body_from_surfs
// Member Type: PUBLIC
// Description: 
// Author     : Steve Owen
// Date       : 9/11/03
//===============================================================================
CubitStatus OCCModifyEngine::create_solid_bodies_from_surfs(DLIList<Surface*> & ref_face_list,
                                          DLIList<BodySM*>& new_bodies,
                                          bool keep_old,
                                          bool heal) const
{
  return CUBIT_SUCCESS;
}

//===============================================================================
// Function   : create_arc_three
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
Curve* OCCModifyEngine::create_arc_three( Point* /*ref_vertex1*/, 
                                            Point* /*ref_vertex2*/,
                                            Point* /*ref_vertex3*/, 
                                            bool /*full*/ )
{ return NULL;

}

//===============================================================================
// Function   : create_arc_three
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
Curve* OCCModifyEngine::create_arc_three( Curve* /*ref_edge1*/, 
                                            Curve* /*ref_edge2*/,
                                            Curve* /*ref_edge3*/, 
                                            bool /*full*/  )
{ return NULL;
}

//===============================================================================
// Function   : create_arc_center_edge
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
Curve* OCCModifyEngine::create_arc_center_edge( Point* /*ref_vertex1*/, 
                                                  Point* /*ref_vertex2*/,
                                                  Point* /*ref_vertex3*/,
                                                  const CubitVector& /*normal*/, 
                                                  double /*radius*/,
                                                  bool /*full*/ ) 
{ 
  return NULL; 
}

CubitStatus 
OCCModifyEngine::create_curve_combine( DLIList<Curve*>& curve_list, 
                                    Curve *&new_curve_ptr )
{
  PRINT_ERROR("Curve combine is not implemented for facet based models\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : get_gqe
// Member Type: PUBLIC
// Description: get the facet geometry query engince instance pointer
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
GeometryQueryEngine *OCCModifyEngine::get_gqe()
{
  return OCCQueryEngine::instance();
}

//===============================================================================
// Function   : is_modify_engine
// Member Type: PUBLIC
// Description: return CUBIT_TRUE if the tb_ptr belongs to this modify engine
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitBoolean OCCModifyEngine::is_modify_engine(const TopologyBridge *tb_ptr) const 
{
  return tb_ptr->get_geometry_query_engine() == OCCQueryEngine::instance();
}

//===============================================================================
// Function   : get_offset_intersections
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine::get_offset_intersections( Curve* /*ref_edge1*/, 
                                                         Curve* /*ref_edge2*/,
                                                         DLIList<CubitVector*>& /*intersection_list*/,
                                                         double /*offset*/,
                                                         CubitBoolean /*ext_first*/ )
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : get_offset_intersections
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine::get_offset_intersections( Curve* /*ref_edge_ptr*/, 
                                                         Surface* /*ref_face_ptr*/,
                                                         DLIList<CubitVector*> & /*intersection_list*/,
                                                         double /*offset*/,
                                                         CubitBoolean /*ext_surf*/ )
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : surface_intersection
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine::surface_intersection( Surface * /*surface1_ptr*/,
                                                     Surface * /*surface2_ptr*/,
                                                     DLIList<Curve*> &/*inter_graph*/,
                                                     const double /*tol*/) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//===============================================================================
// Function   : get_mid_plane
// Member Type: PUBLIC
// Description: 
// Author     : John Fowler
// Date       : 10/02
//===============================================================================
CubitStatus OCCModifyEngine::get_mid_plane( const CubitVector & /*point_1*/,
                                              const CubitVector & /*point_2*/,
                                              const CubitVector & /*point_3*/,
                                              BodySM * /*body_to_trim_to*/,
                                              BodySM *& /*midplane_body*/ ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

CubitStatus OCCModifyEngine::get_spheric_mid_surface( Surface* surface_ptr1,
                                                        Surface* surface_ptr2,
                                                        BodySM* body_to_trim_to,
                                                        BodySM*& midsurface_body ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

CubitStatus OCCModifyEngine::get_conic_mid_surface( Surface* surface_ptr1,
                                                        Surface* surface_ptr2,
                                                        BodySM* body_to_trim_to,
                                                        BodySM*& midsurface_body ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

CubitStatus OCCModifyEngine::get_toric_mid_surface( Surface* surface_ptr1,
                                                        Surface* surface_ptr2,
                                                        BodySM* body_to_trim_to,
                                                        BodySM*& midsurface_body ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//=============================================================================
// Function   : tweak_chamfer
// Member Type: PUBLIC
// Description: Chamfer curves on solid bodies.  The left and right offsets are
//              with respect to the curve direction.  If the given right offset
//              is negative, the left offset is used.  Users can preview to
//              clarify the meaning of left and right.
// Author     : 
// Date       : 
//=============================================================================
CubitStatus OCCModifyEngine::tweak_chamfer( DLIList<Curve*> & /*curve_list*/, 
                                              double /*left_offset*/,
                                              DLIList<BodySM*> & /*new_bodysm_list*/,
                                              double /*right_offset*/,
                                              CubitBoolean /*keep_old_body*/,
                                              CubitBoolean /*preview*/ ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//=============================================================================
// Function   : tweak_chamfer
// Member Type: PUBLIC
// Description: Chamfer vertices on solid or sheet bodies.  On a solid body 
//              there can be up to 3 offsets; on a sheet body up to 2 offsets.
//              The offsets are in the direction of the supplied edges.  If 
//              multiple vertices are supplied, only one offset value is 
//              allowed and the edges are not used.
// Author     : 
// Date       : 
//=============================================================================
CubitStatus
OCCModifyEngine::tweak_chamfer( DLIList<Point*> & /*point_list*/, 
                                  double /*offset1*/,
                                  DLIList<BodySM*> & /*new_bodysm_list*/,
                                  Curve * /*edge1*/,
                                  double /*offset2*/,
                                  Curve * /*edge2*/,
                                  double /*offset3*/,
                                  Curve * /*edge3*/,
                                  CubitBoolean /*keep_old_body*/,
                                  CubitBoolean /*preview*/ ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//=============================================================================
// Function   : tweak_fillet
// Member Type: PUBLIC
// Description: Create a round fillet (or blend) at the given curves on solid 
//              bodies.
// Author     : 
// Date       : 
//=============================================================================
CubitStatus OCCModifyEngine::tweak_fillet( DLIList<Curve*> & /*curve_list*/, 
                                             double /*radius*/,
                                             DLIList<BodySM*> & /*new_bodysm_list*/,
                                             CubitBoolean /*keep_old_body*/,
                                             CubitBoolean /*preview*/ ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//=============================================================================
// Function   : tweak_fillet
// Member Type: PUBLIC
// Description: Create a round fillet (or blend) at the given curves on a solid 
//              body.  The fillet has a variable radius from the start to the
//              end of the curve.
// Author     : 
// Date       : 
//=============================================================================
CubitStatus OCCModifyEngine::tweak_fillet( Curve * /*curve_ptr*/, 
                                             double /*start_radius*/,
                                             double /*end_radius*/,
                                             BodySM *& /*new_bodysm_ptr*/,
                                             CubitBoolean /*keep_old_body*/,
                                             CubitBoolean /*preview*/ ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//=============================================================================
// Function   : tweak_fillet
// Member Type: PUBLIC
// Description: Create a round fillet (or blend) at the given vertices on sheet
//              bodies.
// Author     : 
// Date       : 
//=============================================================================
CubitStatus
OCCModifyEngine::tweak_fillet( DLIList<Point*> & /*ref_vertex_list*/, 
                                 double /*radius*/,
                                 DLIList<BodySM*> & /*new_bodysm_list*/,
                                 CubitBoolean /*keep_old_body*/,
                                 CubitBoolean /*preview*/ ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//=============================================================================
// Function   : tweak_move
// Member Type: PUBLIC
// Description: Tweak specified faces of a volume or volumes along a vector.
// Author     : 
// Date       : 
//=============================================================================
CubitStatus OCCModifyEngine::tweak_move( DLIList<Surface*> & /*surface_list*/, 
                                           const CubitVector & /*delta*/,
                                           DLIList<BodySM*> & /*new_bodysm_list*/, 
                                           CubitBoolean /*keep_old_body*/ ,
                                           CubitBoolean show_preview) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//=============================================================================
// Function   : tweak_move
// Member Type: PUBLIC
// Description: Tweak specified curves of a sheet body along a vector.
// Author     : 
// Date       : 
//=============================================================================
CubitStatus OCCModifyEngine::tweak_move( DLIList<Curve*> & /*curve_list*/,
                                           const CubitVector & /*delta*/,
                                           DLIList<BodySM*> & /*new_bodysm_list*/, 
                                           CubitBoolean /*keep_old_body*/,
                                           CubitBoolean show_preview ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//=============================================================================
// Function   : tweak_offset
// Member Type: PUBLIC
// Description: Tweak specified faces of a volume or volumes by offsetting
//              those faces by the offset distance.
// Author     : 
// Date       : 
//=============================================================================
CubitStatus OCCModifyEngine::tweak_offset( DLIList<Surface*> & /*surface_list*/, 
                                             double /*offset_distance*/,
                                             DLIList<BodySM*> & /*new_bodysm_list*/,
                                             CubitBoolean /*keep_old_body*/,
                                             CubitBoolean show_preview ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//=============================================================================
// Function   : tweak_offset
// Member Type: PUBLIC
// Description: Tweak specified curves of a sheet body or bodies by offsetting
//              those curves by the offset distance.
// Author     : 
// Date       : 
//=============================================================================
CubitStatus OCCModifyEngine::tweak_offset( DLIList<Curve*> & /*curve_list*/,  
                                             double /*offset_distance*/,
                                             DLIList<BodySM*> & /*new_bodysm_list*/,
                                             CubitBoolean /*keep_old_body*/,
                                             CubitBoolean show_preview ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//=============================================================================
// Function   : tweak_remove
// Member Type: PUBLIC
// Description: Function to remove surfaces from a body and then extend the 
//              remaining surfaces to fill the gap or hole.
// Author     : 
// Date       : 
//=============================================================================
CubitStatus OCCModifyEngine::tweak_remove( DLIList<Surface*> & /*surface_list*/,
                                             DLIList<BodySM*> & /*new_bodysm_list*/,
                                             CubitBoolean /*extend_adjoining*/,
                                             CubitBoolean /*keep_surface*/,
                                             CubitBoolean /*keep_old_body*/,
                                             CubitBoolean show_preview ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//=============================================================================
// Function   : tweak_remove
// Member Type: PUBLIC
// Description: Function to remove curves from a sheet body and then extend the 
//              remaining curves or fill the gap or hole.
// Author     : 
// Date       : 
//=============================================================================
CubitStatus OCCModifyEngine::tweak_remove( DLIList<Curve*> & /*curve_list*/,
                                             DLIList<BodySM*> & /*new_bodysm_list*/, 
                                             CubitBoolean /*keep_old_body*/,
                                             CubitBoolean show_preview ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//=============================================================================
// Function   : tweak_target
// Member Type: PUBLIC
// Description: Tweak specified faces of a volume or volumes up to a target 
//              surface.
// Author     : 
// Date       : 
//=============================================================================
CubitStatus OCCModifyEngine::tweak_target( DLIList<Surface*> & /*surface_list*/,
                                           DLIList<Surface*> & ,
                                           DLIList<BodySM*> & /*new_bodysm_list*/,
                                             CubitBoolean /*reverse_flg*/,
                                             CubitBoolean /*keep_old_body*/,
                                             CubitBoolean show_preview ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//=============================================================================
// Function   : tweak_target
// Member Type: PUBLIC
// Description: Tweak specified edges of a surface or set of surfaces (in sheet
//              bodies) up to a target surface.
// Author     : 
// Date       : 
//=============================================================================
CubitStatus OCCModifyEngine::tweak_target( DLIList<Curve*> & /*curve_list*/,
                                           DLIList<Surface*> & /*target_surfs*/,
                                           DLIList<BodySM*> & /*new_bodysm_list*/, 
                                           CubitBoolean ,
                                           CubitBoolean /*keep_old_body*/,
                                           CubitBoolean show_preview ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//=============================================================================
// Function   : tweak_target
// Member Type: PUBLIC
// Description: Tweak specified edges of a sheet body or bodies up to a target
//              curve that is part of a sheet body.  The target is a surface 
//              created by thickening the owning surface of the target curve.
// Author     : 
// Date       : 
//=============================================================================
CubitStatus OCCModifyEngine::tweak_target( DLIList<Curve*> & /*curve_list*/,
                                           DLIList<Curve*> & /*target_curve_ptr*/, 
                                           DLIList<BodySM*> & /*new_bodysm_list*/, 
                                           CubitBoolean,
                                           CubitBoolean /*keep_old_body*/,
                                           CubitBoolean show_preview ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

CubitStatus OCCModifyEngine::remove_curve_slivers( BodySM* /*body*/,
                                                   double /*lengthlimit*/ ) const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

//================================================================================
// Description: Creates a net surface.
// Author     : Tyronne Lim
// Date       : 08/18/03
//================================================================================
CubitStatus OCCModifyEngine::create_net_surface( DLIList<Surface*>& /*ref_face_list*/, 
                                                   BodySM *& /*new_body*/,
                                                   DLIList<DLIList<CubitVector*>*> & /*vec_lists_u*/, 
                                                   DLIList<DLIList<CubitVector*>*> & /*vec_lists_v*/, 
                                                   double /*net_tol*/, 
                                                   CubitBoolean /*heal*/ ) const
{
   PRINT_ERROR("Function not implemented in this engine.\n");
   return CUBIT_FAILURE;
}

//================================================================================
// Description: Creates a net surface.
// Author     : Tyronne Lim
// Date       : 08/18/03
//================================================================================
CubitStatus OCCModifyEngine::create_net_surface( DLIList<Curve*>& /*u_curves*/, 
                                                   DLIList<Curve*>& /*v_curves*/,
                                                   BodySM *& /*new_body*/, 
                                                   double /*net_tol*/, 
                                                   CubitBoolean /*heal*/ ) const
{
   PRINT_ERROR("Function not implemented in this engine.\n");
   return CUBIT_FAILURE;
}

//================================================================================
// Description: Creates an offset surface.
// Author     : Tyronne Lim
// Date       : 08/18/03
//================================================================================
CubitStatus OCCModifyEngine::create_offset_surface( Surface* /*ref_face_ptr*/, 
                                                      BodySM*& /*new_body*/, 
                                                      double /*offset_distance*/ ) const
{
   PRINT_ERROR("Function not implemented in this engine.\n");
   return CUBIT_FAILURE;
}

//================================================================================
// Description: Creates an offset body.
// Author     : Tyronne Lim
// Date       : 08/18/03
//================================================================================
CubitStatus OCCModifyEngine::create_offset_body( BodySM* body_ptr, 
                                                   BodySM*& new_bodysm, 
                                                   double offset_distance ) const
{
  return CUBIT_FAILURE;
}

//================================================================================
// Description: Creates a skin surface.
// Author     : Tyronne Lim
// Date       : 08/18/03
//================================================================================
CubitStatus OCCModifyEngine::create_skin_surface( DLIList<Curve*>& /*curves*/, 
                                                    BodySM*& /*new_body*/ ) const
{
   PRINT_ERROR("Function not implemented in this engine.\n");
   return CUBIT_FAILURE;
}

//================================================================================
// Description: Creates a body from lofting surfaces.
// Author     : Tyronne Lim
// Date       : 08/18/03
//================================================================================
CubitStatus OCCModifyEngine::loft_surfaces( Surface * /*face1*/, 
                                              const double & /*takeoff1*/,
                                              Surface * /*face2*/, 
                                              const double & /*takeoff2*/,
                                              BodySM*& /*new_body*/,
                                              CubitBoolean /*arc_length_option*/, 
                                              CubitBoolean /*twist_option*/,
                                              CubitBoolean /*align_direction*/, 
                                              CubitBoolean /*perpendicular*/,
                                              CubitBoolean /*simplify_option*/ ) const
{
   PRINT_ERROR("Function not implemented in this engine.\n");
   return CUBIT_FAILURE;
}

//================================================================================
// Description: Creates a body by lofting surfaces between bodies.
// Author     : Tyronne Lim
// Date       : 08/18/03
//================================================================================
CubitStatus OCCModifyEngine::loft_surfaces_to_body( Surface * /*face1*/, 
                                                      const double & /*takeoff1*/,
                                                      Surface * /*face2*/, 
                                                      const double & /*takeoff2*/,
                                                      BodySM*& /*new_body*/,
                                                      CubitBoolean /*arc_length_option*/, 
                                                      CubitBoolean /*twist_option*/,
                                                      CubitBoolean /*align_direction*/, 
                                                      CubitBoolean /*perpendicular*/,
                                                      CubitBoolean /*simplify_option*/ ) const
{
   PRINT_ERROR("Function not implemented in this engine.\n");
   return CUBIT_FAILURE;
}
 
//================================================================================
// Description: Creates a surface.
// Author     : Tyronne Lim
// Date       : 08/18/03
//================================================================================
CubitStatus OCCModifyEngine::create_surface( DLIList<CubitVector*>& /*vec_list*/, 
                                               BodySM *& /*new_body*/, 
                                               Surface * /*ref_face_ptr*/,
                                               CubitBoolean /*project_points*/ ) const
{
   PRINT_ERROR("Function not implemented in this engine.\n");
   return CUBIT_FAILURE;
}

//================================================================================
// Description: Creates a weld surface.
// Author     : Tyronne Lim
// Date       : 08/18/03
//================================================================================
CubitStatus OCCModifyEngine::create_weld_surface( CubitVector & /*root*/,
                                                    Surface * /*ref_face1*/, 
                                                    double /*leg1*/, 
                                                    Surface * /*ref_face2*/, 
                                                    double /*leg2*/,
                                                    BodySM *& /*new_body*/ ) const
{
   PRINT_ERROR("Function not implemented in this engine.\n");
   return CUBIT_FAILURE;
}

CubitStatus OCCModifyEngine::webcut_with_sweep_surfaces(
                                 DLIList<BodySM*> &blank_bodies,
                                 DLIList<Surface*> &surfaces,
                                 const CubitVector& sweep_vector,
                                 bool sweep_perp, 
                                 bool through_all,
                                 bool outward,
                                 bool up_to_next, 
                                 Surface *stop_surf, 
                                 Curve *curve_to_sweep_along, 
                                 DLIList<BodySM*> &results_list,
                                 CubitBoolean imprint)
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

CubitStatus OCCModifyEngine::webcut_with_sweep_curves(
                                 DLIList<BodySM*> &blank_bodies,
                                 DLIList<Curve*> &curves,
                                 const CubitVector& sweep_vector,
                                 bool through_all, 
                                 Surface *stop_surf, 
                                 Curve *curve_to_sweep_along, 
                                 DLIList<BodySM*> &results_list,
                                 CubitBoolean imprint)
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

CubitStatus OCCModifyEngine::webcut_with_sweep_surfaces_rotated(
                                 DLIList<BodySM*> &blank_bodies,
                                 DLIList<Surface*> &surfaces,
                                 const CubitVector &point, 
                                 const CubitVector &sweep_axis, 
                                 double angle, 
                                 Surface *stop_surf, 
                                 bool up_to_next, 
                                 DLIList<BodySM*> &results_list,
                                 CubitBoolean imprint)
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

CubitStatus OCCModifyEngine::webcut_with_sweep_curves_rotated(
                                 DLIList<BodySM*> &blank_bodies,
                                 DLIList<Curve*> &curves,
                                 const CubitVector &point, 
                                 const CubitVector &sweep_axis, 
                                 double angle, 
                                 Surface *stop_surf, 
                                 DLIList<BodySM*> &results_list,
                                 CubitBoolean imprint)
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

CubitStatus OCCModifyEngine::scale( BodySM *&body, const CubitVector& factors )
{
  return OCCQueryEngine::instance()->scale( body, factors );
}

CubitStatus OCCModifyEngine::tolerant_imprint( DLIList<BodySM*> &bodies_in,
                                               DLIList<BodySM*> &new_bodies,
                                               DLIList<TopologyBridge*>*,
                                               DLIList<TopologyBridge*>* )  const
{
  PRINT_ERROR("Option not supported for mesh based geometry.\n");
  return CUBIT_FAILURE;
}

// EOF
