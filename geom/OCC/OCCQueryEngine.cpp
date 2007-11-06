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
#include "gp_Ax1.hxx"
#include "gp_Ax2.hxx"
#include "Geom_Surface.hxx"
#include "Geom_Curve.hxx"
#include "BRepBuilderAPI_Transform.hxx"
#include "TColgp_Array1OfPnt.hxx"
#include "Poly_Array1OfTriangle.hxx"
#include "Poly_Triangle.hxx"
#include "Handle_Poly_Triangulation.hxx"
#include "Poly_Polygon3D.hxx"
#include "Handle_Poly_Polygon3D.hxx"
#include "OCCQueryEngine.hpp"
#include "OCCModifyEngine.hpp"
#include "Poly_Triangulation.hxx"
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
#include "GfxDebug.hpp"
#include "KDDTree.hpp"
#include "RTree.hpp"
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
//#include "TopOpeBRep_ShapeIntersector.hxx"
//#include "BRepAdaptor_Curve.hxx"
//#include "TopOpeBRepTool_ShapeTool.hxx"
//#include "BRepPrimAPI_MakePrism.hxx"
//#include "TopOpeBRep_Point2d.hxx"
#include "BRepExtrema_DistShapeShape.hxx"
#include "BRepAlgoAPI_Section.hxx"
#include "BRepBuilderAPI_MakeEdge.hxx"
#include "gp_Lin.hxx"
using namespace NCubitFile;

OCCQueryEngine* OCCQueryEngine::instance_ = NULL;

const int OCCQueryEngine::OCCQE_MAJOR_VERSION = 6;
const int OCCQueryEngine::OCCQE_MINOR_VERSION = 1;
const int OCCQueryEngine::OCCQE_SUBMINOR_VERSION = 0;

TopTools_DataMapOfShapeInteger *OCCQueryEngine::OCCMap = new TopTools_DataMapOfShapeInteger;
TopTools_DataMapOfShapeInteger *OCCQueryEngine::OCCMapr = new TopTools_DataMapOfShapeInteger;
DLIList<TopologyBridge*> *OCCQueryEngine::CGMList = new DLIList<TopologyBridge*>;

std::map<int, TopologyBridge*>* OccToCGM = new std::map<int, TopologyBridge*>;
typedef std::map<int, TopologyBridge*>::value_type valType;
int iTotalTBCreated = 0;
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
// Description: reflect body about an axis
// Author     : sjowen
// Date       : 9/7/01
//================================================================================
CubitStatus OCCQueryEngine::reflect( BodySM *bodysm,
                                     const CubitVector& axis)
{
  OCCBody *body = CAST_TO(bodysm, OCCBody);
  if (!body)
  {
    PRINT_ERROR("Attempt to reflect OCC-based geometry Body.  This body is not OCCBody.");
    return CUBIT_FAILURE;
  }

  body->reflect( axis.x(), axis.y(), axis.z() );

  return CUBIT_SUCCESS;
}


//================================================================================
// Description:  This function queries OCC for the necessary facets
//               information needed in facetting a RefFace.  This
//               information is stored and output in gMem.  The
//               number of triangles, points and facets are also
//               output.
// Author     :  Jane Hu
// Date       :  10/25/07
//================================================================================
CubitStatus OCCQueryEngine::get_graphics( Surface* surface_ptr,
                                          int& number_triangles,
                                          int& number_points,
                                          int& number_facets,
                                          GMem* g_mem,
                                          unsigned short normal_tolerance,
                                          double distance_tolerance,
                                          double max_edge_length) const
{
    // Because this may be unnecessarily called twice,
    // say there is one triangle.
  if (!g_mem)
  {
    number_triangles = 1;
    number_points = 3;
    number_facets = 4;
    return CUBIT_SUCCESS;
  }

  OCCSurface *occ_surface_ptr = CAST_TO(surface_ptr, OCCSurface);
  TopoDS_Face * Topo_Face = occ_surface_ptr->get_TopoDS_Face();
  if (!Topo_Face)
    return CUBIT_FAILURE;

  TopLoc_Location L;
  Handle_Poly_Triangulation facets = BRep_Tool::Triangulation(*Topo_Face, L);

  //if necessary, the face tolerance can be returned. now, no use.
  //double tol = BRep_Tool::Tolerance(*Topo_Face);   

  number_points = facets->NbNodes();
  number_triangles = facets->NbTriangles();
  assert(number_points == 3 * number_triangles);
  number_facets = 4 * number_triangles; 
  
  Poly_Array1OfTriangle triangles(0, number_facets-1);
  triangles.Assign( facets->Triangles() );
  int *facetList =  new int[number_facets];
  //needs to test that N1, N2, N3 index are starting from 0 to number_points-1
  //otherwise needs to update either facetList or gPnts to make consistent.
  //It's possible also that N's starting from 1.
  int minN = 0;
  for (int i = 0; i < triangles.Length(); i++)
  {
     Poly_Triangle triangle = triangles.Value( i );
     int N1, N2, N3;
     triangle.Get(N1, N2, N3); 
     facetList[4 * i] = 3;
     facetList[4 * i + 1] = N1;
     minN = (minN < N1 ? minN : N1);
     facetList[4 * i + 2] = N2;
     minN = (minN < N2 ? minN : N2);
     facetList[4 * i + 3] = N3;
     minN = (minN < N3 ? minN : N3);
  } 
  g_mem->replace_facet_list( facetList, number_facets, number_facets); 

  TColgp_Array1OfPnt points(minN, minN + number_points-1);
  points.Assign(facets->Nodes());
  GPoint *gPnts= new GPoint[number_points + minN];
  for (int i = minN; i < number_points + minN; i ++)
  {
     gp_Pnt gp_pnt = points.Value(i);
     GPoint gPnt;
     gPnt.x = gp_pnt.X();
     gPnt.y = gp_pnt.Y();
     gPnt.z = gp_pnt.Z();
     gPnts[i] = gPnt;
  }
  g_mem->replace_point_list( gPnts, number_points, number_points + minN);

  return CUBIT_SUCCESS;
}

//================================================================================
// Description: This function queries OCC for the edge information
//              needed in facetting a RefEdge.  This information is
//              stored and output in g_mem.
// Author     : Jane Hu
// Date       : 10/26/07
//================================================================================
CubitStatus OCCQueryEngine::get_graphics( Curve* curve_ptr,
                                          int& num_points,
                                          GMem* gMem,
                                          double /*tolerance*/ ) const
{
  //  get the OCCCurve.
  OCCCurve *occ_curve_ptr = CAST_TO(curve_ptr,OCCCurve);
  assert (gMem);

    
  TopoDS_Edge * Topo_Edge = occ_curve_ptr->get_TopoDS_Edge();
  if (!Topo_Edge)
    return CUBIT_FAILURE;

  TopLoc_Location L;
  Handle_Poly_Polygon3D facets = BRep_Tool::Polygon3D(*Topo_Edge, L);

  num_points = facets->NbNodes();
  TColgp_Array1OfPnt points(0, num_points-1);
  points.Assign(facets->Nodes());

  //! Note: If the polygon is closed, the point of closure is 
  //! repeated at the end of its table of nodes. Thus, on a closed 
  //! triangle the function NbNodes returns 4. 
  GPoint *gPnts= new GPoint[num_points];
  for (int i = 0; i < num_points ; i ++)
  {
     gp_Pnt gp_pnt = points.Value(i);
     GPoint gPnt;
     gPnt.x = gp_pnt.X();
     gPnt.y = gp_pnt.Y();
     gPnt.z = gp_pnt.Z();
     gPnts[i] = gPnt;
  }
  gMem->replace_point_list( gPnts, num_points, num_points );
 
  return CUBIT_SUCCESS;
}

//================================================================================
// Description: Given surface and number of point on u and v parametric
//              direction, find the 3-d point locations
// Author     : Jane Hu
// Date       : 10/22/07
//================================================================================
CubitStatus OCCQueryEngine::get_isoparametric_points(Surface* surface,
                                                     int &nu, int &nv,
                                                     GMem*& g_mem) const
{
  OCCSurface* occ_surface = CAST_TO(surface, OCCSurface);
  TopoDS_Face* Tops_face = occ_surface->get_TopoDS_Face();
  TopoDS_Face the_face;
  if (Tops_face == NULL)
  {
    PRINT_ERROR("This surface is not OCCSurface.");
    return CUBIT_FAILURE;
  }

  the_face = *Tops_face;
  Handle_Geom_Surface HGeom_surface = BRep_Tool::Surface(the_face);

  assert (nu > 1 && nv > 1);
  double u1, u2, v1, v2;
  HGeom_surface->Bounds(u1, u2, v1, v2); 
  double interval1  = (u2 - u1)/(nu-1);
  double interval2 = (v2 - v1)/(nv-1);
  
  g_mem = new GMem[nu];
  //nu and nv must be given to calculate the points.
  for (int i = 0; i < nu; i++)
  {
     Handle_Geom_Curve HGeom_curve = HGeom_surface->UIso(u1 + i * interval1); 
     g_mem[i].allocate_polylines(nv-1);
     for (int j = 0; j <  nv; j++)
     {
        gp_Pnt pnt = HGeom_curve->Value(v1 + j * interval2);
        g_mem[i].point_list()[j].x = pnt.X();
        g_mem[i].point_list()[j].y = pnt.Y();
        g_mem[i].point_list()[j].z = pnt.Z();
     }
     g_mem[i].pointListCount = nv;

  }
  return CUBIT_SUCCESS;
}

CubitStatus OCCQueryEngine::get_u_isoparametric_points(Surface* surface,
                                                       double v, int&n,
                                                       GMem*& g_mem) const
{
  OCCSurface* occ_surface = CAST_TO(surface, OCCSurface);
  TopoDS_Face* Tops_face = occ_surface->get_TopoDS_Face();
  TopoDS_Face the_face;
  if (Tops_face == NULL)
  {
    PRINT_ERROR("This surface is not OCCSurface.");
    return CUBIT_FAILURE;
  }

  the_face = *Tops_face;

  Handle_Geom_Surface HGeom_surface = BRep_Tool::Surface(the_face);

  //n must be given to calculate the points.
  assert (n > 1);
  double u1, u2, v1, v2;
  HGeom_surface->Bounds(u1, u2, v1, v2);
  double interval = (u2 - u1)/(n -1); 
  
  Handle_Geom_Curve HGeom_curve = HGeom_surface->VIso(v);
  g_mem = new GMem;
  g_mem->allocate_polylines(n-1);
  for (int j = 0; j < n; j++)
  {
      gp_Pnt pnt = HGeom_curve->Value(u1 + j * interval);
      g_mem->point_list()[j].x = pnt.X();
      g_mem->point_list()[j].y = pnt.Y();
      g_mem->point_list()[j].z = pnt.Z();
   }
   g_mem->pointListCount = n;

   return CUBIT_SUCCESS;
}

CubitStatus OCCQueryEngine::get_v_isoparametric_points(Surface* surface,
                                                       double u, int&n,
                                                       GMem*&g_mem) const
{
  OCCSurface* occ_surface = CAST_TO(surface, OCCSurface);
  TopoDS_Face* Tops_face = occ_surface->get_TopoDS_Face();
  TopoDS_Face the_face;
  if (Tops_face == NULL)
  {
    PRINT_ERROR("This surface is not OCCSurface.");
    return CUBIT_FAILURE;
  }

  the_face = *Tops_face;

  Handle_Geom_Surface HGeom_surface = BRep_Tool::Surface(the_face);

  //n must be given to calculate the points.
  assert (n > 1);
  double u1, u2, v1, v2;
  HGeom_surface->Bounds(u1, u2, v1, v2);
  double interval = (v2 - v1)/(n -1);

  Handle_Geom_Curve HGeom_curve = HGeom_surface->UIso(u);
  g_mem = new GMem;
  g_mem->allocate_polylines(n-1);
  for (int j = 0; j < n; j++)
  {
      gp_Pnt pnt = HGeom_curve->Value(v1 + j * interval);
      g_mem->point_list()[j].x = pnt.X();
      g_mem->point_list()[j].y = pnt.Y();
      g_mem->point_list()[j].z = pnt.Z();
   }
   g_mem->pointListCount = n;

   return CUBIT_SUCCESS;
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
// Description:  Calculate for intersection points between a curve and 
//               a segment defined by two points or
//               between two curves or between a curve and a surface.
// Author     :  Jane Hu
// Date       :  10/15/07
//================================================================================
CubitStatus OCCQueryEngine::get_intersections( Curve* curve, 
                                               CubitVector& point1,
                                               CubitVector& point2,
                                               DLIList<CubitVector*>& intscts,
                                               CubitBoolean bounded,
                                               CubitBoolean closest)
{
  OCCCurve *occ_curve =  CAST_TO(curve, OCCCurve);
  if (occ_curve == NULL)
  {
     PRINT_ERROR("Option not supported for non-occ based geometry.\n");
     return CUBIT_FAILURE;
  }

  OCCPoint *pt1 = new OCCPoint(point1);
  OCCPoint *pt2 = new OCCPoint(point2);
  Curve *curve2 = 
             OCCModifyEngine::instance()->make_Curve(pt1, pt2);
  if (curve2 == NULL)
  {
    PRINT_ERROR( "Unable to create OCC EDGE from points\n" );
    return CUBIT_FAILURE;
  }
  
  OCCCurve *occ_curve2 = CAST_TO(curve2, OCCCurve);
  return get_intersections(occ_curve, occ_curve2, intscts, bounded, closest);
}

CubitStatus OCCQueryEngine::get_intersections( Curve* curve1, 
                                               Curve* curve2,
                                               DLIList<CubitVector*>& intscts,
                                               CubitBoolean bounded,
                                               CubitBoolean closest)
{
  OCCCurve *occ_curve1 =  CAST_TO(curve1, OCCCurve);
  if (occ_curve1 == NULL)
  {
     PRINT_ERROR("Option not supported for non-occ based geometry.\n");
     return CUBIT_FAILURE;
  }

  OCCCurve *occ_curve2 =  CAST_TO(curve2, OCCCurve);
  if (occ_curve2 == NULL)
  {
     PRINT_ERROR("Option not supported for non-occ based geometry.\n");
     return CUBIT_FAILURE;
  }

  //currently, there's no effect on 'closest' argument or bounded.
  BRepExtrema_DistShapeShape distShapeShape(
                                            *(occ_curve1->get_TopoDS_Edge()),
                                            *(occ_curve2->get_TopoDS_Edge()));

  //distShapeShape.Perform();
  if (!distShapeShape.IsDone())
  {
     PRINT_ERROR("Cannot calculate the intersection points for the input curves.\n");
     return CUBIT_FAILURE;
  }
  
  if (distShapeShape.Value() < get_sme_resabs_tolerance())
  {
     int numPnt = distShapeShape.NbSolution();
     for (int i = 1; i <= numPnt; i++)
     {
        gp_Pnt aPoint = distShapeShape.PointOnShape1(i);
     
        CubitVector* cv = new CubitVector(aPoint.X(), aPoint.Y(), aPoint.Z());
        intscts.append(cv);
     }
  }
  return CUBIT_SUCCESS;
}

CubitStatus
OCCQueryEngine::get_intersections( Curve* curve, Surface* surface,
                                   DLIList<CubitVector*>& intscts,
                                   CubitBoolean bounded )
{
  // There's no effect of bounded =  false. 
  OCCCurve *occ_curve =  CAST_TO(curve, OCCCurve);
  if (occ_curve == NULL)
  {
     PRINT_ERROR("Option not supported for non-occ based geometry.\n");
     return CUBIT_FAILURE;
  }

  OCCSurface *occ_surface =  CAST_TO(surface, OCCSurface);
  if (occ_surface == NULL)
  {
     PRINT_ERROR("Option not supported for non-occ based geometry.\n");
     return CUBIT_FAILURE;
  }
   
  //currently, there's no effect on 'closest' argument or bounded.
  BRepExtrema_DistShapeShape distShapeShape(*(occ_curve->get_TopoDS_Edge()),
                                            *(occ_surface->get_TopoDS_Face()));

  //distShapeShape.Perform();
  if (!distShapeShape.IsDone())
  {
     PRINT_ERROR("Cannot calculate the intersection points for the input curve and surface.\n");
     return CUBIT_FAILURE;
  }
  
  if (distShapeShape.Value() < get_sme_resabs_tolerance())
  {
     int numPnt = distShapeShape.NbSolution();
     for (int i = 1; i <= numPnt; i++)
     {
        gp_Pnt aPoint = distShapeShape.PointOnShape1(i);
     
        CubitVector* cv = new CubitVector(aPoint.X(), aPoint.Y(), aPoint.Z());
        intscts.append(cv);
     }
  }
 
  return CUBIT_SUCCESS;
}

//================================================================================
// Description: Find extrema position on an entity list
// Author     : Jane Hu
// Date       : 10/30/07
//================================================================================
CubitStatus
OCCQueryEngine::entity_extrema( DLIList<GeometryEntity*> &ref_entity_list,
                                  const CubitVector *dir1,
                                  const CubitVector *dir2,
                                  const CubitVector *dir3,
                                  CubitVector &extrema,
                                  GeometryEntity *&extrema_entity_ptr )
{
  //in Acis, the api_entity_extrema is used to calculate "possible 
  //self-intersecting sweeping and to align lofting sections"
  PRINT_ERROR("There's no such call in OCC correponding to Acis call."); 
  return CUBIT_FAILURE;
}

//================================================================================
// Description: Find distance between two entities and closest positions.
// Author     : Jane Hu
// Date       : 10/19/07
//================================================================================
CubitStatus
OCCQueryEngine::entity_entity_distance( GeometryEntity *entity1,
                                          GeometryEntity *entity2,
                                          CubitVector &pos1, CubitVector &pos2,
                                          double &distance )
{
  TopoDS_Shape * shape1;
  TopoDS_Shape * shape2;
  if ((shape1 = get_TopoDS_Shape_of_entity(entity1)) == NULL)
  {
    PRINT_ERROR( "problem occured getting OCC entity.\n"
      "       Aborting.\n" );
    return CUBIT_FAILURE;
  }

  if( (shape2 = get_TopoDS_Shape_of_entity( entity2 )) == NULL )
  {
    PRINT_ERROR( "problem occured getting OCC entity.\n"
      "       Aborting.\n");
    return CUBIT_FAILURE;
  }

  BRepExtrema_DistShapeShape distShapeShape(*shape1, *shape2);
  //distShapeShape.Perform();
  
  if (!distShapeShape.IsDone())
  {
    PRINT_ERROR( "problem occured getting distance between OCC entities.\n"
      "       Aborting.\n");
    return CUBIT_FAILURE;
  }

  distance = distShapeShape.Value();
  gp_Pnt pnt1 = distShapeShape.PointOnShape1(1);
  gp_Pnt pnt2 = distShapeShape.PointOnShape2(2);
  pos1 = CubitVector(pnt1.X(), pnt1.Y(), pnt1.Z());
  pos2 = CubitVector(pnt2.X(), pnt2.Y(), pnt2.Z());
  return CUBIT_SUCCESS;
}

TopoDS_Shape* OCCQueryEngine::get_TopoDS_Shape_of_entity(TopologyBridge *entity_ptr)
{
  if (OCCBody *body_ptr = CAST_TO( entity_ptr, OCCBody))
  {
    TopoDS_Shape* theShape = body_ptr->get_TopoDS_Shape();
    if (!theShape)
    {
      PRINT_ERROR("OCCBody without TopoDS_Shape at %s:%d.\n", __FILE__, __LINE__ );
      return NULL;
    }
    return theShape;
  }

  else if (OCCLump * lump_ptr = CAST_TO( entity_ptr,OCCLump))
  {
    TopoDS_Solid * theSolid = lump_ptr->get_TopoDS_Solid();
    if(theSolid)
      return (TopoDS_Shape*) theSolid; 
    else
    {
      PRINT_ERROR("OCCLump without TopoDS_Solid at %s:%d.\n", __FILE__, __LINE__ );
      return NULL;
    }
  }

  else if( OCCSurface *surface_ptr = CAST_TO( entity_ptr, OCCSurface))
  {
    TopoDS_Face *theFace = surface_ptr->get_TopoDS_Face();
    if(!theFace)
    {
      PRINT_ERROR("OCCSurface without TopoDS_Face at %s:%d.\n", __FILE__, __LINE__ );
      return NULL;
    }

    return (TopoDS_Shape*) theFace;
  }

  else if( OCCCurve *curve_ptr = CAST_TO( entity_ptr, OCCCurve))
  {
    TopoDS_Edge *theEdge = curve_ptr->get_TopoDS_Edge();
    if (!theEdge)
    {
      PRINT_ERROR("OCCCurve without TopoDS_Edge at %s:%d.\n", __FILE__, __LINE__ );
      return NULL;
    }

    return (TopoDS_Shape*) theEdge;
  }

  else if( OCCPoint *point_ptr = CAST_TO( entity_ptr, OCCPoint))
  {
    TopoDS_Vertex *thePoint = point_ptr->get_TopoDS_Vertex(); 
    if (!thePoint)
    {
      PRINT_ERROR("OCCPoint without TopoDS_Point at %s:%d.\n", __FILE__, __LINE__ );
      return NULL;
    }

    return (TopoDS_Shape*) thePoint;
  }
  
  PRINT_ERROR("Non-OCC TopologyBridge at %s:%d.\n", __FILE__, __LINE__ );
  return NULL;

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
OCCQueryEngine::import_temp_geom_file(FILE* file_ptr,
                                        const char* /*file_name*/,
                                        const char* file_type,
                                        DLIList<TopologyBridge*> &bridge_list )
{
  //make sure that file_type == "OCC"
  if( !strcmp( file_type,"OCC") )
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
		//surface->add_shell((OCCShell*)(*CGMList)[OCCMap->Find(aShape)]);
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
			//populate_topology_bridge_vertex(*poedge, imported_entities, curve);
		} else {
			coedge = (OCCCoEdge*)(*CGMList)[Map->Find(*poedge)];
		}
		coedge->add_loop((OCCLoop*)occ_to_cgm(aShape));
	}
	return CUBIT_SUCCESS;
}

Point* OCCQueryEngine::populate_topology_bridge_vertex(TopoDS_Shape aShape)
{
	OCCPoint *point;
	TopExp_Explorer Ex;
	for (Ex.Init(aShape, TopAbs_VERTEX); Ex.More(); Ex.Next()) {
		TopoDS_Vertex *povertex = new TopoDS_Vertex;
		*povertex = TopoDS::Vertex(Ex.Current());
		if (!OCCMap->IsBound(*povertex)) 
                {
	 	  printf("Adding vertex\n");
                  iTotalTBCreated++;
	 	  point = new OCCPoint(povertex);
		  CGMList->append(point);
//		  imported_entities.append(point);
		  OCCMap->Bind(*povertex, iTotalTBCreated);
                  OccToCGM->insert(valType(iTotalTBCreated,
                                          (TopologyBridge*)point));
		} 
                else {
                  int i = OCCMap->Find(*povertex);
	 	  point = (OCCPoint*)(OccToCGM->find(i))->second;
		}
	}
	return point;
}

TopologyBridge* OCCQueryEngine::occ_to_cgm(TopoDS_Shape shape)
{
       if(!OCCMap->IsBound(shape))
	 return (TopologyBridge*) NULL;

       int k = OCCMap->Find(shape);
       return (OccToCGM->find(k))->second;
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

  //int num_points, num_edges, num_facets;
  //int num_cfet, num_fet;

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
  /*
  status = restore_facets( file_ptr, file_reader.get_endian(),
                           num_points, num_edges,
                           num_facets, points_array, num_cfet,
                           num_fet, cfet_array, fet_array );
  */
  if( status == CUBIT_FAILURE)
  {
    PRINT_ERROR("Problems restore facets\n");
    return CUBIT_FAILURE;
  }

  //Restore Topology
  /*
  status = restore_topology( file_ptr, file_reader.get_endian(),
                             num_points, points_array,
                             num_cfet, cfet_array, num_fet,
                             fet_array, imported_entities);
  */
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
// Creator       : Jane Hu
//
// Creation Date : 10/29/07
//-------------------------------------------------------------------------
CubitStatus
OCCQueryEngine::delete_solid_model_entities( BodySM* bodysm ) const
{
  OCCBody* fbody = dynamic_cast<OCCBody*>(bodysm);
  if (!fbody)
    return CUBIT_FAILURE;

  TopoDS_Shape* shape = fbody->get_TopoDS_Shape();

  if (!shape)
    return CUBIT_FAILURE;

  //remove the entry from the map
  int k;
  if(OCCMap->IsBound(*shape))
  {
    k = OCCMap->Find(*shape);
    
    if(!OCCMap->UnBind(*shape))
      PRINT_ERROR("The OccBody and TopoDS_Shape pair is not in the map!");

    if(!OccToCGM->erase(k))
      PRINT_ERROR("The OccBody and TopoDS_Shape pair is not in the map!");
  }
  // Remove the links between OCC and Cubit
  //  unhook_ENTITY_from_VGI(shape);

  delete shape;
  delete bodysm;

  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Delete a OCCSurface and child entities.
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 10/29/07
//-------------------------------------------------------------------------
CubitStatus
OCCQueryEngine::delete_solid_model_entities( Surface* surface ) const
{
  OCCSurface* fsurf = dynamic_cast<OCCSurface*>(surface);
  if (!fsurf)
    return CUBIT_FAILURE;

  TopoDS_Face* face = fsurf->get_TopoDS_Face();

  if(!face)
     return CUBIT_FAILURE;

  //remove the entry from the map
  int k;
  if(OCCMap->IsBound(*face))
  {
    k = OCCMap->Find(*face);

    if(!OCCMap->UnBind(*face))
      PRINT_ERROR("The OccSurface and TopoDS_Face pair is not in the map!");

    if(!OccToCGM->erase(k))
      PRINT_ERROR("The OccSurface and TopoDS_Face pair is not in the map!");
  }

  // Remove the links between OCC and Cubit
  //  unhook_ENTITY_from_VGI(face);

  delete face;
  delete surface;
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
  if (!fcurve )
    return CUBIT_FAILURE;

  TopoDS_Edge * edge = fcurve->get_TopoDS_Edge();
  if (!edge)
    return CUBIT_FAILURE;
 
  //remove the entry from the map
  int k;
  if(OCCMap->IsBound(*edge))
  {
    k = OCCMap->Find(*edge);

    if(!OCCMap->UnBind(*edge))
      PRINT_ERROR("The OccCurve and TopoDS_Edge pair is not in the map!");

    if(!OccToCGM->erase(k))
      PRINT_ERROR("The OccCurve and TopoDS_Edge pair is not in the map!");
  }

  // Remove the links between OCC and Cubit
  //unhook_ENTITY_from_VGI(edge);

  delete edge;
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
  OCCPoint* fpoint = dynamic_cast<OCCPoint*>(point);
  if (!fpoint)
    return CUBIT_FAILURE;

  TopoDS_Vertex *vertex = fpoint->get_TopoDS_Vertex();
  if (!vertex)
    return CUBIT_FAILURE;

  //remove the entry from the map
  int k;
  if(OCCMap->IsBound(*vertex))
  {
    k = OCCMap->Find(*vertex);

    if(!OCCMap->UnBind(*vertex))
      PRINT_ERROR("The OccPoint and TopoDS_Vertex pair is not in the map!");

    if(!OccToCGM->erase(k))
      PRINT_ERROR("The OccPoint and TopoDS_Vertex pair is not in the map!");
  }

  // Remove the links between OCC and Cubit
  //unhook_ENTITY_from_VGI(vertex);

  delete vertex;
  delete point;
  return CUBIT_SUCCESS;
}

CubitStatus OCCQueryEngine::fire_ray(BodySM * body,
                                     const CubitVector &start,
                                     const CubitVector &unit,
                                     DLIList<double>& ray_parms,
                                    DLIList<GeometryEntity*> *entity_list) const
{
  CubitStatus status = CUBIT_SUCCESS;

  //- fire a ray at the specified body, returning the entities hit and
  //- the parameters along the ray; return CUBIT_FAILURE if error
  // - line body intersection. 
  gp_Pnt p(start.x(), start.y(), start.z());
  gp_Dir dir(unit.x(), unit.y(), unit.z());
  gp_Lin L(p, dir);
  TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(L); 
  
  OCCBody *occBody = CAST_TO(body, OCCBody);
  if (occBody == NULL)
  {
     PRINT_ERROR("Option not supported for non-occ based geometry.\n");
     return CUBIT_FAILURE;
  }

  BRepExtrema_DistShapeShape distShapeShape(edge,
					*(occBody->get_TopoDS_Shape()));
  //distShapeShape.Perform();
  if (!distShapeShape.IsDone())
  {
     PRINT_ERROR("Cannot calculate the intersection points for the input body.\n");
     return CUBIT_FAILURE;
  }

   if (distShapeShape.Value() < get_sme_resabs_tolerance())
  {
     int numPnt = distShapeShape.NbSolution();
     for (int i = 1; i <= numPnt; i++)
     {
       double para;
       distShapeShape.ParOnEdgeS1(i , para);
       ray_parms.append(para);

       TopoDS_Shape shape = distShapeShape.SupportOnShape2(i);
       int k = OCCMap->Find(shape);
       std::map<int,TopologyBridge*>::iterator it = OccToCGM->find(k);
       TopologyBridge* tb = (*it).second;
       entity_list->append((GeometryEntity*)tb);
     }
  } 
  return status;
}
double OCCQueryEngine::get_sme_resabs_tolerance() const
{
  PRINT_WARNING("OCC doesn't have its standard linear tolerance.\n");
  return 1e-6; 
}
// Gets solid modeler's resolution absolute tolerance

double OCCQueryEngine::set_sme_resabs_tolerance( double )
{
  PRINT_ERROR("OCCQueryEngine::set_sme_resabs_tolerance not yet implemented.\n");
  return 0.0;
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


//=============================================================================
//Function:   create_super_facet_bounding_box(PUBLIC)
//Description: Find the bounding box of a list of BodySMs
//Author: jdfowle
//Date: 12/15/03
//=============================================================================
CubitStatus OCCQueryEngine::create_super_bounding_box(
                                DLIList<BodySM*>& body_list,
                                CubitBox& super_box )
{
BodySM *bodySM;
int i;
CubitStatus status = CUBIT_SUCCESS;

  body_list.reset();
  for ( i = 0; i < body_list.size(); i++ ) {
    bodySM = body_list.get_and_step();  
    OCCBody* occBody = CAST_TO(bodySM, OCCBody);
    super_box |= occBody->get_bounding_box();
  }

  return status;
}

const char* fqe_xform_err = "Transform not implemented for OCC geometry.\n";
CubitStatus OCCQueryEngine::restore_transform( BodySM* body )
{
  OCCBody* facetbod = dynamic_cast<OCCBody*>(body);
  return facetbod ? facetbod->restore( ) : CUBIT_FAILURE;
}
CubitStatus OCCQueryEngine::translate( BodySM* body, const CubitVector& d )
{
  OCCBody* theBody = dynamic_cast<OCCBody*>(body);
  return theBody ? theBody->move( d.x(), d.y(), d.z() ) : CUBIT_FAILURE;
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

//-------------------------------------------------------------------------
// Purpose       : Transform a Solid, Surface, Curve, or Vertex
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 10/23/07
//-------------------------------------------------------------------------
CubitStatus OCCQueryEngine::translate( GeometryEntity* entity,
                                       const CubitVector& v )
{
  TopoDS_Shape * shape;
  if ((shape = get_TopoDS_Shape_of_entity(entity)) == NULL)
  {
    PRINT_ERROR( "problem occured getting OCC entity.\n"
      "       Aborting.\n" );
    return CUBIT_FAILURE;
  }

  gp_Vec aVec(v.x(), v.y(),v.z());
  gp_Trsf aTrsf;  
  aTrsf.SetTranslation(aVec);

  BRepBuilderAPI_Transform aBRepTrsf(*shape, aTrsf);
  return CUBIT_SUCCESS;
}

//a is angular value of rotation in radians
CubitStatus OCCQueryEngine::rotate( GeometryEntity* entity,
                                    const CubitVector& v, double a )
{
  TopoDS_Shape * shape;
  if ((shape = get_TopoDS_Shape_of_entity(entity)) == NULL)
  {
    PRINT_ERROR( "problem occured getting OCC entity.\n"
      "       Aborting.\n" );
    return CUBIT_FAILURE;
  }

  gp_Pnt aOrigin(0,0,0);
  gp_Dir aDir(v.x(), v.y(), v.z());
  gp_Ax1 anAxis(aOrigin, aDir);

  //a is angular value of rotation in radians 
  gp_Trsf aTrsf;
  aTrsf.SetRotation(anAxis, a);

  BRepBuilderAPI_Transform aBRepTrsf(*shape, aTrsf);
  return CUBIT_SUCCESS;
}

CubitStatus OCCQueryEngine::scale( GeometryEntity* entity, double f )
{
  TopoDS_Shape * shape;
  if ((shape = get_TopoDS_Shape_of_entity(entity)) == NULL)
  {
    PRINT_ERROR( "problem occured getting OCC entity.\n"
      "       Aborting.\n" );
    return CUBIT_FAILURE;
  }

  gp_Trsf aTrsf;
  aTrsf.SetScaleFactor(f);

  BRepBuilderAPI_Transform aBRepTrsf(*shape, aTrsf);
  return CUBIT_SUCCESS;
}

CubitStatus OCCQueryEngine::scale( GeometryEntity* , const CubitVector&  )
{
  PRINT_ERROR("non-uniform scaling is not performed on OCC bodies");
  return CUBIT_FAILURE;
}

// like ACIS, here v is the normal of symmetric plane.
CubitStatus OCCQueryEngine::reflect( GeometryEntity* entity, 
                                     const CubitVector&  v)
{
  TopoDS_Shape * shape;
  if ((shape = get_TopoDS_Shape_of_entity(entity)) == NULL)
  {
    PRINT_ERROR( "problem occured getting OCC entity.\n"
      "       Aborting.\n" );
    return CUBIT_FAILURE;
  }
 
  gp_Pnt aOrigin(0,0,0);
  gp_Dir aDir(v.x(), v.y(), v.z());
  gp_Ax2 anAx2(aOrigin, aDir);

  gp_Trsf aTrsf;
  aTrsf.SetMirror(anAx2);

  BRepBuilderAPI_Transform aBRepTrsf(*shape, aTrsf);
  return CUBIT_SUCCESS;
}

//===============================================================================
// Function   : bodies_overlap
// Member Type: PUBLIC
// Description: determine if OCC-based bodies overlap
// Author     : 
// Date       : 10/07
//===============================================================================
CubitBoolean OCCQueryEngine::bodies_overlap (BodySM * body_ptr_1,
                                             BodySM * body_ptr_2 ) const
{
  OCCBody *occ_body1 = CAST_TO(body_ptr_1, OCCBody);
  if (!occ_body1)
  {
    PRINT_ERROR("Can't calculate intersection of non-OCC bodies.");
    return CUBIT_FALSE;
  }
 
  OCCBody *occ_body2 = CAST_TO(body_ptr_2, OCCBody);
  if (!occ_body2)
  {
    PRINT_ERROR("Can't calculate intersection of non-OCC bodies.");
    return CUBIT_FALSE;
  }

  CubitBox box_1 = occ_body1->get_bounding_box();
  CubitBox box_2 = occ_body2->get_bounding_box();
  if ( !box_1.overlap( GEOMETRY_RESABS, box_2 ) )
    return CUBIT_FALSE;

  TopoDS_Shape *shape1 = occ_body1->get_TopoDS_Shape();
  TopoDS_Shape *shape2 = occ_body2->get_TopoDS_Shape(); 
  
  //BRepAlgoAPI_Section calculates intersection between faces only.
  TopExp_Explorer Ex1, Ex2;
  for (Ex1.Init(*shape1, TopAbs_SOLID); Ex1.More(); Ex1.Next())
  {
    TopoDS_Solid *posolid1 =  new TopoDS_Solid;
    *posolid1 = TopoDS::Solid(Ex1.Current());
    OCCLump * lump1 = new OCCLump(posolid1); 
    for (Ex2.Init(*shape2, TopAbs_SOLID); Ex2.More(); Ex2.Next())
    {
       TopoDS_Solid *posolid2 =  new TopoDS_Solid;
       *posolid2 = TopoDS::Solid(Ex2.Current());
       OCCLump * lump2 = new OCCLump(posolid2);
       CubitBoolean is_overlap = volumes_overlap(lump1, lump2);
       if(is_overlap)
	 return CUBIT_TRUE;
    }
  }
  return CUBIT_FALSE;
}

CubitBoolean OCCQueryEngine::volumes_overlap (Lump *lump1, Lump *lump2 ) const
{
  OCCLump *occ_lump1 = CAST_TO(lump1, OCCLump);
  if (!occ_lump1)
  {
    PRINT_ERROR("Can't calculate intersection of non-OCC solids.");
    return CUBIT_FALSE;
  }

  OCCLump *occ_lump2 = CAST_TO(lump2, OCCLump);
  if (!occ_lump2)
  {
    PRINT_ERROR("Can't calculate intersection of non-OCC solids.");
    return CUBIT_FALSE;
  }

  CubitBox box_1 = occ_lump1->bounding_box();
  CubitBox box_2 = occ_lump2->bounding_box();
  if ( !box_1.overlap( GEOMETRY_RESABS, box_2 ) )
    return CUBIT_FALSE;

  TopoDS_Shape *shape1 = (TopoDS_Shape*)(occ_lump1->get_TopoDS_Solid());
  TopoDS_Shape *shape2 = (TopoDS_Shape*)(occ_lump2->get_TopoDS_Solid());
  
  //BRepAlgoAPI_Section calculates intersection between faces only.
  TopExp_Explorer Ex1, Ex2;
  for (Ex1.Init(*shape1, TopAbs_FACE); Ex1.More(); Ex1.Next())  
  {
    for (Ex2.Init(*shape2, TopAbs_FACE); Ex2.More(); Ex2.Next()) 
    {
        BRepAlgoAPI_Section section(Ex1.Current(), Ex2.Current());
        if (section.HasGenerated())
  	  return CUBIT_TRUE;
    }
  }
  return CUBIT_FALSE;
}

//EOF
