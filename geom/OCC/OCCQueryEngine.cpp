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
#include "BRepTools_WireExplorer.hxx"
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
#include "TopoDS_Compound.hxx"
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
#include "TopTools_DataMapOfShapeInteger.hxx"
#include "BRepExtrema_DistShapeShape.hxx"
#include "BRepAlgoAPI_Section.hxx"
#include "BRepBuilderAPI_MakeEdge.hxx"
#include "gp_Lin.hxx"
using namespace NCubitFile;

OCCQueryEngine* OCCQueryEngine::instance_ = NULL;

const int OCCQueryEngine::OCCQE_MAJOR_VERSION = 6;
const int OCCQueryEngine::OCCQE_MINOR_VERSION = 2;
const int OCCQueryEngine::OCCQE_SUBMINOR_VERSION = 0;

TopTools_DataMapOfShapeInteger *OCCQueryEngine::OCCMap = new TopTools_DataMapOfShapeInteger;
TopTools_DataMapOfShapeInteger *OCCQueryEngine::OCCMapr = new TopTools_DataMapOfShapeInteger;
DLIList<TopologyBridge*> *OCCQueryEngine::CGMList = new DLIList<TopologyBridge*>;

std::map<int, TopologyBridge*>* OCCQueryEngine::OccToCGM = new std::map<int, TopologyBridge*>;
typedef std::map<int, TopologyBridge*>::value_type valType;
int OCCQueryEngine::iTotalTBCreated = 0;
CubitBoolean PRINT_RESULT = CUBIT_FALSE;
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
  temp_filename += ".occ";

  if( export_solid_model( ref_entity_list, temp_filename.c_str(), "OCC",
                          cubit_version ) == CUBIT_FAILURE )
  {
    PRINT_ERROR( "Error occured while trying to save temporary OCC_BASED_GEOMETRY file\n");
    return CUBIT_FAILURE;
  }

  int size_after = ref_entity_list.size();

  if( size_before > size_after )
  {
    created_file +=  temp_filename;
    created_file_type += "OCC";
  }
  return CUBIT_SUCCESS;
}

//===========================================================================
//Function Name:export_solid_model
//Member Type:  PUBLIC
//Description:  function called for save/restore to save temporary FACET file
//Author:       Jane Hu
//Date:         11/16/2007
//===========================================================================

CubitStatus OCCQueryEngine::export_solid_model( DLIList<TopologyBridge*>& ref_entity_list,
                                                     const char* file_name,
                                                     const char* file_type,
                                                     const CubitString &,
                                                     const char*)
{
  if( strcmp( file_type, "OCC" ) != 0 )
  {
     //PRINT_ERROR("The specified file type, %s, is not supported!\n", filetype );
     return CUBIT_FAILURE;
  }
 
  DLIList<OCCBody*>    OCC_bodies;
  DLIList<OCCSurface*> OCC_surfaces;
  DLIList<OCCCurve*>   OCC_curves;
  DLIList<OCCPoint*>   OCC_points;

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
      OCC_points.append( pt );

    //if it is a Curve
    else if( OCCCurve* curve = CAST_TO( ref_entity_ptr, OCCCurve) )
      OCC_curves.append( curve );
    
    //if it is a Surface -- I don't think you can ever have a free surface
    //without it being a Body
    else if( OCCSurface* surf = CAST_TO( ref_entity_ptr, OCCSurface) )
      OCC_surfaces.append( surf );
   
    //if it is a Body
    else if( OCCBody* body = CAST_TO( ref_entity_ptr, OCCBody ) )
      OCC_bodies.append( body );

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

  int free_body_count = OCC_bodies.size();
  int free_curve_count = OCC_curves.size();
  int free_point_count = OCC_points.size();
  int free_surface_count = OCC_surfaces.size();

  //if nothing to write out...return
  if( free_body_count == 0 && free_surface_count == 0 && 
      free_curve_count == 0 && free_point_count == 0)
    return CUBIT_SUCCESS;

  //save the facets (geometry info )
  CubitStatus status;

  //write out topology and attributes
  status = write_topology( file_name,
                           OCC_bodies, OCC_surfaces,
                           OCC_curves, OCC_points );
  if( status == CUBIT_FAILURE ) return CUBIT_FAILURE;

  if( free_body_count || free_surface_count || 
      free_curve_count || free_point_count )
      PRINT_INFO( "\nExported:" );

   int flg = 0;

   if( free_body_count )
   {
      if( flg )PRINT_INFO( "         " );else flg=1;
      if( DEBUG_FLAG( 153 ) )
      {
        if( free_body_count == 1 )
           PRINT_INFO( "%4d OCC Body\n", free_body_count );
        else
           PRINT_INFO( "%4d OCC Bodies\n", free_body_count );
      }
   }

   if( free_surface_count )
   {
      if( flg )PRINT_INFO( "         " );else flg=1;
      if( free_surface_count == 1 )
         PRINT_INFO( "%4d OCC Surface\n", free_surface_count );
      else
         PRINT_INFO( "%4d OCC Surface\n", free_surface_count );
   }

   if( free_curve_count )
   {
      if( flg )PRINT_INFO( "         " );else flg=1;
      if( free_curve_count == 1 )
         PRINT_INFO( "%4d OCC Curve\n", free_curve_count );
      else
         PRINT_INFO( "%4d OCC Curves\n", free_curve_count );
   }

   if( free_point_count )
   {
      if( flg )PRINT_INFO( "         " );else flg=1;
      if( free_point_count == 1 )
         PRINT_INFO( "%4d OCC Point\n", free_point_count );
      else
         PRINT_INFO( "%4d OCC Points\n", free_point_count );
   }
   PRINT_INFO( "\n" );

   return CUBIT_SUCCESS;
}

CubitStatus
OCCQueryEngine::write_topology( const char* file_name,
                                DLIList<OCCBody*> &OCC_bodies,
                                DLIList<OCCSurface*> &OCC_surfaces,
                                DLIList<OCCCurve*> &OCC_curves,
                                DLIList<OCCPoint*> &OCC_points )
{

  int i;
  //Create a compound shape to export
  BRep_Builder B;
  TopoDS_Compound Co;
  B.MakeCompound(Co);

  //Add every shape to the compound
  for (i = 0; i < OCC_bodies.size(); i++)
  {
     OCCBody* body = OCC_bodies.get_and_step();
     TopoDS_CompSolid *shape = body->get_TopoDS_Shape();

     //check if this body is build backwards from lump. if so,
     //the body and its CompSolid doesn't have bounded relationship
     //established. In this case, each individual lump of the body 
     // will be exported as TopoDS_Solid without a CompSolid
     if(OCCMap->IsBound(*shape))
       B.Add(Co, *shape);
     else
     {   
        DLIList<Lump*> lumps = body->lumps();
	for(int i = 0; i < lumps.size(); i++)
	{
	  OCCLump *occ_lump = (OCCLump *) lumps.get_and_step();
	  B.Add(Co, *(occ_lump->get_TopoDS_Solid()));
	}
     }
  }

  for (i = 0; i < OCC_surfaces.size(); i++)
  {
     TopoDS_Face *face = OCC_surfaces.get_and_step()->get_TopoDS_Face();
     B.Add(Co, *face);
  }

  for (i = 0; i < OCC_curves.size(); i++)
  {
     TopoDS_Edge *edge = OCC_curves.get_and_step()->get_TopoDS_Edge();
     B.Add(Co, *edge);
  }

  for (i = 0; i < OCC_points.size(); i++)
  {
     TopoDS_Vertex *vertex = OCC_points.get_and_step()->get_TopoDS_Vertex();
     B.Add(Co, *vertex);
  }
 
  char* file = new char[sizeof(file_name)];
  strcpy(file, file_name);
  
  if(!BRepTools::Write(Co, file))
    return CUBIT_FAILURE;
 
  return CUBIT_SUCCESS;
}


CubitStatus
OCCQueryEngine::import_temp_geom_file(FILE* file_ptr,
                                      const char* file_name,
                                      const char* file_type,
                                      DLIList<TopologyBridge*> &bridge_list )
{
  //make sure that file_type == "OCC"
  if( !strcmp( file_type,"OCC") )
    return import_solid_model( file_name, file_type, bridge_list );
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
  PRINT_RESULT = print_results;
  
  imported_entities = populate_topology_bridge(*aShape);
  return CUBIT_SUCCESS;
}

DLIList<TopologyBridge*> OCCQueryEngine::populate_topology_bridge(TopoDS_Shape aShape)
{
	DLIList<TopologyBridge*> tblist;
        // suitable to popolate for a TopoDS_CompSolid or TopoDS_Compound shape.
        TopExp_Explorer Ex;
        for (Ex.Init(aShape, TopAbs_COMPSOLID); Ex.More(); Ex.Next())
          tblist.append(populate_topology_bridge(TopoDS::CompSolid(Ex.Current())));

        for (Ex.Init(aShape, TopAbs_SOLID, TopAbs_COMPSOLID); Ex.More(); Ex.Next())
 	  tblist.append(populate_topology_bridge(TopoDS::Solid(Ex.Current()), CUBIT_TRUE));

	for (Ex.Init(aShape, TopAbs_SHELL, TopAbs_SOLID); Ex.More(); Ex.Next())
          tblist.append(populate_topology_bridge(TopoDS::Shell(Ex.Current())));

        for (Ex.Init(aShape, TopAbs_FACE, TopAbs_SHELL); Ex.More(); Ex.Next())
          tblist.append(populate_topology_bridge(TopoDS::Face(Ex.Current())));

	for (Ex.Init(aShape, TopAbs_WIRE, TopAbs_FACE); Ex.More(); Ex.Next())
          tblist.append(populate_topology_bridge(TopoDS::Wire(Ex.Current())));

	for (Ex.Init(aShape, TopAbs_EDGE, TopAbs_WIRE); Ex.More(); Ex.Next())
          tblist.append(populate_topology_bridge(TopoDS::Edge(Ex.Current())));

	for (Ex.Init(aShape, TopAbs_VERTEX, TopAbs_EDGE); Ex.More(); Ex.Next())
          tblist.append(populate_topology_bridge(TopoDS::Vertex(Ex.Current())));
	return tblist;
}

BodySM* OCCQueryEngine::populate_topology_bridge(TopoDS_CompSolid aShape)
{
        TopoDS_CompSolid *posolid =  new TopoDS_CompSolid;
        *posolid = aShape;
        OCCBody *body;
        if (!OCCMap->IsBound(*posolid))
        {
                if(PRINT_RESULT)
                        PRINT_INFO("Adding Bodies.\n");
                iTotalTBCreated++;
                body = new OCCBody(posolid);
                OCCMap->Bind(*posolid, iTotalTBCreated);
                OccToCGM->insert(valType(iTotalTBCreated,
                                (TopologyBridge*)body));
        }
        else
        {
                int k = OCCMap->Find(*posolid);
                body = (OCCBody*)(OccToCGM->find(k))->second;
        }
        TopExp_Explorer Ex;
        DLIList<Lump*> lumps;
        for (Ex.Init(aShape, TopAbs_SOLID); Ex.More(); Ex.Next())
        {
           Lump* lump = populate_topology_bridge(TopoDS::Solid(Ex.Current()));
  	   lumps.append(lump);
	}
	body->lumps(lumps);
        return body;
}

Lump* OCCQueryEngine::populate_topology_bridge(TopoDS_Solid aShape,
						 CubitBoolean build_body)
{
	//one OCCBody corresponds one OCCLump
	TopoDS_Solid *posolid =  new TopoDS_Solid;
 	*posolid = aShape;
	OCCLump *lump;
        OCCBody *body;
	if (!OCCMap->IsBound(*posolid))
 	{
		if(PRINT_RESULT)
			PRINT_INFO("Adding solids.\n");
                iTotalTBCreated++;
		lump = new OCCLump(posolid);
   		if (build_body)
		{
                   DLIList<Lump*> lumps;
		   lumps.append(lump);
		   body = new OCCBody(lumps);
		}
		OCCMap->Bind(*posolid, iTotalTBCreated);
                OccToCGM->insert(valType(iTotalTBCreated,
				(TopologyBridge*)lump));
	}
	else 
	{
		int k = OCCMap->Find(*posolid);
		lump = (OCCLump*)(OccToCGM->find(k))->second;
	}
	TopExp_Explorer Ex;
        for (Ex.Init(aShape, TopAbs_SHELL); Ex.More(); Ex.Next())
                populate_topology_bridge(TopoDS::Shell(Ex.Current()));
	return lump;
}

OCCShell* OCCQueryEngine::populate_topology_bridge(TopoDS_Shell aShape)
{
	TopoDS_Shell *poshell = new TopoDS_Shell;
	*poshell = aShape;
	OCCShell *shell ;
        if (!OCCMap->IsBound(*poshell))
        {
                if(PRINT_RESULT)
                        PRINT_INFO("Adding shells.\n");
                iTotalTBCreated++;
		shell = new OCCShell(poshell);
                OCCMap->Bind(*poshell, iTotalTBCreated);
                OccToCGM->insert(valType(iTotalTBCreated,
                                (TopologyBridge*)shell));
        }
        else
        {
                int k = OCCMap->Find(*poshell);
                shell = (OCCShell*)(OccToCGM->find(k))->second;
        }

	TopExp_Explorer Ex;
        DLIList<Surface *> surfaces;
        for (Ex.Init(aShape, TopAbs_FACE); Ex.More(); Ex.Next())
	{
           Surface* face = populate_topology_bridge(TopoDS::Face(Ex.Current()));
	   surfaces.append(face);
	}
	shell->add_surfaces(surfaces); 
	return shell;
}

Surface* OCCQueryEngine::populate_topology_bridge(TopoDS_Face aShape)
{
	TopoDS_Face *poface = new TopoDS_Face;
	*poface = aShape;
	OCCSurface *surface;
	if (!OCCMap->IsBound(*poface))
  	{
		if(PRINT_RESULT)
                        PRINT_INFO("Adding faces.\n");
        	iTotalTBCreated++;
		surface = new OCCSurface(poface);

		OCCMap->Bind(*poface, iTotalTBCreated);
  		OccToCGM->insert(valType(iTotalTBCreated,
				(TopologyBridge*)surface));
	} 
	else 
	{
                int k = OCCMap->Find(*poface);
		surface = (OCCSurface*)(OccToCGM->find(k))->second;
	}
	TopExp_Explorer Ex;
        for (Ex.Init(aShape, TopAbs_WIRE); Ex.More(); Ex.Next())
	 	populate_topology_bridge(TopoDS::Wire(Ex.Current()));
	return surface;
}

OCCLoop* OCCQueryEngine::populate_topology_bridge(TopoDS_Wire aShape)
{
	TopoDS_Wire *powire = new TopoDS_Wire;
	*powire = aShape;
	OCCLoop *loop ;
        if (!OCCMap->IsBound(*powire))
        {
                if(PRINT_RESULT)
                        PRINT_INFO("Adding loops.\n");
                iTotalTBCreated++;
        	loop = new OCCLoop(powire);
 		OCCMap->Bind(*powire, iTotalTBCreated);
                OccToCGM->insert(valType(iTotalTBCreated,
                                (TopologyBridge*)loop));
        }
        else
        {
                int k = OCCMap->Find(*powire);
                loop = (OCCLoop*)(OccToCGM->find(k))->second;
        }

	BRepTools_WireExplorer Ex;
        DLIList <OCCCoEdge*> coedges;
	for (Ex.Init(aShape); Ex.More(); Ex.Next())
 	{
	   Curve* curve = populate_topology_bridge(Ex.Current()) ;
           OCCCurve *occ_curve = CAST_TO(curve, OCCCurve);
           TopoDS_Edge *edge = occ_curve->get_TopoDS_Edge( );
  	   OCCCoEdge * coedge = new OCCCoEdge(edge, curve, loop, 
 	   (Ex.Orientation()== TopAbs_FORWARD ? CUBIT_FORWARD : CUBIT_REVERSED));
	   coedges.append(coedge);
	}
	loop->add_coedges(coedges);
	return loop;
}

Curve* OCCQueryEngine::populate_topology_bridge(TopoDS_Edge aShape)
{
        Curve *curve;
        TopoDS_Edge *poedge = new TopoDS_Edge;
	*poedge = aShape;
	if (!OCCMap->IsBound(*poedge)) 
        {
		if(PRINT_RESULT)
                        PRINT_INFO("Adding edges.\n");
                iTotalTBCreated++;
                curve = new OCCCurve(poedge);
		OCCMap->Bind(*poedge, iTotalTBCreated);
                OccToCGM->insert(valType(iTotalTBCreated,
				(TopologyBridge*)curve));
	}
        else 
        {
                int i = OCCMap->Find(*poedge);
                curve = (OCCCurve*)(OccToCGM->find(i))->second;
	}

        TopExp_Explorer Ex;
        for (Ex.Init(aShape, TopAbs_VERTEX); Ex.More(); Ex.Next())
          populate_topology_bridge(TopoDS::Vertex(Ex.Current()));
	return curve;
}

Point* OCCQueryEngine::populate_topology_bridge(TopoDS_Vertex aShape)
{
	OCCPoint *point;
	TopoDS_Vertex *povertex = new TopoDS_Vertex;
	*povertex = aShape;
	if (!OCCMap->IsBound(*povertex)) 
        {
 	  if(PRINT_RESULT)
                        PRINT_INFO("Adding vertices.\n");
          iTotalTBCreated++;
 	  point = new OCCPoint(povertex);
	  OCCMap->Bind(*povertex, iTotalTBCreated);
          OccToCGM->insert(valType(iTotalTBCreated,
                          (TopologyBridge*)point));
	} 
        else 
	{
          int i = OCCMap->Find(*povertex);
 	  point = (OCCPoint*)(OccToCGM->find(i))->second;
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
