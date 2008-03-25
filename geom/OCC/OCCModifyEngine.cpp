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
#include "gp_Ax2.hxx"
#include "gp_Dir.hxx"
#include "gp_Hypr.hxx"
#include "gp_Parab.hxx"
#include "gp_Elips.hxx"
#include "gp_Pln.hxx"
#include "gp_Cylinder.hxx"
#include "gp_Cone.hxx"
#include "gp_Sphere.hxx"
#include "gp_Torus.hxx"
#include "BRepBuilderAPI_MakeShell.hxx"
#include "BRepBuilderAPI_MakeSolid.hxx"
#include "BRepBuilderAPI_MakeWire.hxx"
#include "TopoDS_Shape.hxx"
#include "TColgp_Array1OfPnt.hxx"
#include "GC_MakeArcOfCircle.hxx"
#include "GC_MakeArcOfHyperbola.hxx"
#include "GC_MakeArcOfParabola.hxx"
#include "GC_MakeArcOfEllipse.hxx"
#include "GC_MakeSegment.hxx"
#include "GC_MakeTrimmedCone.hxx"
#include "GC_MakeTrimmedCylinder.hxx"
#include "Geom_BezierCurve.hxx"
#include "Handle_Geom_Plane.hxx"
#include "BRepPrimAPI_MakePrism.hxx"
#include "BRepPrimAPI_MakeCone.hxx"
#include "BRepPrimAPI_MakeTorus.hxx"
#include "GC_MakeEllipse.hxx"
#include "BRepBuilderAPI_MakeEdge.hxx"
#include "BRepAdaptor_Surface.hxx"
#include "BRepBuilderAPI_MakeFace.hxx"
#include "BRepBuilderAPI_Sewing.hxx"
#include "BRepBuilderAPI_Copy.hxx"
#include "BRep_Tool.hxx"
#include "BRep_Builder.hxx"
#include "GProp_GProps.hxx"
#include "BRepGProp.hxx"
#include "TopoDS.hxx"
#include "TopologyBridge.hpp"
#include "BRepAlgoAPI_Fuse.hxx"
#include "BRepAlgoAPI_Cut.hxx"
#include "BRepAlgoAPI_Section.hxx"
#include "BRepPrimAPI_MakeSphere.hxx"
#include "BRepPrimAPI_MakeBox.hxx"
#include "BRepPrimAPI_MakeWedge.hxx"
#include "Handle_Geom_TrimmedCurve.hxx"
#include "Handle_Geom_RectangularTrimmedSurface.hxx"
#include "TopExp_Explorer.hxx"
#include "OCCModifyEngine.hpp"
#include "OCCQueryEngine.hpp"
#include "CubitMessage.hpp"
#include "CubitDefines.h"
#include "TopTools_DataMapOfShapeInteger.hxx"
#include "BRepFeat_SplitShape.hxx"
#include "TopOpeBRep_ShapeIntersector.hxx"
#include "TopTools_ListIteratorOfListOfShape.hxx"
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
 
  BRepBuilderAPI_Copy api_copy(*theEdge);

  TopoDS_Shape newShape = api_copy.ModifiedShape(*theEdge);
 
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
       Point * Pnt = make_Point(*third_point);
       curve = make_Curve(type, point1_ptr, Pnt, mid_points);
    }
    else
    {
       PRINT_ERROR("Cannot create an OCC curve from the given duplicated points.\n");
       return (Curve *)NULL;
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
// Description: make a  spline curve by using the points on surface.
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
  
  if (curve_type != SPLINE_CURVE_TYPE
      && curve_type != STRAIGHT_CURVE_TYPE)
  {
     PRINT_ERROR("Cannot create an OCC curve from the given curve_type.\n"
                 "Candidates are SPLINE_CURVE_TYPE and STRAIGHT_CURVE_TYPE.\n");
     return (Curve *)NULL;
  }

  OCCPoint* occ_point1 = CAST_TO(const_cast<Point*>(point1_ptr), OCCPoint);
  OCCPoint* occ_point2 = CAST_TO(const_cast<Point*>(point2_ptr), OCCPoint);

  if (occ_point1 == NULL || occ_point2 == NULL)
  {
     PRINT_ERROR("Cannot create an OCC curve from the given points.\n"
                 "Possible incompatible geometry engines.\n");
     return (Curve *)NULL;
  }
    
  //project all points on the surface if possible
  OCCSurface* occ_face = NULL;
  if (face_ptr != NULL)
     occ_face = CAST_TO(face_ptr, OCCSurface);
 
  gp_Pnt pt;
  int size = 2+vector_list.size();
  TColgp_Array1OfPnt points(1, size);
  CubitVector* vector = NULL;
  CubitVector closest_location;
  for(int i = 1; i <= size; i++)
  {
     if (i == 1) 
     {
       TopoDS_Vertex *point = occ_point1->get_TopoDS_Vertex();
       pt = BRep_Tool::Pnt(*point);
       vector = new CubitVector(point1_ptr->coordinates());
     }
     else if (i == size)
     {
       TopoDS_Vertex *point = occ_point2->get_TopoDS_Vertex();
       pt = BRep_Tool::Pnt(*point);
       vector = new CubitVector(point2_ptr->coordinates()); 
     } 
     else
     {
       vector = vector_list.get_and_step();
       pt.SetCoord(vector->x(), vector->y(), vector->z());
     } 

     if (occ_face != NULL)
     {
       occ_face->closest_point(*vector, &closest_location);
       pt.SetCoord(closest_location.x(), closest_location.y(), closest_location.z()) ;  	 
     }

     points.SetValue(i, pt); 
  }    
     
  //make curve according to the curve type.
  if(curve_type == SPLINE_CURVE_TYPE)
  {
    Geom_BezierCurve BezierCurve(points);
    Handle(Geom_Curve) curve_ptr(&BezierCurve);
    TopoDS_Edge new_edge = BRepBuilderAPI_MakeEdge(curve_ptr);
    return OCCQueryEngine::instance()->populate_topology_bridge(new_edge); 
  }

  else if(curve_type == STRAIGHT_CURVE_TYPE)
  {
    TColgp_Array1OfPnt two_points(1,2); 
    two_points.SetValue(1, points.Value(1));
    two_points.SetValue(2, points.Value(size));
    Geom_BezierCurve BezierCurve(two_points);
    Handle(Geom_Curve) curve_ptr(&BezierCurve);
    TopoDS_Edge new_edge = BRepBuilderAPI_MakeEdge(curve_ptr);
    return OCCQueryEngine::instance()->populate_topology_bridge(new_edge); 
  }

  return (Curve*) NULL;
}

//===============================================================================
// Function   : make_Curve
// Member Type: PUBLIC
// Description: make a curve
// For STRAIGHT_CURVE_TYPE:
//    intermediate_point_ptr  is not used
//
// For PARABOLA_CURVE_TYPE
//    intermediate_point_ptr is the tip of the parabola
//
// For HYPERBOLA_CURVE_TYPE
//    intermediate_point_ptr is the center of its two foci
//
// For ELLIPSE_CURVE_TYPE
//    intermediate_point_ptr is the center of the ellipse
//    the point who is farther away to the center is the vertex of the ellipse
//    the others point projects to major axis at focus.
//    sense is used to determine which part of the ellipse is required
//
// For ARC_CURVE_TYPE
//    arc passes three points
//
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

  CubitVector v1(point1_ptr->coordinates());
  CubitVector v2(point2_ptr->coordinates());

  gp_Pnt pt1(v1.x(),v1.y(), v1.z());
  gp_Pnt pt2(v2.x(),v2.y(), v2.z());

  CubitVector v3;
  gp_Pnt pt3;
  double tol = OCCQueryEngine::instance()->get_sme_resabs_tolerance();

  Handle(Geom_TrimmedCurve) curve_ptr;
  if(intermediate_point_ptr != NULL)
  {
    v3 = *intermediate_point_ptr;
    pt3.SetCoord(v3.x(),v3.y(), v3.z());
  }

  if (curve_type == STRAIGHT_CURVE_TYPE)
     curve_ptr = GC_MakeSegment(pt1,pt2);

  else if (curve_type == ARC_CURVE_TYPE)
  {
     assert(intermediate_point_ptr != NULL);
     curve_ptr = GC_MakeArcOfCircle(pt1, pt2, pt3);
  }

  else if (curve_type == ELLIPSE_CURVE_TYPE)
  {
     assert(intermediate_point_ptr != NULL);
     
     //calculate for the axis
     double d1 = (v1 - v3).length_squared(); 
     double d2 = (v2 - v3).length_squared();

     CubitVector x = d1 >= d2 ? v1-v3 : v2-v3;
     x.normalize();
     gp_Dir x_dir(x.x(), x.y(), x.z());

     CubitVector N = (v1 - v3) * (v2 - v3); 
     if(N.length_squared() < tol * tol)
     {
       PRINT_ERROR("Cannot create an ellipse curve from the given points.\n"
                 "3 points are in the same line.\n");
       return (Curve *)NULL;
     }
     N.normalize();
     if (sense == CUBIT_REVERSED)
       N = -N;
     gp_Dir N_dir(N.x(), N.y(), N.z());

     gp_Pnt center(v3.x(), v3.y(), v3.z());
     gp_Ax2 axis(center, N_dir, x_dir); 

     //calculate for the major and minor radius.
     double major = d1 >= d2 ? sqrt(d1): sqrt(d2);
     double other_d = d1 >= d2 ? sqrt(d2) : sqrt(d1);
     double c = cos((v1 - v3).interior_angle(v2 - v3)) * other_d;
     double minor = sqrt(major * major - c * c);

     gp_Elips ellipse(axis, major, minor);
     CubitBoolean use_sense = (sense == CUBIT_FORWARD ? CUBIT_TRUE : CUBIT_FALSE); 
     curve_ptr = GC_MakeArcOfEllipse(ellipse, pt1, pt2, use_sense);
  }

  else if(curve_type == PARABOLA_CURVE_TYPE || 
          curve_type == HYPERBOLA_CURVE_TYPE)
  {
    assert(intermediate_point_ptr != NULL);

    //find the directrix and focus of the parabola
    //or the axis, major radius and minor radius of the hyperbola
    CubitVector width_vec = v2 - v1;
    if(width_vec.length_squared() < tol * tol)
    {
       PRINT_ERROR("Cannot create a parabola or hyperbola curve from the given points.\n"
                 "2 end points are the same.\n");
       return (Curve *)NULL;
    }

    CubitVector midpoint_vec = (v1 + v2)/2.0;
    CubitVector height_vec = midpoint_vec - v3;
    gp_Pnt center(v3.x(), v3.y(), v3.z());
 
    if (height_vec.length_squared() < tol * tol)
    { 
       PRINT_ERROR("Cannot create a parabola or hyperbola curve from the given points.\n"
                 "3 points are in the same line.\n");
       return (Curve *)NULL;
    }
    CubitVector x = height_vec;
    x.normalize();
    gp_Dir x_dir(x.x(), x.y(), x.z());
 
    CubitVector N = x * (v2 - v1);  
    if (N.length_squared() < tol * tol)
    {
       PRINT_ERROR("Cannot create a parabola or hyperbola curve from the given points.\n"
                 "3 points are in the same line.\n");
       return (Curve *)NULL;
    }
    N.normalize();
    gp_Dir N_dir(N.x(), N.y(), N.z());

    gp_Ax2 axis(center, N_dir, x_dir);  

    if(curve_type == HYPERBOLA_CURVE_TYPE)
    { 
       //    (focus2) (v3) . (v2)
       //          .   .   . (midpoint = focus1)
       //                  . (v1)
       CubitVector focus2 = 2 * v3 - midpoint_vec;

       //according to the definition of hyperbola,
       //2 * a = length(v2 - focus2)-length(v2 - focus1)

       double major = (v2 - focus2).length()/2.0 - (v2 - midpoint_vec).length()/2.0;

       // if a = 1/2 length major axis, b = 1/2 length minor axis and
       // c = distance center to focus, then a*a + b*b = c*c

       double c_squared = (midpoint_vec - v3).length_squared();
       double minor = sqrt(c_squared  - major*major );
       gp_Hypr hypt(axis, major, minor);
       curve_ptr =
             GC_MakeArcOfHyperbola(hypt, pt1, pt2, CUBIT_TRUE);
    }

    else
    {
       // Find the focus of this parabola.
       // Since for a parabola with its peak at the origin, y = (1/(4*a))*x^2,
       // and since we have restricted this parabola to be symmetric (per the 
       // FastQ method, see the FastQ file getwt.f), we can use the following 
       // relationship to
       // determine "a", the distance the focus lies from the peak on the line
       // formed by the peak and the midpoint of the start and end points`
       double a = width_vec.length_squared()/(16. * height_vec.length()); 
       gp_Parab parab(axis, a);
       curve_ptr =
		GC_MakeArcOfParabola(parab, pt1, pt2, CUBIT_TRUE);
    } 
  }

  else
  {
      PRINT_ERROR("In OCCModifyEngine::make_Curve\n"
                  "       Invalid curve type.\n");
      return (Curve *)NULL;
  }

  TopoDS_Edge new_edge = BRepBuilderAPI_MakeEdge(curve_ptr);
  return OCCQueryEngine::instance()->populate_topology_bridge(new_edge);
}

//===============================================================================
// Function   : make_Surface
//              This function creates a surface given an existing surface, copy.
// Member Type: PUBLIC
// Description: make a surface, OCC allows to create a stand along surface,
//              however the CGM has a design of making all surfaces in a (sheet)
//              body. This will add complexity on all free surface related 
//              calculation and modification, and adding potential bugs too. 
// Author     : Jane Hu
// Date       : 02/08
//===============================================================================
Surface* OCCModifyEngine::make_Surface( Surface * surface_ptr,
                                 CubitBoolean extended_from) const
{
  OCCSurface* occ_surface = CAST_TO(surface_ptr, OCCSurface);
  if (!occ_surface)
  {
     PRINT_ERROR("Cannot create an OCC surface from the given surface.\n"
                 "Possible incompatible geometry engines.\n");
     return (Surface *)NULL;
  }

  //Start of the codes
  double UMax, VMax, UMin, VMin;
  occ_surface->get_param_range_U(UMin, UMax);
  occ_surface->get_param_range_V(VMin, VMax);

  TopoDS_Face *theFace = occ_surface->get_TopoDS_Face();
  if( !theFace)
  {
     PRINT_ERROR("Cannot create an OCC surface from the given surface.\n"
                 "Possible incompatible geometry engines.\n");
     return (Surface *)NULL;
  }

  TopoDS_Face newFace;
  BRepAdaptor_Surface asurface(*theFace);

  CubitBox bounding_box = GeometryQueryTool::instance()->model_bounding_box();
  double const height = 2*(bounding_box.diagonal()).length();
  CubitBox box = occ_surface->bounding_box();
  double ratio = height/(box.diagonal().length());

  double middleU = (UMin + UMax)/2.0;
  double middleV = (VMin + VMax)/2.0;
  double U1 = middleU - (UMax-UMin)/2.0 * ratio;
  double U2 = middleU + (UMax-UMin)/2.0 * ratio;
  double V1 = middleV - (VMax - VMin)/2.0 * ratio;
  double V2 = middleV + (VMax - VMin)/2.0 * ratio;

  if (extended_from == CUBIT_TRUE)
  {
     // We need to get the type of surface.
     GeometryType type = occ_surface->geometry_type();
     if (type  == PLANE_SURFACE_TYPE)
     {
        gp_Pln plane = asurface.Plane();
        newFace = BRepBuilderAPI_MakeFace(plane, U1, U2, V1, V2);
     }
     else if(type == CONE_SURFACE_TYPE)
     {
       //make an infinite cone.
       //Given this lets create another face that is extended from it.
       if(asurface.GetType() == GeomAbs_Cone)
       {
         gp_Cone cone = asurface.Cone();
         gp_Pnt Apex = cone.Apex();
         double semi_angle = cone.SemiAngle();
         gp_Pnt p2;
         double radius2;
         gp_XYZ xyz;
         if (semi_angle > 0)
           xyz = Apex.XYZ() + cone.Position().Direction().XYZ()*height;
             
  	 else
	   xyz = Apex.XYZ() - cone.Position().Direction().XYZ()*height;

         p2.SetXYZ(xyz);
	 radius2 = height * tan(fabs(semi_angle));
         Handle(Geom_RectangularTrimmedSurface) trimmed_cone;
         trimmed_cone = GC_MakeTrimmedCone(Apex, p2, 0, radius2); 
         newFace = BRepBuilderAPI_MakeFace(trimmed_cone);
       }
       else
       {
         gp_Cylinder cylinder = asurface.Cylinder();
         //gp_Pnt pnt = asurface.Value(UMin,V1); 
         double radius = cylinder.Radius();
         gp_Ax1 axis = cylinder.Axis(); 
         Handle(Geom_RectangularTrimmedSurface) trimmed_cyl;
         trimmed_cyl = GC_MakeTrimmedCylinder(axis, radius, height);
         newFace = BRepBuilderAPI_MakeFace(trimmed_cyl);
       } 
     }
     else if(type == SPHERE_SURFACE_TYPE)
     {
       //make a whole sphere.
       gp_Sphere sphere = asurface.Sphere();
       newFace = BRepBuilderAPI_MakeFace(sphere); 
     }
     else if(type == TORUS_SURFACE_TYPE)
     {
       //make a whole torus
       gp_Torus torus = asurface.Torus();
       newFace = BRepBuilderAPI_MakeFace(torus);
     }
     else if(type == SPLINE_SURFACE_TYPE ) 
     {
       //extend the surfaces using the equation if possible.
       Handle(Geom_BezierSurface) bezier = asurface.Bezier();
       newFace = BRepBuilderAPI_MakeFace(bezier, U1, U2, V1, V2);
     }
  }
 
  else
  {
    BRepBuilderAPI_Copy api_copy(*theFace);
    TopoDS_Shape newShape = api_copy.ModifiedShape(*theFace);
    newFace = TopoDS::Face(newShape);
  }
  
  Surface *surface = OCCQueryEngine::instance()->populate_topology_bridge(
                               newFace, CUBIT_TRUE);

  return surface;
}

//===============================================================================
// Function   : make_Surface
// Member Type: PUBLIC
// Description: make a surface of type surface_type, given the list of curves.
//              check edges option is done in GeometryModifyTool level, so 
//              disregard this option.
// Author     : Jane Hu
// Date       : 02/08
//===============================================================================
Surface* OCCModifyEngine::make_Surface( GeometryType surface_type,
                                 DLIList<Curve*>& curve_list,
                                 Surface * old_surface_ptr,
                                 bool check_edges) const
{
  //Create TopoDS_Edge list to make a surface.
  DLIList<DLIList<TopoDS_Edge*>*> topo_edges_loops;
  curve_list.reset() ;
    
  //check no intersections of the TopoDS_Edge's.
  //need to check that no intersection in the middle of the curves, not at
  //vertices or out of boundary.

  int count = 0; //intersection point should be the same as curve_list size.
  for ( int i = 0 ; i < curve_list.size()-1 ; i++ )
  {
     for(int j = i+1; j < curve_list.size(); j ++)
     {
        DLIList<CubitVector*> intscts;
 	CubitBoolean bounded = CUBIT_TRUE;//dummy arg.
	CubitBoolean closest = CUBIT_TRUE;//dummy arg.
        CubitStatus yes_int = 
              OCCQueryEngine::instance()->get_intersections(curve_list[i],
				curve_list[j], intscts, bounded, closest);
        if(yes_int)
        {
           //check intscts point should be vertex or outside boundary.
 	   if (intscts.size() > 2 )  
	   {
	     PRINT_ERROR("In OCCModifyEngine::make_Surface\n"
                 "       Cannot make Surface with intersecting curves.\n");
             return (Surface *)NULL;
           }
           else
           {
             for(int k = 0; k < intscts.size(); k++)
             {
               CubitVector *v = intscts.get_and_step();
	       CubitPointContainment is_on = CAST_TO(curve_list[i],OCCCurve)->
					point_containment(*v);
               if (is_on == CUBIT_PNT_BOUNDARY)
               {
	 	 is_on = CAST_TO(curve_list[j],OCCCurve)->
				point_containment(*v);
		 if (is_on == CUBIT_PNT_BOUNDARY)
                   count++;
               }
               else if(is_on == CUBIT_PNT_INSIDE)
               {
                 PRINT_ERROR("In OCCModifyEngine::make_Surface\n"
                 "       Cannot make Surface with intersecting curves.\n");
                 return (Surface *)NULL;
               }
	     }
	   }
        }
     }
  }
 
  if (count > curve_list.size()) 
  {
      PRINT_ERROR("In OCCModifyEngine::make_Surface\n"
                "       Cannot make Surface with intersecting curves.\n");
      return (Surface *)NULL;
  }

  CubitStatus stat = sort_curves(curve_list, topo_edges_loops); 
  if( stat == CUBIT_FAILURE ) //case of one disconnected curve 
     return (Surface*) NULL;
 
  // Use the topo_edges to make a topo_face
  const TopoDS_Face* topo_face = make_TopoDS_Face(surface_type,
					topo_edges_loops, old_surface_ptr) ;
 
  if(topo_face == NULL)
  {
     PRINT_ERROR("In OCCModifyEngine::make_Surface\n"
                 "       Cannot make Surface object.\n");
     return (Surface *)NULL;
  }

  // make the topology bridges for the face
  TopoDS_Face the_face = *topo_face;
  Surface *surface = OCCQueryEngine::instance()->populate_topology_bridge(
                               the_face, CUBIT_TRUE); 
  return surface ;
}

//===============================================================================
// Function   : sort_curves
// Member Type: PROTECTED
// Description: sort the curves so they are in order and make closed loop 
// Author     : Jane Hu
// Date       : 03/08
//===============================================================================
CubitStatus OCCModifyEngine::sort_curves(DLIList<Curve*> curve_list,
                                         DLIList<DLIList<TopoDS_Edge*>*>& topo_edges_loops)const
{
  topo_edges_loops.clean_out();

  DLIList<TopoDS_Edge*>* topo_edges[curve_list.size()];
  for(int i = 0; i < curve_list.size(); i++)
    topo_edges[i] = new DLIList<TopoDS_Edge*>;

  curve_list.reset() ;
  Curve const* curve_ptr = NULL ;
  OCCCurve* occ_curve = NULL;
  TopoDS_Edge* topo_edge = NULL;

  OCCPoint* start = NULL;
  OCCPoint* end = NULL;
  DLIList<OCCPoint*> point_list;
  double tol = OCCQueryEngine::instance()->get_sme_resabs_tolerance();
  CubitBoolean new_end = CUBIT_TRUE;
  int size = curve_list.size();

  int count = 0;
  for ( int i = 0 ; i < size ; i++ )
  {
     for(int j = 0; j < curve_list.size(); j ++)
     {
        curve_ptr = curve_list.get() ; 
        occ_curve = CAST_TO(const_cast<Curve*>(curve_ptr), OCCCurve);

        if(occ_curve ==  NULL)
        {
           PRINT_ERROR("In OCCModifyEngine::sort_curves\n"
                       "       Got a NULL pointer to OCCCurve\n") ;
           return CUBIT_FAILURE;
        }

        point_list.clean_out();
        occ_curve->get_points(point_list);
        //assert(point_list.size()==2);

        if (i == 0)
        {
          start = point_list.get();
          end = point_list.pop();  
          break;
        }

        if(end->is_equal(*(point_list.get()), tol) ||
                end->is_equal(*(point_list.step_and_get()),tol))
        {
           end = point_list.step_and_get();
           new_end = CUBIT_TRUE;
           break;
        }
        curve_list.step();
     }

     if (new_end)//found next curve 
     {
        topo_edge = occ_curve->get_TopoDS_Edge();
        topo_edges[count]->append(topo_edge);
        curve_list.remove();
        if(start->is_equal( *end, tol))  //formed a closed loop
        {
          i = 0;
          size = curve_list.size() ;
          topo_edges_loops.append(topo_edges[count]);
          count++;
        }
        else
          new_end = CUBIT_FALSE;
     }
     else
     {
        PRINT_ERROR("In OCCModifyEngine::sort_curves\n"
                    "  Curve list  can't  form closed loops    \n") ;
        return CUBIT_FAILURE;
     }
  }

  if( new_end == CUBIT_FALSE ) //case of one disconnected curve
  {
     PRINT_ERROR("In OCCModifyEngine::sort_curves\n"
                 "  Curve list  can't  form closed loops    \n") ;
     return CUBIT_FAILURE;
  }
  return CUBIT_SUCCESS;
} 
//===============================================================================
// Function   : make_TopoDS_Face
// Member Type: PROTECTED
// Description: make a opoDS_Face of type surface_type, given the list of 
//              TopoDS_Edge. the TopoDS_Edge's should be in order in loops.
//              check edges option is done in GeometryModifyTool level, so
//              disregard this option.
// Author     : Jane Hu
// Date       : 02/08
//===============================================================================
const TopoDS_Face* OCCModifyEngine::make_TopoDS_Face(GeometryType surface_type,
			       DLIList<DLIList<TopoDS_Edge*>*> topo_edges_list,
			       Surface * old_surface_ptr)const
{
  // Make sure a supported type of surface is being requested.
  if ( surface_type != PLANE_SURFACE_TYPE  &&
       surface_type != BEST_FIT_SURFACE_TYPE)
  {
      PRINT_WARNING("In OCCGeometryEngine::make_TopoDS_Face\n"
                    "       At this time, cannot make a TopoDS_Face that isn't"
                    " planar or best fit.\n");
      return (TopoDS_Face *)NULL;
  }
 
  // Set the TopoDS_Face pointer, if requested.
  TopoDS_Face *fit_Face = NULL;
  Handle_Geom_Surface S;
  if ( old_surface_ptr != NULL )
  {
      OCCSurface *surf = CAST_TO(old_surface_ptr, OCCSurface );
      fit_Face = surf->get_TopoDS_Face();
      S = BRep_Tool::Surface(*fit_Face);
  }
 
  // Make a wire from the topo_edges.
  // Coincident TopoDS_Vertex will be deleted by OCC.
  if(topo_edges_list.size() == 0)
      return (TopoDS_Face*) NULL;

  DLIList<TopoDS_Wire*> wires;
  GProp_GProps myProps;
  double max_area  = 0.0;
  TopoDS_Wire* out_Wire = NULL;
  TopoDS_Wire test_Wire;

  DLIList<TopoDS_Edge*>* topo_edges; 
  //check and make sure the outer loop is in the first
  for(int i = 0; i < topo_edges_list.size() ; i++)
  {
    topo_edges = topo_edges_list.get_and_step();
    BRepBuilderAPI_MakeWire aWire(*(topo_edges->get()));
    for(int j = 1; j < topo_edges->size(); j++)
      aWire.Add(*(topo_edges->step_and_get()));

    test_Wire = aWire.Wire();
    wires.append(&test_Wire);
   
    if (topo_edges_list.size() == 1)
      break;

    BRepBuilderAPI_MakeFace made_face(test_Wire);

    if (!made_face.IsDone())
    {
       PRINT_ERROR("In OCCModifyEngine::make_TopoDS_Face\n"
                   "   Cannot find the best fit surface for given curves.\n");
       return (TopoDS_Face *)NULL;
    }
    TopoDS_Face test_face = made_face.Face();
    BRepGProp::SurfaceProperties(test_face, myProps); 
    double area = myProps.Mass();
    out_Wire = max_area > area ? out_Wire : &test_Wire;
    max_area = max_area > area ? max_area : area;
  } 

  if (out_Wire)
  {
    wires.remove(out_Wire);
    wires.insert_first(out_Wire);
  }

  //create the TopoDS_Face
  const TopoDS_Face* topo_face = NULL;
  CubitBoolean error = CUBIT_FALSE;

  for(int i = 0; i < topo_edges_list.size() ; i++)
  {
    TopoDS_Wire *the_wire = wires.get_and_step();
    if (i == 0)
    {
      if( old_surface_ptr != NULL )
      {
        BRepBuilderAPI_MakeFace made_face(S, *the_wire);
        if (!made_face.IsDone())
        {
          error = CUBIT_TRUE;
          break;
        }
        topo_face = new TopoDS_Face(made_face.Face());
      }
      else
      {
        CubitBoolean is_planar = (surface_type == PLANE_SURFACE_TYPE ?
				  CUBIT_TRUE : CUBIT_FALSE); 
        BRepBuilderAPI_MakeFace made_face(*the_wire, is_planar);
        if (!made_face.IsDone())
        {
          error = CUBIT_TRUE;
          break;
        }

        topo_face = new TopoDS_Face(made_face.Face());
      }
    }
    else
    {
      BRepBuilderAPI_MakeFace made_face(*topo_face, *the_wire);
      if (!made_face.IsDone())
      {
        error = CUBIT_TRUE;
        break;
      }

      delete topo_face;
      topo_face = new TopoDS_Face(made_face.Face());
    }
  } 

  if(error)
  {
    PRINT_ERROR("In OCCModifyEngine::make_TopoDS_Face\n"
                 "   Cannot find the best fit surface for given curves.\n");
    return (TopoDS_Face *)NULL;
  }

  return topo_face;
}
//===============================================================================
// Function   : make_Lump
// Member Type: PUBLIC
// Description: make a lump of one shell
// Author     : Jane Hu
// Date       : 02/08
//===============================================================================
Lump* OCCModifyEngine::make_Lump( DLIList<Surface*>& surface_list ) const
{
  if (surface_list.size() == 0) 
    return (Lump*) NULL;

  Surface* surface = surface_list.get();
  OCCSurface* occ_surface = CAST_TO(surface, OCCSurface);
  if(!occ_surface)
  {
     PRINT_ERROR("Cannot create an OCC Lump from the given surfaces.\n"
                 "Possible incompatible geometry engines.\n");
     return (Lump *)NULL;
  }
  TopoDS_Face* face_ptr = occ_surface->get_TopoDS_Face();

  for(int i = 1; i < surface_list.size(); i++)
  {
     Surface* surface = surface_list.step_and_get();
     OCCSurface* occ_surface = CAST_TO(surface, OCCSurface);
     if(!occ_surface)
     {
        PRINT_ERROR("Cannot create an OCC Lump from the given surfaces.\n"
                    "Possible incompatible geometry engines.\n");
        return (Lump *)NULL;
     }
     TopoDS_Face* face2_ptr = occ_surface->get_TopoDS_Face();
    
     BRepAlgoAPI_Fuse  aFuse(*face_ptr, *face2_ptr);
     if(aFuse.ErrorStatus() == 1) //The Object is created but Nothing is Done
     {
        PRINT_ERROR("Cannot create an OCC Lump from the given surfaces.\n"
                    "Check surface connectivities.\n");
        return (Lump *)NULL;
     }
     else if (aFuse.ErrorStatus() != 0)
     {
        PRINT_ERROR("Cannot create an OCC Lump from the given surfaces.\n"
                    "OCC internal error.\n");
        return (Lump *)NULL;
     }
     
     TopTools_ListOfShape shapes;
     shapes.Assign( aFuse.Modified(*face_ptr));
     face_ptr = & (TopoDS::Face(shapes.First()));
  }    
  
  Handle_Geom_Surface HGeom_surface = BRep_Tool::Surface(*face_ptr);
  BRepBuilderAPI_MakeShell aMakeShell(HGeom_surface);
  if(!aMakeShell.IsDone())
  {
     PRINT_ERROR("Cannot create an OCC Lump from the given surfaces.\n"
                 "OCC internal error.\n");
     return (Lump *)NULL;
  }

  TopoDS_Shell aShell = aMakeShell.Shell();
  BRepBuilderAPI_MakeSolid aMakeSolid(aShell);
  if (!aMakeSolid.IsDone())
  {
     PRINT_ERROR("Cannot create an OCC Lump from the given surfaces.\n"
                 "OCC internal error.\n");
     return (Lump *)NULL;
  }

  TopoDS_Solid aSolid = aMakeSolid.Solid();

  return 
    OCCQueryEngine::instance()->populate_topology_bridge(aSolid, CUBIT_TRUE); 
}

//===============================================================================
// Function   : make_BodySM
// Member Type: PUBLIC
// Description: make a BodySM from a surface
// Author     : Jane Hu
// Date       : 02/08
//===============================================================================
BodySM* OCCModifyEngine::make_BodySM( Surface *surface ) const
{
  OCCSurface* occ_surface = CAST_TO(surface, OCCSurface);
  if(!occ_surface)
  {
     PRINT_ERROR("Cannot create an OCC body from the given surface.\n"
                 "Possible incompatible geometry engines.\n");
     return (BodySM *)NULL;
  }

  OCCBody* occ_body = occ_surface->my_body();
  if(occ_body)
     return occ_body;

  TopoDS_Face* face = occ_surface->get_TopoDS_Face();
  surface = OCCQueryEngine::instance()->populate_topology_bridge(*face, CUBIT_TRUE);
   
  return CAST_TO(surface, OCCSurface)->my_body();
}



//===============================================================================
// Function   : make_BodySM
// Member Type: PUBLIC
// Description: make a BodySM given a list of Lumps.
// Author     : Jane Hu
// Date       : 02/08
//===============================================================================
BodySM* OCCModifyEngine::make_BodySM( DLIList<Lump*>& lump_list ) const
{
  if (lump_list.size() == 0)
    return (BodySM*) NULL;

  //Create a compsolid shape, save all BodySM's correponding to lump_list
  //for deletion.
  DLIList<BodySM*> bodysm_list;
  TopoDS_CompSolid CS;
  BRep_Builder B;
  B.MakeCompSolid(CS);

  //Add every shape to the CompSolid
  for(int i = 0; i < lump_list.size(); i++)
  {
     Lump* lump = lump_list.get_and_step();
     OCCLump* occ_lump = CAST_TO(lump, OCCLump);
     if(!occ_lump)
     {
        PRINT_ERROR("Cannot create an OCC BodySM from the given lumps.\n"
                    "Possible incompatible geometry engines.\n");
        return (BodySM *)NULL;
     }
     TopoDS_Solid* solid = occ_lump->get_TopoDS_Solid();
     B.Add(CS, *solid);

     BodySM* bodysm_ptr = occ_lump->get_body();
     if(bodysm_ptr == NULL)
          continue;

     bodysm_list.append_unique(bodysm_ptr);
  }
 
  BodySM* bodysm = OCCQueryEngine::instance()->populate_topology_bridge(CS);

  if(bodysm)
  {
     //remove each Lump's body from the BodyList
     for(int i = 0; i < bodysm_list.size(); i++)
     {
        BodySM* bodysm_ptr = bodysm_list.get_and_step();

        OCCQueryEngine::instance()->unhook_BodySM_from_OCC(bodysm_ptr);
     }
  } 
  return bodysm;

}


//===============================================================================
// Function   : sphere
// Member Type: PUBLIC
// Description: build an OCC sphere
// Author     : Jane Hu
// Date       : 03/08
//===============================================================================
BodySM* OCCModifyEngine::sphere(double radius) const
{
  if (radius <= 0)
    return (BodySM*) NULL; 

  TopoDS_Solid S = BRepPrimAPI_MakeSphere(radius);
  
  Lump* lump = OCCQueryEngine::instance()->populate_topology_bridge(S, 
								CUBIT_TRUE);
  if (lump == NULL)
    return (BodySM*)NULL;

  return CAST_TO(lump, OCCLump)->body();
}


//===============================================================================
// Function   : brick
// Member Type: PUBLIC
// Description: build an OCC brick 
// Author     : Jane Hu
// Date       : 03/08
//===============================================================================
BodySM* OCCModifyEngine::brick( double wid, double dep, double hi ) const
{
  if (wid <= 0 || dep <=0 || hi <= 0)
    return (BodySM*)NULL;
  
  TopoDS_Solid S = BRepPrimAPI_MakeBox(wid, dep, hi);

  Lump* lump = OCCQueryEngine::instance()->populate_topology_bridge(S,
								CUBIT_TRUE);

  if (lump == NULL)
    return (BodySM*)NULL;

  return CAST_TO(lump, OCCLump)->body();
}


//===============================================================================
// Function   : brick
// Member Type: PUBLIC
// Description: create an OCC brick given center axes and extension
// Author     : Jane Hu
// Date       : 03/08
//===============================================================================
BodySM* OCCModifyEngine::brick( const CubitVector& center, 
                                const CubitVector axes[3],
                                const CubitVector &extension) const
{
  CubitVector left_corner =  center - 0.5 * extension;
  gp_Pnt left_point(left_corner.x(), left_corner.y(), left_corner.z());
  gp_Dir main_dir(axes[2].x(), axes[2].y(), axes[2].z());
  gp_Dir x_dir(axes[0].x(), axes[0].y(), axes[0].z());
  gp_Ax2 Axis(left_point, main_dir, x_dir);
  TopoDS_Solid S = BRepPrimAPI_MakeBox( Axis, extension.x(), extension.y(),
					extension.z());

  Lump* lump =  OCCQueryEngine::instance()->populate_topology_bridge(S,
								CUBIT_TRUE);
  if (lump == NULL)
    return (BodySM*)NULL;

  return CAST_TO(lump, OCCLump)->body(); 
}

//===============================================================================
// Function   : prism
// Member Type: PUBLIC
// Description: create an OCC prism 
// Author     : Jane Hu  
// Date       : 03/08
//===============================================================================
BodySM* OCCModifyEngine::prism( double height, int sides, double major,
                               double minor) const
{
  //currently OCC only support 4 sided prism
  if (sides != 4)
  {
    PRINT_ERROR("Option not supported for OCC based geometry.\n");
    return (BodySM*) NULL;
  }

  return brick(2 * major, 2 * minor, height); 
}

//===============================================================================
// Function   : pyramid
// Member Type: PUBLIC
// Description: create an OCC pyramid 
// Author     : Jane Hu
// Date       : 03/08
//===============================================================================
BodySM* OCCModifyEngine::pyramid( double height, int sides, double major,
                                 double minor, double top) const
{
  //currently OCC only support 4 sided pyramid
  if (sides != 4)
  {
    PRINT_ERROR("Option not supported for OCC based geometry.\n");
    return (BodySM*) NULL;
  }
  
  TopoDS_Solid S = BRepPrimAPI_MakeWedge( 2 * major, 2* minor, height,
                                        2 * top);

  Lump* lump =  OCCQueryEngine::instance()->populate_topology_bridge(S,
                                                                CUBIT_TRUE);
  if (lump == NULL)
    return (BodySM*)NULL;

  return CAST_TO(lump, OCCLump)->body();
  
}

//===============================================================================
// Function   : cylinder
// Member Type: PUBLIC
// Description: create an OCC cylinder, its base shape can be ellipse with
//		r1, r2 while r3 = 0; or it can be a cone with r1, 
//		r3 while r2 = 0. 
// Author     : Jane Hu
// Date       : 03/08
//===============================================================================
BodySM* OCCModifyEngine::cylinder( double hi, double r1, double r2, double r3 ) const
{
  if(r2 > 0 && r3 > 0)
  {
    PRINT_WARNING("Can not make elliptical cone for OCC engine.\n");
    return (BodySM*) NULL;
  }

  TopoDS_Solid S;
  if(r3 == 0)//elliptical based cylinder
  {
    gp_Pnt center(0.0, 0.0, 0.0);
    gp_Dir main_dir(0.0, 0.0, 1.0);
    gp_Dir x_dir(1.0, 0.0, 0.0);
    gp_Ax2 Axis(center, main_dir, x_dir); 
    Handle(Geom_Curve) curve_ptr = GC_MakeEllipse(Axis, r1, r2); 
    TopoDS_Edge new_edge = BRepBuilderAPI_MakeEdge(curve_ptr);
    BRepBuilderAPI_MakeWire aWire(new_edge);

    TopoDS_Wire test_Wire = aWire.Wire();

    BRepBuilderAPI_MakeFace made_face(test_Wire);

    if (!made_face.IsDone())
    {
       PRINT_ERROR("In OCCModifyEngine::cylinder\n"
                   "   Cannot create elliptical surface for given radii.\n");
       return (BodySM *)NULL;
    }
    TopoDS_Face test_face = made_face.Face();    
    gp_Vec V(0.0, 0.0, hi);
    TopoDS_Shape S1 = BRepPrimAPI_MakePrism(test_face, V);    
    S = TopoDS::Solid(S1);
  }

  else // cone
    S = BRepPrimAPI_MakeCone(r1, r3, hi);

  Lump* lump = OCCQueryEngine::instance()->populate_topology_bridge(S,
                                                                CUBIT_TRUE);

  if (lump == NULL)
  {
    PRINT_ERROR("In OCCModifyEngine::cylinder\n"
                "   Cannot create a cylinder for given radii.\n");
    return (BodySM*)NULL;
  }

  return CAST_TO(lump, OCCLump)->body();
}

//===============================================================================
// Function   : torus
// Member Type: PUBLIC
// Description: create an OCC torus
// Author     : Jane Hu
// Date       : 03/08
//===============================================================================
BodySM* OCCModifyEngine::torus( double r1, double r2 ) const
{
  if (r1 <= 0 || r2 <= 0)
    return (BodySM*) NULL;
 
  TopoDS_Solid S = BRepPrimAPI_MakeTorus(r1, r2);

  Lump* lump = OCCQueryEngine::instance()->populate_topology_bridge(S,
                                                                CUBIT_TRUE);

  if (lump == NULL)
  {
    PRINT_ERROR("In OCCModifyEngine::torus\n"
                "   Cannot create a torus for given radii.\n");
    return (BodySM*)NULL;
  }

  return CAST_TO(lump, OCCLump)->body();
}

//===============================================================================
// Function   : planar_sheet
// Member Type: PUBLIC
// Description: create an OCC planar_sheet with four vectors. 
// Author     : Jane Hu
// Date       : 03/08
//===============================================================================
BodySM* OCCModifyEngine::planar_sheet ( const CubitVector& p1,
                                       const CubitVector& p2,
                                       const CubitVector& p3,
                                       const CubitVector& p4) const
{
  Point* point1 = make_Point(p1);
  Point* point2 = make_Point(p2);
  Point* point3 = make_Point(p3);
  Point* point4 = make_Point(p4);
  Curve * curve1 = make_Curve( point1, point2);
  if (curve1 == NULL)
	return (BodySM*) NULL;
  Curve * curve2 = make_Curve( point2, point3); 
  if (curve2 == NULL)
        return (BodySM*) NULL;
  Curve * curve3 = make_Curve( point3, point4);
  if (curve3 == NULL)
        return (BodySM*) NULL;
  Curve * curve4 = make_Curve( point4, point1);
  if (curve4 == NULL)
        return (BodySM*) NULL;
  DLIList<Curve*> curves;
  curves.append(curve1);
  curves.append(curve2);
  curves.append(curve3);
  curves.append(curve4);
  Surface* surface = make_Surface(PLANE_SURFACE_TYPE, curves);
  if (surface == NULL)
	return (BodySM*) NULL;
  return CAST_TO(surface,OCCSurface)->my_body();
}

//===============================================================================
// Function   : copy_body
// Member Type: PUBLIC
// Description: copy an OCC-based body
// Author     : Jane Hu
// Date       : 03/08
//===============================================================================
BodySM* OCCModifyEngine::copy_body ( BodySM* bodyPtr ) const
{
  OCCBody* occ_body = CAST_TO(bodyPtr, OCCBody);
  if (!occ_body)
  {
     PRINT_ERROR("Cannot create an OCC bodySM from the given bodySM.\n"
                 "Possible incompatible geometry engines.\n");
     return (BodySM *)NULL;
  }

  TopoDS_CompSolid *theCS = occ_body->get_TopoDS_Shape();
  
  if (theCS == NULL) //sheet body
  {
    Surface* surface = make_Surface(occ_body->my_sheet_surface());
    if (surface == NULL)
    {
       PRINT_ERROR("Cannot create an OCC sheet bodySM from the given bodySM.\n");
       return (BodySM *)NULL;
    }

    return CAST_TO(surface,OCCSurface)->my_body();
  }

  if(OCCQueryEngine::instance()->OCCMap->IsBound(*theCS))
  {
    BRepBuilderAPI_Copy api_copy(*theCS);

    TopoDS_Shape newShape = api_copy.ModifiedShape(*theCS);

    TopoDS_CompSolid newCS = TopoDS::CompSolid(newShape);

    return OCCQueryEngine::instance()->populate_topology_bridge(newCS);
  }

  //single lump body
  Lump *lump = occ_body->lumps().get();
  TopoDS_Solid solid = *(CAST_TO(lump, OCCLump)->get_TopoDS_Solid());
  BRepBuilderAPI_Copy api_copy(solid);
  TopoDS_Shape newShape = api_copy.ModifiedShape(solid);
  TopoDS_Solid newSolid = TopoDS::Solid(newShape);
  lump = OCCQueryEngine::instance()->populate_topology_bridge(newSolid,
							CUBIT_TRUE);
  
  return CAST_TO(lump, OCCLump)->body();
}

//===============================================================================
// Function   : stitch_surfs
// Member Type: PUBLIC
// Description: stitch all surfs and try to make a shell body.
// Author     : Jane Hu
// Date       : 03/08
//===============================================================================
CubitStatus OCCModifyEngine::stitch_surfs(
                      DLIList<BodySM*>& surf_bodies,
                      BodySM*& stitched_body) const
{
  stitched_body = NULL;
  if (surf_bodies.size()==0)
    return CUBIT_SUCCESS;

  if (surf_bodies.size()==1)
  {
    stitched_body = surf_bodies.get();
    return CUBIT_SUCCESS;
  }

  DLIList<TopoDS_Face*> faces_to_stitch;
  for (int i = 0; i < surf_bodies.size(); i++)
  {
     BodySM * tool_body = surf_bodies.get_and_step();
     OCCBody* occ_body = CAST_TO(tool_body, OCCBody);
     OCCSurface* surface = occ_body->my_sheet_surface();
     if (surface == NULL)
     {
       PRINT_ERROR("Can't stitch non-sheet bodySM's. \n");
       return CUBIT_FAILURE;
     } 
     TopoDS_Face* topods_face = surface->get_TopoDS_Face();
     if (topods_face != NULL)
       faces_to_stitch.append(topods_face);
  }

  TopoDS_Shape* first_face  = faces_to_stitch.get();
  TopoDS_Face* second_face = NULL;
  TopoDS_Shape fuse;
  for(int i = 1; i < faces_to_stitch.size(); i++)
  {
     second_face = faces_to_stitch[i];
     fuse = BRepAlgoAPI_Fuse(*first_face, *second_face);
     first_face = &fuse;
  }

  TopExp_Explorer Ex;
  int count_face = 0;
  int count_shell = 0;
  for (Ex.Init(fuse, TopAbs_FACE, TopAbs_SHELL); Ex.More(); Ex.Next())
    count_face++;
  for (Ex.Init(fuse, TopAbs_SHELL, TopAbs_SOLID); Ex.More(); Ex.Next())
    count_shell++;

  if (count_face != 1 && count_shell != 1)
  {
     PRINT_ERROR("Can't stitch all surfaces into one BodySM's. \n");
     return CUBIT_FAILURE;
  }
 
  DLIList<TopologyBridge*> tbs = OCCQueryEngine::instance()->
     populate_topology_bridge(fuse);
  OCCBody* body = CAST_TO(tbs.get(), OCCBody);
  if (body)
    stitched_body = body ;

  else
  {
    OCCSurface* face = CAST_TO(tbs.get(), OCCSurface);
    if(face)
      stitched_body = face->my_body();
  }

  return CUBIT_SUCCESS; 
}

//===============================================================================
// Function   : subtract
// Member Type: PUBLIC
// Description: subtract boolean operation on OCC-based bodies
// Author     : Jane Hu
// Date       : 03/08
//===============================================================================
CubitStatus     OCCModifyEngine::subtract(DLIList<BodySM*> &tool_body_list,
                                          DLIList<BodySM*> &from_bodies,
                                          DLIList<BodySM*> &new_bodies,
                                          bool imprint,
                                          bool keep_old) const
{
  //need to implement "imprint" function.
  // copy the bodies in case subtraction has some errors
  DLIList<TopoDS_Shape*> tool_bodies_copy;
  DLIList<TopoDS_Shape*> from_bodies_copy;
  DLIList<CubitBoolean> is_volume;
  for (int i = 0; i < from_bodies.size(); i++)
  {
    BodySM* body = from_bodies.get_and_step();
    OCCBody* occ_body = CAST_TO(body, OCCBody);
    OCCSurface* surface = occ_body->my_sheet_surface();
    OCCShell*   shell = occ_body->shell();
    is_volume.append( CUBIT_TRUE);
    if(surface)
    {
       TopoDS_Face* topo_face = surface->get_TopoDS_Face();
       BRepBuilderAPI_Copy api_copy(*topo_face);
       TopoDS_Shape newShape = api_copy.ModifiedShape(*topo_face);
       TopoDS_Shape* newShape_ptr = new TopoDS_Shape(newShape);
       from_bodies_copy.append(newShape_ptr);
       is_volume.change_to( CUBIT_FALSE);
    }
    else if(shell)
    {
       TopoDS_Shell* topo_shell = shell->get_TopoDS_Shell();
       BRepBuilderAPI_Copy api_copy(*topo_shell);
       TopoDS_Shape newShape = api_copy.ModifiedShape(*topo_shell); 
       TopoDS_Shape* newShape_ptr = new TopoDS_Shape(newShape);
       from_bodies_copy.append(newShape_ptr);
       is_volume.change_to( CUBIT_FALSE);
    }
    else
    {
       DLIList<Lump*> lumps = occ_body->lumps();
       if (lumps.size() > 1)
       {
	 PRINT_ERROR("Can't do boolean operation on CompSolid types. \n");
         return CUBIT_FAILURE;
       }
 
       TopoDS_Solid* solid = CAST_TO(lumps.get(), OCCLump)->get_TopoDS_Solid();
       BRepBuilderAPI_Copy api_copy(*solid);
       TopoDS_Shape newShape = api_copy.ModifiedShape(*solid);
       TopoDS_Shape* newShape_ptr = new TopoDS_Shape(newShape);
       from_bodies_copy.append(newShape_ptr);
    }
  }

  DLIList<CubitBox*> tool_boxes;
  for (int i = 0; i < tool_body_list.size(); i++)
  {
    BodySM* body = tool_body_list.get_and_step();
    OCCBody* occ_body = CAST_TO(body, OCCBody);     
    OCCSurface* surface = occ_body->my_sheet_surface();
    OCCShell*   shell = occ_body->shell();
    if(surface)
    {
       TopoDS_Face* topo_face = surface->get_TopoDS_Face();
       BRepBuilderAPI_Copy api_copy(*topo_face);
       TopoDS_Shape newShape = api_copy.ModifiedShape(*topo_face);
       TopoDS_Shape* newShape_ptr = new TopoDS_Shape(newShape);
       tool_bodies_copy.append(newShape_ptr);
    }
    else if(shell)
    {
       TopoDS_Shell* topo_shell = shell->get_TopoDS_Shell();
       BRepBuilderAPI_Copy api_copy(*topo_shell);
       TopoDS_Shape newShape = api_copy.ModifiedShape(*topo_shell);
       TopoDS_Shape* newShape_ptr = new TopoDS_Shape(newShape);
       tool_bodies_copy.append(newShape_ptr);
    }
    else
    {
       DLIList<Lump*> lumps = occ_body->lumps();
       if (lumps.size() > 1)
       {
         PRINT_ERROR("Can't do boolean operation on CompSolid types. \n");
         return CUBIT_FAILURE;
       }

       TopoDS_Solid* solid = CAST_TO(lumps.get(), OCCLump)->get_TopoDS_Solid();
       BRepBuilderAPI_Copy api_copy(*solid);
       TopoDS_Shape newShape = api_copy.ModifiedShape(*solid);
       TopoDS_Shape* newShape_ptr = new TopoDS_Shape(newShape);
       tool_bodies_copy.append(newShape_ptr);
    }
    CubitBox *tool_box = new CubitBox(occ_body->get_bounding_box());
    tool_boxes.append(tool_box);
  }

  double tol = OCCQueryEngine::instance()->get_sme_resabs_tolerance(); 
  int fraction_remaining = 100;

  // subtract the tool body from each body in the list

  CubitMessage* cmi = CubitMessage::instance(); 
  TopoDS_Shape*  from_shape = from_bodies_copy.get();
  DLIList<TopologyBridge*> tbs;
  for (int i = 0; i < from_bodies_copy.size(); i++)
  {
    BodySM* from_body = from_bodies.get();
    CubitBox box1 = CAST_TO(from_body, OCCBody)->get_bounding_box();
    int count = 0;  //count for not preforming cut
    for(int j = 0; j < tool_body_list.size(); j ++)
    {
      if (cmi->Interrupt())
      {
         PRINT_ERROR("Subtraction interrupted.  Aborting...\n");
         while (tool_boxes.size())
           delete tool_boxes.pop();
         while (tool_bodies_copy.size())
            delete tool_bodies_copy.pop();
         while (from_bodies_copy.size())
            delete from_bodies_copy.pop();
         return CUBIT_FAILURE;
      }
      CubitBox tool_box = *tool_boxes.get_and_step();  
      if(!tool_box.overlap(tol,box1))
      {
        count++;
        continue;
      } 
      TopoDS_Shape* tool_shape = tool_bodies_copy.get_and_step();
      //bodies overlap, proceed with the subtract
      TopoDS_Shape cut_shape = BRepAlgoAPI_Cut(*from_shape, *tool_shape);
 
      //compare to see if the from_shape has gotten cut.
      CubitStatus stat;
      if(is_volume[i])
      {
        GProp_GProps myProps;
        BRepGProp::VolumeProperties(*from_shape, myProps);
        double orig_mass = myProps.Mass();
        BRepGProp::VolumeProperties(cut_shape, myProps);
        double after_mass = myProps.Mass();
        if((-after_mass + orig_mass) <= tol*tol*tol)
        {
          //Add imprint code here 
          if(imprint)
            stat = imprint_toposhapes(from_shape, tool_shape);
          if(!stat)
          {
            PRINT_ERROR("Can't do imprint operation on the body. \n");
            count++;
          }
          continue;
        }
      }
      else
      {
        GProp_GProps myProps;
        BRepGProp::SurfaceProperties(*from_shape, myProps);
        double orig_mass = myProps.Mass();
        BRepGProp::SurfaceProperties(cut_shape, myProps);
        double after_mass = myProps.Mass();
        if((-after_mass + orig_mass) <= tol*tol)
        {
          //Add imprint code here
          if(imprint)
            stat = imprint_toposhapes(from_shape, tool_shape);
          if(!stat)
          {
            PRINT_ERROR("Can't do imprint operation on the body. \n");
            count++;
          }
          continue;
        }
      }
      delete from_shape;
      from_shape = new TopoDS_Shape(cut_shape);
    }

    //ok, we're done wih all cuts, construct new Body'
    if (count < tool_body_list.size() )
      tbs += OCCQueryEngine::instance()->populate_topology_bridge(*from_shape);
    else
    {
      PRINT_INFO("The %d body did not change because cutting tools are not interscting with it.\n", i+1);
      from_bodies.change_to(NULL);
    }
    from_bodies.step();
    from_shape = from_bodies_copy.step_and_get();

    // done with this j iteration; write out count, if necessary
    if (from_bodies.size() * tool_body_list.size() > 1)
    {
       int frac_done = (100 * (i+1)) / (from_bodies.size());
       if ((100 - frac_done) < fraction_remaining)
       {
          fraction_remaining = 100 - frac_done;
          PRINT_INFO("%d\% remaining.\n ", fraction_remaining+1);
       }
    }
  }

  for (int i = 0; i< tbs.size(); i++)
  {
    BodySM* bodysm = CAST_TO(tbs.get_and_step(), BodySM);
    if (bodysm)
      new_bodies.append(bodysm);
  }    

  //ok, we're done wih all cuts, delete unnecessaries. 
  while (tool_boxes.size())
    delete tool_boxes.pop();
  while (tool_bodies_copy.size())
    delete tool_bodies_copy.pop();
  if (!keep_old)
  {
    from_bodies.remove_all_with_value(NULL);
    OCCQueryEngine::instance()->delete_solid_model_entities(from_bodies); 
    OCCQueryEngine::instance()->delete_solid_model_entities( tool_body_list);
  }
  return CUBIT_SUCCESS; 
}

//===============================================================================
// Function   : imprint_toposhapes
// Member Type: PROTECTED
// Description: imprint boolean operation on OCC-based bodies
// Author     : Jane HU
// Date       : 03/08
//===============================================================================
CubitStatus OCCModifyEngine::imprint_toposhapes(TopoDS_Shape*& from_shape, 
                                                TopoDS_Shape* tool_shape)const
{
  TopOpeBRep_ShapeIntersector intersector;
  intersector.InitIntersection(*from_shape, *tool_shape);
  BRepFeat_SplitShape splitor(*from_shape);
  TopTools_ListOfShape list_of_edges;

  int max_edge = 0;
  TopoDS_Shape from_face; 
  for(; intersector.MoreIntersection(); intersector.NextIntersection())
  {
    TopoDS_Shape face1 = intersector.ChangeFacesIntersector().Face(1);
    TopoDS_Shape face2 = intersector.ChangeFacesIntersector().Face(2);
    BRepAlgoAPI_Section section(face1, face2);
    TopTools_ListOfShape temp_list_of_edges;
    temp_list_of_edges.Assign(section.SectionEdges());
    int num_edges = temp_list_of_edges.Extent();
    if (max_edge < num_edges)
    {
      list_of_edges.Assign(temp_list_of_edges);  
      max_edge =  num_edges ;
      from_face = face1;
    }
  }

  TopTools_ListIteratorOfListOfShape Itor;
  Itor.Initialize(list_of_edges);
  DLIList<Curve*> curve_list;
  for(; Itor.More(); Itor.Next())
  {
    TopoDS_Edge edge = TopoDS::Edge(Itor.Value());
    if (max_edge == 1) //surface solid imprint
    {
      splitor.Add(edge, TopoDS::Face(from_face));  
      break;
    }
    Curve* curve = OCCQueryEngine::instance()->populate_topology_bridge(edge);
    curve_list.append(curve);
  }
  if (max_edge > 1)
  {      
    DLIList<DLIList<TopoDS_Edge*>*> edge_lists;
    CubitStatus stat = sort_curves(curve_list, edge_lists); 
    if (!stat)
    {
      PRINT_ERROR("can't do solid solid imprint without a closed loop.\n");
      return CUBIT_FAILURE;
    }  
    assert(edge_lists.size() == 1);
    DLIList<TopoDS_Edge*>* edge_list;
    edge_list = edge_lists.get();
    BRepBuilderAPI_MakeWire myWire;
    for(int i = 0; i < edge_list->size(); i++)
    {
      TopoDS_Edge e = *(edge_list->get_and_step());
      myWire.Add(e); 
    }
    splitor.Add(myWire.Wire(),TopoDS::Face(from_face));
  } 
  splitor.Build();
  if(splitor.IsDone())
    delete from_shape;
  from_shape = new TopoDS_Shape(splitor.Shape());
  /*
  TopExp_Explorer Ex;
  int num_face = 0;
  for (Ex.Init(*from_shape, TopAbs_FACE); Ex.More(); Ex.Next())
  {
     TopoDS_Face face = TopoDS::Face(Ex.Current());
     num_face++;
  }
  */
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
