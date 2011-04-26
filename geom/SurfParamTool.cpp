//- Class: SurfParamTool
//-------------------------------------------------------------------------
// Filename      : SurfParamTool.cpp
//
// Purpose       : This is the generic version of ParamTool when the
//                 geometry engine has a sufficient parameterization.
//                 It uses the tool's existing functions to get transform
//                 between uv and xyz spaces.
//
// Creator       : Christopher Hynes
//
// Creation Date : 7/10/2002
//
// Owner         : Christopher Hynes
//-------------------------------------------------------------------------

#include "SurfParamTool.hpp"
#include "CastTo.hpp"
#include "Surface.hpp"
#include "DLIList.hpp" 

//-------------------------------------------------------------------------
// Function:    SurfParamTool
// Description: constructor
// Author:      chynes
// Date:        7/10/2002
//-------------------------------------------------------------------------
SurfParamTool::SurfParamTool(Surface *surf) 
{

	//- update private variables
	refSurf = surf;

}

//-------------------------------------------------------------------------
// Function:    SurfParamTool
// Description: deconstructor
// Author:      chynes
// Date:        7/10/2002
//-------------------------------------------------------------------------
SurfParamTool::~SurfParamTool() {}

//===================================================================================
// Function: set_up_space (Public)
// Description: sets up space of flattening
// Author: chynes
// Date: 7/10/02
//===================================================================================
CubitStatus SurfParamTool::set_up_space(void) {
	
	CubitStatus rv = CUBIT_SUCCESS;

	
	return rv; 
}

//===================================================================================
// Function: transform_to_uv (Public)
// Description: same as title, the local sizing will be returned in the z coord 
// Author: chynes
// Date: 7/10/02
//===================================================================================
CubitStatus SurfParamTool::transform_to_uv(const CubitVector &xyz_location, CubitVector &uv_location) 
{
	double u,v;

	CubitStatus rv = refSurf->u_v_from_position(xyz_location, u, v);
	uv_location.set(u,v,1.0);

	CubitVector du, dv;
	uv_derivitives(u,v,du,dv);

	return rv;
}

//===================================================================================
// Function: transform_to_xyz (Public)
// Description: same as title
// Author: chynes
// Date: 7/10/02
//===================================================================================
CubitStatus SurfParamTool::transform_to_xyz(CubitVector &xyz_location, const CubitVector &uv_location) 
{
	xyz_location = refSurf->position_from_u_v(uv_location.x(), uv_location.y());

	return CUBIT_SUCCESS;
}

CubitStatus SurfParamTool::uv_derivitives( double u_param, double v_param, CubitVector &du, CubitVector &dv )
{
	return refSurf->uv_derivitives (u_param, v_param, du, dv);
}

//void intersect2DLines(double x0, double y0, double x1, double y1, 
//					  double x2, double y2, double x3, double y3,
//					  double &xc, double &yc)
////
//// intersect two lines in 2D.  First line defined by points 0 and 1, second by points 2 and 3
////
//{
//	double u = ( (x3-x2)*(y0-y2) - (y3-y2)*(x0-x2) ) / ( (y3-y2)*(x1-x0) - (x3-x2)*(y1-y0) );
//
//	xc = x0 + u * (x1-x0);
//	yc = y0 + u * (y1-y0);
//}
//
//	
//
//CubitStatus SurfParamTool::circumcenter(double u0, double v0, 
//						 double u1, double v1,
//						 double u2, double v2,
//						 CubitVector &center)
////
//// Calculates the center of the center of the circumellipse for the three points in parameter space
//// this (u,v) should map back to the circumcenter of the circle of the three points in three-space.
////
//{
//
//	// first, lets calculate the circumcenter without the gradients to make sure the system works
//
//	// our method for finding a circumcenter utilizes the fact that the perpendicalar bisectors of each side 
//	// of the triangle intersect at the circumcenter
//
//	// find the three bisectors of the sides
//
//	double u3 = u0 + (u1 - u0) * 0.5;	// midpoint of side 0-1
//	double v3 = v0 + (v1 - v0) * 0.5;
//
//	double u4 = u3 + v0 - v1;	// point on the perp bisector of side 0-1
//	double v4 = v3 + u1 - u0;
//
//	double u5 = u1 + (u2 - u1) * 0.5;	// midpoint of side 1-2
//	double v5 = v1 + (v2 - v1) * 0.5;
//
//	double u6 = u5 + v1 - v2;	// point on the perp bisector of side 1-2
//	double v6 = v5 + u2 - u1;
//
//	double u7 = u0 + (u2 - u0) * 0.5;	// midpoint of side 0-2
//	double v7 = v0 + (v2 - v0) * 0.5;
//
//	double u8 = u7 + v0 - v2;
//	double v8 = v7 + u2 - u0;
//
//	// intersect two bisectors to give a cirumcenter
//
//	double xc, yc;
//	intersect2DLines(u3,v3,u4,v4,u5,v5,u6,v6,xc,yc);
//	double xc2, yc2;
//	intersect2DLines(u3,v3,u4,v4,u7,v7,u8,v8,xc2,yc2);
//	double xc3, yc3;
//	intersect2DLines(u5,v5,u6,v6,u7,v7,u8,v8,xc3,yc3);
//
//
//
//
//	//// get derivitives of each point
//
//	//CubitVector du0, du1, du2, dv0, dv1, dv2, du, dv;
//	//uv_derivitives (u0, v0, du0, dv0);
//	//uv_derivitives (u1, v1, du1, dv1);
//	//uv_derivitives (u2, v2, du2, dv2);
//
//	//// for now, average derivities first
//
//	//du = (du0 + du1 + du2) / 3.0; n
//	//dv = (dv0 + dv1 + dv2) / 3.0;
//
//	//// Calculate the rotation angle to the principal axis
//
//	//double alpha = atan2( 2.0 * ( du.x()*dv.x() + du.y()*dv.y() ), (du.length_squared() - dv.length_squared())) * 0.5;
//
//	//// Calculate principal derivitives
//
//	//double cosa = acos(alpha);
//	//double sina = asin(alpha);
//	//CubitVector dup = cosa * du + sina * dv;
//	//CubitVector dvp = -sina * du + cosa * dv;
//
//	////  
//
//
//	return CUBIT_SUCCESS;
//}




//EOF

//- Class: TestParamTool
//-------------------------------------------------------------------------
// Filename      : TestParamTool.cpp
//
// Purpose       : This is the generic version of ParamTool when the
//                 geometry engine has a sufficient parameterization.
//                 It uses the tool's existing functions to get transform
//                 between uv and xyz spaces.
//
// Creator       : Christopher Hynes
//
// Creation Date : 7/10/2002
//
// Owner         : Christopher Hynes
//-------------------------------------------------------------------------

#include "CastTo.hpp"
#include "Surface.hpp"
#include "DLIList.hpp" 

//-------------------------------------------------------------------------
// Function:    TestParamTool
// Description: constructor
// Author:      chynes
// Date:        7/10/2002
//-------------------------------------------------------------------------
TestParamTool::TestParamTool() 
{

	//- update private variables
//	refSurf = surf;

  uRange = 1.0;
  vRange = 4.0;
  
  xMin = -0.5;
  yMin = -0.5;

  xMax = 0.5;
  yMax = 0.5;


  zDepth = 0.5;  // constant z coordinate of test surface

}

//-------------------------------------------------------------------------
// Function:    TestParamTool
// Description: deconstructor
// Author:      chynes
// Date:        7/10/2002
//-------------------------------------------------------------------------
TestParamTool::~TestParamTool() {}

//===================================================================================
// Function: set_up_space (Public)
// Description: sets up space of flattening
// Author: chynes
// Date: 7/10/02
//===================================================================================
CubitStatus TestParamTool::set_up_space(void)
{
	
	CubitStatus rv = CUBIT_SUCCESS;

	
	return rv; 
}

//===================================================================================
// Function: transform_to_uv (Public)
// Description: same as title, the local sizing will be returned in the z coord 
// Author: chynes
// Date: 7/10/02
//===================================================================================
CubitStatus TestParamTool::transform_to_uv(const CubitVector &xyz_location, CubitVector &uv_location) 
{
  // calculate the u,v from the x,y

  double u = (xyz_location.x()-xMin)/(xMax-xMin) * uRange;
  double v = (xyz_location.y()-yMin)/(yMax-yMin) * vRange;
	uv_location.set(u,v,1.0);

	return CUBIT_SUCCESS;
}

//===================================================================================
// Function: transform_to_xyz (Public)
// Description: same as title
// Author: chynes
// Date: 7/10/02
//===================================================================================
CubitStatus TestParamTool::transform_to_xyz(CubitVector &xyz_location, const CubitVector &uv_location) 
{
  double u = (xyz_location.x()-xMin)/(xMax-xMin) * uRange;
  // get x,y from u,v

  double x = uv_location.x()/uRange + xMin;
  double y = uv_location.y()/vRange + yMin;

  xyz_location.set(x,y,zDepth);

	return CUBIT_SUCCESS;
}

CubitStatus TestParamTool::uv_derivitives( double u_param, double v_param, CubitVector &du, CubitVector &dv )
{
  du.set(1.0/uRange,0,0);
  dv.set(0,1.0/vRange,0);
  return CUBIT_SUCCESS;
}
