//- Class: PeriodicParamTool
//-------------------------------------------------------------------------
// Filename      : PeriodicParamTool.cpp
//
// Purpose       : Handles surfaces that have a periodic parameterization
//		           Uses, but modifies, the underlying engine parameterization
//
//
// Creator       : Ray J. Meyers
//
// Creation Date : 12/15/2008
//
// Owner         : Ray J. Meyers
//-------------------------------------------------------------------------

#include "PeriodicParamTool.hpp"
//#include "CastTo.hpp"
#include "Surface.hpp"
//#include "DLIList.hpp"

//-------------------------------------------------------------------------
// Function:    PeriodicParamTool
// Description: constructor
// Author:      chynes
// Date:        7/10/2002
//-------------------------------------------------------------------------
PeriodicParamTool::PeriodicParamTool(Surface *surf) 
{
	//- update private variables
	refSurf = surf;
	uOffset = 0.0;
	vOffset = 0.0;
  uPeriod = 0.0;
  vPeriod = 0.0;
  mirrorSurface = false;

}

//-------------------------------------------------------------------------
// Function:    PeriodicParamTool
// Description: deconstructor
// Author:      chynes
// Date:        7/10/2002
//-------------------------------------------------------------------------
PeriodicParamTool::~PeriodicParamTool() {}

//===================================================================================
// Function: set_up_space (Public)
// Description: sets up space of flattening
// Author: chynes
// Date: 7/10/02
//===================================================================================
CubitStatus PeriodicParamTool::set_up_space(double u_period, double v_period, double u_offset, double v_offset)
{
	// store the u and periods 

	uPeriod = u_period;
	vPeriod = v_period;
	uOffset = u_offset;
	vOffset = v_offset;

	
	CubitStatus rv = CUBIT_SUCCESS;

	return rv; 
}

//===================================================================================
// Function: transform_to_uv (Public)
// Description: same as title, the local sizing will be returned in the z coord 
// Author: chynes
// Date: 7/10/02
//===================================================================================
CubitStatus PeriodicParamTool::transform_to_uv(const CubitVector &xyz_location, CubitVector &uv_location) 
{
	 
	double u,v;

	CubitStatus rv = refSurf->u_v_from_position(xyz_location, u, v);

	// offset values to avoid parameter discontinuity

	if (uPeriod && u < uOffset)
	{
		u += uPeriod;
	}
  
  // mirror surface if required to get correct loop orientation
  if (mirrorSurface)
  {
    u = -u;
  }

	if (vPeriod && v < vOffset)
	{
		v += vPeriod;
	}

	uv_location.set(u,v,1.0);


	return rv;
}

//===================================================================================
// Function: transform_to_xyz (Public)
// Description: same as title
// Author: chynes
// Date: 7/10/02
//===================================================================================
CubitStatus PeriodicParamTool::transform_to_xyz(CubitVector &xyz_location, const CubitVector &uv_location) 
{
	double u = uv_location.x();	
  if (mirrorSurface)
  {
    u = -u;
  }
	if (u > uPeriod)
	{
		u = u - uPeriod;
	}
	double v = uv_location.y();
	if (v > vPeriod)
	{
		v = v - vPeriod;
	}
	xyz_location = refSurf->position_from_u_v(u,v);


	return CUBIT_SUCCESS;
}

void PeriodicParamTool::mirror_surface(bool true_false)
{
  mirrorSurface = true_false;
  if (mirrorSurface)
    PRINT_INFO("Loops appear backwards, mirroring surface...\n");
 
}

CubitStatus PeriodicParamTool::uv_derivitives(double u_param, double v_param, 
                                              CubitVector &du, CubitVector &dv)
{
  if (mirrorSurface)
    u_param = -u_param;
  if (u_param > uPeriod)
    u_param = u_param-uPeriod;
  if (v_param > vPeriod)
    v_param = v_param-vPeriod;
  refSurf->uv_derivitives(u_param, v_param, du, dv);
  return CUBIT_SUCCESS;
}
//EOF
