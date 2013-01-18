//- Class: PeriodicParamTool
//-------------------------------------------------------------------------
// Filename      : PeriodicParamTool.cpp
//
// Purpose       : Handles surfaces that have a periodic parameterization
//		           Uses, but modifies, the underlying engine parameterization
//
// Creator       : Ray J. Meyers
//
// Creation Date : 12/15/2008
//
// Owner         : Ray J. Meyers
//-------------------------------------------------------------------------

#ifndef PERIODIC_PARAM_TOOL_HPP
#define PERIODIC_PARAM_TOOL_HPP

#include "ParamTool.hpp"
#include "CubitGeomConfigure.h"

class Surface;
class CubitVector;
class CubitNode;
template <class X> class DLIList;

class CUBIT_GEOM_EXPORT PeriodicParamTool : public ParamTool
{
public:

	//- constructor
	PeriodicParamTool(Surface *surf);

	//- destructor
	~PeriodicParamTool();
	
	CubitStatus set_up_space(double u_period, double v_period, double u_offset, double v_offset);

	CubitStatus transform_to_uv(const CubitVector &xyz_location, CubitVector &uv_location);

	CubitStatus transform_to_xyz(CubitVector &xyz_location, const CubitVector &uv_location);

  void mirror_surface( bool true_false);

  CubitStatus uv_derivitives( double u_param, double v_param, 
		CubitVector &du, CubitVector &dv );

private:

	Surface *refSurf;
	double uPeriod;
	double vPeriod;
	double uOffset;
	double vOffset;
  bool mirrorSurface;

	//- reference surface
};

#endif // ACIS_PARAM_TOOL_HPP

