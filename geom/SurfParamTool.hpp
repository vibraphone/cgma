//- Class: SurfParamTool
//-------------------------------------------------------------------------
// Filename      : SurfParamTool.hpp
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

#ifndef SURF_PARAM_TOOL_HPP
#define SURF_PARAM_TOOL_HPP

#include "ParamTool.hpp"
#include "CubitGeomConfigure.h"

class Surface;
class CubitVector;
template <class X> class DLIList;

//! Transforms between a surface's uv parameter and xyz space. 
class CUBIT_GEOM_EXPORT SurfParamTool : public ParamTool
{
public:

	//- constructor
	SurfParamTool(Surface *surf);

	//- deconstructor
	~SurfParamTool();
	
	CubitStatus set_up_space(void);

	CubitStatus transform_to_uv(const CubitVector &xyz_location, CubitVector &uv_location);

	CubitStatus transform_to_xyz(CubitVector &xyz_location, const CubitVector &uv_location);

	CubitStatus uv_derivitives( double u_param, double v_param, 
                                    CubitVector &du, CubitVector &dv );

	static CubitStatus circumcenter(double u0, double v0, 
                                        double u1, double v1,
                                        double u2, double v2,
                                        CubitVector &center);  // temp function to be moved elsewhere

private:

	Surface *refSurf;
	//- reference surface
};

class Surface;
class CubitVector;
template <class X> class DLIList;

//! Transforms between a surface's uv parameter and xyz space. 
class CUBIT_GEOM_EXPORT TestParamTool : public ParamTool
{
public:

	//- constructor
	TestParamTool();

	//- deconstructor
	~TestParamTool();
	
	CubitStatus set_up_space(void);

	CubitStatus transform_to_uv(const CubitVector &xyz_location, CubitVector &uv_location);

	CubitStatus transform_to_xyz(CubitVector &xyz_location, const CubitVector &uv_location);

	CubitStatus uv_derivitives( double u_param, double v_param, 
                                    CubitVector &du, CubitVector &dv );

	static CubitStatus circumcenter(double u0, double v0, 
                                        double u1, double v1,
                                        double u2, double v2,
                                        CubitVector &center);  // temp function to be moved elsewhere

private:

//	Surface *refSurf;
	//- reference surface

  double uRange, vRange;
  double zDepth;

  double xMin, yMin, xMax, yMax;

};

#endif // ACIS_PARAM_TOOL_HPP

