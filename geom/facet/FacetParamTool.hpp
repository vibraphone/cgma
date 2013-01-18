//-----------------------------------------------------------------------------
//
//   File:FacetParamTool.hpp
//   
//   Purpose: interface for 3D-2D parameterization algorithms for facets
//
//
//-----------------------------------------------------------------------------



#ifndef FACET_PARAM_TOOL_HPP
#define FACET_PARAM_TOOL_HPP

#include "ParamTool.hpp"
class Surface;
class FacetSurface;
class CubitFacet;
class CubitVector;

class FacetParamTool : public ParamTool
{

public:
  FacetParamTool(int numn, int nume, double* nodes, int* tri);
       //~constructor
  FacetParamTool(Surface *surf);
	//- constructor

  ~FacetParamTool();
       //~deconstructor

  CubitStatus set_up_space(void);

  CubitStatus transform_to_uv(const CubitVector &xyz_location, CubitVector &uv_location);

  CubitStatus transform_to_xyz(CubitVector &xyz_location, const CubitVector &uv_location);

  CubitStatus locate_point_in_uv(FacetSurface *surf, const CubitVector &the_point, CubitFacet *&tri_ptr);

  CubitStatus exhaustive_locate_point_in_uv(FacetSurface *surf, const CubitVector &the_point, CubitFacet *&tri_ptr); 

private:

	Surface *refSurf;
	//- reference surface
  
};

#endif

