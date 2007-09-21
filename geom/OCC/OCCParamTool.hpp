//-----------------------------------------------------------------------------
//
//   File:OCCParamTool.hpp
//   
//   Purpose: interface for 3D-2D parameterization algorithms for facets
//
//
//-----------------------------------------------------------------------------



#ifndef FACET_PARAM_TOOL_HPP
#define FACET_PARAM_TOOL_HPP

#include "ParamTool.hpp"
class Surface;
class OCCSurface;
class CubitFacet;
class CubitVector;

class OCCParamTool : public ParamTool
{

public:
  OCCParamTool(int numn, int nume, double* nodes, int* tri);
       //~constructor
  OCCParamTool(Surface *surf);
	//- constructor

  ~OCCParamTool();
       //~deconstructor

  CubitStatus set_up_space(void);

  CubitStatus transform_to_uv(CubitVector &xyz_location, CubitVector &uv_location);

  CubitStatus transform_to_xyz(CubitVector &xyz_location, CubitVector &uv_location);

  CubitStatus locate_point_in_uv(OCCSurface *surf, CubitVector &the_point, CubitFacet *&tri_ptr);

  CubitStatus exhaustive_locate_point_in_uv(OCCSurface *surf, CubitVector &the_point, CubitFacet *&tri_ptr); 

#ifdef BOYD14
  CubitStatus export_facets(int numn, int numf, double *points, int *facets);
  // - debug function for testing inputs to roadkill

  int flatten();
      //accessor function to the parameterization algorithms

  double* get_uvs_sizing(double& ratio, double*& sizings);
      //after calling flatten, returns an array of the u,v parameters, accessing
      //d_flat_mesh
#endif
private:

	Surface *refSurf;
	//- reference surface
  
};

#endif

