#undef NDEBUG
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <list>
#include <map>
#include <string>

#include "GeometryModifyTool.hpp"
#include "GeometryQueryTool.hpp"
#include "GMem.hpp"
#include "OCCModifyEngine.hpp"
#include "OCCQueryEngine.hpp"
#include "RefVertex.hpp"
#include "RefEdge.hpp"
#include "RefFace.hpp"
#include "RefVolume.hpp"
#include "InitCGMA.hpp"
#include "CubitCompat.hpp"

#define STRINGIFY(S) XSTRINGIFY(S)
#define XSTRINGIFY(S) #S
#ifndef SRCDIR
#  define SRCDIR "."
#endif

using std::list;
using std::pair;
using std::set;
using std::map;
using std::vector;

int main( int argc, char** argv ) {

  CubitStatus status = InitCGMA::initialize_cgma("OCC");
  if (CUBIT_SUCCESS != status) return 1;
  GeometryQueryTool *gti = GeometryQueryTool::instance();
   
  // Read in the geometry from files specified on the command line
  const char* filename = STRINGIFY(SRCDIR) "/LeverArm.brep";
  status = CubitCompat_import_solid_model(filename, "OCC");
  if (status != CUBIT_SUCCESS) {
    PRINT_ERROR("Problems reading geometry file %s.\n", filename);
    return 1;
  }
 
  DLIList<RefEdge*> my_curvs;
  DLIList<RefFace*> my_surfs;
  gti->ref_edges(my_curvs);
  gti->ref_faces(my_surfs);

  // ***** Discretize the curves ***** //

  map< RefEdge*, list<double> > discrete_curves;

  for(int i = 0; i < my_curvs.size(); ++i) {
    
    list<double> sample_params;
    RefEdge* this_curv = my_curvs.get_and_step();
    //Insert nine arbitrary points on each curves
    double delta = this_curv->end_param() - this_curv->start_param();
    for(int j = 1; j < 10; ++j)
      sample_params.push_back( this_curv->start_param() + (static_cast<double>(j) * delta / 10.) );
    
    discrete_curves.insert( std::make_pair(this_curv, sample_params) );

  }

  //Pick one of the surfaces arbitrarily and attempt to find the u-v coordinates
  //of the points discretizing the boundary of this surface.

  {

    RefFace* this_surf = my_surfs.next(7);
    DLIList<RefVertex*> surf_points;
    DLIList<RefEdge*> surf_curves;
    this_surf->ref_vertices(surf_points);
    this_surf->ref_edges(surf_curves);
    
    CubitVector coord0, coord1;
    double u, v;

    double max_u, min_u, max_v, min_v;
    this_surf->get_param_range_U(min_u, max_u);
    this_surf->get_param_range_V(min_v, max_v);

    //Deals with curve end points.
    
    for(int i = 0; i < surf_points.size(); ++i) {

      RefVertex* this_point = surf_points.get_and_step();
      
      u = v = -1.e300;
      coord0 = this_point->coordinates();
      this_surf->u_v_from_position(coord0, u, v, &coord1);

      if( u < min_u || v < min_v ) {
	printf("Failed to find a parameter with valid range for an apex point\n");
	printf("  point coords:     %e %e %e\n", coord0.x(), coord0.y(), coord0.z());
	printf("  projected coords: %e %e %e\n", coord1.x(), coord1.y(), coord1.z());
	printf("  params returned:  %e %e\n\n", u, v);
        return 1;
      }

    }
    
    //Deals with discretization points on curves.

    for(int i = 0; i < surf_curves.size(); ++i) {
      
      RefEdge* this_curv = surf_curves.get_and_step(); 
     
      assert(discrete_curves.count(this_curv) == 1);
      list<double> curve_params = discrete_curves.find(this_curv)->second;
      list<double>::iterator it = curve_params.begin(), it_end = curve_params.end(); 

      for( ; it != it_end; ++it) {

	u = v = -1.e300;
	this_curv->position_from_u(*it, coord0);
	this_surf->u_v_from_position(coord0, u, v, &coord1);

	if( u < min_u || v < min_v ) {
	  printf("Failed to find a parameter with valid range for a curve point\n");
	  printf("  point coords:     %e %e %e\n", coord0.x(), coord0.y(), coord0.z());
	  printf("  projected coords: %e %e %e\n", coord1.x(), coord1.y(), coord1.z());
	  printf("  params returned:  %e %e\n\n", u, v);
          return 1;
	}

      }

    }

  }    
    
  printf("done...\n");

  return(0);

}
