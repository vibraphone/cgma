//-------------------------------------------------------------------------
// Filename      : AutoMidsurfaceTool.cpp
//
// Purpose       : 
//   Create a midsurface of a body/volume given a thickness range. If no thickness range is given
//   then make a educated guess at the thickness (using something like Volume/Surface Area).
//   
//   Create Midsurface Volume <id_list> auto [<lower_tol> <upper_tol>]
//
// Creator       : Sam Showman
//
// Creation Date : 05/10/2008
//-------------------------------------------------------------------------
#ifndef AutoMidsurfaceTool_HPP
#define AutoMidsurfaceTool_HPP

#include "CubitDefines.h"
#include "CubitVector.hpp"
#include "CubitPlane.hpp"

class GeometryModifyEngine;
class Body;
class RefFace;
class BodySM;
class Surface;
class Curve;
class RefEdge;
template <class X> class DLIList;

class AutoMidsurfaceTool
{
public:

    AutoMidsurfaceTool();
    ~AutoMidsurfaceTool(){}

    CubitStatus midsurface(
	DLIList<Body*> &body_list_in,
    DLIList<BodySM*> &body_list_out,
    DLIList<Body*> &old_bodies_midsurfaced,
    DLIList<double> &thickness_out,
    double lower_tol = CUBIT_DBL_MAX,
    double upper_tol = CUBIT_DBL_MAX,
    CubitBoolean delete_midsurfaced = CUBIT_FALSE,
    CubitBoolean preview = CUBIT_FALSE);
    //- automatically midsurfaces a volume based on surface pairs and surface area
    //- body_list_in - list of bodies to midsurface
	//- body_list_out - result bodies
    //- lower_tol - lower tolerance
    //- upper_tol - upper tolerance
    //- delete_midsurfaced - delete the midsurfaced solid
    //- transp_midsurfaced - make the midsurfaced solid transparent
	//- preview - preview the results

private:
	CubitBoolean curve_in_surface(Curve *curve_in, Surface *surf_in);

	CubitStatus get_planar_mid_surface(
		RefFace* ref_face1,
		RefFace* ref_face2,
		BodySM* body_sm_to_trim_to,
		BodySM*& midsurface_body_sm,
		GeometryModifyEngine *gme_ptr );

	CubitStatus find_offset_pair_patches(
        DLIList<RefFace*> pairs_list_0,
        DLIList<RefFace*> pairs_list_1,	
        DLIList<RefFace*>& red_faces,
        DLIList<RefFace*>& yellow_faces,
        DLIList<double>& offset_distances);

    CubitStatus random_loc_on_surface( Surface* face_ptr, CubitVector &loc );

    CubitBoolean check_surf_pairs(double min_thick, double max_thick,
        DLIList<RefFace*> check_list, Body* body_in );
};

#endif

