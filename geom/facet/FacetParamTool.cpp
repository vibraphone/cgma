//===================================================================================
//
// File: FacetParamTool.cpp
//
// Descripiton: interface for 3D-2D parameterization algorithms
//
//===================================================================================


#include "FacetParamTool.hpp"
#include "CastTo.hpp"
#include "Surface.hpp"
#include "FacetSurface.hpp"
#include "FacetEvalTool.hpp"
#include "FacetorTool.hpp"
#include "CubitPoint.hpp"
#include "CubitPointData.hpp"
#include "CubitFacet.hpp"
#include "CubitFacetData.hpp"
#include "CubitBox.hpp"
#include "TDFacetBoundaryPoint.hpp"

#ifdef ROADKILL
#include "ParamManager.hpp"
#endif

#define DETERM(p1,q1,p2,q2,p3,q3)\
     ((q3)*((p2)-(p1)) + (q2)*((p1)-(p3)) + (q1)*((p3)-(p2)))
#define INSIDE_TOL 1.0e-6

//===================================================================================
// Function: FacetParamTool (Public)
// Description: Constructor
// Author: chynes
// Date: 1/02
//===================================================================================
FacetParamTool::FacetParamTool(int numn, int nume, double* nodes, int* tri) 
{}


//-------------------------------------------------------------------------
// Function:    FacetParamTool
// Description: constructor
// Author:      chynes
// Date:        7/10/2002
//-------------------------------------------------------------------------
FacetParamTool::FacetParamTool(Surface *surf) 
{
	refSurf = surf;
}

//===================================================================================
// Function: ~FacetParamTool
// Description: Deconstructor
// Author: chynes
// Date: 1/02
//===================================================================================
FacetParamTool::~FacetParamTool() {}

//===================================================================================
// Function: set_up_space (Public)
// Description: sets up space of flattening, flattens, and then updates surface
// Author: chynes
// Date: 7/10/02
//===================================================================================
CubitStatus FacetParamTool::set_up_space() 
{ 
#ifdef ROADKILL
	FacetEvalTool *fetool = CAST_TO(refSurf, FacetSurface)->get_eval_tool();

	if ((fetool->loops())->size() != 1)
	{
		PRINT_WARNING("Cannot parameterize surface.  Multiple loops detected\n");
		return CUBIT_FAILURE;
	}

	DLIList<CubitPoint *> my_point_list;
	fetool->get_points(my_point_list);
	DLIList<CubitFacet *> my_facet_list;
	fetool->get_facets(my_facet_list);

  // make arrays out of the points and facets
	double *points = new double [3 * my_point_list.size()];
	int *facets = new int [3 * my_facet_list.size()];
	if (!points || !facets)
	{
		PRINT_ERROR("Could not define parameterization for surface (out of memory)\n");
		return CUBIT_FAILURE;
	}
	int ii;
	CubitPoint *pt;
	my_point_list.reset();
	for (ii=0; ii<my_point_list.size(); ii++)
	{
		pt = my_point_list.get_and_step();
		points[ii*3] = pt->x();
		points[ii*3+1] = pt->y();
		points[ii*3+2] = pt->z();
		pt->marked(ii);
	}
	CubitFacet *facet;
	CubitPoint *pts[3];
	for (ii=0; ii<my_facet_list.size(); ii++)
	{
		facet = my_facet_list.get_and_step();    
		facet->points( pts[0], pts[1], pts[2] );
		facets[ii*3]   = pts[0]->marked();
		facets[ii*3+1] = pts[1]->marked();
		facets[ii*3+2] = pts[2]->marked();
	}

	//debug, export facets
//	export_facets(my_point_list.size(), my_facet_list.size(), points, facets);
//	PRINT_ERROR("Debugging Failure.  Please remove export_facets from FacetParamTool\n");
//	return CUBIT_FAILURE;

	// do the parameterization
	ParamManager parammanager( my_point_list.size(), my_facet_list.size(),
							points, facets );
	if(!parammanager.flatten())
	{
		PRINT_ERROR("Surface Parameterizer Failed\n");
		return CUBIT_FAILURE;
	}
	else
	{
		double ratio;
		double *sizes;
		parammanager.build_sizing();
		double *uv = parammanager.get_uvs_sizing( ratio, sizes ); 
		
		// update the points with uvs values

		TDFacetBoundaryPoint *td_fbp;
		CubitBoolean on_internal_boundary;
		my_point_list.reset();
		for (ii=0; ii<my_point_list.size(); ii++)
		{
			pt = my_point_list.get_and_step();
			DLIList <CubitFacet *> facet_list;
			pt->facets_on_surf( fetool->tool_id(), facet_list, on_internal_boundary );
			if (on_internal_boundary)
			{
				td_fbp = TDFacetBoundaryPoint::get_facet_boundary_point( pt );
				if (!td_fbp)
				{
					TDFacetBoundaryPoint::add_facet_boundary_point( pt );
					td_fbp = TDFacetBoundaryPoint::get_facet_boundary_point( pt );
					td_fbp->add_surf_facets( facet_list );
					td_fbp->set_uvs( facet_list.get(), uv[ii*2], uv[ii*2+1], sizes[ii]);
				}
				else
				{
					if (td_fbp->set_uvs( facet_list.get(),
						uv[ii*2], uv[ii*2+1], sizes[ii] ) != CUBIT_SUCCESS)
					{
						td_fbp->add_surf_facets( facet_list );
						td_fbp->set_uvs( facet_list.get(), uv[ii*2], uv[ii*2+1], sizes[ii]);
					}
				}
			}
			else
			{
				pt->set_uvs( uv[ii*2], uv[ii*2+1], sizes[ii]);
			}
		}
	
		PRINT_INFO("Surface Parameterization succeeded\n");
		delete [] sizes;
		delete [] uv;
    
	}
	
	// clean up

	delete [] points;
	delete [] facets;
	return CUBIT_SUCCESS;

#else
	return CUBIT_FAILURE;
#endif

}

//===================================================================================
// Function: transform_to_uv (Public)
// Description: same as title, set_up_space must have been already called
// the local sizing will be returned in the z coord 
// of the uv location vector
// Author: chynes
// Date: 7/10/02
//===================================================================================
CubitStatus FacetParamTool::transform_to_uv(const CubitVector &xyz_location, CubitVector &uv_location) 
{
	DLIList<CubitFacet *> facet_list; 
	FacetEvalTool *fetool;
	CubitFacet *tri_ptr;
	CubitBoolean outside = CUBIT_FALSE;
	CubitVector closest_point;
	CubitVector area_coord;
	double u, v, s;
	double compare_tol;

	// find best compare_tol
	fetool = CAST_TO(refSurf, FacetSurface)->get_eval_tool();
	compare_tol = 1e-3*(fetool->bounding_box().diagonal().length());

	// locate point
	CubitStatus rv = CUBIT_SUCCESS;
	fetool->get_facets(facet_list);
	rv = FacetEvalTool::project_to_facets(facet_list, 
										  tri_ptr, 
										  0, 
										  compare_tol, 
										  xyz_location, 
										  CUBIT_FALSE, 
										  &outside, 
										  &closest_point, 
										  NULL);
	// determine barycentric coordinates for in facet
	if(rv == CUBIT_SUCCESS) 
	{
		FacetEvalTool::facet_area_coordinate(tri_ptr, closest_point, area_coord);

		// extract data
		double u0, u1, u2, v0, v1, v2, s0, s1, s2;
		CubitPoint *p0, *p1, *p2;
		tri_ptr->points(p0, p1, p2);
		p0->get_uvs(tri_ptr, u0, v0, s0);
		p1->get_uvs(tri_ptr, u1, v1, s1);
		p2->get_uvs(tri_ptr, u2, v2, s2);

		// determine u coordinate
		u = (area_coord.x()*u0) + (area_coord.y()*u1) + (area_coord.z()*u2);

		// determine v coordinate
		v = (area_coord.x()*v0) + (area_coord.y()*v1) + (area_coord.z()*v2);

		// determine sizing
		s = (area_coord.x()*s0) + (area_coord.y()*s1) + (area_coord.z()*s2);

		uv_location.set(u,v,s);
	}

		return rv;	
}

//===================================================================================
// Function: transform_to_xyz (Public)
// Description: same as title
// Author: chynes
// Date: 7/10/02
//===================================================================================
CubitStatus FacetParamTool::transform_to_xyz(CubitVector &xyz_location, const CubitVector &uv_location) 
{
	CubitStatus rv = CUBIT_SUCCESS;
	CubitFacet *tri_ptr;
	FacetSurface *facet_surf = CAST_TO(refSurf, FacetSurface);
	CubitVector area_coord;
	double x, y, z;

	//locate point
	rv = locate_point_in_uv(facet_surf, uv_location, tri_ptr);

	if(rv == CUBIT_SUCCESS)
	{
		// create uv facet
		CubitPoint* uvpoint_array[3];
		DLIList<CubitPoint *> point_list;
		CubitPoint *pt_ptr;
		double u,v;
		int ii;
		tri_ptr->points(point_list);
		for(ii = 0; ii<3; ii++)
		{
			pt_ptr = point_list.get_and_step();
			pt_ptr->get_uv(tri_ptr, u, v);
			uvpoint_array[ii] = (CubitPoint *) new CubitPointData(u, v, 0.0);
		}
		CubitFacet *uvfacet_ptr = (CubitFacet *) new CubitFacetData( uvpoint_array[0], 
																	 uvpoint_array[1], 
																	 uvpoint_array[2] );

		// determine barycentric coordinates of uv point in uv facet
		FacetEvalTool::facet_area_coordinate(uvfacet_ptr, uv_location, area_coord);

		//delete created objects
		delete uvfacet_ptr;
		for(ii = 0; ii<3; ii++)
		{
			pt_ptr = uvpoint_array[ii];
			delete pt_ptr;
		}

		CubitPoint *p0, *p1, *p2;
		CubitVector coord0, coord1, coord2;
		tri_ptr->points(p0,p1,p2);
		coord0 = p0->coordinates();
		coord1 = p1->coordinates();
		coord2 = p2->coordinates();

		//determine x coordinate
		x =   (area_coord.x()*coord0.x())
			+ (area_coord.y()*coord1.x())
			+ (area_coord.z()*coord2.x());
		//determine y coordinate
		y =   (area_coord.x()*coord0.y())
			+ (area_coord.y()*coord1.y())
			+ (area_coord.z()*coord2.y());
		//determine z coordinate
		z =   (area_coord.x()*coord0.z())
			+ (area_coord.y()*coord1.z())
			+ (area_coord.z()*coord2.z());

		xyz_location.set(x,y,z);

	}

	return rv;
}


//===================================================================================
// Function: locate_point_in_uv (Public)
// Description: same as title (the_point must have coords (u,v,0.0))
// Author: chynes
// Date: 7/10/02
//===================================================================================
CubitStatus FacetParamTool::
locate_point_in_uv(FacetSurface *surf, const CubitVector &the_point, CubitFacet *&tri_ptr) 
{
	CubitStatus rv = CUBIT_SUCCESS;

	DLIList<CubitFacet *> facet_list;
	int tool_id = surf->get_eval_tool()->tool_id();
	surf->get_eval_tool()->get_facets(facet_list);
	tri_ptr = facet_list.get();

	// loop until we find something

	CubitPoint *p0, *p1, *p2;
	double aa, bb, cc;
	double u0, u1, u2, v0, v1, v2;
	CubitBoolean found = CUBIT_FALSE;
	while(!found && rv == CUBIT_SUCCESS)
	{
		tri_ptr->points(p0, p1, p2);
		p0->get_uv(tri_ptr, u0, v0);
		p1->get_uv(tri_ptr, u1, v1);
		p2->get_uv(tri_ptr, u2, v2);
		aa = DETERM(the_point.x(), the_point.y(), u1, v1, u2, v2);
		bb = DETERM(u0, v0, the_point.x(), the_point.y(), u2, v2);
		cc = DETERM(u0, v0, u1, v1, the_point.x(), the_point.y());
		if (aa > -INSIDE_TOL &&
			bb > -INSIDE_TOL &&
			cc > -INSIDE_TOL)
		{
			found = CUBIT_TRUE;  // this is the one
		}
		else
		{
			// set up to check the next logical neighbor
			if (aa <= bb && aa <= cc) {
				int edge_index = 1;
				tri_ptr = tri_ptr->adjacent( edge_index, &tool_id );
			}
			else if (bb <= aa && bb <= cc) {
				int edge_index = 2;
				tri_ptr = tri_ptr->adjacent( edge_index, &tool_id );
			}
			else {
				int edge_index = 0;
				tri_ptr = tri_ptr->adjacent( edge_index, &tool_id);
			}
			// check to see if we've left the triangulation
      
			if (tri_ptr == NULL)
			{
				rv = exhaustive_locate_point_in_uv(surf, the_point, tri_ptr );
				found = CUBIT_TRUE;
			}
		}
	}

	return rv;
}
//===================================================================================
// Function: exhaustive_locate_point_in_uv (Public)
// Description: same as title
// Author: chynes
// Date: 7/10/02
//===================================================================================
CubitStatus FacetParamTool::
exhaustive_locate_point_in_uv(FacetSurface *surf, const CubitVector &the_point, CubitFacet *&tri_ptr) 
{
	CubitStatus rv = CUBIT_SUCCESS;

	DLIList<CubitFacet *> facet_list;
	surf->get_eval_tool()->get_facets(facet_list);

	// loop until we find something

	CubitPoint *p0, *p1, *p2;
	int ii;
	double aa, bb, cc;
	double u0, u1, u2, v0, v1, v2;
	CubitBoolean found = CUBIT_FALSE;
	for(ii = 0; ii < facet_list.size() && !found; ii++)
	{
		tri_ptr = facet_list.get_and_step();
		tri_ptr->points(p0, p1, p2);
		p0->get_uv(tri_ptr, u0, v0);
		p1->get_uv(tri_ptr, u1, v1);
		p2->get_uv(tri_ptr, u2, v2);
		aa = DETERM(the_point.x(), the_point.y(), u1, v1, u2, v2);
		bb = DETERM(u0, v0, the_point.x(), the_point.y(), u2, v2);
		cc = DETERM(u0, v0, u1, v1, the_point.x(), the_point.y());
		if (aa > -INSIDE_TOL &&
			bb > -INSIDE_TOL &&
			cc > -INSIDE_TOL)
		{
			found = CUBIT_TRUE;  // this is the one
		}
	}
	if (!found)
	{
		rv = CUBIT_FAILURE;
		tri_ptr = NULL;
	}
	return rv;
}

//===================================================================================
// Function: flatten
// Description: accessor function to the parameterization algorithms
// Author: chynes
// Date: 1/02
//====================================================================================


//EOF
