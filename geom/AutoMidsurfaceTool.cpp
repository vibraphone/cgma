//-------------------------------------------------------------------------
// Filename      : AutoMidsurfaceTool.cpp
//
// Purpose       : 
//   Create a midsurface of a body/volume given a thickness range. If no thickness range is given
//   then make a educated guess at the thickness (using something like Volume/Surface Area).
//   
//   Create Midsurface Volume <id_list> auto [<min_thickness> <max_thickness>]
//
// Creator       : Sam Showman
//
// Creation Date : 05/10/2008
//-------------------------------------------------------------------------
#include "AutoMidsurfaceTool.hpp"
#include "GeometryModifyEngine.hpp"
#include "GeometryQueryTool.hpp"
#include "BodySM.hpp"
#include "Body.hpp"
#include "GeometryModifyTool.hpp"
#include "RefFace.hpp"
#include "Surface.hpp"
#include "Curve.hpp"
#include "RefEntity.hpp"
#include "SurfaceOverlapTool.hpp"
#include "CubitPlane.hpp"
#include "GfxPreview.hpp"
#include "GfxDebug.hpp"
#include "RefVertex.hpp"
#include "ProgressTool.hpp"
#include "AppUtil.hpp"
#include "BoundingBoxTool.hpp"
#include "GeomMeasureTool.hpp"

// Constructor - nothing going on here
AutoMidsurfaceTool::AutoMidsurfaceTool()
{
}

// The main midsurface function
// lower_tol, upper_tol and preview are optional
CubitStatus AutoMidsurfaceTool::midsurface(
	DLIList<Body*> &body_list_in,
	DLIList<BodySM*> &body_list_out,
    DLIList<Body*> &old_bodies_midsurfaced,
    DLIList<double> &thickness_out,
	double lower_limit,
	double upper_limit,
    CubitBoolean delete_midsurfaced,
    CubitBoolean preview)
{
    if(lower_limit == CUBIT_DBL_MAX)// no limit set
        lower_limit = -CUBIT_DBL_MAX;
    double lower_tol = CUBIT_DBL_MAX;
    double upper_tol = CUBIT_DBL_MAX;
    const double auto_thickness_margin = 0.05; // if the user wants to automatically find the search
	                                          // thickness then this var give the search margin around the
	                                          // guess
    ProgressTool* prog_tool = 0;
    if(body_list_in.size()>5)
        prog_tool = AppUtil::instance()->progress_tool();

	// At lease one body must be provided
	if(body_list_in.size() < 1)
	{
		PRINT_ERROR( "No bodies given for midsurfacing\n" );
		return CUBIT_FAILURE;
	}

	// The surfaceOverlapTool is persistent so we need to save the 
	// max_gap and such to restore them at the end or if we error out
	// save current settings
	double max_gap_save = SurfaceOverlapTool::instance()->get_gap_max();
	double min_gap_save = SurfaceOverlapTool::instance()->get_gap_min();
	int normal_type = SurfaceOverlapTool::instance()->get_normal_type();
	CubitBoolean cubit_bool_save = SurfaceOverlapTool::instance()->get_check_within_bodies();
    CubitBoolean skip_facing_surfaces = SurfaceOverlapTool::instance()->get_skip_facing_surfaces();

    // we want to only find overlap within a body
    SurfaceOverlapTool::instance()->set_check_within_bodies(CUBIT_TRUE);
    // 1=any, 2=opposite, 3=same  - we want to find only the overlaps that normals
    // pointing in the opposite directions
    SurfaceOverlapTool::instance()->set_normal_type(2);
    // Don't pickup surfaces that face each other
    SurfaceOverlapTool::instance()->set_skip_facing_surfaces(CUBIT_TRUE);

    // list of bodies that fail to midsurface
    DLIList<Body*> failing_bodies; 

    GeometryModifyEngine* gme = 0;
    GeometryQueryEngine* gqe = 0;

	// loop over every body and try to create midsurface(s)
	int i = 0;
	CubitStatus return_status = CUBIT_FAILURE;

    if(prog_tool)
        prog_tool->start(0,body_list_in.size());

    for(i = body_list_in.size();i--;)
    {
        if(prog_tool)
            prog_tool->step();

		Body* cur_body = body_list_in[i];
		if(!cur_body)
			continue;

		BodySM* body_sm = cur_body->get_body_sm_ptr();
		if(!body_sm)
			continue;

        if(cur_body->is_sheet_body())
        {
            PRINT_INFO("Body %d is a sheet body.\n",cur_body->id());
            continue;
        }

		// Grab the geometrymodify and geometryquery engines to use later
		gqe = cur_body->get_geometry_query_engine();
		gme = GeometryModifyTool::instance()->get_engine(body_sm);

		if(!gqe || !gme)
			continue;

		// Here are the steps to finding/creating the midsurface
		// 1. If the user did not give a thickness range to search then
		//    make an educated guess at the proper thickness range. The assumption
		//    is that the midsurface is a square and the thickness is constant.
		//    The resulting equation is a third order polynomial that is solved using
		//    a few newton iterations. The initial thickness guess is Volume/Area
		// 2. Using the given search distances use the SurfaceOverlapTool to find
		//    surface pairs.
        // 3. If there is only one surface pair then use the existing midsurface commands
		// 4. Find if the surface pairs represent two surface patches
		// 5. If there are only two surface patches try to offset one of the patches
		// 6. (this step is commented out for now) - If 5 fails or there are more than
        //           two surface patches then try the following:
		//         - Use the manual midsurface creation function to create midsurfaces for each
		//           pair of surfaces.
		//         - Unite all of the created midsurfaces together
		//         - remove any surfaces that have a curve touching a surface pair 
		//         - Regularize the resulting body
		// 7. Done

		{
			PRINT_DEBUG_198("AUTOMATICALLY calculating search range\n");
			DLIList<RefVolume*> vol_list;
			cur_body->ref_volumes(vol_list);
			double total_vol = 0;
            double total_vol_bb = 0;
			for(int vol_cnt = 0; vol_cnt < vol_list.size(); vol_cnt++)
			{
                CubitVector cg;
                double temp_volume;
                vol_list[vol_cnt]->mass_properties(cg,temp_volume);
                CubitBox vol_bb = vol_list[vol_cnt]->bounding_box();
                total_vol += temp_volume;
                total_vol_bb += vol_bb.x_range()*vol_bb.y_range()*vol_bb.z_range();
            }

            if(total_vol<0 || total_vol > total_vol_bb)
            {
                PRINT_INFO("Could not midsurface Body %d - try healing the body.\n",cur_body->id());
                failing_bodies.append(cur_body);
                continue;
            }
                
			PRINT_DEBUG_198("Volume of %f\n",total_vol);

			DLIList<RefFace*> face_list;
			cur_body->ref_faces(face_list);
			double total_surf = 0;
			for(int surf_cnt = 0; surf_cnt < face_list.size(); surf_cnt++)
				total_surf += face_list[surf_cnt]->area();
			PRINT_DEBUG_198("Area of %f\n",total_surf);

			double t_g = total_vol/(total_surf/2.0);
            double initial_guess = t_g;
			PRINT_DEBUG_198("Initial guess of thickness %f\n",t_g);
			// use a newton solver to get a more accurate estimate the thickness of the volume
			for(int n_i = 0;n_i<100;n_i++)
			{
				double tol_newton = GEOMETRY_RESABS;
				double t_gn = t_g + tol_newton;
				double f_prime = ((2.0*total_vol + sqrt(total_vol*t_g*t_g*t_g)*4.0 - total_surf*t_g)
					-(2.0*total_vol + sqrt(total_vol*t_gn*t_gn*t_gn)*4.0 - total_surf*t_gn))/
					(t_g-t_gn);

				// avoid divide by zero
				if(fabs(f_prime)<tol_newton)
					break;

				double t_old = t_g;
				t_g = t_g - (2.0*total_vol + sqrt(total_vol*t_g*t_g*t_g)*4.0 - total_surf*t_g)/f_prime;
				
				PRINT_DEBUG_198("Guess %d Thickness %f\n",n_i,t_g);
				if(fabs(t_g-t_old)<tol_newton)
				{
					PRINT_DEBUG_198("Converged with thickness of %f in %d steps\n",t_g,n_i);
					break;
				}
                if(t_g<0.0)
                {
					PRINT_DEBUG_198("thickness less than zero setting back to initial guess\n");
                    t_g = fabs(initial_guess);
					break;
                }
			}
			upper_tol = t_g + t_g*auto_thickness_margin;
			lower_tol = t_g - t_g*auto_thickness_margin;
            upper_tol = upper_tol <= upper_limit?upper_tol:upper_limit;
            lower_tol = lower_tol >= lower_limit?lower_tol:lower_limit;

			PRINT_DEBUG_198("Guessing a thickness of %f to %f\n",lower_tol,upper_tol);
		}

		// set the lower and upper search distances
		SurfaceOverlapTool::instance()->set_gap_max(upper_tol);
		SurfaceOverlapTool::instance()->set_gap_min(lower_tol);

		DLIList<RefFace*> ref_face_list,list1,list2;
        DLIList<RefEntity*> faces_to_draw; 
        cur_body->ref_faces(ref_face_list);
        // find the surface pairs
        SurfaceOverlapTool::instance()->find_overlapping_surfaces(ref_face_list,list1,list2,faces_to_draw);

        int tweak_iters = 4;
        for(int tweak = 0;tweak<tweak_iters;tweak++)
        {
            // if we didn't find anything then the part may be long and selender so grow the search thickness
            if(list1.size()==0 && list2.size() == 0)
            {
                if(tweak == tweak_iters-1 && lower_limit != -CUBIT_DBL_MAX && upper_limit != CUBIT_DBL_MAX)
                {
                    // on the last try use the user defined limits
                    lower_tol = lower_limit;
                    upper_tol = upper_limit;
                }
                else
                {
                    lower_tol = (upper_tol + lower_tol)/2.0;
                    upper_tol += lower_tol*auto_thickness_margin*2;
                    upper_tol = upper_tol <= upper_limit?upper_tol:upper_limit;
                    lower_tol = lower_tol >= lower_limit?lower_tol:lower_limit;
                }

                PRINT_DEBUG_198("Guessing again with thickness of %f to %f\n",lower_tol,upper_tol);
                SurfaceOverlapTool::instance()->set_gap_max(upper_tol);
                SurfaceOverlapTool::instance()->set_gap_min(lower_tol);
                SurfaceOverlapTool::instance()->find_overlapping_surfaces(ref_face_list,list1,list2,faces_to_draw);
            }

            DLIList<RefFace*> check_list;
            check_list += list1;
            check_list += list2;

            if(check_list.size() == 0 )
                continue;

            // make sure the pairs will match the solid within 10% or so
            if(!check_surf_pairs(lower_tol,upper_tol,check_list,cur_body))
            {
                list1.clean_out();
                list2.clean_out();
                continue;
            }
            break;
        }

        if(list1.size() != list2.size())
        {
            PRINT_INFO("Could not find workable surface pairs for Body %d - try using the Sheet Offset command. \n",cur_body->id());
            failing_bodies.append(cur_body);
            continue;
        }
		else if(list1.size() == 0 || list2.size() == 0)
		{
				PRINT_INFO("No surface pairs found for Body %d - try changing the search range\n",cur_body->id());
			failing_bodies.append(cur_body);
			continue;
		}

		// get the first pair and see if there are only two patches
		DLIList<RefFace*> red_faces;
		red_faces.append(list1[0]);
		DLIList<RefFace*> yellow_faces;
		yellow_faces.append(list2[0]);
		DLIList<RefFace*> paired_faces;
		paired_faces += list1;
		paired_faces += list2;
		paired_faces.uniquify_unordered();

		// red surfaces
		while(1)
		{
			int start_cnt = red_faces.size();
			DLIList<RefEdge*> red_edges;
			int j = 0;
			for(j =0;j<red_faces.size();j++)
				red_faces[j]->ref_edges(red_edges);
			red_edges.uniquify_unordered();
			for(j =0;j<red_edges.size();j++)
				red_edges[j]->ref_faces(red_faces);
			red_faces.uniquify_unordered();
			red_faces.intersect_unordered(paired_faces);
			if(start_cnt == red_faces.size())
				break;
		}

		// yellow surfaces
		while(1)
		{
			int start_cnt = yellow_faces.size();
			DLIList<RefEdge*> yellow_edges;
			int j = 0;
			for(j =0;j<yellow_faces.size();j++)
				yellow_faces[j]->ref_edges(yellow_edges);
			yellow_edges.uniquify_unordered();
			for(j =0;j<yellow_edges.size();j++)
				yellow_edges[j]->ref_faces(yellow_faces);
			yellow_faces.uniquify_unordered();
			yellow_faces.intersect_unordered(paired_faces);
			if(start_cnt == yellow_faces.size())
				break;
		}

        DLIList<BodySM*> results;
        bool midsurface_done = false;

        if(DEBUG_FLAG(198))
        {
            int j = 0;
            PRINT_INFO("Trying surface offset to create the mid_surface\n");
            PRINT_INFO("Red surface ");
            for(j = 0;j < red_faces.size();j++)
            {
                GfxDebug::draw_ref_face(red_faces[j],CUBIT_RED);
                PRINT_INFO("%d ",red_faces[j]->id());
            }

            PRINT_INFO("\nYellow surface ");
            for(j = 0;j < yellow_faces.size();j++)
            {
                GfxDebug::draw_ref_face(yellow_faces[j],CUBIT_YELLOW);
                PRINT_INFO("%d ",yellow_faces[j]->id());
            }

            PRINT_INFO("\n");
        }

        // first check to see if we can use the simple midsurface functions
        if(red_faces.size() == 1 && yellow_faces.size() == 1 &&
            paired_faces.size() == red_faces.size() + yellow_faces.size()) 
        {
            RefFace* face_1 = red_faces[0];
            RefFace* face_2 = yellow_faces[0];
            midsurface_done = false;

            if(face_1->geometry_type() == face_2->geometry_type())
            {
                Surface* surf_1 = face_1->get_surface_ptr();
                Surface* surf_2 = face_2->get_surface_ptr();
                BodySM* result_body;
                // grab the distance between surfaces
                CubitVector temp_vec0;
                CubitVector temp_vec1;
                double temp_dist = 0;
                gqe->entity_entity_distance(
                    face_1->get_surface_ptr(),
                    face_2->get_surface_ptr(),
                    temp_vec0,temp_vec1,temp_dist);

                switch(face_1->geometry_type())
                {
                case CONE_SURFACE_TYPE:
                    if(gme->get_conic_mid_surface(surf_1,surf_2,body_sm,result_body) == CUBIT_SUCCESS)
                    {
                        midsurface_done = true;
                        results.append(result_body);
                        thickness_out.append(fabs(temp_dist));
                    }
                    break;
                case PLANE_SURFACE_TYPE:
                    if(get_planar_mid_surface(face_1,face_2,body_sm,result_body,gme) == CUBIT_SUCCESS)
                    {
                        midsurface_done = true;
                        results.append(result_body);
                        thickness_out.append(fabs(temp_dist));
                    }
                    break;
                case SPHERE_SURFACE_TYPE:
                    if(gme->get_spheric_mid_surface(surf_1,surf_2,body_sm,result_body) == CUBIT_SUCCESS)
                    {
                        midsurface_done = true;
                        results.append(result_body);
                        thickness_out.append(fabs(temp_dist));
                    }
                    break;
                case TORUS_SURFACE_TYPE:
                    if(gme->get_toric_mid_surface(surf_1,surf_2,body_sm,result_body) == CUBIT_SUCCESS)
                    {
                        midsurface_done = true;
                        results.append(result_body);
                        thickness_out.append(fabs(temp_dist));
                    }
                    break;
                case CYLINDER_SURFACE_TYPE:
                    if(gme->get_conic_mid_surface(surf_1,surf_2,body_sm,result_body) == CUBIT_SUCCESS)
                    {
                        midsurface_done = true;
                        results.append(result_body);
                        thickness_out.append(fabs(temp_dist));
                    }
                    break;
                default:
                    break;
                }
            }
        }

		if(!midsurface_done &&
            paired_faces.size() == red_faces.size() + yellow_faces.size()) // just do the offset
		{
			int j = 0;
            DLIList<double> offset_distances;
            for(j = 0;j<list1.size();j++)
            {
                CubitVector temp_vec0;
                CubitVector temp_vec1;
				double temp_dist = 0;
				if(!gqe->entity_entity_distance(
					list1[j]->get_surface_ptr(),
					list2[j]->get_surface_ptr(),
					temp_vec0,temp_vec1,temp_dist))
				{
					break;
				}
				offset_distances.append(-temp_dist*.5);
			}

			DLIList<Surface*> red_surfs;
			for(j = 0;j<red_faces.size();j++)
				red_surfs.append(red_faces[j]->get_surface_ptr());

			DLIList<Surface*> yellow_surfs;
			for(j = 0;j<yellow_faces.size();j++)
				yellow_surfs.append(yellow_faces[j]->get_surface_ptr());

            // all of the surfaces are offset the same distance
			double offset_distance = offset_distances[0];
            bool old_error_flag = GET_ERROR_FLAG();
            SET_ERROR_FLAG(false); // don't throw any gme errors
			if( gme->create_offset_sheet(red_surfs,offset_distance,
				NULL,NULL,results))
            {
                midsurface_done = true;
                for(j = 0;j<results.size();j++) // for every body add a thickness
                    thickness_out.append(fabs(offset_distance*2.));
            }
            else if( gme->create_offset_sheet(yellow_surfs,offset_distance,
				NULL,NULL,results)) // try the other direction
            {
                midsurface_done = true;
                for(j = 0;j<results.size();j++) // for every body add a thickness
                    thickness_out.append(fabs(offset_distance*2.));
            }
            else
            {
                PRINT_INFO("Could not create midsurface for Body %d - try using the surface offset command\n",cur_body->id());
                failing_bodies.append(cur_body);
            }
            SET_ERROR_FLAG(old_error_flag); // turn errors back on
        }
        
        if(!midsurface_done && paired_faces.size() != red_faces.size() + yellow_faces.size())
        {
            PRINT_INFO("Could not find workable surface pairs for Body %d - try changing the search range or \n"
                "        using the Sheet Offset command.\n",cur_body->id());
        }

      /*if(!midsurface_done)
        {
			if(DEBUG_FLAG(198))
				PRINT_INFO("Trying the extend, unite, and trim method\n");

			// okay now remove duplicate pairs and unsupported pairs
			DLIList<Surface*> surf_list1;
			DLIList<Surface*> surf_list2;
			bool delete_and_exit = false;
			for(int j = 0;j<list1.size();j++)
			{
				RefFace* face_1 = list1[j];
				RefFace* face_2 = list2[j];

				if(DEBUG_FLAG(198))
				{
					PRINT_INFO("Red surface ");
					GfxDebug::draw_ref_face(face_1,CUBIT_RED);
					PRINT_INFO("%d ",face_1->id());

					PRINT_INFO("\nYellow surface ");
					GfxDebug::draw_ref_face(face_2,CUBIT_YELLOW);
					PRINT_INFO("%d ",face_2->id());

					PRINT_INFO("\n");
				}

				if(face_1->geometry_type() != 	face_2->geometry_type())
					continue;

				Surface* surf_1 = face_1->get_surface_ptr();
				surf_list1.append(surf_1);
				Surface* surf_2 = face_2->get_surface_ptr();
				surf_list2.append(surf_2);
				BodySM* result_body;
				switch(face_1->geometry_type())
				{
				case CONE_SURFACE_TYPE:
					if(gme->get_conic_mid_surface(surf_1,surf_2,body_sm,result_body) == CUBIT_SUCCESS)
						results.append(result_body);
					else
					    delete_and_exit = true;
					break;
				case PLANE_SURFACE_TYPE:
					if(get_planar_mid_surface(face_1,face_2,body_sm,result_body,gme) == CUBIT_SUCCESS)
						results.append(result_body);
					else
					    delete_and_exit = true;
					break;
				case SPHERE_SURFACE_TYPE:
					if(gme->get_spheric_mid_surface(surf_1,surf_2,body_sm,result_body) == CUBIT_SUCCESS)
						results.append(result_body);
					else
					    delete_and_exit = true;
					break;
				case TORUS_SURFACE_TYPE:
					if(gme->get_toric_mid_surface(surf_1,surf_2,body_sm,result_body) == CUBIT_SUCCESS)
						results.append(result_body);
					else
					    delete_and_exit = true;
					break;
				case CYLINDER_SURFACE_TYPE:
					if(gme->get_conic_mid_surface(surf_1,surf_2,body_sm,result_body) == CUBIT_SUCCESS)
						results.append(result_body);
					else
					    delete_and_exit = true;
					break;
				default:
					delete_and_exit = true;
					break;
				}

				if(delete_and_exit)
				{
					PRINT_WARNING("Failed to pair surface %d with surface %d\n",face_1->id(),face_2->id());
					break;
				}
			}

			if(delete_and_exit)
			{
				failing_bodies.append(cur_body);
    			gqe->delete_solid_model_entities(results);
				continue;
			}

			DLIList<BodySM*> unite_results;
			if(results.size()>1)
			{
				bool reg_result = GeometryModifyTool::instance()->boolean_regularize();
				GeometryModifyTool::instance()->boolean_regularize(true);
				if(gme->unite(results,unite_results)== CUBIT_SUCCESS)
				{
					// if the unite works just add them to the result list
					results = unite_results;
				}
				else
				{
					// clean up the created surfaces and move on to the next
					// body
					failing_bodies.append(cur_body);
					gqe->delete_solid_model_entities(results);
					GeometryModifyTool::instance()->boolean_regularize(reg_result);
					continue;
				}

				GeometryModifyTool::instance()->boolean_regularize(reg_result);
			}

			// trim the hanging surfaces 
			DLIList<Surface*> paired_surfs;
			paired_surfs += surf_list1;
			paired_surfs += surf_list2;

			DLIList<Curve*> all_curves;

			int k = 0;
			for(k = 0;k<results.size();k++)
				results[k]->curves(all_curves);

			all_curves.uniquify_unordered();

			DLIList<Surface*> remove_surfs;
			for(k = 0;k<all_curves.size();k++)
				for(int m = 0;m<paired_surfs.size();m++)
					if(curve_in_surface(all_curves[k],paired_surfs[m]))
						all_curves[k]->surfaces(remove_surfs);

			remove_surfs.uniquify_unordered();

			body_list_out += results;
			DLIList<BodySM*> tweak_results;
			if(gme->tweak_remove(remove_surfs,tweak_results,CUBIT_FALSE))
			{
				results = tweak_results;
			}
			else
			{
				// clean up the created surfaces and move on to the next
				// body
				failing_bodies.append(cur_body);
				gqe->delete_solid_model_entities(results);
				continue;
			}

			DLIList<BodySM*> regularize_results;
			// regularize the results
			for(k = 0;k < results.size();k++)
			{
				BodySM* new_body = 0;
				if(gme->regularize_body(results[k],new_body))
					regularize_results.append(new_body);
				else if(DEBUG_FLAG(198))
					PRINT_INFO("Regularize failure\n");
			}
            results = regularize_results;
        }*/

        if(!midsurface_done)
        {
           failing_bodies.append(cur_body);
           continue;
        }

        old_bodies_midsurfaced.append(cur_body);

        if(delete_midsurfaced && !preview)
            GeometryQueryTool::instance()->delete_Body(cur_body);

        return_status = CUBIT_SUCCESS;
        body_list_out += results;
    }

    if(prog_tool)
        prog_tool->end();

    PRINT_INFO("Successfully midsurface %d of %d bodies\n",body_list_out.size(),body_list_in.size());
    if(preview)
    {
        for(int k = 0;k<body_list_out.size();k++)
        {
            DLIList<Surface*> preview_surfaces;
            body_list_out[k]->surfaces(preview_surfaces);
            for(int p = 0;p<preview_surfaces.size();p++)
                GfxPreview::draw_surface_facets_shaded(preview_surfaces[p],CUBIT_BLUE);
        }
        GfxPreview::flush();
        if(gqe)
            gqe->delete_solid_model_entities(body_list_out);
        body_list_out.clean_out();
    }

	if(failing_bodies.size() > 0)
	{
        PRINT_INFO("\n");
		PRINT_INFO("Failed to midsurface Body ");
		for(i = 0;i<failing_bodies.size();i++)
			PRINT_INFO("%d ",failing_bodies[i]->id());
		PRINT_INFO("\n");
	}

	if(DEBUG_FLAG(198))
		GfxDebug::flush();

	SurfaceOverlapTool::instance()->set_check_within_bodies(cubit_bool_save);
	SurfaceOverlapTool::instance()->set_gap_max(max_gap_save);
	SurfaceOverlapTool::instance()->set_normal_type(normal_type);
	SurfaceOverlapTool::instance()->set_gap_min(min_gap_save);
    SurfaceOverlapTool::instance()->set_skip_facing_surfaces(skip_facing_surfaces);

	return return_status;
}


CubitStatus AutoMidsurfaceTool::get_planar_mid_surface( RefFace* ref_face1,
                                                        RefFace* ref_face2,
                                                        BodySM* body_sm_to_trim_to,
                                                        BodySM*& midsurface_body_sm,
                                                        GeometryModifyEngine *gme_ptr )
{
	CubitVector normal_1, normal_2, point_1, point_2, point_3;
	CubitPlane plane_1, plane_2;
	CubitVector p_mid, n_mid;

	point_1 = ref_face1->center_point();
	point_2 = ref_face2->center_point();

	normal_1 = ref_face1->normal_at(point_1);
	normal_2 = ref_face2->normal_at(point_2);

	plane_1 = CubitPlane(normal_1,point_1);
	plane_2 = CubitPlane(normal_2,point_2);

	if(point_1 == point_2)
	{
		PRINT_ERROR( "In GeometryModifyTool:: get_planar_mid_surface\n"
			"       Since both surfaces share the same point, the midsurface is not well-defined\n");
		return CUBIT_FAILURE;
	}
	else
	{
		CubitVector temp1 = point_2;
		temp1 = plane_1.project(temp1);
		temp1 -= point_2;
		if ( temp1.length_squared() < GEOMETRY_RESABS*GEOMETRY_RESABS )
		{
			PRINT_ERROR("In GeometryModifyTool:: get_planar_mid_surface\n"
				"       Since both planes are the same, the midsurface is not well-defined.\n");
			return CUBIT_FAILURE;
		}
	}

	if ( ( normal_1.about_equal( normal_2 ) ) || ( (-normal_1).about_equal( normal_2 ) ) )
	{
		p_mid = (point_1+point_2)/2;
		n_mid = plane_1.normal();
	}
	else
	{
		CubitVector direction_of_line;
		plane_1.intersect(plane_2,p_mid,direction_of_line);
		direction_of_line.normalize();

		// Find if point_1 and point_2 are on the line of intersection
		// If they are, then the mid-plane is not well-defined
		CubitVector p1 = point_1-p_mid;
		CubitVector p2 = point_2-p_mid;
		p1.normalize();
		p2.normalize();

		if(p1==direction_of_line || p1==-direction_of_line)
		{
			PRINT_ERROR("In GeometryModifyTool:: get_planar_mid_surface\n"
				"       P1 is on the line of intersection.\n");
			return CUBIT_FAILURE;
		}

		if(p2==direction_of_line || p2==-direction_of_line)
		{
			PRINT_ERROR("In GeometryModifyTool:: get_planar_mid_surface\n"
				"       P2 is on the line of intersection.\n");
			return CUBIT_FAILURE;
		}

		CubitVector v1 = p1 - (p1%direction_of_line)*direction_of_line;
		v1.normalize();

		CubitVector v2 = p2 - (p2%direction_of_line)*direction_of_line;
		v2.normalize();

		n_mid = v1 - v2;
		n_mid.normalize();
	}

	CubitPlane mid_plane(n_mid, p_mid);
	point_1 = p_mid;

	//find three points that will define the infinite plane from the
	//mid plane.through the point in any direction just not along the
	//normal direction
	CubitVector Xdir(1,0,0), Ydir(0,1,0);
	CubitVector direction1;

	if ( ( ! n_mid.about_equal( Xdir ) ) && ( ! (-n_mid).about_equal( Xdir ) ) )
		direction1 = Xdir + n_mid;
	else
		direction1 = Ydir + n_mid;

	point_2 = p_mid + direction1;
	point_2 = mid_plane.project(point_2);

	direction1 = point_2-point_1;
	CubitVector direction2 = direction1*n_mid;
	point_3 = point_1 + direction2;

	CubitStatus ret = gme_ptr->get_mid_plane(point_1, point_2, point_3,
		body_sm_to_trim_to, midsurface_body_sm );
	return ret;
}

CubitBoolean AutoMidsurfaceTool::curve_in_surface(Curve *curve_in, Surface *surf_in)
{
	CubitVector loc_0;
	CubitVector loc_1;
	CubitVector loc_2;

	curve_in->position_from_fraction(0.1,loc_0);
	curve_in->position_from_fraction(0.5,loc_1);
	curve_in->position_from_fraction(0.9,loc_2);

	GeometryQueryEngine* gqe = surf_in->get_geometry_query_engine();
	double tol = gqe->get_sme_resabs_tolerance();
	CubitVector cl_pnt_0;
	CubitVector cl_pnt_1;
	CubitVector cl_pnt_2;
	surf_in->closest_point(loc_0,&cl_pnt_0);
	surf_in->closest_point(loc_1,&cl_pnt_1);
	surf_in->closest_point(loc_2,&cl_pnt_2);

	if(cl_pnt_0.distance_between(loc_0)<tol &&
		cl_pnt_1.distance_between(loc_1)<tol &&
		cl_pnt_2.distance_between(loc_2)<tol)
	{
		return CUBIT_TRUE;
	}

	return CUBIT_FALSE;
}

CubitStatus AutoMidsurfaceTool::find_offset_pair_patches(
	DLIList<RefFace*> pairs_list_0,
	DLIList<RefFace*> pairs_list_1,	
	DLIList<RefFace*>& red_faces,
	DLIList<RefFace*>& yellow_faces,
	DLIList<double>& offset_distances)
{
	return CUBIT_FAILURE;
}

CubitStatus AutoMidsurfaceTool::random_loc_on_surface( Surface* face_ptr, CubitVector &loc )
{
  GMem g_mem;
  GeometryQueryEngine* gqe = face_ptr->get_geometry_query_engine();
  unsigned short norm_tol = 10;
  double dist_tol = -1.0;
  gqe->get_graphics( face_ptr, &g_mem, norm_tol, dist_tol );

  if(g_mem.fListCount < 1)
  {
    // Decrease tolerance and try again (we can get this for small features)
    norm_tol /= 2;
    gqe->get_graphics( face_ptr, &g_mem, norm_tol, dist_tol);
  }

  if(g_mem.fListCount < 1)
  {
    // Lets give up
    PRINT_ERROR( "Unable to find location on a surface\n" );
    return CUBIT_FAILURE;
  }

  // Use the first triangle
  GPoint p[3];
  GPoint* plist = g_mem.point_list();
  int* facet_list = g_mem.facet_list();
  int c = 0;

  p[0] = plist[facet_list[++c]];
  p[2] = plist[facet_list[++c]];
  p[1] = plist[facet_list[++c]];

  // Get centroid
  CubitVector p1( p[0].x, p[0].y, p[0].z );
  CubitVector p2( p[2].x, p[2].y, p[2].z );
  CubitVector p3( p[1].x, p[1].y, p[1].z );

  CubitVector center = (p1 + p2 + p3)/3.0;

  face_ptr->closest_point_trimmed(center,loc);

  return CUBIT_SUCCESS;
}

CubitBoolean AutoMidsurfaceTool::check_surf_pairs(double min_thick, double max_thick,
                                                 DLIList<RefFace*> check_list, Body* body_in )
{
    double total_area = 0.0;
    DLIList<RefVolume*> vol_list;
    body_in->ref_volumes(vol_list);
    double total_vol = 0;
    for(int vol_cnt = 0; vol_cnt < vol_list.size(); vol_cnt++)
    {
        CubitVector cg;
        double temp_volume;
        vol_list[vol_cnt]->mass_properties(cg,temp_volume);
        total_vol += temp_volume;
    }

    for(int i = 0;i<check_list.size();i++)
        total_area += check_list[i]->area();

    total_area/=2.0;

    if(min_thick*total_area < total_vol && max_thick*total_area > total_vol)
        return CUBIT_TRUE;

    return CUBIT_FALSE;
}
