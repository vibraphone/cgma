#include "SplitCompositeSurfaceTool.hpp"
#include "GeometryModifyTool.hpp"
#include "RefFace.hpp"
#include "CompositeSurface.hpp"
#include "CompositeCurve.hpp"
#include "Curve.hpp"
#include "Point.hpp"
#include "GeometryModifyEngine.hpp"
#include "GeometryQueryTool.hpp"

SplitCompositeSurfaceTool *SplitCompositeSurfaceTool::instance_ = NULL;

SplitCompositeSurfaceTool *SplitCompositeSurfaceTool::instance()
{
  if (instance_ == NULL) 
    instance_ = new SplitCompositeSurfaceTool();
                                                                                
  return instance_;
}

CubitStatus SplitCompositeSurfaceTool::split_surface(RefFace *ref_face_ptr,
                                                     DLIList<CubitVector*> &locations,
                                                     DLIList<DLIList<CubitVector*>*> &vec_lists,
                                                     CubitBoolean preview_flg,
                                                     CubitBoolean create_ref_edges_flg,
                                                     CubitBoolean clear_previous_previews)
{
   get_additional_split_points(ref_face_ptr->get_surface_ptr(), vec_lists);

   return GeometryModifyTool::instance()->split_surface(ref_face_ptr,
      locations, vec_lists, preview_flg, create_ref_edges_flg, clear_previous_previews);
}

CubitStatus SplitCompositeSurfaceTool::split_surface(DLIList<RefFace*> &ref_face_list,
                                                     DLIList<CubitVector*> &locations,
                                                     DLIList<DLIList<DLIList<CubitVector*>*>*> &list_of_vec_lists,
                                                     CubitBoolean preview_flg,
                                                     CubitBoolean create_ref_edges_flg,
                                                     CubitBoolean clear_previous_previews)
{
   for( int jj = ref_face_list.size() ; jj > 0 ; jj--)
   {
      DLIList<DLIList<CubitVector*>*> vec_lists = *( list_of_vec_lists.get_and_step() );
      
      get_additional_split_points(ref_face_list.get_and_step()->get_surface_ptr(), vec_lists );
   }

   return GeometryModifyTool::instance()->split_surface(ref_face_list,
      locations, list_of_vec_lists, preview_flg, create_ref_edges_flg, clear_previous_previews);
}

void SplitCompositeSurfaceTool::find_faces_for_pos(CubitVector &pos, DLIList<Surface*> surf_list,
                        CubitPointContainment &containment,
                        DLIList<Surface*> &out_list)
{
  int j;

  for(j=surf_list.size(); j>0; j--)
  {
    Surface *surf_ptr = surf_list.get_and_step();
    containment = surf_ptr->point_containment(pos);
    if(containment == CUBIT_PNT_INSIDE || containment == CUBIT_PNT_BOUNDARY)
      out_list.append_unique(surf_ptr);
  }
}

void SplitCompositeSurfaceTool::get_additional_split_points(Surface *surf,
                            DLIList<DLIList<CubitVector*>*> &vec_lists)
{
  int i, j, k, num_surfs;
  DLIList<Surface*> surf_list;

  CompositeSurface *cs = dynamic_cast<CompositeSurface*>(surf);
  if(cs && (num_surfs = cs->num_surfs()) > 1)
  {
    for(i=0; i<num_surfs; i++)
    {
      Surface *surf = cs->get_surface(i);
      if(surf)
        surf_list.append(surf);
    }

    DLIList<Curve*> hidden_curves;
    cs->get_hidden_curves(hidden_curves);

    for(i=vec_lists.size(); i--;)
    {
      // Remove any hidden curves that already have points on them.
      DLIList<Curve*> hidden_curves_without_pts = hidden_curves;
      DLIList<CubitVector*> *vec_list = vec_lists.get_and_step();
      for(j=vec_list->size(); j--;)
      {
        CubitVector *cur_vec = vec_list->get_and_step();
        DLIList<Curve*> crvs_to_remove;
        for(k=hidden_curves_without_pts.size(); k>0; k--)
        {
          Curve *crv = hidden_curves_without_pts.get_and_step();
          CubitVector foot;
          crv->closest_point_trimmed(*cur_vec, foot);
          if(foot.about_equal(*cur_vec))
            crvs_to_remove.append_unique(crv);
        }
        hidden_curves_without_pts -= crvs_to_remove;
      }

      // Now we should only have hidden curves with no points on them.
      if(hidden_curves_without_pts.size() > 0)
      {
        vec_list->reset();
        CubitVector *cur_vec = vec_list->get();
        CubitVector *last_vec = vec_list->prev();
        while(cur_vec != last_vec)
        {
          vec_list->move_to(cur_vec);
          CubitVector *next_vec = vec_list->next();
          DLIList<CubitVector*> tmp_list;
          DLIList<double> vals;
          tmp_list.append(cur_vec);
          tmp_list.append(next_vec);
          vals.append(0.0);
          vals.append(1.0);
          Curve *hidden_curve = hidden_curves_without_pts.get();
          CompositeCurve *cc = dynamic_cast<CompositeCurve*>(hidden_curve);
          if(cc)
            hidden_curve = cc->get_curve(0);
          GeometryModifyEngine *gme = GeometryModifyTool::instance()->
            get_engine(hidden_curve);
          Point *pt1 = gme->make_Point(*cur_vec);
          Point *pt2 = gme->make_Point(*last_vec);
          Curve *tmp_crv = gme->make_Curve(STRAIGHT_CURVE_TYPE, pt1, pt2, NULL, CUBIT_FORWARD);
          if(tmp_crv)
          {
            double arc_length = tmp_crv->get_arc_length();
            for(j=hidden_curves_without_pts.size(); j--;)
            {
              CubitVector pos1, pos2;
              double dist;
              Curve *crv = hidden_curves_without_pts.get_and_step();
              GeometryQueryTool::instance()->entity_entity_distance(crv, tmp_crv, pos1, pos2, dist);
              CubitVector v1 = *cur_vec - pos1;
              CubitVector v2 = *last_vec - pos1;
              v1.normalize();
              v2.normalize();
              if(v1 % v2 < .3)
              {
                CubitVector *new_vec = new CubitVector(pos1);
                v1 = (pt1->coordinates() - pos2);
                double percent = v1.length()/arc_length;
                if(percent > 1.0)
                  percent = 1.0;
                if(percent < 0.0)
                  percent = 0.0;
                tmp_list.reset();
                vals.reset();
                while(percent > vals.get())
                {
                  vals.step();
                  tmp_list.step();
                }
                vals.back();
                tmp_list.back();
                vals.insert(percent);
                tmp_list.insert(new_vec);
              }
            }
            delete tmp_crv;
          }
          delete pt1;
          delete pt2;
          if(tmp_list.size() > 2)
          {
            vec_list->move_to(cur_vec);
            tmp_list.reset();
            tmp_list.step();
            for(j=tmp_list.size(); j>2; j--)
            {
              vec_list->insert(tmp_list.get());
              tmp_list.step();
            }
          }
          cur_vec = next_vec;
        }
      }
    }
  }
}

// Passed in direction should be normalized by calling function.
void SplitCompositeSurfaceTool::find_face_with_non_zero_param_dir(DLIList<Surface*> &surf_list, 
                                        CubitVector &dir,
                                        CubitVector &pos,
                                        Surface *&ret_surf,
                                        double &du, double &dv,
                                        double &step)
{
  int i;
  ret_surf = NULL;
  for(i=surf_list.size(); i>0 && !ret_surf; i--)
  {
    Surface *surf = surf_list.get_and_step();
    CubitPointContainment cont = surf->point_containment(pos);
    if(cont == CUBIT_PNT_INSIDE)
    {
      double ulow, uhigh, vlow, vhigh;
      surf->param_dir(dir, pos, du, dv);
      surf->get_param_range_U(ulow, uhigh);
      surf->get_param_range_V(vlow, vhigh);
      double u_len = fabs(ulow-uhigh);
      double v_len = fabs(vlow-vhigh);
      double par_tol;
      if(u_len > v_len)
        par_tol = v_len/100.0;
      else
        par_tol = u_len/100.0;
      par_tol *= par_tol;
      double len = du*du + dv*dv;
      if(len > par_tol)
      {
        step = sqrt(u_len*u_len + v_len*v_len)/20.0;
        ret_surf = surf;
      }
    }
    else if(cont == CUBIT_PNT_BOUNDARY)
    {
      double ulow, uhigh, vlow, vhigh;
      surf->param_dir(dir, pos, du, dv);
      surf->get_param_range_U(ulow, uhigh);
      surf->get_param_range_V(vlow, vhigh);
      double u_len = fabs(ulow-uhigh);
      double v_len = fabs(vlow-vhigh);
      double par_tol;
      if(u_len > v_len)
        par_tol = v_len/100.0;
      else
        par_tol = u_len/100.0;
      par_tol *= par_tol;
      double len = du*du + dv*dv;
      if(len > par_tol)
      {
        double u, v;
        surf->u_v_from_position(pos, u, v);
        CubitVector new_pos = surf->position_from_u_v(u+du, v+dv);
        CubitPointContainment new_cont = surf->point_containment(new_pos);
        if(new_cont == CUBIT_PNT_INSIDE)
        {
          step = sqrt(u_len*u_len + v_len*v_len)/20.0;
          ret_surf = surf;
        }
      }
    }
  }
}
