#include <assert.h>

#include "CompositeEngine.hpp"
#include "PartitionEngine.hpp"

#include "DLIList.hpp"
#include "TDUniqueId.hpp"
#include "CubitTransformMatrix.hpp"
#include "RTree.hpp"

#include "CompositePoint.hpp"
#include "CompositeCurve.hpp"
#include "CompositeCoEdge.hpp"
#include "CompositeLoop.hpp"
#include "CompositeSurface.hpp"
#include "CompositeShell.hpp"
#include "CompositeLump.hpp"
#include "CompositeBody.hpp"
#include "GfxPreview.hpp"

#include "PartitionPoint.hpp"
#include "SegmentedCurve.hpp"

#include "CADefines.hpp"

#include "VGLoopTool.hpp"
#include "Body.hpp"
#include "LumpSM.hpp"

#include "GeometryQueryTool.hpp"
#include "AppUtil.hpp"

typedef VGLoopTool<CompositeSurface,
                   CompositeLoop,
                   CompositeCoEdge,
                   CompositeCurve,
                   CompositePoint> CompLoopTool;

CompositeEngine* CompositeEngine::instance_ = NULL; 


CompositeEngine::~CompositeEngine()
{
  GeometryQueryTool::instance()->unregister_intermediate_engine(this);
}

void CompositeEngine::delete_instance()
{
  if( NULL != instance_ )
  {
    delete instance_;
    instance_ = NULL;
  }
}

CompositeEngine& CompositeEngine::instance()
{
  if( instance_ == NULL )
  {
    instance_ = new CompositeEngine();
    assert( instance != NULL );
  }
  
  return *instance_;
}

// This function goes through all of the curves and vertices and
// removes any named attributes specified below.
void CompositeEngine::remove_imprint_attributes_after_modify( DLIList<BodySM*> &old_sms,
                                DLIList<BodySM*> &new_sms )
{
  int k, m, q, w, g, b, s, t;
  CubitString name("IMPRINT_PREEXISTING");
  std::vector<CubitString> string_list;
  string_list.push_back( name );
  CubitSimpleAttrib geom_attrib( &string_list, 0, 0 );

  DLIList<BodySM*> all_sms = old_sms;
  all_sms += new_sms;
  DLIList<TopologyBridge*> top_bridges;
  CAST_LIST_TO_PARENT(all_sms, top_bridges);
  for(k=top_bridges.size(); k--;)
  {
    TopologyBridge *cur_body = top_bridges.get_and_step();
    DLIList<TopologyBridge*> lumps;
    cur_body->get_children_virt(lumps);
    for(m=lumps.size(); m--;)
    {
      TopologyBridge *cur_lump = lumps.get_and_step();
      DLIList<TopologyBridge*> shells;
      cur_lump->get_children_virt(shells);
      for(q=shells.size(); q--;)
      {
        TopologyBridge *cur_shell = shells.get_and_step();
        DLIList<TopologyBridge*> surfaces;
        cur_shell->get_children_virt(surfaces);
        for(w=surfaces.size(); w--;)
        {
          TopologyBridge *cur_surface = surfaces.get_and_step();
          DLIList<TopologyBridge*> loops;
          cur_surface->get_children_virt(loops);
          for(g=loops.size(); g--;)
          {
            TopologyBridge *cur_loop = loops.get_and_step();
            DLIList<TopologyBridge*> coedges;
            cur_loop->get_children_virt(coedges);
            for(b=coedges.size(); b--;)
            {
              TopologyBridge *cur_coedge = coedges.get_and_step();
              DLIList<TopologyBridge*> curves;
              cur_coedge->get_children_virt(curves);
              for(s=curves.size(); s--;)
              {
                TopologyBridge *cur_curve = curves.get_and_step();
                DLIList<CubitSimpleAttrib> list;
                cur_curve->get_simple_attribute("IMPRINT_PREEXISTING",list);
                if(list.size() != 0)
                  cur_curve->remove_simple_attribute_virt(list.get());
                DLIList<TopologyBridge*> pts;
                cur_curve->get_children_virt(pts);
                for(t=pts.size(); t--;)
                {
                  TopologyBridge *cur_pt = pts.get_and_step();
                  list.clean_out();
                  cur_pt->get_simple_attribute("IMPRINT_PREEXISTING",list);
                  if(list.size() != 0)
                    cur_pt->remove_simple_attribute_virt(list.get());
                }
              }
            }
          }
        }
      }
    }
  }
}

void CompositeEngine::push_imprint_attributes_before_modify
                     ( DLIList<BodySM*> &bodies)
{
}

void CompositeEngine::push_named_attributes_to_curves_and_points
//                     ( DLIList<BodySM*> &bodies, char *name_in)
                     ( DLIList<TopologyBridge*> &in_list, const char *name_in)
{
  int i/*, k, m, q, w, g, b, s, t*/;
  CubitString name(name_in);
  std::vector<CubitString> string_list;
  string_list.push_back( name );
  CubitSimpleAttrib attrib( &string_list, 0, 0 );

  for(i=in_list.size(); i>0; i--)
  {
    TopologyBridge *tb = in_list.get_and_step();
    if(dynamic_cast<BodySM*>(tb))
    {
      DLIList<TopologyBridge*> lumps;
      tb->get_children_virt(lumps);
      push_named_attributes_to_curves_and_points(lumps, name_in);
    }
    else if(dynamic_cast<Lump*>(tb))
    {
      DLIList<TopologyBridge*> shells;
      tb->get_children_virt(shells);
      push_named_attributes_to_curves_and_points(shells, name_in);
    }
    else if(dynamic_cast<ShellSM*>(tb))
    {
      DLIList<TopologyBridge*> surfs;
      tb->get_children_virt(surfs);
      push_named_attributes_to_curves_and_points(surfs, name_in);
    }
    else if(dynamic_cast<Surface*>(tb))
    {
      DLIList<TopologyBridge*> loops;
      tb->get_children_virt(loops);
      push_named_attributes_to_curves_and_points(loops, name_in);
    }
    else if(dynamic_cast<LoopSM*>(tb))
    {
      DLIList<TopologyBridge*> coedges;
      tb->get_children_virt(coedges);
      push_named_attributes_to_curves_and_points(coedges, name_in);
    }
    else if(dynamic_cast<CoEdgeSM*>(tb))
    {
      DLIList<TopologyBridge*> curves;
      tb->get_children_virt(curves);
      push_named_attributes_to_curves_and_points(curves, name_in);
    }
    else if(dynamic_cast<Curve*>(tb))
    {
      append_attrib( tb, attrib );
      DLIList<TopologyBridge*> points;
      tb->get_children_virt(points);
      push_named_attributes_to_curves_and_points(points, name_in);
    }
    else if(dynamic_cast<TBPoint*>(tb))
    {
      append_attrib( tb, attrib );
    }
  }

/*


  DLIList<TopologyBridge*> top_bridges;
  CAST_LIST_TO_PARENT(bodies, top_bridges);
  for(k=top_bridges.size(); k--;)
  {
    TopologyBridge *cur_body = top_bridges.get_and_step();
    DLIList<TopologyBridge*> lumps;
    cur_body->get_children_virt(lumps);
    for(m=lumps.size(); m--;)
    {
      TopologyBridge *cur_lump = lumps.get_and_step();
      DLIList<TopologyBridge*> shells;
      cur_lump->get_children_virt(shells);
      for(q=shells.size(); q--;)
      {
        TopologyBridge *cur_shell = shells.get_and_step();
        DLIList<TopologyBridge*> surfaces;
        cur_shell->get_children_virt(surfaces);
        for(w=surfaces.size(); w--;)
        {
          TopologyBridge *cur_surface = surfaces.get_and_step();
          DLIList<TopologyBridge*> loops;
          cur_surface->get_children_virt(loops);
          for(g=loops.size(); g--;)
          {
            TopologyBridge *cur_loop = loops.get_and_step();
            DLIList<TopologyBridge*> coedges;
            cur_loop->get_children_virt(coedges);
            for(b=coedges.size(); b--;)
            {
              TopologyBridge *cur_coedge = coedges.get_and_step();
              DLIList<TopologyBridge*> curves;
              cur_coedge->get_children_virt(curves);
              for(s=curves.size(); s--;)
              {
                TopologyBridge *cur_curve = curves.get_and_step();
                append_attrib( cur_curve, &attrib );
                DLIList<TopologyBridge*> pts;
                cur_curve->get_children_virt(pts);
                for(t=pts.size(); t--;)
                {
                  TopologyBridge *cur_pt = pts.get_and_step();
                  append_attrib( cur_pt, &attrib );
                }
              }
            }
          }
        }
      }
    }
  }
  */
}

void CompositeEngine::get_all_curves_and_points(DLIList<TopologyBridge*> &tb_list,
                                                DLIList<Curve*> &curves,
                                                DLIList<TBPoint*> &points)
{
  int i;
  Curve *crv;
  TBPoint *pt;
  for(i=tb_list.size(); i>0; i--)
  {
    TopologyBridge *tb = tb_list.get_and_step();
    if(dynamic_cast<BodySM*>(tb))
    {
      DLIList<TopologyBridge*> lumps;
      tb->get_children_virt(lumps);
      get_all_curves_and_points(lumps, curves, points);
    }
    else if(dynamic_cast<Lump*>(tb))
    {
      DLIList<TopologyBridge*> shells;
      tb->get_children_virt(shells);
      get_all_curves_and_points(shells, curves, points);
    }
    else if(dynamic_cast<ShellSM*>(tb))
    {
      DLIList<TopologyBridge*> surfs;
      tb->get_children_virt(surfs);
      get_all_curves_and_points(surfs, curves, points);
    }
    else if(dynamic_cast<Surface*>(tb))
    {
      DLIList<TopologyBridge*> loops;
      tb->get_children_virt(loops);
      get_all_curves_and_points(loops, curves, points);
    }
    else if(dynamic_cast<LoopSM*>(tb))
    {
      DLIList<TopologyBridge*> coedges;
      tb->get_children_virt(coedges);
      get_all_curves_and_points(coedges, curves, points);
    }
    else if(dynamic_cast<CoEdgeSM*>(tb))
    {
      DLIList<TopologyBridge*> tmp_curves;
      tb->get_children_virt(tmp_curves);
      get_all_curves_and_points(tmp_curves, curves, points);
    }
    else if((crv = dynamic_cast<Curve*>(tb)))
    {
      curves.append(crv);
      DLIList<TopologyBridge*> tmp_points;
      tb->get_children_virt(tmp_points);
      get_all_curves_and_points(tmp_points, curves, points);
    }
    else if((pt = dynamic_cast<TBPoint*>(tb)))
    {
      points.append(pt);
    }
  }
}

// Function to apply/remove COMPOSITE_GEOM attributes as necessary based
// on imprinting.
void CompositeEngine::attribute_after_imprinting(DLIList<TopologyBridge*> &tb_list,
                                                 DLIList<Body*> &old_bodies)
{
  DLIList<TopologyBridge*> all_bridges = tb_list;
  int i, j, k;
  for(k = old_bodies.size(); k>0; k--)
  {
    Body *body = old_bodies.get_and_step();
    TopologyBridge *tb = body->bridge_manager()->topology_bridge();
    if(tb)
      all_bridges.append_unique(tb);
  }

  DLIList<Curve*> all_curves;
  DLIList<TBPoint*> all_points;
  get_all_curves_and_points(all_bridges, all_curves, all_points);
  all_curves.uniquify_ordered();
  all_points.uniquify_ordered();

  double geom_factor = GeometryQueryTool::get_geometry_factor();
  double merge_tol = geom_factor*GEOMETRY_RESABS;

  AbstractTree<TBPoint*> *pt_tree = new RTree<TBPoint*>(merge_tol);
  AbstractTree<Curve*> *crv_tree = new RTree<Curve*>(merge_tol);

  DLIList<Curve*> all_curves_with_composite_att;
  DLIList<TBPoint*> all_points_with_composite_att;
  for(k=all_curves.size(); k>0; k--)
  {
    Curve *cur_curve = all_curves.get_and_step();
    crv_tree->add(cur_curve);
    DLIList<CubitSimpleAttrib> list;
    cur_curve->get_simple_attribute("COMPOSITE_GEOM",list);
    if(list.size() > 0)
      all_curves_with_composite_att.append(cur_curve);
  }
  for(k=all_points.size(); k>0; k--)
  {
    TBPoint *cur_point = all_points.get_and_step();
    pt_tree->add(cur_point);
    DLIList<CubitSimpleAttrib> list;
    cur_point->get_simple_attribute("COMPOSITE_GEOM",list);
    if(list.size() > 0)
      all_points_with_composite_att.append(cur_point);
  }

  DLIList<CubitSimpleAttrib> list;
  while(all_points_with_composite_att.size())
  {
    DLIList<TBPoint*> other_pts;
    DLIList<BodySM*> other_bodies;
    DLIList<double> other_distances;

    // For the given pt we will look for "coincident" pts (those within merge tol)
    // and categorize them as either having or not having a composite att.
    TBPoint *cur_pt = all_points_with_composite_att.extract();
    pt_tree->remove(cur_pt);

    BodySM *cur_body = cur_pt->bodysm();
    DLIList<TBPoint*> coincident_pts_with_composite_att, coincident_pts_without_composite_att;
    DLIList<TBPoint*> close_pts;
    CubitBox bbox = cur_pt->bounding_box();
    pt_tree->find(bbox, close_pts);

    // Only keep the closest pt from each body.
    for(j=close_pts.size(); j>0; j--)
    {
      TBPoint *other_pt = close_pts.get_and_step();
      BodySM *other_body = other_pt->bodysm();
      // Don't keep anything that is in the same body as the current pt.
      if(other_body != cur_body)
      {
        double cur_dist_sq = cur_pt->coordinates().distance_between_squared(other_pt->coordinates());
        if(other_bodies.move_to(other_body))
        {
          int list_index = other_bodies.get_index();
          other_distances.reset();
          other_distances.step(list_index);
          double prev_dist_sq = other_distances.get();
          if(cur_dist_sq < prev_dist_sq)
          {
            other_distances.change_to(cur_dist_sq);
            other_pts.reset();
            other_pts.step(list_index);
            other_pts.change_to(other_pt);
          }
        }
        else
        {
          other_bodies.append(other_body);
          other_pts.append(other_pt);
          other_distances.append(cur_dist_sq);
        }
      }
    }
    // Make sure our current pt is added to a list.
    coincident_pts_with_composite_att.append(cur_pt);
    // Classify the coincident pts as either having or not
    // having a composite att.
    for(j=other_pts.size(); j>0; j--)
    {
      TBPoint *pt = other_pts.get_and_step();
      list.clean_out();
      pt->get_simple_attribute("COMPOSITE_GEOM",list);
      if(list.size() > 0)
      {
        coincident_pts_with_composite_att.append(pt);
        if(all_points_with_composite_att.move_to(pt))
          all_points_with_composite_att.extract();
      }
      else
        coincident_pts_without_composite_att.append(pt);
    }

    // If we have found at least one other pt coincident with the current point...
    if(coincident_pts_with_composite_att.size() > 1 ||
      coincident_pts_without_composite_att.size() > 0)
    {
      // If there is at least one pt without a composite att that is an imprinter we
      // will remove all composite atts from coincident pts
      bool found = false;
      for(j=coincident_pts_without_composite_att.size(); j>0 && !found; j--)
      {
        TBPoint *tmp_pt = coincident_pts_without_composite_att.get_and_step();
        list.clean_out();
        tmp_pt->get_simple_attribute("IMPRINTER",list);
        if(list.size() > 0)
          found = true;
      }
      if(found)
      {
        // Remove all composite atts.
        for(j=coincident_pts_with_composite_att.size(); j>0; j--)
        {
          TBPoint *tmp_pt = coincident_pts_with_composite_att.get_and_step();
          list.clean_out();
          tmp_pt->get_simple_attribute("COMPOSITE_GEOM",list);
          if(list.size() > 0)
            tmp_pt->remove_simple_attribute_virt(list.get());
        }
      }
      else
      {
        // There were no imprinter points that didn't have composite atts.  
        // Next we will look for imprinter points with composite atts.  These
        // may have resulted in a new point.  If there is a non composite att
        // point that doesn't have an ORIGINAL att we will know it is new
        // from the imprinter composite att point and we know to put a composite
        // att on it.
        found = false;
        for(j=coincident_pts_with_composite_att.size(); j>0 && !found; j--)
        {
          TBPoint *tmp_pt = coincident_pts_with_composite_att.get_and_step();
          list.clean_out();
          tmp_pt->get_simple_attribute("IMPRINTER",list);
          if(list.size() > 0)
            found = true;
        }
        if(found)
        {
          // Now put a composite att on any point that doesn't have one.
          for(j=coincident_pts_without_composite_att.size(); j>0; j--)
          {
            TBPoint *tmp_pt = coincident_pts_without_composite_att.get_and_step();
            list.clean_out();
            tmp_pt->get_simple_attribute("ORIGINAL", list);
            if(list.size() == 0)
            {
              // The point was not in the original model and therefore was created by 
              // the imprint of a pt with a composite att.  We need to put a composite
              // att on it.
              list.clean_out();
              coincident_pts_with_composite_att.get()->get_simple_attribute("COMPOSITE_GEOM",list);
              tmp_pt->append_simple_attribute_virt(list.get());
            }
          }
        }
      }
    }

    for(i=coincident_pts_with_composite_att.size(); i>0; i--)
    {
      TBPoint *pt = coincident_pts_with_composite_att.get_and_step();
      list.clean_out();
      pt->get_simple_attribute("IMPRINTER",list);
      if(list.size() > 0)
        pt->remove_simple_attribute_virt(list.get());
      list.clean_out();
      pt->get_simple_attribute("ORIGINAL",list);
      if(list.size() > 0)
        pt->remove_simple_attribute_virt(list.get());
    }
    for(i=coincident_pts_without_composite_att.size(); i>0; i--)
    {
      TBPoint *pt = coincident_pts_without_composite_att.get_and_step();
      list.clean_out();
      pt->get_simple_attribute("IMPRINTER",list);
      if(list.size() > 0)
        pt->remove_simple_attribute_virt(list.get());
      list.clean_out();
      pt->get_simple_attribute("ORIGINAL",list);
      if(list.size() > 0)
        pt->remove_simple_attribute_virt(list.get());
    }
  }
  delete pt_tree;

  CubitSense rel_sense;
  while(all_curves_with_composite_att.size())
  {
    DLIList<Curve*> other_crvs;
    DLIList<BodySM*> other_bodies;
    DLIList<double> other_distances;

    Curve *cur_crv = all_curves_with_composite_att.extract();
    crv_tree->remove(cur_crv);

    BodySM *cur_body = cur_crv->bodysm();
    DLIList<Curve*> coincident_crvs_with_composite_att, coincident_crvs_without_composite_att;
    DLIList<Curve*> close_crvs;
    CubitBox bbox = cur_crv->bounding_box();
    crv_tree->find(bbox, close_crvs);

    for(j=close_crvs.size(); j>0; j--)
    {
      Curve *other_crv = close_crvs.get_and_step();
      BodySM *other_body = other_crv->bodysm();
      // Only consider curves from other bodies.
      if(cur_body != other_body)
      {
        if(this->about_spatially_equal(cur_crv, other_crv, rel_sense, geom_factor, 0))
        {
          CubitVector pos1, pos2;
          double cur_dist;
          cur_crv->get_geometry_query_engine()->entity_entity_distance(
            cur_crv, other_crv, pos1, pos2, cur_dist );
          if(other_bodies.move_to(other_body))
          {
            int list_index = other_bodies.get_index();
            other_distances.reset();
            other_distances.step(list_index);
            double prev_dist = other_distances.get();
            if(cur_dist < prev_dist)
            {
              other_distances.change_to(cur_dist);
              other_crvs.reset();
              other_crvs.step(list_index);
              other_crvs.change_to(other_crv);
            }
          }
          else
          {
            other_bodies.append(other_body);
            other_crvs.append(other_crv);
            other_distances.append(cur_dist);
          }
        }
      }
      coincident_crvs_with_composite_att.append(cur_crv);
      for(j=other_crvs.size(); j>0; j--)
      {
        Curve *crv = other_crvs.get_and_step();
        list.clean_out();
        crv->get_simple_attribute("COMPOSITE_GEOM", list);
        if(list.size() > 0)
        {
          coincident_crvs_with_composite_att.append(other_crv);
          if(all_curves_with_composite_att.move_to(other_crv))
            all_curves_with_composite_att.extract();
        }
        else
          coincident_crvs_without_composite_att.append(other_crv);
      }
    }

    // If we have found at least one other crv coincident with the current crv...
    if(coincident_crvs_with_composite_att.size() > 1 ||
      coincident_crvs_without_composite_att.size() > 0)
    {
      // If there is at least one curve without a composite att that is an imprinter we
      // will remove all composite atts from coincident curves
      bool found = false;
      for(j=coincident_crvs_without_composite_att.size(); j>0 && !found; j--)
      {
        Curve *tmp_crv = coincident_crvs_without_composite_att.get_and_step();
        list.clean_out();
        tmp_crv->get_simple_attribute("IMPRINTER",list);
        if(list.size() > 0)
          found = true;
      }
      if(found)
      {
        // Remove all composite atts.
        for(j=coincident_crvs_with_composite_att.size(); j>0; j--)
        {
          Curve *tmp_crv = coincident_crvs_with_composite_att.get_and_step();
          list.clean_out();
          tmp_crv->get_simple_attribute("COMPOSITE_GEOM",list);
          if(list.size() > 0)
            tmp_crv->remove_simple_attribute_virt(list.get());
        }
      }
      else
      {
        // There were no imprinter crvs that didn't have composite atts.  
        // Next we will look for imprinter crvs with composite atts.  These
        // may have resulted in a new crv.  If there is a non composite att
        // crv that doesn't have an ORIGINAL att we will know it is new
        // from the imprinter composite att crv and we know to put a composite
        // att on it.
        found = false;
        for(j=coincident_crvs_with_composite_att.size(); j>0 && !found; j--)
        {
          Curve *tmp_crv = coincident_crvs_with_composite_att.get_and_step();
          list.clean_out();
          tmp_crv->get_simple_attribute("IMPRINTER",list);
          if(list.size() > 0)
            found = true;
        }
        if(found)
        {
          // Now put a composite att on any crv that doesn't have one.
          for(j=coincident_crvs_without_composite_att.size(); j>0; j--)
          {
            Curve *tmp_crv = coincident_crvs_without_composite_att.get_and_step();
            list.clean_out();
            tmp_crv->get_simple_attribute("ORIGINAL", list);
            if(list.size() == 0)
            {
              // The crv was not in the original model and therefore was created by 
              // the imprint of a crv with a composite att.  We need to put a composite
              // att on it.
              list.clean_out();
              coincident_crvs_with_composite_att.get()->get_simple_attribute("COMPOSITE_GEOM",list);
              tmp_crv->append_simple_attribute_virt(list.get());
            }
          }
        }
      }
    }

    for(i=coincident_crvs_with_composite_att.size(); i>0; i--)
    {
      Curve *crv = coincident_crvs_with_composite_att.get_and_step();
      list.clean_out();
      crv->get_simple_attribute("IMPRINTER",list);
      if(list.size() > 0)
        crv->remove_simple_attribute_virt(list.get());
      list.clean_out();
      crv->get_simple_attribute("ORIGINAL",list);
      if(list.size() > 0)
        crv->remove_simple_attribute_virt(list.get());
    }
    for(i=coincident_crvs_without_composite_att.size(); i>0; i--)
    {
      Curve *crv = coincident_crvs_without_composite_att.get_and_step();
      list.clean_out();
      crv->get_simple_attribute("IMPRINTER",list);
      if(list.size() > 0)
        crv->remove_simple_attribute_virt(list.get());
      list.clean_out();
      crv->get_simple_attribute("ORIGINAL",list);
      if(list.size() > 0)
        crv->remove_simple_attribute_virt(list.get());
    }
  }
  delete crv_tree;
  for(i=all_curves.size(); i>0; i--)
  {
    Curve *cur_curve = all_curves.get_and_step();
    list.clean_out();
    cur_curve->get_simple_attribute("IMPRINTER",list);
    if(list.size() > 0)
      cur_curve->remove_simple_attribute_virt(list.get());
    list.clean_out();
    cur_curve->get_simple_attribute("ORIGINAL",list);
    if(list.size() > 0)
      cur_curve->remove_simple_attribute_virt(list.get());
  }
  for(i=all_points.size(); i>0; i--)
  {
    TBPoint *cur_point = all_points.get_and_step();
    list.clean_out();
    cur_point->get_simple_attribute("IMPRINTER",list);
    if(list.size() > 0)
      cur_point->remove_simple_attribute_virt(list.get());
    list.clean_out();
    cur_point->get_simple_attribute("ORIGINAL",list);
    if(list.size() > 0)
      cur_point->remove_simple_attribute_virt(list.get());
  }
}

void CompositeEngine::process_curves_after_imprint(Curve *att_cur, 
                                                   Curve *other_cur,
                                                   DLIList<BodySM*> &new_sms)
{
  DLIList<CubitSimpleAttrib> list;

  if(att_cur == other_cur)
  {
    // This case will happen when we have manually added one of the existing
    // curves on the surface to be imprinted to the "new_ENTITIES" list.  We
    // do this in cases where the curve to imprint on the surface exactly 
    // falls on one of the existing curves.  In this case the ACIS face 
    // doesn't get new curves created but we need to consider the 
    // curve on the face as new because it may have been hidden in a 
    // composite surface and needs to be reintroduced.  So, in this
    // case we will remove the attribute so that the curve is no longer
    // hidden.
    att_cur->get_simple_attribute("COMPOSITE_GEOM",list);
    att_cur->remove_simple_attribute_virt(list.get());
  }
  else
  {
    other_cur->get_simple_attribute("IMPRINT_PREEXISTING",list);
    if(list.size() == 0)
    {
      // This is a new bridge created by imprinting this hidden bridge.  In 
      // this case we need to add a COMPOSITE_GEOM attribute to the new bridge
      // so we don't see a resulting imprinted bridge from the hidden bridge.
      list.clean_out();
      att_cur->get_simple_attribute("COMPOSITE_GEOM",list);
      other_cur->append_simple_attribute_virt(list.get());
    }
    else
    {
      // This bridge existed before the imprint operation.  Therefore it
      // could also have a COMPOSITE_GEOM attribute on it.  Check this.
      list.clean_out();
      other_cur->get_simple_attribute("COMPOSITE_GEOM",list);
      if(list.size() == 0)
      {
        // It doesn't have a COMPOSITE_GEOM attribute so we need to remove
        // the COMPOSITE_GEOM from att_bridge because the hidden nature gets
        // wiped out by the imprinting process.
        att_cur->get_simple_attribute("COMPOSITE_GEOM",list);
        att_cur->remove_simple_attribute_virt(list.get());
        TBOwner *bridge_owner = att_cur->owner();
        CompositeCurve *cc_bridge_owner = dynamic_cast<CompositeCurve*>(bridge_owner);
        if(cc_bridge_owner)
        {
          TBOwner *cc_owner = cc_bridge_owner->owner();
          HiddenEntitySet *hes = dynamic_cast<HiddenEntitySet*>(cc_owner);
          if(hes)
          {
            CompositeSurface *cs = dynamic_cast<CompositeSurface*>(hes->owner());
            if(cs)
              cs->HadBridgeRemoved = 1;
              // This is currently how we are notifying the owning CompositeSurface
              // that it needs to be deactivated and rebuilt.  It really has 
              // nothing to do with the bridge being removed though.  Bad.
          }
        }
      }
      else
      {
        // This bridge was also hidden so do nothing.
      }
    }
  }
}

void CompositeEngine::process_points_after_imprint(TBPoint *att_pt, 
                                                   TBPoint *other_pt,
                                                   DLIList<BodySM*> &new_sms)
{
  int i;
  DLIList<CubitSimpleAttrib> list;

  if(att_pt == other_pt)
  {
    // This case will happen when we have manually added one of the existing
    // pts on the surface to be imprinted to the "new_ENTITIES" list.  We
    // do this in cases where the pt to imprint on the surface exactly 
    // falls on one of the existing pts.  In this case the ACIS face 
    // doesn't get new pts created but we need to consider the 
    // pt on the face as new because it may have been hidden in a 
    // composite curve and needs to be reintroduced.  So, in this
    // case we will remove the attribute so that the pt is no longer
    // hidden.
    att_pt->get_simple_attribute("COMPOSITE_GEOM",list);
    att_pt->remove_simple_attribute_virt(list.get());
  }
  else
  {
    other_pt->get_simple_attribute("IMPRINT_PREEXISTING",list);
    if(list.size() == 0)
    {
      // This is a new bridge created by imprinting this hidden bridge.  In 
      // this case we need to add a COMPOSITE_GEOM attribute to the new bridge
      // if possible so we don't see a resulting imprinted bridge from the hidden bridge.
      int num_visible_curves = 0;
      DLIList<TopologyBridge*> curves;
      att_pt->get_parents_virt(curves);
      for(i=curves.size(); i--;)
      {
        list.clean_out();
        TopologyBridge *c = curves.get_and_step();
        c->get_simple_attribute("COMPOSITE_GEOM", list);
        if(list.size() == 0)
          num_visible_curves++;
      }
      if(num_visible_curves > 2)
      {
        list.clean_out();
        att_pt->get_simple_attribute("COMPOSITE_GEOM",list);
        att_pt->remove_simple_attribute_virt(list.get());
        TBOwner *bridge_owner = att_pt->owner();
        CompositePoint *cp_bridge_owner = dynamic_cast<CompositePoint*>(bridge_owner);
        if(cp_bridge_owner)
        {
          TBOwner *cp_owner = cp_bridge_owner->owner();
          HiddenEntitySet *hes = dynamic_cast<HiddenEntitySet*>(cp_owner);
          if(hes)
          {
            CompositeCurve *cc = dynamic_cast<CompositeCurve*>(hes->owner());
            if(cc)
              cc->HadBridgeRemoved = 1;
              // This is currently how we are notifying the owning CompositeCurve
              // that it needs to be deactivated and rebuilt.  It really has 
              // nothing to do with the bridge being removed though.  Bad.
          }
        }
      }
      else
      {
        list.clean_out();
        att_pt->get_simple_attribute("COMPOSITE_GEOM",list);
        other_pt->append_simple_attribute_virt(list.get());
      }
    }
    else
    {
      // This bridge existed before the imprint operation.  Therefore it
      // could also have a COMPOSITE_GEOM attribute on it.  Check this.
      list.clean_out();
      other_pt->get_simple_attribute("COMPOSITE_GEOM",list);
      if(list.size() == 0)
      {
        // It doesn't have a COMPOSITE_GEOM attribute so we need to remove
        // the COMPOSITE_GEOM from att_bridge because the hidden nature gets
        // wiped out by the imprinting process.
        att_pt->get_simple_attribute("COMPOSITE_GEOM",list);
        att_pt->remove_simple_attribute_virt(list.get());
        TBOwner *bridge_owner = att_pt->owner();
        CompositePoint *cp_bridge_owner = dynamic_cast<CompositePoint*>(bridge_owner);
        if(cp_bridge_owner)
        {
          TBOwner *cp_owner = cp_bridge_owner->owner();
          HiddenEntitySet *hes = dynamic_cast<HiddenEntitySet*>(cp_owner);
          if(hes)
          {
            CompositeCurve *cc = dynamic_cast<CompositeCurve*>(hes->owner());
            if(cc)
              cc->HadBridgeRemoved = 1;
              // This is currently how we are notifying the owning CompositeCurve
              // that it needs to be deactivated and rebuilt.  It really has 
              // nothing to do with the bridge being removed though.  Bad.
          }
        }
      }
      else
      {
        // This bridge was also hidden so do nothing.
      }
    }
  }
}

// This is a copy of the function in MergeTool with the difference that it
// accepts a layer flag to dictate at which level the topology is traversed
// (solid modeler level or virtual level).
CubitBoolean CompositeEngine::about_spatially_equal( Curve *curve_1, Curve *curve_2,
                                               CubitSense &relative_sense, 
                                               double tolerance_factor,
                                               int layer)
{
  if( curve_1 == curve_2 )
    return CUBIT_TRUE;

  relative_sense = CUBIT_FORWARD;
  const double ONE_THIRD = 1.0/3.0;

  // Find the point 1/3 along curve_1 
  CubitVector test_point_1, test_point_2;
  if( curve_1->position_from_fraction( ONE_THIRD, test_point_1 ) != CUBIT_SUCCESS )
    return CUBIT_FALSE;

  // See if the 1/3 point on curve_1 lies on curve_2
  if ( curve_2->closest_point_trimmed(test_point_1, test_point_2)
       != CUBIT_SUCCESS )
  {
    return CUBIT_FALSE;
  }

  if ( GeometryQueryTool::instance()->
       about_spatially_equal(test_point_1, test_point_2, tolerance_factor )
       != CUBIT_SUCCESS )
  {
    return CUBIT_FALSE;
  }
  
  CubitVector tangent_1, tangent_2;
  if( curve_1->closest_point(test_point_2, test_point_1, &tangent_1) != CUBIT_SUCCESS )
    return CUBIT_FALSE;
  if( curve_2->closest_point(test_point_1, test_point_2, &tangent_2) != CUBIT_SUCCESS )
    return CUBIT_FALSE;

  //If one of the curves is zero-length, it will have a zero
  //tangent vector.
  double len_product = tangent_1.length() * tangent_2.length();
  if( len_product > CUBIT_DBL_MIN )
  {
    double dot_product = (tangent_1 % tangent_2);
    if (dot_product < 0)
    relative_sense = CUBIT_REVERSED;
  }
  else
  {
    //If one of the tangents is zero-length, one of the curves had
    //better be as well.
    assert( (curve_1->measure() * curve_2->measure()) < CUBIT_RESABS );
  }

  //compare the start and end vertices to be spatially equal.
  DLIList<TBPoint*> curve_1_points(2), curve_2_points(2);
  DLIList<TopologyBridge*> c1pts, c2pts;
  curve_1->get_children(c1pts, false, layer);
  curve_2->get_children(c2pts, false, layer);
  CAST_LIST( c1pts, curve_1_points, TBPoint );
  CAST_LIST( c2pts, curve_2_points, TBPoint );

  if( curve_1->bridge_sense() == CUBIT_REVERSED )
    curve_1_points.reverse();
  if( curve_2->bridge_sense() == CUBIT_REVERSED )
    curve_2_points.reverse();

  TBPoint* curve_1_start = curve_1_points.get(); 
  curve_1_points.last();
  TBPoint* curve_1_end =  curve_1_points.get();

  TBPoint* curve_2_start = curve_2_points.get(); 
  curve_2_points.last();
  TBPoint* curve_2_end =  curve_2_points.get();

  if (relative_sense == CUBIT_REVERSED)
    std::swap(curve_2_start, curve_2_end);

  if (curve_1_start == curve_1_end ||curve_2_start == curve_2_end)
  {
    CubitVector c1start = curve_1_start->coordinates();
    CubitVector c2start = curve_2_start->coordinates();
    if ((curve_1_start != curve_1_end) ||
        (curve_2_start != curve_2_end) ||
        !GeometryQueryTool::instance()->about_spatially_equal(c1start, c2start, tolerance_factor))
      return CUBIT_FALSE;
  }
  else
  {
    CubitVector c1start = curve_1_start->coordinates();
    CubitVector c1end = curve_1_end->coordinates();
    CubitVector c2start = curve_2_start->coordinates();
    CubitVector c2end = curve_2_end->coordinates();
    if ((curve_1_start == curve_2_end) ||
        (curve_1_end == curve_2_start) ||
        !GeometryQueryTool::instance()->about_spatially_equal(c1start, c2start, tolerance_factor ) ||
        !GeometryQueryTool::instance()->about_spatially_equal(c1end, c2end, tolerance_factor ))
      return CUBIT_FALSE;
  }

  return CUBIT_TRUE;

}

// This function will try to determine if virtual topology bridges have
// been modified and if so will deactivate them so that they can be 
// rebuilt later using the COMPOSITE_GEOM attributes on the underlying
// solid model topology.
void CompositeEngine::remove_modified(DLIList<Surface*> &surfaces,  
                                      DLIList<Curve*> &curves, 
                                      DLIList<TBPoint*> &points)
{
  clean_out_deactivated_geometry();

  int i, j, k, m, n, w;
  int something_changed = 1;
  DLIList<TBPoint*> already_deactivated_points;
  DLIList<Curve*> already_deactivated_curves;
  DLIList<Surface*> already_deactivated_surfs;

  while(something_changed)
  {
    something_changed = 0;

    DLIList<TBPoint*> deactivated_points;
    DLIList<Curve*> deactivated_curves;
    DLIList<Surface*> deactivated_surfs;

    // Look for composite points that are out of date.
    for(w=points.size(); w--;)
    {
      CompositePoint *p = dynamic_cast<CompositePoint*>(points.get_and_step());
      if(p && !already_deactivated_points.is_in_list(p))
        deactivated_points.append(p);
    }
    deactivated_points.uniquify_ordered();

    // Look for composite curves that are out of date.
    for(w=curves.size(); w--;)
    {
      Curve *current_curve = curves.get_and_step();
      CompositeCurve *cur = dynamic_cast<CompositeCurve*>(current_curve);
      if(cur && !already_deactivated_curves.is_in_list(cur))
        deactivated_curves.append(cur);
    }
    deactivated_curves.uniquify_ordered();

    // Look for composite surfaces that are out of date.
    for(w=surfaces.size(); w--;)
    {
      CompositeSurface* csurf = dynamic_cast<CompositeSurface*> (surfaces.get_and_step());
      if (csurf && !already_deactivated_surfs.is_in_list(csurf))
        deactivated_surfs.append(csurf);
    }
    deactivated_surfs.uniquify_ordered();

    already_deactivated_points += deactivated_points;
    already_deactivated_curves += deactivated_curves;
    already_deactivated_surfs += deactivated_surfs;

    something_changed += deactivated_surfs.size() + deactivated_curves.size() +
      deactivated_points.size();

    // Now actually deactivate the out of date composite surfs.
    for(j=deactivated_surfs.size(); j--;)
    {
      CompositeSurface *csurf = dynamic_cast<CompositeSurface*>(deactivated_surfs.get_and_step());

      // We have to also deactivate the boundary curves.  When we deactivate
      // the CompositeSurface it removes all of the CompositeCoEdges associated
      // with it.  However, it doesn't deactivate the composite curves associated
      // with the composite coedges.  Therefore you can end up with a regular
      // CoEdge pointing to a CompositeCurve and if the CompositeCurve has more
      // than 1 curve in it later calls to replace_surface (which will in turn
      // call replace_curve) will fail.
      DLIList<Curve*> boundary_curves;
      csurf->curves(boundary_curves);
      for (k=boundary_curves.size(); k--; )
      {
        CompositeCurve* c = dynamic_cast<CompositeCurve*>(boundary_curves.get_and_step());
        assert(NULL != c);
        deactivated_curves.append_unique(c);
        already_deactivated_curves.append_unique(c);

        DLIList<TBPoint*> boundary_pts;
        c->points(boundary_pts);
        for (int e=boundary_pts.size(); e--; )
        {
          CompositePoint* p = dynamic_cast<CompositePoint*>(boundary_pts.get_and_step());
          deactivated_points.append_unique(p);
          already_deactivated_points.append_unique(p);
          notify_deactivated(p);
        }

        notify_deactivated(c);
      }

      notify_deactivated(csurf);

      DLIList<Curve*> hidden;
      csurf->get_hidden_curves(hidden);
      for (k=hidden.size(); k--; )
      {
        CompositeCurve* hcurve = dynamic_cast<CompositeCurve*>(hidden.pop());
        assert(NULL != hcurve);

        deactivated_curves.append_unique(hcurve);
        already_deactivated_curves.append_unique(hcurve);
        notify_deactivated(hcurve);

        if(hcurve->num_curves() == 1)
        {
          Curve *c = hcurve->get_curve(0);
          DLIList<TopologyBridge*> end_pts;
          c->get_children(end_pts, false, 0);
          for(m=end_pts.size(); m--;)
          {
            TBPoint *cur_p = dynamic_cast<TBPoint*>(end_pts.get_and_step());
            if(cur_p)
            {
              CompositePoint* cp = dynamic_cast<CompositePoint*>(cur_p->owner());
              if(cp)
                cur_p = (TBPoint*)cp;
              TBOwner *own = cur_p->owner();
              HiddenEntitySet *hes = dynamic_cast<HiddenEntitySet*>(own);
              if(hes)
              {
                CompositeCurve *cc = dynamic_cast<CompositeCurve*>(hes->owner());
                if(cc)
                {
                  deactivated_curves.append_unique(cc);
                  already_deactivated_curves.append_unique(cc);
                  notify_deactivated(cc);

                  DLIList<TBPoint*> hidden_pts;
                  cc->get_hidden_points(hidden_pts);
                  for (n=hidden_pts.size(); n--; )
                  {
                    CompositePoint *hpoint = dynamic_cast<CompositePoint*>(hidden_pts.pop());
                    assert(NULL != hpoint);
                    deactivated_points.append_unique(hpoint);
                    already_deactivated_points.append_unique(hpoint);
                    notify_deactivated(hpoint);
                  }
                }
              }
            }
          }
        }
      }
    }

    // Now actually deactivate the out of date composite curves.
    for(j=deactivated_curves.size(); j--;)
    {
      CompositeCurve *ccurve = dynamic_cast<CompositeCurve*>(deactivated_curves.get_and_step());

      DLIList<TBPoint*> boundary_pts;
      ccurve->points(boundary_pts);
      for (k=boundary_pts.size(); k--; )
      {
        CompositePoint* p = dynamic_cast<CompositePoint*>(boundary_pts.get_and_step());
        deactivated_points.append_unique(p);
        already_deactivated_points.append_unique(p);
        notify_deactivated(p);
      }

      notify_deactivated(ccurve);

      int j;
      DLIList<TBPoint*> hidden;
      ccurve->get_hidden_points(hidden);
      for (j=hidden.size(); j--; )
      {
        CompositePoint *hpoint = dynamic_cast<CompositePoint*>(hidden.pop());
        assert(NULL != hpoint);
        deactivated_points.append_unique(hpoint);
        already_deactivated_points.append_unique(hpoint);
        notify_deactivated(hpoint);
      }
    }

    // Now actually deactivate the out of date composite points.
    for(j=deactivated_points.size(); j--;)
    {
      CompositePoint* cpoint = dynamic_cast<CompositePoint*> (deactivated_points.pop());
      notify_deactivated(cpoint);
    }
  }

  int remove_point_atts = 1;
  if(remove_point_atts)
  {
    // Remove any COMPOSITE_GEOM attributes on points that
    // have a valence of more than two (real curves - hidden curves).
    for(i=points.size(); i>0; i--)
    {
      TBPoint *pt = points.get_and_step();
      CompositePoint *cp = dynamic_cast<CompositePoint*>(pt);
      if(cp)
        pt = cp->get_point();
      DLIList<CubitSimpleAttrib> attribs;
      pt->get_simple_attribute("COMPOSITE_GEOM", attribs);
      if(attribs.size() > 0)
      {
        DLIList<TopologyBridge*> tmp_curves;
        pt->get_parents_virt(tmp_curves);
        int num_curves = 0;
        for(j=tmp_curves.size(); j>0; j--)
        {
          TopologyBridge *crv = tmp_curves.get_and_step();
          DLIList<CubitSimpleAttrib> attribs;
          crv->get_simple_attribute("COMPOSITE_GEOM", attribs);
          if(attribs.size() == 0)
            num_curves++;
        }
        if(num_curves != 2)
        {
          for(j=attribs.size(); j>0; j--)
          {
            const CubitSimpleAttrib &csa = attribs.get_and_step();
            pt->remove_simple_attribute_virt(csa);
          }
        }
      }
    }
  }
}

//-------------------------------------------------------------------------
// Purpose       : Constructor
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 08/25/03
//-------------------------------------------------------------------------
CompositeEngine::CompositeEngine()
{
  CubitStatus result = GeometryQueryTool::instance()->
    register_intermediate_engine( this );
  assert(result == CUBIT_SUCCESS);
}

//-------------------------------------------------------------------------
// Purpose       : Get composite-level query results from
//                 query of underlying topology.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/13/02
//-------------------------------------------------------------------------
/*
void CompositeEngine::fix_up_query_results( DLIList<TopologyBridge*>& list,
                                            bool keep_hidden )
{
  DLIList<TopologyBridge*> tmp_list;
  int i;
  
  list.reset();
  for( i = list.size(); i--; )
  {
    TopologyBridge* bridge = list.get();
    TBOwnerSet* set = dynamic_cast<TBOwnerSet*>(bridge->owner());
    if( set )
    {
      tmp_list.clean_out();
      set->get_owners( tmp_list );
      if( tmp_list.size() )
      {
        tmp_list.reset();
        list.change_to( tmp_list.get_and_step() );
        for( int j = 1; j < tmp_list.size(); j++ )
          list.insert( tmp_list.get_and_step() );
      }
      else
      {
        list.change_to(0);
      }
    }
    list.step();
  }
  
  list.remove_all_with_value( 0 );
  
  list.reset();
  for( i = list.size(); i--; )
  {
    TopologyBridge* bridge = list.get();
    TopologyBridge *next = 0;
    while( next = dynamic_cast<TopologyBridge*>(bridge->owner()) )
      bridge = next;
    
    if( !keep_hidden && dynamic_cast<HiddenEntitySet*>(bridge->owner()) )
    {
      list.change_to(0);
    }
    else if( list.get() != bridge )
    {
      if( list.is_in_list( bridge ) )
        list.change_to(0);
      else
        list.change_to(bridge);
    }
    
    list.step();
  }
  
  list.remove_all_with_value( 0 );
}
*/


//-------------------------------------------------------------------------
// Purpose       : Combine two curves
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/13/02
//-------------------------------------------------------------------------
CompositeCurve* CompositeEngine::composite( Curve* keep_curve,
                                            Curve* dead_curve,
                                            TBPoint* keep_point,
                                            bool remove_partition )
{
  if( keep_curve == dead_curve )
  {
    PRINT_ERROR("Cannot remove vertex from single-vertex curve.\n");
    return 0;
  }
  
  CompositeCurve* result = 0;
  CompositeCurve* ckeep = dynamic_cast<CompositeCurve*>(keep_curve);
  CompositeCurve* cdead = dynamic_cast<CompositeCurve*>(dead_curve);
  bool replaced1 = false;
  bool replaced2 = false;
  
  if( !ckeep )
  {
    ckeep = replace_curve( keep_curve );
    replaced1 = true;
  }
  if( !cdead )
  {
    cdead = replace_curve( dead_curve );
    replaced2 = true;
  }
  
  CompositePoint* comppoint = dynamic_cast<CompositePoint*>(keep_point);
  if( keep_point && !comppoint )
  {
    comppoint = dynamic_cast<CompositePoint*>(keep_point->owner());
    assert( comppoint!= NULL );
  }
  
  if( !ckeep || !cdead || 
      !(result = combine(ckeep, cdead, comppoint, remove_partition)) )
  {
    if( replaced1 && ckeep )
    {
      Curve* s = remove_composite( ckeep );
      assert( s != 0 );
    }
    if( replaced2 && cdead )
    {
      Curve* s = remove_composite( cdead );
      assert( s != 0 );
    }
  }
  
  return result;
}

//-------------------------------------------------------------------------
// Purpose       : Replace a "point" curve with a composite
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 11/26/02
//-------------------------------------------------------------------------
CompositePoint* CompositeEngine::replace_point( TBPoint* point )
{
  assert( !dynamic_cast<CompositePoint*>(point) );
  return new CompositePoint( point );
}

//-------------------------------------------------------------------------
// Purpose       : Replace a "real" curve with a composite
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/13/02
//-------------------------------------------------------------------------
CompositeCurve* CompositeEngine::replace_curve( Curve* curve )
{
  DLIList<TopologyBridge*> points, coedges;
  //curve->get_children_virt( points );
  //fix_up_query_results(points);
  curve->get_children( points, true, COMPOSITE_LAYER );

  if(points.size() > 2)
    return 0;
    
  points.reset();
  TBPoint* start_pt = dynamic_cast<TBPoint*>(points.get_and_step());
  TBPoint*   end_pt = dynamic_cast<TBPoint*>(points.get_and_step());

  CompositeCurve* composite = dynamic_cast<CompositeCurve*>(curve);
  
  if(!composite)
    composite = new CompositeCurve( curve );
  else
  {
    if(composite->num_curves() > 1)
    {
      PRINT_ERROR("Error replacing existing composite curve with more than one underlying curve\n");
      return 0;
    }
  }

  CompositePoint* start = dynamic_cast<CompositePoint*>(start_pt);
  if( !start ) 
    start = replace_point( start_pt );
  CompositePoint* end = 0;
  if ( end_pt == start_pt )
    end = start;
  else if ( !(end = dynamic_cast<CompositePoint*>(end_pt)) )
    end = replace_point( end_pt );
  assert( start && end );

  composite->start_point( start );
  composite->  end_point(   end );

  DLIList<TopologyBridge*> existing_composite_coedges;
  if(dynamic_cast<CompositeCurve*>(curve))
  {
    dynamic_cast<CompositeCurve*>(curve)->get_curve(0)->get_parents_virt(coedges);
    curve->get_parents_virt(existing_composite_coedges);
  }
  else
  {
    curve->get_parents_virt( coedges );
  }
  for( int i = coedges.size(); i--; )
  {
    CoEdgeSM* coedge = dynamic_cast<CoEdgeSM*>(coedges.get_and_step());
    assert(coedge);
    CompositeCoEdge* ccoedge = NULL;
    for(int h=existing_composite_coedges.size(); h>0; h--)
    {
      CompositeCoEdge *temp = dynamic_cast<CompositeCoEdge*>(existing_composite_coedges.get_and_step());
      if(temp->get_coedge(0) == coedge)
      {
        ccoedge = temp;
        h = 0;
      }
    }
    if(!ccoedge)
    {
      ccoedge = new CompositeCoEdge( coedge );
      if( composite->get_sense(0) == CUBIT_REVERSED )
        ccoedge->reverse();
      
      assert( ccoedge->get_curve() == 0 );
      composite->add( ccoedge );
    }
  }
  
  return composite;
}


//-------------------------------------------------------------------------
// Purpose       : Replace a "real" surface with a composite
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/13/02
//-------------------------------------------------------------------------
CompositeSurface* CompositeEngine::replace_surface( Surface* surface )
{
  DLIList<TopologyBridge*> loops, coedges, curves;
  int i, j;
  if( dynamic_cast<CompositeSurface*>(surface) )
    return 0;
  
  CompositeSurface* compsurf = new CompositeSurface( surface );
  //surface->get_children_virt( loops );
  //fix_up_query_results(loops);
  surface->get_children( loops, false, COMPOSITE_LAYER );
  loops.reset();
  for( i = loops.size();i--; )
  {
    TopologyBridge* loop_bridge = loops.get_and_step();
    CompositeLoop* comploop = new CompositeLoop();
    compsurf->add( comploop );
    if( loop_bridge->owner() )
    {
      loop_bridge->owner()->swap_bridge( loop_bridge, comploop, false );
      loop_bridge->owner(0);
    }
    //compsurf->hidden_entities().hide( loop_bridge );
    
    coedges.clean_out();
    //loop_bridge->get_children_virt( coedges );
    //fix_up_query_results(coedges);
    loop_bridge->get_children( coedges, false, COMPOSITE_LAYER );
    
    coedges.reset();
    CompositeCoEdge* prev = 0;
    for( j = coedges.size(); j--; )
    {
      CoEdgeSM* coedge = dynamic_cast<CoEdgeSM*>(coedges.get_and_step());
      CompositeCoEdge* comp_coedge = dynamic_cast<CompositeCoEdge*>(coedge);
      if( !comp_coedge )
      {
        curves.clean_out();
        //coedge->get_children_virt( curves );
        //fix_up_query_results(curves);
        coedge->get_children(curves, false, COMPOSITE_LAYER);
        assert( curves.size() == 1 );
        replace_curve( dynamic_cast<Curve*>(curves.get()) );
        comp_coedge = dynamic_cast<CompositeCoEdge*>(coedge->owner());
        assert( comp_coedge!= NULL );
      }
      else
      {
        /*
        PRINT_INFO("\nStart: %lf %lf %lf", comp_coedge->start_point()->coordinates().x(),
          comp_coedge->start_point()->coordinates().y(), comp_coedge->start_point()->coordinates().z());
        PRINT_INFO("\nEnd: %lf %lf %lf", comp_coedge->end_point()->coordinates().x(),
          comp_coedge->end_point()->coordinates().y(), comp_coedge->end_point()->coordinates().z());
       */
      }
      comploop->insert_after( comp_coedge, prev );
      prev = comp_coedge;
    }
  }
  
  return compsurf;
}
      


//-------------------------------------------------------------------------
// Purpose       : Replace a "real" lump with a composite
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/19/02
//-------------------------------------------------------------------------
CompositeLump* CompositeEngine::replace_lump( Lump* lump )
{
  DLIList<TopologyBridge*> shells, surfaces;
  int i, j;
  if( dynamic_cast<CompositeLump*>(lump) )
    return 0;
  
  CompositeLump* complump = new CompositeLump( lump );
  //lump->get_children_virt( shells );
  //fix_up_query_results(shells);
  lump->get_children(shells, false, COMPOSITE_LAYER);
  shells.reset();
  for( i = shells.size();i--; )
  {
    CompositeShell* compshell = new CompositeShell();
    complump->add( compshell );
    TopologyBridge* shell_bridge = shells.get_and_step();
    if( shell_bridge->owner() )
    {
      shell_bridge->owner()->swap_bridge(shell_bridge, compshell, false);
      shell_bridge->owner(0);
    }
    //complump->hidden_entities().hide( shell_bridge );
    
    
    surfaces.clean_out();
    //shell_bridge->get_children_virt( surfaces );
    //fix_up_query_results(surfaces);
    shell_bridge->get_children(surfaces, false, COMPOSITE_LAYER);
    
    surfaces.reset();
    for( j = surfaces.size(); j--; )
    {
      Surface* surface = dynamic_cast<Surface*>(surfaces.get_and_step());
      CompositeSurface* compsurf = dynamic_cast<CompositeSurface*>(surface);
      if( !compsurf )
      {
        compsurf = replace_surface( surface );
        assert( compsurf!= NULL );
      }

      compshell->add( compsurf, CUBIT_FORWARD );
    }
  }
  
  return complump;
}
      


//-------------------------------------------------------------------------
// Purpose       : Replace a "real" bodysm with a composite
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/19/02
//-------------------------------------------------------------------------
CompositeBody* CompositeEngine::replace_body( BodySM* body )
{
  DLIList<TopologyBridge*> lumps;
  int i;
  if( dynamic_cast<CompositeBody*>(body) )
    return 0;
  
  CompositeBody* compbody = new CompositeBody( );
  compbody->add(body);
  body->get_children( lumps, false, COMPOSITE_LAYER );
  lumps.reset();
  for( i = lumps.size();i--; )
  {
    Lump* lump = dynamic_cast<Lump*>(lumps.get_and_step());
    CompositeLump* complump = dynamic_cast<CompositeLump*>(lump);
    if( !complump )
      complump = replace_lump( lump );
    compbody->add( complump );
  }
  
  return compbody;
}

TBPoint* CompositeEngine::remove_composite( CompositePoint* composite )
{
  assert( composite->next_curve() == 0 );
  assert( !dynamic_cast<HiddenEntitySet*>(composite->owner()) );
  TBPoint* result = composite->get_point();
  delete composite;
  return result;
}



//-------------------------------------------------------------------------
// Purpose       : Reomve a 1-curve composite
//
// Special Notes : Inverse of replace_curve(..)
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/13/02
//-------------------------------------------------------------------------
Curve* CompositeEngine::remove_composite( CompositeCurve* composite )
{
  assert( composite->num_curves() == 1 && 
          !composite->has_parent_composite_surface() &&
          !dynamic_cast<HiddenEntitySet*>(composite->owner()) &&
          !composite->is_stitched() );

  CompositeCoEdge* coedge = 0;
  while( (coedge = composite->first_coedge()) )
  {
    assert( !coedge->get_loop() && coedge->num_coedges() == 1 );
    composite->remove( coedge );
    CoEdgeSM* real_coedge = coedge->get_coedge(0);
    coedge->remove_coedge( 0 );
    if( coedge->owner() )
      coedge->owner()->swap_bridge( coedge, real_coedge, false );
    delete coedge;
  }
  
  CompositePoint* sp = composite->start_point();
  CompositePoint* ep = composite->end_point();
  Curve* curve = composite->get_curve( 0 );
  bool reversed = composite->get_sense(0) == CUBIT_REVERSED;
  composite->remove_curve(0);
  if( composite->owner() )
    composite->owner()->swap_bridge( composite, curve, reversed );
 
  composite->start_point(0);
  composite->end_point(0);
  delete composite;
  
  if( ! sp->next_curve() )
    remove_composite(sp);
  if( ep != sp && !ep->next_curve() )
    remove_composite(ep);
  
  // we must notify the graphics of the modify from "real" to virtual -- KGM
  // I realize that this is not where the other notifies are completed but there
  // is no knowledge of the change later on. 
  CubitObservable* observable = dynamic_cast<CubitObservable*>(curve->topology_entity());
  if (observable)
  {
      AppUtil::instance()->send_event(observable, GEOMETRY_MODIFIED);
  }

  return curve;
}


//-------------------------------------------------------------------------
// Purpose       : Remove a 1-surface composite
//
// Special Notes : inverse of replace_surface(..)
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/13/02
//-------------------------------------------------------------------------
Surface* CompositeEngine::remove_composite( CompositeSurface* composite )
{
  assert( ! composite->has_hidden_entities() );
  assert( composite->next_co_surface(0) == 0 );
  CompositeLoop* loop;
  CompositeCoEdge* coedge;
  CompositeCurve* curve;
  
  while( (loop = composite->first_loop()) )
  {
    composite->remove( loop );
    while( (coedge = loop->first_coedge()) )
    {
      loop->remove(coedge);
      curve = coedge->get_curve();
      if( curve->num_curves() == 1 &&
          !curve->has_parent_composite_surface() &&
          !curve->is_stitched()) 
        remove_composite(curve);
    }
    delete loop;
  }
  
  Surface* surface = composite->get_surface( 0 );
  bool reversed = composite->get_sense(0) == CUBIT_REVERSED;
  composite->remove_surface(0);
  if( composite->owner() )
    composite->owner()->swap_bridge( composite, surface, reversed );
  delete composite;
  
  return surface;
}

//-------------------------------------------------------------------------
// Purpose       : Remove a 1-lump composite
//
// Special Notes : inverse of replace_lump(..)
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/19/02
//-------------------------------------------------------------------------
Lump* CompositeEngine::remove_composite( CompositeLump* composite )
{
  assert( composite->num_lumps() == 1 );
  assert( composite->get_body() == 0 );
  
  for( CompositeShell* shell = composite->first_shell();
       shell != 0;
       shell = composite->next_shell( shell ) )
  {
    composite->remove( shell );
    while( shell->first_co_surf() )
    {
      CompositeCoSurf* cos = shell->first_co_surf();
      CompositeSurface* surface = cos->get_surface();
      shell->remove( cos );
      surface->remove( cos );
      delete cos;
      
      if( surface->next_co_surface() == 0 &&
          !surface->has_hidden_entities() )
        remove_composite( surface );
    }
    delete shell;
  }
  
  Lump* result = composite->get_lump( 0 );
  composite->remove_bridge(result);
  if( composite->owner() )
    composite->owner()->swap_bridge( composite, result, false );
  delete composite;
  
  return result;
}


//-------------------------------------------------------------------------
// Purpose       : Remove a 1-body composite
//
// Special Notes : inverse of replace_body(..)
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/19/02
//-------------------------------------------------------------------------
BodySM* CompositeEngine::remove_composite( CompositeBody* composite )
{
  assert( composite->num_bodies() == 1 );
  
  while( CompositeLump* lump = composite->next_lump() )
  {
    composite->remove( lump );
    remove_composite( lump );
  }
  
  BodySM* result = composite->get_body( 0 );
  composite->remove_bridge(result);
  if( composite->owner() )
    composite->owner()->swap_bridge( composite, result, false );
  delete composite;
  
  return result;
}

  
//-------------------------------------------------------------------------
// Purpose       : Combine two curves into a composite
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/13/02
//-------------------------------------------------------------------------
CompositeCurve* CompositeEngine::combine( CompositeCurve* keep, 
                                          CompositeCurve* dead,
                                          CompositePoint* keep_point,
                                          bool remove_partitions )
{
  if( keep->start_point() == keep->end_point() ||
      dead->start_point() == dead->end_point() )
    return 0;
  
    // find the point to remove by compositing.
  CompositePoint* common = keep->common_point( dead );
  if (!common)
    return 0;
    
  if (keep_point && common == keep_point)
  {
    common = keep->other_point(common);
    if (dead->other_point(common) != keep_point)
      return 0;
  }
  
  bool prepend = (common == keep->start_point());
  bool reverse = ((common == dead->end_point()) != prepend);
  
    // Order coedges such that each coedge from the 
    // dead curve is grouped with the corresponding
    // coedge from the curve to be kept.
  DLIList<CompositeCoEdge*> dead_coedges, sorted_coedges;
  CompositeCoEdge* coedge;
  for( coedge = dead->first_coedge();
       coedge != 0;
       coedge = dead->next_coedge( coedge ) )
    dead_coedges.append( coedge );
  
  int keep_count = 0;
  CompositeCoEdge *keep_coedge = NULL;
  for( keep_coedge = keep->first_coedge();
       keep_coedge != 0;
       keep_coedge = keep->next_coedge( keep_coedge ) )
  {
    keep_count++;
    LoopSM* keep_loop = keep_coedge->get_parent_loop();
    CubitSense keep_sense = keep_coedge->sense();
    
    for( int j = dead_coedges.size(); j--;  )
    {
      CompositeCoEdge* dead_coedge = dead_coedges.step_and_get();
      if( dead_coedge->get_parent_loop() == keep_loop &&
         (reverse == (keep_sense != dead_coedge->sense()) ) )
      {
        dead_coedges.extract();
        sorted_coedges.append( dead_coedge );
        break;
      }
    }
  }
  
    // If not all coedges were paired up, can't proceed.
  if( sorted_coedges.size() != keep_count )
    return 0;
  
    // If dead is reversed wrt keep, reverse it such that
    // both curves have the same relative sense.
  if( reverse )
  {
    dead->reverse();
    for( coedge = dead->first_coedge();
         coedge != 0;
         coedge = dead->next_coedge( coedge) )
      coedge->reverse();
  }
  
    // Combine the CoEdges
  sorted_coedges.reset();
  for ( keep_coedge = keep->first_coedge();
        keep_coedge != 0;
        keep_coedge = keep->next_coedge( keep_coedge ) )
  {
    CompositeCoEdge* dead_coedge = sorted_coedges.get_and_step();
    keep_coedge->combine( dead_coedge, prepend );
    
    if ( dead_coedge->owner() )
      dead_coedge->owner()->remove_bridge( dead_coedge );
    delete dead_coedge;
  }
  
  
    // hide point
  keep->hidden_entities().hide( common );
  
    // combine curves
  keep->combine( dead, prepend );
  if (prepend)
    keep->start_point( dead->other_point( common ) );
  else
    keep->end_point( dead->other_point( common ) );
  
  if (dead->owner())
    dead->owner()->remove_bridge(dead);
  dead->start_point(0);
  dead->end_point(0);
  delete dead;
  
  if( remove_partitions )
    remove_partition_point( common );
  
  return keep;
}


//-------------------------------------------------------------------------
// Purpose       : Restore hidden point 
//
// Special Notes : May create point-curve if point is hidden by a surface.
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/27/04
//-------------------------------------------------------------------------
CubitStatus CompositeEngine::restore_point( TBPoint* point )
{
  HiddenEntitySet* set;
  
  CompositePoint* comp = dynamic_cast<CompositePoint*>(point);
  if (!comp)
    comp = dynamic_cast<CompositePoint*>(point->owner());
  if (!comp)
    return CUBIT_FAILURE;
  
  set = dynamic_cast<HiddenEntitySet*>(comp->owner());
  if (!set)
    return CUBIT_FAILURE;
  
  if (dynamic_cast<CompositeCurve*>(set->owner()))
  {
    point = comp->get_point();
    if (!restore_point_in_curve(comp))
      return CUBIT_FAILURE;
    
      // check if we still have a composite point
      // if this was the only point splitting the curve,
      // then the composite point will have been destroyed.
    comp = dynamic_cast<CompositePoint*>(point->owner());
    if (!comp)
      return CUBIT_SUCCESS;
    
      // If the point is no longer hidden, we're done.
    set = dynamic_cast<HiddenEntitySet*>(comp->owner());
    if (!set)
      return CUBIT_SUCCESS;
  }
  
    // Need to create point-curve?
  if (dynamic_cast<CompositeSurface*>(set->owner()))
  {
    if (!restore_point_in_surface(comp))
      return CUBIT_FAILURE;
    
      // If the point is no longer hidden, we're done.
    set = dynamic_cast<HiddenEntitySet*>(comp->owner());
    if (!set)
      return CUBIT_SUCCESS;
  }
  
    // Something went wrong.
  return CUBIT_FAILURE;
}
    


//-------------------------------------------------------------------------
// Purpose       : Restore a point removed when a curve was composited
//                 (split a composite point)
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/13/02
//-------------------------------------------------------------------------
CubitStatus CompositeEngine::restore_point_in_curve( TBPoint* point )
{
    // find CompositePoint
  CompositePoint* comppoint = dynamic_cast<CompositePoint*>(point);
  if( !comppoint )
    comppoint = dynamic_cast<CompositePoint*>(point->owner());
  if( !comppoint )
    return CUBIT_FAILURE;

  // find CompositeCurve.
  HiddenEntitySet* hs = dynamic_cast<HiddenEntitySet*>(comppoint->owner());
  if( !hs )
    return CUBIT_FAILURE;
  CompositeCurve* compcurve = dynamic_cast<CompositeCurve*>(hs->owner());
  if( !compcurve )
    return CUBIT_FAILURE;
  
    // find which curves in the composite contain the point
  DLIList<TopologyBridge*> points;
  int index = -1;
  for( int i = 0; i < compcurve->num_curves() && index < 0; i++ )
  {
    Curve* curve = compcurve->get_curve(i);
    points.clean_out();
    curve->get_children( points, true, COMPOSITE_LAYER );
    for( int j = points.size(); j--; )
    {
      TBPoint* pt = dynamic_cast<TBPoint*>(points.get_and_step());
      assert( pt!= NULL );
      if( pt == comppoint )
      {
      	index = i;
      	break;
      }
    }
  }
  
  assert( index >= 0 );
  Curve *curve1, *curve2;
  return split( compcurve, index, curve1, curve2 );
}


  
  
  

//-------------------------------------------------------------------------
// Purpose       : Split a composite at the index-th point (after the
//                 index-th underlying curve)
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/13/02
//-------------------------------------------------------------------------
CubitStatus CompositeEngine::split( CompositeCurve* curve, int index,
                                    Curve*& result1, Curve*& result2 )
{
  result1 = result2 = 0;
  
    // Split the composite geometry
  CompositeCurve* new_curve = curve->split( curve->get_curve(index) );
  if( !new_curve )
    return CUBIT_FAILURE;
  HiddenEntitySet* ownerSet = dynamic_cast<HiddenEntitySet*>(curve->owner());
  if( ownerSet )
    ownerSet->hide(new_curve);
  
    // Split owning CoEdges
  for( CompositeCoEdge* old_coe = curve->first_coedge();
       old_coe != 0;
       old_coe = curve->next_coedge( old_coe ) )
  {
    CompositeCoEdge* new_coe = old_coe->split( index );
    assert( new_coe != NULL );
    new_curve->add( new_coe );
    if( ownerSet )
      ownerSet->hide(new_coe);
    
    if( old_coe->get_loop() )
    	old_coe->get_loop()->insert_after( new_coe, old_coe );
  }
  
    // update end point of original curve
  DLIList<TopologyBridge*> children;
  int last = curve->num_curves() - 1;
  Curve* last_curve = curve->get_curve(last);
  children.clean_out();
  last_curve->get_children( children, true, COMPOSITE_LAYER );
  assert( children.size() == 2 );
  children.reset();
  if( curve->get_sense(last) == CUBIT_FORWARD )
    children.step();
  CompositePoint* endpt = dynamic_cast<CompositePoint*>(children.get());
  assert( !!endpt );
  CompositePoint* oldendpt = curve->end_point();
  curve->end_point( endpt );
  
    // attach start and end point to new curve
  children.clean_out();
  Curve* first_curve = new_curve->get_curve(0);
  first_curve->get_children( children, true, COMPOSITE_LAYER );
  assert( children.size() == 2 );
  children.reset();
  if( new_curve->get_sense(0) == CUBIT_REVERSED )
    children.step();
  CompositePoint* startpt = dynamic_cast<CompositePoint*>(children.get());
  assert( startpt!= NULL );
  new_curve->start_point( startpt );
  new_curve->end_point( oldendpt );

    // unhide restored (split) point
  HiddenEntitySet& old_hidden = curve->hidden_entities();
  old_hidden.restore( curve->end_point() );
  old_hidden.restore( new_curve->start_point() );
  if( ownerSet )
  {
    ownerSet->hide(curve->end_point());
    if(new_curve->start_point() != curve->end_point() )
      ownerSet->hide(new_curve->start_point());
  }
  
    // move other hidden points from old curve to new curve
    // as needed
  for( int j = 0; j < new_curve->num_curves(); j++ )
  {
    Curve* r_curve = new_curve->get_curve(j);
    children.clean_out();
    r_curve->get_children( children, true, COMPOSITE_LAYER );
    TBPoint* sp = dynamic_cast<TBPoint*>(children.get_and_step() );
    TBPoint* ep = dynamic_cast<TBPoint*>(children.get_and_step() );
    assert( sp && ep && sp != ep );
    if( old_hidden.restore( sp ) )
      new_curve->hidden_entities().hide( sp );
    if( old_hidden.restore( ep ) )
      new_curve->hidden_entities().hide( ep );
  }

    // Remove composite geometry if a 'composite' of one curve
  result1 = curve;
  result2 = new_curve;
  if( curve->num_curves() == 1 && !ownerSet &&
      !curve->has_parent_composite_surface() )
    result1 = remove_composite( curve );
  if( new_curve->num_curves() == 1 && !ownerSet && 
      !new_curve->has_parent_composite_surface() )
    result2 = remove_composite( new_curve );
    
  return CUBIT_SUCCESS;
}


//-------------------------------------------------------------------------
// Purpose       : Create a composite surface by removing a curve
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/13/02
//-------------------------------------------------------------------------
CompositeSurface* CompositeEngine::remove_curve( Curve* dead_curve,
                                                 bool remove_partitions,
                                                 Surface* survivor )
{
  DLIList<TopologyBridge*> temp, coedges, loops, surfaces, curves,
                           surf1_shells, surf2_shells;
  
  if( dynamic_cast<CompositeCurve*>(dead_curve->owner()) )
    dead_curve = dynamic_cast<Curve*>(dead_curve->owner());
  
    // If this is a null-geometry composite curve (a composite
    // point-curve without a corresponding real point-curve),
    // just destroy it and hide the defining point.
  CompositeCurve* compcurve = dynamic_cast<CompositeCurve*>(dead_curve);
  if (compcurve && compcurve->num_curves() == 0)
  {
    CompositeSurface* surf = compcurve->first_coedge()->get_loop()->get_surface();
    CompositePoint* geom_pt = destroy_point_curve(compcurve);
    surf->hidden_entities().hide(geom_pt);
    return surf;
  }
  
    // Get two coedges to remove, and put curve(s) in 'curves'
  dead_curve->get_parents_virt( coedges );
  if( coedges.size() != 2 )  
    return 0;
  
  coedges.reset();
  TopologyBridge* coedge1 = coedges.get_and_step();
  TopologyBridge* coedge2 = coedges.get_and_step();
  
    // Get surfaces corresponding to loops
  loops.clean_out();
  coedge1->get_parents_virt( loops );
  assert( loops.size() == 1 );
  surfaces.clean_out();
  loops.get()->get_parents_virt( surfaces );
  assert( surfaces.size() == 1 );
  Surface* surface1 = dynamic_cast<Surface*>(surfaces.get());
  
  loops.clean_out();
  coedge2->get_parents_virt( loops );
  assert( loops.size() == 1 );
  surfaces.clean_out();
  loops.get()->get_parents_virt( surfaces );
  assert( surfaces.size() == 1 );
  Surface* surface2 = dynamic_cast<Surface*>(surfaces.get());
  
  if (surface2 == survivor)
    std::swap(surface1,surface2);
  
    // Make sure surfaces have same parent shells
  surface1->get_parents_virt( surf1_shells );
  surface2->get_parents_virt( surf2_shells );
  if( surf1_shells != surf2_shells )
    return 0;
  
    // get/make composites of surfaces
  bool created_surf1 = false, created_surf2 = false;
  CompositeSurface* compsurf1 = dynamic_cast<CompositeSurface*>(surface1);
  if( ! compsurf1 )
  {
    compsurf1 = replace_surface( surface1 );
    created_surf1 = true;
  }
  
  CompositeSurface* compsurf2 = 0;
  if( surface2 == surface1 )
  {
    compsurf2 = compsurf1;
  }
  else
  {
    compsurf2 = dynamic_cast<CompositeSurface*>(surface2);
    if( !compsurf2 )
    {
      compsurf2 = replace_surface( surface2 );
      created_surf2 = true;
    }
  }
  assert( compsurf1 && compsurf2 );
  
    // get composite(s) of curve(s)
  compcurve = dynamic_cast<CompositeCurve*>(dead_curve);
  if( !compcurve) 
    compcurve = dynamic_cast<CompositeCurve*>(dead_curve->owner());
  assert( !!compcurve );
  
    // Get composite coedges
  CompositeCoEdge *compcoedge1 = dynamic_cast<CompositeCoEdge*>(coedge1);
  if (!compcoedge1)
    compcoedge1 = dynamic_cast<CompositeCoEdge*>(coedge1->owner());
  CompositeCoEdge *compcoedge2 = dynamic_cast<CompositeCoEdge*>(coedge2);
  if (!compcoedge2)
    compcoedge2 = dynamic_cast<CompositeCoEdge*>(coedge2->owner());
  
    // combine the composite surfaces (not loops yet.)
  if( compsurf1 != compsurf2 )
  {
      // If second surface has opposite sense, fix it so
      // that both have the same sense.
    if( compcoedge2->sense() == compcoedge1->sense() )
    {
      bool b_reverse_coedges = true; // automatically handle reversing the coedges for the loop
      compsurf2->reverse_sense();
      for( CompositeLoop* loop = compsurf2->first_loop();
           loop != 0;
           loop = compsurf2->next_loop( loop ) )
        loop->reverse(b_reverse_coedges);
    }
    
      // Move all loops from old surface to new
    while( CompositeLoop* loop = compsurf2->first_loop() )
    {
      compsurf2->remove( loop );
      compsurf1->add( loop );
    }
    
      // Combine the set of underlying, real surfaces
    compsurf1->combine( compsurf2 );
    if( compsurf2->owner() )
      compsurf2->owner()->remove_bridge( compsurf2 );
    delete compsurf2;
    compsurf2 = 0;
  }
  
  CompositeLoop* loop1 = compcoedge1->get_loop();
  CompositeLoop* loop2 = compcoedge2->get_loop();
    // same loop => remove sipe || split loop
  if( loop1 == loop2 ) 
  {
      // split loop
    if( loop1->next_coedge( compcoedge1 ) != compcoedge2 &&
        loop1->prev_coedge( compcoedge1 ) != compcoedge2 )
    {
      loop2 = new CompositeLoop();
      CompositeCoEdge* coedge = loop1->next_coedge( compcoedge1 );
      CompositeCoEdge* prev = 0;
      while( coedge != compcoedge2 )
      {
      	loop1->remove( coedge );
      	loop2->insert_after( coedge, prev );
      	prev = coedge;
      	coedge = loop1->next_coedge( compcoedge1 );
      }
      compsurf1->add( loop2 );
    }
      
    loop1->remove( compcoedge1 );
    loop1->remove( compcoedge2 );
    if( loop1->num_coedges() == 0 )
    {
      compsurf1->remove( loop1 );
      delete loop1;
    }
  }
    // stitch loops
  else 
  {
    CompositeCoEdge* coedge;
    
      // insert coedges
    while( loop2->num_coedges() > 1 )  // all but the dead one
    {
      coedge = loop2->next_coedge( compcoedge2 );
      loop2->remove( coedge );
      loop1->insert_before( coedge, compcoedge1 );
    }
    loop1->remove( compcoedge1 );
    loop2->remove( compcoedge2 );
    assert( loop2->num_coedges() == 0 );
    if( loop2->get_surface() )
      loop2->get_surface()->remove( loop2 );
    delete loop2;
    
      // If loops only had one coedge each, then we
      // just removed a hole.
    if( loop1->num_coedges() == 0 )
    {
      if( loop1->get_surface() )
        loop1->get_surface()->remove( loop1 );
      delete loop1;
    }
  }
  
  compsurf1->hidden_entities().hide( compcoedge1 );
  compsurf1->hidden_entities().hide( compcoedge2 );
  
    // clean up dead curve(s) and points
  curves.clean_out();
  compcurve->start_point()->get_parents_virt( curves );
  if( curves.size() == 1 )
    compsurf1->hidden_entities().hide( compcurve->start_point() );
  curves.clean_out();
  compcurve->end_point()->get_parents_virt( curves );
  if( curves.size() == 1 )
    compsurf1->hidden_entities().hide( compcurve->end_point() );
  compsurf1->hidden_entities().hide( compcurve );
  
  if( remove_partitions )
    remove_partition_curves( compcurve );
  
  return compsurf1;
}


static void cme_hide_surface( HiddenEntitySet& set,
                              CompositeSurface* surface )
{
  set.hide( surface );
  CompositeLoop* loop = 0;
  while (NULL != (loop = surface->next_loop( loop )))
  {
    set.hide( loop );
    CompositeCoEdge* coedge = loop->first_coedge();
    do {
      set.hide( coedge );
      
      bool hide_curve = true;
      CompositeCurve* curve = coedge->get_curve();
      CompositeCoEdge* curve_coedge = 0;
      while ((curve_coedge = curve->next_coedge( curve_coedge )))
        if (curve_coedge->owner() != &set)
          { hide_curve = false; break; }
      
      if (hide_curve)
      {
        set.hide( curve );
        int num_pts = 1 + (curve->start_point() == curve->end_point());
        for (int i = 0; i < num_pts; i++)
        {
          CompositePoint* pnt = i ? curve->end_point() : curve->start_point();
          bool hide_pnt = true;
          CompositeCurve* pnt_curve = 0;
          while ((pnt_curve = pnt->next_curve( pnt_curve )))
            if (pnt_curve->owner() != &set)
              { hide_pnt = false; break; }
          
          if (hide_pnt)
            set.hide( pnt );
        }
      }
    } while((coedge = loop->next_coedge(coedge)) != loop->first_coedge());
  }
}

static void cme_unhide_surface( CompositeSurface* surf )
{
  HiddenEntitySet* set = dynamic_cast<HiddenEntitySet*>(surf->owner());
  assert( !!set );
  set->restore( surf );
  CompositeLoop* loop = 0;
  while (NULL != (loop = surf->next_loop( loop )))
  {
    set->restore( loop );
    CompositeCoEdge* coedge = loop->first_coedge();
    do {
      if (coedge->owner() == set)
      {
        CompositeCurve* curve = coedge->get_curve();
        set->restore( coedge );
        set->restore( curve );
        if (curve->start_point()->owner() == set)
          set->restore( curve->start_point() );
        if (curve->end_point()->owner() == set)
          set->restore( curve->end_point() );
      }
    } while ((coedge = loop->next_coedge(coedge)) != loop->first_coedge());
  }
}
  


//-------------------------------------------------------------------------
// Purpose       : Create composite volume by removing a surface
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/10/04
//-------------------------------------------------------------------------
CompositeLump* CompositeEngine::remove_surface( Surface* dead_surf,
                                                Surface* stitch_partner,
                                                bool /*remove_partitions*/  )
{
  bool okay = true;
  Surface* sptr;
  if ((sptr = dynamic_cast<Surface*>(dead_surf->owner())))
    dead_surf = sptr;
  if ((sptr = dynamic_cast<Surface*>(stitch_partner->owner())))
    stitch_partner = sptr;
  
    // Get composite surfaces
  bool replaced_1 = false, replaced_2 = false;
  CompositeSurface *surf1 = 0, *surf2 = 0;
  if (NULL == (surf1 = dynamic_cast<CompositeSurface*>(dead_surf)))
  {
    surf1 = replace_surface( dead_surf );
    replaced_1 = true;
  }
  if (NULL == stitch_partner)
  {
    surf2 = surf1;
  }
  else if (NULL == (surf2 = dynamic_cast<CompositeSurface*>(stitch_partner)))
  {
    surf2 = replace_surface( stitch_partner );
    replaced_2 = true;
  }  
  
    // Check if surfaces are stitched
  if (surf2 != surf1)
  {
    std::set<CompositeLoop*> used_loops;
    CompositeLoop* loop = 0;
    while (NULL != (loop = surf1->next_loop(loop)))
    {
      CompositeCoEdge* coedge1 = loop->first_coedge();
      CompositeCurve* curve = coedge1->get_curve();
      CompositeCoEdge* coedge2 = 0;
      while (NULL != (coedge2 = curve->next_coedge(coedge2)))
        if (coedge2->get_loop() &&
            coedge2->get_loop()->get_surface() == surf2)
          break;
      
      if (!coedge2) 
      {
        break;
        okay = false;
      }
      
      CompositeLoop* loop2 = coedge2->get_loop();
      if (used_loops.insert( loop2 ).second == false)
      {
        okay = false;
        break;
      }
      
      CompositeCoEdge* coedge = coedge1;
      while ((coedge = loop->next_coedge(coedge)) != coedge1)
      {
        curve = coedge->get_curve();
        coedge2 = 0;
        while (NULL != (coedge2 = curve->next_coedge(coedge2)))
          if (coedge2->get_loop() == loop2)
            break;
        
        if (!coedge2)
        {
          okay = false;
          break;
        }
      }
    }
    surf1->stitch( surf2 );
    
    if (!okay)
    {
      if (replaced_1)
        remove_composite( surf1 );
      if (replaced_2)
        remove_composite( surf2 );
      return 0;
    }
  }
  
    // Replace volumes/bodies with composites, combine bodies if necessary
  bool combined_bodies = false;
  bool replaced_body = false;
  DLIList<BodySM*> body_list;
  surf1->bodysms( body_list );
  if (body_list.size() != 1)
    return 0;
  
  BodySM* body1 = body_list.pop();
  BodySM* body2 = 0;
  CompositeBody* body = 0;
  
  if (surf2)
  {
    surf2->bodysms( body_list );
    if (body_list.size() != 1)
      return 0;
    body2 = body_list.pop();
  }
  
  if (body2 && body2 != body1)
  {
    body = combine_bodies( body1, body2 );
    combined_bodies = true;
  }
  else
  {
    if (NULL == (body = dynamic_cast<CompositeBody*>(body1)) &&
        NULL == (body = dynamic_cast<CompositeBody*>(body1->owner())))
    {
      body = replace_body( body1 );
      replaced_body = true;
    }
  }
  
  if (!body)
    return 0;
  
    // Get CoSurfaces to combine must be exactly two of them.
  CompositeCoSurf *cosurf1 = 0, *cosurf2 = 0;
  cosurf1 = surf1->next_co_surface(NULL);
  cosurf2 = surf1->next_co_surface(cosurf1);
  if (!cosurf1)
  {
    okay = false;
  }
  else if (cosurf2)
  {
    if (surf2 || surf1->next_co_surface(cosurf2))
      okay = false;
  }
  else if (surf2)
  {
    cosurf2 = surf2->next_co_surface(NULL);
    if (!cosurf2 || surf2->next_co_surface(cosurf2))
      okay = false;
  }
  if (!okay)
  {
//    if (replaced_body)
//      restore_body( body );
    return 0;
  }
  
  
    // combine the composite lumps (not shells yet.)
  CompositeLump* lump1 = cosurf1->get_shell()->get_lump();
  CompositeLump* lump2 = cosurf2->get_shell()->get_lump();
  if( lump1 != lump2 )
  {
    
      // Move all shells from old lump to new
    while( CompositeShell* shell = lump2->first_shell() )
    {
      lump2->remove( shell );
      lump1->add( shell );
    }
    
      // Combine the set of underlying, real lumps
    lump1->combine( lump2 );
    if( lump2->owner() )
      lump2->owner()->remove_bridge( lump2 );
    body->remove(lump2);
    delete lump2;
    lump2 = 0;
  }
  
    // remove surface(s)
  CompositeShell* shell1 = cosurf1->get_shell();
  CompositeShell* shell2 = cosurf2->get_shell();
  cme_hide_surface( lump1->hidden_entities(), surf1 );
  if (surf2)
    cme_hide_surface( lump1->hidden_entities(), surf2 );
  delete cosurf1;
  delete cosurf2;

  if (shell1 != shell2)
  {
    while (CompositeCoSurf* cosurf = shell2->next_co_surf(0))
    {
      shell2->remove( cosurf );
      shell1->add( cosurf );
    }
    lump1->remove( shell2 );
    delete shell2;
  }
  else // if (shell1 == shell2)
  {
    shell2 = split_shell( shell1 );
    if (shell2)
      lump1->add( shell2 );
  }
  
  return lump1;
}


//-------------------------------------------------------------------------
// Purpose       : Test if a CompositeShell needs to be split, and
//                 split it
//
// Special Notes : Copied from PartitionEngine::split_shell
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/11/04
//-------------------------------------------------------------------------
CompositeShell* CompositeEngine::split_shell( CompositeShell* shell_to_split )
{
  std::map<CompositeSurface*, int> marks;

    // Make sure all cosurface surface marks are cleared
  CompositeCoSurf* cosurf = 0;
  while( (cosurf = shell_to_split->next_co_surf( cosurf )) )
    marks[cosurf->get_surface()] = 0;
  
    // Identify non-manifold surfaces, marking them either
    // with a 2 if they can be used to split the volume
    // (if they are part of a connected patch for which the
    // bounadary of that patch intersects the volume boundary
    // at all curves) or a 3 if they are other non-manifold
    // surfaces.  This will get a bit tricky if there are
    // non-manifold surfaces hanging off of the patch of
    // split surfaces.
  
    // First for all non-manifold surfaces, if the surface has
    // a curve that is not shared with any other surface, mark
    // it with a 3.  Otherwise mark it with a 2.
  cosurf = 0;
  DLIList<CompositeSurface*> surf_stack;
  while ( (cosurf = shell_to_split->next_co_surf(cosurf)) )
  {
    CompositeSurface* surf = cosurf->get_surface();
      // If we haven't done this surface yet and it is non-manifold
    if ( !marks[surf] && surf->find_next(cosurf) )
    {
      marks[surf] = 2;
      bool no_free_curve = true;
      CompositeLoop* loop = 0;
      while( no_free_curve && (loop = surf->next_loop(loop)) ) 
      {
        CompositeCoEdge* coedge = loop->first_coedge();
        do 
        {
          CompositeCurve* curve = coedge->get_curve();
            // If the curve has more than one coedge, it 
            // is not a free curve (this also accounts for
            // the case where the curve is a non-manifold 
            // curve on the surface interioir -- e.g. a sipe)
          if ( !curve->next_coedge(coedge) && 
               curve->next_coedge(0) == coedge )
          {
            no_free_curve = false;
            break;
          }
          coedge = loop->next_coedge(coedge);
        } while( coedge != loop->first_coedge() );
      }
      
      if( !no_free_curve )
      {
        marks[surf] = 3;
        surf_stack.append( surf );
      }
    }
  }
  
    // Now for each surface we marked with a three, traverse
    // and mark adjacent surfaces until we come to a curve
    // connected to more that two surfaces.
  while( surf_stack.size() ) 
  {
    CompositeSurface* surf = surf_stack.pop();
    CompositeLoop* loop = 0;
    while ( (loop = surf->next_loop(loop)) )
    {
      CompositeCoEdge* coedge = loop->first_coedge();
      do 
      {
        CompositeCurve* curve = coedge->get_curve();
        int split_count = 0;
        int boundary_count = 0;
        CompositeCoEdge* curve_coe = 0;
        while ( (curve_coe = curve->next_coedge(curve_coe) ) != NULL )
        {
          CompositeSurface* curve_surf = curve_coe->get_loop()->get_surface();
          switch ( marks[curve_surf] ) 
          {
            case 0: boundary_count++; break;
            case 2: split_count++;    break;
          }
        }
        
        if ( split_count == 1 && !boundary_count )
        {
          curve_coe = 0;
          while ( (curve_coe = curve->next_coedge(curve_coe) ) != NULL )
          {
            CompositeSurface* curve_surf = curve_coe->get_loop()->get_surface();
            if ( marks[curve_surf] == 2 )
            {
              marks[curve_surf] = 3;
              surf_stack.append(curve_surf);
            }
          }
        }
        
        coedge = loop->next_coedge( coedge );
      } while( coedge != loop->first_coedge() );
    }
  }
  
    // Now build a new shell by traversing cofaces, marking
    // each with that will go in a new shell with a 1.
  
    // Start with any cosurf that does not have a free
    // non-manifold surface (marked with a 3).  We'll handle
    // free non-manifold surfaces later.
  std::set<CompositeCoSurf*> marked_cosurfs;
  DLIList<CompositeCoSurf*> cosurf_stack;
  cosurf = 0;
  while ( (cosurf = shell_to_split->next_co_surf(cosurf)) )
    if ( marks[cosurf->get_surface()] != 3 )
      break;
  if ( cosurf )
  {
    marked_cosurfs.insert(cosurf);
    cosurf_stack.append( cosurf );
  }
  
    // Traverse over adjacent cosurfaces, marking them with a 1
  while (cosurf_stack.size())
  {
    cosurf = cosurf_stack.pop();
    CompositeSurface* surf = cosurf->get_surface();
    CompositeLoop* loop = 0;
    while ( (loop = surf->next_loop(loop)) )
    {
      CompositeCoEdge* coedge = loop->first_coedge();
      do
      {
        CompositeCurve* curve = coedge->get_curve();
        CompositeCoEdge* curve_coe = 0;
        CompositeCoSurf *boundary_cosurf = 0, *split_cosurf = 0;
        int split_cosurf_count = 0;
        while ( (curve_coe = curve->next_coedge(curve_coe)) )
        {
          if ( curve_coe == coedge )
            continue;
          
          bool same_coe_sense = curve_coe->sense() == coedge->sense();
          CompositeSurface* curve_surf = curve_coe->get_loop()->get_surface();
          CompositeCoSurf* curve_cosurf = 0;
          while ( (curve_cosurf = curve_surf->next_co_surface(curve_cosurf)) )
          {
            if ( curve_cosurf->get_shell() != shell_to_split )
              continue;
            
            bool same_cos_sense = curve_cosurf->sense() == cosurf->sense();
            if ( same_cos_sense == same_coe_sense )
              continue;
            
              // Always choose split surface first if we
              // found one
            if ( marks[curve_cosurf->get_surface()] == 2 ) {
              split_cosurf_count++;
              split_cosurf = curve_cosurf;
            }
            
              // Skip other non-manifold surfaces.  We'll
              // handle those later.
            else if( marks[curve_cosurf->get_surface()] != 3 )
              boundary_cosurf = curve_cosurf;
          }
        }
        
        CompositeCoSurf* next_cosurf = split_cosurf ? split_cosurf : boundary_cosurf;
        if ( marked_cosurfs.find(next_cosurf) == marked_cosurfs.end()
             && split_cosurf_count < 2 )
        {
          marked_cosurfs.insert(next_cosurf);
          cosurf_stack.append(next_cosurf);
        }
      
        coedge = loop->next_coedge(coedge);
      } while( coedge != loop->first_coedge() );
    } // end while (loop)
  } // end while (cosurf_stack.size())
    
  
    // build lists of cosurfaces, one for each shell and
    // one of other non-manifold surfaces
  DLIList<CompositeCoSurf*> marked_list, unmarked_list, other_list;
  while( (cosurf = shell_to_split->next_co_surf(0)) )
  {
    shell_to_split->remove(cosurf);
    if ( marks[cosurf->get_surface()] == 3 )
      other_list.append( cosurf );
    else if( marked_cosurfs.find(cosurf) != marked_cosurfs.end() )
      marked_list.append( cosurf );
    else
      unmarked_list.append( cosurf );
  }
  
    // If one of marked_list or unmarked_list is empty,
    // we can't split the shell yet.  Put cofaces back in
    // shell and exit.
  if ( !marked_list.size() || !unmarked_list.size() )
  {
    marked_list += unmarked_list;
    marked_list += other_list;
    marked_list.reverse();
    while ( marked_list.size() )
    {
      cosurf = marked_list.pop();
      shell_to_split->add( cosurf );
    }
    return 0;
  }
  
    // Put unmarked list back in old shell
  unmarked_list.reverse();
  while ( unmarked_list.size() )
  {
    cosurf = unmarked_list.pop();
    shell_to_split->add(cosurf);
  }
  
    // Put marked list in new shell
  CompositeShell* new_shell = new CompositeShell;
  marked_list.reverse();
  while ( marked_list.size() )
  {
    cosurf = marked_list.pop();
    new_shell->add(cosurf);
  }
  
    // Now sort out other non-manifold surfaces
  
    // Clear marks and get list of surface from cosurfaces
  surf_stack.clean_out();
  while( other_list.size() )
  {
    cosurf = other_list.pop();
    CompositeSurface* surf = cosurf->get_surface();
    if ( marks[surf] )
      surf_stack.append(surf);
  }
  
  insert_nonmanifold_surfaces( surf_stack, shell_to_split, new_shell );
  return new_shell;
}
 
//-------------------------------------------------------------------------
// Purpose       : After a shell is split, determine which of the two
//                 resulting shells each non-manifold surface belongs in.
//
// Special Notes : Copied from PartitionEngine::insert_nonmanifold_surfaces
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/11/04
//-------------------------------------------------------------------------
void CompositeEngine::insert_nonmanifold_surfaces( 
                                   DLIList<CompositeSurface*>& surf_stack,
                                   CompositeShell* shell1,
                                   CompositeShell* shell2 )
{
  DLIList<CompositeSurface*> known_list(surf_stack.size()), 
                           unknown_list(surf_stack.size());
  
  CompositeCoSurf* cosurf;
  CompositeSurface* surf;
  CompositeLoop* loop;
  CompositeCoEdge *coedge, *curve_coe;
  CompositeCurve* curve;
  
    // Loop until we've placed all the surfaces in one
    // shell or the other.
  while ( surf_stack.size() )
  {
    bool did_some = false;
    
      // Put any surfaces for which we immediately
      // know the shell into the appropriate shell.
      // Put others in known_list or unknown_list
      // depending on if we can determine which shell
      // they go in using a geometric comparison.
    known_list.clean_out();
    unknown_list.clean_out();

      // Take all surfaces out of stack in this loop.
      // We might put some back in after the loop (thus
      // the outer loop.)
    while ( surf_stack.size() )
    {
      surf = surf_stack.pop();
      
      CompositeShell* known_shell = 0;
      bool found_shell = false;
      loop = 0;
      while ( (loop = surf->next_loop(loop)) )
      {
        coedge = loop->first_coedge();
        do
        {
          curve = coedge->get_curve();
          curve_coe = 0;
          while ( (curve_coe = curve->next_coedge(curve_coe)) )
          {
            CompositeSurface* surf = curve_coe->get_loop()->get_surface();
            cosurf = 0;
            while ( (cosurf = surf->next_co_surface( cosurf )) )
            {
              if( cosurf->get_shell() == shell1 ||
                  cosurf->get_shell() == shell2 )
              {
                found_shell = true;
                if ( known_shell && known_shell != cosurf->get_shell() )
                  known_shell = 0;
                else
                  known_shell = cosurf->get_shell();
              }
            }
          }
        
          coedge = loop->next_coedge(coedge);
        } while( coedge != loop->first_coedge() );
      } // end while(loop)
        
        // This surface does not intersect the shell at
        // any curve.  We can't do this one yet.
      if ( !found_shell )
      {
        unknown_list.append( surf );
        continue;
      }
      
        // This surface intersected both shells at some
        // curves, but did not have a curve that intersected
        // only one shell.  We can do this one geometricly
        // if we have to.
      if ( !known_shell )
      {
        known_list.append( surf );
        continue;
      }
      
        // If we got this far, then the surface had at least
        // one curve that intersected only one of the shells.
        // We know it goes in that shell.
      did_some = true;
      CompositeCoSurf* cosurf1 = surf->find_first((CompositeShell*)0);
      CompositeCoSurf* cosurf2 = surf->find_next(cosurf1);
      known_shell->add(cosurf1);
      known_shell->add(cosurf2);
    
    } // while(surf_stack)  -- the inside one
    
      // Unknown_list always goes back in surf_stack to
      // try again.
    surf_stack += unknown_list;
    
      // If we did some surfaces, then put the rest back
      // in surf_stack and try again.  If they intersect
      // one of the surfaces we did place in this iteration,
      // we can avoid needing to do geometric checks.
    if ( did_some )
    {
      surf_stack += known_list;
      continue;
    }
    
      // If known_list is empty, somethings wrong (we're
      // going to loop forever.)  Abort the loop and try
      // to recover as best we can.
    if( !known_list.size() )
      break;
    
      // choose a single surface in do a geometric comparison
      // for, and put the rest back in surf_stack
    surf = known_list.pop();
    surf_stack += known_list;
    
    bool in_shell = false;
    if ( ! inside_shell( shell2, surf, in_shell ) )
    {
        // if inside_shell failed, abort.
      surf_stack.append(surf);
      break;
    }
    
    CompositeShell* shell = in_shell ? shell2 : shell1;
    shell->add(surf->find_first((CompositeShell*)0));
    shell->add(surf->find_first((CompositeShell*)0));
  }
  
    // something went wrong
  if( surf_stack.size() )
  {
    PRINT_ERROR("Internal error splitting volume at %s:%d\n"
                "Topology may be invalid.  Please report this.\n",
                __FILE__, __LINE__);
    while( surf_stack.size() ) 
    {
      CompositeSurface* surf = surf_stack.pop();
      cosurf = 0;
      while( (cosurf = surf->next_co_surface(cosurf)) )
        if( !cosurf->get_shell() )
          shell1->add(cosurf);
    }
  }
}


CubitStatus CompositeEngine::inside_shell( CompositeShell* const shell,
                                           CompositeSurface* const surf,
                                           bool& result )
{
    // Find the curve and coedge at which the nonmanifold surface
    // intersects the shells
  CompositeLoop* loop = 0;
  CompositeCoEdge *nonman_coedge = 0;
  while ( !nonman_coedge && (loop = surf->next_loop(loop)) )
  {
    CompositeCoEdge* loop_coedge = loop->first_coedge();
    do 
    {
        // Check if this curve is the curve of intersection
        // Iterate through curve coedges.
      CompositeCurve* curve = loop_coedge->get_curve();
      CompositeCoEdge* curve_coedge = 0;
      while ( (curve_coedge = curve->next_coedge(curve_coedge)) )
      {
        
        CompositeSurface* coedge_surf = curve_coedge->get_loop()->get_surface();
        if( coedge_surf->find_first(shell) )
        {
          nonman_coedge = curve_coedge;
          break;
        }
      }
    
      loop_coedge = loop->next_coedge(loop_coedge);
    } while( !nonman_coedge && loop_coedge != loop->first_coedge() );
  }
  
  if ( !nonman_coedge ) // bad input!
    return CUBIT_FAILURE;

    // There must exist two surfaces in the shell that are manifold 
    // in the shell and that are adjacent to the curve
  CompositeCurve* common_curve = nonman_coedge->get_curve();
  CompositeCoSurf *cosurf1 = 0, *cosurf2 = 0, *cosurf;
  CompositeCoEdge *coedge1 = 0, *coedge2 = 0, *coedge = 0;
  while ( (coedge = common_curve->next_coedge(coedge)) )
  {
    CompositeSurface* surf = coedge->get_loop()->get_surface();
    if ( (cosurf = surf->find_first(shell)) && !surf->find_next(cosurf) )
    {
      if( coedge1 ) {
        coedge2 = coedge;
        cosurf2 = cosurf;
      } else {
        coedge1 = coedge;
        cosurf1 = cosurf;
      }
    }
  }
  if ( !coedge1 || !coedge2 )
    return CUBIT_FAILURE;
  
    // Evaluate normals at midpoint of curve
  CubitVector base, tangent, point;
  double u = (common_curve->start_param()+common_curve->end_param())/2.0;
  common_curve->position_from_u( u, base );
  common_curve->closest_point( base, point, &tangent );
  tangent.normalize();
  
  CubitVector normal1, normal2, normal;
  surf->closest_point( base, 0, &normal );
  cosurf1->get_surface()->closest_point( base, 0, &normal1 );
  cosurf2->get_surface()->closest_point( base, 0, &normal2 );
  
    // Try to handle tangencies
  bool fix1 = (normal1 * normal).length_squared() < CUBIT_RESABS;
  bool fix2 = (normal2 * normal).length_squared() < CUBIT_RESABS;
  if ( fix1 || fix2 )
  {
    CubitVector dir = tangent * normal;
    double len = dir.length();
    assert(len > GEOMETRY_RESABS);
    dir /= len;
    if ( nonman_coedge->sense() == CUBIT_FORWARD )
      dir = -dir;
      
    CubitVector diag = surf->bounding_box().diagonal();
    if ( fix1 )
    {
      CubitVector d = cosurf1->get_surface()->bounding_box().diagonal();
      if ( diag.x() < d.x() ) diag.x(d.x());
      if ( diag.y() < d.y() ) diag.y(d.y());
      if ( diag.z() < d.z() ) diag.z(d.z());
    }
    if( fix2 )
    {
      CubitVector d = cosurf2->get_surface()->bounding_box().diagonal();
      if ( diag.x() < d.x() ) diag.x(d.x());
      if ( diag.y() < d.y() ) diag.y(d.y());
      if ( diag.z() < d.z() ) diag.z(d.z());
    }
    
    double step = 1e-3 * fabs( diag % dir );
    if ( step < 2*GEOMETRY_RESABS )
      step = 2*GEOMETRY_RESABS;
    
    for ( int i = 0; i < 1000; i++ )
    {
      point = base + i * step * dir;
      surf->closest_point( base, 0, &normal );
      if( fix1 )
        cosurf1->get_surface()->closest_point( base, 0, &normal1 );
      if( fix2 )
        cosurf2->get_surface()->closest_point( base, 0, &normal2 );
      
      bool done1 = !fix1 || (normal1 * normal).length_squared() > CUBIT_RESABS;
      bool done2 = !fix2 || (normal2 * normal).length_squared() > CUBIT_RESABS;
      if ( done1 && done2 )
      {
        fix1 = fix2 = false;
        break;
      }
    }
  }
  
  if ( fix1 || fix2 )
  {
    PRINT_ERROR("Failed to adjust for surface tangencies.\n"
                "This is a BUG.  %s:%d\n", __FILE__, __LINE__ );
    return CUBIT_FAILURE;
  }
  
  if ( nonman_coedge->sense() == CUBIT_FORWARD )
    normal = -normal;
  if ( cosurf1->sense() == coedge1->sense() )
    normal1 = -normal1;
  if ( cosurf2->sense() == coedge2->sense() )
    normal2 = -normal2;
  
  result = tangent.vector_angle(normal1, normal ) <=
           tangent.vector_angle(normal1, normal2);
  return CUBIT_SUCCESS;
}

  



//-------------------------------------------------------------------------
// Purpose       : Destroy a point-curve
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/28/04
//-------------------------------------------------------------------------
CompositePoint* CompositeEngine::destroy_point_curve( CompositeCurve* curve )
{
  CompositePoint* geom_pt = curve->start_point();
  CompositeCoEdge* coedge = curve->first_coedge();
  CompositeLoop* loop = coedge->get_loop();
  CompositeSurface* surf = loop->get_surface();

    // Make sure topology looks like a point curve
  if (curve->num_curves()        >  0       ||
      curve->end_point()         != geom_pt ||
      curve->next_coedge(coedge) != NULL    ||
      loop->next_coedge(coedge)  != coedge   )
  {
    assert( geom_pt == curve->start_point() );
    assert( curve->next_coedge(coedge) == NULL );
    assert( loop->next_coedge(coedge) == coedge );
    return 0;
  }

    // Destroy point-curve and corresponding topology
  curve->remove(coedge);
  loop->remove(coedge);
  delete coedge;

  surf->remove(loop);
  delete loop;

  curve->start_point(0);
  curve->end_point(0);
  delete curve;

  return geom_pt;
}
  


CubitStatus CompositeEngine::restore_curve( Curve* curve )
{
  CompositeSurface* new_surf = 0;
  
    // find compostie surface owning the curve
  TopologyBridge* bridge = curve;
  HiddenEntitySet* owner_set = NULL;
  while( bridge && !(owner_set = dynamic_cast<HiddenEntitySet*>(bridge->owner())) )
  {
    bridge = dynamic_cast<TopologyBridge*>(bridge->owner());
  }
  
  if( !owner_set )
    return CUBIT_FAILURE;
  
  CompositeSurface* compsurf = dynamic_cast<CompositeSurface*>(owner_set->owner());
  if( !compsurf )
    return CUBIT_FAILURE;
  
  CompositeCurve* compcurve = dynamic_cast<CompositeCurve*>(curve);
  if( !compcurve )
    compcurve = dynamic_cast<CompositeCurve*>(curve->owner());
  if( !compcurve )
    return CUBIT_FAILURE;
  
  CompositePoint *start_point = compcurve->start_point();
  CompositePoint *end_point   = compcurve->end_point();
  CompositeCurve *stitch_partner = 0; // find_stitch( compcurve );
  
    // Clean up any null-geometry point-curves attached to the
    // start and/or end points.
  CompositeCurve* itor = start_point->next_curve(NULL);
  while (itor)
  {
    CompositeCurve* next = start_point->next_curve(itor);
    if (itor->num_curves() == 0)
      destroy_point_curve(itor);
    itor = next;
  }
  itor = end_point->next_curve(NULL);
  while (itor)
  {
    CompositeCurve* next = end_point->next_curve(itor);
    if (itor->num_curves() == 0)
      destroy_point_curve(itor);
    itor = next;
  }
  
    // Check if end points are hidden by a CompositeCurve, and
    // if so, split the composite curve so that the end point(s)
    // are no longer hidden.
  HiddenEntitySet* hs;
  if( (hs = dynamic_cast<HiddenEntitySet*>(start_point->owner())) 
      && dynamic_cast<CompositeCurve*>(hs->owner())
      && !restore_point_in_curve( start_point ) )
    return CUBIT_FAILURE;
  if( end_point != start_point
      && (hs = dynamic_cast<HiddenEntitySet*>(end_point->owner()))
      && dynamic_cast<CompositeCurve*>(hs->owner())
      && !restore_point_in_curve( end_point ) )
    return CUBIT_FAILURE;
  
    // Find which loop(s) to insert CoEdge in, and where in
    // the loop(s) to insert it.  For the start point of the
    // curve, start_loop is the loop containing that point and
    // start_prev_coedge and start_next_coedge are the previous
    // and next coedges in the loop at that point, respectively.
    // Similarly, end_loop, end_prev_coedge and end_next_coedge
    // for the end point.
  CompositeLoop *start_loop = 0, *end_loop = 0;
  CompositeCoEdge *start_prev_coedge = 0, *start_next_coedge = 0;
  CompositeCoEdge *end_prev_coedge = 0, *end_next_coedge = 0;
  
  if ( ! find_coedges( compsurf, compcurve, start_point, 
                       start_prev_coedge, start_next_coedge ) )
    return CUBIT_FAILURE;
  if( start_prev_coedge )
  {
    start_loop = start_prev_coedge->get_loop();
    assert(start_next_coedge && 
           start_next_coedge->get_loop() == start_loop );
  }

  if ( ! find_coedges( compsurf, compcurve, end_point, 
                       end_prev_coedge, end_next_coedge ) )
    return CUBIT_FAILURE;
  if( end_prev_coedge )
  {
    end_loop = end_prev_coedge->get_loop();
    assert(end_next_coedge &&
           end_next_coedge->get_loop()== end_loop );
  }
                                     
    // must be all or none
  assert( !start_loop || (start_prev_coedge && start_next_coedge) );
  assert( !end_loop || (end_prev_coedge && end_next_coedge) );
    // closed curve?
  assert( (start_point != end_point) || 
          (start_prev_coedge == end_prev_coedge && start_next_coedge == end_next_coedge) );
  
    // Find coedges, and un-hide coedegs, curve, and
    // end points.  
  CompositeCoEdge* coedge1 = compcurve->first_coedge();
  while (coedge1 && coedge1->owner() != &compsurf->hidden_entities())
    coedge1 = compcurve->next_coedge(coedge1);
  assert( coedge1 );
  CompositeCoEdge* coedge2 = compcurve->next_coedge( coedge1 );
  while (coedge2 && coedge2->owner() != &compsurf->hidden_entities())
    coedge2 = compcurve->next_coedge(coedge2);
  if (!coedge2)
  {
    DLIList<CompositeCurve*> stitched;
    compcurve->get_stitched( stitched );
    stitched.remove( compcurve );
    while (!coedge2 && stitched.size())
    {
      CompositeCurve* other = stitched.pop();
      coedge2 = other->first_coedge();
      while (coedge2 && coedge2->owner() != &compsurf->hidden_entities())
        coedge2 = other->next_coedge(coedge2);
    }
  }
  assert (coedge2 || compcurve->geometry_type() == POINT_CURVE_TYPE);
  compsurf->hidden_entities().restore( coedge1 );
  if (coedge2) // no coedge2 for point-curves
    compsurf->hidden_entities().restore( coedge2 );

  compsurf->hidden_entities().restore( compcurve );
  if( compcurve->start_point()->owner() == &(compsurf->hidden_entities()) )
    compsurf->hidden_entities().restore( compcurve->start_point() );
  if( compcurve->end_point()->owner() == &(compsurf->hidden_entities()) )
    compsurf->hidden_entities().restore( compcurve->end_point() );
  if( stitch_partner )
    compsurf->hidden_entities().restore( stitch_partner );

    // If neither point intersected a loop (topologically), then
    // create a new loop containing the curve.
  if( !start_loop && !end_loop )
  {
      // hole
    if ( compcurve->start_point() == compcurve->end_point() )
    {
      start_loop = new CompositeLoop();
      start_loop->insert_after(coedge1,0);
      compsurf->add(start_loop);
      if (coedge2)  // no coedge2 for point-curves
      {
        end_loop = new CompositeLoop();
        end_loop->insert_after(coedge2,0);
        compsurf->add(end_loop);
        if ( CompLoopTool::loop_angle_metric(coedge1) >
             CompLoopTool::loop_angle_metric(coedge2) )
        {
          compsurf->remove(end_loop);
          new_surf = split_surface(compsurf, start_loop, end_loop);
        }
        else
        {
          compsurf->remove(start_loop);
          new_surf = split_surface(compsurf, end_loop, start_loop);
        }
      }
    }
      // hardline
    else
    {
      CompositeLoop* new_loop = new CompositeLoop();
      new_loop->insert_after( coedge1, 0 );
      new_loop->insert_after( coedge2, coedge1 );
      compsurf->add( new_loop );
    }
  }

    // If only one of the end points intersected a loop, then
    // create a sipe in that loop.
  else if( !start_loop || !end_loop )
  {
    CompositeCoEdge* prev;
    CompositeLoop* loop;
    if( start_loop )
    {
      loop = start_loop;
      prev = start_prev_coedge;
    }
    else 
    {
      loop = end_loop;
      prev = end_prev_coedge;
    }
    
    if( coedge1->start_point() == prev->end_point() )
    {
      loop->insert_after( coedge1, prev );
      loop->insert_after( coedge2, coedge1 );
    }
    else
    {
      assert( coedge2->start_point() == prev->end_point() );
      loop->insert_after( coedge2, prev );
      loop->insert_after( coedge1, coedge2 );
    }
  }
    
    // If the end points of the curve intersected different
    // loops, combine the loops such that the curve becomes
    // a "bridge" between them.
  else if( start_loop != end_loop )
  {
    CompositeCoEdge* prev = start_prev_coedge;
    CompositeCoEdge* coedge = end_next_coedge;
    CompositeCoEdge* next = 0;
    
      // Which of the two coedges for the curve we are
      // restoring do we want to begin with (and store
      // the other as other_coedge).
    CompositeCoEdge* other_coedge;
    if( coedge1->start_point() == prev->end_point() )
    {
      start_loop->insert_after( coedge1, prev );
      prev = coedge1;
      other_coedge = coedge2;
    }
    else if( coedge2->start_point() == prev->end_point() )
    {
      start_loop->insert_after( coedge2, prev );
      prev = coedge2;
      other_coedge = coedge1;
    }
    else assert( 0 );
    
    while( end_loop->first_coedge() ) // while loop has coedges
    {
      next = end_loop->next_coedge( coedge );
      end_loop->remove( coedge );
      start_loop->insert_after( coedge, prev );
      prev = coedge;
      coedge = next;
    }
    
      // The other coedge for the curve we are restoring...
    start_loop->insert_after( other_coedge, prev );
    
    assert( end_loop->num_coedges() == 0 );
    compsurf->remove( end_loop );
    delete end_loop;
  }
    
    // If both end points of the curve intersected the same
    // loop, then split the loop (and the composite surface)
    // into two.
  else
  {
    assert( start_loop == end_loop );

      // Special case:
      // Hole intersecting original loop at one point
      // (when loop is split, oringinal loop has all
      //  original coedges, and new loop has only the
      //  curve we are restoring)  Figure out which
      // coedge belongs in the hole.
    if( coedge1->start_point() == coedge1->end_point() )
    {
      assert( start_next_coedge == end_next_coedge );
      assert( start_prev_coedge == end_prev_coedge );
      
      CubitVector prev_tan, coe1_tan, coe2_tan, norm, junk;
      CubitVector point = coedge1->start_point()->coordinates();
      start_prev_coedge->get_curve()->closest_point( point, junk, &prev_tan );
      if( start_prev_coedge->sense() == CUBIT_FORWARD ) // yes, forward!!!
      	prev_tan *= -1.0;
      coedge1->get_curve()->closest_point( point, junk, &coe1_tan );
      coe2_tan = coe1_tan;
      if( coedge1->sense() == CUBIT_REVERSED )
      	coe1_tan *= -1.0;
      if( coedge2->sense() == CUBIT_REVERSED )
      	coe2_tan *= -1.0;
      compsurf->closest_point( point, 0, &norm );
      	
      double angle1 = norm.vector_angle( prev_tan, coe1_tan );
      double angle2 = norm.vector_angle( prev_tan, coe2_tan );
      
      if( angle2 < angle1 )
      {
      	CompositeCoEdge* temp = coedge2;
        coedge2 = coedge1;
        coedge1 = temp;
      }
    }
    
      // Normal case (not a hole)
    else
    {
        // Make sure coedge1 is the reverse one
      if( coedge1->sense() == CUBIT_FORWARD )
        std::swap(coedge1,coedge2);
    }
    
    end_loop = new CompositeLoop();
    start_loop->insert_after( coedge1, end_prev_coedge );
    end_loop->insert_after( coedge2, 0 );
    
    CompositeCoEdge* coedge = end_next_coedge;
    CompositeCoEdge* prev = coedge2;
    while( coedge != start_next_coedge )
    {
      CompositeCoEdge* next = start_loop->next_coedge( coedge );
      start_loop->remove( coedge );
      end_loop->insert_after( coedge, prev );
      prev = coedge;
      coedge = next;
    }
    
    new_surf = split_surface(  compsurf, start_loop, end_loop );
    
  }
  
  if( new_surf )
  {
    if( compsurf->next_co_surface() )
    {
      CompositeCoSurf* cos = 0;
      while( (cos = compsurf->next_co_surface(cos)) )
        cos->get_shell()->add( new_surf, cos->sense() );
    }
    
    if( ! new_surf->has_hidden_entities() &&
        ! new_surf->next_co_surface(0) )
      remove_composite( new_surf );
  }
  
  if( ! compsurf->has_hidden_entities() &&
      ! compsurf->next_co_surface(0) )
    remove_composite( compsurf );
  
  return CUBIT_SUCCESS;
}


//-------------------------------------------------------------------------
// Purpose       : Given a point that intersects a loop in the passed 
//                 surface, find the location in that loop at which the
//                 passed curve should be inserted.  If there are only two
//                 coedges in the loop at the passed point, then the answer
//                 is obvious.  If there is a sipe at the point, then use
//                 topological information of underlying real surfaces to
//                 determine the two coedges between which the new curve
//                 should be inserted.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/19/03
//-------------------------------------------------------------------------
CubitStatus CompositeEngine::find_coedges( CompositeSurface* surface,
                                           CompositeCurve* curve,
                                           CompositePoint* point,
                                           CompositeCoEdge*& previous,
                                           CompositeCoEdge*& next )
{
  const char* const bad_loop_message = "Internal error: Invalid loop. (%s:%d)\n";
  
    // Find list of all coedges around passed point
    // and in passed surface.
  DLIList<CompositeCoEdge*> point_coedges;
  CompositeCurve* pt_curve = 0;
  while ( (pt_curve = point->next_curve(pt_curve)) )
  {
    CompositeCoEdge* coedge = 0;
    while ( (coedge = pt_curve->next_coedge(coedge)) )
    {
      if (coedge->get_loop() && coedge->get_loop()->get_surface() == surface)
        point_coedges.append(coedge);
    }
  }
  
    // Point is at end of a sipe/hardline
  if ( point_coedges.size() == 0 )
  {
    previous = next = 0;
    return CUBIT_SUCCESS;
  }
  
    // One coedge - closed curve
  if ( point_coedges.size() == 1 &&
       point_coedges.get()->start_point() ==
       point_coedges.get()->end_point() )
  {
    previous = next = point_coedges.get();
    return CUBIT_SUCCESS;
  }
  
    // Broken loop?
  if ( point_coedges.size() % 2 != 0 )
  {
    PRINT_ERROR(bad_loop_message,__FILE__,__LINE__);
    assert(0);
    return CUBIT_FAILURE;
  }
  
    // If only two coedges, then we're done
  if ( point_coedges.size() == 2 )
  {
    previous = point_coedges.get();
    next = point_coedges.next();
    
    if ( previous->start_point() == point )
      std::swap(previous, next);
    
    return CUBIT_SUCCESS;
  }
  
    // Find previous/next coedges by using order of
    // coedges about point in underlying surfaces.
    // Coedges must occur in clock-wise order about point.
  
    // Get the real curve at the point
  int curve_index = point == curve->start_point() ? 0 : curve->num_curves() - 1;
  Curve* real_curve = curve->get_curve(curve_index);
  
    // Get two coedges in composite
  DLIList<TopologyBridge*> coedge_bridges, loop_bridges(1), surface_bridges(1);
  real_curve->get_parents_virt(coedge_bridges);
  DLIList<CoEdgeSM*> curve_coedges;
  while (coedge_bridges.size())
  {
    TopologyBridge* coe_bridge = coedge_bridges.pop();
    
    loop_bridges.clean_out();
    surface_bridges.clean_out();
    
    coe_bridge->get_parents_virt(loop_bridges);
    assert(loop_bridges.size() == 1);
    
    loop_bridges.get()->get_parents_virt(surface_bridges);
    assert(surface_bridges.size() == 1);
    
    if ( surface->contains_bridge(surface_bridges.get()) )
    {
      curve_coedges.append(dynamic_cast<CoEdgeSM*>(coe_bridge));
    }
  }
  if ( curve_coedges.size() != 2 )
  {
    PRINT_ERROR(bad_loop_message,__FILE__,__LINE__);
    assert(0);
    return CUBIT_FAILURE;
  }
  
    // Assign coedges such that prev_coe_sm goes "out of"
    // the point and next_coe_sm goes "into" the point.
    // (Unless surface sense is reversed in composite surface.)
  curve_coedges.reset();
  CoEdgeSM *prev_coe_sm = curve_coedges.get();
  CoEdgeSM *next_coe_sm = curve_coedges.next();
  CompositeCoEdge* prev_owner = dynamic_cast<CompositeCoEdge*>(prev_coe_sm->owner());
  CompositeCoEdge* next_owner = dynamic_cast<CompositeCoEdge*>(next_coe_sm->owner());
  assert( prev_owner && next_owner );
  if ( prev_owner->start_point() == point )
  {
    assert(next_owner->end_point() == point);
  }
  else
  {
    assert(next_owner->start_point() == point);
    assert(prev_owner->end_point() == point);
    std::swap(prev_coe_sm, next_coe_sm);
  }
    
    // Now iterate through real coedges around point to
    // find the next/prev composite coedges.
  CompositeCoEdge* prev_result 
    = find_next_point_coedge( surface, prev_coe_sm, point, point_coedges );
  CompositeCoEdge* next_result 
    = find_next_point_coedge( surface, next_coe_sm, point, point_coedges );
  
  if( !prev_result || !next_result )
  {
    PRINT_ERROR("Internal error: Invalid composite. (%s:%d)\n",__FILE__,__LINE__);
    return CUBIT_FAILURE;
  }

  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Find next coedge by traversing real coedges on
//                 surfaces hidden by composite.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/19/03
//-------------------------------------------------------------------------
CompositeCoEdge* CompositeEngine::find_next_point_coedge( 
                                     CompositeSurface* const compsurf,
                                     CoEdgeSM* const first_coedge,
                                     CompositePoint* point,
                                     DLIList<CompositeCoEdge*>& point_coedges )
{
  DLIList<TopologyBridge*> loop_bridges(1), coedge_bridges, curve_pts(2);
  CoEdgeSM* coedge = first_coedge;
  
  coedge_bridges.clean_out();
  coedge->get_children(coedge_bridges,true,COMPOSITE_LAYER-1);
  assert(coedge_bridges.size() == 1);
  Curve* curvesm = dynamic_cast<Curve*>(coedge_bridges.get());
 
  do
  {
      // Get loop from coedge
    loop_bridges.clean_out();
    coedge->get_parents_virt(loop_bridges);
    assert(loop_bridges.size() == 1);
    LoopSM* loop_sm = dynamic_cast<LoopSM*>(loop_bridges.get());
    
      // Make sure we're still inside the composite.
      // It's a bug if this check fails.
    loop_bridges.clean_out();
    loop_sm->get_parents_virt(loop_bridges);
    if ( !compsurf->contains_bridge(loop_bridges.get()) )
      { assert(0); break; }
    
      // Get direction of curve wrt point : forward if the
      // curve ends at the point or reverse if it begins at
      // the point.
    curve_pts.clean_out();
    curvesm->get_children(curve_pts,true,COMPOSITE_LAYER-1);
    curve_pts.reset();
    CubitSense curve_sense;
    if ( curve_pts.next()->owner() == point )
      curve_sense = CUBIT_FORWARD;
    else if( curve_pts.get()->owner() == point )
      curve_sense = CUBIT_REVERSED;
    else
      { assert(0); break; } //bug - shouldn't happen
      
      // Get the next coedge in the loop. If the coedge
      // begins at the point, then the adjacent coedge
      // is the previous one in the loop.  Otherwise
      // it is the next one in the loop.
    coedge_bridges.clean_out();
    loop_sm->get_children(coedge_bridges,true,COMPOSITE_LAYER-1);
    coedge_bridges.move_to(coedge);
    assert(coedge_bridges.get() == coedge);
    if ( coedge->sense() == curve_sense )
      coedge = dynamic_cast<CoEdgeSM*>(coedge_bridges.next());
    else
      coedge = dynamic_cast<CoEdgeSM*>(coedge_bridges.prev());

      // Are we done?  This is the real end condition.
      // Any other return path is an error.
    CompositeCoEdge* result = dynamic_cast<CompositeCoEdge*>(coedge->owner());
    assert(!!result);
    if ( point_coedges.is_in_list(result) )
      return result;

      // get the curve for the current coedge
    coedge_bridges.clean_out();
    coedge->get_children(coedge_bridges,true,COMPOSITE_LAYER-1);
    assert(coedge_bridges.size() == 1);
    curvesm = dynamic_cast<Curve*>(coedge_bridges.get());

      // Get the other coedge on the curve
    coedge_bridges.clean_out();
    curvesm->get_parents_virt(coedge_bridges);
    if ( coedge_bridges.size() != 2 )
      break;  // Reached boundary of composite
    coedge_bridges.move_to( coedge );
    assert(coedge_bridges.get() == coedge );
    coedge = dynamic_cast<CoEdgeSM*>(coedge_bridges.next());
  
  } while( coedge != first_coedge );
  
  return 0;
}
  
  



//-------------------------------------------------------------------------
// Purpose       : 
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 
//-------------------------------------------------------------------------
CompositeSurface* CompositeEngine::split_surface( 
                                 CompositeSurface* surf_to_split,
                                 CompositeLoop* loop_on_orig,
                                 CompositeLoop* loop_on_new )
{
  int i;
  
  DLIList<CoEdgeSM*> front_list;
  DLIList<TopologyBridge*> bridge_list, loop_list, coedge_list;
  DLIList<TopologyBridge*> hidden_loops_to_move;
  DLIList<Surface*> surfs_to_move;
  
    // Do advancing front across surfaces underlying composite surface
    // to find which must be moved to the new composite.
    
    // Begin with the coedges in loop_on_new.
  CompositeCoEdge* coedge = loop_on_new->first_coedge();
  do
  { 
    for( i = 0; i < coedge->num_coedges(); i++ )
      front_list.append( coedge->get_coedge(i) );
    coedge = loop_on_new->next_coedge( coedge );
  } while( coedge != loop_on_new->first_coedge() );
  
    // Loop until front_list is empty
  while( front_list.size() )
  {
    CoEdgeSM* coe_real = front_list.pop();
    
      // Loop and Surface from CoEdge
    bridge_list.clean_out();
    coe_real->get_parents_virt( bridge_list );
    assert( bridge_list.size() == 1 );
    TopologyBridge* loop = bridge_list.get();
    bridge_list.clean_out();
    loop->get_parents_virt( bridge_list );
    if( bridge_list.size() != 1 )
    {
      assert( loop == loop_on_new );
      assert( bridge_list.size() == 0 );
      continue;
    }
    
    Surface* surf = dynamic_cast<Surface*>(bridge_list.get());
    assert( surf!= NULL  );
    
    if( surf->owner() != surf_to_split )
      continue;
      
      // Add surface to list of Surfaces to move to new composite
    surfs_to_move.append_unique( surf );  
    
      // Get list of all surface coedges (coedge_list)
      // and all loops to hidden_loops_to_move
    loop_list.clean_out();
    surf->get_children_virt( loop_list );
    for( i = loop_list.size(); i--; )
    {
      TopologyBridge* loop_bridge = loop_list.get_and_step();
      if( hidden_loops_to_move.is_in_list( loop_bridge ) )
        continue;
      hidden_loops_to_move.append( loop_bridge );
      bridge_list.clean_out();
      loop_bridge->get_children( bridge_list, true, COMPOSITE_LAYER-1 );
      coedge_list += bridge_list;
    }
    
      // For each coedge on the surface (all coedges in coedge_list),
      // get the curve and search for any other coedges that are 
      // parents of that curve and children of a surface hidden by
      // the composite.  Add them to front_list.
    while( coedge_list.size() )
    {
      TopologyBridge* coe_bridge = coedge_list.pop();
      bridge_list.clean_out();
      coe_bridge->get_children_virt( bridge_list );
      assert( bridge_list.size() == 1 );
      TopologyBridge* curve_bridge = bridge_list.get();
      bridge_list.clean_out();
      curve_bridge->get_parents_virt( bridge_list );
      
      while( bridge_list.size() )
      {
      	TopologyBridge* other_coe = bridge_list.pop();
	      if( other_coe == coe_bridge ) 
	        continue;
	
	      CompositeCoEdge* compcoe = dynamic_cast<CompositeCoEdge*>(other_coe);
	      if( !compcoe )
	      {
	        compcoe = dynamic_cast<CompositeCoEdge*>(other_coe->owner());
	      }
	      if( !compcoe || compcoe->get_loop() != loop_on_orig )
	        front_list.append( dynamic_cast<CoEdgeSM*>( other_coe ) );
      }
    }
  }
  
    // Check if surface needs to be split.  A new, closed loop
    // does not always indicate the surface needs to be split.
    // See PR#2140 for an example.
  if (surfs_to_move.size() == surf_to_split->num_surfs())
  {
    surf_to_split->add( loop_on_new );
    return 0;
  }
  
    // Split the composite (pass in indices of surfaces that 
    // should be in new composite rather than old.)
  VGArray<int> index_array( surfs_to_move.size() );
  for( i = 0; i < surfs_to_move.size(); i++ )
    index_array[i] = surf_to_split->index_of( surfs_to_move.get_and_step() );
  CompositeSurface* new_surf = surf_to_split->split( index_array );
  if( !new_surf )
  {
    surf_to_split->add( loop_on_new );
    return 0;
  }
  
  assert( new_surf->num_surfs() && surf_to_split->num_surfs() );
  new_surf->add( loop_on_new );

    // Move hidden coedges from old surface's hidden set to
    // new surface.
  HiddenEntitySet* new_set = &(new_surf->hidden_entities());
  HiddenEntitySet* old_set = &(surf_to_split->hidden_entities());
  for ( i = 0; i < new_surf->num_surfs(); i++ )
  {
    loop_list.clean_out();
    new_surf->get_surface(i)->get_children_virt(loop_list);
    for ( int j = loop_list.size(); j--; )
    {
      coedge_list.clean_out();
      loop_list.get_and_step()->get_children( coedge_list, true, COMPOSITE_LAYER );
      for ( int k = coedge_list.size(); k--; )
      {
        TopologyBridge* coedge_ptr = coedge_list.get_and_step();
        if ( coedge_ptr->owner() == old_set )
        {
          old_set->restore(coedge_ptr);
          new_set->hide(coedge_ptr);
        }
      }
    }
  }
  

    // Move hidden curves from old surface's hidden set to
    // new surface.
  DLIList<Curve*> hidden_curves;
  old_set->hidden_curves( hidden_curves );
  for( i = hidden_curves.size(); i--; )
  {
    Curve* curve = hidden_curves.get_and_step();
    CompositeCurve* comp_curve = dynamic_cast<CompositeCurve*>(curve);
    assert(!!comp_curve);
    
    bool all_new = true;
    bool all_old = true;
    CompositeCoEdge* coedge_ptr = 0;
    while( (coedge_ptr = comp_curve->next_coedge(coedge_ptr)) )
    {
      if( coedge_ptr->owner() != new_set )
        all_new = false;
      if( coedge_ptr->owner() != old_set )
        all_old = false;
    }
    if( ! all_old )
      old_set->restore( curve );
    if( all_new )
      new_set->hide( curve );
  }

    // Move hidden points from old surface's hidden set to
    // new surface.
  DLIList<TBPoint*> hidden_points;
  old_set->hidden_points( hidden_points );
  for( i = hidden_points.size(); i--; )
  {
    TBPoint* point = hidden_points.get_and_step();
    CompositePoint* comp_pt = dynamic_cast<CompositePoint*>(point);
    assert(!!comp_pt);
    
    CompositeCurve* curve = 0;
    bool all_new = true;
    bool all_old = true;
    while ( (curve = comp_pt->next_curve(curve)) )
    {
      if ( curve->owner() != new_set )
        all_new = false;
      if ( curve->owner() != old_set )
        all_old = false;
    }
    
    if( ! all_old )
      old_set->restore( point );
    if( all_new )
      new_set->hide( point );
  }
  


    // figure out which visible loops need to be moved
  DLIList<CompositeLoop*> loops_to_move;
  for( CompositeLoop* loop = surf_to_split->first_loop();
       loop != 0;
       loop = loop->next_loop() )
  {
    CompositeCoEdge* a_coedge = loop->first_coedge();
    CoEdgeSM* real_coe = a_coedge->get_coedge(0);
    assert( real_coe!= NULL  );
    bridge_list.clean_out();
    real_coe->get_parents_virt( bridge_list );
    assert( bridge_list.size() == 1 );
    TopologyBridge* loop_bridge = bridge_list.get();
    bridge_list.clean_out();
    loop_bridge->get_parents_virt( bridge_list );
    Surface* surf = dynamic_cast<Surface*>(bridge_list.get());
    assert( surf!= NULL  );
    if( surf->owner() == new_surf )
      loops_to_move.append( loop );
  }
  
    // move visible loops to new surface
  while( loops_to_move.size() )
  {
    CompositeLoop* loop = loops_to_move.pop();
    surf_to_split->remove( loop );
    new_surf->add( loop );
  }
  
  assert( loop_on_orig->get_surface() == surf_to_split );  
  return new_surf;
}


    
//-------------------------------------------------------------------------
// Purpose       : Composite Curves
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 
//-------------------------------------------------------------------------
CompositeCurve* CompositeEngine::remove_point( TBPoint* dead_point,
                                               bool remove_partition,
                                               Curve* survivor )
{ 
  if( dynamic_cast<CompositePoint*>(dead_point->owner()) )
    dead_point = dynamic_cast<TBPoint*>(dead_point->owner());
 
  DLIList<TopologyBridge*> query_results;
  dead_point->get_parents_virt( query_results );
  if( query_results.size() != 2 )
    return 0;
  
  Curve* curve1 = dynamic_cast<Curve*>(query_results.get_and_step());
  Curve* curve2 = dynamic_cast<Curve*>(query_results.get_and_step());
  assert( curve1 && curve2 );
  
  query_results.clean_out();
  curve1->get_children( query_results, COMPOSITE_LAYER );
  assert( query_results.size() == 2 && query_results.move_to( dead_point ) );
  query_results.move_to( dead_point );
  TBPoint* other1 = dynamic_cast<TBPoint*>(query_results.step_and_get());
  
  query_results.clean_out();
  curve2->get_children( query_results, COMPOSITE_LAYER );
  assert( query_results.size() == 2 && query_results.move_to( dead_point ) );
  query_results.move_to( dead_point );
  TBPoint* other2 = dynamic_cast<TBPoint*>(query_results.step_and_get());
  
    // If curves form a closed two-curve loop, we need to 
    // specify which point to keep.
  TBPoint* keep = other1 == other2 ? other1 : 0;
  if (survivor == curve2) std::swap(curve1, curve2);
  return composite( curve1, curve2, keep, remove_partition );
}

//-------------------------------------------------------------------------
// Purpose       : Remove curve partitions beneath a composite curve
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/11/03
//-------------------------------------------------------------------------
CubitStatus CompositeEngine::remove_partition_point( CompositePoint* comp )
{
  PartitionPoint* pt = dynamic_cast<PartitionPoint*>(comp->get_point());
  if( !pt || pt->real_point() )
    return CUBIT_SUCCESS;
  
  if( comp->next_curve() )
    return CUBIT_FAILURE;
  
  if ( !PartitionEngine::instance().remove_point( pt ) )
    return CUBIT_FAILURE;
  
  assert( !comp->get_point() );
  clean_out_deactivated_geometry();
  
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Remove surface partitions beneath a composite surface
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/11/03
//-------------------------------------------------------------------------
CubitStatus CompositeEngine::remove_partition_curves( CompositeCurve* curve )
{
  if( ! dynamic_cast<HiddenEntitySet*>(curve->owner()) )
    { assert(0); return CUBIT_FAILURE; }
  
  DLIList<CompositeCurve*> hidden_curves;
  hidden_curves.append( curve );
    
    // Split composite such that each partition curve to be removed
    // is owned by a seperate composite curve.  Put the partitition
    // curves in dead_curves and their owning composite curves in
    // dead_composites.
  int i;
  CubitStatus result = CUBIT_SUCCESS;
  DLIList<PartitionCurve*> dead_curves;
  while( hidden_curves.size() )
  {
    CompositeCurve* comp = hidden_curves.pop();
    
      // If this curve has a single underlying curve,
      // we're done with it.  Add the underlying curve
      // to the dead curve list if it is a PartitionCurve
      // and move on to the next composite curve.
    if( comp->num_curves() == 1 )
    {
      SegmentedCurve* segcurve = dynamic_cast<SegmentedCurve*>(comp->get_curve(0));
      if( segcurve )
      {
        dead_curves.append(segcurve);
      }
      continue;
    }
    
      // Search for a partition curve in the composite
    for ( i = 1; i < comp->num_curves(); i++ )
      if( dynamic_cast<SegmentedCurve*>(comp->get_curve(i)) )
        break;
        
      // Composite doesn't contain any partition curves.
      // Move on to the next composite curve.
    if( i == comp->num_curves() )
      continue;
    
      // Split composite curve at the partition curve. 
    Curve* r1, *r2;
    if ( !split( comp, i, r1, r2 ) )
    {
      result = CUBIT_FAILURE;
      continue;
    }
    
      // Add the resulting composites back onto hidden_curves.
      // Continue processing them until all partitions are in their
      // own composite, or the composite contains no more partitions.
    CompositeCurve* comp1 = dynamic_cast<CompositeCurve*>(r1);
    CompositeCurve* comp2 = dynamic_cast<CompositeCurve*>(r2);
    assert(comp1 && comp2);
    hidden_curves.append(comp1);
    hidden_curves.append(comp2);
  }
  
  
    // Now actually remove the partition curves
  while ( dead_curves.size() )
    if ( ! PartitionEngine::instance().remove_curve( dead_curves.pop() ) )
      result = CUBIT_FAILURE;
  
    // Last, delete the remaining composite curves
  clean_out_deactivated_geometry();
  
  return result;
}

//-------------------------------------------------------------------------
// Purpose       : Create a point-curve given a point hidden by a composite
//                 surface.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 02/27/04
//-------------------------------------------------------------------------
CompositeCurve* CompositeEngine::restore_point_in_surface( TBPoint* pt )
{
    // If real point is hidden by a composite surface, it should
    // have already been replaced by a composite point. Get the
    // composite point.
  CompositePoint* point = dynamic_cast<CompositePoint*>(pt);
  if (!point && !(point = dynamic_cast<CompositePoint*>(pt->owner())))
    return NULL;
  
    // Get the surface hiding the point.
  HiddenEntitySet* hidden_set = dynamic_cast<HiddenEntitySet*>(point->owner());
  if (!hidden_set)
    return NULL;
    
  CompositeSurface* surface = dynamic_cast<CompositeSurface*>(hidden_set->owner());
  if (!surface)
    return NULL;
  
    // Check if the point is already owned by a point curve.
  CompositeCurve* curve = 0;
  while ( (curve = point->next_curve(curve)) )
    if (curve->geometry_type() == POINT_CURVE_TYPE)
      return restore_curve(curve) ? curve : NULL;
  
    // Construct null-geometry composite curve
  surface->hidden_entities().restore(point);
  curve = new CompositeCurve(point);
  CompositeCoEdge* coedge = new CompositeCoEdge(curve);
  CompositeLoop* loop = new CompositeLoop();
  loop->insert_after( coedge, NULL );
  surface->add(loop);
  return curve;
}
  
//-------------------------------------------------------------------------
// Purpose       : Create a new cosurface with appropriate sense
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/27/04
//-------------------------------------------------------------------------
static CompositeCoSurf* cme_create_cosurf( CompositeLump* lump,
                                           CompositeSurface* surf )
{
  Surface* real_surf = surf->get_surface(0);
  
  DLIList<TopologyBridge*> shells, lumps(1);
  real_surf->get_parents( shells );
  while (shells.size())
  {
    TopologyBridge* shell = shells.pop();
    lumps.clean_out();
    shell->get_parents( lumps );
    assert( lumps.size() == 1 );
    TopologyBridge* real_lump = lumps.pop();
    if (real_lump->owner() == lump)
    {
      CubitSense sense = real_surf->get_shell_sense( dynamic_cast<ShellSM*>(shell) );
      if (sense == CUBIT_UNKNOWN)
        sense = CUBIT_FORWARD;
      else if (surf->get_sense(0) == CUBIT_REVERSED)
        sense = (sense == CUBIT_REVERSED) ? CUBIT_FORWARD : CUBIT_REVERSED;
    
      CompositeCoSurf* result = new CompositeCoSurf(sense);
      surf->add( result );
      return result;
    }
  }
  
  return NULL;
}

//-------------------------------------------------------------------------
// Purpose       : Restore a surface hidden by a composite volume
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/27/04
//-------------------------------------------------------------------------
CubitStatus CompositeEngine::restore_surface( Surface* surf,
                                              Surface*& stitch_partner )
{ 
    // find composite lump owning the curve
  TopologyBridge* bridge = surf;
  HiddenEntitySet* owner_set = NULL;
  while( bridge && !(owner_set = dynamic_cast<HiddenEntitySet*>(bridge->owner())) )
  {
    bridge = dynamic_cast<TopologyBridge*>(bridge->owner());
  }
  
  if( !owner_set )
    return CUBIT_FAILURE;
  
  CompositeLump* lump = dynamic_cast<CompositeLump*>(owner_set->owner());
  if( !lump )
    return CUBIT_FAILURE;
  
  CompositeSurface* surf1 = dynamic_cast<CompositeSurface*>(surf);
  if( !surf1 )
    surf1 = dynamic_cast<CompositeSurface*>(surf1->owner());
  if( !surf1 )
    return CUBIT_FAILURE;
  
    // Check if and curves are hidden by a CompositeSurface, and
    // if so, split the composite surface so that the curve(s)
    // are no longer hidden.
  HiddenEntitySet* hs;
  CompositeLoop* loop;
  CompositeCoEdge* coedge;
  for (loop = surf1->first_loop(); loop; loop = surf1->next_loop(loop))
  {
    coedge = loop->first_coedge();
    do {
      CompositeCurve* curv = coedge->get_curve();
      if ( (hs = dynamic_cast<HiddenEntitySet*>(curv->owner()))
        && dynamic_cast<CompositeSurface*>(hs->owner())
        && !restore_curve(curv))
        return CUBIT_FAILURE;
    } while ((coedge = loop->next_coedge(coedge)) != loop->first_coedge());
  } 
  
    // Un-hide surf's and children
  CompositeSurface* surf2 = surf1->unstitch();
  cme_unhide_surface( surf1 );
  if (surf2 != surf)
    cme_unhide_surface( surf2 );
  
    // Find affected shells.
  bool all_connected = true;
  DLIList<CompositeShell*> modified_shell_list;
  DLIList<CompositeCurve*> curve_list;
  for (loop = surf1->first_loop(); loop; loop = surf1->next_loop(loop))
  {
    coedge = loop->first_coedge();
    do {
      CompositeCurve* curv = coedge->get_curve();
      curv->get_stitched( curve_list );
      bool found_adj_surf = false;
      while (curve_list.size())
      {
        curv = curve_list.pop();
        CompositeCoEdge* curve_coedge = 0;
        while ((curve_coedge = curv->next_coedge(curve_coedge)))
        {
          if (!curve_coedge->get_loop())
            continue;

          CompositeSurface* adj_surf = curve_coedge->get_loop()->get_surface();
          if (adj_surf == surf1 || adj_surf == surf2)
            continue;
            
          CompositeCoSurf* cosurf = 0;
          while ((cosurf = adj_surf->next_co_surface( cosurf )))
          {
            if (cosurf->get_shell()->get_lump() == lump)
            {
              found_adj_surf = true;
              modified_shell_list.append_unique( cosurf->get_shell() );
            }
          }
        }
      }
      
      if (!found_adj_surf)
        all_connected = false;
        
    } while ((coedge = loop->next_coedge( coedge )) != loop->first_coedge());
  }
  
  
    // Create Co-Surfacs
  CompositeCoSurf *cosurf2, *cosurf1 = cme_create_cosurf( lump, surf1 );
  if (surf1 == surf2)
  {
    CubitSense sense = cosurf1->sense() == CUBIT_FORWARD ?
                       CUBIT_REVERSED : CUBIT_FORWARD;
    cosurf2 = new CompositeCoSurf( sense );
    surf2->add( cosurf2 );
  }
  else
  {
    cosurf2 = cme_create_cosurf( lump, surf2 );
  }


  CompositeShell* shell_to_split = 0;
  if (modified_shell_list.size() == 0)
  {
    // No adjacent shells -- create void
    CompositeShell* shell_to_split = new CompositeShell();
    lump->add( shell_to_split );
    
    // If surface is closed, set all_connected to true
    // to indicate that the volume must be split (to create
    // the void)
    loop = 0;
    all_connected = true;
    while ((loop = surf1->next_loop(loop)))
    {
      coedge = loop->first_coedge();
      do {
        CompositeCurve* curv = coedge->get_curve();
        CompositeCoEdge* curv_coedge = 0;
        bool closed = false;
        while ((curv_coedge = curv->next_coedge(curv_coedge)))
        {
          if (curv_coedge != coedge &&
              curv_coedge->get_loop() &&
              (curv_coedge->get_loop()->get_surface() == surf1 ||
               curv_coedge->get_loop()->get_surface() == surf2))
          {
            closed = true;
            break;
          }
        }
        if (!closed)
          all_connected = false;
      } while ((coedge = loop->next_coedge(coedge)) != loop->first_coedge());
    }
  }
  else
  {
    // otherwise combine all connected shells
    modified_shell_list.reverse();
    shell_to_split = modified_shell_list.pop();
    while (modified_shell_list.size())
    {
      CompositeShell* shell = modified_shell_list.pop();
      while (CompositeCoSurf* cosurf = shell->first_co_surf())
      {
        shell->remove( cosurf );
        shell_to_split->add( cosurf );
      }
      lump->remove( shell );
      delete shell;
    }
  }
  
  shell_to_split->add( cosurf1 );
  shell_to_split->add( cosurf2 );
  
  stitch_partner = surf1 == surf2 ? NULL : surf2;
  if (all_connected)
  {
    CompositeLump* new_lump = split_lump( shell_to_split );
    lump->get_body()->add( new_lump );
    //if (!new_lump->has_hidden_entities())
    //  remove_composite( new_lump );
  }
//  if (!lump->has_hidden_entities())
//    remove_composite( lump );
  
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Test if a CompositeShell needs to be split, and if
//                 so, split it and the owning lump.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/27/04
//-------------------------------------------------------------------------
CompositeLump* CompositeEngine::split_lump( CompositeShell* shell_to_split )
{
  int i, j;
  
    // split the shell
  CompositeShell* new_shell = split_shell( shell_to_split );
  if (!new_shell)
    return shell_to_split->get_lump();
    
    // Get list of shells to move to new lump
  shell_to_split->get_lump()->add( new_shell );
  DLIList<CompositeShell*> shell_list, shells_to_move;
  shell_list.append( new_shell );
  shells_to_move.append( new_shell );
  while (shell_list.size())
  {
    CompositeShell* shell = shell_list.pop();
    CompositeCoSurf* cosurf = 0;
    while ((cosurf = shell->next_co_surf( cosurf )))
    {
      CompositeSurface* surf = cosurf->get_surface();
      surf = surf->get_stitch_partner();
      if (!surf)
        continue;
      
      CompositeCoSurf* surf_cosurf = 0;
      while ((surf_cosurf = surf->next_co_surface( surf_cosurf )))
      {
        CompositeShell* surf_shell = surf_cosurf->get_shell();
        if (surf_shell == shell_to_split ||
            surf_shell->get_lump() != shell_to_split->get_lump() ||
            shells_to_move.is_in_list( surf_shell ))
          continue;
        
        shell_list.append( surf_shell );
        shells_to_move.append( surf_shell );
      }
    }
  }
  
    // Get list of real lumps defining composite lump that
    // are to be moved to the new composite lump
  DLIList<Lump*> lumps_to_move;
  DLIList<TopologyBridge*> shells, lumps;
  for (i = shells_to_move.size(); i--; )
  {
    CompositeShell* shell = shells_to_move.get_and_step();
    CompositeCoSurf* cosurf = 0;
    while ((cosurf = shell->next_co_surf( cosurf )))
    {
      CompositeSurface* surf = cosurf->get_surface();
      for (j = 0; j < surf->num_surfs(); j++)
      {
        Surface* real_surf = surf->get_surface(j);
        real_surf->get_parents_virt( shells );
        while (shells.size())
        {
          shells.pop()->get_parents_virt( lumps );
          assert( lumps.size() == 1 );
          TopologyBridge* lump = lumps.pop();
          if (lump->owner() == shell_to_split->get_lump())
            lumps_to_move.append_unique( dynamic_cast<Lump*>(lump) );
        }
      }
    }
  }
  
  if (lumps_to_move.size() == shell_to_split->get_lump()->num_lumps())
    return shell_to_split->get_lump();
  
    // Split composite lump and move shells to new lump
  VGArray<int> vol_indices( lumps_to_move.size() );
  lumps_to_move.reset();
  for (i = lumps_to_move.size(); i--; )
    vol_indices[i] = shell_to_split->get_lump()->index_of( lumps_to_move.next(i) );
  CompositeLump* new_lump = shell_to_split->get_lump()->split( vol_indices );
  new_lump->add( new_shell );
  shells_to_move.reverse();
  while (shells_to_move.size())
  {
    CompositeShell* shell = shells_to_move.pop();
    shell_to_split->get_lump()->remove( shell );
    new_lump->add( new_shell );
  }
  
    // Move any hidden entities from old composite to new
  HiddenEntitySet* old_set = &shell_to_split->get_lump()->hidden_entities();
  DLIList<TopologyBridge*> surfs;
  while (lumps_to_move.size())
  {
    Lump* lump = lumps_to_move.pop();
    lump->get_children( shells, true, COMPOSITE_LAYER-1 );
    while (shells.size())
    {
      TopologyBridge* shell = shells.pop();
      shell->get_children( surfs, true, COMPOSITE_LAYER );
      while (surfs.size())
      {
        CompositeSurface* surf = dynamic_cast<CompositeSurface*>(surfs.pop());
        if (surf->owner() != old_set)
          continue;
          
        cme_unhide_surface( surf );
        cme_hide_surface( new_lump->hidden_entities(), surf );
      }
    }
  }

  return new_lump;
}
  


//-------------------------------------------------------------------------
// Purpose       : Recreate composites stored in attributes
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/18/02
//-------------------------------------------------------------------------
CubitStatus CompositeEngine::import_geometry( DLIList<TopologyBridge*>& imported_geom )
{
  CubitStatus result = CUBIT_SUCCESS;
  int i;
  
  DLIList<BodySM*>  bodies;
  DLIList<Lump*>    lumps,    temp_lumps;
  DLIList<Surface*> surfaces, temp_surfaces;
  DLIList<Curve*>   curves,   temp_curves;
  DLIList<TBPoint*>  points,   temp_points;
  
  CAST_LIST( imported_geom, bodies,   BodySM  );
  CAST_LIST( imported_geom, lumps,    Lump    );
  CAST_LIST( imported_geom, surfaces, Surface );
  CAST_LIST( imported_geom, curves,   Curve   );
  CAST_LIST( imported_geom, points,   TBPoint   );

  i = bodies.size() + lumps.size() + surfaces.size() + curves.size() + points.size();
  if( i != imported_geom.size() )
  {
    PRINT_WARNING("CompositeEngine: unknown geometry types encountered during import.\n");
  }
  
    // get points from free curves
  for( i = curves.size(); i--; )
  {
    Curve* ptr = curves.get_and_step();
    temp_points.clean_out();
    ptr->points( temp_points );
    points += temp_points;
  }
  
    // get points and curves from free surfaces
  for( i = surfaces.size(); i--; )
  {
    Surface* ptr = surfaces.get_and_step();
    
    temp_points.clean_out();
    temp_curves.clean_out();
    
    ptr->points( temp_points );
    ptr->curves( temp_curves );
    
    points += temp_points;
    curves += temp_curves;
  }
  
    // get child geometry from free lumps
    // (are there ever any free lumps?)
  for( i = lumps.size(); i--; )
  {
    Lump* ptr = lumps.get_and_step();
    
    temp_points.clean_out();
    temp_curves.clean_out();
    temp_surfaces.clean_out();
    
    ptr->points( temp_points );
    ptr->curves( temp_curves );
    ptr->surfaces( temp_surfaces );
    
    points += temp_points;
    curves += temp_curves;
    surfaces += temp_surfaces;
  }
  
    // get child geometry from bodies
  for( i = bodies.size(); i--; )
  {
    BodySM* ptr = bodies.get_and_step();
    
    temp_points.clean_out();
    temp_curves.clean_out();
    temp_surfaces.clean_out();
    temp_lumps.clean_out();
    
    ptr->points( temp_points );
    ptr->curves( temp_curves );
    ptr->surfaces( temp_surfaces );
    ptr->lumps( temp_lumps );
    
    points += temp_points;
    curves += temp_curves;
    surfaces += temp_surfaces;
    lumps += temp_lumps;
  }

 
  // assert that lists do not contain duplicates
  //assert( (i = points.size(),   points.uniquify_ordered(),   i == points.size()  ) );
  //assert( (i = curves.size(),   curves.uniquify_ordered(),   i == curves.size()  ) );  
  points.uniquify_ordered(); // can have duplicates in some non-manifold cases
  curves.uniquify_ordered(); // can have duplicates in some non-manifold cases
  assert( (i = surfaces.size(), surfaces.uniquify_ordered(), i == surfaces.size()) );
  assert( (i = lumps.size(),    lumps.uniquify_ordered(),    i == lumps.size()   ) );
  
  if ( CGMApp::instance()->attrib_manager()->auto_actuate_flag(CA_COMPOSITE_VG) &&
       CGMApp::instance()->attrib_manager()->auto_read_flag(CA_COMPOSITE_VG) )
  {

    // Loop through the surfaces and if there is a 
    // TOPLOGY_BRIDGE_ID attribute get the id and set
    // it in the toplogy bridge for the surface.
    surfaces.reset();
    for(i=surfaces.size(); i--;)
    {
      Surface *s = surfaces.get_and_step();
      GeometryEntity *ge = dynamic_cast<GeometryEntity*>(s);
      TopologyBridge *tb = dynamic_cast<TopologyBridge*>(s);
      if(ge && tb)
      {
        DLIList<CubitSimpleAttrib> list;
        ge->get_simple_attribute("TOPOLOGY_BRIDGE_ID",list);
        list.reset();
        for(int j = list.size(); j--;)
        {
          const CubitSimpleAttrib& attrib = list.get_and_step();
          assert(attrib.int_data_list().size() == 1);
          int id = attrib.int_data_list()[0];
          ge->set_saved_id(id);
          
          std::vector<CubitString> names = attrib.string_data_list();
          if(!names.empty())
            names.erase(names.begin());

          ge->set_saved_names( names );

          tb->remove_simple_attribute_virt(attrib);
        }
      }
    }

    // Need to do curves here too.
    curves.reset();
    for(i=curves.size(); i--;)
    {
      Curve *c = curves.get_and_step();
      GeometryEntity *ge = dynamic_cast<GeometryEntity*>(c);
      TopologyBridge *tb = dynamic_cast<TopologyBridge*>(c);
      if(ge && tb)
      {
        DLIList<CubitSimpleAttrib> list;
        ge->get_simple_attribute("TOPOLOGY_BRIDGE_ID",list);
        list.reset();
        for(int j = list.size(); j--;)
        {
          const CubitSimpleAttrib& attrib = list.get_and_step();
          assert(attrib.int_data_list().size() == 1);
          int id = attrib.int_data_list()[0];
          ge->set_saved_id(id);
          
          std::vector<CubitString> names = attrib.string_data_list();
          if(!names.empty())
            names.erase(names.begin());
          ge->set_saved_names( names );

          tb->remove_simple_attribute_virt(attrib);
        }
      }
    }
    
    // create composites top-down
    if( !create_composites( bodies   ) ) result = CUBIT_FAILURE;
    if( !create_composites( surfaces ) ) result = CUBIT_FAILURE;
    if( !create_composites( curves   ) ) result = CUBIT_FAILURE;
    if( !create_composites( points   ) ) result = CUBIT_FAILURE;   
    
    
    for( i = bodies.size(); i--; )
    {
      BodySM* ptr = bodies.get_and_step();
      
      temp_curves.clean_out();
      temp_surfaces.clean_out();
      
      ptr->curves( temp_curves );
      ptr->surfaces( temp_surfaces );      
   
      //for each curve or surface, if it is not a composite,
      //and has a COMPOSIT_DATA_ATTRIB_NAME, 
      //which is an ENTITY_NAME, convert it into an ENTITY_NAME attribute

      for( int k=temp_curves.size(); k--; )
      {
        Curve *tmp_curve =  temp_curves.get_and_step();

        DLIList<CubitSimpleAttrib> list;
        tmp_curve->get_simple_attribute("COMPOSITE_ATTRIB",list);
        for( int j = list.size(); j--; )
        {
          CubitSimpleAttrib tmp_attrib = list.get_and_step();
          std::vector<CubitString> string_list = tmp_attrib.string_data_list();
          if( string_list[1] == "ENTITY_NAME" && !dynamic_cast<CompositeCurve*>( tmp_curve ) )
          {         
            //convert the attribute into an ENTITY_NAME attribute                        
            tmp_curve->remove_simple_attribute_virt( tmp_attrib );
            tmp_attrib.string_data_list().erase( tmp_attrib.string_data_list().begin() );            
            tmp_curve->append_simple_attribute_virt( tmp_attrib );   
          }
        }                
      }
      
      for( int k=temp_surfaces.size(); k--; )
      {
        Surface *tmp_surf = temp_surfaces.get_and_step();

        DLIList<CubitSimpleAttrib> list;
        tmp_surf->get_simple_attribute("COMPOSITE_ATTRIB",list);
        for( int j = list.size(); j--; )
        {
          CubitSimpleAttrib tmp_attrib = list.get_and_step();
          std::vector<CubitString> string_list = tmp_attrib.string_data_list();
          if( string_list[1] == "ENTITY_NAME" && !dynamic_cast<CompositeSurface*>( tmp_surf ) )
          {       
            //convert the attribute into an ENTITY_NAME attribute            
            tmp_surf->remove_simple_attribute_virt( tmp_attrib );
            tmp_attrib.string_data_list().erase( tmp_attrib.string_data_list().begin() );
            tmp_surf->append_simple_attribute_virt( tmp_attrib );
          }
        }                
      }
    }   
  }

  for ( i = bodies.size(); i--; )
    strip_attributes( bodies.get_and_step() );
  for ( i = surfaces.size(); i--; )
    strip_attributes( surfaces.get_and_step() );
  for ( i = curves.size(); i--; )
    strip_attributes( curves.get_and_step() );
  for ( i = points.size(); i--; )
    strip_attributes( points.get_and_step() );
  
    // update imported_geom list to contain composites rather than
    // entities used to create the composites
  for( i = imported_geom.size(); i--; )
  {
    TopologyBridge* bridge = imported_geom.step_and_get();
    TopologyBridge* virt = dynamic_cast<TopologyBridge*>(bridge->owner());
    if( virt )
      imported_geom.change_to( virt );
  }
    
/* DON'T DO THIS.  IT CHANGES THE ORDER OF THE IMPORT LIST,
   CAUSING IDS TO CHANGE.  FURTHER, IT RANDOMIZES THE IDS. 
   
    // composites were put in the list once for each underlying
    // entity.  need to remove the duplicates.
  imported_geom.uniquify_unordered();
*/
    // that's all folks
  return result;
} 

bool CompositeEngine::is_composite(TBOwner *bridge_owner)
{
  bool ret = false;
  if(bridge_owner)
  {
    if(dynamic_cast<CompositeBody*>(bridge_owner) ||
       dynamic_cast<CompositeLump*>(bridge_owner) ||
       dynamic_cast<CompositeSurface*>(bridge_owner) ||
       dynamic_cast<CompositeCurve*>(bridge_owner) ||
       dynamic_cast<CompositeCoEdge*>(bridge_owner) ||
       dynamic_cast<CompositePoint*>(bridge_owner))
    {
      ret = true;
    }
  }
  return ret;
}

bool CompositeEngine::is_composite(TopologyBridge *bridge)
{
  bool ret = false;
  if(bridge)
  {
    if(dynamic_cast<CompositeBody*>(bridge) ||
       dynamic_cast<CompositeLump*>(bridge) ||
       dynamic_cast<CompositeSurface*>(bridge) ||
       dynamic_cast<CompositeCurve*>(bridge) ||
       dynamic_cast<CompositeCoEdge*>(bridge) ||
       dynamic_cast<CompositePoint*>(bridge))
    {
      ret = true;
    }
  }
  return ret;
}

bool CompositeEngine::is_partition(TBOwner *bridge_owner)
{
  return false;
}

void CompositeEngine::remove_attributes( DLIList<TopologyBridge*> &bridge_list )
{
  //remove attributes off of underlying/hidden entities
  DLIList<Curve*> curve_list, temp_curves;
  DLIList<Surface*> surface_list, temp_surfaces;
  DLIList<Lump*> lump_list;
  DLIList<BodySM*> body_list;

  CAST_LIST( bridge_list, curve_list, Curve );
  CAST_LIST( bridge_list, surface_list, Surface );
  CAST_LIST( bridge_list, lump_list, Lump );
  CAST_LIST( bridge_list, body_list, BodySM );

  int i;
  for( i = surface_list.size(); i--; )
  {
    temp_curves.clean_out();
    surface_list.get_and_step()->curves( temp_curves );
    curve_list += temp_curves;
  }

  for( i = lump_list.size(); i--; )
  {
    Lump* lump = lump_list.get_and_step();

    temp_curves.clean_out();
    lump->curves( temp_curves );
    curve_list += temp_curves;

    temp_surfaces.clean_out();
    lump->surfaces( temp_surfaces );
    surface_list += temp_surfaces;
  }
  for( i = body_list.size(); i--; )
  {
    BodySM* body = body_list.get_and_step();

    temp_curves.clean_out();
    body->curves( temp_curves );
    curve_list += temp_curves;

    temp_surfaces.clean_out();
    body->surfaces( temp_surfaces );
    surface_list += temp_surfaces;
  }

//  DLIList<CompositeCurve*> ccurve_list;
//  DLIList<CompositeSurface*> csurf_list;

//  CAST_LIST( curve_list, ccurve_list, CompositeCurve );
//  CAST_LIST( surface_list, csurf_list, CompositeSurface );

//  ccurve_list.uniquify_unordered();
//  csurf_list.uniquify_unordered();
  curve_list.uniquify_unordered();
  surface_list.uniquify_unordered();

  int j,k;
  for( i = curve_list.size(); i--; )
  {
    Curve *tmp_curve = curve_list.get_and_step();
    strip_attributes(tmp_curve);
    CompositeCurve* tmp_comp_curve = dynamic_cast<CompositeCurve*>(tmp_curve);
    if(tmp_comp_curve)
    {
      for( j = 0; j < tmp_comp_curve->num_curves(); j++ )
      {
        strip_attributes( tmp_comp_curve->get_curve(j) );

        //remove attributes off underlying points too
        DLIList<TBPoint*> hidden_points;
        tmp_comp_curve->get_hidden_points( hidden_points );
        for( k=hidden_points.size(); k--; )
          strip_attributes( hidden_points.get_and_step() );
      }
    }
  }
  for( i = surface_list.size(); i--; )
  {
    Surface *tmp_surf = surface_list.get_and_step();
    strip_attributes(tmp_surf);
    CompositeSurface *tmp_comp_surf = dynamic_cast<CompositeSurface*>(tmp_surf);
    if(tmp_comp_surf)
    {
      for( int j = 0; j < tmp_comp_surf->num_surfs(); j++ )
        strip_attributes( tmp_comp_surf->get_surface(j) );
    }
  }

  //remove attrigutes off of bridges passed in
  for( i=bridge_list.size(); i--; )
    strip_attributes( bridge_list.get_and_step() );
}

void CompositeEngine::remove_attributes_from_unmodifed_virtual( DLIList<TopologyBridge*> &bridge_list )
{
  //remove attributes off of underlying/hidden entities
  DLIList<Curve*> curve_list, temp_curves;
  DLIList<Surface*> surface_list, temp_surfaces;
  DLIList<Lump*> lump_list;
  DLIList<BodySM*> body_list;

  CAST_LIST( bridge_list, curve_list, Curve );
  CAST_LIST( bridge_list, surface_list, Surface );
  CAST_LIST( bridge_list, lump_list, Lump );
  CAST_LIST( bridge_list, body_list, BodySM );

  int i;
  for( i = surface_list.size(); i--; )
  {
    temp_curves.clean_out();
    surface_list.get_and_step()->curves( temp_curves );
    curve_list += temp_curves;
  }

  for( i = lump_list.size(); i--; )
  {
    Lump* lump = lump_list.get_and_step();

    temp_curves.clean_out();
    lump->curves( temp_curves );
    curve_list += temp_curves;

    temp_surfaces.clean_out();
    lump->surfaces( temp_surfaces );
    surface_list += temp_surfaces;
  }
  for( i = body_list.size(); i--; )
  {
    BodySM* body = body_list.get_and_step();

    temp_curves.clean_out();
    body->curves( temp_curves );
    curve_list += temp_curves;

    temp_surfaces.clean_out();
    body->surfaces( temp_surfaces );
    surface_list += temp_surfaces;
  }

  DLIList<CompositeCurve*> ccurve_list;
  DLIList<CompositeSurface*> csurf_list;

  CAST_LIST( curve_list, ccurve_list, CompositeCurve );
  CAST_LIST( surface_list, csurf_list, CompositeSurface );

  ccurve_list.uniquify_unordered();
  csurf_list.uniquify_unordered();

  int k;
  for( i = ccurve_list.size(); i--; )
  {
    CompositeCurve* tmp_comp_curve = ccurve_list.get_and_step();
    /*
    for( j = 0; j < tmp_comp_curve->num_curves(); j++ )
    {
//      strip_attributes( tmp_comp_curve->get_curve(j) );

      //remove attributes off underlying points too
    }
    */
    DLIList<TBPoint*> hidden_points;
    tmp_comp_curve->get_hidden_points( hidden_points );
    for( k=hidden_points.size(); k--; )
      strip_attributes( hidden_points.get_and_step() );
  }
  for( i = csurf_list.size(); i--; )
  {
    CompositeSurface *tmp_comp_surf = csurf_list.get_and_step();
    /*
    for( int j = 0; j < tmp_comp_surf->num_surfs(); j++ )
      strip_attributes( tmp_comp_surf->get_surface(j) );
      */
    DLIList<Curve*> hidden_curves;
    tmp_comp_surf->get_hidden_curves( hidden_curves );
    for( k=hidden_curves.size(); k--; )
      strip_attributes( hidden_curves.get_and_step() );
  }

  //remove attrigutes off of bridges passed in
  /*
  for( i=bridge_list.size(); i--; )
    strip_attributes( bridge_list.get_and_step() );
    */
}

void CompositeEngine::strip_attributes( TopologyBridge* bridge )
{
  const char* const attrib_names[] = { "COMPOSITE_GEOM",
                                       "COMPOSITE_STITCH",
                                       "COMPOSITE_SENSE",
                                       "COMPOSITE_ATTRIB",
                                       "COMPOSITE_NULLGEOM",
                                       "TOPOLOGY_BRIDGE_ID",
                                       0 };
  
  DLIList<CubitSimpleAttrib> list;
  for( int i = 0; attrib_names[i]; i++ )
  {
    bridge->get_simple_attribute( attrib_names[i], list );
    while( list.size() )
    {
      bridge->remove_simple_attribute_virt(list.pop());
    }
  }
}

      
//-------------------------------------------------------------------------
// Purpose       : Retreive the CSA with the specified name if it exists
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/18/02
//-------------------------------------------------------------------------
CubitSimpleAttrib
CompositeEngine::find_attribute_by_name( TopologyBridge* bridge, 
                                         const CubitString name )
{
  CubitSimpleAttrib result;
  DLIList<CubitSimpleAttrib> attrib_list;
  bridge->get_simple_attribute( name, attrib_list );
  attrib_list.reset();
  if ( attrib_list.size() )
    result = attrib_list.extract();
    
  return result;
}  
      
//-------------------------------------------------------------------------
// Purpose       : Create composite surfaces
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/18/02
//-------------------------------------------------------------------------
CubitStatus CompositeEngine::create_composites( DLIList<Curve*>& list )
{
  CubitStatus result = CUBIT_FAILURE;
  
  for( int i = list.size(); i--; )
  {
    Curve* curve = list.get_and_step();
    CompositeCurve* ccurve;
    ccurve = dynamic_cast<CompositeCurve*>(curve->owner());
    if( ccurve )
      curve = ccurve;
    CubitSimpleAttrib attrib = find_attribute_by_name( curve, "COMPOSITE_GEOM" );
    if( !attrib.isEmpty() )
    {
      curve->remove_simple_attribute_virt( attrib );
      
/*
      Curve* stitch_partner = 0;
      if( attrib->int_data_list()->size() )
      {
        assert(attrib->int_data_list()->size() == 1 );
        attrib->int_data_list()->reset();
        int stitch_partner_id = *attrib->int_data_list()->get();
        TopologyBridge* bridge = stitchMap[stitch_partner_id];
        if( !bridge )
        {
          stitchMap[stitch_partner_id] = curve;
          continue;
        }
        
        stitch_partner = dynamic_cast<Curve*>(bridge);
        assert( !!stitch_partner );
      }
      
      CompositeCurve* ccurve;
      ccurve = dynamic_cast<CompositeCurve*>(curve->owner());
      if( ccurve )
        curve = ccurve;
      if( stitch_partner )
      {
        ccurve = dynamic_cast<CompositeCurve*>(stitch_partner->owner());
        if( ccurve )
          stitch_partner = ccurve;
      }
      
      CompositeSurface* surf = remove_curve( curve, stitch_partner );
      if( !surf )
      {
        PRINT_ERROR("Error creating composite surface.  Could not remove curve.\n");
        result = CUBIT_FAILURE;
      }
*/
      CompositeSurface* surf = remove_curve( curve );
      if( !surf )
      {
        PRINT_ERROR("Error creating composite surface.  Could not remove curve.\n");
        result = CUBIT_FAILURE;
      }
      else
      {
        // Tell the composite surface which surfaces to ignore during
        // evaluation.
        for(int j=0; j<surf->num_surfs(); j++)
        {
          Surface *srf = surf->get_surface(j);
          CubitSimpleAttrib ignore_attrib = find_attribute_by_name( srf, "COMPOSITE_IGNORE" );
          if( !ignore_attrib.isEmpty() )
          {
            surf->ignore_surface(srf);
            srf->remove_simple_attribute_virt( ignore_attrib );
          }
        }
        surf->read_attributes();
      }
    }
  }
  return result;
}
      
  
//-------------------------------------------------------------------------
// Purpose       : Create composite curves
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/18/02
//-------------------------------------------------------------------------
CubitStatus CompositeEngine::create_composites( DLIList<TBPoint*>& list )
{
  CubitStatus result = CUBIT_FAILURE;
  
    // Restore composite curves.
  int i;
  for( i = list.size(); i--; )
  {
    TBPoint* point = list.get_and_step();
    CompositePoint* cpoint;
    cpoint = dynamic_cast<CompositePoint*>(point->owner());
    if( cpoint )
      point = cpoint;
    CubitSimpleAttrib attrib = find_attribute_by_name( point, "COMPOSITE_GEOM" );
    if( !attrib.isEmpty() )
    {
      point->remove_simple_attribute_virt( attrib );
/*      
      TBPoint* stitch_partner = 0;
      if( attrib->int_data_list()->size() )
      {
        assert(attrib->int_data_list()->size() == 1 );
        attrib->int_data_list()->reset();
        int stitch_partner_id = *attrib->int_data_list()->get();
        TopologyBridge* bridge = stitchMap[stitch_partner_id];
        if( !bridge )
        {
          stitchMap[stitch_partner_id] = point;
          continue;
        }
        
        stitch_partner = dynamic_cast<TBPoint*>(bridge);
        assert( !!stitch_partner );
      }
*/      
/*        
      if( stitch_partner )
      {
        cpoint = dynamic_cast<CompositePoint*>(stitch_partner->owner());
        if( cpoint )
          stitch_partner = cpoint;
      }
*/        
      CompositeCurve* curve = remove_point( point /*, stitch_partner*/ );
      if( !curve )
      {
        PRINT_ERROR("Error creating composite curve.  Could not remove point.\n");
        result = CUBIT_FAILURE;
      }
      else
      {
        curve->read_attributes();
      }
    }
  }
  
    // Restore point curves
  for( i = list.size(); i--; )
  {
    TBPoint* point = list.get_and_step();
    CompositePoint* cpoint;
    cpoint = dynamic_cast<CompositePoint*>(point->owner());
    if( cpoint )
      point = cpoint;
    CubitSimpleAttrib attrib = find_attribute_by_name( point, "COMPOSITE_NULLGEOM" );
    if( !attrib.isEmpty() )
    {
      point->remove_simple_attribute_virt( attrib );
   
      CompositeCurve* curve = restore_point_in_surface( point );
      if( !curve )
      {
        PRINT_ERROR("Error creating composite curve.  Could not rrestore point-curve.\n");
        result = CUBIT_FAILURE;
      }
      else
      {
        curve->read_attributes();
      }
    }
  }
  
  return result;
}
      

//-------------------------------------------------------------------------
// Purpose       : Create composite lumps
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/18/02
//-------------------------------------------------------------------------
CubitStatus CompositeEngine::create_composites( DLIList<Surface*>&  )
  { return CUBIT_SUCCESS; }
  
//-------------------------------------------------------------------------
// Purpose       : Create composite bodies
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/18/02
//-------------------------------------------------------------------------
CubitStatus CompositeEngine::create_composites( DLIList<BodySM*>&  )
{ 
  return CUBIT_SUCCESS; 
}

//-------------------------------------------------------------------------
// Purpose       : Save composite geometry
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/18/02
//-------------------------------------------------------------------------
CubitStatus CompositeEngine::export_geometry( DLIList<TopologyBridge*>& list )
{
  int i;
  CubitStatus result = CUBIT_SUCCESS;

  if ( CGMApp::instance()->attrib_manager()->auto_update_flag(CA_COMPOSITE_VG) &&
       CGMApp::instance()->attrib_manager()->auto_write_flag(CA_COMPOSITE_VG) )
  {

    DLIList<Curve*> curve_list, temp_curves;
    DLIList<Surface*> surface_list, temp_surfaces;
    DLIList<Lump*> lump_list, temp_lumps;
    DLIList<BodySM*> body_list;

    CAST_LIST( list, curve_list, Curve );
    CAST_LIST( list, surface_list, Surface );
    CAST_LIST( list, lump_list, Lump );
    CAST_LIST( list, body_list, BodySM );

    for( i = surface_list.size(); i--; )
    {
      Surface* surf = surface_list.get_and_step();
      temp_curves.clean_out();
      surf->curves( temp_curves );
      curve_list += temp_curves;
    }
    for( i = lump_list.size(); i--; )
    {
      Lump* lump = lump_list.get_and_step();

      temp_curves.clean_out();
      lump->curves( temp_curves );
      curve_list += temp_curves;

      temp_surfaces.clean_out();
      lump->surfaces( temp_surfaces );
      surface_list += temp_surfaces;
    }
    for( i = body_list.size(); i--; )
    {
      BodySM* body = body_list.get_and_step();

      temp_curves.clean_out();
      body->curves( temp_curves );
      curve_list += temp_curves;

      temp_surfaces.clean_out();
      body->surfaces( temp_surfaces );
      surface_list += temp_surfaces;

      temp_lumps.clean_out();
      body->lumps( temp_lumps );
      lump_list += temp_lumps;
    }

    DLIList<CompositeCurve*> ccurve_list;
    DLIList<CompositeSurface*> csurf_list;
    DLIList<CompositeLump*> clump_list;
    DLIList<CompositeBody*> cbody_list;

    CAST_LIST( curve_list, ccurve_list, CompositeCurve );
    CAST_LIST( surface_list, csurf_list, CompositeSurface );
  //  CAST_LIST( lump_list, clump_list, CompositeLump );
  //  CAST_LIST( body_list, cbody_list, CompositeBody );

    ccurve_list.uniquify_unordered();
    csurf_list.uniquify_unordered();


   // Add an attribute to each topology bridge specifying what its id is.  Upon restoring
   // the file these attributes will be read so that the topology id can be stored with the
   // topology bridge.  Normally, the ids would get set when creating the associated RefEntity
   // and the topology bridge id would be set to be that of the RefEntity.  However, for the
   // underlying entities of a composite RefEntities are not created so the ids never get
   // set when restoring.  Thus we are doing it manually here.  Only doing curves and
   // surfaces right now.
   csurf_list.reset();
   for(i=csurf_list.size(); i--;)
   {
     CompositeSurface *cs = csurf_list.get_and_step();
     int number_surfs = cs->num_surfs();
     for(int j=0; j<number_surfs; ++j)
     {
       Surface *s = cs->get_surface(j);
       GeometryEntity *ge = dynamic_cast<GeometryEntity*>(s);
       TopologyBridge *tb = dynamic_cast<TopologyBridge*>(s);
       if(ge && tb)
       {
         std::vector<int> int_list;
         int_list.push_back(ge->get_saved_id() );

         //save out the names
         std::vector<CubitString> names;
         ge->get_saved_names( names );
         names.insert(names.begin(), CubitString("TOPOLOGY_BRIDGE_ID"));
         CubitSimpleAttrib geom_attrib( &names, 0, &int_list );
         append_attrib(tb, geom_attrib);
       }
     }
   }

   ccurve_list.reset();
   for(i=ccurve_list.size(); i--;)
   {
     CompositeCurve *cc = ccurve_list.get_and_step();
     int number_curves = cc->num_curves();
     for(int j=0; j<number_curves; ++j)
     {
       Curve *c = cc->get_curve(j);
       GeometryEntity *ge = dynamic_cast<GeometryEntity*>(c);
       TopologyBridge *tb = dynamic_cast<TopologyBridge*>(c);
       if(ge && tb)
       {
         std::vector<int> int_list;
         int_list.push_back( ge->get_saved_id() );

         //save out the names
         std::vector<CubitString> names;
         ge->get_saved_names( names );
         names.insert(names.begin(), CubitString("TOPOLOGY_BRIDGE_ID"));

         CubitSimpleAttrib geom_attrib( &names, 0, &int_list );
         append_attrib(tb, geom_attrib);
       }
     }
   }

    for( i = ccurve_list.size(); i--; )
      save( ccurve_list.get_and_step() );
    for( i = csurf_list.size(); i--; )
      save( csurf_list.get_and_step() );
    for( i = clump_list.size(); i--; )
      save( clump_list.get_and_step() );
    for( i = cbody_list.size(); i--; )
      save( cbody_list.get_and_step() );
  }
  
  list.reset();
  DLIList<TopologyBridge*> underlying;
  for( i = list.size(); i--; )
  {
    TopologyBridge* ptr = list.step_and_get();
    underlying.clean_out();
    if( CompositeCurve* curve = dynamic_cast<CompositeCurve*>(ptr) )
      for( int j = 0; j < curve->num_curves(); j++ )
        underlying.append( curve->get_curve(j) );
    else if( CompositeSurface* surf = dynamic_cast<CompositeSurface*>(ptr) )
      for( int j = 0; j < surf->num_surfs(); j++ )
        underlying.append( surf->get_surface(j) );
/*
    else if( CompositeLump* lump = dynamic_cast<CompositeLump*>)(ptr) )
      for( int j = 0; j < lump->num_entities(); j++ )
        underlying.append( lump->lump(j) );
    else if( CompositeBody* body = dynamic_cast<CompositeBody*>)(ptr) )
      for( int j = 0; j < body->num_entities(); j++ )
        underlying.append( body->body(j) );
*/
    if( underlying.size() )
    {
      list.change_to( underlying.pop() );
      list += underlying;
    }
  }
  
  return result;
}


void CompositeEngine::append_attrib( TopologyBridge* tb, 
                                     const CubitSimpleAttrib& csa )
{
  CubitSimpleAttrib old = find_attribute_by_name( tb,
                           csa.character_type() );
  if( !old.isEmpty() )
  {
    tb->remove_simple_attribute_virt( old );
  }
  tb->append_simple_attribute_virt( csa );
}


CubitStatus CompositeEngine::save( CompositePoint* point )
{
  DLIList<CompositePoint*> points;
  point->get_stitched( points );
  if( points.size() > 1 )
  {
    CubitString name("COMPOSITE_STITCH");
    std::vector<CubitString> string_list;
    string_list.push_back(name);
    CubitSimpleAttrib geom_attrib( &string_list, 0, 0 );
    int uid = TDUniqueId::generate_unique_id();
    geom_attrib.int_data_list().push_back( uid );
    
    for( int i = points.size(); i--; )
      append_attrib( points.step_and_get()->get_point(), geom_attrib );
  }
  return CUBIT_SUCCESS;
}
    

CubitStatus CompositeEngine::save( CompositeCurve* curve  )
{
  DLIList<TBPoint*> hidden_points;
  std::vector<CubitString> string_list;

  if (curve->num_curves() == 0) // point-curve
  {
    assert(hidden_points.size() == 0);
    CubitString ptname("COMPOSITE_NULLGEOM");
    string_list.clear();
    string_list.push_back( ptname );
    CubitSimpleAttrib null_geom_attrib( &string_list, 0, 0 );
    CompositePoint* pt = curve->start_point();
    assert(curve->end_point() == pt);
    append_attrib( pt, null_geom_attrib );
    curve->write_attributes();
    return CUBIT_SUCCESS;
  }
  
  ;
  string_list.push_back( CubitString("COMPOSITE_GEOM") );
  CubitSimpleAttrib geom_attrib( &string_list, 0, 0 );

  int i;

  curve->write_attributes();
  curve->get_hidden_points( hidden_points );
  for( i = hidden_points.size(); i--; )
  {
    TBPoint* point = hidden_points.get_and_step();
    if( CompositePoint* cpoint = dynamic_cast<CompositePoint*>(point) )
      save(cpoint);
    
    append_attrib( point, geom_attrib );
  }


  if( curve->is_stitched() && curve == curve->primary_stitched_curve() )
  {
    int stitch_uid = TDUniqueId::generate_unique_id();
    DLIList<CompositeCurve*> curve_list;
    curve->get_stitched( curve_list );

    string_list[0] = "COMPOSITE_STITCH";
    CubitSimpleAttrib stitch_attrib( &string_list, 0, 0 );
    stitch_attrib.int_data_list().push_back(stitch_uid);
    for( i = curve_list.size(); i--; )
      append_attrib( curve_list.step_and_get(), stitch_attrib );
    stitch_attrib.int_data_list().clear();

    if( curve_list.move_to( curve ) )
      curve_list.extract();

    for( i = curve_list.size(); i--; )
      save( curve_list.step_and_get() );
  }

  
  string_list[0] = "COMPOSITE_SENSE";
  CubitSimpleAttrib sense_attrib( &string_list, 0, 0 );
  
  for( i = 0; i < curve->num_curves(); i++ )
  {
    CubitSimpleAttrib old = find_attribute_by_name( curve->get_curve(i), "COMPOSITE_SENSE" );
    if( !old.isEmpty() )
    {
      if( curve->get_sense(i) == CUBIT_FORWARD )
        curve->get_curve(i)->remove_simple_attribute_virt( old );
    }
    else if( curve->get_sense(i) == CUBIT_REVERSED )
    {
      curve->get_curve(i)->append_simple_attribute_virt( sense_attrib );
    }
  }

  return CUBIT_SUCCESS;
}

   
CubitStatus CompositeEngine::save( CompositeSurface* surf )
{ 
  DLIList<Curve*> comp_curves;
  surf->hidden_entities().hidden_curves( comp_curves );
  surf->write_attributes();
  
  std::vector<CubitString> string_list;
  string_list.push_back( CubitString("COMPOSITE_GEOM") );
  CubitSimpleAttrib geom_attrib( &string_list, 0, 0 );

  int i, j;

  for( i = comp_curves.size(); i--; )
  {
    CompositeCurve* ccurve = 
      dynamic_cast<CompositeCurve*>(comp_curves.get_and_step());
    assert(ccurve!= NULL );
    
    for( j = 0; j < ccurve->num_curves(); j++ )
    {      
      Curve* tb = ccurve->get_curve(j);
      append_attrib( tb, geom_attrib );
    }  
  }
  

  if( surf->get_stitch_partner() )
  {
    CubitSimpleAttrib old =
      find_attribute_by_name( surf->get_stitch_partner(), "COMPOSITE_STITCH" );
    if( !old.isEmpty() )
    {
      append_attrib( surf, old );
    }
    else
    {
      int stitch_uid = TDUniqueId::generate_unique_id();
      string_list[0] = "COMPOSITE_STITCH";
      CubitSimpleAttrib stitch_attrib( &string_list, 0, 0);
      stitch_attrib.int_data_list().push_back( stitch_uid );
      append_attrib( surf, stitch_attrib );
      stitch_attrib.int_data_list().clear();
    }
  }

  
  string_list[0] = "COMPOSITE_SENSE";
  CubitSimpleAttrib sense_attrib( &string_list, 0, 0 );
  
  for( i = 0; i < surf->num_surfs(); i++ )
  {
    CubitSimpleAttrib old = find_attribute_by_name( surf->get_surface(i), "COMPOSITE_SENSE" );
    if( !old.isEmpty() )
    {
      if( surf->get_sense(i) == CUBIT_FORWARD )
        surf->get_surface(i)->remove_simple_attribute_virt( old );
    }
    else if( surf->get_sense(i) == CUBIT_REVERSED )
    {
      surf->get_surface(i)->append_simple_attribute_virt( sense_attrib );
    }
  }

  string_list[0] = "COMPOSITE_IGNORE";
  CubitSimpleAttrib ignore_attrib( &string_list, 0, 0 );

  DLIList<Surface*> srfs;
  surf->get_ignored_surfs(srfs);
  for( i = 0; i < surf->num_surfs(); i++ )
  {
    Surface *srf = surf->get_surface(i);
    if(srfs.is_in_list(srf))
      srf->append_simple_attribute_virt(ignore_attrib);
  }


  return CUBIT_SUCCESS;
}

CubitStatus CompositeEngine::save( CompositeLump* lump )
{
  DLIList<Surface*> comp_surfs;
  lump->hidden_entities().hidden_surfaces( comp_surfs );
  
  CubitString name("COMPOSITE_GEOM");
  std::vector<CubitString> string_list;
  string_list.push_back( name );
  CubitSimpleAttrib geom_attrib( &string_list, 0, 0 );

  int i, j;
  
  for( i = comp_surfs.size(); i--; )
  {
    CompositeSurface* csurf = 
      dynamic_cast<CompositeSurface*>(comp_surfs.get_and_step());
    assert(csurf!= NULL );
    
    for( j = 0; j < csurf->num_surfs(); j++ )
    {
      append_attrib( csurf, geom_attrib );
      geom_attrib.int_data_list().clear();
    }  
  }
  
  return CUBIT_SUCCESS;
}

CubitStatus CompositeEngine::save( CompositeBody* )
{ return CUBIT_SUCCESS; }


CompositePoint* CompositeEngine::stitch_points( TBPoint* pt1, TBPoint* pt2 )
{
  CompositePoint* cp1 = dynamic_cast<CompositePoint*>(pt1);
  CompositePoint* cp2 = dynamic_cast<CompositePoint*>(pt2);
  if( !cp1 ) cp1 = replace_point( pt1 );
  if( !cp2 ) cp2 = replace_point( pt2 );
  
  while( CompositeCurve* curve = cp2->next_curve() )
  {
    if( curve->start_point() == cp2 )
      curve->start_point( cp1 );
    if( curve->end_point() == cp2 )
      curve->end_point( cp1 );
  }
  
  cp1->stitch( cp2 );
  return cp1;
}

CompositeCurve* CompositeEngine::stitch_curves( Curve* curve1, Curve* curve2 )
{
  CompositeCurve* c1 = dynamic_cast<CompositeCurve*>(curve1);
  CompositeCurve* c2 = dynamic_cast<CompositeCurve*>(curve2);
  if( !c1 ) c1 = replace_curve( curve1 );
  if( !c2 ) c2 = replace_curve( curve2 );
  
  bool reversed = false, fail = false;
  if( c1->start_point() == c2->end_point() )
    reversed = true;
  else if( c1->start_point() != c2->start_point() )
    fail = true;
    
  if (c1->start_point() == c1->end_point())
  {
    CubitVector pt, tan1, tan2, junk;
    c1->position_from_fraction( 0.5, pt );
    c1->closest_point( pt, junk, &tan1 );
    c2->closest_point( pt, junk, &tan2 );
    reversed = (tan1 % tan2) < 0.0;
  }
  
  if( fail || 
      (reversed && c1->end_point() != c2->start_point()) ||
      (!reversed && c1->end_point() != c2->end_point()) )
  {
    if( c1->num_curves() == 1 )
      remove_composite( c1 );
    if( c2->num_curves() == 1 )
      remove_composite( c2 );
    return 0;
  }
  
  //while( CompositeCoEdge* coedge = c2->first_coedge() )
  //{
    //c2->remove( coedge );
    //c1->add( coedge );
    //if( reversed )
   //   coedge->reverse();
  //}
  CompositeCoEdge* coedge = 0;
  if (reversed) while ((coedge = c2->next_coedge( coedge )))
    coedge->reverse();
    

  c1->stitch( c2 );
  return c1;
}

CompositeSurface* 
CompositeEngine::stitch_surfaces( Surface* /*surf1*/, Surface* /*surf2*/ )
{
  return 0;
/*
  CompositeSurface* cs1 = dynamic_cast<CompositeSurface*>(surf1);
  CompositeSurface* cs2 = dynamic_cast<CompositeSurface*>(surf2);
  
    // Can only stitch pairs of surfaces. Fail if either 
    // surface is already stitched with another surface.
  if( (cs1 && cs1->get_stitch_partner()) ||
      (cs2 && cs2->get_stitch_partner())  )
    return 0;
  
  if( ! cs1 ) cs1 = replace_surface( surf1 );
  if( ! cs2 ) cs2 = replace_surface( surf2 );
  
    // special case - topological sphere
  if( !cs1->first_loop() && !cs2->first_loop() )
  {
    if( cs1->stitch( cs2 ) )
      return cs1;  // success
    
    if( !cs1->has_hidden_entities() )
      remove_composite( cs1 );
    if( !cs2->has_hidden_entities() )
      remove_composite( cs2 );
    return 0; // fail
  }
  
    // make sure all curves are merged
  CompositeLoop* loop = 0;
  bool reversed = false;
  bool failed = false;
  while( !failed && (loop = cs1->next_loop(loop)) )
  {
    CompositeCoEdge* coedge = loop->first_coedge();
    do
    {
      if( ! coedge->get_curve()->find_coedge( cs2 ) )
      {
        failed = true;
        break;
      }
      coedge = coedge->next();
    }
    while( coedge != loop->first_coedge() );
  }
  
  if( failed )
  {
    if( ! cs1->has_hidden_entities() )
      remove_composite( cs1 );
    if( ! cs2->has_hidden_entities() )
      remove_composite( cs2 );
  
    return 0;
  }
  
  while( CompositeCoSurf* cosurf = cs2->next_co_surface() )
  {
    cs2->remove( cosurf );
    cs1->add( cosurf );
    if( reversed )
    {
      CubitSense sense = cosurf->sense();
      sense = sense == CUBIT_REVERSED ? CUBIT_FORWARD : CUBIT_REVERSED;
      cosurf->sense( sense );
    }
  }
  cs1->stitch( cs2 );

  return cs1;
*/
}

TBPoint* CompositeEngine::insert_point(CompositeCurve* curve, double u)
{
  int i, index;
  double param;
  TBPoint *result = 0;
  CubitVector position;
  DLIList<TBPoint*> pts;
  double tolsqr = GEOMETRY_RESABS*GEOMETRY_RESABS;
  TBPoint *start, *end;
  DLIList<TopologyBridge*> pt_list;
  Curve *rcurve = 0;

  // Get the curve in the CompositeCurve on which the
  // parameter u lies.
  if( ! curve->curve_param( u, param, index ) )
    return 0;
  rcurve = curve->get_curve(index);

  // Get the position corresponding with u.
  if( !rcurve->position_from_u(param,position) )
    return 0;
  
  // Get the CompositeCurve's hidden points.  We will
  // first see if the position to insert is on top of 
  // one of the hidden points.
  curve->get_hidden_points(pts);
  for(i=pts.size(); i--;)
  {
    TBPoint *cur_pt = pts.get_and_step();

    // Don't use GEOMETRY_RESABS if we can get a
    // value from the solid modeling engine.
    double tmp_tol = tolsqr;
    CompositePoint *cp = dynamic_cast<CompositePoint*>(cur_pt);
    GeometryQueryEngine *gqe = NULL;
    if(cp)
    {
      TBPoint *real_pt = cp->get_point();
      gqe = real_pt->get_geometry_query_engine();
    }
    else
      gqe = cur_pt->get_geometry_query_engine();
    if(gqe)
    {
      double tmp_tol = gqe->get_sme_resabs_tolerance();
      tmp_tol *= tmp_tol;
      if(tmp_tol > tolsqr)
        tolsqr = tmp_tol;
    }
    if( (cur_pt->coordinates() - position).length_squared() < tmp_tol )
    {
      result = cur_pt;
      i = 0;
    }
  }

  // If the insert position wasn't one of the hidden 
  // points look at the ends of the underlying curves
  // to see if it is one of the end points.
  if(!result)
  {
    rcurve->get_children_virt(pt_list);
    assert(pt_list.size());
    pt_list.reset();
    start = dynamic_cast<TBPoint*>(pt_list.get());
    end = dynamic_cast<TBPoint*>(pt_list.next());
    assert( start && end );
    
    if( (start->coordinates() - position).length_squared() < tolsqr )
    {
      result = start;
    }
    else if( (end->coordinates() - position).length_squared() < tolsqr )
    {
      result = end;
    }
    else
    {
      // If the insert position wasn't one of the end points go 
      // ahead and insert a point into the curve.
      result = PartitionEngine::instance().insert_point(rcurve,param);
      if( !result )
        return 0;
    }
  }
    
  // Work with a real point here because functions called from
  // restore_point_in_curve() may remove the composite all together
  // leaving us with a stale CompositePoint pointer at this level.
  CompositePoint *comp_pt = dynamic_cast<CompositePoint*>(result);
  if(comp_pt)
    result = comp_pt->get_point();

  // Finally, undo the composite.
  if ( ! restore_point_in_curve( result ) )
    return 0;
    
  TBPoint* comp = dynamic_cast<CompositePoint*>(result->owner());
  return comp ? comp : result;
}

TBPoint* CompositeEngine::insert_point_curve( CompositeSurface* surf,
                                            const CubitVector& pos )
{
    // find closest surface
  int i = surf->closest_underlying_surface( pos );
  Surface *dummy_surf, *real_surf = surf->get_surface(i);
  
  TBPoint* pt = PartitionEngine::instance().insert_point_curve( real_surf, pos, dummy_surf );
  if ( !pt )
    return 0;
  
  if (!restore_point(pt))
    return 0;
  
  return dynamic_cast<CompositePoint*>(pt->owner());
}

int CompositeEngine::is_hidden(TopologyBridge *tb)
{
  int ret = 0;
  int done = 0;
  TopologyBridge *cur_bridge = tb;
  while(!done)
  {
    TBOwner *tbowner = cur_bridge->owner();
    if(dynamic_cast<BridgeManager*>(tbowner))
    {
      done = 1;
      ret = 0;
    }
    else if(dynamic_cast<HiddenEntitySet*>(tbowner))
    {
      done = 1;
      ret = 1;
    }
    cur_bridge = dynamic_cast<TopologyBridge*>(tbowner);
    if(!cur_bridge)
      done = 1;
  }
  return ret;
}

CubitStatus CompositeEngine::insert_curve( DLIList<Surface*>& surfaces,
                                           DLIList<CubitVector*>& polyline,
                                           DLIList<Surface*>& new_surfaces,
                                           DLIList<Curve*>& new_curves )
{
  int i;
  
  DLIList<Surface*> surfs_to_reverse;
  DLIList<Surface*> real_surfs;
  DLIList<TBPoint*> points;
  surfaces.reset();

  for ( i = surfaces.size(); i--; )
  {
    Surface* surf = surfaces.get_and_step();
    CompositeSurface* comp = dynamic_cast<CompositeSurface*>(surf);
    if ( !comp )
    {
      real_surfs.append(surf);
      continue;
    }
    
    // Partitioning needs all of the surfaces involved in the partition
    // to have facets oriented in the same direction w.r.t. the volume (in or out).
    // For surface facets that will have to be flipped we also want to reverse
    // the sense of the PartitionSurfaces to be consistent.  Therefore, we
    // create a list here of surfaces/facets that need to be flipped.  We
    // will use the sense info from the CompositeSurface to determine this.
    for ( int j = 0; j < comp->num_surfs(); j++ )
    {
      real_surfs.append( comp->get_surface( j ) );
      if(comp->get_sense(j) == CUBIT_REVERSED)
        surfs_to_reverse.append(comp->get_surface(j));
    }
      
      // append all real points visible on composite surfaces
      // to list.
    CompositeLoop* loop = 0;
    while( (loop = comp->next_loop(loop) ) != NULL )
    {
      CompositeCoEdge* coe = loop->first_coedge();
      do {
        points.append(coe->get_curve()->start_point()->get_point());
        points.append(coe->get_curve()->end_point()->get_point());
        coe = coe->next();
      } while( coe != loop->first_coedge() );
    }
  }
  
  CubitStatus s = PartitionEngine::instance().
    insert_curve( real_surfs, polyline, new_surfaces, new_curves, NULL,
    surfs_to_reverse.size() ? &surfs_to_reverse : NULL);
  if ( !s )
    return CUBIT_FAILURE;
  
  new_curves.last();
  for ( i = new_curves.size(); i--; )
  {
    Curve* curve = new_curves.step_and_get();
    CompositeCurve* ccurve = dynamic_cast<CompositeCurve*>(curve->owner());
    if ( !ccurve ) continue;
    
    restore_curve( ccurve );
    ccurve = dynamic_cast<CompositeCurve*>(curve->owner());
    if ( ccurve )
      new_curves.change_to(ccurve);
  }
  
    // The remainder of this function tries to composite over
    // 2-valence vertics in the chain of resulting curves.
    // These vertices are the result of intersections with 
    // hidden curves.  If there is only one result curve, then
    // skip the rest.
  if (new_curves.size() == 1)
  {
    clean_out_deactivated_geometry();
    return CUBIT_SUCCESS;
  }
  
  
    // Remove any 2-valence vertices in chain of resulting curves.
    // These are intersections with hidden curves.
  new_curves.reset();
  Curve* prev_curve = new_curves.get();
  DLIList<TopologyBridge*> prev_pts(2), next_pts(2), pt_curves;
  DLIList<TBPoint*> dead_points;
  prev_curve->get_children( prev_pts, false, COMPOSITE_LAYER );
  for ( i = new_curves.size(); i--; )
  {
    Curve* curve = new_curves.step_and_get();
    curve->get_children( next_pts, false, COMPOSITE_LAYER );
    while( prev_pts.size() )
    {
      TBPoint* point = dynamic_cast<TBPoint*>(prev_pts.pop());
      if( !next_pts.move_to( point ) )
        continue;
      next_pts.extract();
      
      pt_curves.clean_out();
      point->get_parents( pt_curves );
      if ( pt_curves.size() != 2 || points.is_in_list(point) )
        continue;
       
      if ( CompositePoint* comp = dynamic_cast<CompositePoint*>(point) )
        if ( points.is_in_list(comp->get_point()) )
          continue;
       
      dead_points.append( point );
      prev_pts.clean_out();
    } 
       
    prev_pts = next_pts;
    next_pts.clean_out();
  }

      
  while( dead_points.size() )
  {
    TBPoint* point = dead_points.pop();
    if ( CompositePoint* comp = dynamic_cast<CompositePoint*>(point->owner()) )
      point = comp;
      
    pt_curves.clean_out();
    point->get_parents( pt_curves );
    if ( pt_curves.size() == 1 )
      continue;
    assert(pt_curves.size() == 2);
    
    CompositeCurve* keep = remove_point(point);
    if (!keep)
      continue;
    
    bool add = true;
    if ( pt_curves.move_to(keep) )
    {
      pt_curves.extract();
      add = false;
    }
    
    while( pt_curves.size() )
    {
      Curve* stale_ptr = reinterpret_cast<Curve*>(pt_curves.pop());
      if ( new_curves.move_to( stale_ptr ) )
      {
        if ( add )
        {
          new_curves.change_to( keep );
          add = false;
        }
        else
        {
          new_curves.remove();
        }
      }
    }
    if ( add )
      new_curves.append( keep );
  }
 
  clean_out_deactivated_geometry();
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Update for destroted underlying topology
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/14/04
//-------------------------------------------------------------------------
void CompositeEngine::notify_deactivated (CompositeBody* body)
{
  while (CompositeLump* vol_ptr = body->next_lump(NULL))
  {
    body->remove( vol_ptr );
  }
  deactivatedList.append_unique( body );
}

//-------------------------------------------------------------------------
// Purpose       : Update for destroted underlying topology
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/14/04
//-------------------------------------------------------------------------
void CompositeEngine::notify_deactivated( CompositeLump* vol )
{
  while (CompositeShell* shell = vol->next_shell(NULL))
  {
    while (CompositeCoSurf* cosurf = shell->first_co_surf())
    {
      shell->remove( cosurf );
      CompositeSurface* surf = cosurf->get_surface();
      if (surf)
        surf->remove( cosurf );
      delete cosurf;
    }
    vol->remove( shell );
    deactivatedList.append_unique( shell );
  }
  deactivatedList.append_unique( vol );
}

  
//-------------------------------------------------------------------------
// Purpose       : Update for destroted underlying topology
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/31/03
//-------------------------------------------------------------------------
void CompositeEngine::notify_deactivated (CompositeSurface* surface)
{
  while (CompositeCoSurf* cosurf_ptr = surface->next_co_surface(NULL))
  {
    if (cosurf_ptr->get_shell())
      cosurf_ptr->get_shell()->remove(cosurf_ptr);
    surface->remove(cosurf_ptr);
    delete cosurf_ptr;
  }
  
  while (CompositeLoop* loop_ptr = surface->first_loop())
  {
    while (CompositeCoEdge* coedge_ptr = loop_ptr->first_coedge())
    {
      CompositeCurve* curve = coedge_ptr->get_curve();
      if (curve)
      {
        curve->remove(coedge_ptr);
          // Clean up point-curves with surface...
        if (curve->num_curves() == 0)
          notify_deactivated(curve);
      }
      loop_ptr->remove(coedge_ptr);
      deactivatedList.append_unique(coedge_ptr);
    }
    surface->remove(loop_ptr);
    deactivatedList.append_unique(loop_ptr);
  }
  deactivatedList.append_unique(surface);
}

void CompositeEngine::notify_deactivated (CompositeCurve* curve)
{
  while (CompositeCoEdge* coedge_ptr = curve->next_coedge(NULL))
  {
    if (coedge_ptr->get_loop())
      coedge_ptr->get_loop()->remove(coedge_ptr);
    curve->remove(coedge_ptr);
    deactivatedList.append_unique(coedge_ptr);
  }
  
  if (curve->start_point())
    curve->start_point(0);
  if (curve->end_point())
    curve->end_point(0);
    
  deactivatedList.append_unique(curve);
}

void CompositeEngine::notify_deactivated (CompositePoint* point)
{
  while (CompositeCurve* curve_ptr = point->next_curve(NULL))
  {
    if (point == curve_ptr->start_point())
      curve_ptr->start_point(0);
    if (point == curve_ptr->end_point())
      curve_ptr->end_point(0);
  }
  
  deactivatedList.append_unique(point);
}

void CompositeEngine::clean_out_deactivated_geometry()
{
  while (deactivatedList.size())
    delete deactivatedList.pop();
}

        
CubitStatus CompositeEngine::translate( CompositeBody* body, 
                                        const CubitVector& delta )
{
  assert (body && delta.length_squared());
  return CUBIT_FAILURE;
}

CubitStatus CompositeEngine::rotate( CompositeBody* body,
                                     const CubitVector& axis,
                                     double degrees )
{
  assert( body && axis.length_squared() && degrees );
  return CUBIT_FAILURE;
}

CubitStatus CompositeEngine::scale( CompositeBody* body,
                                    const CubitVector& factors )
{
  assert( body && factors.x() && factors.y() && factors.z() );
  return CUBIT_FAILURE;
}

CubitStatus CompositeEngine::reflect( CompositeBody* body,
                                      const CubitVector& axis )

{
  assert( body && axis.length() );
  return CUBIT_FAILURE;
}

CubitStatus CompositeEngine::restore_transform( CompositeBody* body )
{
  assert( 0 != body );
  return CUBIT_FAILURE;
}


CubitStatus CompositeEngine::translate( CompositeSurface* surf, 
                                        const CubitVector& delta )
{
  assert (surf && delta.length_squared());
  return CUBIT_FAILURE;
}

CubitStatus CompositeEngine::rotate( CompositeSurface* surf,
                                     const CubitVector& axis,
                                     double degrees )
{
  assert( surf && axis.length_squared() && degrees );
  return CUBIT_FAILURE;
}

CubitStatus CompositeEngine::scale( CompositeSurface* surf,
                                    const CubitVector& factors )
{
  assert( surf && factors.x() && factors.y() && factors.z() );
  return CUBIT_FAILURE;
}

CubitStatus CompositeEngine::reflect( CompositeSurface* surf,
                                      const CubitVector& axis )

{
  assert( surf && axis.length() );
  return CUBIT_FAILURE;
}


CubitStatus CompositeEngine::translate( CompositeCurve* curve, 
                                        const CubitVector& delta )
{
  assert (curve && delta.length_squared());
  return CUBIT_FAILURE;
}

CubitStatus CompositeEngine::rotate( CompositeCurve* curve,
                                     const CubitVector& axis,
                                     double degrees )
{
  assert( curve && axis.length_squared() && degrees );
  return CUBIT_FAILURE;
}

CubitStatus CompositeEngine::scale( CompositeCurve* curve,
                                    const CubitVector& factors )
{
  assert( curve && factors.x() && factors.y() && factors.z() );
  return CUBIT_FAILURE;
}

CubitStatus CompositeEngine::reflect( CompositeCurve* curve,
                                      const CubitVector& axis )

{
  assert( curve && axis.length() );
  return CUBIT_FAILURE;
}

CubitStatus CompositeEngine::notify_transform( TopologyBridge* bridge,
                                               const CubitTransformMatrix& xform )
{
  int i;
  DLIList<TopologyBridge*> lumps, shells, surfaces, tmp_list;
  
  if (BodySM* body_sm = dynamic_cast<BodySM*>(bridge))
    body_sm->get_children( lumps, true, COMPOSITE_LAYER);
  else if (Lump* lump = dynamic_cast<Lump*>(bridge))
    lumps.append(lump);
  else if (ShellSM* shell = dynamic_cast<ShellSM*>(bridge))
    shells.append(shell);
  else if (Surface* surface = dynamic_cast<Surface*>(bridge))
    surfaces.append(surface);
  else  
    return CUBIT_SUCCESS;
  
  for (i = lumps.size(); i--; )
  {
    TopologyBridge* lump = lumps.step_and_get();
    CompositeLump* comp_lump = dynamic_cast<CompositeLump*>(lump);
    if (comp_lump)
      comp_lump->update();
  }

  while (lumps.size())
  {
    TopologyBridge* lump = lumps.pop();
    tmp_list.clean_out();
    lump->get_children( tmp_list, true, COMPOSITE_LAYER );
    shells += tmp_list;
  }
  
  while (shells.size())
  {
    TopologyBridge* shell = shells.pop();
    tmp_list.clean_out();
    shell->get_children( tmp_list, true, COMPOSITE_LAYER );
    surfaces += tmp_list;
  }
  
  for (i = surfaces.size(); i--; )
  {
    TopologyBridge* surface = surfaces.step_and_get();
    CompositeSurface* composite = dynamic_cast<CompositeSurface*>(surface);
    if (composite)
    {
      composite->notify_transformed();
      composite->update();
    }
    DLIList<Curve*> curves;
    surface->curves(curves);
    int j;
    for (j = 0; j < curves.size(); j++)
    {
       CompositeCurve* comp = 
          dynamic_cast<CompositeCurve*>(curves.get_and_step());
       if (comp)
          comp->update();
    }
  }
  
  double det = xform.sub_matrix( 3, 3 ).determinant();
  bool reverse = det < 0.0;
  if (!reverse)
    return CUBIT_SUCCESS;

  // if we get here, the transform was a reflection so we need to make sure all
  // loops and coedges get reversed
  DLIList<TopologyBridge*> loops, tmp_loops, coedges;
  while (surfaces.size())
  {
    surfaces.pop()->get_children( tmp_loops, true, COMPOSITE_LAYER );
    loops += tmp_loops;
    while (tmp_loops.size())
    {
      tmp_loops.pop()->get_children( tmp_list, true, COMPOSITE_LAYER );
      coedges += tmp_list;
      tmp_list.clean_out();
    }
 }

  // reverse loops and coedges separately, since composite coedges in non-composite
  // loops also need to be reversed
  bool b_reverse_coedges = false;  
  while (loops.size())
  {
    CompositeLoop* loop = dynamic_cast<CompositeLoop*>(loops.pop());
    if (loop)
      loop->reverse(b_reverse_coedges);
  }

  // reverse the coedges
  while (coedges.size())
  {
    CompositeCoEdge* comp = dynamic_cast<CompositeCoEdge*>(coedges.pop());
    if (comp)
    {
      comp->reverse();
    }
  }

  return CUBIT_SUCCESS;
}


//-------------------------------------------------------------------------
// Purpose       : Combine bodies
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/11/04
//-------------------------------------------------------------------------
CompositeBody* CompositeEngine::combine_bodies( BodySM* body1, BodySM* body2 )
{
  CompositeBody* comp1 = dynamic_cast<CompositeBody*>(body1);
  if (!comp1)
    comp1 = dynamic_cast<CompositeBody*>(body1->owner());
  if (!comp1)
    comp1 = replace_body( body1 );
  
  CompositeBody* comp2 = dynamic_cast<CompositeBody*>(body2);
  if (!comp2)
    comp2 = dynamic_cast<CompositeBody*>(body2->owner());
  if (!comp2)
    comp2 = replace_body( body2 );
  
  assert(comp1 && comp2);
  if (comp1 == comp2)
    return comp1;
    
  comp1->combine( comp2 );
  
  while (CompositeLump* lump = comp2->next_lump(0))
  {
    comp2->remove( lump );
    comp1->add( lump );
  }
  
  delete comp2;
  return comp1;
}

void CompositeEngine::get_tbs_with_bridge_manager_as_owner( TopologyBridge *source_bridge, 
                                                            DLIList<TopologyBridge*> &tbs )
{
  TBOwner* tb_owner = source_bridge->owner();
  CompositeSurface *comp_surf = CAST_TO(tb_owner, CompositeSurface );
  if( comp_surf )  
  { 
    if( comp_surf->bridge_manager() )
    {
      tbs.append( comp_surf );
      return;
    }
  } 

  CompositeCurve *comp_curve = CAST_TO(tb_owner, CompositeCurve );
  if( comp_curve )  
  { 
    if( comp_curve->bridge_manager() )
    {
      tbs.append( comp_curve );
      return;
    }
  } 

  CompositePoint *comp_pt = CAST_TO(tb_owner, CompositePoint );
  if( comp_pt )  
  { 
    if( comp_pt->bridge_manager() )
    {
      tbs.append( comp_pt );
      return;
    }
  } 

  HiddenEntitySet *hidden_ent_set = CAST_TO(tb_owner, HiddenEntitySet );
  if( hidden_ent_set )
  {
    TBOwner *owner = hidden_ent_set->owner();
    CompositeCurve *comp_curve = CAST_TO(owner, CompositeCurve );
    if( comp_curve && comp_curve->bridge_manager() )
    {
      tbs.append( comp_curve );
      return;
    }

    CompositeSurface *comp_surf = CAST_TO(owner, CompositeSurface );
    if( comp_surf )
      return;

    assert(0);
  }
  else
  {
    CompositePoint *comp_pt = CAST_TO(tb_owner, CompositePoint);   

    if( comp_pt )
    {
      TBOwner *other_owner = comp_pt->owner();

      HiddenEntitySet *hidden_ent_set = CAST_TO(other_owner, HiddenEntitySet );
      if( hidden_ent_set )
      {
        TBOwner *owner = hidden_ent_set->owner();
        CompositeCurve *comp_curve = CAST_TO(owner, CompositeCurve );
        if( comp_curve && comp_curve->bridge_manager() )
        {
          tbs.append( comp_curve );
          return;
        }       

        CompositeSurface *comp_surf = CAST_TO(owner, CompositeSurface );
        if( comp_surf && comp_surf->bridge_manager() )
        {
          tbs.append( comp_surf );
          return;
        }       
      }
    }

    CompositeCurve *comp_curve = CAST_TO(tb_owner, CompositeCurve);   

    if( comp_curve )
    {
      TBOwner *other_owner = comp_curve->owner();

      HiddenEntitySet *hidden_ent_set = CAST_TO(other_owner, HiddenEntitySet );
      if( hidden_ent_set )
      {
        TBOwner *owner = hidden_ent_set->owner();
        CompositeSurface *comp_surf = CAST_TO(owner, CompositeSurface );
        if( comp_surf->bridge_manager() )
        {
          tbs.append( comp_surf );
          return;
        }        
      }
    }
  }

    return;
}


