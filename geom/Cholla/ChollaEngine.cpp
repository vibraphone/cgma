//- Class:       ChollaEngine
//- Description: Creates the topology for a given geometry described by facets.
//- Owner:       Steven J. Owen
//- Checked by:
//- Version:
//#include <list>
#include <set>
#include <map>
#include <vector>
#include "CubitDefines.h"
#include "ChollaEngine.hpp"
#include "DLIList.hpp"
#include "TDGeomFacet.hpp"
#include "CastTo.hpp"
#include "ChollaSkinTool.hpp"
#include "ChollaVolume.hpp"
#include "ChollaSurface.hpp"
#include "ChollaCurve.hpp"
#include "ChollaPoint.hpp"
#include "CubitFacet.hpp"
#include "CubitFacetData.hpp"
#include "CubitFacetEdge.hpp"
#include "CubitFacetEdgeData.hpp"
#include "CubitPoint.hpp"
#include "CubitPointData.hpp"
#include "FacetEntity.hpp"
#include "FacetEvalTool.hpp"
#include "FacetDataUtil.hpp"
#include "debug.hpp"
#include "TDFacetBoundaryEdge.hpp"
#include "TDFacetBoundaryPoint.hpp"
#include "CurveFacetEvalTool.hpp"
#include "Cholla.h"
#include "GfxDebug.hpp"
#include "TDFacetboolData.hpp"
#include "GMem.hpp"

//============================================================================
//Function:  ChollaEngine (PUBLIC) (constructor)
//============================================================================
ChollaEngine::ChollaEngine()
{
  hashCurveArray = NULL;
  hashCurveSize = 0;
  hashPointArray = NULL;
  hashPointSize = 0;
}

//============================================================================
//Function:  ChollaEngine (PUBLIC) (constructor)
//============================================================================
ChollaEngine::ChollaEngine(DLIList<FacetEntity*> &face_list,
                           DLIList<FacetEntity*> &edge_list,
                           DLIList<FacetEntity*> &point_list )
{
  faceList = face_list;
  edgeList = edge_list;
  pointList = point_list;
  set_up_tool_datas();
  hashCurveArray = NULL;
  hashCurveSize = 0;
  hashPointArray = NULL;
  hashPointSize = 0;
  doFlip = CUBIT_FALSE;
}

//============================================================================
//Function:  ChollaEngine (PUBLIC) (constructor)
//============================================================================
ChollaEngine::ChollaEngine(DLIList<CubitFacet*>     &facet_list,
                           DLIList<CubitFacetEdge*> &edge_list,
                           DLIList<CubitPoint*>     &point_list )
{
  CAST_LIST(facet_list, faceList, FacetEntity);
  CAST_LIST(edge_list, edgeList, FacetEntity);
  CAST_LIST(point_list, pointList, FacetEntity);
  set_up_tool_datas();
  hashCurveArray = NULL;
  hashCurveSize = 0;
  hashPointArray = NULL;
  hashPointSize = 0;
  doFlip = CUBIT_FALSE;
}

//============================================================================
//Function:  ChollaEngine (PUBLIC) (constructor)
//Notes: This case is used only when the cholla entities have been generated 
//       directly rather than using create_geometry
//       Does not use the TDGeomFacet Tooldatas
//============================================================================
ChollaEngine::ChollaEngine(DLIList<CubitFacet*>     &facet_list,
                           DLIList<CubitFacetEdge*> &edge_list,
                           DLIList<CubitPoint*>     &point_list,
                           DLIList<ChollaVolume *>  &cholla_volumes,
                           DLIList<ChollaSurface *> &cholla_surfaces,
                           DLIList<ChollaCurve *>   &cholla_curves,
                           DLIList<ChollaPoint *>   &cholla_points )
{
  CAST_LIST(facet_list, faceList, FacetEntity);
  CAST_LIST(edge_list, edgeList, FacetEntity);
  CAST_LIST(point_list, pointList, FacetEntity);
  //set_up_tool_datas();  TDGeomFacet tooldatas should have already been added

  hashCurveArray = NULL;
  hashCurveSize = 0;
  hashPointArray = NULL;
  hashPointSize = 0;
  doFlip = CUBIT_FALSE;
  
  chollaVolumeList = cholla_volumes;
  chollaSurfaceList = cholla_surfaces;
  chollaCurveList = cholla_curves;
  chollaPointList = cholla_points;
}
                          

//============================================================================
//Function:  set_up_tool_datas
//============================================================================
void ChollaEngine::set_up_tool_datas( )
{
  int ii;
  FacetEntity *fe_ptr;
  for (ii=0; ii<faceList.size(); ii++)
  {
    fe_ptr = faceList.get_and_step();
    TDGeomFacet::add_geom_facet( fe_ptr, -1 );
  }
  for (ii=0; ii<pointList.size(); ii++)
  {
    fe_ptr = pointList.get_and_step();
    TDGeomFacet::add_geom_facet( fe_ptr, -1 );
  }
}

CubitStatus ChollaEngine::build_eval_tools()
{
  CubitBoolean use_feature_angle = CUBIT_FALSE; 
  double min_dot = 0.0;                    
  int interp_order = 0;                 
  CubitBoolean smooth_non_manifold = CUBIT_FALSE;
  CubitBoolean split_surfaces = CUBIT_FALSE;

  return build_eval_tools( chollaSurfaceList,
                           chollaCurveList,
                           interp_order, use_feature_angle, 
                           min_dot, smooth_non_manifold,
                           split_surfaces );
}

//============================================================================
//Function:  delete_tool_datas
//============================================================================
void ChollaEngine::delete_tool_datas( )
{
  int ii;
  FacetEntity *fe_ptr;
  for (ii=0; ii<faceList.size(); ii++)
  {
    fe_ptr = faceList.get_and_step();
    fe_ptr->delete_TD( &TDGeomFacet::is_geom_facet );
  }
  for (ii=0; ii<edgeList.size(); ii++)
  {
    fe_ptr = edgeList.get_and_step();
    fe_ptr->delete_TD( &TDGeomFacet::is_geom_facet );
  }
  for (ii=0; ii<pointList.size(); ii++)
  {
    fe_ptr = pointList.get_and_step();
    fe_ptr->delete_TD( &TDGeomFacet::is_geom_facet );
  }
}

//============================================================================
//Function:  delete_eval_tools
//============================================================================
void ChollaEngine::delete_eval_tools()
{

  delete_tool_datas();
  int ii;
  ChollaSurface *cs_ptr = NULL;
  FacetEvalTool *fe_tool_ptr;
  for (ii = chollaSurfaceList.size(); ii > 0; ii-- )
  {
    cs_ptr = chollaSurfaceList.get_and_step();
    fe_tool_ptr = cs_ptr->get_eval_tool();
    if (fe_tool_ptr)
      delete fe_tool_ptr;
  }
  ChollaCurve *cc_ptr = NULL;
  CurveFacetEvalTool *ce_tool_ptr;
  for (ii = chollaCurveList.size(); ii > 0; ii-- )
  {
    cc_ptr = chollaCurveList.get_and_step();
    ce_tool_ptr = cc_ptr->get_eval_tool();
    if (ce_tool_ptr)
      delete ce_tool_ptr;
  }
  faceList.clean_out();
  edgeList.clean_out();
  pointList.clean_out();
}

//============================================================================
//Function:  delete_eval_tools_but_not_facets
//============================================================================
void ChollaEngine::delete_eval_tools_but_not_facets()
{

  delete_tool_datas();
  int ii;
  ChollaSurface *cs_ptr = NULL;
  FacetEvalTool *fe_tool_ptr;
  for (ii = chollaSurfaceList.size(); ii > 0; ii-- )
  {
    cs_ptr = chollaSurfaceList.get_and_step();
    fe_tool_ptr = cs_ptr->get_eval_tool();
    if (fe_tool_ptr)
    {
      DLIList<CubitFacet *>facets;
      fe_tool_ptr->remove_facets(facets);
      delete fe_tool_ptr;
    }
  }
  ChollaCurve *cc_ptr = NULL;
  CurveFacetEvalTool *ce_tool_ptr;
  for (ii = chollaCurveList.size(); ii > 0; ii-- )
  {
    cc_ptr = chollaCurveList.get_and_step();
    ce_tool_ptr = cc_ptr->get_eval_tool();
    if (ce_tool_ptr)
    {
      delete ce_tool_ptr;
    }
  }
  faceList.clean_out();
  edgeList.clean_out();
  pointList.clean_out();
}

//============================================================================
//Function:  ~ChollaEngine (PUBLIC) (destructor)
//============================================================================
ChollaEngine::~ChollaEngine()
{
}

//============================================================================
//Function:  delete (PUBLIC)
//Description: removes all cholla entities stored with the cholla engine
//             usually called before destructor
//============================================================================
void ChollaEngine::delete_me()
{
  delete_tool_datas();
  int ii;
    //clean up any data remaining.
  for (ii = chollaPointList.size(); ii > 0; ii-- )
    delete chollaPointList.remove();
  for (ii = chollaCurveList.size(); ii > 0; ii-- )
    delete chollaCurveList.remove();
  for (ii = chollaSurfaceList.size(); ii > 0; ii-- )
    delete chollaSurfaceList.remove();
  
}

//==================================================================================
//Function:  create_geometry (PUBLIC)
//Description:  Interface function to acutally create the topology
//              for the given mesh.
//==================================================================================
CubitStatus ChollaEngine::create_geometry(
  CubitBoolean use_feature_angle,   // use an angle to define where to break surfaces
  double angle,                     // the feature angle
  int interp_order,                 // 0=linear, 4=b-spline patches
  CubitBoolean smooth_non_manifold, // check for continuity accross >2 valence edges
  CubitBoolean split_surfaces)      // create new FacetEntities to split surfaces at
                                    // at features.  Otherwise - don't mess with the
                                    // FacetEntities. (creates TDs instead)

{
  if (0)
  {
    dump("cyl.cholla", angle);
  }

  //- convert feature angle to a dot product

  double min_dot = 0.0;
  if (use_feature_angle)
  {
    if (angle > 180.0) angle = 180.0;
    if (angle < 0.0) angle = 0.0;
    double rad_angle = (180.0 - angle) * CUBIT_PI / 180.0;
    min_dot = cos( rad_angle );
  }

  //- create one facet surface to start with

  ChollaSkinTool c_skin_tool;
  ChollaSurface *cholla_surface_ptr = NULL;
  CubitStatus stat = c_skin_tool.skin_2d(faceList, cholla_surface_ptr);
  if ( stat != CUBIT_SUCCESS )
    return stat;
  if ( cholla_surface_ptr )
    chollaSurfaceList.append(cholla_surface_ptr);

  // before building the surfaces, orient the facets so they are consistent
  
  DLIList<CubitFacet *>total_facet_list;
  CAST_LIST( faceList, total_facet_list, CubitFacet );
  //stat = check_all_facet_orientations( total_facet_list, doFlip );
  if (stat!=CUBIT_SUCCESS)
    return stat;
  
  // generate the topology from the facets

  stat = create_volume_boundaries( chollaSurfaceList,
                                   use_feature_angle, min_dot, 
                                   split_surfaces );
  if ( stat == CUBIT_FAILURE )
    return stat;

  stat = create_surface_boundaries( chollaSurfaceList, chollaCurveList,
                                    use_feature_angle, min_dot );
  if ( stat == CUBIT_FAILURE )
    return stat;

  stat = create_curve_boundaries( chollaCurveList, chollaPointList );
  if ( stat == CUBIT_FAILURE )
    return stat;

    // Okay.  Now we are ready to actually build the geometry.
  int mydebug = 0;
  if (mydebug)
    print_me();

  stat = build_eval_tools( chollaSurfaceList,
                           chollaCurveList,
                           interp_order, use_feature_angle, 
                           min_dot, smooth_non_manifold,
                           split_surfaces );

  return stat;
}

//=============================================================================
//Function:  create_volume_boundaries (PRIVATE)
//Description: creates the surfaces based on the sideset and element block
//             information
//Author: sjowen
//Date: 10/17/00
//=============================================================================
CubitStatus ChollaEngine::create_volume_boundaries( 
  DLIList<ChollaSurface*> &cholla_surface_sheets,  // output global list of surfaces
  CubitBoolean use_feature_angle,   // define surfaces based on feature angle
  double min_dot,      // minimum dot product between face normals
  CubitBoolean split_surfaces )  // create new FacetEntities to split surfaces at
                                 // at features.  Otherwise - don't mess with the
                                 // FacetEntities. (creates TDs instead)
{
  CubitStatus rv = CUBIT_SUCCESS;

  // Split these surfaces so that there is only one set of continuous
  // faces per surface 

  int ii;
  int num_surfaces = cholla_surface_sheets.size();
  cholla_surface_sheets.reset();
  for ( ii = num_surfaces; ii > 0; ii-- )
  {

    ChollaSurface *chsurf_ptr = cholla_surface_sheets.get_and_step();

    DLIList<FacetEntity*> surf_facet_list;
    FacetEntity *sfacet;    
    chsurf_ptr->get_facets(surf_facet_list);
    DLIList<CubitFacet*> facet_list;
    CubitFacet *facet;
    CubitFacetEdge *fedge;
    int kk, mm;
    int mydebug = 0;
    if(mydebug){
      GfxDebug::clear();
    }
    DLIList<CubitFacetEdge *> feature_edge_list;
    for ( kk = surf_facet_list.size(); kk > 0; kk-- ) {
      sfacet = surf_facet_list.get_and_step();
      facet_list.clean_out();
      sfacet->facets(facet_list);
      facet = facet_list.get();
      TDFacetboolData* tdf = TDFacetboolData::get(facet);
      if ( tdf ) {
        int *cptr;
        cptr = tdf->get_edge_indices( (bool) facet->is_backwards());
        for ( mm = 0; mm < 3; mm++ ) {
          if ( cptr[mm] != 0 ) {
            fedge = facet->edge((2+mm)%3);
            if(mydebug){
              if(facet->is_backwards()){
                facet->edge((1+mm)%3)->debug_draw(CUBIT_GREEN);
                facet->edge((mm)%3)->debug_draw(CUBIT_RED);
              
                fedge->debug_draw(CUBIT_BLUE);
              }
              else
                fedge->debug_draw(CUBIT_WHITE);
            }
            
            TDGeomFacet::add_geom_facet(fedge, -1);
            fedge->set_as_feature();
            feature_edge_list.append( fedge );
          }            
        }
      }
    }
    if(mydebug){
      GfxDebug::mouse_xforms();
    }
    // make a list of feature edges

    rv = chsurf_ptr->add_preexisting_feature_edges( feature_edge_list );
		if (rv != CUBIT_SUCCESS)
			return rv;
		
    if (use_feature_angle)
    {
      rv = chsurf_ptr->feature_angle( min_dot, feature_edge_list );
    }
    else
    {
      rv = chsurf_ptr->non_manifold_edges( feature_edge_list );
    }
    if (rv != CUBIT_SUCCESS)
      return rv;

     // crack the surface at the feature edges.  create new edges and
     // points so the facet representation is discontinuous.

    rv = make_features( feature_edge_list, split_surfaces );
    if (rv != CUBIT_SUCCESS)
       return rv;

    // split up the surface

    rv = chsurf_ptr->split_surface( cholla_surface_sheets );
    if (rv != CUBIT_SUCCESS)
       return rv;
  }

  // Clean up any edges that do not form complete loops as a result of 
  // feature angle.  For this implementation we will allow features to exist
  // within the surface without defining a complete loop. -- so the next
  // piece of code is never executed.

  CubitBoolean use_complete_loops_only = CUBIT_FALSE;
  if (use_feature_angle && use_complete_loops_only)
  {
    for ( ii = cholla_surface_sheets.size(); ii > 0; ii-- )
    {
      ChollaSurface *chsurf_ptr = cholla_surface_sheets.get_and_step();
      chsurf_ptr->clean_features();
    }
  }

  // Now that we've broken everything into surfaces, update the surface IDs
  // on the boundary facet tool data

  if (!split_surfaces)
  {
    for ( ii = cholla_surface_sheets.size(); ii > 0; ii-- )
    {
      ChollaSurface *chsurf_ptr = cholla_surface_sheets.get_and_step();
      chsurf_ptr->update_boundary_tool_data();
    }
  }
  return rv;
}


//=============================================================================
//Function:  create_surface_boundaries (PRIVATE)
//Description: creates the curves based on the facet surface information
//Author: sjowen
//Date: 12/3/00
//=============================================================================
CubitStatus ChollaEngine::create_surface_boundaries(
  DLIList<ChollaSurface*> &cholla_surface_list, // global list of surfaces
  DLIList<ChollaCurve*> &cholla_curve_list,     // output global list of curves
  CubitBoolean use_feature_angle,
  double min_dot )
{
  CubitStatus stat = CUBIT_SUCCESS;

    // determine the boundaries for each surface.  One curve per surface

  int ii;
  for ( ii = cholla_surface_list.size(); ii > 0; ii-- )
  {
    ChollaSurface *chsurf_ptr = cholla_surface_list.get_and_step();
    DLIList<ChollaCurve*> chcurve_list;
    chsurf_ptr->get_curves( chcurve_list );
    if (chcurve_list.size() == 0)
    {
      DLIList<FacetEntity*> facet_list;
      chsurf_ptr->get_facets(facet_list);
      ChollaSkinTool c_skin_tool;
      stat = c_skin_tool.skin_2d(facet_list, chsurf_ptr);
      if ( stat != CUBIT_SUCCESS )
        return stat;
    }
  }

  // create a hash list of curves - to speed up classification

  stat = init_hash_curves();
  if (stat != CUBIT_SUCCESS)
  {
    delete_hash_curves();
    return stat;
  }

  // loop through each of the edges on the surfaces
  // Determine which curve it is a part of.
  // Create a new ChollaCurve for each curve
  // Curves are created wherever there is a unique set of associated
  // surfaces 

  int jj, kk;
  for ( ii = cholla_surface_list.size(); ii > 0; ii-- )
  {

    ChollaSurface *chsurf_ptr = cholla_surface_list.get_and_step();
    DLIList<ChollaCurve*> chcurve_list;
    chsurf_ptr->get_curves( chcurve_list );

    // curently there should only be one curve list per surface

    for (jj=chcurve_list.size(); jj>0; jj--)
    {
      ChollaCurve *chcurv_ptr = chcurve_list.get_and_step();
      DLIList<FacetEntity*> facet_list =  chcurv_ptr->get_facet_list();
      FacetEntity *edge_ptr;
      for ( kk = 0; kk < facet_list.size(); kk++)
      {
        edge_ptr = facet_list.get_and_step();
        stat = classify_edge( edge_ptr, cholla_curve_list, chsurf_ptr );
        if (stat != CUBIT_SUCCESS) 
          return stat;
      }

      // delete this ChollaCurve - it should have been replaced by one
      // or more curves bounding this surface

      chsurf_ptr->remove_curve( chcurv_ptr );
      delete chcurv_ptr;
    }
  }
  delete_hash_curves();

  // Split these curves so that there is only one string of continuous
  // edges per curve (it will also order the edges and set the start
  // and end nodes for each curve) 

  int num_curves = cholla_curve_list.size();

  cholla_curve_list.reset();
  for ( ii = num_curves; ii > 0 && stat == CUBIT_SUCCESS; ii-- )
  {
    ChollaCurve *chcurv_ptr = cholla_curve_list.get();

    // if necessary mark nodes that will serve as feature breaks (vertices 
    // will be generatedat them)

    if (use_feature_angle)
    {
      stat = chcurv_ptr->feature_angle( min_dot );
      if (stat != CUBIT_SUCCESS)
         return stat;
    }

    // split the curve based on various criteria

    stat = chcurv_ptr->split_curve( cholla_curve_list );

    // delete this curve (new curves were created in split_curve)

    delete chcurv_ptr;
    cholla_curve_list.change_to(NULL);
    cholla_curve_list.step();
  }
  cholla_curve_list.remove_all_with_value(NULL);

  // update the point->curve associativity

  for (ii=0; ii<cholla_curve_list.size(); ii++)
  {
    ChollaCurve *chcurv_ptr = cholla_curve_list.get_and_step();
    CubitPoint *start_ptr, *end_ptr;
    chcurv_ptr->get_ends( start_ptr, end_ptr );
    TDGeomFacet *td = TDGeomFacet::get_geom_facet( start_ptr );
    td->add_cholla_curve( chcurv_ptr );
    td = TDGeomFacet::get_geom_facet( end_ptr );
    td->add_cholla_curve( chcurv_ptr );
  }

  return stat;
}

//=============================================================================
//Function:  init_hash_curves (PRIVATE)
//Description: create a hash array of all curves.  They are hashed based on the
//             smallest id of any surface attached to the curve 
//Author: sjowen
//Date: 3/7/01
//=============================================================================
CubitStatus ChollaEngine::init_hash_curves( )
{
  /* === find the next highest prime number */

  int num_surfs = chollaSurfaceList.size();
  int i;
  hashCurveSize = num_surfs;
  if (num_surfs < 30) hashCurveSize = 31;
  else 
  {
    i=2;
    while (i<hashCurveSize*0.5 + 1) {
      if (hashCurveSize % i == 0) {
        i=2;
        hashCurveSize++;
      }
      else {
        i++;
      }
    }
  }
  hashCurveArray = new DLIList<ChollaCurve*>[hashCurveSize];

  int key = 0;
  ChollaCurve *chcurv_ptr;
  DLIList<ChollaSurface*> *fsm_list_ptr; 
  for (i=0; i<chollaCurveList.size(); i++)
  {
    chcurv_ptr = chollaCurveList.get_and_step();
    fsm_list_ptr = chcurv_ptr->get_surface_list_ptr();
    key = get_curve_hash_key( fsm_list_ptr );
    hashCurveArray[key].append( chcurv_ptr );   
  }
  return CUBIT_SUCCESS;
}

//=============================================================================
//Function:  delete_hash_curves (PRIVATE)
//Description: delete the hash curve stuff
//Author: sjowen
//Date: 3/7/01
//=============================================================================
void ChollaEngine::delete_hash_curves( )
{
  if (hashCurveArray)
   delete [] hashCurveArray;
  hashCurveArray = NULL;
  hashCurveSize = 0;
}

//=============================================================================
//Function:  get_curve_hash_key (PRIVATE)
//Description: 
//Author: sjowen
//Date: 3/7/01
//=============================================================================
int ChollaEngine::get_curve_hash_key(
  DLIList<ChollaSurface*> *fsm_list_ptr )
{
  int key, j;
  ChollaSurface *chsurf_ptr;
  if (fsm_list_ptr->size() == 0)
  {
    key = 0;
  }
  else
  {
    key = INT_MAX;
    for (j=0; j<fsm_list_ptr->size(); j++)
    {
      chsurf_ptr = fsm_list_ptr->get_and_step();
      if (chsurf_ptr->get_id() < key)
        key = chsurf_ptr->get_id();
    }
  }
  key = key % hashCurveSize;
  return key;
}


//=============================================================================
//Function:  classify_edge (PRIVATE)
//Description: sorts a edge into its correct curve based on its associated
//             surfaces and sidesets/nodesets.  Creates a new block curve if 
//             necessary.
//Author: sjowen
//Date: 12/3/00
//=============================================================================
CubitStatus ChollaEngine::classify_edge(
  FacetEntity *edge_ptr,       // the edge we are classifying
  DLIList<ChollaCurve*> &cholla_curve_list,  // add to one of these
  ChollaSurface *cholla_surf_mesh_ptr )   // the current surface
{
  CubitStatus rv = CUBIT_SUCCESS;

  // see if we have already classified this edge (from another surface)

  TDGeomFacet *td_gm_edge = TDGeomFacet::get_geom_facet(edge_ptr);
  if (td_gm_edge->get_hit_flag() != 0)
    return rv;
  td_gm_edge->set_hit_flag(1);

  // get the surfaces adjacent to this edge

  DLIList<ChollaSurface*> this_chsurf_list;
  td_gm_edge->get_cholla_surfs( this_chsurf_list );
  int this_num_adj = this_chsurf_list.size();
   
  // see if the surfaces defined on this face match any 
  // of the existing block curves

  DLIList<ChollaSurface*> *chsurf_list_ptr;
  ChollaCurve *chcurv_ptr = NULL;
  int found = 0;
  int key = get_curve_hash_key( &this_chsurf_list );
  int ii;
  for (ii=0; ii<hashCurveArray[key].size() && !found; ii++)
  {
    chcurv_ptr = hashCurveArray[key].get();
    // the first one checked should be the same as the last one checked (don't
    // use get_and_step here)  This should speed things up

    // check if surfaces are the same

    chsurf_list_ptr = chcurv_ptr->get_surface_list_ptr( );
    int num_adj = chsurf_list_ptr->size();
    if (num_adj == this_num_adj)
    {
      found = 1;
      int jj, kk;
      for (jj=chsurf_list_ptr->size(); jj>0 && found; jj--)
      {
        int same_surf = 0;
        ChollaSurface *chsurf_ptr = chsurf_list_ptr->get_and_step();
        for(kk=this_chsurf_list.size(); kk>0 && !same_surf; kk--)
        {
          ChollaSurface *this_chsurf_ptr = this_chsurf_list.get_and_step();
          if (this_chsurf_ptr == chsurf_ptr)
          {
            same_surf = 1;
          }
        }
        if (!same_surf)
          found = 0;
      }
    }
    if (!found)
      hashCurveArray[key].step();
  }

  // if the unique set of surfaces that this edge is associated 
  // with was found to already exist for a facet curve -- add the
  // edge to the block curve

  if (found)
  {

    // add the edge to the block curve mesh (make sure it is only added once)

    //int was_added = 
    chcurv_ptr->add_facet_unique( edge_ptr );

    // add the curve to the surface (if needed)

    cholla_surf_mesh_ptr->add_curve_unique( chcurv_ptr );

    // add the curve to the edge
  
    td_gm_edge->add_cholla_curve( chcurv_ptr );
  }

  // if the unique set of surfaces that this edge is associated
  // with is not found, then create a new facet curve and add the edge to it

  else
  {

    // create it and update surface and nodeset info

    int block_id = td_gm_edge->get_block_id();
    ChollaCurve *new_chcurv_ptr = new ChollaCurve( block_id ); 
    for (int mm=0; mm<this_num_adj; mm++)
    {
      new_chcurv_ptr->add_surface( this_chsurf_list.get_and_step() );
    }

    // add the edge

    new_chcurv_ptr->add_facet( edge_ptr );

    // update the surface with this new curve

    cholla_surf_mesh_ptr->add_curve( new_chcurv_ptr );

    // add the new curve to the global list

    cholla_curve_list.append( new_chcurv_ptr ); 

    // add the curve to the edge
  
    td_gm_edge->add_cholla_curve( new_chcurv_ptr );

    // add the new curve to the hash table

    hashCurveArray[key].append( new_chcurv_ptr );
  }
  return rv;
}

//=============================================================================
//Function:  create_curve_boundaries (PRIVATE)
//Description: creates the points based on the nodeset and curve information
//Author: sjowen
//Date: 12/6/00
//=============================================================================
CubitStatus ChollaEngine::create_curve_boundaries(
  DLIList<ChollaCurve*> &cholla_curve_list,   // global list of curves
  DLIList<ChollaPoint*> &cholla_point_list )  // output global list of points
{

  CubitStatus stat = CUBIT_SUCCESS;

  // hash the points for speed

  stat = init_hash_points();
  if (stat != CUBIT_SUCCESS)
  {
    delete_hash_points();
    return stat;
  }

  // loop through each of the end nodes on the curves
  // Determine which point it is a part of.
  // Create a new ChollaPoint for each point

  //int mydebug = 0;
  int ii, kk;
  for ( ii = cholla_curve_list.size(); ii > 0; ii-- )
  {
    ChollaCurve *chcurv_ptr = cholla_curve_list.get_and_step();
    CubitPoint *point_ptr[2];
    chcurv_ptr->get_ends( point_ptr[0], point_ptr[1] );
    for ( kk = 0; kk < 2; kk++)
    {
      stat = classify_point( point_ptr[kk], cholla_point_list, chcurv_ptr );
      if (stat != CUBIT_SUCCESS)
      {
        delete_hash_points();
        return stat;
      }
    }
  }
  delete_hash_points();
  return stat;
}

//=============================================================================
//Function:  classify_point (PRIVATE)
//Description: sorts a point into its correct point based on its associated
//             curve.  Creates a new block point if necessary
//Author: sjowen
//Date: 12/6/00
//=============================================================================
CubitStatus ChollaEngine::classify_point(
  CubitPoint *point_ptr,                     // the node to classify
  DLIList<ChollaPoint*> &cholla_point_list,   // global list of points
  ChollaCurve *chcurv_ptr )                 // curve that the end point is on
{
  int ii;
  int found = 0;
  TDGeomFacet *td_node = TDGeomFacet::get_geom_facet( point_ptr );
  ChollaPoint *chpnt_ptr = NULL;
  DLIList<ChollaCurve*> fcm_list;
  td_node->get_cholla_curves(fcm_list);
  int key = get_point_hash_key( &fcm_list );
  for (ii = 0; ii < hashPointArray[key].size() && !found; ii++)
  {
    chpnt_ptr = hashPointArray[key].get();
    FacetEntity *this_point_ptr = chpnt_ptr->get_facets();
    if (this_point_ptr == point_ptr)
    {
      found = 1;
    }
    else
    {
      hashPointArray[key].step();
    }
  } 

  if (found)
  {
    chpnt_ptr->add_curve( chcurv_ptr );
    chcurv_ptr->add_point( chpnt_ptr );
  }
  else
  {
    ChollaPoint *new_chpnt_ptr = new ChollaPoint(); 
    new_chpnt_ptr->add_facet( point_ptr );
    new_chpnt_ptr->add_curve( chcurv_ptr );
    chcurv_ptr->add_point( new_chpnt_ptr );
    cholla_point_list.append( new_chpnt_ptr );
    hashPointArray[key].append( new_chpnt_ptr );
  }
  return CUBIT_SUCCESS;
}

//=============================================================================
//Function:  init_hash_points (PRIVATE)
//Description: create a hash array of all points.  They are hashed based on the
//             smallest id of any curve attached to the curve 
//Author: sjowen
//Date: 3/7/01
//=============================================================================
CubitStatus ChollaEngine::init_hash_points( )
{
  /* === find the next highest prime number */

  int num_curves = chollaCurveList.size();
  int i;
  hashPointSize = num_curves;
  if (num_curves < 30) hashPointSize = 31;
  else 
  {
    i=2;
    while (i<hashPointSize*0.5 + 1) {
      if (hashPointSize % i == 0) {
        i=2;
        hashPointSize++;
      }
      else {
        i++;
      }
    }
  }
  hashPointArray = new DLIList<ChollaPoint*>[hashPointSize];

  int key = 0;
  ChollaPoint *chpnt_ptr;
  DLIList<ChollaCurve*> *fcm_list_ptr; 
  for (i=0; i<chollaPointList.size(); i++)
  {
    chpnt_ptr = chollaPointList.get_and_step();
    fcm_list_ptr = chpnt_ptr->get_curve_list_ptr();
    key = get_point_hash_key( fcm_list_ptr );
    hashPointArray[key].append( chpnt_ptr );   
  }
  return CUBIT_SUCCESS;
}

//=============================================================================
//Function:  delete_hash_points (PRIVATE)
//Description: delete the hash point stuff
//Author: sjowen
//Date: 3/7/01
//=============================================================================
void ChollaEngine::delete_hash_points( )
{
  if (hashPointArray)
   delete [] hashPointArray;
  hashPointArray = NULL;
  hashPointSize = 0;
}

//=============================================================================
//Function:  get_point_hash_key (PRIVATE)
//Description: 
//Author: sjowen
//Date: 3/7/01
//=============================================================================
int ChollaEngine::get_point_hash_key(
  DLIList<ChollaCurve*> *fcm_list_ptr )
{
  int key, j;
  ChollaCurve *chcurv_ptr;
  if (fcm_list_ptr->size() == 0)
  {
    key = 0;
  }
  else
  {
    key = INT_MAX;
    for (j=0; j<fcm_list_ptr->size(); j++)
    {
      chcurv_ptr = fcm_list_ptr->get_and_step();
      if (chcurv_ptr->get_id() < key)
        key = chcurv_ptr->get_id();
    }
  }
  key = key % hashPointSize;
  return key;
}

//=============================================================================
//Function:  facet_dimension (PRIVATE)
//Description: returns the dimension of a mesh entity.
//=============================================================================
int ChollaEngine::facet_dimension(FacetEntity *facet_ptr)
{
  if (CAST_TO(facet_ptr, CubitFacet) != NULL)
    return 2;
  if (CAST_TO(facet_ptr, CubitFacetEdge) != NULL)
    return 1;
  if (CAST_TO(facet_ptr, CubitPoint) != NULL)
    return 0;
  return -1;
}


//===============================================================================
//Function:  build_eval_tools (PRIVATE)
//Description:  build the CUBIT geometry based on the Facet entity class lists
//===============================================================================
CubitStatus ChollaEngine::build_eval_tools(
  DLIList<ChollaSurface*> &cholla_surface_list,
  DLIList<ChollaCurve*> &cholla_curve_list,
  int interp_order,
  CubitBoolean use_feature_angle,
  double min_dot,
  CubitBoolean smooth_non_manifold,
  CubitBoolean split_surfaces )
{
  CubitStatus stat = CUBIT_SUCCESS;  

  if (stat == CUBIT_SUCCESS)
    stat = build_surface_and_curve_eval_tools( cholla_surface_list,
                                   interp_order, min_dot );  
  if (stat == CUBIT_SUCCESS)
    stat = build_curve_eval_tools( cholla_curve_list, interp_order );

  if (stat == CUBIT_SUCCESS)
    if (interp_order == 4)
      stat = clean_geometry( smooth_non_manifold, split_surfaces,
                           use_feature_angle, min_dot, cholla_curve_list );
  return stat;
}


//===============================================================================
//Function:  build_curve_eval_tools (PRIVATE)
//Description:  From the cholla curve list, create CurveFacetEvalTools
//===============================================================================
CubitStatus ChollaEngine::build_curve_eval_tools(
  DLIList<ChollaCurve*> &cholla_curve_list,
  int interp_order )
{
  CubitStatus stat = CUBIT_SUCCESS;
  int kk;
  interp_order = interp_order;
  for ( kk = cholla_curve_list.size(); kk > 0; kk-- )
  {
    ChollaCurve *chcurv_ptr = cholla_curve_list.get_and_step();
    CurveFacetEvalTool *curv_eval_tool_ptr = chcurv_ptr->get_eval_tool();
    if (curv_eval_tool_ptr == NULL)
    {
      CubitPoint *start_point, *end_point;
      chcurv_ptr->get_ends( start_point, end_point );
      
      // if this is a curve without a parent surface then handle it 
      // differently.  (Curves with parents use the surface to evaluate to
      // With only a curve, it must evaluate to the curve)
      // Note the CurveFacetEvalTool for curves that have parent surfaces
      // are computed in build_surface_eval_tools
      
      DLIList<ChollaSurface*> chsurf_list = chcurv_ptr->get_surfaces();
      if (chsurf_list.size() == 0)
      {
        DLIList<FacetEntity*> facet_list;
        DLIList<CubitPoint*> point_list;  // needs to be filled in
        facet_list = chcurv_ptr->get_facet_list();
        DLIList<CubitFacetEdge*> edge_list;
        CAST_LIST( facet_list, edge_list, CubitFacetEdge );
        if (stat != CUBIT_SUCCESS)
          return stat;     

        curv_eval_tool_ptr = new CurveFacetEvalTool;
        stat = curv_eval_tool_ptr->initialize( edge_list, point_list );

        if ( stat != CUBIT_SUCCESS )
        {
            return stat;
        }
        chcurv_ptr->assign_eval_tool(curv_eval_tool_ptr);
      }
      else
      {       
        //fix up the orientation of the Cholla curve
        stat = chcurv_ptr->order_edges();
        if ( stat != CUBIT_SUCCESS )
        {         
          return stat;
        }
        
        DLIList<FacetEntity*> facet_ents;
        facet_ents = chcurv_ptr->get_facet_list();
        DLIList<CubitFacetEdge*> edge_facets;
        CAST_LIST(facet_ents, edge_facets, CubitFacetEdge );
        DLIList<CubitPoint*> ordered_points;
        
        for( int k=0; k<edge_facets.size(); k++ )
        {
          CubitFacetEdge *tmp_edge = edge_facets.get_and_step();          
          if( k==0 )
          {
            ordered_points.append( tmp_edge->point(0) );            
            ordered_points.append( tmp_edge->point(1) );            
          }
          else        
            ordered_points.append( tmp_edge->point(1) );                      
        }     

        FacetEvalTool *surf_eval_tool = chsurf_list.get()->get_eval_tool();
        CurveFacetEvalTool *curv_eval_tool_ptr = new CurveFacetEvalTool;       
        edge_facets.reset();        
        stat = curv_eval_tool_ptr->initialize( edge_facets, ordered_points, surf_eval_tool );          

        if ( stat != CUBIT_SUCCESS )
        {
          return stat;
        }
        chcurv_ptr->assign_eval_tool( curv_eval_tool_ptr );           
      }
    }   
  }
  return stat;
}

//===============================================================================
//Function:  build_surface_eval_tools (PRIVATE)
//Description:  From the facet surface list, create the FacetEvalTools
//===============================================================================
CubitStatus ChollaEngine::build_surface_and_curve_eval_tools(
  DLIList<ChollaSurface*> &cholla_surface_list,
  int interp_order,
  double min_dot)
{
  CubitStatus stat = CUBIT_SUCCESS;
  int ii, kk;

  // make sure the facet flags have been reset

  for ( kk = cholla_surface_list.size(); kk > 0; kk-- )
  {
    ChollaSurface *chsurf_ptr = cholla_surface_list.get_and_step();
    chsurf_ptr->reset_facet_flags();
  }

  // now loop through surfaces and create them

  int mydebug = 0;
  for ( kk = cholla_surface_list.size(); kk > 0; kk-- )
  {
    ChollaSurface *chsurf_ptr = cholla_surface_list.get_and_step();
    DLIList<FacetEntity*> facet_entity_list;    
    chsurf_ptr->get_facets(facet_entity_list);
    DLIList<CubitFacet*> facet_list;
    CAST_LIST( facet_entity_list, facet_list, CubitFacet );

    FacetEvalTool *eval_tool_ptr = new FacetEvalTool();
    eval_tool_ptr->replace_facets(facet_list);

    chsurf_ptr->assign_eval_tool(eval_tool_ptr);

    // go through each of this surface's curves and create CurveFacetEvalTools

    DLIList<ChollaCurve *> chcurv_list;
    chsurf_ptr->get_curves( chcurv_list );
    for (ii=0; ii<chcurv_list.size(); ii++)
    {
      ChollaCurve *chcurv_ptr = chcurv_list.get_and_step();   

      if (chcurv_ptr->get_eval_tool() == NULL)
      {
        //fix up the orientation of the Cholla curve
        stat = chcurv_ptr->order_edges();
        if ( stat != CUBIT_SUCCESS )
        {         
          PRINT_ERROR("Problems ordering edges!!!!!\n");          
          return stat;
        }        
        
        DLIList<FacetEntity*> facet_ents;
        facet_ents = chcurv_ptr->get_facet_list();
        DLIList<CubitFacetEdge*> edge_facets;
        CAST_LIST(facet_ents, edge_facets, CubitFacetEdge );
        DLIList<CubitPoint*> ordered_points;
        
        for( int k=0; k<edge_facets.size(); k++ )
        {
          CubitFacetEdge *tmp_edge = edge_facets.get_and_step();          
          if( k==0 )
          {
            ordered_points.append( tmp_edge->point(0) );            
            ordered_points.append( tmp_edge->point(1) );            
          }
          else        
            ordered_points.append( tmp_edge->point(1) );                      
        }     

        CurveFacetEvalTool *curv_eval_tool_ptr = new CurveFacetEvalTool;       
        edge_facets.reset();        
        stat = curv_eval_tool_ptr->initialize( edge_facets, ordered_points, eval_tool_ptr );          

        if ( stat != CUBIT_SUCCESS )
        {
          return stat;
        }
        chcurv_ptr->assign_eval_tool( curv_eval_tool_ptr );     
      }
    }
  }
  return stat;
}

CubitStatus ChollaEngine::rebuild_surface_and_curve_eval_tools(
  DLIList<ChollaSurface*> &cholla_surface_list,
  int interp_order,
  double min_dot )
{
  CubitStatus stat = CUBIT_SUCCESS;
  int ii, kk;

  // make sure the facet flags have been reset
  for ( kk = cholla_surface_list.size(); kk > 0; kk-- )
  {
    ChollaSurface *chsurf_ptr = cholla_surface_list.get_and_step();
    chsurf_ptr->reset_facet_flags();    
  }

  // get unique list of curves
  DLIList<ChollaCurve *> all_chcurves;
  for ( kk = cholla_surface_list.size(); kk > 0; kk-- )
  {
    ChollaSurface *chsurf_ptr = cholla_surface_list.get_and_step();
    DLIList<ChollaCurve *> chcurv_list;
    chsurf_ptr->get_curves( chcurv_list );
    all_chcurves += chcurv_list;    
  }
  all_chcurves.uniquify_ordered();

  // now loop through surfaces and create them

  int mydebug = 0;
  for ( kk = cholla_surface_list.size(); kk > 0; kk-- )
  {
    ChollaSurface *chsurf_ptr = cholla_surface_list.get_and_step();
    DLIList<FacetEntity*> facet_entity_list;
    DLIList<CubitPoint*> point_list;
    chsurf_ptr->get_points(point_list);
    chsurf_ptr->get_facets(facet_entity_list);
    DLIList<CubitFacet*> facet_list;
    CAST_LIST( facet_entity_list, facet_list, CubitFacet );

    FacetEvalTool *eval_tool_ptr = chsurf_ptr->get_eval_tool();
    if (NULL != eval_tool_ptr)
    {
      DLIList<CubitFacet*> facets;
      eval_tool_ptr->remove_facets(facets);
      delete eval_tool_ptr;
    }

    eval_tool_ptr = new FacetEvalTool(facet_list, point_list, interp_order, min_dot);
    //eval_tool_ptr = new FacetEvalTool();
    //eval_tool_ptr->replace_facets(facet_list);

    chsurf_ptr->assign_eval_tool(eval_tool_ptr);
  }

  // go through each of the curves and create CurveFacetEvalTools
  for (ii=0; ii<all_chcurves.size(); ii++)
  {
    ChollaCurve *chcurv_ptr = all_chcurves.get_and_step();    
    
    CubitStatus rv = chcurv_ptr->order_edges();

    CurveFacetEvalTool* curve_eval = chcurv_ptr->get_eval_tool();
    if (NULL != curve_eval)
    {
      DLIList<CubitFacetEdge*> eval_facets;
      curve_eval->remove_facets(eval_facets);
      delete curve_eval;
    }

    DLIList<FacetEntity*> facet_ents;
    facet_ents = chcurv_ptr->get_facet_list();
    DLIList<CubitFacetEdge*> edge_facets;
    CAST_LIST(facet_ents, edge_facets, CubitFacetEdge );
    DLIList<CubitPoint*> ordered_points;

    for( int k=0; k<edge_facets.size(); k++ )
    {
      CubitFacetEdge *tmp_edge = edge_facets.get_and_step();          
      if( k==0 )
      {
        ordered_points.append( tmp_edge->point(0) );            
        ordered_points.append( tmp_edge->point(1) );            
      }
      else        
        ordered_points.append( tmp_edge->point(1) );                      
    }     

    CurveFacetEvalTool *curv_eval_tool_ptr = new CurveFacetEvalTool;       
    edge_facets.reset();
    ordered_points.reset();
    stat = curv_eval_tool_ptr->initialize( edge_facets, ordered_points);          

    if ( stat != CUBIT_SUCCESS )
    {
      return stat;
    }
    chcurv_ptr->assign_eval_tool( curv_eval_tool_ptr );       
  }

  return stat; 
}

//===============================================================================
//Function:  determine_curve_orientation (static PUBLIC)
//Description:  Determine orientation of the curve with respect to the surface
//===============================================================================
CubitStatus ChollaEngine::determine_curve_orientation( 
  ChollaSurface *chsurf_ptr,
  ChollaCurve *chcurv_ptr,
  CubitSense & orientation )
{
  // get the first mesh edge on the curve. use it to determine the 
  // orientation of the curve with respect to the surface

  DLIList<FacetEntity*> facet_list = chcurv_ptr->get_facet_list();
  facet_list.reset();
  FacetEntity *facet_ptr = facet_list.get();
  CubitFacetEdge *edge_ptr = CAST_TO( facet_ptr, CubitFacetEdge );
  if (!edge_ptr)
    return CUBIT_FAILURE; // facet edges are not on the list - we shouldn't be here

  // get the adjacent face on this surface

  DLIList<FacetEntity*> adj_face_list;
  edge_ptr->get_parents( adj_face_list );

  DLIList<FacetEntity*> face_ptrs;
  int jj;
  for ( jj = 0; jj < adj_face_list.size(); jj++ )
  {
    FacetEntity* face_ptr = adj_face_list.get_and_step();    
    TDGeomFacet *td_gm_face = TDGeomFacet::get_geom_facet(face_ptr); 
    DLIList<ChollaSurface*> chsurf_list;
    if(td_gm_face)
    {
      td_gm_face->get_cholla_surfs( chsurf_list );
      if (chsurf_list.size())
      {
        ChollaSurface *face_chsurf_ptr = chsurf_list.get();
        if(face_chsurf_ptr == chsurf_ptr)
        {
          face_ptrs.append(face_ptr);
        }
      }
    }
  }
  if(!face_ptrs.size())
    return CUBIT_FAILURE;  // didn't find an adj face on the surface ???
    

  // if there are two, we'll just say the orientation is unknown
  // this happens when there is a non-manifold, callers are responsible for dealing with those
  if(face_ptrs.size() >= 2)
  {
    orientation = CUBIT_UNKNOWN;
    return CUBIT_SUCCESS;
  }

  FacetEntity* face_ptr = face_ptrs[0];

  // determine orientation of nodes on this mesh face

  CubitPoint *start_ptr, *end_ptr;
  chcurv_ptr->get_ends( start_ptr, end_ptr );
  end_ptr = edge_ptr->other_point( start_ptr );
  if (end_ptr == NULL)    
  {  
    return CUBIT_FAILURE;  // the edge list may not be ordered correctly??
  }
  CubitPoint *points[3];
  CubitFacet *tri_ptr = CAST_TO( face_ptr, CubitFacet );
  tri_ptr->points( points[0], points[1], points[2] );
  
  int found = 0;
  for ( jj = 0; jj < 3 && !found; jj++ )
  {
    if (points[jj] == start_ptr)
    {
      int next_jj = (jj + 1)%3;
      int prev_jj = (jj + 2)%3;
      if(points[next_jj] == end_ptr)
      {
        found = 1;
        orientation = CUBIT_FORWARD;
      }
      else if(points[prev_jj] == end_ptr)
      {
        found = 1;
        orientation = CUBIT_REVERSED;
      }
    } 
  } 
  if (!found)
    return CUBIT_FAILURE;  // couldn't determine the orientation

  return CUBIT_SUCCESS;
}

//===============================================================================
//Function:  clean_geometry (PRIVATE)
//Description:  fix the edge control points and the normals so they are conforming
//              (or non-conforming) accross curves
//===============================================================================
CubitStatus ChollaEngine::clean_geometry( 
  CubitBoolean smooth_non_manifold,
  CubitBoolean split_surfaces,
  CubitBoolean use_feature_angle,
  double mindot,
  DLIList <ChollaCurve *> &cholla_curve_list )
{
  int iedge;
  ChollaCurve *chcurv_ptr;
  FacetEntity *fedge_ptr;
  CubitFacetEdge *edge_ptr;
  DLIList <FacetEntity *> facet_list; 
  DLIList<CubitFacetEdge *> feature_edge_list;
  int icurve;
  for (icurve=0; icurve<cholla_curve_list.size(); icurve++)
  {
    chcurv_ptr = cholla_curve_list.get_and_step();
    facet_list.clean_out();
    facet_list = chcurv_ptr->get_facet_list( );
    for (iedge=0; iedge<facet_list.size(); iedge++)
    {
      fedge_ptr = facet_list.get_and_step();
      edge_ptr = CAST_TO( fedge_ptr, CubitFacetEdge );
      feature_edge_list.append( edge_ptr );
    }
  }

  return fix_geometry( smooth_non_manifold,split_surfaces, 
                       use_feature_angle,mindot, feature_edge_list );
}

//=============================================================================
//Function:  check_all_facet_orientations (PUBLIC)
//Description: check for consistent orientation of the facets and set the 
//             isBackwards flag on the facets if necessary  
//Author: sjowen
//Date: 8/14/01
//=============================================================================
CubitStatus ChollaEngine::check_all_facet_orientations(
  DLIList<CubitFacet*> &facet_list,
  CubitBoolean do_flip)
{
  CubitStatus stat = CUBIT_SUCCESS;
  //return CUBIT_SUCCESS;
  // mark facets

  int ii;
  for (ii=0; ii<facet_list.size(); ii++)
    facet_list.get_and_step()->marked(1);

  int mydebug = 0;
  if (mydebug)
  {
    dcolor(CUBIT_YELLOW);
    dfldraw(facet_list);
    dview();
  }

  // recursively loop through all facets making sure we are oriented the same
  int nfacets = 0;
  CubitFacet *facet_ptr;
  while ( nfacets < facet_list.size() )
  {
    facet_ptr = facet_list.get_and_step();
    if (facet_ptr->marked() == 1)
    {
      stat = check_facet_orientation( facet_ptr, do_flip, nfacets, mydebug );
      if (stat!= CUBIT_SUCCESS)
        return stat;
    }
  }
  return stat;
}


//=============================================================================
//Function:  check_facet_orientation (PUBLIC)
//Description: check for consistent orientation of the facets and change
//             orientation if necessary  (is NOT recursive function)
//Notes:  If do_flip is CUBIT_TRUE, then the facet topology is actually
//        changed.  Otherwise the is_backwards flag is set
//Author: sjowen
//Date: 6/3/01
//=============================================================================
CubitStatus ChollaEngine::check_facet_orientation(
  CubitFacet *start_facet, 
  CubitBoolean do_flip,
  int &nfacets, int mydebug )
{

  int jj;
  CubitFacet *adjfacet, *facet;
  CubitPoint *cpt0, *cpt1, *adjpt;
  DLIList<CubitFacet *>temp_flist;
  DLIList<CubitFacet *>adj_facet_list;
  temp_flist.append(start_facet);
  if (mydebug)
  {
    dcolor(CUBIT_RED);
    dfdraw(start_facet);
  }
  nfacets++;
  start_facet->marked( 0 );
  int adjidx; 
  while(temp_flist.size())
  {
    facet = temp_flist.pop();
    if (facet->marked() == 0)
    {
      for (jj=0; jj<3; jj++)
      {
        cpt0 = facet->point(jj);
        cpt1 = facet->point((jj+1)%3);
        adj_facet_list.clean_out();
        facet->shared_facets(cpt0, cpt1, adj_facet_list);
        if(adj_facet_list.size()==1)
        {
          adjfacet = adj_facet_list.get();
          if (adjfacet->marked() == 1)
          {
            adjidx = adjfacet->point_index( cpt1 ); 
            if (adjidx == -1)
            {
              return CUBIT_FAILURE;
            }
            adjpt = adjfacet->point((adjidx+1)%3);
            if (adjpt != cpt0)
            {        
              // the orientation is not consistent -- flip the orientation of the 
              // adjacent facet

              // an assertion here means that the facets really don't share an
              // edge.  Should not happen!
              
              if (adjfacet->point((adjidx+2)%3) != cpt0)
              {
                return CUBIT_FAILURE;
              }
       
              if (do_flip)
              {
                adjfacet->flip();
              }
              else
              {
                if(adjfacet->is_backwards())
                  adjfacet->is_backwards( 0 );
                else
                  adjfacet->is_backwards( 1 );
              }
              if(mydebug){
                adjfacet->debug_draw(CUBIT_WHITE);
              }
            }

            // put its neighbor on the list to be processed later

            if (adjfacet->marked() == 1)
            {
              temp_flist.append( adjfacet );
              adjfacet->marked( 0 );
              if(mydebug)
              {
                dfdraw(adjfacet);
              }
              nfacets++;
            }
          }
        }
      }
    }
  }
  if(mydebug)
    GfxDebug::mouse_xforms();
  return CUBIT_SUCCESS;
}


#ifdef ALPHA_CABLE
//=============================================================================
//Function:  split_surface_at_edges (PUBLIC)
//Description: Given a list of edges that divides a surface into two, this
// function splits the surface and modifies the engine to represent the
// new geometry.
//Author: mbrewer
//Date: 4/22/2004
//=============================================================================
CubitStatus ChollaEngine::split_surface_at_edges(
  ChollaSurface* owning_surface,
  DLIList<CubitFacetEdge *> &feature_edge_list,
  DLIList<ChollaCurve*> &cholla_curves,
  DLIList<ChollaSurface*> &cholla_surfaces)
{
  int i=0;
  DLIList<ChollaSurface*> original_surfaces = chollaSurfaceList;
  DLIList<ChollaCurve*> original_curves = chollaCurveList;
  //mark edges as feature
  CubitFacetEdge* edge_ptr=NULL;
  for(i=0;i<feature_edge_list.size();++i){
    edge_ptr=feature_edge_list.get_and_step();
    TDGeomFacet::add_geom_facet(edge_ptr, -1); 
    TDGeomFacet* td_gm_edge = TDGeomFacet::get_geom_facet(edge_ptr);
    edge_ptr->set_as_feature();
  }
   //Michael temp
  DLIList<FacetEntity*> *surf_facet_list=owning_surface->get_facet_list_ptr();
  int surface_id = owning_surface->get_id();
  for(i=0;i<surf_facet_list->size();++i){
    FacetEntity* face_ptr = surf_facet_list->get_and_step();
    TDGeomFacet *td_gm_face = TDGeomFacet::get_geom_facet(face_ptr);
    td_gm_face->set_hit_flag(surface_id);
  }
  
  CubitStatus rv = make_features( feature_edge_list, CUBIT_FALSE );
  if (rv != CUBIT_SUCCESS)
    return rv;
  
    // split up the surface
  rv = owning_surface->split_surface( cholla_surfaces );
  if (rv != CUBIT_SUCCESS)
    return rv;

  for ( i = cholla_surfaces.size(); i > 0; i-- )
  {
    ChollaSurface *chsurf_ptr = cholla_surfaces.get_and_step();
    chsurf_ptr->update_boundary_tool_data();
    if(chsurf_ptr != owning_surface){
      chollaSurfaceList.append(chsurf_ptr);
     
    }
  }
  rv = create_surface_boundaries( chollaSurfaceList, chollaCurveList,
                                  CUBIT_FALSE, 0.0);
  if ( !rv )
    return rv;

  rv = create_curve_boundaries( chollaCurveList, chollaPointList );
  if ( !rv )
    return rv;

  rv = build_eval_tools( chollaSurfaceList,
                         chollaCurveList,
                         0, CUBIT_FALSE, 
                         0, CUBIT_FALSE, CUBIT_FALSE);
  
  
  cholla_surfaces = chollaSurfaceList;
  cholla_curves = chollaCurveList;

  cholla_surfaces -= original_surfaces;
  cholla_curves -= original_curves;
  
  return rv;
}
#endif//cable

//=============================================================================
//Function:  make_features (PRIVATE)
//Description: add additional edges and points where feature angle requires
//             a break in the surface
//Author: sjowen
//Date:  4/31/01
//=============================================================================
CubitStatus ChollaEngine::make_features(
  DLIList<CubitFacetEdge *> &feature_edge_list,
  CubitBoolean split_surfaces )
{
  CubitStatus rv = CUBIT_SUCCESS;
  int ii, jj;
  CubitFacetEdge *edge_ptr;
  CubitPoint *point_ptr;
  int mydebug = 0;
  if(mydebug)
    GfxDebug::clear();
  for (ii=0; ii<feature_edge_list.size(); ii++)
  {
    edge_ptr = feature_edge_list.get_and_step();
    if(mydebug)
      edge_ptr->debug_draw(CUBIT_YELLOW);
    for (jj=0; jj<2; jj++)
    {
      point_ptr = edge_ptr->point( jj );
      point_ptr->marked( 0 );
    }
    edge_ptr->set_flag( 1 );  // mark edges as being feature edges
  }

  if (mydebug)
  {
    GfxDebug::mouse_xforms();
    dcolor(CUBIT_WHITE);
    deldraw(feature_edge_list);
    dview();
  }

  for (ii=0; ii<feature_edge_list.size(); ii++)
  {
    edge_ptr = feature_edge_list.get_and_step();
    for (jj=0; jj<2; jj++)
    {
      point_ptr = edge_ptr->point( jj );
      if (0 == point_ptr->marked())
      {
        point_ptr->marked(1);
        if (split_surfaces)
          rv = crack_open_point( point_ptr, mydebug );
        else
          rv = insert_discontinuity_at_point( point_ptr );
        if (rv != CUBIT_SUCCESS)
          return rv;
      }
    }
  }

  for (ii=0; ii<feature_edge_list.size(); ii++)
  {
    feature_edge_list.get_and_step()->set_flag( 0 );
  }
  return rv;
}

//=============================================================================
//Function:  crack_open_point (PRIVATE)
//Description: split the surface facets at a point if required.  Create
//             new points and edges as needed
//Author: sjowen
//Date:  4/31/01
//=============================================================================
CubitStatus ChollaEngine::crack_open_point(
  CubitPoint *point_ptr,  
  int mydebug )
{
  CubitStatus rv = CUBIT_SUCCESS;

  // get facets at this point and mark them all ready for recursion function
  // get_facets_at_point

  int ii, jj;
  DLIList<CubitFacet *> facet_list;
  point_ptr->facets( facet_list );
  CubitFacet *facet_ptr, *adj_facet_ptr;
  for (ii=0; ii<facet_list.size(); ii++)
  {
    facet_ptr = facet_list.get_and_step();
    facet_ptr->marked( 0 );
  }

  // go through and find groups of adjacent facets that will form part of
  // new surfaces

  int totfacets = facet_list.size();
  int numfacets = 0;
  while(numfacets < totfacets)
  {
    facet_list.clean_out();
    point_ptr->facets( facet_list );
    facet_ptr = facet_list.get();
    DLIList<CubitFacet *> adj_facet_list;
    DLIList<CubitFacetEdge *>feature_edge_list;

    // make a list of all facets adjacet this point that are bounded
    // by feature edges.  Start with facet_ptr and recurse until
    // we bump into feature edges at the point.

    rv = get_facets_at_point( point_ptr, facet_ptr, adj_facet_list, 
                              feature_edge_list );
    numfacets += adj_facet_list.size();

    // if the number of facets for this feature are the same as the number of
    // facets at this point, then we are done with this point.  Otherwise
    // create new point and edges 

    if (numfacets < totfacets)
    {
      // create a new point and update the facets

      CubitPoint *new_point_ptr = 
         new CubitPointData( point_ptr->coordinates() );
      TDGeomFacet::add_geom_facet( new_point_ptr , -1);
      for (ii=0; ii<adj_facet_list.size(); ii++)
      {
        adj_facet_ptr = adj_facet_list.get_and_step();

        // facet modifications can only occur with native CUBIT facet data
        // if we are using someone else's data (ie. Sierra) then fail.
        // (Should not be using split_surfaces = TRUE)

        CubitFacetData *cfd_ptr = CAST_TO( adj_facet_ptr, CubitFacetData );
        if (!cfd_ptr)
          return CUBIT_FAILURE;
        new_point_ptr->add_facet( adj_facet_ptr );
        point_ptr->remove_facet( adj_facet_ptr );
        int pidx = adj_facet_ptr->point_index( point_ptr );
        cfd_ptr->set_point( new_point_ptr, pidx );
      }

      // The existing and the new points in some cases have to maintain 
      // the same normal...  so they need to keep track of
      // each other.  They will likely wind up on different surfaces
      // and we don't want to go search for them later

      TDGeomFacet *td_existing_point = TDGeomFacet::get_geom_facet( point_ptr );
      TDGeomFacet *td_new_point = TDGeomFacet::get_geom_facet( new_point_ptr );
      DLIList<CubitPoint *> partner_point_list;

      // add the existing point's partner list to the new point

      CubitPoint *partner_point_ptr;
      TDGeomFacet *td_partner;
      td_existing_point->get_partner_points( partner_point_list );
      for (ii=0; ii<partner_point_list.size(); ii++)
      {
        partner_point_ptr = partner_point_list.get_and_step();
        td_partner = TDGeomFacet::get_geom_facet( partner_point_ptr );
          //added mbrewer:  the check for a pointer here shouldn't
          // be necessary, but I have a case where it fails.  I don't
          // see where the logic is wrong.
        if(!td_partner){
          TDGeomFacet::add_geom_facet( partner_point_ptr , -1);
          td_partner = TDGeomFacet::get_geom_facet( partner_point_ptr );
        }
        td_partner->add_partner_point( new_point_ptr );
        td_new_point->add_partner_point( partner_point_ptr );
      }

      // Add each other as a partner

      td_existing_point->add_partner_point( new_point_ptr );
      td_new_point->add_partner_point( point_ptr );

      // create new edges at the point

      CubitFacetEdgeData *cfed_ptr;
      for (ii=0; ii<feature_edge_list.size(); ii++)
      {
        CubitFacetEdge *edge_ptr = feature_edge_list.get_and_step();
        cfed_ptr = CAST_TO( edge_ptr, CubitFacetEdgeData );
        if (!cfed_ptr)
          return CUBIT_FAILURE;
        rv = CUBIT_FAILURE;
        for (jj=0; jj<adj_facet_list.size() && rv == CUBIT_FAILURE; jj++)
        {
          rv = cfed_ptr->remove_facet( adj_facet_list.get_and_step() );
        }
        if (rv != CUBIT_SUCCESS)
          return rv;
        CubitPoint *p1 = edge_ptr->other_point( point_ptr );

        // create a new edge - the edge constructor will update
        // the facet and its edge use.

        CubitFacetEdge *new_edge_ptr =
          new CubitFacetEdgeData( new_point_ptr, p1 );
        TDGeomFacet::add_geom_facet( new_edge_ptr, -1 );
        new_edge_ptr->set_as_feature();
        if (mydebug)
        {
          dcolor(CUBIT_RED);
          dedraw(new_edge_ptr);
        }

        // The existing and the new edges have to maintain the same 
        // control point locations, so they need to keep track of
        // each other.  They will likely wind up on different surfaces
        // and we don't want to go search for them later

        TDGeomFacet *td_existing_edge = TDGeomFacet::get_geom_facet( edge_ptr );
        TDGeomFacet *td_new_edge = TDGeomFacet::get_geom_facet( new_edge_ptr );
        DLIList<CubitFacetEdge *> partner_edge_list;

        // add the existing edge's partner list to the new edge

        CubitFacetEdge *partner_edge_ptr;
        td_existing_edge->get_partner_edges( partner_edge_list );
        for (jj=0; jj<partner_edge_list.size(); jj++)
        {
          partner_edge_ptr = partner_edge_list.get_and_step();
          td_partner = TDGeomFacet::get_geom_facet( partner_edge_ptr );
          td_partner->add_partner_edge( new_edge_ptr );
          td_new_edge->add_partner_edge( partner_edge_ptr );
        }

        // Add each other as a partner

        td_existing_edge->add_partner_edge( new_edge_ptr );
        td_new_edge->add_partner_edge( edge_ptr );
      }

      // update the other edges with the new point

      for (ii=0; ii<adj_facet_list.size(); ii++)
      {
        adj_facet_ptr = adj_facet_list.get_and_step();
        for (jj=0; jj<3; jj++)
        {
          CubitFacetEdge *edge_ptr = adj_facet_ptr->edge( jj );
          cfed_ptr = CAST_TO( edge_ptr, CubitFacetEdgeData );
          if (!cfed_ptr)
            return CUBIT_FAILURE;
          if (cfed_ptr->point( 0 ) == point_ptr)
            cfed_ptr->set_point( new_point_ptr, 0 );
          if (cfed_ptr->point( 1 ) == point_ptr)
            cfed_ptr->set_point( new_point_ptr, 1 );
        }
      }
    }
  }

  return rv;
}

//=============================================================================
//Function:  insert_discontinuity_at_point (PRIVATE)
//Description: does the same as crack_open_point, but does not create new 
//             facet entities at surface boundaries.  Instead, it creates 
//             facet boundary tooldatas to hold the additional control point 
//             and normal data
//Author: sjowen
//Date:  4/31/01
//=============================================================================
CubitStatus ChollaEngine::insert_discontinuity_at_point(
  CubitPoint *point_ptr )
{
  CubitStatus rv = CUBIT_SUCCESS;

  // get facets at this point and mark them all ready for recursion function
  // get_facets_at_point

  int ii;
  DLIList<CubitFacet *> facet_list;
  point_ptr->facets( facet_list );
  CubitFacet *facet_ptr = NULL;
  for (ii=0; ii<facet_list.size(); ii++)
  {
    facet_ptr = facet_list.get_and_step();
    facet_ptr->marked( 0 );
  }

  int mydebug=0;
  if (mydebug)
  {
    dfldraw(facet_list);
    dview();
    mydebug = mydebug;
  }

  // go through and find groups of adjacent facets that will form part of
  // new surfaces

  TDFacetBoundaryPoint *td_fbp = 
    TDFacetBoundaryPoint::get_facet_boundary_point( point_ptr );
  int totfacets = facet_list.size();
  int numfacets = 0;
  while(numfacets < totfacets)
  {
    ii = 0;
    int found = 0;
    while(!found && ii < totfacets )
    {
      ii++;
      facet_ptr = facet_list.get_and_step();
      if (facet_ptr != NULL)
        found = 1;
    }
    assert(found);
    DLIList<CubitFacet *> adj_facet_list;
    DLIList<CubitFacetEdge *>feature_edge_list;

    // make a list of all facets adjacet this point that are bounded
    // by feature edges.  Start with facet_ptr and recurse until
    // we bump into feature edges at the point.

    rv = get_facets_at_point( point_ptr, facet_ptr, adj_facet_list, 
                              feature_edge_list );

    if (rv != CUBIT_SUCCESS)
      return rv;
    numfacets += adj_facet_list.size();
    for (ii=0; ii<adj_facet_list.size(); ii++)
    {
      facet_list.move_to( adj_facet_list.get_and_step() );
      facet_list.change_to( NULL );
    }

    // create a new boundary point tooldata if needed and add the facets
    // associated with this surface

    if (td_fbp == NULL)
    {
      TDFacetBoundaryPoint::add_facet_boundary_point( point_ptr );
      td_fbp = TDFacetBoundaryPoint::get_facet_boundary_point( point_ptr );
    }
    td_fbp->add_surf_facets( adj_facet_list );

    // create new boundary edge tooldatas at the point
    for (ii=0; ii<feature_edge_list.size(); ii++)
    {
      CubitFacetEdge *edge_ptr = feature_edge_list.get_and_step();

      TDFacetBoundaryEdge *td_fbe = 
        TDFacetBoundaryEdge::get_facet_boundary_edge( edge_ptr );
      if (td_fbe == NULL)
      {
        TDFacetBoundaryEdge::add_facet_boundary_edge( edge_ptr );
        td_fbe = TDFacetBoundaryEdge::get_facet_boundary_edge( edge_ptr );
      }
      td_fbe->add_surf_facet( adj_facet_list );
    }
  }

  return rv;
}

//=============================================================================
//Function:  get_facets_at_point (PRIVATE)
//Description: make a list of all faces adjacent to a single point that are 
//             bounded by feature edges - recursive function
//Author: sjowen
//Date:  5/2/01
//=============================================================================
CubitStatus ChollaEngine::get_facets_at_point(
  CubitPoint *point_ptr,
  CubitFacet *facet_ptr,
  DLIList<CubitFacet *> &facet_list,
  DLIList<CubitFacetEdge *> &feature_edge_list)
{
  CubitStatus rv = CUBIT_SUCCESS;
  CubitFacetEdge *edge_ptr;
  CubitFacet *adj_facet;
  facet_list.append(facet_ptr);
  facet_ptr->marked( 1 );

  for (int ii=0; ii<3; ii++)
  {
    edge_ptr = facet_ptr->edge( ii );
    if (edge_ptr->contains( point_ptr ))
    {
      if (!edge_ptr->is_feature())
      {
        if (edge_ptr->num_adj_facets() >2)
          return CUBIT_FAILURE;
        adj_facet = edge_ptr->other_facet( facet_ptr );
        if (adj_facet == NULL)
          return CUBIT_FAILURE;
        if (adj_facet->marked() == 0)
        {
          rv = get_facets_at_point( point_ptr, adj_facet, 
                                    facet_list, feature_edge_list );
          if (rv != CUBIT_SUCCESS)
            return rv;
        }
      }
      else
      {
        if (edge_ptr->num_adj_facets() > 1)
        {
          feature_edge_list.append_unique( edge_ptr );
        }
      }
    }
  }
  return rv;
}

//===============================================================================
//Function:  fix_geometry (PUBLIC)
//Description:  fix the edge control points and the normals so they are conforming
//              (or non-conforming) accross curves
//===============================================================================
CubitStatus ChollaEngine::fix_geometry(
  CubitBoolean smooth_non_manifold,
  CubitBoolean split_surfaces,
  CubitBoolean use_feature_angle,
  double mindot,
  DLIList <CubitFacetEdge *> &feature_edge_list )
{
  CubitStatus stat = CUBIT_SUCCESS;
  use_feature_angle = use_feature_angle;

  // unmark all the edges and points

  CubitFacetEdge *edge_ptr; 
  int iedge;
  int flag = (smooth_non_manifold) ? 0 : 1;
  for (iedge = 0; iedge < feature_edge_list.size(); iedge++)
  {
    edge_ptr = feature_edge_list.get_and_step();
    edge_ptr->set_flag( flag );
    edge_ptr->point( 0 )->marked( flag );
    edge_ptr->point( 1 )->marked( flag );
  }

  // for non-manifold edges, find sets of adjacent facets that 
  // we can smooth across.  If adjacent facets do not satisfy the
  // feature angle criteria, then they should have C1 continuity
  // across the edge.  To do this, the normals at the edge's
  // points on partner edges will need to be the same.  This next
  // piece of code potentially alters the point normals.  As a result
  // edge and facet control points will need to be updated.

  DLIList<CubitFacet *>update_facet_list;
  if (smooth_non_manifold)
  {
    DLIList <CubitPoint *> changed_points;
    for (iedge = 0; iedge < feature_edge_list.size(); iedge++)
    {
      edge_ptr = feature_edge_list.get_and_step();
      if (edge_ptr->get_flag() == 1)
        continue;
      edge_ptr->set_flag( 1 );
      if (split_surfaces)
        stat = fix_split_non_manifold_edge( edge_ptr, mindot, 
                                            changed_points );
      else
        stat = fix_non_manifold_edge( edge_ptr, mindot, 
                                      changed_points );
      if (stat != CUBIT_SUCCESS)
        return stat;
    }
    if (changed_points.size() > 0)
      stat = update_edges_at_points( split_surfaces,
                                     changed_points, 
                                     update_facet_list,
                                     mindot );
    if (stat != CUBIT_SUCCESS)
      return stat;
  }

  // For all cases we need to ensure that partner edges share the same
  // control, point locations.  This piece of code goes through and
  // computes an average location for control points on partner edges
  // and sets them back on the edge

  for (iedge = 0; iedge < feature_edge_list.size(); iedge++)
  {
    edge_ptr = feature_edge_list.get_and_step();
    if (edge_ptr->get_flag() == 2)
      continue;
    edge_ptr->set_flag( 2 );
    if (split_surfaces)
      stat = fix_split_edge( edge_ptr, update_facet_list );
    else
      stat = fix_edge( edge_ptr, update_facet_list );
    if (stat != CUBIT_SUCCESS)
      return stat;
  }

  // since we have potentially modified the point normals and edge 
  // control points, we need to fix the facet control points (interior
  // to the facet) The update_facet_list contains a non-unique list of
  // facets that will need to be recomputed.

  int ii;
  CubitFacet *facet_ptr;
  for (ii=0; ii<update_facet_list.size(); ii++)
    update_facet_list.get_and_step()->marked( 0 );

  for (ii=0; ii<update_facet_list.size(); ii++)
  {
    facet_ptr = update_facet_list.get_and_step();
    if (facet_ptr->marked() == 0)
    {
      facet_ptr->marked( 1 );
      stat = FacetEvalTool::init_bezier_facet( facet_ptr );
      if (stat != CUBIT_SUCCESS)
        return stat;
    }
  }

  for (ii=0; ii<update_facet_list.size(); ii++)
    update_facet_list.get_and_step()->marked( 0 );

  return stat;
}

//===============================================================================
//Function:  update_edges_at_points (PRIVATE)
//Description: recompute the edge control points at all facets sharing
//             points in the list
//Date: 05/01
//===============================================================================
CubitStatus ChollaEngine::update_edges_at_points(
  CubitBoolean split_surfaces,        // we are splitting facets to create geometry
  DLIList<CubitPoint *> &point_list,  // points whose edges will be updated
  DLIList<CubitFacet *> &facet_update_list,  // append facets adjacent edges that
                                             // will be updated later
  double mindot )
{
  split_surfaces = split_surfaces;
  CubitStatus stat = CUBIT_SUCCESS;
  CubitPoint *point_ptr;
  CubitFacetEdge *edge_ptr;
  DLIList <CubitFacetEdge *> edge_list;
  int ii;
  for (ii=0; ii<point_list.size(); ii++)
  {
    point_ptr = point_list.get_and_step();
    point_ptr->edges( edge_list );
    point_ptr->facets( facet_update_list );
  }
  for (ii=0; ii<edge_list.size(); ii++)
  {
    edge_ptr = edge_list.get_and_step();
    edge_ptr->set_flag( 0 );
  }
  for (ii=0; ii<edge_list.size(); ii++)
  {
    edge_ptr = edge_list.get_and_step();
    if (edge_ptr->get_flag() == 0)
    {
      edge_ptr->set_flag( 1 );
      stat = FacetEvalTool::init_bezier_edge( edge_ptr, mindot );
    }
  }
  return stat;
}

//===============================================================================
//Function:  fix_edge (PRIVATE)
//Description: Enforce the condition where all partner edges to this edge
//             share the same control points.  Compute the average of the
//             partner control points and set on the edge
//Date: 05/01
//===============================================================================
CubitStatus ChollaEngine::fix_edge( CubitFacetEdge *edge_ptr,
                                   DLIList<CubitFacet *> &update_facet_list)
{
  if (edge_ptr->num_adj_facets() == 1)
    return CUBIT_SUCCESS;
  TDFacetBoundaryEdge *td_fbe =
    TDFacetBoundaryEdge::get_facet_boundary_edge( edge_ptr );
  if (!td_fbe)
    return CUBIT_FAILURE;

  edge_ptr->facets( update_facet_list );
  CubitStatus stat = td_fbe->merge_control_points();
  return stat;
}

//===============================================================================
//Function:  fix_split_edge (PRIVATE)
//Description: Enforce the condition where all partner edges to this edge
//             share the same control points.  Compute the average of the
//             partner control points and set on the edge
//Date: 05/01
//===============================================================================
CubitStatus ChollaEngine::fix_split_edge( CubitFacetEdge *edge_ptr,
                                   DLIList<CubitFacet *> &update_facet_list)
{

  // we only need to do this if the edge has any partners.  Otherwise return

  TDGeomFacet *td = TDGeomFacet::get_geom_facet( edge_ptr );
  if (td->num_partner_edges() == 0)
    return CUBIT_SUCCESS;
 
  // get the control points for this edge (oriented with
  // respect to point(0) on the edge)--- This means that
  // the rest of the partners to this edge must also
  // provide control points with respect to the same orientation

  CubitVector control_points[3];
  CubitPoint *start_ptr = edge_ptr->point(0);
  edge_ptr->get_control_points( start_ptr, control_points );

  // get all this edge's partners

  int ii;
  CubitVector partner_control_points[3];
  CubitFacetEdge *partner_edge_ptr;
  CubitPoint *partner_start_ptr;
  DLIList <CubitFacetEdge *> partner_edge_list;
  td->get_partner_edges( partner_edge_list );
  for (ii=0; ii<partner_edge_list.size(); ii++)
  {
    partner_edge_ptr = partner_edge_list.get_and_step();
    partner_edge_ptr->set_flag( 2 );

    // determine the same orientation as the first edge

    partner_start_ptr = partner_edge_ptr->point(0);
    TDGeomFacet *td_point = TDGeomFacet::get_geom_facet( partner_start_ptr );
    if (!td_point->is_partner( start_ptr ))
    {
      partner_start_ptr = partner_edge_ptr->point(1);
      td_point = TDGeomFacet::get_geom_facet( partner_start_ptr );
      if (!td_point->is_partner( start_ptr ))
      {
        PRINT_ERROR("Topology error generating facet geometry\n");
        return CUBIT_FAILURE;
      }
    }

    // get the partner's control points and add to the current
    // totals

    partner_edge_ptr->get_control_points( partner_start_ptr, 
                                      partner_control_points );
    control_points[0] += partner_control_points[0];
    control_points[1] += partner_control_points[1];
    control_points[2] += partner_control_points[2];

    // add this edge's adjacent facets to the update facet list

    partner_edge_ptr->facets( update_facet_list );
  }
  edge_ptr->facets( update_facet_list );

  // average the control point locations

  int nedges = partner_edge_list.size() + 1;
  control_points[0] /= (double)nedges;
  control_points[1] /= (double)nedges;
  control_points[2] /= (double)nedges;

  // set the new control points back on the edges

  edge_ptr->set_control_points( start_ptr, control_points );

  TDGeomFacet *td_point;
  for (ii=0; ii<partner_edge_list.size(); ii++)
  {
    partner_edge_ptr = partner_edge_list.get_and_step();

    // determine the same orientation as the first edge

    partner_start_ptr = partner_edge_ptr->point(0);
    td_point = TDGeomFacet::get_geom_facet( partner_start_ptr );
    if (!td_point->is_partner( start_ptr ))
    {
      partner_start_ptr = partner_edge_ptr->point(1);
      td_point = TDGeomFacet::get_geom_facet( partner_start_ptr );
      if (!td_point->is_partner( start_ptr ))
      {
        PRINT_ERROR("Topology error generating facet geometry\n");
        return CUBIT_FAILURE;
      }
    }
    partner_edge_ptr->set_control_points( partner_start_ptr,
                                           control_points );
  }
  return CUBIT_SUCCESS;
}

//===============================================================================
//Function:  fix_non_split_manifold_edge (PRIVATE)
//Description: Enforce the condition where partner edges on non-manifold edges
//             with >2 adjacent facets that do not satisfy feature criteria
//             will have the same normal accross the edge
//Date: 05/01
//===============================================================================
CubitStatus ChollaEngine::fix_split_non_manifold_edge(
  CubitFacetEdge *edge_ptr,
  double min_dot,
  DLIList <CubitPoint *> &changed_points )
{
  // first check if this edge has 3 or more adjacent facets.  Return now if not

  CubitStatus stat = CUBIT_SUCCESS;
  TDGeomFacet *td = TDGeomFacet::get_geom_facet( edge_ptr );
  int nedges = td->num_partner_edges() + 1;
  if (nedges < 3)
    return stat;

  // make an array of partner edges

  int *partner0 = new int [nedges];
  int *partner1 = new int [nedges];
  double *partner_dot = new double [nedges];
  int num_partners = 0;
  CubitFacetEdge **edges = new CubitFacetEdge * [nedges];
  if (!edges || !partner0 || !partner1 || !partner_dot)
    return CUBIT_FAILURE;
  edges[0] = edge_ptr;
  DLIList<CubitFacetEdge *> partner_list;
  td->get_partner_edges( partner_list );
  int ii;
  for (ii=0; ii< partner_list.size(); ii++)
  {
    edges[ii+1] = partner_list.get_and_step();
    edges[ii+1]->set_flag( 1 );
  }

  
  // check each facet against every other facet at this edge to
  // see if we need to enforce continuity.  There may be several
  // cases where the feature angle criteria is (not) met.  Find
  // the pair that has the greatest dot product between them.
  // If more than 1 pair has the same dot product and meets criteria
  // then it arbitrarily uses the first one.

  int jj, kk;
  for (ii=0; ii<nedges-1; ii++)
  {
    CubitFacet *facet_ptr = edges[ii]->adj_facet( 0 );
    CubitVector normal = facet_ptr->normal();
    normal.normalize();
    double best_dot = -1.0;
    for (jj=ii+1; jj<nedges; jj++)
    {
      CubitFacet *afacet_ptr = edges[jj]->adj_facet( 0 );
      CubitVector anormal = afacet_ptr->normal();
      anormal.normalize();
      double dot = fabs(normal % anormal);
      if (dot > min_dot && dot > best_dot)
      {

        // check to see if this partner pair has not already
        // been used in another partnership.  If so, then check
        // to see if this match is better.  Otherwise, ignore
        // this match.

        int found = 0;
        for (kk=0; kk<num_partners; kk++)
        {
          if (partner0[kk] == ii || partner0[kk] == jj ||
              partner1[kk] == ii || partner1[kk] == jj)
          {
            found = 1;
            if (dot > partner_dot[kk])
            {
              partner_dot[kk] = dot;
              partner0[kk] = ii;
              partner1[kk] = jj;
            }
          }
        }
        if (!found)
        {
          partner_dot[num_partners] = dot;
          partner0[num_partners] = ii;
          partner1[num_partners] = jj;
          num_partners++;
        }
        best_dot = dot;
      }
    }
  }

  // we found a partner facet(s) that needs to have C1 continuity
  // across this facet.  Merge the normals on the edge so they
  // are the same

  for (ii=0; ii<num_partners; ii++)
  {
    CubitFacetEdge *edge0 = edges[partner0[ii]];
    CubitFacetEdge *edge1 = edges[partner1[ii]];
    for (jj=0; jj<2; jj++)
    {
      CubitPoint *pt0 = edge0->point(jj);

      // find its point partner on the partner edge

      CubitPoint *pt1 = edge1->point( 0 );
      td = TDGeomFacet::get_geom_facet( pt0 );
      if (!td->is_partner( pt1 ))
      {
        pt1 = edge1->point( 1 );
        if (!td->is_partner( pt1 ))
        {
          PRINT_ERROR("Internal error defining continuity across non-manifold edges\n");
          delete [] edges;
          return CUBIT_FAILURE;
        }
      }
      if (pt0->marked() == 0 || pt1->marked() == 0)
      {
        pt0->marked( 1 );
        pt1->marked( 1 );
        changed_points.append( pt0 );
        changed_points.append( pt1 );
        stat = merge_normals( pt0, pt1 );
        if (stat != CUBIT_SUCCESS)
          return stat;
      }
    }
  }
  
  delete [] edges;
  delete [] partner0;
  delete [] partner1;
  delete [] partner_dot;

  return stat;
}

//===============================================================================
//Function:  fix_non_manifold_edge (PRIVATE)
//Description: Enforce the condition where partner edges on non-manifold edges
//             with >2 adjacent facets that do not satisfy feature criteria
//             will have the same normal accross the edge
//Date: 05/01
//===============================================================================
CubitStatus ChollaEngine::fix_non_manifold_edge(
  CubitFacetEdge *edge_ptr,
  double min_dot,
  DLIList <CubitPoint *> &changed_points )
{
  // first check if this edge has 3 or more adjacent facets.  Return now if not

  CubitStatus stat = CUBIT_SUCCESS;
  int nfacets = edge_ptr->num_adj_facets();
  if (nfacets < 3)
    return stat;

  TDFacetBoundaryEdge *td_fbe =
    TDFacetBoundaryEdge::get_facet_boundary_edge( edge_ptr );
  if (td_fbe == NULL)
  {
    return CUBIT_FAILURE;  // there is supposed to be a tool data here!
  }

  // make an array of partner edges

  int *partner0 = new int [nfacets];
  int *partner1 = new int [nfacets];
  double *partner_dot = new double [nfacets];
  int num_partners = 0;
  CubitFacet **facets = new CubitFacet * [nfacets];
  if (!facets || !partner0 || !partner1 || !partner_dot)
    return CUBIT_FAILURE;
  DLIList <CubitFacet *>adj_facets;
  edge_ptr->facets( adj_facets );
  int ii;
  for (ii=0; ii< adj_facets.size(); ii++)
  {
    facets[ii] = adj_facets.get_and_step();
    facets[ii]->marked( 1 );
  }

  
  // check each facet against every other facet at this edge to
  // see if we need to enforce continuity.  There may be several
  // cases where the feature angle criteria is (not) met.  Find
  // the pair that has the greatest dot product between them.
  // If more than 1 pair has the same dot product and meets criteria
  // then it arbitrarily uses the first one.

  int jj, kk;
  for (ii=0; ii<nfacets-1; ii++)
  {
    CubitFacet *facet_ptr = facets[ii];
    CubitVector normal = facet_ptr->normal();
    normal.normalize();
    double best_dot = -1.0;
    for (jj=ii+1; jj<nfacets; jj++)
    {
      CubitFacet *afacet_ptr = facets[jj];
      CubitVector anormal = afacet_ptr->normal();
      anormal.normalize();
      double dot = fabs(normal % anormal);
      if (dot > min_dot && dot > best_dot)
      {

        // check to see if this partner pair has not already
        // been used in another partnership.  If so, then check
        // to see if this match is better.  Otherwise, ignore
        // this match.

        int found = 0;
        for (kk=0; kk<num_partners; kk++)
        {
          if (partner0[kk] == ii || partner0[kk] == jj ||
              partner1[kk] == ii || partner1[kk] == jj)
          {
            found = 1;
            if (dot > partner_dot[kk])
            {
              partner_dot[kk] = dot;
              partner0[kk] = ii;
              partner1[kk] = jj;
            }
          }
        }
        if (!found)
        {
          partner_dot[num_partners] = dot;
          partner0[num_partners] = ii;
          partner1[num_partners] = jj;
          num_partners++;
        }
        best_dot = dot;
      }
    }
  }

  // we found a partner facet(s) that needs to have C1 continuity
  // across this facet.  Merge the normals on the edge so they
  // are the same

  if (num_partners > 0)
  {
    CubitPoint *pt0 = edge_ptr->point(0);
    CubitPoint *pt1 = edge_ptr->point(1);
    TDFacetBoundaryPoint *td_fbp0 =
      TDFacetBoundaryPoint::get_facet_boundary_point( pt0 );
    TDFacetBoundaryPoint *td_fbp1 =
      TDFacetBoundaryPoint::get_facet_boundary_point( pt1 );
    if (!td_fbp0 || !td_fbp1)
      return CUBIT_FAILURE;

    for (ii=0; ii<num_partners; ii++)
    {
      td_fbp0->merge_normals( facets[partner0[ii]], facets[partner1[ii]] );
      td_fbp1->merge_normals( facets[partner0[ii]], facets[partner1[ii]] );
    }
    pt0->marked( 1 );
    pt1->marked( 1 );
    changed_points.append( pt0 );
    changed_points.append( pt1 );
  }
  
  delete [] facets;
  delete [] partner0;
  delete [] partner1;
  delete [] partner_dot;

  return stat;
}

//===============================================================================
//Function:  merge_normals (PRIVATE)
//Description: Enforce continuity across facets by merging the normals at two
//             points.  Assumes the points are coincident.
//Date: 05/01
//===============================================================================
CubitStatus ChollaEngine::merge_normals(
  CubitPoint *pt0,
  CubitPoint *pt1)
{

  DLIList<CubitFacet*> adj_facet_list0;
  pt0->facets(adj_facet_list0);
  if (adj_facet_list0.size() == 0)
    return CUBIT_FAILURE;

  CubitVector avg_normal(0.0e0, 0.0e0, 0.0e0);
  double totangle = 0.0e0;
  double angle;
  CubitFacet *facet;
  int j;

  // weight the normal by the spanning angle at the point

  for (j = 0; j < adj_facet_list0.size(); j++)
  {
    facet = adj_facet_list0.get_and_step();
    angle = facet->angle( pt0 );
    facet->weight( angle );
    totangle += angle;
  }
  DLIList<CubitFacet*> adj_facet_list1;
  pt1->facets(adj_facet_list1);
  if (adj_facet_list1.size() == 0)
    return CUBIT_FAILURE;
  for (j = 0; j < adj_facet_list1.size(); j++)
  {
    facet = adj_facet_list1.get_and_step();
    angle = facet->angle( pt1 );
    facet->weight( angle );
    totangle += angle;
  }

  // computed weighted normal

  CubitVector normal;
  for (j = 0; j < adj_facet_list0.size(); j++)
  {
    facet = adj_facet_list0.get_and_step();
    normal = facet->normal();
    normal.normalize();
    avg_normal += (facet->weight() / totangle) * normal;
  }
  
  // orientation of facets may be opposite on opposing surfaces.
  // Check for this case and correct of necessary
  
  CubitVector norm0, norm1;
  norm0 = pt0->normal();
  norm0.normalize();
  norm1 = pt1->normal();
  norm1.normalize();
  double dot = norm0 % norm1;
  double sign = 1.0;
  if (dot < 0.0)
    sign = -1.0;

  for (j = 0; j < adj_facet_list1.size(); j++)
  {
    facet = adj_facet_list1.get_and_step();
    normal = sign * facet->normal();
    normal.normalize();
    avg_normal += (facet->weight() / totangle) * normal;
  }

  // set the new normal on both points

  avg_normal.normalize();
  pt0->normal(avg_normal);
  CubitVector temp_vector = sign * avg_normal;
  pt1->normal(temp_vector);
  double coefd = -(pt0->coordinates()%avg_normal);
  pt0->d_coef( coefd );
  pt1->d_coef( sign * coefd );

  return CUBIT_SUCCESS;
}

//===============================================================================
//Function:  print_me (PRIVATE)
//Description: debug
//===============================================================================
void ChollaEngine::print_me()
{
  FILE *fp = fopen( "fbg.out", "w" );
  int ii, jj, kk;
  for (ii=0; ii<chollaSurfaceList.size(); ii++)
  {
    ChollaSurface *chsurf_ptr = chollaSurfaceList.get_and_step();
    fprintf( fp, "*********** Surface %d ***************\n", chsurf_ptr->get_id());
    DLIList<ChollaCurve *> fcm_list;
    chsurf_ptr->get_curves( fcm_list );
    for (jj=0; jj<fcm_list.size(); jj++)
    {
      ChollaCurve *chcurv_ptr = fcm_list.get_and_step();
      fprintf(fp, "   Curve %d\n", chcurv_ptr->get_id() );
      //CubitVector start, end;
      DLIList<ChollaPoint *> fpm_list = chcurv_ptr->get_points();
      for (kk=0; kk<fpm_list.size(); kk++)
      {
        ChollaPoint *chpnt_ptr = fpm_list.get_and_step();
        FacetEntity *facet_point = chpnt_ptr->get_facets();
        if (facet_point)
        {
          CubitPoint *point_ptr = CAST_TO(facet_point, CubitPoint);
          CubitVector pt = point_ptr->coordinates();
          fprintf(fp, "      (%d) x = %8.4f y = %8.4f z = %8.4f\n",
                  chpnt_ptr->get_id(), pt.x(), pt.y(), pt.z());
        }
      }
    }
  }
  fclose(fp);
}

//===============================================================================
//Function:  dump (PUBLIC)
//Description: debug
//===============================================================================
void ChollaEngine::dump( const char *filename, double angle )
{
  int include_results = 0;
  int num_face = faceList.size();
  int num_edge = edgeList.size();
  int num_vert = pointList.size();
  int *face_edge = new int [3*num_face];
  int *edge_vert = new int [2*num_edge];
  double *vert = new double [3*num_vert];

  int ii;
  CubitVector loc;
  CubitPoint *point_ptr;
  pointList.reset();
  for (ii=0; ii<num_vert; ii++)
  {
    point_ptr = dynamic_cast<CubitPoint *>(pointList.get_and_step());
    point_ptr->set_id( ii );
    loc = point_ptr->coordinates();
    vert[ii*3] = loc.x();
    vert[ii*3+1] = loc.y();
    vert[ii*3+2] = loc.z();
  }

  CubitFacetEdge *edge_ptr;
  edgeList.reset();
  for (ii=0; ii<num_edge; ii++)
  {
    edge_ptr = dynamic_cast<CubitFacetEdge *>(edgeList.get_and_step());
    edge_ptr->set_id( ii );
    edge_vert[ii*2] = edge_ptr->point( 0 )->id();
    edge_vert[ii*2+1] = edge_ptr->point( 1 )->id();
  }

  CubitFacet *facet_ptr;
  faceList.reset();
  for(ii=0; ii<num_face; ii++)
  {
    facet_ptr = dynamic_cast<CubitFacet *>(faceList.get_and_step());
    face_edge[ii*3] = facet_ptr->edge( 0 )->id();
    face_edge[ii*3+1] = facet_ptr->edge( 1 )->id();
    face_edge[ii*3+2] = facet_ptr->edge( 2 )->id();
  }

  double *edge_control_points = new double [3*3*num_edge];
  double *tri_control_points = new double [3*6*num_face];
  double *quad_control_points = NULL;

  char mesh_filename[256];
  sprintf(mesh_filename, "%s.mesh", filename);
  dumpMesh(mesh_filename, include_results, angle, num_face, 0, num_edge,
           num_vert, face_edge, NULL, edge_vert, vert, edge_control_points, 
           tri_control_points, quad_control_points);

  constructBezier( angle, num_face, 0, num_edge, num_vert, face_edge, NULL,
                   edge_vert, vert, edge_control_points, tri_control_points,
                   quad_control_points );

  char result_filename[256];
  sprintf(result_filename, "%s.results", filename);
  dumpResults(result_filename, num_edge, num_face, 0, 
              edge_control_points, tri_control_points, quad_control_points );

  include_results = 1;
  char full_filename[256];
  sprintf(full_filename, "%s.full", filename);
  dumpMesh(full_filename, include_results, angle, num_face, 0, num_edge,
           num_vert, face_edge, NULL, edge_vert, vert, edge_control_points, 
           tri_control_points, quad_control_points);

  delete [] face_edge;
  delete [] edge_vert;
  delete [] vert;
  delete [] edge_control_points;
  delete [] tri_control_points;
}

//===============================================================================
//Function:  mark_features (PUBLIC)
//Description: Calling this before ChollaEngine create geometry / topology
// allows you to define features points. This can be used in conjunction
// with the feature angle option or exclusively (original intention).
//===============================================================================
void ChollaEngine::mark_features (DLIList<CubitPoint*> &feature_points)
{
  int ii;
  TDGeomFacet *td_gm = NULL;
  CubitPoint *curr_pt = NULL;

    // loop through supplied CubitPoint*
    // set the hit flag on each one
  
  for (ii=0; ii < feature_points.size(); ++ii)
  {
    curr_pt = feature_points.get_and_step();
    TDGeomFacet::add_geom_facet(curr_pt, -1);
    td_gm = TDGeomFacet::get_geom_facet(curr_pt);
    td_gm->set_hit_flag(1);
    curr_pt->set_as_feature();
  }
}

//===============================================================================
//Function:  mark_features (PUBLIC)
//Description: Calling this before ChollaEngine create geometry / topology
// allows you to define features points. This can be used in conjunction
// with the feature angle option or exclusively (original intention).
//===============================================================================
void ChollaEngine::mark_features (DLIList<CubitFacetEdge*> &feature_edges)
{
  int ii;
  CubitFacetEdge *curr_edge = NULL;
	
	// loop through supplied CubitFaceEdge*
	// set the hit flag on each one

  for (ii=0; ii < feature_edges.size(); ++ii)
  {
    curr_edge = feature_edges.get_and_step();
    curr_edge->set_as_feature();
  }
}

  
static CubitStatus create_tri_facets(int *face_list,int current_position,
                                     CubitPoint **point_list,
                                     DLIList <CubitFacet*> &facet_list)
{
  int step = face_list[current_position];
    //This function is hard coded to work for splitting 4, 5, and 6
    //faceted points into triangles.
  if ( step > 6 || step < 4 )
  {
    PRINT_ERROR("Trying to split facets with wrong function.\n");
    return CUBIT_FAILURE;
  }
    //Basically get the points in a list first.
  int ii = current_position + 1, jj;
  DLIList <CubitPoint*> temp_points;
  
  for( jj = ii; jj < ii+step; jj++ )
  {
    int index_to_point = face_list[jj];
    assert(index_to_point >= 0 );
    temp_points.append(point_list[index_to_point]);
  }
    //Now cycle through the points and creat facets from them.
  CubitPoint *point_1, *point_2, *point_3;
    //This definatly could be improved to get more optimal triangles...
    //Create the first two facets in the normal way.  This is the
    //same for 4, 5 and 6 points.
  for ( jj = 0; jj < 2; jj++ )
  {
    point_1 = temp_points.get_and_step();
    point_2 = temp_points.get_and_step();
    point_3 = temp_points.get();
      //If these are colinear we are screwed...
    CubitFacet *new_facet = (CubitFacet *)
      new CubitFacetData(point_1, point_2, point_3);
    facet_list.append(new_facet);
  }
    //create the third facet for 5 points.
  if ( step == 5 )
  {
      //Now create the last facet with the first point, last point and
      //third point.
    temp_points.reset();
    point_1 = temp_points.prev();
    point_2 = temp_points.get();
    point_3 = temp_points.next(2);
    CubitFacet *new_facet = (CubitFacet *)
      new CubitFacetData(point_1, point_2, point_3);
    facet_list.append(new_facet);
  }
    //create the third and fourth for 6 points.
  else if ( step == 6 )
  {
      //get the current point (should be position 5/6)
    point_1 = temp_points.get();
      //get the next point. (should be position 6/6)
    point_2 = temp_points.next();
      //get the prev(2) point. ( should be position 3/6)
    point_3 = temp_points.prev(2);
    CubitFacet *new_facet = (CubitFacet *)
      new CubitFacetData(point_1, point_2, point_3);
    facet_list.append(new_facet);
      // get the next point. (should be 6/6)
    point_1 = temp_points.next();
      // get the next 2 point ( should be 1/6)
    point_2 = temp_points.next(2);
      //point3 stays the same (position 3/6).
    new_facet = (CubitFacet *)
      new CubitFacetData(point_1, point_2, point_3);
    facet_list.append(new_facet);
  }
  return CUBIT_SUCCESS;
}

//===================================================================================
// Description: given facets in the form of a GMem data structure create new
//              CubitFacet and CubitPoint classes and return in lists
// Notes: 
// Author: 
// Date: 
//===================================================================================
CubitStatus ChollaEngine::get_facets(GMem& gMem, DLIList<CubitFacet*> &facet_list, 
                                     DLIList<CubitPoint*> &dl_point_list,
                                     bool insert_null_facets )
{
  if(gMem.fListCount == 0)
    return CUBIT_FAILURE;

  GPoint* plist = gMem.point_list();
  CubitPoint **point_list = (CubitPoint **)
    new CubitPointData *[gMem.pointListCount];
  double x, y, z;
  int i,j;
  for (i = 0; i<gMem.pointListCount; i++) 
  {
    x = plist[i].x;
    y = plist[i].y;
    z = plist[i].z;
    CubitPoint *new_point = (CubitPoint *) new CubitPointData(x, y, z);
    point_list[i] = new_point;
    dl_point_list.append(new_point);
  } 
  int *face_list = gMem.facet_list();
  int step;
  DLIList<CubitPoint*> temp_points;
  CubitFacet *new_facet;
  int ii, index_to_point;
  CubitStatus stat1;
  for ( i = 0; i < gMem.fListCount; i++ )
  {
    step = face_list[i];
      //Now get the points for this facet.
    switch (step)
    {
      case 3:
        ii = i + 1;
        temp_points.clean_out();
        for( j = ii; j < ii+step; j++ )
        {
          index_to_point = face_list[j];
          assert(index_to_point >= 0 && index_to_point <= gMem.pointListCount );
          temp_points.append(point_list[index_to_point]);
        }
        new_facet = (CubitFacet *)
          new CubitFacetData(temp_points.get(), 
                             temp_points.next(),
                             temp_points.next(2) );
        facet_list.append(new_facet);
        i = i+step;
        break;
      case 4:
      case 5:
      case 6:
        stat1 = create_tri_facets(face_list,i, point_list,
                                  facet_list);
        if ( stat1 != CUBIT_SUCCESS )
        {
          PRINT_ERROR("Facets were not triangular and couldn't be split\n");
          return stat1;
        }
        i = i+step;
        break;
      default:
          //We currently can't handle this...
          //Eventually we will have to split this polygon into more facets.
        PRINT_ERROR("Facets were not triangular, and couldn't be made so.\n");
        PRINT_ERROR("Surface has: %d points forming a facet.\n", step);
        if( insert_null_facets )
          facet_list.append( NULL );
        return CUBIT_FAILURE;
    }
  }
  delete [] point_list;
  return CUBIT_SUCCESS;
}

//=============================================================================
//Function: collapses a given curve to a given point_to_keep
//          
//Description:
//Note: the hashing and other data structures are not updated during collapse.
//Author: william roshan quadros
//Date:  12/15/08
//=============================================================================
CubitStatus ChollaEngine::collapse_curve( ChollaCurve *cholla_curve, ChollaPoint *point_to_keep )
{
  // Make sure there are no underlying facet edges connected with this curve
  if( cholla_curve->get_facet_list().size() )
    return CUBIT_FAILURE;
  
  // Disassociate from cholla_surfaces 
  DLIList<ChollaSurface *> &ref_surfaces = cholla_curve->get_surfaces();
  int i;
  ChollaSurface *cholla_surface;
  for( i = 0; i < ref_surfaces.size(); i++ )
  {
    cholla_surface = ref_surfaces.get_and_step();
    cholla_surface->remove_curve( cholla_curve );
  }
  ref_surfaces.clean_out();

  // Disassociate from cholla point_to_delete 
  // Make sure curve has at least one point
  DLIList<ChollaPoint *> &cholla_pnts = cholla_curve->get_points();
  if( 0 == cholla_pnts.size()  )
    return CUBIT_FAILURE;

    // free eval tool
  if( cholla_curve->get_eval_tool() )
  {
    delete cholla_curve->get_eval_tool(); 
  }

  // find the point to delete
  ChollaPoint *first_pnt = cholla_pnts.get_and_step();
  ChollaPoint *last_pnt = cholla_pnts.get();
  ChollaPoint *point_to_del = NULL;
  switch( cholla_pnts.size() )
  {

  case 2:    
    if( point_to_keep == last_pnt )
    {
      point_to_del = first_pnt;
    }
    else
    {
      if( point_to_keep == first_pnt )
      {
        point_to_del = last_pnt;
      }
    }

    // disassociate end points from cholla_curve
    if( point_to_del )
    {
      point_to_del->remove_curve( cholla_curve );
      cholla_pnts.remove( point_to_del );
    }
    if( point_to_keep )
    {
      point_to_keep->remove_curve( cholla_curve );
      cholla_pnts.remove( point_to_keep );
    }

    // Merge two points 
    merge_two_points( point_to_keep, point_to_del );

    break;

  case 1:
    // disassociate end points from cholla_curve
    point_to_del = first_pnt;
    point_to_del->remove_curve( cholla_curve );
    cholla_pnts.remove( point_to_del );
    
    if( 0 == point_to_del->get_curve_list_ptr()->size() )
    {
      // remove point_to_del from ChollaEngine
      remove_point( point_to_del );

      // free point_to_del
      delete point_to_del;
    }    
    break;

  default:
    return CUBIT_FAILURE;
  }


  // Disassociate curve from chollaEngine
  this->remove_curve( cholla_curve );

  // free the cholla_curve
  delete cholla_curve;

  return CUBIT_SUCCESS;
}

  // Merge two points 
CubitStatus ChollaEngine::merge_two_points( ChollaPoint *point_to_keep, ChollaPoint *point_to_del )
{
  // move cholla curves from point_to_del to point_to_keep
  DLIList<ChollaCurve *> &del_pnt_curves = point_to_del->get_curves();
  int i;
  ChollaCurve *cholla_curve;
  for( i = 0; i < del_pnt_curves.size(); i++ )
  {
    cholla_curve = del_pnt_curves.get_and_step();
    
    cholla_curve->remove_point( point_to_del );
    cholla_curve->add_point( point_to_keep );
    point_to_keep->add_curve( cholla_curve );
  }
  del_pnt_curves.clean_out();

  // if point_to_keep facet is null use point_to_del's facet
  if( NULL == point_to_keep->get_facets() && point_to_del->get_facets() )
  {
    point_to_keep->add_facet( point_to_del->get_facets() );
  }
  
  // what we should do with data member id?  should we move if point_to_keep->get_id() !=  0?

  // remove point_to_del from ChollaEngine
  remove_point( point_to_del );

  // free point_to_del
  delete point_to_del;

  return CUBIT_SUCCESS;
}



CubitStatus ChollaEngine::disassociate_surface( ChollaSurface *cholla_surf )
{
   // remove surface from volume
  DLIList<ChollaVolume *> cholla_volumes;
  cholla_surf->get_volumes(cholla_volumes);
  int i;
  ChollaVolume *cholla_volume;
  for (i=0; i<cholla_volumes.size(); i++)
  {
    cholla_volume = cholla_volumes.get_and_step();
    cholla_volume->remove_surface( cholla_surf );
  }
  
   // Disassociate curves from suface
  DLIList< ChollaCurve *> cholla_curves;
  cholla_surf->get_curves( cholla_curves );
  
  ChollaCurve *cholla_curve;
  for( i = 0; i < cholla_curves.size(); i++ )
  {
    cholla_curve = cholla_curves.get_and_step();
    
    cholla_curve->remove_surface( cholla_surf );
    cholla_surf->remove_curve( cholla_curve );
  }
  cholla_curves.clean_out();

  // Remove facet tool if available ?
  if( cholla_surf->get_eval_tool() )
  {
#ifdef _DEBUG
    //PRINT_INFO("WARNING: delete ch_surf->facet_eval_tool safely \n");//   
#endif
    //delete cholla_surf->get_eval_tool();
    // set eval tool to NULL
  }
  
  // remove from ChollaEngine
  remove_surface( cholla_surf );

  return CUBIT_SUCCESS;

}
//=============================================================================
//Function: collapses a given surface 
//Description: currently it disassociates surface with its curves and ChollaEngine 
//Note: make sure there are no underlying facets and model is automatically watertight when input cholla_surf is deleted
//Author: william roshan quadros
//Date:  12/15/08
//=============================================================================

CubitStatus ChollaEngine::collapse_surface( ChollaSurface *cholla_surf )
{
  // Make sure there are no facets 
  if( cholla_surf->get_facet_list().size() )
     return CUBIT_FAILURE; 

  disassociate_surface( cholla_surf );
  
  // free cholla surface
  delete cholla_surf;

  return CUBIT_SUCCESS;
}

CubitStatus ChollaEngine::remove_facet_entity( CubitFacet *facet, ChollaSurface *cholla_surf )
{
  this->faceList.remove( facet );

  if( cholla_surf )
    cholla_surf->remove_facet( facet );

  return CUBIT_SUCCESS;
}
  
CubitStatus ChollaEngine::remove_facet_entity( CubitFacetEdge *facet_edge, ChollaCurve *cholla_curve )
{
  this->edgeList.remove( facet_edge );

  if( cholla_curve )
    cholla_curve->remove_facet( facet_edge );

  return CUBIT_SUCCESS;
}

CubitStatus ChollaEngine::remove_facet_entity( CubitPoint *facet_pnt, ChollaPoint *cholla_point )
{
  this->pointList.remove( facet_pnt );

  if( cholla_point )
    cholla_point->remove_facet(/* facet_pnt */);

  return CUBIT_SUCCESS;
}

CubitStatus ChollaEngine::remove_facet_entity( CubitFacet *facet, std::set<ChollaEntity *> &cholla_surfs )
{
  std::set<ChollaEntity *>::iterator set_it;
  for( set_it = cholla_surfs.begin(); set_it != cholla_surfs.end(); set_it++ )
  {
    if( dynamic_cast< ChollaSurface *>(*set_it) )
    {
      static_cast<ChollaSurface *>(*set_it)->remove_facet( facet );
    }
  }

  this->faceList.remove( facet );

  return CUBIT_SUCCESS;
}
  
CubitStatus ChollaEngine::remove_facet_entity( CubitFacetEdge *facet_edge, std::set<ChollaEntity *> &cholla_curves )
{
  
  std::set<ChollaEntity *>::iterator set_it;
  ChollaEntity *ptr_entity;
  ChollaCurve *ptr_curve;
  for( set_it = cholla_curves.begin(); set_it != cholla_curves.end(); set_it++ )
  {
    ptr_entity = *set_it;
    ptr_curve = dynamic_cast<ChollaCurve *>( ptr_entity );
    if( dynamic_cast<ChollaCurve*>(*set_it) )
    { 
      dynamic_cast<ChollaCurve*>(*set_it)->remove_facet( facet_edge );
    }
  }

  this->edgeList.remove( facet_edge );

  return CUBIT_SUCCESS;
}

CubitStatus ChollaEngine::replace_facet_entity( CubitFacetEdge *remove_edge, CubitFacetEdge *replace_edge, std::set<ChollaEntity *> cholla_curves )
{
  std::set<ChollaEntity *>::iterator set_it; 
  for( set_it = cholla_curves.begin(); set_it != cholla_curves.end(); set_it++ )
  {
    if( dynamic_cast<ChollaCurve*>(*set_it) )
    { 
      static_cast<ChollaCurve*>(*set_it)->replace_facet( remove_edge, replace_edge );
    }
  }

  return CUBIT_SUCCESS;
}


CubitStatus ChollaEngine::remove_facet_entity( CubitPoint *facet_pnt, std::set<ChollaEntity *> &cholla_points )
{
  this->pointList.remove( facet_pnt );

  std::set<ChollaEntity *>::iterator set_it;
  for( set_it = cholla_points.begin(); set_it != cholla_points.end(); set_it++ )
  {
    if( dynamic_cast<ChollaPoint*>(*set_it) )
    { 
      static_cast<ChollaPoint*>(*set_it)->remove_facet( facet_pnt );
    }
  }

  return CUBIT_SUCCESS;
}

/*
// disassoicate ch_entity from ch_points, ch_curve, ch_surf, and ch_engine
CubitStatus ChollaEngine::disassociate_curve( ChollaCurve *ch_curve )
{
  
  // disassociate from cholla points
  ch_curve->disassociate_from_points();
  
  // disassociate from cholla surface
  ch_curve->disassociate_from_surfaces();

  // disassociate from cholla engine
  remove_curve( ch_curve );

  return CUBIT_SUCCESS;
}
*/

CubitStatus ChollaEngine::create_surface( int block_id,
                                          ChollaSurface *&new_ch_surf )
{
  new_ch_surf = new ChollaSurface( block_id );
  
  // add it to cholla engine
  chollaSurfaceList.append( new_ch_surf ); 

  return CUBIT_SUCCESS;
}

 // Create new ch_curve and associate with end ch_points.  Add ch_curve to chollaEngine
CubitStatus ChollaEngine::create_curve(int block_id, 
                                       ChollaPoint *new_ch_pnt0, 
                                       ChollaPoint *new_ch_pnt1, 
                                       ChollaCurve *&new_ch_curve
                                    )
{

  
  new_ch_curve = new ChollaCurve( block_id );

  // associate with cholla points
  new_ch_curve->add_point( new_ch_pnt0 );
  new_ch_pnt0->add_curve( new_ch_curve );
  new_ch_curve->add_point( new_ch_pnt1 );
  new_ch_pnt1->add_curve( new_ch_curve );


  // update start and end facet points
  new_ch_curve->set_start( CAST_TO( new_ch_pnt0->get_facets(), CubitPointData) );
  new_ch_curve->set_end( CAST_TO( new_ch_pnt1->get_facets(), CubitPointData) );

  // add it to cholla engine
  chollaCurveList.append( new_ch_curve ); 

  return CUBIT_SUCCESS;
}

CubitStatus ChollaEngine::create_point( CubitPoint *pnt, ChollaPoint * &new_ch_pnt )
{
  new_ch_pnt = new ChollaPoint();
  new_ch_pnt->add_facet( pnt );
  double *coordinate = new double[3];
  CubitVector coord_vect = pnt->coordinates();
  coordinate[0] = coord_vect[0];
  coordinate[1] = coord_vect[1];
  coordinate[2] = coord_vect[2];
  new_ch_pnt->assign_geometric_point( (void*) coordinate );

  return CUBIT_SUCCESS;
}


// disassoicate ch_curve from ch_points, ch_surfs, and ch_engine
CubitStatus ChollaEngine::disassociate_curve( ChollaCurve *ch_curve, bool disassociate_with_vert, bool disassociate_with_surf, bool disassociate_with_model )
{
  // remove ch_curve from ch_point
  int i;
  if( disassociate_with_vert )
  {
    for( i = 0; i < ch_curve->get_points().size(); i++ )
    {
      ch_curve->get_points().get_and_step()->remove_curve( ch_curve );
    }
    // remove ch_points
    ch_curve->get_points().clean_out();
  }

  if( disassociate_with_surf )
  {
    // remove ch_curve in ch_surface
    ChollaSurface *ptr_ch_surf;   
    for( i = 0; i < ch_curve->get_surfaces().size(); i++ )
    {
      ptr_ch_surf = ch_curve->get_surfaces().get_and_step();
      ptr_ch_surf->remove_curve( ch_curve );
    }
    // remove ch_surface
    ch_curve->get_surfaces().clean_out();
  }

  if( disassociate_with_model )
  {
    // remove from ch_engie
    this->remove_curve( ch_curve );
  }

  return CUBIT_SUCCESS;
}

//=============================================================================
//Function: detach_volumes 
//Description: Create independent manifold volumes from the non-manifold set 
//Note: 
//Author:sjowen
//Date:  09/10/09
//=============================================================================
CubitStatus ChollaEngine::detach_volumes()
{
  CubitStatus rv = CUBIT_SUCCESS;
 
  // define maps between the original entities and their copies so we can 
  // use them for detaching the facets
  
  std::map<ChollaSurface *, ChollaSurface *> surf_map;
  std::multimap<ChollaCurve *, ChollaCurve *> curve_map;
  std::multimap<ChollaPoint *, ChollaPoint *> point_map; 
  
  for( int i=chollaVolumeList.size(); i--; )
  {
    ChollaVolume *chvol = chollaVolumeList.get_and_step();
    DLIList<ChollaSurface*> ch_surfaces;
    chvol->get_surfaces( ch_surfaces );    

    std::map<ChollaCurve *, ChollaCurve *> tmp_curve_map;
    std::map<ChollaPoint *, ChollaPoint *> tmp_point_map; 
    
    DLIList<ChollaSurface*> surfaces_to_detach;
    DLIList<ChollaCurve*> curves_to_detach;
    for( int j=ch_surfaces.size(); j--; )
    {
      ChollaSurface *chsurf_ptr = ch_surfaces.get_and_step();     

      // look for any surfaces that are associated with exactly two volumes
      if (chsurf_ptr->num_volumes() == 2)      
        surfaces_to_detach.append( chsurf_ptr );   
      else
      {
        //look for any curves that are attached to two volumes
        DLIList<ChollaCurve*> ch_curves;
        chsurf_ptr->get_curves( ch_curves );
        for( int k=ch_curves.size(); k--; )
        {
          ChollaCurve *tmp_curve = ch_curves.get_and_step();
          if( tmp_curve->num_volumes() > 1 )
            curves_to_detach.append( tmp_curve );
        }
      }
    }   

    //filter out curves that are in surfs in 'surfaces_to_detach'
    //since they will get detached when the parent surface gets detached
    for( int j=curves_to_detach.size(); j--; )
    {
      ChollaCurve *tmp_curve = curves_to_detach.get();
      DLIList<ChollaSurface*> tmp_surfs = tmp_curve->get_surfaces();

      for( int k=tmp_surfs.size(); k--; )
      {
        ChollaSurface *tmp_surf = tmp_surfs.get_and_step();
        if( surfaces_to_detach.is_in_list( tmp_surf ) )
        {
          curves_to_detach.change_to(NULL);
          break;
        }
      }
      curves_to_detach.step();    
    }  
    
    curves_to_detach.uniquify_unordered();
    curves_to_detach.remove_all_with_value( NULL );

    if( surfaces_to_detach.size() || curves_to_detach.size() )
    {
      rv = detach_surfaces( surfaces_to_detach, curves_to_detach, chvol, surf_map, tmp_curve_map, tmp_point_map );
      if (rv != CUBIT_SUCCESS)
        return CUBIT_FAILURE;
      
      //for all detached surfaces, detach all the facets and create new ones for
      //the detached surfaces
      rv = detach_facets(surfaces_to_detach, curves_to_detach, chvol, surf_map, tmp_curve_map, tmp_point_map);
      if (rv != CUBIT_SUCCESS)
        return rv;

      // set the cubitpoints at the correct start/end locations on the cholla curves   
      rv = set_curve_endpoints(tmp_point_map, tmp_curve_map );
      if (rv != CUBIT_SUCCESS)
        return rv;
    }

    //append the map to the multi-map
    std::map<ChollaCurve*, ChollaCurve*>::iterator iter = tmp_curve_map.begin();
    for( ; iter!=tmp_curve_map.end(); iter++)
      curve_map.insert( std::pair<ChollaCurve*,ChollaCurve*>(iter->first, iter->second) );

    //append the map to the multi-map
    std::map<ChollaPoint*, ChollaPoint*>::iterator pt_iter = tmp_point_map.begin();
    for( ; pt_iter!=tmp_point_map.end(); pt_iter++)
      point_map.insert( std::pair<ChollaPoint*,ChollaPoint*>(pt_iter->first, pt_iter->second) );   
  } 
  
  // create the evaluation tools on the new entities
  
  CubitBoolean use_feature_angle = CUBIT_FALSE; 
  double min_dot = 0.0;                    
  int interp_order = 0;                 
  CubitBoolean smooth_non_manifold = CUBIT_FALSE;
  CubitBoolean split_surfaces = CUBIT_FALSE;
  DLIList<ChollaSurface *> cholla_surfaces;
  DLIList<ChollaCurve *> cholla_curves;
  
  // get lists of new surfaces and curves from the maps
  
  std::map<ChollaSurface *, ChollaSurface *>::iterator smap_it;
  for (smap_it=surf_map.begin(); smap_it != surf_map.end(); smap_it++)
  {
    cholla_surfaces.append(smap_it->second);
  }
  
  std::multimap<ChollaCurve *, ChollaCurve *>::iterator cmap_it;
  for (cmap_it=curve_map.begin(); cmap_it != curve_map.end(); cmap_it++)
  {
    cholla_curves.append(cmap_it->second);
  }
  
  // rebuild the new curves so the edges are oriented correctly
  
  for (int ii=0; ii<cholla_curves.size(); ii++)
  {
    ChollaCurve *cholla_curve = cholla_curves.get_and_step();
    DLIList<ChollaPoint *>  cholla_points = cholla_curve->get_points();    
    int periodic = 0;
    if (cholla_points.size() == 1)
      periodic = 1;
    else
      assert(cholla_points.size() == 2);  // assuming we have exactly two end points so far
    CubitPoint *start_point, *end_point;
    cholla_curve->get_ends(start_point, end_point);        
    int max_edges = cholla_curve->num_edges();
    cholla_curve->clean_out_edges();  // edges will be added in build_curve_from_edges 
    rv = cholla_curve->build_curve_from_edges( start_point, periodic, max_edges, NULL, cholla_curve );
    if (rv != CUBIT_SUCCESS)
      return rv;
  } 
  
  // replace the facets on the existing eval tools on the curves and surfaces bounding the interface.
  // for simplicity, just replace them on all the existing eval tools on 
  
  for (int icurv = 0; icurv < chollaCurveList.size(); icurv++)
  {
    ChollaCurve *chcurv_ptr = chollaCurveList.get_and_step();
    CurveFacetEvalTool *ceval = chcurv_ptr->get_eval_tool();
    if (ceval != NULL)
    {
      DLIList<FacetEntity *> fedges = chcurv_ptr->get_facet_list();
      DLIList<CubitFacetEdge *> edges;
      CAST_LIST( fedges, edges, CubitFacetEdge );
      ceval->replace_facets(edges);
    }
  }
  
  for (int isurf=0; isurf<chollaSurfaceList.size(); isurf++)
  {
    ChollaSurface *chsurf_ptr = chollaSurfaceList.get_and_step();
    FacetEvalTool *seval = chsurf_ptr->get_eval_tool();
    if (seval != NULL)
    {
      DLIList<FacetEntity *> ffacets = chsurf_ptr->get_facet_list();
      DLIList<CubitFacet *> facets;
      CAST_LIST( ffacets, facets, CubitFacet );
      seval->replace_facets(facets);
    }
  }
 
  
  // build the eval tools
      
  rv = build_eval_tools(cholla_surfaces,
                        cholla_curves,
                        interp_order, use_feature_angle, 
                        min_dot, smooth_non_manifold,
                        split_surfaces );
    
  return rv;
}

//=============================================================================
//Function: detach_surface 
//Description: given a non-manifold surface in a cholla model, create a copy
//             and update child entities
//Notes:  assumes that chsurf_ptr has exactly 2 adjacent volumes
//Author:sjowen
//Date:  09/18/09
//=============================================================================
CubitStatus ChollaEngine::detach_surfaces(DLIList<ChollaSurface*> &chsurfs,
                                          DLIList<ChollaCurve*> &chcurves,
                                         ChollaVolume *vol_getting_new_surfs,
                                         std::map<ChollaSurface *, ChollaSurface *> &surf_map,
                                         std::map<ChollaCurve *, ChollaCurve *> &curve_map,
                                         std::map<ChollaPoint *, ChollaPoint *> &point_map)
{
  CubitStatus rv = CUBIT_SUCCESS;
  
  // detach the surface from its volumes

  DLIList<ChollaCurve*> all_curves_in_surfs;
  DLIList<ChollaPoint*> all_points_in_surfs;

  std::multimap<ChollaCurve*, ChollaSurface*> orig_curve_to_new_surf_map;
  
  for( int i=chsurfs.size(); i--; )
  {
    ChollaSurface *chsurf_ptr = chsurfs.get_and_step();
    DLIList<ChollaVolume *> chvol_list;
    chsurf_ptr->get_volumes(chvol_list);
    ChollaVolume *chvol1_ptr = chvol_list.get_and_step();
    ChollaVolume *chvol2_ptr = chvol_list.get();

    if( chvol2_ptr != vol_getting_new_surfs )      
    {
      chvol1_ptr = chvol2_ptr;
      chvol2_ptr = vol_getting_new_surfs;      
    }

    assert(chvol1_ptr != NULL);  
    assert(chvol2_ptr != NULL);
    assert(chvol1_ptr != chvol2_ptr);
    
    // create a copy of the non-manifold surface and attach it to volume 2

    ChollaSurface *newchsurf_ptr = new ChollaSurface( chsurf_ptr->get_block_id() );
    chollaSurfaceList.append(newchsurf_ptr);
    surf_map.insert(std::pair<ChollaSurface *, ChollaSurface *>(chsurf_ptr, newchsurf_ptr));

    chvol2_ptr->remove_surface(chsurf_ptr);
    chsurf_ptr->remove_volume(chvol2_ptr);
    chvol2_ptr->add_surface(newchsurf_ptr);
    newchsurf_ptr->add_volume(chvol2_ptr);
    
    chsurf_ptr->set_merge_partner(newchsurf_ptr);
    newchsurf_ptr->set_merge_partner(chsurf_ptr);

    // detach the curves  
    DLIList<ChollaCurve *> chcurv_list;
    chsurf_ptr->get_curves(chcurv_list);    

    for( int k=chcurv_list.size(); k--; )
    {
      ChollaCurve *chcurve = chcurv_list.get_and_step();
      orig_curve_to_new_surf_map.insert( std::pair<ChollaCurve*, ChollaSurface*>( chcurve, newchsurf_ptr ) );
    }

    all_curves_in_surfs += chcurv_list;

    DLIList<ChollaPoint*> chpt_list;
    chsurf_ptr->get_vertices( chpt_list );
    all_points_in_surfs += chpt_list;
  }

  all_curves_in_surfs.uniquify_unordered();  
  
  for(int icurv = 0; icurv < all_curves_in_surfs.size(); icurv++)
  {    
    ChollaCurve *chcurv_ptr = all_curves_in_surfs.get_and_step();

    std::multimap<ChollaCurve*, ChollaSurface*>::iterator iter, upper_iter;
    DLIList<ChollaSurface*> new_surfs;

    iter = orig_curve_to_new_surf_map.find( chcurv_ptr );

    assert( iter != orig_curve_to_new_surf_map.end() );

    upper_iter = orig_curve_to_new_surf_map.upper_bound( iter->first );

    for(; iter!=upper_iter; ++iter )
      new_surfs.append( iter->second );

    rv = detach_curve( chcurv_ptr, new_surfs, vol_getting_new_surfs, curve_map, point_map );
    if (rv != CUBIT_SUCCESS)
      return rv;  
  }

  for( int i=chcurves.size(); i--; )
  {
    DLIList<ChollaPoint*> tmp_pts = chcurves.get_and_step()->get_points();
    all_points_in_surfs += tmp_pts;
  }

  rv = detach_curves( chcurves, vol_getting_new_surfs, curve_map, point_map );
  if( rv != CUBIT_SUCCESS )
    return rv;

  //add the points in the cholla     
  all_points_in_surfs.uniquify_unordered();
  for (int ipt = 0; ipt < all_points_in_surfs.size(); ipt++)
  {
    ChollaPoint *chpt_ptr = all_points_in_surfs.get_and_step();
    
    rv = detach_point( chpt_ptr, vol_getting_new_surfs, point_map, curve_map );
    if (rv != CUBIT_SUCCESS)
      return rv;
  }   
  
  return rv;
}

CubitStatus ChollaEngine::detach_curves( DLIList<ChollaCurve*> &curves,
                                         ChollaVolume *detaching_volume,                                         
                                         std::map<ChollaCurve *, ChollaCurve *> &curve_map,
                                         std::map<ChollaPoint *, ChollaPoint *> &point_map )
{
  CubitStatus rv = CUBIT_SUCCESS;

  DLIList<ChollaSurface*> tmp_surfs;
  DLIList<ChollaPoint*> points_to_detach;
  for( int i=curves.size(); i--; )
  {    
    ChollaCurve *tmp_curve = curves.get_and_step();    

    //get the surfaces in 'detaching_volume'
    tmp_surfs.clean_out();
    tmp_surfs = tmp_curve->get_surfaces();
    for( int k=0; k<tmp_surfs.size(); k++ )
    {
      if( !tmp_surfs[k]->is_in_volume( detaching_volume) )
        tmp_surfs[k] = NULL;      
    }
    tmp_surfs.remove_all_with_value(NULL);

    rv = detach_curve( tmp_curve, tmp_surfs, detaching_volume, curve_map, point_map );
    if (rv != CUBIT_SUCCESS)
      return rv;

    DLIList<ChollaPoint*> tmp_pts = tmp_curve->get_points();
    points_to_detach += tmp_pts;
  }

  //go through the map, setting up merge partners
  std::map<ChollaCurve*, ChollaCurve*>::iterator iter = curve_map.begin();
  for( ; iter != curve_map.end(); iter++ )
  {
    ChollaCurve *ch_curv = iter->first;
    ChollaCurve *new_ch_curv = iter->second;

    ch_curv->set_merge_partner( new_ch_curv );
    new_ch_curv->set_merge_partner( ch_curv );
  }
  return rv;
}                                          

//=============================================================================
//Function: detach_curve 
//Description: given a non-manifold curve in a cholla model, create a copy
//             and update child entities. updates the curve_map
//Author:sjowen
//Date:  09/18/09
//=============================================================================
CubitStatus ChollaEngine::detach_curve(ChollaCurve *chcurv_ptr,
                                       DLIList<ChollaSurface*> &new_surfs,
                                       ChollaVolume *chvol2_ptr,
                                       std::map<ChollaCurve*, ChollaCurve*> &curve_map,
                                       std::map<ChollaPoint*, ChollaPoint*> &point_map )
{   
  // create a copy of the curve on the surface and add it to volume 2  

  ChollaCurve *newchcurv_ptr = new ChollaCurve( chcurv_ptr->get_block_id() );
  chollaCurveList.append( newchcurv_ptr );
  curve_map.insert(std::pair<ChollaCurve *, ChollaCurve *>(chcurv_ptr, newchcurv_ptr));
  
  for( int i=new_surfs.size(); i--; )
  {
    ChollaSurface *newchsurf_ptr = new_surfs.get_and_step();
    newchsurf_ptr->add_curve(newchcurv_ptr);
    newchcurv_ptr->add_surface(newchsurf_ptr);    
  }

  // any surfaces attached to this chcurv that are in vol2 are removed and the newchcurv is added
  DLIList<ChollaSurface *> chcsurf_list = chcurv_ptr->get_surfaces();
  for (int icsurf = 0; icsurf < chcsurf_list.size(); icsurf++)
  {
    ChollaSurface *chcsurf_ptr = chcsurf_list.get_and_step();

    //if the adjacent surfaces are in both volumes, don't do anything
    if(chcsurf_ptr->is_in_volume(chvol2_ptr) )      
    {
      chcurv_ptr->remove_surface(chcsurf_ptr);
      chcsurf_ptr->remove_curve(chcurv_ptr);
      newchcurv_ptr->add_surface(chcsurf_ptr);
      chcsurf_ptr->add_curve_unique(newchcurv_ptr);
    }  
  } 

  return CUBIT_SUCCESS;
}

//=============================================================================
//Function: detach_point 
//Description: given a non-manifold point in a cholla model, create a copy
//             and update connectivity. updates the point_map
//Author:sjowen
//Date:  09/18/09
//=============================================================================
CubitStatus ChollaEngine::detach_point(ChollaPoint *chpt_ptr,
                                       ChollaVolume *chvol2_ptr,
                                       std::map<ChollaPoint *, ChollaPoint *> &point_map,
                                       std::map<ChollaCurve *, ChollaCurve *> &curve_map )                                       
{ 
  
  
  ChollaPoint *newchpt_ptr = new ChollaPoint();
  chollaPointList.append( newchpt_ptr );
  point_map.insert(std::pair<ChollaPoint*, ChollaPoint*>(chpt_ptr, newchpt_ptr));  
     
  //find the curves that it is in
  DLIList<ChollaCurve *> chptcurv_list = chpt_ptr->get_curves();
  for (int iptcurv = 0; iptcurv < chptcurv_list.size(); iptcurv++)
  {
    ChollaCurve *chptcurv_ptr = chptcurv_list.get_and_step();

    // for curves that were copied (on the interface), add the new chollapoint to the new curve

    std::map<ChollaCurve *, ChollaCurve *>::iterator cmap_it;
    cmap_it = curve_map.find(chptcurv_ptr);
    if (cmap_it != curve_map.end())
    {
      ChollaCurve *new_chcurv_ptr = cmap_it->second;      
      new_chcurv_ptr->add_point(newchpt_ptr);
      newchpt_ptr->add_curve(new_chcurv_ptr);
      
      chpt_ptr->set_merge_partner( newchpt_ptr );
      newchpt_ptr->set_merge_partner( chpt_ptr );
    }
    else
    {
      // remove curve in vol 2 (not on interface) from the original point and add it to the new point
      
      if (chptcurv_ptr->is_in_volume(chvol2_ptr))
      {
        chpt_ptr->remove_curve(chptcurv_ptr);
        chptcurv_ptr->remove_point(chpt_ptr);
        newchpt_ptr->add_curve(chptcurv_ptr);
        chptcurv_ptr->add_point(newchpt_ptr);

        chpt_ptr->set_merge_partner( newchpt_ptr );
        newchpt_ptr->set_merge_partner( chpt_ptr );
      }
    }
  }

  return CUBIT_SUCCESS;
}
  

//=============================================================================
//Function: detach_facets 
//Description: detach the individual facets to create a manifold representation
//Note: makes a copy of the facets on the merged surface and its child entities
//      and reconnects them locally
//Author:sjowen
//Date:  09/10/09
//=============================================================================
CubitStatus ChollaEngine::detach_facets(DLIList<ChollaSurface*> &chsurfs,
                                        DLIList<ChollaCurve*> &chcurves,
                                        ChollaVolume *chvol,
                                        std::map<ChollaSurface *, ChollaSurface *> &surf_map,
                                        std::map<ChollaCurve *, ChollaCurve *> &curve_map,
                                        std::map<ChollaPoint *, ChollaPoint *> &point_map)
{  
  
  CubitStatus rv = CUBIT_SUCCESS;
  
  std::vector<CubitPoint *> new_points;
  std::vector<CubitFacetEdge *> new_edges;
  
  rv = copy_facets_at_interface( chsurfs, chcurves, new_points, new_edges,
                                 surf_map, curve_map, point_map );
  if (rv != CUBIT_SUCCESS)
    return rv;
  
  rv = connect_facets_at_interface( chsurfs, chcurves, chvol, new_points, new_edges );
  if (rv != CUBIT_SUCCESS)
    return rv;
  
  return rv;
}
  
//=============================================================================
//Function: copy_facets_at_interface 
//Description: 
//Author:sjowen
//Date:  09/18/09
//=============================================================================
CubitStatus ChollaEngine::copy_facets_at_interface(DLIList<ChollaSurface*> &chsurfs,
                                                   DLIList<ChollaCurve*> &chcurves,
                                                   std::vector<CubitPoint *> &new_points,
                                                   std::vector<CubitFacetEdge *> &new_edges,                                                   
                                                   std::map<ChollaSurface *, ChollaSurface *> &surf_map,
                                                   std::map<ChollaCurve *, ChollaCurve *> &curve_map,
                                                   std::map<ChollaPoint *, ChollaPoint *> &point_map)

{
  CubitStatus rv = CUBIT_SUCCESS;
  
  // first set the marked flags on all facets entities on this surface to -1
  
  int ifacet;
  DLIList<FacetEntity *> facet_list;
  for( int i=chsurfs.size(); i--; )
  {
    ChollaSurface *chsurf_ptr = chsurfs.get_and_step();
    chsurf_ptr->get_facets( facet_list );
  }

  for( int i=chcurves.size(); i--; )
  {
    DLIList<FacetEntity*> tmp_facets = chcurves.get_and_step()->get_facet_list();
    facet_list += tmp_facets;
  }

  FacetDataUtil::mark_facets(facet_list, FACET_ENTITY_UNINITIALIZED);
  
  // create a copy of each of the facet entities on the surface.  The marked flag in the facet will
  // keep track of the new entity created.  It will be a location in the new_points or new_edges
  // array.  Once we are finished with creating points and edges, we can create the new facets
  // given the references in the marked flags.
  
  // create new points
  
  rv = copy_points_at_interface(facet_list, new_points, surf_map, curve_map, point_map);
  if (rv != CUBIT_SUCCESS)
    return rv;
  
  // create new edges
  
  rv = copy_edges_at_interface(facet_list, new_points, new_edges, surf_map, curve_map, point_map);
  if (rv != CUBIT_SUCCESS)
    return rv;
  
  // create new facets
  
  DLIList<CubitFacet *> new_facets;
  for (ifacet = 0; ifacet<facet_list.size(); ifacet++)
  {
    FacetEntity *facet_ptr = facet_list.get_and_step();
    CubitFacet *cfacet_ptr = dynamic_cast<CubitFacet *> (facet_ptr);
    if( cfacet_ptr )
    {
      CubitFacetEdge *fedge, *newfedges[3];
      CubitFacet *newcfacet_ptr;
      for (int ii=0; ii<3; ii++)
      {
        fedge = cfacet_ptr->edge(ii);
        int idx = fedge->marked();
        newfedges[ii] = new_edges[idx];
      }
      newcfacet_ptr = (CubitFacet *) new CubitFacetData(newfedges[0], newfedges[1], newfedges[2]);
      new_facets.append( newcfacet_ptr );
      FacetEntity *newfacet_ptr = dynamic_cast<FacetEntity *> (newcfacet_ptr);
      set_new_facet_owners( 2, facet_ptr, newfacet_ptr, surf_map, curve_map, point_map );
    } 
  }
  
  // make sure facets are oriented consistently on new volume
  if( new_facets.size() )
  {
    CubitFacet *start_facet = new_facets.get(); 
    CubitBoolean do_flip = CUBIT_TRUE;
    int nfacets;
    int mydebug = 0;
    rv = check_facet_orientation(start_facet, do_flip, nfacets, mydebug );
  }
  
  return rv;
}


CubitStatus ChollaEngine::detach_facet_edges(DLIList<ChollaCurve*> &chcurves,
                                             ChollaVolume *detaching_volume,                                             
                                             std::map<ChollaCurve *, ChollaCurve *> &curve_map,
                                             std::map<ChollaPoint *, ChollaPoint *> &point_map )
{  
  
  CubitStatus rv = CUBIT_SUCCESS;
  
  std::vector<CubitPoint *> new_points;
  std::vector<CubitFacetEdge *> new_edges;
  
  rv = copy_facet_edges_at_interface( chcurves, new_points, new_edges, curve_map, point_map );
  if (rv != CUBIT_SUCCESS)
    return rv;

  for (int i=chcurves.size(); i--; )
  {
    ChollaCurve *chcurv_ptr = chcurves.get_and_step();    
    rv = connect_points_at_interface( chcurv_ptr, detaching_volume, new_points );
    if (rv != CUBIT_SUCCESS)
      return rv;
    
    rv = connect_edges_at_interface( chcurv_ptr, detaching_volume, new_edges );
    if (rv != CUBIT_SUCCESS)
      return rv;
  }   
  
  return rv;
}


CubitStatus ChollaEngine::copy_facet_edges_at_interface(DLIList<ChollaCurve*> &chcurves,
                                                   std::vector<CubitPoint *> &new_points,
                                                   std::vector<CubitFacetEdge *> &new_edges,                                                                                                      
                                                   std::map<ChollaCurve *, ChollaCurve *> &curve_map,
                                                   std::map<ChollaPoint *, ChollaPoint *> &point_map)

{
  CubitStatus rv = CUBIT_SUCCESS;
  
  // first set the marked flags on all facets entities to -1  
  DLIList<FacetEntity *> facet_list;
  for( int i=chcurves.size(); i--; )
  {
    ChollaCurve *chcurv_ptr = chcurves.get_and_step();
    DLIList<FacetEntity*> tmp_facet_list = chcurv_ptr->get_facet_list();
    facet_list += tmp_facet_list;
  }

  FacetDataUtil::mark_facets(facet_list, FACET_ENTITY_UNINITIALIZED);
  
  // create a copy of each of the facet entities on the surface.  The marked flag in the facet will
  // keep track of the new entity created.  It will be a location in the new_points or new_edges
  // array.  Once we are finished with creating points and edges, we can create the new facets
  // given the references in the marked flags.
  
  // create new points

  std::map<ChollaSurface *, ChollaSurface *> dummy_surf_map;
  
  rv = copy_points_at_interface(facet_list, new_points, dummy_surf_map, curve_map, point_map);
  if (rv != CUBIT_SUCCESS)
    return rv;
  
  // create new edges
  
  rv = copy_edges_at_interface(facet_list, new_points, new_edges, dummy_surf_map, curve_map, point_map);
  if (rv != CUBIT_SUCCESS)
    return rv;  
  
  return rv;
}

//=============================================================================
//Function: copy_points_at_interface 
//Description: copy the points at the interface        
//Author:sjowen
//Date:  09/18/09
//=============================================================================
CubitStatus ChollaEngine::copy_points_at_interface(DLIList<FacetEntity *> &facet_list,
                                                   std::vector<CubitPoint *> &new_points,                                                   
                                                   std::map<ChollaSurface *, ChollaSurface *> &surf_map,
                                                   std::map<ChollaCurve *, ChollaCurve *> &curve_map,
                                                   std::map<ChollaPoint *, ChollaPoint *> &point_map)
{
  int iploc = 0;
  
  FacetEntity *fe_ptr, *newfe_ptr;
  CubitPoint *point_ptr, *newpoint_ptr;

  for (int ifacet = 0; ifacet<facet_list.size(); ifacet++)
  {
    FacetEntity *facet_ptr = facet_list.get_and_step();
    CubitFacet *cfacet_ptr = dynamic_cast<CubitFacet*>(facet_ptr);
    if( cfacet_ptr )
    {     
      for (int ii=0; ii<3; ii++)
      {
        point_ptr = cfacet_ptr->point(ii);        
        if (point_ptr->marked() == FACET_ENTITY_UNINITIALIZED)
        {
          newpoint_ptr = (CubitPoint *) new CubitPointData( point_ptr->x(), point_ptr->y(), point_ptr->z() );
          new_points.push_back(newpoint_ptr);
          point_ptr->marked(iploc++);
          fe_ptr = dynamic_cast<FacetEntity *> (point_ptr);
          newfe_ptr = dynamic_cast<FacetEntity *> (newpoint_ptr);
          set_new_facet_owners(0, fe_ptr, newfe_ptr, surf_map, curve_map, point_map );
        }   
      } 
    } 
    else
    {
      CubitFacetEdge *cfacet_edge_ptr = dynamic_cast<CubitFacetEdge*>(facet_ptr);      
      for (int ii=0; ii<2; ii++)
      {
        point_ptr = cfacet_edge_ptr->point(ii);
        if (point_ptr->marked() == FACET_ENTITY_UNINITIALIZED)
        {
          newpoint_ptr = (CubitPoint *) new CubitPointData( point_ptr->x(), point_ptr->y(), point_ptr->z() );
          new_points.push_back(newpoint_ptr);
          point_ptr->marked(iploc++);
          fe_ptr = dynamic_cast<FacetEntity *> (point_ptr);
          newfe_ptr = dynamic_cast<FacetEntity *> (newpoint_ptr);
          set_new_facet_owners(0, fe_ptr, newfe_ptr, surf_map, curve_map, point_map );
        }   
      } 
    }
  }
  
  return CUBIT_SUCCESS;
}

//=============================================================================
//Function: copy_edges_at_interface 
//Description: copy the edges at the interface        
//Author:sjowen
//Date:  09/18/09
//=============================================================================
CubitStatus ChollaEngine::copy_edges_at_interface(DLIList<FacetEntity *> &facet_list,
                                                  std::vector<CubitPoint *> &new_points,
                                                  std::vector<CubitFacetEdge *> &new_edges,                                                  
                                                  std::map<ChollaSurface *, ChollaSurface *> &surf_map,
                                                  std::map<ChollaCurve *, ChollaCurve *> &curve_map,
                                                  std::map<ChollaPoint *, ChollaPoint *> &point_map)
{
  int ieloc = 0;
  CubitFacetEdge *edge_ptr, *newedge_ptr;
  for (int ifacet = 0; ifacet<facet_list.size(); ifacet++)
  {
    FacetEntity *facet_ptr = facet_list.get_and_step();
    CubitFacet *cfacet_ptr = dynamic_cast<CubitFacet *> (facet_ptr);
    if( cfacet_ptr )
    {      
      for (int ii=0; ii<3; ii++)
      {
        edge_ptr = cfacet_ptr->edge(ii);
        if (edge_ptr->marked() == FACET_ENTITY_UNINITIALIZED)
        {
          CubitPoint *p0 = edge_ptr->point( 0 );
          CubitPoint *p1 = edge_ptr->point( 1 );
          int idx0 = p0->marked();
          int idx1 = p1->marked();
          CubitPoint *newp0 = new_points[idx0];
          CubitPoint *newp1 = new_points[idx1];
          newedge_ptr = (CubitFacetEdge *) new CubitFacetEdgeData( newp0, newp1 );
          new_edges.push_back(newedge_ptr);
          edge_ptr->marked(ieloc++);
          FacetEntity *fe_ptr = dynamic_cast<FacetEntity *> (edge_ptr);
          FacetEntity *newfe_ptr = dynamic_cast<FacetEntity *> (newedge_ptr);
          set_new_facet_owners( 1, fe_ptr, newfe_ptr, surf_map, curve_map, point_map );
        }  
      } 
    }
    else
    {
      CubitFacetEdge *edge_ptr = dynamic_cast<CubitFacetEdge*>(facet_ptr);
      if (edge_ptr->marked() == FACET_ENTITY_UNINITIALIZED)
      {
        CubitPoint *p0 = edge_ptr->point( 0 );
        CubitPoint *p1 = edge_ptr->point( 1 );
        int idx0 = p0->marked();
        int idx1 = p1->marked();
        CubitPoint *newp0 = new_points[idx0];
        CubitPoint *newp1 = new_points[idx1];
        newedge_ptr = (CubitFacetEdge *) new CubitFacetEdgeData( newp0, newp1 );
        new_edges.push_back(newedge_ptr);
        edge_ptr->marked(ieloc++);
        FacetEntity *fe_ptr = dynamic_cast<FacetEntity *> (edge_ptr);
        FacetEntity *newfe_ptr = dynamic_cast<FacetEntity *> (newedge_ptr);
        set_new_facet_owners( 1, fe_ptr, newfe_ptr, surf_map, curve_map, point_map );
      }
    }
  }
  return CUBIT_SUCCESS;
}

//=============================================================================
//Function: connect_facets_at_interface 
//Description: detach the facets from original points and edges and reattach to new copy         
//Author:sjowen
//Date:  09/18/09
//=============================================================================
CubitStatus ChollaEngine::connect_facets_at_interface(DLIList<ChollaSurface*> &chsurfs,
                                                      DLIList<ChollaCurve*> &chcurves,
                                                      ChollaVolume *chvol_ptr,
                                                      std::vector<CubitPoint *> &new_points,
                                                      std::vector<CubitFacetEdge *> &new_edges)
{
  CubitStatus rv = CUBIT_SUCCESS;  
  DLIList<ChollaCurve *> chcurv_list;
  
  for( int i=chsurfs.size(); i--; )
  {
    ChollaSurface *tmp_surf = chsurfs.get_and_step();    
    DLIList<ChollaCurve*> tmp_curves;
    tmp_surf->get_curves( tmp_curves );
    chcurv_list += tmp_curves;
  }

  chcurv_list += chcurves;
  chcurv_list.uniquify_unordered();

  for (int icrv = 0; icrv < chcurv_list.size(); icrv++)
  {
    ChollaCurve *chcurv_ptr = chcurv_list.get_and_step();    
    rv = connect_points_at_interface( chcurv_ptr, chvol_ptr, new_points );
    if (rv != CUBIT_SUCCESS)
      return rv;
    
    rv = connect_edges_at_interface( chcurv_ptr, chvol_ptr, new_edges );
    if (rv != CUBIT_SUCCESS)
      return rv;
  } 

  return rv;
}

//=============================================================================
//Function: connect_points_at_interface 
//Description: find adjacent CubitFacetEdges on volume 2 that needs updating 
//Notes: chcurv_ptr is a curve at the interface on the original volume
//       chvol_ptr is the second volume (new (copied) entities belong to vol 2)
//       This functions gets all edges attached to the original CubitPoints on the 
//       original curve 'chcurv_ptr'.  If any edge is on volume 2 (the volume we are 
//       splitting off) we need to update the edges on that volume to contain the
//       the new CubitPoint.
//Author:sjowen
//Date:  09/18/09
//=============================================================================
CubitStatus ChollaEngine::connect_points_at_interface(ChollaCurve *chcurv_ptr,
                                                      ChollaVolume *chvol_ptr,
                                                      std::vector<CubitPoint *> &new_points)
{
  DLIList<CubitPoint *> cp_list;
  chcurv_ptr->get_facet_points(cp_list, CUBIT_TRUE);
  for (int ip = 0; ip<cp_list.size(); ip++)
  {
    CubitPoint *cp_ptr = cp_list.get_and_step();
    CubitPoint *newcp_ptr = new_points[cp_ptr->marked()];
    
    // set the point into edges that are on volume 2.
    // Note that there is no direct reference from points to edges in our data structure 
    // so no need to add/remove the edge from the point
    
    DLIList<CubitFacetEdge *> pedge_list;
    cp_ptr->edges(pedge_list);
    for (int iedge=0; iedge<pedge_list.size(); iedge++)
    {
      CubitFacetEdge *edge_ptr = pedge_list.get_and_step();
      TDGeomFacet *td_geom = TDGeomFacet::get_geom_facet( edge_ptr );
      if (td_geom->is_in_volume( chvol_ptr ))
      {
        CubitPoint *p0 = edge_ptr->point(0);
        CubitPoint *p1 = edge_ptr->point(1);
        CubitFacetEdgeData *cfed_ptr = dynamic_cast<CubitFacetEdgeData *> (edge_ptr);
        if (p0 == cp_ptr)
          cfed_ptr->set_point( newcp_ptr, 0 );
        else if (p1 == cp_ptr)
          cfed_ptr->set_point( newcp_ptr, 1 );
      }
    }    
    
    // remove facets in volume 2 from the point and add the to the new point
    // update the point reference on the facet
    
    DLIList<CubitFacet *> pfacet_list;
    cp_ptr->facets(pfacet_list);
    for (int ifacet = 0; ifacet < pfacet_list.size(); ifacet++)
    {
      CubitFacet *facet_ptr = pfacet_list.get_and_step();
      TDGeomFacet *td_geom = TDGeomFacet::get_geom_facet( facet_ptr );
      if (td_geom->is_in_volume( chvol_ptr ))
      {
        cp_ptr->remove_facet( facet_ptr );        
        newcp_ptr->add_facet( facet_ptr );
        int ptidx = facet_ptr->point_index(cp_ptr);
        CubitFacetData *cfd_ptr = dynamic_cast<CubitFacetData *> (facet_ptr);
        cfd_ptr->set_point(newcp_ptr, ptidx);
      }
    }
  }
  return CUBIT_SUCCESS;
}

//=============================================================================
//Function: connect_edges_at_interface 
//Description: detach the facet edges from original facets and reattach to new copy         
//Notes:    chcurv_ptr is a curve at the interface on the original volume
//       chvol_ptr is the second volume (new (copied) entities belong to vol 2)
//       This functions gets all CubitFacets attached to any CubitFacetEdges on the  
//       original curve 'chcurv_ptr'.  If any facet is on volume 2 (the volume we are 
//       splitting off) we need to update the edges on that facet on the new volume
//       to contain the new CubitFacetEdge.
//Author:sjowen
//Date:  09/18/09
//=============================================================================
CubitStatus ChollaEngine::connect_edges_at_interface(ChollaCurve *chcurv_ptr,
                                                     ChollaVolume *chvol_ptr,
                                                     std::vector<CubitFacetEdge *> &new_edges)
{
  DLIList<FacetEntity *> fe_list = chcurv_ptr->get_facet_list();
  for (int ie=0; ie<fe_list.size(); ie++)
  {
    FacetEntity *fe_ptr = fe_list.get_and_step();
    CubitFacetEdge *edge_ptr = dynamic_cast<CubitFacetEdge *> (fe_ptr);
    assert(edge_ptr != NULL);
    DLIList<CubitFacet *> efacet_list;
    edge_ptr->facets(efacet_list);
    for (int ifacet=0; ifacet<efacet_list.size(); ifacet++)
    {
      CubitFacet *facet_ptr = efacet_list.get_and_step();
      TDGeomFacet *td_geom = TDGeomFacet::get_geom_facet( facet_ptr );
      if (td_geom->is_in_volume( chvol_ptr ))
      {
        edge_ptr->remove_facet( facet_ptr );
        CubitFacetEdge *newedge_ptr = new_edges[edge_ptr->marked()];
        newedge_ptr->add_facet( facet_ptr );
        int edidx = facet_ptr->edge_index(edge_ptr);
        CubitFacetData *cfd_ptr = dynamic_cast<CubitFacetData *> (facet_ptr);
        cfd_ptr->edge(newedge_ptr, edidx);
      }
    }
  }
  return CUBIT_SUCCESS;
}

//=============================================================================
//Function: set_new_facet_owners 
//Description: update the ownenrship of the new detached facet entity based upon
//             the map set up in detach_volumes
//            
//Author:sjowen
//Date:  09/14/09
//=============================================================================
CubitStatus ChollaEngine::set_new_facet_owners(int type, //0, 1, or 2 based on dimension of facet entity
                                               FacetEntity *fe_ptr, FacetEntity *newfe_ptr,                                                
                                               std::map<ChollaSurface *, ChollaSurface *> &surf_map,
                                               std::map<ChollaCurve *, ChollaCurve *> &curve_map,
                                               std::map<ChollaPoint *, ChollaPoint *> &point_map )
{
  
  // The tooldata on the original facet entity should tell us what cholla entity it belongs
  // to.  Using the map, we can then determine its new cholla entity partner.  With the
  // new merge partner, set the ownership of the new facet.  Note that this manages one-to-many
  // ownership cases where a facet lies on (or is owned by) any number of geometric entities
  
  TDGeomFacet *td_geom = TDGeomFacet::get_geom_facet( fe_ptr );
  TDGeomFacet::add_geom_facet(newfe_ptr, td_geom->get_block_id());
  
  // the original facet entity is owned by one or more surfaces
  
  DLIList<ChollaSurface *> surf_list;
  td_geom->get_cholla_surfs( surf_list );  
  for (int jj=0; jj<surf_list.size(); jj++)
  {
    ChollaSurface *surf_ptr = surf_list.get_and_step();
    std::map<ChollaSurface *, ChollaSurface *>::iterator map_it;
    map_it = surf_map.find( surf_ptr );
    assert(map_it != surf_map.end());
    ChollaSurface *newsurf_ptr = map_it->second;    
    TDGeomFacet *newtdgeom = TDGeomFacet::get_geom_facet( newfe_ptr );
    newtdgeom->add_cholla_surf( newsurf_ptr );
    
    // add this facet entity to the cholla surface only if this is a facet
    
    if (type == 2)
      newsurf_ptr->add_facet(newfe_ptr);
  }

  
  // the original facet entity is owned by one or more curves
  
  DLIList<ChollaCurve *> curv_list;
  td_geom->get_cholla_curves( curv_list );
  for (int jj=0; jj<curv_list.size(); jj++)
  {
    ChollaCurve *curv_ptr = curv_list.get_and_step();
    std::map<ChollaCurve *, ChollaCurve *>::iterator map_it, upper_iter;
    map_it = curve_map.find( curv_ptr );
    assert(map_it != curve_map.end());       
    ChollaCurve *newcurv_ptr = NULL;
    newcurv_ptr = map_it->second;
    
    TDGeomFacet *newtdgeom = TDGeomFacet::get_geom_facet( newfe_ptr );
    newtdgeom->add_cholla_curve( newcurv_ptr );
    
    // add this facet entity to the cholla curve only if this is a facetedge
    
    if (type == 1)
      newcurv_ptr->add_facet(newfe_ptr);      
  }
  
  // the original facet entity is owned by one or more points (vertices)
  
  DLIList<ChollaPoint *> point_list;
  td_geom->get_cholla_points( point_list );
  for (int jj=0; jj<point_list.size(); jj++)
  {
    ChollaPoint *point_ptr = point_list.get_and_step();
    std::map<ChollaPoint *, ChollaPoint *>::iterator map_it, upper_iter;
    map_it = point_map.find( point_ptr );
    assert(map_it != point_map.end());
    ChollaPoint *newchpoint_ptr = NULL;
    newchpoint_ptr = map_it->second;
  
    TDGeomFacet *newtdgeom = TDGeomFacet::get_geom_facet( newfe_ptr );
    newtdgeom->add_cholla_point( newchpoint_ptr );
    
    // add this facet entity to the cholla point only if this is a cubitpoint
    
    if (type == 0)
      newchpoint_ptr->add_facet(newfe_ptr);
  }
  return CUBIT_SUCCESS;

}

//=============================================================================
//Function: verify_points_to_curves 
//Description: verify the connectivity between points and curves           
//Author:sjowen
//Date:  09/16/09
//=============================================================================
CubitStatus ChollaEngine::verify_points_to_curves()
{
  for (int ii=0; ii<chollaPointList.size(); ii++)
  {
    ChollaPoint *chpt = chollaPointList.get_and_step();
    if (!chpt->verify_curves())
    {
      PRINT_ERROR("ChollaPoint %d not associated with one of its curves.\n", chpt->get_id());
      return CUBIT_FAILURE;
    }
  }
  for (int ii=0; ii<chollaCurveList.size(); ii++)
  {
    ChollaCurve *chcrv = chollaCurveList.get_and_step();
    if (!chcrv->verify_points())
    {
      PRINT_ERROR("ChollaCurve %d not associated with one of its points.\n", chcrv->get_id());
      return CUBIT_FAILURE;
    }
  }
  return CUBIT_SUCCESS;
}

CubitStatus ChollaEngine::set_curve_endpoints(std::map<ChollaPoint *, ChollaPoint *> &point_map,
                                              std::map<ChollaCurve *, ChollaCurve *> &curve_map )
                                              
{   
  DLIList<ChollaPoint*> original_chpts;
  DLIList<ChollaPoint*> new_chpts;
  DLIList<ChollaCurve*> original_chcurves;
  DLIList<ChollaCurve*> new_chcurves;  

  std::map<ChollaCurve*, ChollaCurve*>::iterator curv_iter = curve_map.begin();
  for(; curv_iter!=curve_map.end(); curv_iter++ )
  {
    ChollaCurve *orig_curve = curv_iter->first;
    ChollaCurve *new_curve = curv_iter->second;

    //curves
    original_chcurves.append( orig_curve );
    new_chcurves.append( new_curve );
    
    //points
    DLIList<ChollaPoint*> tmp_pts = orig_curve->get_points();
    original_chpts += tmp_pts;

    tmp_pts.clean_out();
    tmp_pts = new_curve->get_points();
    new_chpts += tmp_pts;
  }

  original_chpts.uniquify_unordered();
  original_chcurves.uniquify_unordered();
  new_chpts.uniquify_unordered();
  new_chcurves.uniquify_unordered();  
  
  std::map<ChollaPoint*,ChollaPoint*>::iterator pt_iter;
  for( int i=original_chpts.size(); i--; )
  {
    ChollaPoint *orig_pt = original_chpts.get_and_step();

    pt_iter = point_map.find( orig_pt );
    assert( pt_iter != point_map.end() );

    ChollaPoint *new_pt = pt_iter->second;    
    
    CubitPoint *new_cp = dynamic_cast<CubitPoint *>(new_pt->get_facets());
    CubitPoint *orig_cp = dynamic_cast<CubitPoint *>(orig_pt->get_facets());

    //get the curves attached to it that are in the original surface
    DLIList<ChollaCurve*> adj_curves = orig_pt->get_curves();

    for( int j=adj_curves.size(); j--; )
    {
      ChollaCurve *orig_curve = adj_curves.get_and_step();

      curv_iter = curve_map.find( orig_curve );

      //use the original curve's start/end pts as a guide to set the new curve's start/end pts
      if( curv_iter != curve_map.end() )
      {
        //if there is a corresponding new curve, get it
        ChollaCurve *new_curve_in_new_surf = curv_iter->second;

        CubitPoint *start, *end;
        orig_curve->get_ends(start, end);        
        assert(orig_cp != NULL);
        assert(start!= NULL);
        assert(end != NULL);

        if ( orig_cp == start )
        {
          new_curve_in_new_surf->set_start(new_cp);
        }        
        if ( orig_cp == end )
        {
          new_curve_in_new_surf->set_end(new_cp);
        }        
      }   
    }
  }

  //for each new point
  for( int i=new_chpts.size(); i--; )
  {
    ChollaPoint *new_point = new_chpts.get_and_step();

    DLIList<ChollaCurve*> curves_in_pt = new_point->get_curves();

    for( int k=curves_in_pt.size(); k--; )
    {
      ChollaCurve *curve_in_pt = curves_in_pt.get_and_step();

      //get adjacent curves that are not new (not in 'new_curves')
      if( !original_chcurves.is_in_list( curve_in_pt ) )
      {
        DLIList<ChollaPoint *> chpts_on_curve = curve_in_pt->get_points();        
        
        // one point on the curve assumes a periodic curve.  set both ends the same

        if (chpts_on_curve.size() == 1)
        {
          ChollaPoint *chpt_on_crv = chpts_on_curve.get();
          CubitPoint *cp = dynamic_cast<CubitPoint *> (chpt_on_crv->get_facets());
          curve_in_pt->set_start(cp);
          curve_in_pt->set_end(cp);
        }

        // standard case. one point of curve is on the interface and one is not.
        // In this case, one of the points has been replaced with a new point.
        // so one of the end point pointers are out of date.  Determine which one
        // and then set it.

        else if (chpts_on_curve.size() == 2)
        {
          ChollaPoint *chpt1_on_crv = chpts_on_curve.get_and_step();
          ChollaPoint *chpt2_on_crv = chpts_on_curve.get();
          CubitPoint *cp1 = dynamic_cast<CubitPoint *> (chpt1_on_crv->get_facets());
          CubitPoint *cp2 = dynamic_cast<CubitPoint *> (chpt2_on_crv->get_facets());
          CubitPoint *curstart, *curend;
          curve_in_pt->get_ends(curstart, curend);
          assert(curstart != NULL);
          assert(curend != NULL);

          if (curstart == cp1)
          {
            curve_in_pt->set_end(cp2);
          }
          else if (curstart == cp2)
          {
            curve_in_pt->set_end(cp1);
          }
          else if (curend == cp2)
          {
            curve_in_pt->set_start(cp1);
          }
          else if (curend == cp1)
          {
            curve_in_pt->set_start(cp2);
          }
          else
          {
            //do this by proximity
            double dist1 = curstart->coordinates().distance_between( cp1->coordinates() );
            double dist2 = curstart->coordinates().distance_between( cp2->coordinates() );

            if( dist1 < dist2 )
              curve_in_pt->set_start( cp1 );
            else if(dist2 < dist1 )
              curve_in_pt->set_start( cp2 );
            else
              assert(0);
          }
        }
        else
          assert(0);
      }
    } 
  }
  return CUBIT_SUCCESS;
}

// EOF
