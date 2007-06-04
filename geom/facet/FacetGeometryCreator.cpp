//- Class:       FacetGeometryCreator
//- Description: Creates the topology for a given geometry described by facets.
//- Owner:       Steven J. Owen
//- Checked by:
//- Version:
#include "CubitDefines.h"
#include "FacetGeometryCreator.hpp"
#include "FacetQueryEngine.hpp"
#include "DLIList.hpp"
#include "TDGeomFacet.hpp"
#include "CastTo.hpp"
#include "FacetSkinTool.hpp"
#include "FacetSurfaceMesh.hpp"
#include "FacetCurveMesh.hpp"
#include "FacetPointMesh.hpp"
#include "Body.hpp"
#include "Lump.hpp"
#include "ShellSM.hpp"
#include "Surface.hpp"
#include "Curve.hpp"
#include "GeometryQueryTool.hpp"
#include "GeometryModifyTool.hpp"
#include "PointSM.hpp"
#include "FacetShell.hpp"
#include "FacetSurface.hpp"
#include "CubitFacet.hpp"
#include "CubitFacetEdge.hpp"
#include "CubitPoint.hpp"
#include "FacetEntity.hpp"
#include "FacetEvalTool.hpp"
#include "debug.hpp"
#include "CurveSM.hpp"
#include "TDFacetBoundaryEdge.hpp"
#include "TDFacetBoundaryPoint.hpp"


//============================================================================
//Function:  FacetGeometryCreator (PUBLIC) (constructor)
//============================================================================
FacetGeometryCreator::FacetGeometryCreator()
{
  hashCurveArray = NULL;
  hashCurveSize = 0;
  hashPointArray = NULL;
  hashPointSize = 0;
}

//============================================================================
//Function:  FacetGeometryCreator (PUBLIC) (constructor)
//============================================================================
FacetGeometryCreator::FacetGeometryCreator(DLIList<FacetEntity*> &face_list,
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
}

//============================================================================
//Function:  set_up_tool_datas
//============================================================================
void FacetGeometryCreator::set_up_tool_datas( )
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

//============================================================================
//Function:  delete_tool_datas
//============================================================================
void FacetGeometryCreator::delete_tool_datas( )
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

//==================================================================================
//Function:  ~FacetGeometryCreator (PUBLIC) (destructor)
//==================================================================================
FacetGeometryCreator::~FacetGeometryCreator()
{
  delete_tool_datas();
  int ii;
    //clean up any data remaining.
  for (ii = facetPointList.size(); ii > 0; ii-- )
    delete facetPointList.remove();
  for (ii = facetCurveList.size(); ii > 0; ii-- )
    delete facetCurveList.remove();
  for (ii = facetSurfaceList.size(); ii > 0; ii-- )
    delete facetSurfaceList.remove();
  
}

//=============================================================================
//Function:  create_volume_boundaries (PRIVATE)
//Description: creates the surfaces based on the sideset and element block
//             information
//Author: sjowen
//Date: 10/17/00
//=============================================================================
CubitStatus FacetGeometryCreator::create_volume_boundaries( 
  DLIList<FacetSurfaceMesh*> &facet_surface_sheets,  // output global list of surfaces
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
  int num_surfaces = facet_surface_sheets.size();
  facet_surface_sheets.reset();
  for ( ii = num_surfaces; ii > 0; ii-- )
  {

    FacetSurfaceMesh *fsm_ptr = facet_surface_sheets.get_and_step();
    
    // make a list of feature edges
    
    DLIList<CubitFacetEdge *> feature_edge_list;
    if (use_feature_angle)
    {
      rv = fsm_ptr->feature_angle( min_dot, feature_edge_list );
    }
    else
    {
      rv = fsm_ptr->non_manifold_edges( feature_edge_list );
    }
    if (rv != CUBIT_SUCCESS)
      return rv;

     // crack the surface at the feature edges.  create new edges and
     // points so the facet representation is discontinuous.

    rv = FacetQueryEngine::make_features( feature_edge_list, split_surfaces );
    if (rv != CUBIT_SUCCESS)
       return rv;

    // split up the surface

    rv = fsm_ptr->split_surface( facet_surface_sheets );
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
    for ( ii = facet_surface_sheets.size(); ii > 0; ii-- )
    {
      FacetSurfaceMesh *fsm_ptr = facet_surface_sheets.get_and_step();
      fsm_ptr->clean_features();
    }
  }

  // Now that we've broken everything into surfaces, update the surface IDs
  // on the boundary facet tool data

  if (!split_surfaces)
  {
    for ( ii = facet_surface_sheets.size(); ii > 0; ii-- )
    {
      FacetSurfaceMesh *fsm_ptr = facet_surface_sheets.get_and_step();
      fsm_ptr->update_boundary_tool_data();
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
CubitStatus FacetGeometryCreator::create_surface_boundaries(
  DLIList<FacetSurfaceMesh*> &facet_surface_list, // global list of surfaces
  DLIList<FacetCurveMesh*> &facet_curve_list,     // output global list of curves
  CubitBoolean use_feature_angle,
  double min_dot )
{
  CubitStatus stat = CUBIT_SUCCESS;

    // determine the boundaries for each surface.  One curve per surface

  int ii;
  for ( ii = facet_surface_list.size(); ii > 0; ii-- )
  {
    FacetSurfaceMesh *fsm_ptr = facet_surface_list.get_and_step();
    DLIList<FacetCurveMesh*> fcurve_list;
    fsm_ptr->get_curves( fcurve_list );
    if (fcurve_list.size() == 0)
    {
      DLIList<FacetEntity*> facet_list;
      fsm_ptr->get_facets(facet_list);
      FacetSkinTool f_skin_tool;
      stat = f_skin_tool.skin_2d(facet_list, fsm_ptr);
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
  // Create a new FacetCurveMesh for each curve
  // Curves are created wherever there is a unique set of associated
  // surfaces 

  int jj, kk;
  for ( ii = facet_surface_list.size(); ii > 0; ii-- )
  {

    FacetSurfaceMesh *fsm_ptr = facet_surface_list.get_and_step();
    DLIList<FacetCurveMesh*> fcurve_list;
    fsm_ptr->get_curves( fcurve_list );

    // curently there should only be one curve list per surface

    for (jj=fcurve_list.size(); jj>0; jj--)
    {
      FacetCurveMesh *fcm_ptr = fcurve_list.get_and_step();
      DLIList<FacetEntity*> facet_list =  fcm_ptr->get_facets();
      FacetEntity *edge_ptr;
      for ( kk = 0; kk < facet_list.size(); kk++)
      {
        edge_ptr = facet_list.get_and_step();
        stat = classify_edge( edge_ptr, facet_curve_list, fsm_ptr );
        if (stat != CUBIT_SUCCESS) 
          return stat;
      }

      // delete this FacetCurveMesh - it should have been replaced by one
      // or more curves bounding this surface

      fsm_ptr->remove_curve( fcm_ptr );
    }
  }
  delete_hash_curves();

  // Split these curves so that there is only one string of continuous
  // edges per curve (it will also order the edges and set the start
  // and end nodes for each curve) 

  int num_curves = facet_curve_list.size();

  facet_curve_list.reset();
  for ( ii = num_curves; ii > 0 && stat == CUBIT_SUCCESS; ii-- )
  {
    FacetCurveMesh *fcm_ptr = facet_curve_list.get();

    // if necessary mark nodes that will serve as feature breaks (vertices 
    // will be generatedat them)

    if (use_feature_angle)
    {
      stat = fcm_ptr->feature_angle( min_dot );
      if (stat != CUBIT_SUCCESS)
         return stat;
    }

    // split the curve based on various criteria

    stat = fcm_ptr->split_curve( facet_curve_list );

    // delete this curve (new curves were created in split_curve)

    delete fcm_ptr;
    facet_curve_list.change_to(NULL);
    facet_curve_list.step();
  }
  facet_curve_list.remove_all_with_value(NULL);

  // update the point->curve associativity

  for (ii=0; ii<facet_curve_list.size(); ii++)
  {
    FacetCurveMesh *fcm_ptr = facet_curve_list.get_and_step();
    CubitPoint *start_ptr, *end_ptr;
    fcm_ptr->get_ends( start_ptr, end_ptr );
    TDGeomFacet *td = TDGeomFacet::get_geom_facet( start_ptr );
    td->add_facet_curve( fcm_ptr );
    td = TDGeomFacet::get_geom_facet( end_ptr );
    td->add_facet_curve( fcm_ptr );
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
CubitStatus FacetGeometryCreator::init_hash_curves( )
{
  /* === find the next highest prime number */

  int num_surfs = facetSurfaceList.size();
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
  hashCurveArray = new DLIList<FacetCurveMesh*>[hashCurveSize];

  int key = 0;
  FacetCurveMesh *fcm_ptr;
  DLIList<FacetSurfaceMesh*> *fsm_list_ptr; 
  for (i=0; i<facetCurveList.size(); i++)
  {
    fcm_ptr = facetCurveList.get_and_step();
    fsm_list_ptr = fcm_ptr->get_surface_list_ptr();
    key = get_curve_hash_key( fsm_list_ptr );
    hashCurveArray[key].append( fcm_ptr );   
  }
  return CUBIT_SUCCESS;
}

//=============================================================================
//Function:  delete_hash_curves (PRIVATE)
//Description: delete the hash curve stuff
//Author: sjowen
//Date: 3/7/01
//=============================================================================
void FacetGeometryCreator::delete_hash_curves( )
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
int FacetGeometryCreator::get_curve_hash_key(
  DLIList<FacetSurfaceMesh*> *fsm_list_ptr )
{
  int key, j;
  FacetSurfaceMesh *fsm_ptr;
  if (fsm_list_ptr->size() == 0)
  {
    key = 0;
  }
  else
  {
    key = INT_MAX;
    for (j=0; j<fsm_list_ptr->size(); j++)
    {
      fsm_ptr = fsm_list_ptr->get_and_step();
      if (fsm_ptr->get_id() < key)
        key = fsm_ptr->get_id();
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
CubitStatus FacetGeometryCreator::classify_edge(
  FacetEntity *edge_ptr,       // the edge we are classifying
  DLIList<FacetCurveMesh*> &facet_curve_list,  // add to one of these
  FacetSurfaceMesh *facet_surf_mesh_ptr )   // the current surface
{
  CubitStatus rv = CUBIT_SUCCESS;

  // see if we have already classified this edge (from another surface)

  TDGeomFacet *td_gm_edge = TDGeomFacet::get_geom_facet(edge_ptr);
  if (td_gm_edge->get_hit_flag() != 0)
    return rv;
  td_gm_edge->set_hit_flag(1);

  // get the surfaces adjacent to this edge

  DLIList<FacetSurfaceMesh*> this_fsurf_list;
  td_gm_edge->get_facet_surfs( this_fsurf_list );
  int this_num_adj = this_fsurf_list.size();

  // see if the surfaces defined on this face match any 
  // of the existing block curves

  DLIList<FacetSurfaceMesh*> *fsurf_list_ptr;
  FacetCurveMesh *fcm_ptr = NULL;
  int found = 0;
  int key = get_curve_hash_key( &this_fsurf_list );
  int ii;
  for (ii=0; ii<hashCurveArray[key].size() && !found; ii++)
  {
    fcm_ptr = hashCurveArray[key].get();
    // the first one checked should be the same as the last one checked (don't
    // use get_and_step here)  This should speed things up

    // check if surfaces are the same

    fsurf_list_ptr = fcm_ptr->get_surface_list_ptr( );
    int num_adj = fsurf_list_ptr->size();
    if (num_adj == this_num_adj)
    {
      found = 1;
      int jj, kk;
      for (jj=fsurf_list_ptr->size(); jj>0 && found; jj--)
      {
        int same_surf = 0;
        FacetSurfaceMesh *fsm_ptr = fsurf_list_ptr->get_and_step();
        for(kk=this_fsurf_list.size(); kk>0 && !same_surf; kk--)
        {
          FacetSurfaceMesh *this_fsm_ptr = this_fsurf_list.get_and_step();
          if (this_fsm_ptr == fsm_ptr)
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

    int was_added = fcm_ptr->add_facet_unique( edge_ptr );

    // add the curve to the surface (if needed)

    facet_surf_mesh_ptr->add_curve_unique( fcm_ptr );

    // add the curve to the edge
  
    td_gm_edge->add_facet_curve( fcm_ptr );
  }

  // if the unique set of surfaces that this edge is associated
  // with is not found, then create a new facet curve and add the edge to it

  else
  {

    // create it and update surface and nodeset info

    int block_id = td_gm_edge->get_block_id();
    FacetCurveMesh *new_fcm_ptr = new FacetCurveMesh( block_id ); 
    for (int mm=0; mm<this_num_adj; mm++)
    {
      new_fcm_ptr->add_surface( this_fsurf_list.get_and_step() );
    }

    // add the edge

    new_fcm_ptr->add_facet( edge_ptr );

    // update the surface with this new curve

    facet_surf_mesh_ptr->add_curve( new_fcm_ptr );

    // add the new curve to the global list

    facet_curve_list.append( new_fcm_ptr ); 

    // add the curve to the edge
  
    td_gm_edge->add_facet_curve( new_fcm_ptr );

    // add the new curve to the hash table

    hashCurveArray[key].append( new_fcm_ptr );
  }
  return rv;
}

//=============================================================================
//Function:  create_curve_boundaries (PRIVATE)
//Description: creates the points based on the nodeset and curve information
//Author: sjowen
//Date: 12/6/00
//=============================================================================
CubitStatus FacetGeometryCreator::create_curve_boundaries(
  DLIList<FacetCurveMesh*> &facet_curve_list,   // global list of curves
  DLIList<FacetPointMesh*> &facet_point_list )  // output global list of points
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
  // Create a new FacetPointMesh for each point

#ifdef BOYD17 
  int mydebug = 0;
#endif 
  int ii, kk;
  for ( ii = facet_curve_list.size(); ii > 0; ii-- )
  {
    FacetCurveMesh *fcm_ptr = facet_curve_list.get_and_step();
    CubitPoint *point_ptr[2];
    fcm_ptr->get_ends( point_ptr[0], point_ptr[1] );
    for ( kk = 0; kk < 2; kk++)
    {
      stat = classify_point( point_ptr[kk], facet_point_list, fcm_ptr );
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
//Function:  classify_node (PRIVATE)
//Description: sorts a node into its correct point based on its associated
//             curve.  Creates a new block point if necessary
//Author: sjowen
//Date: 12/6/00
//=============================================================================
CubitStatus FacetGeometryCreator::classify_point(
  CubitPoint *point_ptr,                     // the node to classify
  DLIList<FacetPointMesh*> &facet_point_list,   // global list of points
  FacetCurveMesh *fcm_ptr )                 // curve that the end point is on
{
  int ii;
  int found = 0;
  TDGeomFacet *td_node = TDGeomFacet::get_geom_facet( point_ptr );
  FacetPointMesh *fpm_ptr;
  DLIList<FacetCurveMesh*> fcm_list;
  td_node->get_facet_curves(fcm_list);
  int key = get_point_hash_key( &fcm_list );
  for (ii = 0; ii < hashPointArray[key].size() && !found; ii++)
  {
    fpm_ptr = hashPointArray[key].get();
    FacetEntity *this_point_ptr = fpm_ptr->get_facets();
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
    fpm_ptr->add_curve( fcm_ptr );
    fcm_ptr->add_point( fpm_ptr );
  }
  else
  {
    FacetPointMesh *new_fpm_ptr = new FacetPointMesh(); 
    new_fpm_ptr->add_facet( point_ptr );
    new_fpm_ptr->add_curve( fcm_ptr );
    fcm_ptr->add_point( new_fpm_ptr );
    facet_point_list.append( new_fpm_ptr );
    hashPointArray[key].append( new_fpm_ptr );
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
CubitStatus FacetGeometryCreator::init_hash_points( )
{
  /* === find the next highest prime number */

  int num_curves = facetCurveList.size();
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
  hashPointArray = new DLIList<FacetPointMesh*>[hashPointSize];

  int key = 0;
  FacetPointMesh *fpm_ptr;
  DLIList<FacetCurveMesh*> *fcm_list_ptr; 
  for (i=0; i<facetPointList.size(); i++)
  {
    fpm_ptr = facetPointList.get_and_step();
    fcm_list_ptr = fpm_ptr->get_curve_list_ptr();
    key = get_point_hash_key( fcm_list_ptr );
    hashPointArray[key].append( fpm_ptr );   
  }
  return CUBIT_SUCCESS;
}

//=============================================================================
//Function:  delete_hash_points (PRIVATE)
//Description: delete the hash point stuff
//Author: sjowen
//Date: 3/7/01
//=============================================================================
void FacetGeometryCreator::delete_hash_points( )
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
int FacetGeometryCreator::get_point_hash_key(
  DLIList<FacetCurveMesh*> *fcm_list_ptr )
{
  int key, j;
  FacetCurveMesh *fcm_ptr;
  if (fcm_list_ptr->size() == 0)
  {
    key = 0;
  }
  else
  {
    key = INT_MAX;
    for (j=0; j<fcm_list_ptr->size(); j++)
    {
      fcm_ptr = fcm_list_ptr->get_and_step();
      if (fcm_ptr->get_id() < key)
        key = fcm_ptr->get_id();
    }
  }
  key = key % hashPointSize;
  return key;
}


//===============================================================================
//Function:  build_geometry (PRIVATE)
//Description:  build the CUBIT geometry based on the Facet entity class lists
//===============================================================================
CubitStatus FacetGeometryCreator::build_geometry(
  DLIList<FacetSurfaceMesh*> &facet_surface_list,
  DLIList<FacetCurveMesh*> &facet_curve_list,
  DLIList<FacetPointMesh*> &facet_point_list,
  int interp_order,
  CubitBoolean use_feature_angle,
  double min_dot,
  CubitBoolean smooth_non_manifold,
  CubitBoolean split_surfaces )
{
  CubitStatus stat = CUBIT_SUCCESS;

  stat = build_point_geometry( facet_point_list );
  
  if (stat == CUBIT_SUCCESS)
    stat = build_curve_geometry( facet_curve_list );

  if (stat == CUBIT_SUCCESS)
    stat = build_surface_geometry( facet_surface_list,
                                   interp_order, min_dot );
  if (stat == CUBIT_SUCCESS)
    if (interp_order == 4)
      stat = clean_geometry( smooth_non_manifold, split_surfaces,
                           use_feature_angle, min_dot, facet_curve_list );
  return stat;
}

//===============================================================================
//Function:  build_point_geometry (PRIVATE)
//Description:  From the facet point list, create geometric points for each
//===============================================================================
CubitStatus FacetGeometryCreator::build_point_geometry(
  DLIList<FacetPointMesh*> &facet_point_list )
{
  CubitStatus stat = CUBIT_SUCCESS;
  int kk;
#ifdef BOYD17 
  int mydebug = 0;
#endif
  for ( kk = facet_point_list.size(); kk > 0; kk-- )
  {
    FacetPointMesh *fpm_ptr = facet_point_list.get_and_step();
    Point *point_ptr = fpm_ptr->get_geometric_point();
    if (point_ptr == NULL)
    {
      FacetEntity *facet_ptr = fpm_ptr->get_facets();
      CubitPoint *cp_ptr = CAST_TO( facet_ptr, CubitPoint );
      PointSM *point = NULL;
      stat = GeometryModifyTool::instance()->make_facet_point( cp_ptr, point );
      if ( point == NULL || stat != CUBIT_SUCCESS )
      {
        PRINT_ERROR("Problems building mesh based points.\n");
        return stat;
      }
      point_ptr = (Point *)point;
      fpm_ptr->assign_geometric_point(point_ptr);
    }
  }
  return stat;
}

//===============================================================================
//Function:  build_curve_geometry (PRIVATE)
//Description:  From the facet curve list, create geometric curves for each
//===============================================================================
CubitStatus FacetGeometryCreator::build_curve_geometry(
  DLIList<FacetCurveMesh*> &facet_curve_list )
{
  CubitStatus stat = CUBIT_SUCCESS;
  int kk;
  for ( kk = facet_curve_list.size(); kk > 0; kk-- )
  {
    FacetCurveMesh *fcm_ptr = facet_curve_list.get_and_step();
    Curve *curv_ptr = fcm_ptr->get_geometric_curve();
    if (curv_ptr == NULL)
    {
      DLIList<FacetPointMesh*> fpoint_list = fcm_ptr->get_points( );
      FacetPointMesh *fpm0_ptr = fpoint_list.get_and_step();
      FacetPointMesh *fpm1_ptr = fpoint_list.get_and_step();
      CubitPoint *start_point, *end_point;
      fcm_ptr->get_ends( start_point, end_point );
      if (fpm0_ptr->get_facets() != start_point)
      {
        FacetPointMesh *temp_ptr;
        temp_ptr = fpm0_ptr;
        fpm0_ptr = fpm1_ptr;
        fpm1_ptr = temp_ptr;
      }
      PointSM *start_ptr = (PointSM *) fpm0_ptr->get_geometric_point();
      PointSM *end_ptr = (PointSM *) fpm1_ptr->get_geometric_point();
      
      // if this is a curve without a parent surface then handle it 
      // differently.  (Curves with parents use the surface to evaluate to
      // With only a curve, it must evaluate to the curve)
      
      DLIList<FacetSurfaceMesh*> fsm_list = fcm_ptr->get_surfaces();
      if (fsm_list.size() == 0)
      {
        DLIList<FacetEntity*> facet_list;
        DLIList<CubitPoint*> point_list;  // needs to be filled in
        facet_list = fcm_ptr->get_facets();
        DLIList<CubitFacetEdge*> edge_list;
        CAST_LIST( facet_list, edge_list, CubitFacetEdge );
        if (stat != CUBIT_SUCCESS)
          return stat;
        stat = FacetQueryEngine::instance()->make_facet_curve( start_ptr, end_ptr,
                                                                  edge_list, point_list,
                                                                  curv_ptr );
        if (stat == CUBIT_SUCCESS)
        {
          CurveSM *curvsm_ptr = CAST_TO( curv_ptr, CurveSM );
          GeometryQueryTool::instance()->make_free_RefEdge( curvsm_ptr );
        }
      }
      else
      {
        stat = FacetQueryEngine::instance()->make_facet_curve( start_ptr, end_ptr,
                                                                  curv_ptr );
        if ( curv_ptr == NULL || stat != CUBIT_SUCCESS )
        {
          PRINT_ERROR("Problems building mesh based curves.\n");
          return stat;
        }
      }
      fcm_ptr->assign_geometric_curve(curv_ptr);
    }
  }
  return stat;
}

//===============================================================================
//Function:  build_loop_geometry (PRIVATE)
//Description:  From the block curve list of a surface, create geometric loops 
//===============================================================================
CubitStatus FacetGeometryCreator::build_loop_geometry(
  DLIList<FacetCurveMesh*> &facet_curve_list,   // curves on a surface
  FacetSurfaceMesh *fsm_ptr,                    // the surface
  DLIList<LoopSM*> &loop_list,                  // append to this loop list
  int &ncurves )                                // number of curves used
{

  // find the first unused curve on the list.  Use that as the starting curve
  // in the loop

  ncurves = 0;
  CubitStatus stat = CUBIT_SUCCESS;
  int ii, jj;
  int found = 0;
  FacetCurveMesh *fcm_ptr = NULL;
  for(ii=0; ii<facet_curve_list.size() && !found; ii++)
  {
    fcm_ptr = facet_curve_list.get_and_step();
    if (fcm_ptr->get_flag() == 0)
      found = 1;
  }
  if (!found)
    return CUBIT_FAILURE;  // all the curves have already been used in this surface

  // get the first mesh edge on the first curve. use it to determine the 
  // orientation of the curve with respect to the surface

  DLIList<FacetEntity*> facet_list = fcm_ptr->get_facets();
  FacetEntity *facet_ptr = facet_list.get();
  CubitFacetEdge *edge_ptr = CAST_TO( facet_ptr, CubitFacetEdge );
  if (!edge_ptr)
    return CUBIT_FAILURE; // facet edges are not on the list ????

  // get the adjacent face on this surface

  DLIList<FacetEntity*> adj_face_list;
  FacetEntity *face_ptr = NULL;
  edge_ptr->get_parents( adj_face_list );
  found = 0;
  for ( jj = 0; jj < adj_face_list.size() && !found; jj++ )
  {
    face_ptr = adj_face_list.get_and_step();
    TDGeomFacet *td_gm_face = TDGeomFacet::get_geom_facet(face_ptr); 
    DLIList<FacetSurfaceMesh*> fsurf_list;
    td_gm_face->get_facet_surfs( fsurf_list );
    FacetSurfaceMesh *face_fsm_ptr = fsurf_list.get();
    if(face_fsm_ptr == fsm_ptr)
    {
      found = 1;
    }
  }
  if (!found)
    return CUBIT_FAILURE;  // didn't find an adj face on the surface ???
    
  // determine orientation of nodes on this mesh face

  CubitPoint *start_ptr, *end_ptr;
  fcm_ptr->get_ends( start_ptr, end_ptr );
  end_ptr = edge_ptr->other_point( start_ptr );
  if (end_ptr == NULL)
    return CUBIT_FAILURE;  // the edge list may not be ordered correctly??
  CubitPoint *points[3];
  CubitFacet *tri_ptr = CAST_TO( face_ptr, CubitFacet );
  tri_ptr->points( points[0], points[1], points[2] );

  found = 0;
  CubitSense orientation;
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
  fcm_ptr->get_ends( start_ptr, end_ptr );

  // create a new coedge

  CoEdgeSM *coedge_ptr = NULL;
  Curve *curv_ptr = fcm_ptr->get_geometric_curve();
  stat = GeometryModifyTool::instance()->make_facet_coedge( (CurveSM*)curv_ptr, orientation,
                                                      coedge_ptr );
  if ( coedge_ptr == NULL || stat != CUBIT_SUCCESS )
  {
    PRINT_ERROR("Problems building mesh based coedges.\n");
    return stat;
  }

  // start a list of coedges for this loop

  DLIList<CoEdgeSM*> coedge_list;
  coedge_list.append( coedge_ptr );
  fcm_ptr->set_flag(1);
  FacetCurveMesh *last_fcm_ptr = fcm_ptr;

  // loop through and determine the rest of the coedges based on the
  // orientation of the first curve in the loop
   
  fcm_ptr->get_ends( start_ptr, end_ptr );
  FacetEntity *n0_ptr, *n1_ptr;
  if (orientation == CUBIT_FORWARD)
  {
    n0_ptr = start_ptr;
    n1_ptr = end_ptr;
  }
  else
  {
    n1_ptr = start_ptr;
    n0_ptr = end_ptr;
  }

  while (n0_ptr != n1_ptr)
  {

    // find the next curve in the loop

    int ntimes = 0;
    found = 0;
    while(!found)
    {
      for ( jj = 0; jj < facet_curve_list.size() && !found; jj++)
      {
        fcm_ptr = facet_curve_list.get_and_step();
        if (last_fcm_ptr != fcm_ptr && (fcm_ptr->get_flag() == 0 || ntimes > 0))
        {
          fcm_ptr->get_ends( start_ptr, end_ptr );
          if (start_ptr == n1_ptr)
          {
            found = 1;
            orientation = CUBIT_FORWARD;
            n1_ptr = end_ptr;
          }
          else if (end_ptr == n1_ptr)
          {
            found = 1;
            orientation = CUBIT_REVERSED;
            n1_ptr = start_ptr;
          }
        }
      }
      if (!found)
      {

        // in this special case - we have gone through the whole list without
        // finding the next curve on the loop.  In the first pass we have assumed
        // that each curve is only used once.  This may not be the case.  Go
        // back and see of we can find a curve that has already been used that
        // will work.  If we fail the second time through, then we fail!

        ntimes++;
        if (ntimes > 1)
          // return CUBIT_FAILURE;  // didn't find the next curve in the loop
          // DEBUG ONLY -- THIS NEEDS TO BE FIXED
          goto create_loop;
      }
    }

    // create a new coedge

    coedge_ptr = NULL;
    curv_ptr = fcm_ptr->get_geometric_curve();
    stat = GeometryModifyTool::instance()->make_facet_coedge( (CurveSM*)curv_ptr, orientation,
                                                         coedge_ptr );
    if ( coedge_ptr == NULL || stat != CUBIT_SUCCESS )
    {
      PRINT_ERROR("Problems building mesh based coedges.\n");
      return stat;
    }
    coedge_list.append( coedge_ptr );
    if (coedge_list.size() > 2*facet_curve_list.size())
      return CUBIT_FAILURE;  // more curves in the loop than their are curves ????
    fcm_ptr->set_flag(1);
    last_fcm_ptr = fcm_ptr;
  } 

  // create the loop

create_loop:
  LoopSM *loop_ptr = NULL;
  stat = GeometryModifyTool::instance()->make_facet_loop( coedge_list, loop_ptr );
  if ( coedge_ptr == NULL || stat != CUBIT_SUCCESS )
  {
    PRINT_ERROR("Problems building mesh based loops.\n");
    return stat;
  }
  loop_list.append( loop_ptr );
  ncurves = coedge_list.size();
      
  return stat;
}

//===============================================================================
//Function:  build_surface_geometry (PRIVATE)
//Description:  From the facet surface list, create geometric surface,
//              loops and coedges for each surface in the list
//===============================================================================
CubitStatus FacetGeometryCreator::build_surface_geometry(
  DLIList<FacetSurfaceMesh*> &facet_surface_list,
  int interp_order,
  double min_dot )
{
  CubitStatus stat = CUBIT_SUCCESS;
  int ii, kk;

  // make sure the facet flags have been reset

  FacetCurveMesh *fcm_ptr;
  for ( kk = facet_surface_list.size(); kk > 0; kk-- )
  {
    FacetSurfaceMesh *fsm_ptr = facet_surface_list.get_and_step();
    fsm_ptr->reset_facet_flags();
  }

  // now loop through surfaces and create them

  int mydebug = 0;
  for ( kk = facet_surface_list.size(); kk > 0; kk-- )
  {
    int num_curves = 0;
    FacetSurfaceMesh *fsm_ptr = facet_surface_list.get_and_step();

    DLIList<FacetCurveMesh*> fcurve_list;
    fsm_ptr->get_curves( fcurve_list );

    // init flags for curves on this surface

    for (ii=0; ii< fcurve_list.size(); ii++)
    {
      fcm_ptr = fcurve_list.get_and_step();
      fcm_ptr->set_flag(0);
      if (mydebug)
      {
        fcm_ptr->draw();
      }
    }

    // generate loops - keep going until we have used all the curves in the list

    DLIList<LoopSM*> loop_list;
    while (num_curves < fcurve_list.size())
    {
      int ncurves = 0;
      stat = build_loop_geometry( fcurve_list, fsm_ptr, loop_list, ncurves );
      if (stat != CUBIT_SUCCESS)
      {
        PRINT_ERROR("Can't build geometric loop from facets\n");
        return stat;
      }
      num_curves += ncurves;
    }

    // create the surface

    Surface *surf_ptr = fsm_ptr->get_geometric_surface();
    if (surf_ptr == NULL)
    {
      DLIList<FacetEntity*> facet_entity_list;
      DLIList<CubitPoint*> point_list;
      fsm_ptr->get_points(point_list);
      fsm_ptr->get_facets(facet_entity_list);
      DLIList<CubitFacet*> facet_list;
      CAST_LIST( facet_entity_list, facet_list, CubitFacet );
      stat = GeometryModifyTool::instance()->make_facet_surface(facet_list,
                                                                point_list,
                                                                loop_list,
                                                                interp_order,
                                                                min_dot,
                                                                surf_ptr);
      if ( surf_ptr == NULL || stat != CUBIT_SUCCESS )
      {
        PRINT_ERROR("Problems building mesh based surfaces.\n");
        return stat;
      }
      fsm_ptr->assign_geometric_surface(surf_ptr);
    }
  }
  return stat;
}

//===============================================================================
//Function:  clean_geometry (PRIVATE)
//Description:  fix the edge control points and the normals so they are conforming
//              (or non-conforming) accross curves
//===============================================================================
CubitStatus FacetGeometryCreator::clean_geometry( 
  CubitBoolean smooth_non_manifold,
  CubitBoolean split_surfaces,
  CubitBoolean use_feature_angle,
  double mindot,
  DLIList <FacetCurveMesh *> &facet_curve_list )
{
  int iedge;
  FacetCurveMesh *fcm_ptr;
  FacetEntity *fedge_ptr;
  CubitFacetEdge *edge_ptr;
  DLIList <FacetEntity *> facet_list; 
  DLIList<CubitFacetEdge *> feature_edge_list;
  int icurve;
  for (icurve=0; icurve<facet_curve_list.size(); icurve++)
  {
    fcm_ptr = facet_curve_list.get_and_step();
    facet_list.clean_out();
    facet_list = fcm_ptr->get_facets( );
    for (iedge=0; iedge<facet_list.size(); iedge++)
    {
      fedge_ptr = facet_list.get_and_step();
      edge_ptr = CAST_TO( fedge_ptr, CubitFacetEdge );
      feature_edge_list.append( edge_ptr );
    }
  }

  return FacetQueryEngine::instance()->fix_geometry( smooth_non_manifold,
                                                        split_surfaces, 
                                                        use_feature_angle, 
                                                        mindot, feature_edge_list );
}

//===============================================================================
//Function:  print_me (PRIVATE)
//Description: debug
//===============================================================================
void FacetGeometryCreator::print_me()
{
  FILE *fp = fopen( "fbg.out", "w" );
  int ii, jj, kk;
  for (ii=0; ii<facetSurfaceList.size(); ii++)
  {
    FacetSurfaceMesh *fsm_ptr = facetSurfaceList.get_and_step();
    fprintf( fp, "*********** Surface %d ***************\n", fsm_ptr->get_id());
    DLIList<FacetCurveMesh *> fcm_list;
    fsm_ptr->get_curves( fcm_list );
    for (jj=0; jj<fcm_list.size(); jj++)
    {
      FacetCurveMesh *fcm_ptr = fcm_list.get_and_step();
      fprintf(fp, "   Curve %d\n", fcm_ptr->get_id() );
#ifdef BOYD17 
      CubitVector start, end;
#endif
      DLIList<FacetPointMesh *> fpm_list = fcm_ptr->get_points();
      for (kk=0; kk<fpm_list.size(); kk++)
      {
        FacetPointMesh *fpm_ptr = fpm_list.get_and_step();
        FacetEntity *facet_point = fpm_ptr->get_facets();
        if (facet_point)
        {
          CubitPoint *point_ptr = CAST_TO(facet_point, CubitPoint);
          CubitVector pt = point_ptr->coordinates();
          fprintf(fp, "      (%d) x = %8.4lf y = %8.4lf z = %8.4lf\n",
                  fpm_ptr->get_id(), pt.x(), pt.y(), pt.z());
        }
      }
    }
  }
  fclose(fp);
}

