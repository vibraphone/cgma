//- Class:       ChollaSurface
//- Description: Temporary class for constructing the facet-based geometry
//-              
//- Owner:       Steven J. Owen
//- Checked by:
//- Version:
#include "ChollaSurface.hpp"
#include "ChollaCurve.hpp"
#include "TDGeomFacet.hpp"
#include "CastTo.hpp"
#include "CubitPoint.hpp"
#include "CubitFacet.hpp"
#include "CubitFacetData.hpp"
#include "CubitFacetEdge.hpp"
#include "CubitPointData.hpp"
#include "CubitFacetEdgeData.hpp"
#include "debug.hpp"
#include "TDFacetBoundaryPoint.hpp"
#include "TDFacetBoundaryEdge.hpp"
#include "GfxDebug.hpp"
#include "FacetEvalTool.hpp"
#include "FacetDataUtil.hpp"

//===============================================================================
//Function:  ChollaSurface (PUBLIC) (constructor)
//===============================================================================
ChollaSurface::ChollaSurface(int block_id)
{
  static int count = 100;
  id = count++;
  myFlag = CUBIT_FALSE;
  mySurface = NULL;
  myEvalTool = NULL;
  blockId = block_id;
  myMergePartner = NULL;
}
//===============================================================================
//Function:  ~ChollaSurface (PUBLIC) (destructor)
//===============================================================================
ChollaSurface::~ChollaSurface()
{
}

//===========================================================================
//Function Name: check_faceting
//
//Member Type:  PRIVATE
//Descriptoin:  check the edge/face orientations  
//===========================================================================
void ChollaSurface::check_faceting()
{

}


//=============================================================================
//Function:  split_surface (PUBLIC)
//Description: split this surface into multiple ChollaSurface where there are
//             discontinuous faces.  
//Author: sjowen
//Date: 12/22/00
//=============================================================================
CubitStatus ChollaSurface::split_surface( 
  DLIList<ChollaSurface*> &facet_surface_list 
  )
{
  DLIList<ChollaSurface*> new_surface_list;
  CubitStatus stat = CUBIT_SUCCESS;

  // Go through the surfaceElemList and pull faces off one by one as we
  // determine which surface it belongs to.  Continue until we have depleted
  // the list

  int jj;
  int mydebug = 0;
  int num_surfs_created = 0;
  while( surfaceElemList.size() > 0)
  {

    // start with the first face and create a list of all elements 
    // attached to the face

    DLIList<FacetEntity*> face_list;
    FacetEntity *start_face_ptr = surfaceElemList.get_and_step();
    stat = get_adj_facets( start_face_ptr, face_list, mydebug );
    if (stat != CUBIT_SUCCESS)
      return stat;

    // if we have the same number of faces on the face_list as we do
    // on the surfaceElemList, then we are done.  This surface is
    // not multiply connected.  Oherwise continue...

    if (face_list.size() == surfaceElemList.size())
    {

      // if this surface had a curve already defined (its a 2D topology
      // defined in skin2D) then its no longer valid if the surface was split
      // (for 3D the curves aren't defined until later)

      if (num_surfs_created > 0 && curveList.size() > 0)
      {
        ChollaCurve *fcm_ptr = curveList.get();  // there should only be 1
        curveList.remove();
        fcm_ptr->remove_td_associativity( this );
        delete fcm_ptr;
      }
      return CUBIT_SUCCESS;
    }

    // create a new surface to hold the face info

    num_surfs_created++;
    ChollaSurface *fsm_ptr = new ChollaSurface( blockId );
    facet_surface_list.append( fsm_ptr );

    // update the geometric curve pointer

    fsm_ptr->assign_geometric_surface( NULL );

    // add the faces to this surface and update the surface 
    // pointers in the face tool data

      // surfaceElemList: nullify the items then compress the list
      // afterwards, instead of removing the items one by one, for
      // speed.
    for (jj=face_list.size(); jj > 0; jj--)
    {
      FacetEntity *face_ptr = face_list.get_and_step();
      TDGeomFacet *td_gm_face = TDGeomFacet::get_geom_facet(face_ptr); 
      td_gm_face->remove_cholla_surf( this );
      td_gm_face->add_cholla_surf( fsm_ptr );
      surfaceElemList.move_to_nearby( face_ptr );
      surfaceElemList.extract();
      fsm_ptr->add_facet( face_ptr );
    }
  }

  return CUBIT_SUCCESS;
}


#if 0
//=============================================================================
//Function:  get_adj_facets (PRIVATE)
//Description: recursive funstion that creates a list of all faces connected
//             the passed in face  
//Author: sjowen
//Date: 12/22/00
//=============================================================================
CubitStatus ChollaSurface::get_adj_facets( 
  FacetEntity *start_face_ptr,
  DLIList<FacetEntity*> &face_list,
  int mydebug)
{
  CubitStatus stat = CUBIT_SUCCESS;

  // add this face to the list

  if (mydebug)
  {
    //start_face_ptr->draw( CUBIT_RED );
    //CDrawingTool::instance()->flush();
  }
  face_list.append( start_face_ptr );
  TDGeomFacet *td_gm_face = TDGeomFacet::get_geom_facet(start_face_ptr);
  td_gm_face->set_hit_flag( 0 );

  // loop through its edges

  CubitFacetEdge *edge_ptr;
  FacetEntity *face_ptr;
  DLIList<CubitFacetEdge*> edge_list;
  start_face_ptr->edges( edge_list );
  int ii;
  for (ii=0; ii<edge_list.size(); ii++)
  {
    edge_ptr = edge_list.get_and_step();
    
    // edges that already have a tool data defined are the result
    // of a feature angle.  Don't traverse past a feature angle edge
    
    TDGeomFacet *td_gm_edge = TDGeomFacet::get_geom_facet(edge_ptr);
    if (td_gm_edge == NULL)
    {
      DLIList<FacetEntity*> adj_face_list;
      edge_ptr->get_parents( adj_face_list );

      // keep traversing only if there are two adjacent faces to this edge,
      // otherwise, this is a boundary

      if (adj_face_list.size() == 2)
      {
        face_ptr = adj_face_list.get_and_step();
        if (face_ptr == start_face_ptr)
          face_ptr = adj_face_list.get();

        // go to its neighbor if it is part of the surface

        td_gm_face = TDGeomFacet::get_geom_facet(face_ptr); 
        if (td_gm_face->get_hit_flag() == id)
        {
          stat = get_adj_facets( face_ptr, face_list, mydebug );
          if (stat != CUBIT_SUCCESS)
            return stat;
        }
      }
    }
  }
  return stat;
}
#endif

//=============================================================================
//Function:  get_adj_facets (PRIVATE)
//Description: non recursive function that creates a list of all facets connected
//             to the passed in facet  
//Author: sjowen
//Date: 12/22/00
//=============================================================================
CubitStatus ChollaSurface::get_adj_facets( 
  FacetEntity *start_face_ptr,
  DLIList<FacetEntity*> &face_list,
  int mydebug,
  bool bound_check,
  bool feature_edge_check)
{
  //int found = 0;
  int ii;
  CubitStatus stat = CUBIT_SUCCESS;
  DLIList<FacetEntity*> temp_list;
  FacetEntity *face_ptr = NULL;
  FacetEntity *adj_face_ptr = NULL;
  DLIList<CubitFacetEdge *>edge_list;
  CubitFacetEdge *edge_ptr = NULL;
  DLIList<FacetEntity *>adj_face_list;

  if (mydebug)
  {
    for(ii=0; ii<surfaceElemList.size(); ii++)
    {
      face_ptr = surfaceElemList.get_and_step();
      TDGeomFacet *td_gm_face = TDGeomFacet::get_geom_facet(face_ptr);
      PRINT_INFO("%d ", td_gm_face->get_hit_flag());
      if (ii%10 == 0)
      {
        PRINT_INFO("\n");
      }
    }
  }

  // add this face to the list

  temp_list.append( start_face_ptr );
  TDGeomFacet *td_gm_face = TDGeomFacet::get_geom_facet(start_face_ptr);
  td_gm_face->set_hit_flag( 0 );

  while (temp_list.size())
  {
    face_ptr = temp_list.pop();
    td_gm_face = TDGeomFacet::get_geom_facet(face_ptr); 
    if (td_gm_face->get_hit_flag() == 0)
    {
      face_list.append( face_ptr );
      if (mydebug)
      {
        face_ptr->debug_draw( CUBIT_RED );
        GfxDebug::flush();
      }
      edge_list.clean_out();
      face_ptr->edges( edge_list );
      for (ii=0; ii<edge_list.size(); ii++)
      {
        edge_ptr = edge_list.get_and_step();

        // edges that already have a tool data defined are the result
        // of a feature angle.  Don't traverse past a feature angle edge
    
        TDGeomFacet *td_gm_edge = TDGeomFacet::get_geom_facet(edge_ptr);
        if (td_gm_edge == NULL || !feature_edge_check )
        {
          adj_face_list.clean_out();
          edge_ptr->get_parents( adj_face_list );

          // keep traversing only if there are two adjacent faces to this edge,
          // otherwise, this is a boundary

          if (adj_face_list.size() != 2)
          {
            continue;
          }

          if( bound_check )
          {
            TDFacetBoundaryEdge *td_facet_bnd_edge = TDFacetBoundaryEdge::get_facet_boundary_edge( edge_ptr );
            if( td_facet_bnd_edge )
              continue;
          }
            
          adj_face_ptr = adj_face_list.get_and_step();
          if (adj_face_ptr == face_ptr)
            adj_face_ptr = adj_face_list.get();

          // go to its neighbor if it is part of the surface

          td_gm_face = TDGeomFacet::get_geom_facet(adj_face_ptr); 
          if (td_gm_face->get_hit_flag() == id)
          {
            temp_list.append( adj_face_ptr );
            td_gm_face->set_hit_flag( 0 );
          }
          
        }
      }
    }
  }
  return stat;
}



//=============================================================================
//Function:  feature_angle (PRIVATE)
//Description: mark all edges that exceed the specified feature angle
//             min_dot is the minimum dot product between adjacent face normals
//Author: sjowen
//Date: 12/22/00
//=============================================================================
CubitStatus ChollaSurface::feature_angle( 
  double min_dot,
  DLIList<CubitFacetEdge *> &feature_edge_list)
{
  //CubitStatus stat = CUBIT_SUCCESS;
  int ii, jj;

  // compute face normals
  int mydebug = 0;
  double dot;
  CubitVector face_normal, adj_face_normal;
  FacetEntity *face_ptr, *adj_face_ptr;
  TDGeomFacet *td_gm_face;
  CubitFacet *tri_ptr;
  for (ii=0; ii<surfaceElemList.size(); ii++)
  {
    face_ptr = surfaceElemList.get_and_step();
    td_gm_face = TDGeomFacet::get_geom_facet(face_ptr);
    tri_ptr = CAST_TO( face_ptr, CubitFacet );
    face_normal = tri_ptr->normal( );
    face_normal.normalize();
    td_gm_face->set_normal( face_normal );
  }

  // check adjacencies and compute the dot product between them
  // where dot product is less than the min_dot, create a tool data
  // on the edge
  
  if(mydebug)
    GfxDebug::clear();
  for (ii=0; ii<surfaceElemList.size(); ii++)
  {
    face_ptr = surfaceElemList.get_and_step();
    CubitFacet* curr_facet = CAST_TO(face_ptr, CubitFacet);
    CubitFacet* temp_facet = NULL;
    double curr_area=curr_facet->area();
    td_gm_face = TDGeomFacet::get_geom_facet(face_ptr);
    face_normal = td_gm_face->get_normal();
    DLIList<CubitFacetEdge*> edge_list;
    face_ptr->edges( edge_list );
    for (jj=0; jj<edge_list.size(); jj++)
    {
      CubitFacetEdge *edge_ptr = edge_list.get_and_step();
      TDGeomFacet *td_gm_edge = TDGeomFacet::get_geom_facet(edge_ptr);
      if (!td_gm_edge)
      {
        DLIList<FacetEntity*> adj_face_list;
        edge_ptr->get_parents( adj_face_list );

        // it has to be an internal edge - ignore boundaries
        if (adj_face_list.size() == 2)
        {
          adj_face_ptr = adj_face_list.get_and_step();
          
          if (adj_face_ptr == face_ptr)
            adj_face_ptr = adj_face_list.get();          
          td_gm_face = TDGeomFacet::get_geom_facet(adj_face_ptr);

          // make sure the adj face is on the same surface

          if (td_gm_face->get_hit_flag() == id)
          {

            // test the dot product between normals

            adj_face_normal = td_gm_face->get_normal();
            temp_facet = CAST_TO(adj_face_ptr, CubitFacet);
            double adj_area=temp_facet->area();
              //mbrewer:: ensure no NULL
            dot = adj_face_normal % face_normal;
            bool add_an_edge = true;
            if(mydebug){
              PRINT_INFO("This area %f, other %f\n",curr_area,adj_area);
            }
            if(curr_area<(1.e-10*adj_area)){
              DLIList<CubitFacetEdge*> edge_list_tmp;
              CubitFacetEdge* edge_ptr_tmp=NULL;
              curr_facet->edges(edge_list_tmp);
              double edge_length = edge_ptr->length();
              int k = 0;
              if(edge_list_tmp.size()>0)
                add_an_edge=false;
              for(k=edge_list_tmp.size();k>0;k--){
                edge_ptr_tmp=edge_list_tmp.get_and_step();
                if(edge_ptr_tmp != edge_ptr &&
                   edge_ptr_tmp->length() > edge_length){
                  add_an_edge=true;
                }
              }
            }
            if(adj_area<(1.e-10*curr_area))
            {
              DLIList<CubitFacetEdge*> edge_list_tmp;
              CubitFacetEdge* edge_ptr_tmp=NULL;
              temp_facet->edges(edge_list_tmp);
              double edge_length = edge_ptr->length();
              int k = 0;
              if(edge_list_tmp.size()>0)
                add_an_edge=false;
              for(k=edge_list_tmp.size();k>0;k--){
                edge_ptr_tmp=edge_list_tmp.get_and_step();
                if(edge_ptr_tmp != edge_ptr &&
                   edge_ptr_tmp->length() > edge_length){
                  add_an_edge=true;
                }
              }
            }
            if (dot <= min_dot && add_an_edge )
            {
              if(mydebug){
                edge_ptr->debug_draw(CUBIT_MAGENTA);
              }
              TDGeomFacet::add_geom_facet(edge_ptr, -1); 
              td_gm_edge = TDGeomFacet::get_geom_facet(edge_ptr);
              edge_ptr->set_as_feature();
              feature_edge_list.append( edge_ptr );
            }
          }
        }

        // non-manifold edges (edges with more than 2 adj facets)
        // must be features

        else if (adj_face_list.size() > 2)
        {
          TDGeomFacet::add_geom_facet(edge_ptr, -1); 
          td_gm_edge = TDGeomFacet::get_geom_facet(edge_ptr);
          edge_ptr->set_as_feature();
          feature_edge_list.append( edge_ptr );
        }
      }
    }
  }
  if(mydebug){
    GfxDebug::mouse_xforms();
  }
  return CUBIT_SUCCESS;
}

//=============================================================================
//Function:  add_preexisting_feature_edges (PRIVATE)
//Description: edges that were marked previously in function ChollaEngine::mark_features
//             are added to the feature edge list
//Author: sjowen
//Date: 01/08
//=============================================================================
CubitStatus ChollaSurface::add_preexisting_feature_edges( 
	DLIList<CubitFacetEdge *> &feature_edge_list)
{
	DLIList<CubitFacetEdge *> edge_list;
	DLIList<CubitFacet *> facet_list;
	CAST_LIST( surfaceElemList, facet_list, CubitFacet );
	FacetDataUtil::get_edges( facet_list, edge_list );
	int iedge;
	CubitFacetEdge *edge;
	for (iedge=0; iedge < edge_list.size(); iedge++)
	{
		edge = edge_list.get_and_step();
		if (edge->is_feature())
			feature_edge_list.append(edge);
	}
	return CUBIT_SUCCESS;
}
	

//=============================================================================
//Function:  non_manifold_edges (PRIVATE)
//Description: mark all edges that are non-manifold (have more than 2 adj
//             facets)
//Author: sjowen
//Date: 5/01
//=============================================================================
CubitStatus ChollaSurface::non_manifold_edges( 
  DLIList<CubitFacetEdge *> &feature_edge_list)
{
  //CubitStatus stat = CUBIT_SUCCESS;
  int ii, jj;
  DLIList<CubitFacetEdge*> edge_list;
  FacetEntity *face_ptr;
  TDGeomFacet *td_gm_face;
  for (ii=0; ii<surfaceElemList.size(); ii++)
  {
    face_ptr = surfaceElemList.get_and_step();
    td_gm_face = TDGeomFacet::get_geom_facet(face_ptr);
    edge_list.clean_out();
    face_ptr->edges( edge_list );
    for (jj=0; jj<edge_list.size(); jj++)
    {
      CubitFacetEdge *edge_ptr = edge_list.get_and_step();
      TDGeomFacet *td_gm_edge = TDGeomFacet::get_geom_facet(edge_ptr);
      if (!td_gm_edge)
      {
        DLIList<FacetEntity*> adj_face_list;
        edge_ptr->get_parents( adj_face_list );

        // non-manifold edges (edges with more than 2 adj facets)
        // must be features

        if (adj_face_list.size() > 2 || edge_ptr->is_feature())
        {
          TDGeomFacet::add_geom_facet(edge_ptr, -1); 
          td_gm_edge = TDGeomFacet::get_geom_facet(edge_ptr);
          edge_ptr->set_as_feature();
          feature_edge_list.append( edge_ptr );
        }
      }
    }
  }
  return CUBIT_SUCCESS;
}



//=============================================================================
//Function:  clean_features (PRIVATE)
//Description: clean up edges that do not form complete loops as a result 
//             of feature angle
//Author: sjowen
//Date: 12/22/00
//=============================================================================
CubitStatus ChollaSurface::clean_features( )
{
  int ii, jj;
  for (ii=0; ii<surfaceElemList.size(); ii++)
  {
    FacetEntity *face_ptr = surfaceElemList.get_and_step();
    DLIList<CubitFacetEdge*> edge_list;
    face_ptr->edges( edge_list );
    for (jj=0; jj<edge_list.size(); jj++)
    {
      CubitFacetEdge *edge_ptr = edge_list.get_and_step();
      TDGeomFacet *td_gm_edge = TDGeomFacet::get_geom_facet(edge_ptr);
      if (td_gm_edge != NULL)
      {

        // make sure this edge actually forms the boundary between two surfaces

        DLIList<FacetEntity*> adj_face_list;
        edge_ptr->get_parents( adj_face_list );

        // it has to be an internal edge - ignore boundaries

        if (adj_face_list.size() == 2)
        {
          FacetEntity *adj_face_ptr = adj_face_list.get_and_step();
          if (adj_face_ptr == face_ptr)
            adj_face_ptr = adj_face_list.get();          
          TDGeomFacet *td_gm_adjface = TDGeomFacet::get_geom_facet(adj_face_ptr);
          if (td_gm_adjface != NULL)
          {
            DLIList<ChollaSurface*> surf_list;
            td_gm_adjface->get_cholla_surfs( surf_list );
            ChollaSurface *adj_surf = surf_list.get();  // must be only 1 surf
            if (adj_surf == this)
            {

              // the same surface is on both sides of this edge - Therefore it
              // does not form a complete loop.  Delete the tool data from
              // this edge to indicate that it is not at a boundary.

              edge_ptr->delete_TD( &TDGeomFacet::is_geom_facet );
            }
          }
        }
      }
    }
  }
  return CUBIT_SUCCESS;
}

//=============================================================================
//Function:  init_hit_flags (PUBLIC)
//Description: initialize the hit flags to the block surface id - used for 
//             traversing the faces
//Author: sjowen
//Date: 01/22/01
//=============================================================================
void ChollaSurface::init_hit_flags()
{
  int ii;
  for (ii=0; ii<surfaceElemList.size(); ii++)
  {
    FacetEntity *face_ptr = surfaceElemList.get_and_step();
    TDGeomFacet *td_gm = TDGeomFacet::get_geom_facet(face_ptr);
    td_gm->set_hit_flag( id );
  }
}

//=============================================================================
//Function:  get_points (PUBLIC)
//Description: compile a list of points on the surface 
//Author: sjowen
//Date: 4/01
//=============================================================================
void ChollaSurface::get_points( DLIList<CubitPoint *> &point_list )
{
  int ii;
  CubitPoint *point[3];
  CubitFacet *facet;
  FacetEntity *facet_entity;
  for (ii=0; ii<surfaceElemList.size(); ii++)
  {
    facet_entity = surfaceElemList.get_and_step();
    facet = CAST_TO( facet_entity, CubitFacet );
    facet->points( point[0], point[1], point[2] );
    point[0]->marked(0);
    point[1]->marked(0);
    point[2]->marked(0);
  }

  int jj;
  for (ii=0; ii<surfaceElemList.size(); ii++)
  {
    facet_entity = surfaceElemList.get_and_step();
    facet = CAST_TO( facet_entity, CubitFacet );
    facet->points( point[0], point[1], point[2] );
    for (jj=0; jj<3; jj++)
    {
      if (point[jj]->marked() == 0)
      {
        point[jj]->marked(1);
        point_list.append(point[jj]);
      }
    }
  }
}

//=============================================================================
//Function:  update_boundary_tool_data (PUBLIC)
//Description: update the surface IDs on the boundary facet tooldatas
//Author: sjowen
//Date: 5/01
//=============================================================================
CubitStatus ChollaSurface::update_boundary_tool_data()
{
  int ii, jj;
  CubitPoint *point_ptr;
  CubitFacetEdge *edge_ptr;
  CubitFacet *facet_ptr;
  FacetEntity *facet_entity;

  TDFacetBoundaryPoint *td_fbp;
  TDFacetBoundaryEdge *td_fbe;
  for (ii=0; ii<surfaceElemList.size(); ii++)
  {
    facet_entity = surfaceElemList.get_and_step();
    facet_ptr = CAST_TO( facet_entity, CubitFacet );
    for (jj=0; jj<3; jj++)
    {
      point_ptr = facet_ptr->point(jj);
      td_fbp = TDFacetBoundaryPoint::get_facet_boundary_point( point_ptr );
      if (td_fbp != NULL)
      {
        td_fbp->set_surf_id( facet_ptr, id );
      }

      edge_ptr = facet_ptr->edge(jj);
      td_fbe = TDFacetBoundaryEdge::get_facet_boundary_edge( edge_ptr );
      if (td_fbe != NULL)
      {
        td_fbe->set_surf_id( facet_ptr, id );
      }
    }
  }

  return CUBIT_SUCCESS;
}

//=============================================================================
//Function:  reset_facet_flags (PUBLIC)
//Description: 
//Author: sjowen
//Date: 6/01
//=============================================================================
void ChollaSurface::reset_facet_flags()
{
  CubitFacet *facet_ptr;
  FacetEntity *fe_ptr;
  int ii;
  for (ii=0; ii<surfaceElemList.size(); ii++)
  {
    fe_ptr = surfaceElemList.get_and_step();
    facet_ptr = CAST_TO( fe_ptr, CubitFacet );
    facet_ptr->marked( 0 );
  }
}

//=============================================================================
//Function:  is_adjacent (PUBLIC)
//Description: determine if this surface is adjacent to the given surface
//Author: sjowen
//Date: 02/11/2004
//=============================================================================
CubitBoolean ChollaSurface::is_adjacent( ChollaSurface *other_surf )
{
  DLIList<ChollaSurface *> *adjsurf_list_ptr;
  ChollaCurve *curv;
  ChollaSurface *surf;
  int icurv, ii;
  for(icurv=0; icurv< curveList.size(); icurv++)
  {
    curv = curveList.get_and_step();
    adjsurf_list_ptr = curv->get_surface_list_ptr();
    for(ii=0; ii<adjsurf_list_ptr->size(); ii++)
    {
      surf = adjsurf_list_ptr->get_and_step();
      if (surf != this && surf == other_surf)
      {
        return CUBIT_TRUE;
      }
    }
  }
  return CUBIT_FALSE; 
}

//=============================================================================
//Function:  get_loop_edges (PUBLIC)
//Description: return the ordered list of edges on the boundary of this surface 
//Author: sjowen
//Date: 02/23/2004
//=============================================================================
DLIList<DLIList<CubitFacetEdge *>*> *ChollaSurface::get_loop_edges(  )
{
  assert(myEvalTool != NULL);
  return myEvalTool->loops();
}


static int icolor = 0;
void ChollaSurface::debug_draw()
{
  icolor++;
  icolor = (icolor%15)+1;
  dcolor(icolor);
  dldraw(surfaceElemList);
}

//=============================================================================
//Function:  flip_facets (PUBLIC)
//Description: invert all facets on this surface 
//Author: sjowen
//Date: 09/10/09
//=============================================================================
void ChollaSurface::flip_facets()
{
  FacetEntity *facet_entity;
  CubitFacet *facet_ptr;
  for (int ii=0; ii<surfaceElemList.size(); ii++)
  {
    facet_entity = surfaceElemList.get_and_step();
    facet_ptr = dynamic_cast<CubitFacet *> (facet_entity);
    assert( facet_ptr != NULL );
    facet_ptr->flip();
  }
}

//=============================================================================
//Function:  is_in_volume (PUBLIC)
//Description: return whether this surface is in a particular volume 
//Author: sjowen
//Date: 09/11/09
//=============================================================================
CubitBoolean ChollaSurface::is_in_volume( ChollaVolume *chvol_ptr )
{
  ChollaVolume *mychvol_ptr;
  for (int ii=0; ii<volList.size(); ii++)
  {
    mychvol_ptr = volList.get_and_step();
    if (mychvol_ptr == chvol_ptr)
      return CUBIT_TRUE;
  }
  return CUBIT_FALSE;
}

//=============================================================================
//Function:  get_vertices (PUBLIC)
//Description: get the list of ChollaPoints on this surface
//Author: sjowen
//Date: 09/11/09
//=============================================================================
void ChollaSurface::get_vertices( DLIList<ChollaPoint *> &chpt_list )
{
  chpt_list.clean_out();
  ChollaCurve *chcurv_ptr;
  for (int ii=0; ii<curveList.size(); ii++)
  {
    chcurv_ptr = curveList.get_and_step();
    DLIList<ChollaPoint *> chc_pts = chcurv_ptr->get_points(); 
    chpt_list += chc_pts;
  }
  chpt_list.uniquify_unordered();
}


//EOF

