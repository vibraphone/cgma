//- Class:       ChollaSkinTool
//- Description: Creates a skin defined by facets and points, depending on the
//-              dimensionality of the desired skin.
//- Owner:       Steve Owen
//- Checked by:
//- Version:
#include "ChollaSkinTool.hpp"
#include "FacetEntity.hpp"
#include "TDGeomFacet.hpp"
#include "CastTo.hpp"
#include "ChollaSurface.hpp"
#include "ChollaCurve.hpp"
#include "CubitFacetEdge.hpp"
#include "TDFacetBoundaryEdge.hpp"
#include "debug.hpp"
#include "GfxDebug.hpp"

static int mydebug = 0;
//==================================================================================
//Function:  ChollaSkinTool (PUBLIC) (constructor)
//==================================================================================
ChollaSkinTool::ChollaSkinTool()
{
}
//==================================================================================
//Function:  ~ChollaSkinTool (PUBLIC) (destructor)
//==================================================================================
ChollaSkinTool::~ChollaSkinTool()
{
}

//==================================================================================
//Function:  skin_1d  (PUBLIC)
//Description:  creates a skin of the given facet entities.
//==================================================================================
CubitStatus ChollaSkinTool::skin_1d(DLIList<FacetEntity*> &facet_list,
                                   ChollaCurve *&facet_curve_mesh_ptr)
{
  CubitStatus rv = CUBIT_SUCCESS;
  
  // create a ChollaCurve if we have to (only if this is a 1D facet)
  
  if (!facet_curve_mesh_ptr)
  {
    FacetEntity *edge_ptr = facet_list.get();
    TDGeomFacet *td_gm = TDGeomFacet::get_geom_facet( edge_ptr );
    facet_curve_mesh_ptr = new ChollaCurve( td_gm->get_block_id() );
    
    // associate all of the tooldata on the faces of this surf with the 
    // new ChollaCurve
    
    int ii;
    for (ii=0; ii<facet_list.size(); ii++)
    {
      edge_ptr = facet_list.get_and_step();
      facet_curve_mesh_ptr->add_facet( edge_ptr );
      td_gm = TDGeomFacet::get_geom_facet( edge_ptr );
      td_gm->add_cholla_curve( facet_curve_mesh_ptr );
    }
  }
  
  // Note: the start and end points of this curve will be defined in 
  // ChollaCurve::split_curve.  The BlockPointMesh objects at these points
  // will be defined in MeshGeometryCreator::classify_node
  
  return rv;
}

//==================================================================================
//Function:  skin_2d  (PUBLIC)
//Description:  creates a skin of the given mesh entities.
//==================================================================================
CubitStatus ChollaSkinTool::skin_2d(DLIList<FacetEntity*> &facet_list,
                                   ChollaSurface *&facet_surface_mesh_ptr)
{
  if(facet_list.size() == 0){
    return CUBIT_SUCCESS;
  }
  
  int debugflag=0;
  if (debugflag)
  {
    dcolor(CUBIT_YELLOW);
    dldraw( facet_list );
    dview();
    dcolor(CUBIT_RED);
  }
  CubitStatus rv = CUBIT_SUCCESS;
  int block_id;

  // create a ChollaSurface if we have to (only if this is a 2D model)

  FacetEntity *face_ptr = facet_list.get();
  TDGeomFacet *td_gm = TDGeomFacet::get_geom_facet( face_ptr );
  block_id = td_gm->get_block_id();
  if (!facet_surface_mesh_ptr)
  {
    facet_surface_mesh_ptr = new ChollaSurface(block_id);

    // associate all of the tooldata on the faces of this surf with the new ChollaSurface

    int ii;
    for (ii=0; ii<facet_list.size(); ii++)
    {
      face_ptr = facet_list.get_and_step();
      facet_surface_mesh_ptr->add_facet( face_ptr );
      td_gm = TDGeomFacet::get_geom_facet( face_ptr );
      td_gm->add_cholla_surf( facet_surface_mesh_ptr );
      td_gm->set_hit_flag( facet_surface_mesh_ptr->get_id() );
    }
  }

  // create a single ChollaCurve for this surface (assumes one loop of edges)

  ChollaCurve *fcm_ptr = new ChollaCurve( block_id );
  facet_surface_mesh_ptr->add_curve( fcm_ptr );
  fcm_ptr->add_surface( facet_surface_mesh_ptr );

  // loop through all faces on this surface searching for the boundary edges

  int jj, kk, ll;
  for ( kk = 0; kk < facet_list.size(); kk++)
  {
    face_ptr = facet_list.get_and_step();
    DLIList<CubitFacetEdge*> edge_list;
    face_ptr->edges( edge_list );

    // loop through each edge on this face searching for boundary edges

    for (jj=edge_list.size(); jj > 0; jj--)
    {
      CubitFacetEdge *edge_ptr = edge_list.get_and_step();

      // check if this edge has already been processed from an adjacent surface.
      // If it has, then tool data would have already been defined at this edge
      // and by definition would be at the boundary (only tooldatas on edges
      // at the boundary ofa surface are created)

      TDGeomFacet *td_gm_edge = TDGeomFacet::get_geom_facet(edge_ptr);
      int on_boundary = 0;
      if (td_gm_edge != NULL)
      { 
        on_boundary = 1;

        // check for internal C-zero edge

        if (edge_ptr->num_adj_facets() == 2)
        {
          TDFacetBoundaryEdge *td_fbe = 
            TDFacetBoundaryEdge::get_facet_boundary_edge( edge_ptr );
          if (td_fbe != NULL)
          {
            if (td_fbe->is_internal_edge())
            {
              on_boundary = 0;
            }
          }
        } 
      }

      // check for general case where no tool data yet defined

      else 
      {
        DLIList<FacetEntity*> adj_face_list;

        // check the adjacent faces to this edge.  If only one adjacent face, then
        // it is on the boundary.  If more than one face, then the other face(s)
        // must be associated with a surface other than facet_surface_mesh_ptr
        // in order to be on the boundary

        edge_ptr->get_parents( adj_face_list );
        if (adj_face_list.size() <= 1)
        {
          on_boundary = 1;
        }
        else
        {
          for (ll=adj_face_list.size(); ll> 0 && !on_boundary; ll--)
          {
            FacetEntity *adj_face_ptr = adj_face_list.get_and_step();
            if (adj_face_ptr != face_ptr)
            {
              TDGeomFacet *td_gm_adjface = TDGeomFacet::get_geom_facet(adj_face_ptr);
              DLIList<ChollaSurface*> surf_list;
              td_gm_adjface->get_cholla_surfs( surf_list );

              // if it doesn't have an associated surface yet, then it is 
              // a neighboring surface that has not been defined yet (this
              // should only occur for the 2D case)
                       
              if (surf_list.size() == 0)
              {
                on_boundary = 1;
              }
              else
              {

                // there should only be one surface associated with
                // each face - otherwise we've screwed up somewhere  

                assert ( surf_list.size() == 1 );

                // if the surface is not the same as the current surface
                // then we are at the boundary

                ChollaSurface *check_bsm_ptr = surf_list.get();
                if (facet_surface_mesh_ptr != check_bsm_ptr)
                {
                  on_boundary = 1;
                }
              }
            }
          }
        }
      }
      if (on_boundary)
      {
        // create a tool data if needed

        if (td_gm_edge == NULL)
        {
          TDGeomFacet::add_geom_facet(edge_ptr, -1); 
          td_gm_edge = TDGeomFacet::get_geom_facet( edge_ptr );
          edge_ptr->set_as_feature();
        }

        // add the pointer to this surface onto the edge tool data

        td_gm_edge->add_cholla_surf( facet_surface_mesh_ptr );

        // add this edge to the curve

        fcm_ptr->add_facet( edge_ptr );
        if (mydebug)
          dedraw(edge_ptr);
      }
    }
  }

  return rv;
}

       
              
  

      
