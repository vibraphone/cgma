//-------------------------------------------------------------------------
// Filename		 : SyncTool.cpp
//
// Purpose		 : A class to facilitate synchronizing a solid modeler with
//               Cubit.	   
//
// Special Notes :
//
// Creator       : Derek Quam
//
// Creation Date : 07/03/2007
//
// Owner         : ???
//-------------------------------------------------------------------------
#include "SyncTool.hpp"
#include "Model.hpp"
#include "CompositeCurve.hpp"
#include "castto.hpp"
#include "CompositeSurface.hpp"
#include "RefFace.hpp"
#include "Gmem.hpp"
#include "ProeSurface.hpp"

//================================================================================
// Description: Notify the CGM that the geometry of the surfaces has been modified.
// Author     : Derek Quam
// Date       : 07/06/07
//================================================================================
CubitStatus SyncTool::notify_updated_surfaces()
{
  DLIList<RefEntity *> surf_list;
  Model::instance()->ref_entity_list("surface",surf_list, CUBIT_FALSE);
  for (int i = 0; i < surf_list.size(); i++)
  {
    surf_list.get()->notify_all_observers(GEOMETRY_MODIFIED);
    RefFace* curr_refface = CAST_TO(surf_list.get_and_step(),RefFace);
    if( curr_refface ){
      Surface* curr_surf = curr_refface->get_surface_ptr();
      CompositeSurface* curr_comp_surf = CAST_TO( curr_surf,CompositeSurface);
      if( curr_comp_surf ){
        int num_surfs = curr_comp_surf->num_surfs();
        curr_comp_surf->update();
        //GMem gmem;
        for (int j = 0; j < num_surfs; j++)
          {
          curr_comp_surf->notify_topology_modified( CAST_TO( curr_comp_surf->get_surface(j), TopologyBridge ) );
          ProeSurface* curr_proe_surf = CAST_TO( curr_comp_surf->get_surface(j),ProeSurface);
          if( curr_proe_surf )
            curr_proe_surf->tessellate_surface();
          }
        //curr_comp_surf->get_graphics(gmem);
        curr_comp_surf->recache_graphics();
        }
      }
  }
  return CUBIT_SUCCESS;
}

//================================================================================
// Description: Notify the CGM that the geometry of the curves has been modified.
// Author     : Derek Quam
// Date       : 07/06/07
//================================================================================
CubitStatus SyncTool::notify_updated_curves()
{
  DLIList<RefEntity *> curve_list;
  Model::instance()->ref_entity_list("curve", curve_list, CUBIT_FALSE);
  for (int i = 0; i < curve_list.size(); i++)
  {
    curve_list.get()->notify_all_observers(GEOMETRY_MODIFIED);
    RefEdge* curr_refedge = CAST_TO(curve_list.get_and_step(),RefEdge);
    if( curr_refedge ){
       Curve* curr_curve = curr_refedge->get_curve_ptr();
       CompositeCurve* curr_comp_curve = CAST_TO(curr_curve, CompositeCurve);
       if( curr_comp_curve ){
         curr_comp_curve->update();
         }
      }
  }
  return CUBIT_SUCCESS;
}
//================================================================================
// Description: Update any Virtual Geometry to reflect any updates.
// Author     : Alex Hays
// Date       : 07/19/07
//================================================================================
CubitStatus SyncTool::notify_virtual_curves()
{
  //Indiscriminantly update all Virtual Geometry contained in the Cubit model

  //Tell Virtual (composite) Curves to update cached length data
  DLIList<RefEntity *> curve_list;
  Model::instance()->ref_entity_list("curve", curve_list, CUBIT_FALSE);
  for (int i = 0; i < curve_list.size(); i++)
  {
    RefEdge* curr_refedge = CAST_TO(curve_list.get_and_step(),RefEdge);
    if( curr_refedge ){
       Curve* curr_curve = curr_refedge->get_curve_ptr();
       CompositeCurve* curr_comp_curve = CAST_TO(curr_curve, CompositeCurve);
       if( curr_comp_curve ){
         curr_comp_curve->update();
         }
      }
  }
  
  return CUBIT_SUCCESS;
}
//================================================================================
// Description: Update any Virtual Geometry to reflect any updates.
// Author     : Alex Hays
// Date       : 07/19/07
//================================================================================
CubitStatus SyncTool::notify_virtual_surfaces()
{
  //Indiscriminantly update all Virtual Geometry contained in the Cubit model

  //Tell Virtual (composite) Surfaces to update all their contained surfaces
  DLIList<RefEntity *> surf_list;
  Model::instance()->ref_entity_list("surface",surf_list, CUBIT_FALSE);
  for (int i = 0; i < surf_list.size(); i++)
  {
    RefFace* curr_refface = CAST_TO(surf_list.get_and_step(),RefFace);
    if( curr_refface ){
      Surface* curr_surf = curr_refface->get_surface_ptr();
      CompositeSurface* curr_comp_surf = CAST_TO( curr_surf,CompositeSurface);
      if( curr_comp_surf ){
        int num_surfs = curr_comp_surf->num_surfs();
        for (int j = 0; j < num_surfs; j++)
          {
          curr_comp_surf->notify_topology_modified( CAST_TO( curr_comp_surf->get_surface(j), TopologyBridge ) );

          }
        }
      }
  }
  return CUBIT_SUCCESS;
}