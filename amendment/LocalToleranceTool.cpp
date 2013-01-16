//-------------------------------------------------------------------------
// Filename      : LocalToleranceTool.cpp
//
// Purpose       : The tool used to calculate local tolerance at RefEntities (Vert, Curve, Surface, Volume)
//
// Special Notes :
//
// Creator       : William Roshan Quadros
//
// Creation Date : 11/16/2010
//
//-------------------------------------------------------------------------

#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <sstream>
#include <iostream>

#include "LocalToleranceTool.hpp"

#include "CubitDefines.h"
#include "GeometryDefines.h"
#include "GeometryQueryTool.hpp"
#include "MergeTool.hpp"

#include "RefEntityFactory.hpp"
#include "BodySM.hpp"
#include "BasicTopologyEntity.hpp"

#include "RefEntity.hpp"
#include "RefVertex.hpp"
#include "RefEdge.hpp"
#include "RefFace.hpp"
#include "RefVolume.hpp"


// Instance of static members
LocalToleranceTool* LocalToleranceTool::instance_ = 0;

//===================================================================================
// Description: constructor
// Notes:  
// Author: roshan
// Date: 11/16/2010
//===================================================================================
LocalToleranceTool::LocalToleranceTool()
{

}

//===================================================================================
// Description: destructor
// Notes:  
// Author: roshan
// Date: 11/16/2010
//===================================================================================
LocalToleranceTool::~LocalToleranceTool()
{

}

//-------------------------------------------------------------------------
// Purpose       : Controls access and creation of the sole instance of this
//                 class.
//
// Special Notes :
//
// Creator       : William Roshan Quadros
//
// Creation Date : 11/16/2010
//-------------------------------------------------------------------------

LocalToleranceTool* LocalToleranceTool::instance(void)
{
   // Check to see if we have created an instance of the class
   // If not, proceed to create one.

  //if (LocalToleranceTool::instance_ == 0)
  // {
  //   LocalToleranceTool::instance_ = new LocalToleranceTool() ;

  //    
  // }

   // Return the a pointer to the instance of the class.
  return NULL; //LocalToleranceTool::instance_ ;
}

//-------------------------------------------------------------------------
// Purpose       : Deletes instance variable
//
// Special Notes :
//
// Creator       : William Roshan Quadros
//
// Creation Date : 11/16/2010
//-------------------------------------------------------------------------
void LocalToleranceTool::delete_instance()
{
  if( NULL != instance_ )
  {
    delete instance_;
    instance_ = NULL;
  }
}

//-------------------------------------------------------------------------
// Purpose       : calculates the local tolerances at all ref_entities of input bodies
//
// Special Notes :
//
// Creator       : William Roshan Quadros
//
// Creation Date : 11/16/2010
//-------------------------------------------------------------------------
bool LocalToleranceTool::calculate_local_tolerances( DLIList<BodySM*> body_sm_list )
{
  bool rv = true;
  const double default_vert_local_tol = 0.001;
  const double default_edge_local_tol = 0.01;
  const double default_face_local_tol = 0.02;
  const double default_vol_local_tol = 0.03;
  int i;

  // get all bodies
  DLIList<Body *> bodies;
  RefEntityFactory::instance()->bodies( bodies );

  // find all ref_vert in the bodies
  DLIList<RefVertex *> ref_verts;
  GeometryQueryTool::instance()->ref_vertices( ref_verts );

  // for now set default value in all ref_vert
  RefVertex *ref_vert;
  for( i = 0; i < ref_verts.size(); i++ )
  {
    ref_vert = ref_verts.get_and_step();
    ref_vert->local_tolerance( default_vert_local_tol );
  }

  // find all ref_edge in the bodies
  DLIList<RefEdge *> ref_edges;
  GeometryQueryTool::instance()->ref_edges( ref_edges );

  // for now set default value in all ref_edges
  RefEdge *ref_edge;
  for( i = 0; i < ref_edges.size(); i++ )
  {
    ref_edge = ref_edges.get_and_step();
    ref_edge->local_tolerance( default_edge_local_tol );
  }

  // find all ref_face in the bodies
  DLIList<RefFace *> ref_faces;
  GeometryQueryTool::instance()->ref_faces( ref_faces );

  // for now set default value in all ref_faces
  RefFace *ref_face;
  for( i = 0; i < ref_faces.size(); i++ )
  {
    ref_face = ref_faces.get_and_step();
    ref_face->local_tolerance( default_face_local_tol );
  }

  // find all ref_vol in the bodies
  DLIList<RefVolume *> ref_vols;
  GeometryQueryTool::instance()->ref_volumes( ref_vols );

  // for now set default value in all ref_vols
  RefVolume *ref_vol;
  for( i = 0; i < ref_vols.size(); i++ )
  {
    ref_vol = ref_vols.get_and_step();
    ref_vol->local_tolerance( default_vol_local_tol );
  }
  
  return rv;
}

//-------------------------------------------------------------------------
// Purpose       : calculates the local tolerances automatically at all the ref_entities of input bodies
//
// Special Notes :
//
// Creator       : William Roshan Quadros
//
// Creation Date : 11/16/2010
//-------------------------------------------------------------------------
bool LocalToleranceTool::calculate_local_tolerances_automatically( DLIList<BodySM*> body_sm_list )
{
  bool rv = true;

  // find mergable surfaces
  double geom_factor = GeometryQueryTool::get_geometry_factor();
  double merge_tolerance = geom_factor*GEOMETRY_RESABS;
  DLIList< DLIList<Surface*>*> lists_of_mergeable_surfaces;
  MergeTool::instance()->find_only_mergeable_surfaces ( body_sm_list, lists_of_mergeable_surfaces, merge_tolerance );


  /* Pseudocode for future use
  // Map holding overlapping surface pairs
  overlapping_surfs_map;
  // Add all surfs to spatial tree
  SpatialTree.add(all_surfs);
  // Create array of all surfs sorted by surface area (smallest to largest)
  all_surf_array.sort();
  // Loop through progressively larger merge tolerances
  for(cur_tol=.00001; cur_tol < .1; cur_tol *= 10.0)
  {
    // Each time array size may change because surfs have been removed.
    int num_surfs = all_surf_array.size(), i;
    for(i=0; i<num_surfs; i++)
    {
      cur_surf = all_surf_array[i];
      // Get surfs that are close to the current surf.
      close_surfs = SpatialTree.get_close_surfs(cur_surf);
      // Keep a running total of overlap area so we know when current surf
      // is overlapping completely with other surfs.
      double overlap_area = 0.0;
      int j;
      // Loop over all of the close surfs and look for overlaps.
      for(j=close_surfs.size(); j>0; j--)
      {
        cur_close_surf = close_surfs.get_and_step();
        // If there is already a map entry between these two surfs don't check them again.
        if(!overlapping_surfs_map.entry_exists(cur_surf, cur_close_surf))
        {
          // Measure overlap for these two surfs at current tolerance
          cur_overlap = measure_overlap(cur_surf, cur_close_surf, cur_tol);
          if(cur_overlap > 0.0)
          {
            // add entry for these two surfs (at the current tol)
            overlapping_surfs_map.add_entry(cur_surf, cur_close_surf, cur_tol);
            // Update the overlapping area for cur_surf
            overlap_area += cur_overlap;
            // If overlapping area for cur_surf is all of cur_surf then we are done processing it.
            if(overlap_area == cur_surf.surface_area())
            {
              // Remove cur_surf from all_surf_array
              // Remove cur_surf from SpatialTree (if this helps with efficiency)
              // jump out of processing of cur_surf
              j=0;
            }
          }
        }
      }
    }
  }
  */
  return rv;
}

//-------------------------------------------------------------------------
// Purpose       : debugging function to print local tolerances at all ref_entities of the input bodies
//
// Special Notes :
//
// Creator       : William Roshan Quadros
//
// Creation Date : 11/16/2010
//-------------------------------------------------------------------------
bool LocalToleranceTool::print_local_tolerances( DLIList<BodySM*> body_sm_list )
{
  bool rv = true;
  int i;

  // get all bodies
  DLIList<Body *> bodies;
  RefEntityFactory::instance()->bodies( bodies );

  // find all ref_vert in the bodies
  DLIList<RefVertex *> ref_verts;
  GeometryQueryTool::instance()->ref_vertices( ref_verts );

  // for now set default value in all ref_vert
  RefVertex *ref_vert;
  for( i = 0; i < ref_verts.size(); i++ )
  {
    ref_vert = ref_verts.get_and_step();
    PRINT_INFO( " vert id %d, local_tol %lf \n", ref_vert->id(), ref_vert->local_tolerance() );
  }

  // find all ref_edge in the bodies
  DLIList<RefEdge *> ref_edges;
  GeometryQueryTool::instance()->ref_edges( ref_edges );

  // for now set default value in all ref_edges
  RefEdge *ref_edge;
  for( i = 0; i < ref_edges.size(); i++ )
  {
    ref_edge = ref_edges.get_and_step();
    PRINT_INFO( " edge id %d, local_tol %lf \n", ref_edge->id(), ref_edge->local_tolerance() );
  }

  // find all ref_face in the bodies
  DLIList<RefFace *> ref_faces;
  GeometryQueryTool::instance()->ref_faces( ref_faces );

  // for now set default value in all ref_faces
  RefFace *ref_face;
  for( i = 0; i < ref_faces.size(); i++ )
  {
    ref_face = ref_faces.get_and_step();
    PRINT_INFO( " face id %d, local_tol %lf \n", ref_face->id(), ref_face->local_tolerance() );
  }

  // find all ref_vol in the bodies
  DLIList<RefVolume *> ref_vols;
  GeometryQueryTool::instance()->ref_volumes( ref_vols );

  // for now set default value in all ref_vols
  RefVolume *ref_vol;
  for( i = 0; i < ref_vols.size(); i++ )
  {
    ref_vol = ref_vols.get_and_step();
    PRINT_INFO( " vol id %d, local_tol %lf \n", ref_vol->id(), ref_vol->local_tolerance() );
  }
  
  return rv;
}

