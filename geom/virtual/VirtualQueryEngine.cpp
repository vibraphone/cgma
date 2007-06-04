//-------------------------------------------------------------------------
// Filename      : VirtualQueryEngine.cc 
//
// Purpose       : This file contains the implementation of the class 
//                 VirtualQueryEngine. 
//
// Special Notes : There are several functions in this class which are 
//                 inhereted from GeometryQueryEngine that are not
//                 implemented or applicable in this class.  These
//                 functions, if called, will print an error message
//                 and return the most appropriate error condition available.
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 08/01/97 
//-------------------------------------------------------------------------

// ********** BEGIN STANDARD INCLUDES      **********
#include <assert.h>
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN CUBIT INCLUDES         **********

#include "VirtualQueryEngine.hpp"
#include "CompositeTool.hpp"
#include "PartitionTool.hpp"
#include "CompositeEngine.hpp"
#include "PartitionEngine.hpp"
#include "GeometryUtil.hpp"

#include "CubitDefines.h"
#include "GeometryDefines.h"
#include "CubitVector.hpp"

#include "CastTo.hpp"
#include "GMem.hpp"

#include "ModelQueryEngine.hpp"
#include "ModelEntity.hpp"

#include "Body.hpp"
#include "CoVolume.hpp"
#include "RefVolume.hpp"
#include "RefFace.hpp"
#include "Loop.hpp"
#include "CoEdge.hpp"
#include "RefEdge.hpp"
#include "Chain.hpp"
#include "CoVertex.hpp"
#include "RefVertex.hpp"
#include "RefEntityFactory.hpp"
#include "CoFace.hpp"


#include "CompositePoint.hpp"
#include "CompositeCurve.hpp"
#include "CompositeSurface.hpp"
#include "CompositeLump.hpp"
#include "CompositeBody.hpp"
#include "PartitionBody.hpp"
#include "PartitionLump.hpp"
#include "PartitionSurface.hpp"
#include "PartitionCurve.hpp"


#include "Shell.hpp"
#include "GeometryQueryTool.hpp"
#include "CubitPointData.hpp"
#include "CubitFacetData.hpp"

#include "CACompositeVG.hpp"
#include "CAPartitionVG.hpp"
#include "CAVirtualVG.hpp"

// ********** END CUBIT INCLUDES           **********

// ********** BEGIN STATIC DECLARATIONS    **********
VirtualQueryEngine* VirtualQueryEngine::instance_ = NULL; 
//- A pointer to the single instance of this class


const int VirtualQueryEngine::VGE_MAJOR_VERSION = 10;
const int VirtualQueryEngine::VGE_MINOR_VERSION = 0;
const int VirtualQueryEngine::VGE_SUBMINOR_VERSION = 0;

const double VirtualQueryEngine::param_epsilon_fraction = 1000000;
//- This value is used to determine epsilon for
//- compensating for rounding error in parameter
//- conversion between virtaul geometry and real
//- geometry.

// ********** END STATIC DECLARATIONS      **********

// ********** BEGIN PUBLIC FUNCTIONS       **********

//-------------------------------------------------------------------------
// Purpose	 : Destructor
//
// Special Notes :
//
// Creator	 : Jason Kraftcheck
//
// Creation Date : 08/01/97
//-------------------------------------------------------------------------
VirtualQueryEngine::~VirtualQueryEngine( )
{
//  destroy_all_virtual_entities();
//	if( virtual_entity_list_.size() )
//		PRINT_WARNING("Not all virtual geometry has been removed.  "
//		              "This is a bug (memory leak).  Please report it.\n");
//	assert( virtual_entity_list_.size() == 0 );
  instance_ = NULL;
}


//-------------------------------------------------------------------------
// Purpose	 : Return a description of this ModelingEngine
//
// Special Notes :
//
// Creator	 : Jason Kraftcheck
//
// Creation Date : 08/01/97
//-------------------------------------------------------------------------
int VirtualQueryEngine::get_major_version()
{
  return VGE_MAJOR_VERSION;
}

int VirtualQueryEngine::get_minor_version()
{
  return VGE_MINOR_VERSION;
}

int VirtualQueryEngine::get_subminor_version()
{
  return VGE_SUBMINOR_VERSION;
}

CubitString VirtualQueryEngine::get_engine_version_string()
{
  CubitString return_this(
    "VirtualQueryEngine v1.2 by Jason Kraftcheck.\n" );
  return_this += "I-CARVE Lab., University of Wisconsin - Madison,  August 1, 1997.\n";
  return_this += "This ModelingEngine performs the necessary topological ";
  return_this += "operations for working with VirtualEntities, including ";
  return_this += "objects of type CompositeEntity, PartitionEntity, and ";
  return_this += "ParasiteEntity.  ComposinteModelingEngine also ";
  return_this += "provides the interface for creating, destroying, and ";
  return_this += "modifying VirtualEntities.\n";
  
  return return_this;
}

CubitStatus VirtualQueryEngine::transform_vec_position(
  CubitVector const& ,
  BodySM *,
  CubitVector & ) const
{
  return CUBIT_FAILURE;
}
//-------------------------------------------------------------------------
// Purpose       : This function determines faceting information for use
//		             by HOOPS.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 
//-------------------------------------------------------------------------
CubitStatus VirtualQueryEngine::get_graphics(
  Surface* surf_ptr, int& number_triangles,
  int& number_points, int& number_facets,
  GMem* gMem,
  unsigned short normal_tol,
  double dist_tol,
  double max_edge_length ) const
{
  CompositeSurface* cs_ptr = CAST_TO(surf_ptr, CompositeSurface );
  PartitionSurface* ps_ptr = CAST_TO(surf_ptr, PartitionSurface );
  if( cs_ptr != NULL )
  {
    return get_composite_surface_facetting( cs_ptr,
                                            number_triangles, number_points,
                                            number_facets, gMem,
                                            normal_tol, dist_tol, max_edge_length );
  }
  else if( ps_ptr != NULL )
  {
    return get_partition_surface_facetting( ps_ptr,
                                            number_triangles, number_points,
                                            number_facets, gMem,
                                            normal_tol, dist_tol, max_edge_length );
  }
  else
  {
    PRINT_INFO("VirtualQueryEngine::get_graphics_facets"
               "( s%p, %d, %d, %d, %p, %d, %f ) = CUBIT_FAILURE\n",
               surf_ptr, number_triangles, 
               number_points, number_facets, gMem, normal_tol, dist_tol);
    return CUBIT_FAILURE;
  }
}

//-----------------------------------------------------------------------
// Purpose	 : This function returns the necessary faceting information
//		   to display the curve.
//
// Special Notes :
//
// Creator	 : Jason Kraftcheck
//
// Creation Date : 10/10/97
//-----------------------------------------------------------------------
CubitStatus VirtualQueryEngine::get_graphics( Curve* curve_ptr,
                                                 int& num_points,
                                                 GMem* gMem,
                                                 double ) const
{
  CompositeCurve* ccurve = CAST_TO( curve_ptr, CompositeCurve );
  PartitionCurve* pcurve = CAST_TO( curve_ptr, PartitionCurve );
  
  CubitStatus result;
  
  if( ccurve != NULL )
  {
    result = get_composite_curve_facetting( ccurve, num_points, gMem );
  }
  else if( pcurve != NULL )
  {
    result = get_partition_curve_facetting( pcurve, num_points, gMem );
  }
  else
  {
    PRINT_ERROR("Unknown curve type or NULL curve passed to "
                "VirtualQueryEngine::get_graphics_edges()\n");
    return CUBIT_FAILURE;
  }
  
  return result;
}

// This is the template for an unimplemented version of this function.
CubitStatus VirtualQueryEngine::get_isoparametric_points(Surface*,
                                                            int& nu, int& nv,
                                                            GMem*& ) const
{
  nu = nv = 0;
  return CUBIT_FAILURE;
}

// This is the template for an unimplemented version of this function.
CubitStatus VirtualQueryEngine::get_u_isoparametric_points(Surface*,
                                                            double , int& ,
                                                            GMem*& ) const
{
  return CUBIT_FAILURE;
}

// This is the template for an unimplemented version of this function.
CubitStatus VirtualQueryEngine::get_v_isoparametric_points(Surface*,
                                                            double , int&,
                                                            GMem*&) const
{
  return CUBIT_FAILURE;
}

//-----------------------------------------------------------------------
// Purpose	 : This function constructs the edge facetting information
//		   of a CompositeCurve for use by HOOPS.
//
// Special Notes :
//
// Creator	 : Jason Kraftcheck
//
// Creation Date : 10/10/97
//-----------------------------------------------------------------------
CubitStatus VirtualQueryEngine::get_composite_curve_facetting(
  CompositeCurve* ccurve_ptr,
  int& num_points,
  GMem* gMem ) const
{
    // Just return if gMem is NULL
  if (gMem == NULL)
  {
    num_points = 0;
    return CUBIT_SUCCESS;
  }
  
  GMem current_gmem; 
  int point_count = 0; //Current total of points returned
    // from underlying edges
  
  int next_index = 0; //index of next insert spot in passed GMem object.
  
    //Get points for each underlying curve
  current_gmem.clean_out();
  for( int i = 0; i < ccurve_ptr->num_curves(); i++ )
  {
      // Get the next curve and its GME
    Curve* curve_ptr = ccurve_ptr->get_curve(i);
    GeometryQueryEngine* GQE_ptr = 
      curve_ptr->get_geometry_query_engine();
    
      // Get the GME to facet the curve
    int this_point_count = 0;
    CubitStatus current_status = 
      GQE_ptr->get_graphics (curve_ptr, this_point_count, &current_gmem);
    if( current_status == CUBIT_FAILURE )
      return CUBIT_FAILURE;
      
    // Make sure the passed in gmem is big enough
    gMem->allocate_more_polylines(this_point_count-1);
      //assert( gMem->pts_size >= point_count );
    
      // If this isn't the first curve, we don't need to
      // repeat one of the end nodes.
    if (point_count != 0)
      this_point_count--;
      // Add to the total number of points so far
    point_count += this_point_count;
    
      //Copy points to passed GMem, in the correct order
      //for the sense of the underlying edge wrt the
      //composite.
    CubitSense sense = ccurve_ptr->get_sense( i );
    
    if( sense == CUBIT_FORWARD )
    {
      if (next_index == 0)
        for(int j = 0; j < point_count; j++)
        {
          gMem->point_list()[j] = current_gmem.point_list()[j];
        }
      else
          // cur_gmem has one more point than we will access,
          // so last index is point_count-next_index.
        for(int j = next_index; j < point_count; j++)
        {
          gMem->point_list()[j] = current_gmem.point_list()[j - next_index + 1];
        }
    }
    else //( sense == CUBIT_REVERSED )
    {
        // All but the first curve will not use one of its points.
        // This is taken into account by decrementing point_count
        // earlier.
      for(int j = next_index; j < point_count; j++)
      {
        gMem->point_list()[j] = current_gmem.point_list()[point_count-j-1];
      }
    }
    next_index = point_count;
    gMem->pointListCount = point_count;
  }
  
  num_points = point_count;
  gMem->pointListCount = point_count;
  return CUBIT_SUCCESS;
}

//-----------------------------------------------------------------------
// Purpose	 : This function constructs the edge facetting information
//		   of a PartitionCurve for use by HOOPS.
//
// Special Notes : It is currently not possible to split a face which
//                 has more than one loop (holes) because we cannot
//                 determine which partition the holes should lie on with
//                 the current implementation of ParasiteCurve.
//
// Creator       : Jason Kraftcheck 
//
// Creation Date : 10/10/97
//-----------------------------------------------------------------------
CubitStatus VirtualQueryEngine::get_partition_curve_facetting(
  PartitionCurve* pcurve_ptr, int& num_steps, GMem* gMem ) const
{
  assert( gMem != NULL );
  CubitStatus result = pcurve_ptr->get_graphics( *gMem );
  num_steps = gMem->pointListCount;
  return result;
}


//-------------------------------------------------------------------------
// Purpose       : get Hoops facetting information for a CompositeSurface.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/17/98
//-------------------------------------------------------------------------
CubitStatus VirtualQueryEngine::get_composite_surface_facetting(
  CompositeSurface* surf_ptr,
  int& number_triangles,
  int& number_points,
  int& number_facets,
  GMem* gMem,
  unsigned short normal_tol,
  double dist_tol,
  double longest_edge ) const
{
  if( gMem == NULL )
    return CUBIT_FAILURE;
    
  if (surf_ptr->get_graphics(*gMem))
  {
    number_triangles = gMem->fListCount / 4;
    number_points = gMem->pointListCount;
    number_facets = gMem->fListCount;
    return CUBIT_SUCCESS;
  }
    
  number_triangles = 0;
  number_points = 0;
  number_facets = 0;

    // Get the underlying Surfaces
  GeometryQueryEngine* gqe_ptr;
  GMem *gMem_list = new GMem[ surf_ptr->num_surfs() ];
  
    // Get the facets of each Surface
  int i;
  for(i = 0; i < surf_ptr->num_surfs(); i++)
  {
    int tri_count, pt_count, facet_count;
    Surface* surface_ptr = surf_ptr->get_surface(i);
    
    gqe_ptr = surface_ptr->get_geometry_query_engine();
    
    tri_count = pt_count = facet_count = 0;
    gqe_ptr->get_graphics( surface_ptr, tri_count, 
                                  pt_count, facet_count, &(gMem_list[i]),
                                  normal_tol, dist_tol, longest_edge );
    
    number_triangles += tri_count;
    number_facets += gMem_list[i].fListCount;
    number_points += gMem_list[i].pointListCount;
  }
  
  int point_offset = 0;
  int face_offset = 0;
  gMem->replace_point_list(new GPoint[number_points],
                           number_points, number_points);
  gMem->replace_facet_list(new int[number_facets],
                           number_facets, number_facets);
  
  for( i = 0; i < surf_ptr->num_surfs(); i++ )
  {
    int j;
    for(  j = 0; j < gMem_list[i].pointListCount; j++ )
    {
      assert((point_offset + j) < gMem->pointListCount );
      gMem->point_list()[j + point_offset].x = gMem_list[i].point_list()[j].x;
      gMem->point_list()[j + point_offset].y = gMem_list[i].point_list()[j].y;
      gMem->point_list()[j + point_offset].z = gMem_list[i].point_list()[j].z;
    }
    
    int k;
    for( k = 0; k < gMem_list[i].fListCount;)
    {
      int count = gMem_list[i].facet_list()[k];
      assert( k + count < gMem->fListCount );
      
      gMem->facet_list()[k + face_offset] = count;
      k++;
      for( int l = 0; l < count; l++, k++ )
      {
        assert( k + face_offset < gMem->fListCount );
        assert( k < gMem_list[i].fListCount );
        gMem->facet_list()[k + face_offset] = gMem_list[i].facet_list()[k] 
          + point_offset;
      }
    }
    
    point_offset += j;
    face_offset += k;
  }
  
  gMem->consolidate_points(10*GEOMETRY_RESABS);
  number_points = gMem->pointListCount;
//     // Debug!!!
//   int q;
//   for (q = 0; q < bte_list.size(); q++)
//   {
//     PRINT_INFO("Printing Points for Surface %d.\n", q);
//     int r;
//     for (r = 0; r < gMem_list[q].pointListCount; r++)
//     {
//       PRINT_INFO("%g %g %g\n",
//                  gMem_list[q].point_list()[r].x,
//                  gMem_list[q].point_list()[r].y,
//                  gMem_list[q].point_list()[r].z);
//     }
//     PRINT_INFO("Printing FList for Surface %d.\n", q);
//     for (r = 0; r < gMem_list[q].fListCount; r++)
//     {
//       PRINT_INFO("%d %s",
//                  gMem_list[q].facet_list()[r],
//                  r % 4 ? "" : "\n");
//     }
//   }
//   PRINT_INFO("Composite points:\n");
//   for (q = 0; q < gMem->pointListCount; q++)
//   {
//     PRINT_INFO("%g %g %g\n",
//                gMem->point_list()[q].x,
//                gMem->point_list()[q].y,
//                gMem->point_list()[q].z);
//   }
//   PRINT_INFO("Composite facet list:\n");
//   for (q = 0; q < gMem->fListCount; q++)
//   {
//     PRINT_INFO("%d %s",
//                gMem->facet_list()[q],
//                q % 4 ? "" : "\n");
//   }
  
  delete [] gMem_list;	
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : get hoops facets for a partition surface
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck 
//
// Creation Date : 06/18/98
//-------------------------------------------------------------------------
CubitStatus VirtualQueryEngine::get_partition_surface_facetting(
  PartitionSurface* surf_ptr,
  int& num_triangles,
  int& num_points,
  int& num_facets,
  GMem* gMem,
  unsigned short ,
  double,
  double ) const
{
  int i, j;
  
  DLIList<CubitFacetData*> surf_facets;
  DLIList<CubitPoint*> surf_points;
  
    // get surface facets
  surf_ptr->get_facet_data( surf_facets );

    // get list of points from list of facets
  for( i = surf_facets.size(); i--; )
  {
    CubitFacetData* facet = surf_facets.step_and_get();
    for( j = 0; j< 3; j++ )
      facet->point(j)->marked(1);
  }
  for( i = surf_facets.size(); i--; )
  {
    CubitFacetData* facet = surf_facets.step_and_get();
    for( j = 0; j< 3; j++ )
    {
      if( facet->point(j)->marked() )
      {
        facet->point(j)->marked(0);
        surf_points.append( facet->point(j) );
      }
    }
  }

  if( ! gMem ) return CUBIT_SUCCESS;

    // allocate storage for facet data
  gMem->allocate_tri( surf_facets.size() );
  gMem->fListCount = surf_facets.size() * 4;
  gMem->pointListCount = surf_points.size();
    
    // set return values and if gMem is NULL return
  num_triangles = surf_facets.size();
  num_facets = gMem->fListCount;
  num_points = gMem->pointListCount;
  
    // put points in GMem and mark each point with
    // its index in the array
  surf_points.reset();
  for( i = 0; i < gMem->pointListCount; i++ )
  {
    GPoint* gpt = gMem->point_list() + i;
    CubitPoint* pt = surf_points.get_and_step();
    pt->marked(i);
    gpt->x = (float)(pt->coordinates().x() );
    gpt->y = (float)(pt->coordinates().y() );
    gpt->z = (float)(pt->coordinates().z() );
  }
  
    // put each facet in the GMem.  retreive
    // point indices from marks.
  surf_facets.reset();
  int* gfacet_ptr = gMem->facet_list();
  for( i = 0; i < surf_facets.size(); i++ )
  {
    CubitFacetData* facet = surf_facets.get_and_step();
    *gfacet_ptr = 3;
    gfacet_ptr++;
    for( j = 0; j < 3; j++ )
    {
      *gfacet_ptr = facet->point(j)->marked();
      gfacet_ptr++;
    }
  }
  
    // clear point marks
  for( i = surf_points.size(); i--; )
    surf_points.step_and_get()->marked(0);
  
  return CUBIT_SUCCESS;
}

// ********** END PUBLIC FUNCTIONS         **********

// ********** BEGIN PROTECTED FUNCTIONS    **********

//-------------------------------------------------------------------------
// Purpose       : The constructor of the VirtualQueryEngine class.
//
// Special Notes : 
//
// Creator       : Wes Gill
//
// Creation Date : 10/5/01
//-------------------------------------------------------------------------

VirtualQueryEngine::VirtualQueryEngine()
{
  assert( !instance_);
    //Add the VirtualQueryEngine to the GeometryQueryTool's gqeList and
    //set the VGE to the GQT's default engine.
   //GeometryQueryTool::instance()->add_gqe(this);
   //GeometryQueryTool::instance()->set_default_engine(this); 
  CompositeEngine::instance();
  PartitionEngine::instance();
}


//-------------------------------------------------------------------------
// Purpose       : get child virtual geometry
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck 
//
// Creation Date : 08/11/98
//-------------------------------------------------------------------------
void VirtualQueryEngine::get_VEs( TopologyEntity* te_ptr,
                                     DLIList<TopologyBridge*>& ve_list, 
                                     CubitBoolean visible,
                                     const CubitBoolean children_too)
{
  DLIList<TopologyBridge*> child_list;
  Body*      body_ptr   = CAST_TO( te_ptr, Body );
  CoVolume*  covol_ptr  = CAST_TO( te_ptr, CoVolume );
  RefVolume* vol_ptr    = CAST_TO( te_ptr, RefVolume );
  Shell*     shell_ptr  = CAST_TO( te_ptr, Shell );
  CoFace*    coface_ptr = CAST_TO( te_ptr, CoFace );
  RefFace*   face_ptr   = CAST_TO( te_ptr, RefFace );
  Loop*      loop_ptr   = CAST_TO( te_ptr, Loop );
  CoEdge*    coedge_ptr = CAST_TO( te_ptr, CoEdge );
  RefEdge*   edge_ptr   = CAST_TO( te_ptr, RefEdge );
  Chain*     chain_ptr  = CAST_TO( te_ptr, Chain );
  CoVertex*  covex_ptr  = CAST_TO( te_ptr, CoVertex );
  RefVertex* vex_ptr    = CAST_TO( te_ptr, RefVertex );
  
  if( (body_ptr != NULL) || (covol_ptr != NULL) )
  {
    if (children_too) {
      DLIList<RefVolume*> vol_list;
      te_ptr->ref_volumes( vol_list );
      for( int i = vol_list.size(); i > 0; i-- )
      {
        child_list.clean_out();
        get_VEs( vol_list.get_and_step(), child_list, visible,
                 children_too);
        ve_list.merge_unique( child_list );
      }
    }
  }
  else if( vol_ptr != NULL )
  {
    get_VEs( vol_ptr, ve_list, visible, children_too );
  }
  else if( (shell_ptr != NULL) || (coface_ptr != NULL) )
  {
    if (children_too) {
      DLIList<RefFace*> face_list;
      te_ptr->ref_faces( face_list );
      for( int i = face_list.size(); i> 0; i-- )
      {
        child_list.clean_out();
        get_VEs( face_list.get_and_step(), child_list, visible,
                 children_too);
        ve_list.merge_unique( child_list );
      }
    }
  }
  else if( face_ptr != NULL )
  {
    get_VEs( face_ptr, ve_list, visible, children_too );
  }
  else if( (loop_ptr != NULL) || (coedge_ptr != NULL) )
  {
    if (children_too) {
      DLIList<RefEdge*> edge_list;
      te_ptr->ref_edges( edge_list );
      for( int i = edge_list.size(); i > 0; i-- )
      {
        child_list.clean_out();
        get_VEs( edge_list.get_and_step(), child_list, visible,
                 children_too);
        ve_list.merge_unique( child_list );
      }
    }
  }
  else if( edge_ptr != NULL )
  {
    get_VEs( edge_ptr, ve_list, visible, children_too );
  }
  else if( (chain_ptr != NULL) || (covex_ptr != NULL ) )
  {
    if (children_too) {
      DLIList<RefVertex*> vertex_list;
      te_ptr->ref_vertices( vertex_list );
      for( int i = vertex_list.size(); i > 0; i-- )
      {
        child_list.clean_out();
        get_VEs( vertex_list.get_and_step(), child_list, visible,
                 children_too);
        ve_list.merge_unique( child_list );
      }
    }
  }
  else if( vex_ptr != NULL )
  {
    get_VEs( vex_ptr, ve_list, visible, children_too );
  }
  else
	{
	  PRINT_ERROR("Unknown type of TopologyEntity passed to "
		            "VirtualQueryEngine::get_VEs(..).  Passed pointer"
								" may be stale.\n");
	}
}

void VirtualQueryEngine::get_VEs( RefVolume* volume_ptr,
                                     DLIList<TopologyBridge*>& ve_list,
                                     CubitBoolean visible,
                                     const CubitBoolean children_too )
{
#ifdef BOYD17
  DLIList<BasicTopologyEntity*> tmp_list, bte_list;
#endif  
  DLIList<BasicTopologyEntity*> bte_list;
  int i;
  
  DLIList<TopologyBridge*> tb_list;
  volume_ptr->bridge_manager()->get_bridge_list( tb_list );
  
  for( i = tb_list.size(); i--; )
  {  
    TopologyBridge* tb_ptr = tb_list.get_and_step();
    PartitionEntity* ve_ptr = dynamic_cast<PartitionEntity*>(tb_ptr);
    CompositeLump* ce_ptr = dynamic_cast<CompositeLump*>(tb_ptr);
    if( ve_ptr || ce_ptr )
    {
      ve_list.append( tb_ptr );
    }
  }
  
  if (!children_too) return;
  DLIList<TopologyBridge*> child_list;
  
  DLIList<RefFace*> face_list;
  volume_ptr->ref_faces( face_list );
  DLIList<BasicTopologyEntity*> temp_face_list;
  CAST_LIST_TO_PARENT( face_list, temp_face_list );
  bte_list.merge_unique( temp_face_list );
//  bte_list.merge_unique( face_list );
  for( i = bte_list.size();i > 0;i--)
  {
    child_list.clean_out();
    get_VEs( bte_list.get_and_step(), child_list, visible );
    ve_list.merge_unique( child_list );
  }
}

void VirtualQueryEngine::get_VEs( RefFace* face_ptr,
                                     DLIList<TopologyBridge*>& ve_list,
                                     CubitBoolean visible,
                                     const CubitBoolean children_too )
{
#ifdef BOYD17
  DLIList<BasicTopologyEntity*> tmp_list, bte_list;
#endif 
  DLIList<BasicTopologyEntity*> bte_list;
  DLIList<TopologyBridge*> tb_list;
  int i;
  
  face_ptr->bridge_manager()->get_bridge_list( tb_list );
  tb_list.reset();
  
  for( i = tb_list.size(); i--; )
  {
  
    TopologyBridge* ve_ptr = tb_list.get_and_step();
    if( dynamic_cast<PartitionEntity*>(ve_ptr) ||
        dynamic_cast<CompositeSurface*>(ve_ptr) )
    {
      ve_list.append( ve_ptr );
    }
  }
  
  
  if (!children_too) return;

  DLIList<TopologyBridge*> child_list;
  DLIList<RefEdge*> edge_list;
  face_ptr->ref_edges( edge_list );
  DLIList<BasicTopologyEntity*> temp_edge_list;
  CAST_LIST_TO_PARENT( edge_list, temp_edge_list );
  bte_list.merge_unique( temp_edge_list );

  for( i = bte_list.size();i > 0;i--)
  {
    child_list.clean_out();
    get_VEs( bte_list.get_and_step(), child_list, visible );
    ve_list.merge_unique( child_list );
  }
}

void VirtualQueryEngine::get_VEs( RefEdge* edge_ptr,
                                     DLIList<TopologyBridge*>& ve_list,
                                     CubitBoolean visible,
                                     const CubitBoolean children_too )
{
#ifdef BOYD17
  DLIList<BasicTopologyEntity*> tmp_list, bte_list;
#endif
  DLIList<BasicTopologyEntity*> bte_list;
  DLIList<TopologyBridge*> tb_list;
  int i;
  
  edge_ptr->bridge_manager()->get_bridge_list( tb_list );
  tb_list.reset();
  
  for( i = tb_list.size(); i--; )
  {
  
    TopologyBridge* ve_ptr = tb_list.get_and_step();
    if( dynamic_cast<PartitionCurve*>(ve_ptr) ||
        dynamic_cast<CompositeCurve*>(ve_ptr) )
    {
      ve_list.append( ve_ptr );
    }
  }
  
  if (!children_too) return;
  DLIList<TopologyBridge*> child_list;
  DLIList<RefVertex*> vex_list;
  edge_ptr->ref_vertices( vex_list );
  DLIList<BasicTopologyEntity*> temp_vex_list;
  CAST_LIST_TO_PARENT( vex_list, temp_vex_list );
  bte_list.merge_unique( temp_vex_list );

  for( i = bte_list.size();i > 0;i--)
  {
    child_list.clean_out();
    get_VEs( bte_list.get_and_step(), child_list, visible );
    ve_list.merge_unique( child_list );
  }
}

void VirtualQueryEngine::get_VEs( RefVertex* vex_ptr,
                                     DLIList<TopologyBridge*>& ve_list,
                                     CubitBoolean ,
                                     const CubitBoolean )
{
  DLIList<TopologyBridge*> tb_list;
  vex_ptr->bridge_manager()->get_bridge_list( tb_list );
  tb_list.reset();
  for( int i = tb_list.size(); i--; )
  {
    TopologyBridge* tb = tb_list.get_and_step();
    if( dynamic_cast<PartitionPoint*>(tb) ||
        dynamic_cast<CompositePoint*>(tb) )
      ve_list.append( tb );
  }
}


//-------------------------------------------------------------------------
// Purpose       : check if the entity is virtual or has virtual entities. 
//
// Special Notes : all is_virtual entity, has_virtual = true
//                 only body, refvolume, refface can be non-virtual but has
//                 virtual.
//
// Creator       : Jane Hu
//
// Creation Date : 10/01/04
//-------------------------------------------------------------------------
CubitBoolean VirtualQueryEngine::has_virtual(TopologyEntity* te_ptr)
{
  if (te_ptr == NULL)
    return CUBIT_FALSE;
  
  if (is_virtual(te_ptr))
    return CUBIT_TRUE;
  
    // check if children has partition or composite entities
  int i;
  Body      * body_ptr = CAST_TO (te_ptr, Body);
  RefVolume * vol_ptr  = CAST_TO (te_ptr, RefVolume);
  RefFace   * face_ptr = CAST_TO (te_ptr, RefFace);
  CubitBoolean has_vir = CUBIT_FALSE;
  
  if (body_ptr != NULL)
  {
    DLIList<RefVolume*> vol_list;
    te_ptr->ref_volumes(vol_list);
    for (i = vol_list.size(); i > 0; i--)
    {
      if (has_virtual(vol_list.get_and_step()))
      {
        has_vir = CUBIT_TRUE;
        break;
      }
    }
  }
                                                                                
  if (vol_ptr != NULL)
  {
    DLIList<RefFace*> face_list;
    te_ptr->ref_faces(face_list);
    for (i = face_list.size(); i > 0; i --)
    {
       if (has_virtual(face_list.get_and_step()))
       { 
         has_vir = CUBIT_TRUE;
         break;
       }
    }
  }
                                                                                
  if (face_ptr != NULL)
  {
    DLIList<RefEdge *> edge_list;
    te_ptr->ref_edges(edge_list);
    for (i = edge_list.size(); i > 0; i --)
    {
       if (is_virtual(edge_list.get_and_step()))
       {
         has_vir = CUBIT_TRUE;
         break;
       }
    }
  }

  return has_vir;
} 

CubitBoolean VirtualQueryEngine::is_virtual(TopologyEntity *entity,
                                               const CubitBoolean children_too)
//Check if a TopologyEntity has any topologyBridge that is composite or 
//partition type.
{
  if ( entity == NULL )
      return CUBIT_FALSE;
  
  if (CAST_TO(entity, Body) !=NULL)
  {
     int i;
     DLIList<TopologyBridge*> list;
     entity->bridge_manager()->get_bridge_list(list);
     TopologyBridge * tb;
     for (i = list.size(); i > 0; i--)
     {
        tb =  list.get_and_step();
        if (CAST_TO(tb, PartitionBody) != NULL ||
            CAST_TO(tb, CompositeBody) != NULL)
           return CUBIT_TRUE;
      }
  }

  DLIList<TopologyBridge*> ve_list;
  get_VEs(entity, ve_list, CUBIT_TRUE, children_too);
  return (ve_list.size() ? CUBIT_TRUE : CUBIT_FALSE);
}

CubitBoolean VirtualQueryEngine::is_virtual(DLIList<Body*> &entity_list,
                                               const CubitBoolean children_too)
{
  if (!children_too) return CUBIT_FALSE;
                                                                                
  int i;
  for (i = entity_list.size(); i > 0; i--)
    if (is_virtual(entity_list.get_and_step(), children_too)) return CUBIT_TRUE;                                                                                
  return CUBIT_FALSE;
}

CubitBoolean VirtualQueryEngine::is_virtual(DLIList<RefEntity*> &entity_list,
                                               const CubitBoolean children_too)
{
  int i;
  for (i = entity_list.size(); i > 0; i--) {
    RefEntity *ref_entity = entity_list.get_and_step();
    if (is_virtual(CAST_TO(ref_entity, TopologyEntity), children_too)) return CUBIT_TRUE;
  }
                                                                                
  return CUBIT_FALSE;
}

//-------------------------------------------------------------------------
// Purpose	 : Sort a DLList of RefEdges to be used in the construction
//		   of a CompositeCurve.
//
// Special Notes : Sorted using topological information
//
// Creator	 : Jason Kraftcheck
//
// Creation Date : 08/04/97
//-------------------------------------------------------------------------
CubitStatus VirtualQueryEngine::sort_edges( 
  DLIList<RefEdge*>& edge_list ) const
{
  DLIList<RefEdge*> sorted_edge_list;
  RefEdge *first, *current = NULL, *prev;
  CubitStatus return_value = CUBIT_SUCCESS;
  
  if( edge_list.size( ) <= 0 )
  {
    return CUBIT_FAILURE;
  }
  
  if( edge_list.size() == 1 )
  {
    return CUBIT_SUCCESS;
  }
  
    //Note: this part badly needs some optimization.  Given the array-based
    // implementation of the DLList class, this is a very inefficient
    // approach to sorting the edges.  Look into a different seach
    // algorithm, and possible switching functions in DLList.
  
    //add first Curve
  edge_list.reset( );
  first = prev = edge_list.remove( );
  sorted_edge_list.append( first );
  
    //sort the rest of the Curves
  int i;
  for( i = 0; i < edge_list.size( ); i++ ) //loop for each remaining RefEdge
  {
    edge_list.reset( );
      //look through list for next RefEdge
    for( int j = 0; j < edge_list.size( ); j++ )
    {
      current = edge_list.get( );
      if( prev->common_ref_vertex( current ) != NULL )
      {
        edge_list.remove( );
        sorted_edge_list.append( current );
        i--; //edge_list.size() decreases with removal of edge
        prev = current;
        break;
      }
      else if( first->common_ref_vertex( current ) != NULL )
      {
        edge_list.remove( );
        sorted_edge_list.insert_first( current );
        first = current;
        i--; //edge_list.size() decreases with removal of edge
        break;
      }
      else
      {
        edge_list.step( );
      }
    }//end for(j)
  }//end for(i)
  
    //check if all the RefEdges got sorted
  if( edge_list.size( ) > 0 )
  {
    return_value =  CUBIT_FAILURE;
    
      //Sort the rest of the edges anyway, for use by
      //other functions than composite edge creation.
    
    sort_edges( edge_list );
    edge_list.reset( );
    for( i = 0; i < edge_list.size(); i++ )
    {
      current = edge_list.get_and_step( );
      sorted_edge_list.append( current );
    }
  }
  else
  {
    return_value = CUBIT_SUCCESS;
  }
  
  edge_list.clean_out( );
  edge_list = sorted_edge_list;
  edge_list.reset( );
  return return_value;
}



// ********** BEGIN PRIVATE FUNCTIONS      **********

//-------------------------------------------------------------------------
// Purpose: Display a default error message for certain inherited functions
//-------------------------------------------------------------------------
void VirtualQueryEngine::default_error_message( 
  const char callers_name[] ) const
{
  PRINT_ERROR("A call was made to:\n");
  PRINT_ERROR("VirtualQueryEngine::");
  PRINT_ERROR(callers_name);
  PRINT_ERROR("\n");
  PRINT_ERROR("Although this function is inhereted from GeometryQueryEngine\n");
  PRINT_ERROR("it is not applicable or implememtend in VirtualQueryEngine.\n");
  PRINT_ERROR("\nTHIS IS A BUG.  This function should NEVER be called.\n");
  PRINT_ERROR("\nThis function exists only because it is pure virtual in\n");
  PRINT_ERROR("the abstract class GeometryQueryEngine, and must be\n");
  PRINT_ERROR("defined in VirtualQueryEngine for the class to be\n");
  PRINT_ERROR("instantiated.\n");
}

// ********** END PRIVATE FUNCTIONS        **********

// ********** BEGIN HELPER CLASSES         **********
// ********** END HELPER CLASSES           **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

// ********** BEGIN STATIC FUNCTIONS       **********
// ********** END STATIC FUNCTIONS         **********


	
CubitStatus 
VirtualQueryEngine::get_intersections( Curve* , Curve* ,
                                         DLIList<CubitVector*>& /*intersection_list*/,
                                         CubitBoolean /*bounded*/,
                                         CubitBoolean /*closest*/ )
{
  return CUBIT_FAILURE;
}

CubitStatus
VirtualQueryEngine::get_intersections(Curve* curve, CubitVector &point1,
                                         CubitVector &point2,
                                         DLIList<CubitVector*>& intersection_list,
                                         CubitBoolean bounded,
                                         CubitBoolean closest )
{
  DLIList<TopologyBridge*> curve_list;
  get_underlying_curves( curve, curve_list );
  int i;
  for (i = 0; i < curve_list.size(); i++)
  {
    // Get the next curve and its GME
    Curve* curve_ptr = CAST_TO(curve_list.get_and_step(), Curve);
    GeometryQueryEngine* GQE_ptr =
      curve_ptr->get_geometry_query_engine();
    GQE_ptr->get_intersections(curve_ptr, point1, point2, intersection_list,
                               bounded, closest);
  }
  return CUBIT_SUCCESS;
}

CubitStatus 
VirtualQueryEngine::get_intersections( Curve*, Surface*,
                                          DLIList<CubitVector*>& /* intersection_list */,
                                          CubitBoolean /* bounded */)
{
   return CUBIT_FAILURE;
}

//-------------------------------------------------------------------------
// Purpose       : This function get underlying curves for virtual curves
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date :
//-------------------------------------------------------------------------
CubitStatus VirtualQueryEngine::get_underlying_curves(Curve * curve_ptr,
                                        DLIList<TopologyBridge*>& curve_list)
{
   assert (curve_ptr);
   CompositeCurve *comp_curve = CAST_TO(curve_ptr, CompositeCurve);
   PartitionCurve *part_curve = CAST_TO(curve_ptr, PartitionCurve);
   if ( comp_curve )
   {
     int i;
     for (i = 0; i < comp_curve->num_curves(); i ++)
     {
       part_curve = CAST_TO(comp_curve->get_curve(i), PartitionCurve);
       if (!part_curve)
         curve_list.append_unique(CAST_TO(comp_curve->get_curve(i),
                                  TopologyBridge));
       else
       {
         if(dynamic_cast<Curve*>(part_curve->partitioned_entity()))
         {
           curve_list.append_unique(CAST_TO(part_curve->partitioned_entity(),
                                    TopologyBridge));
         }
       }
     }
   }
   else if ( part_curve )
     curve_list.append(CAST_TO(part_curve->partitioned_entity(),
                       TopologyBridge));
   return CUBIT_SUCCESS;
}

CubitStatus VirtualQueryEngine::get_underlying_surfaces(Surface * surf_ptr,
                                        DLIList<TopologyBridge*>& surf_list)
{
   assert (surf_ptr);
   CompositeSurface *comp_surf = CAST_TO(surf_ptr, CompositeSurface);
   PartitionSurface *part_surf = CAST_TO(surf_ptr, PartitionSurface);
   if ( comp_surf )
   {
     int i;
     for (i = 0; i < comp_surf->num_surfs(); i ++)
     {
       part_surf = CAST_TO(comp_surf->get_surface(i), PartitionSurface);
       if (!part_surf)
         surf_list.append_unique(CAST_TO(comp_surf->get_surface(i),
                                  TopologyBridge));
       else
       {
         if(dynamic_cast<Surface*>(part_surf->partitioned_entity()))
         {
           surf_list.append_unique(CAST_TO(part_surf->partitioned_entity(),
                                    TopologyBridge));
         }
       }
     }
   }
   else if ( part_surf )
     surf_list.append(CAST_TO(part_surf->partitioned_entity(),
                       TopologyBridge));
   return CUBIT_SUCCESS;
}

//================================================================================
// Description: Find extrema position on an entity list
// Author     : 
// Date       : 
//================================================================================
CubitStatus 
VirtualQueryEngine::entity_extrema( DLIList<GeometryEntity*> &, 
                                       const CubitVector *, 
                                       const CubitVector *,
                                       const CubitVector *, 
                                       CubitVector &,
                                       GeometryEntity *& )
{
  PRINT_ERROR("Entity extrema calculation not yet supported for virtual geometry.\n");
  return CUBIT_FAILURE;
}
 
//================================================================================
// Description: Find distance between two entities and closest positions.
// Author     : 
// Date       : 
//================================================================================
CubitStatus
VirtualQueryEngine::entity_entity_distance( GeometryEntity *ge1,
                                           GeometryEntity *ge2,
                                           CubitVector &p1, CubitVector &p2,
                                           double &distance )
{
  CompositeCurve *cc1 = dynamic_cast<CompositeCurve*>(ge1);
  CompositeCurve *cc2 = dynamic_cast<CompositeCurve*>(ge2);
  CompositeSurface *cs1 = dynamic_cast<CompositeSurface*>(ge1);
  CompositeSurface *cs2 = dynamic_cast<CompositeSurface*>(ge2);

  DLIList<GeometryEntity*> geometry_entities_1;
  DLIList<GeometryEntity*> geometry_entities_2;
  DLIList<GeometryQueryEngine*> gqes_1;
  DLIList<GeometryQueryEngine*> gqes_2;

  int i, j;
  int still_ok = 1;

  if(cc1)
  {
    for(i=cc1->num_curves()-1; i>-1 && still_ok; i--)
    {
      GeometryEntity *ge = cc1->get_curve(i);
      if(dynamic_cast<PartitionEntity*>(ge))
        still_ok = 0;
      GeometryQueryEngine *gqe = ge->get_geometry_query_engine();
      geometry_entities_1.append(ge);
      gqes_1.append(gqe);
    }
  }
  else if(cs1)
  {
    for(i=cs1->num_surfs()-1; i>-1 && still_ok; i--)
    {
      GeometryEntity *ge = cs1->get_surface(i);
      if(dynamic_cast<PartitionEntity*>(ge))
        still_ok = 0;
      GeometryQueryEngine *gqe = ge->get_geometry_query_engine();
      geometry_entities_1.append(ge);
      gqes_1.append(gqe);
    }
  }
  else
  {
    if(dynamic_cast<PartitionEntity*>(ge1))
      still_ok = 0;
    else if(dynamic_cast<Curve*>(ge1) || dynamic_cast<Surface*>(ge1))
    {
      GeometryQueryEngine *gqe = ge1->get_geometry_query_engine();
      geometry_entities_1.append(ge1);
      gqes_1.append(gqe);
    }
    else
    {
      PRINT_ERROR("Entity-entity distance not yet supported for virtual geometry.\n");
      return CUBIT_FAILURE;
    }
  }
  if(cc2)
  {
    for(i=cc2->num_curves()-1; i>-1 && still_ok; i--)
    {
      GeometryEntity *ge = cc2->get_curve(i);
      if(dynamic_cast<PartitionEntity*>(ge))
        still_ok = 0;
      GeometryQueryEngine *gqe = ge->get_geometry_query_engine();
      geometry_entities_2.append(ge);
      gqes_2.append(gqe);
    }
  }
  else if(cs2)
  {
    for(i=cs2->num_surfs()-1; i>-1 && still_ok; i--)
    {
      GeometryEntity *ge = cs2->get_surface(i);
      if(dynamic_cast<PartitionEntity*>(ge))
        still_ok = 0;
      GeometryQueryEngine *gqe = ge->get_geometry_query_engine();
      geometry_entities_2.append(ge);
      gqes_2.append(gqe);
    }
  }
  else
  {
    if(dynamic_cast<PartitionEntity*>(ge2))
      still_ok = 0;
    else if(dynamic_cast<Curve*>(ge2) || dynamic_cast<Surface*>(ge2))
    {
      GeometryQueryEngine *gqe = ge2->get_geometry_query_engine();
      geometry_entities_2.append(ge2);
      gqes_2.append(gqe);
    }
    else
    {
      PRINT_ERROR("Entity-entity distance not yet supported for virtual geometry.\n");
      return CUBIT_FAILURE;
    }
  }

  if(still_ok)
  {
    double smallest_distance = CUBIT_DBL_MAX;
    for(i=geometry_entities_1.size(); i--;)
    {
      GeometryEntity *ge1 = geometry_entities_1.get_and_step();
      GeometryQueryEngine *gqe1 = gqes_1.get_and_step();
      for(j=geometry_entities_2.size(); j--;)
      {
        GeometryEntity *ge2 = geometry_entities_2.get_and_step();
        GeometryQueryEngine *gqe2 = gqes_2.get_and_step();

        CubitVector cur_pos1, cur_pos2;
        double cur_distance;

        gqe1->entity_entity_distance(ge1, ge2, cur_pos1, cur_pos2, cur_distance);

        if(cur_distance < smallest_distance)
        {
          smallest_distance = cur_distance;
          distance = cur_distance;
          p1 = cur_pos1;
          p2 = cur_pos2;
        }
      }
    }
  }
  else
  {
    PRINT_ERROR("Entity-entity distance not yet supported for virtual geometry.\n");
    return CUBIT_FAILURE;
  }

  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Remove all possible virtual geometry.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/28/00
//-------------------------------------------------------------------------
void VirtualQueryEngine::remove_virtual_geometry( RefEntity* entity_ptr,
                                             CubitBoolean all_children )
{
  if ( Body* body_ptr = dynamic_cast<Body*>(entity_ptr) )
    remove_virtual_geometry(body_ptr, all_children);
  else if( RefVolume* vol_ptr = dynamic_cast<RefVolume*>(entity_ptr) )
    remove_virtual_geometry(vol_ptr, all_children);
  else if( RefFace* face_ptr = dynamic_cast<RefFace*>(entity_ptr) )
    remove_virtual_geometry(face_ptr);
}

void VirtualQueryEngine::remove_virtual_geometry( Body* body_ptr,
                                             bool all_children )
{
  int i;
  DLIList<RefVolume*> volumes, part_results;
  body_ptr->ref_volumes(volumes);
  for ( i = volumes.size(); i--; ) {
    RefVolume* vol = volumes.step_and_get();
    if ( ! dynamic_cast<PartitionLump*>(vol->get_lump_ptr()) )
      volumes.change_to(0);
  }
  volumes.remove_all_with_value(0);
  if ( volumes.size() )
    PartitionTool::instance()->unpartitionAll(volumes, part_results);
  
  if ( !all_children )
    return;
    
  volumes.clean_out();
  body_ptr->ref_volumes(volumes);
  for ( i = volumes.size(); i--; )
    remove_virtual_geometry( volumes.step_and_get(), true );
}

void VirtualQueryEngine::remove_virtual_geometry( RefVolume* vol_ptr,
                                             bool all_children )
{
  int i;
  DLIList<RefFace*> faces, part_results;
  vol_ptr->ref_faces( faces );
  for ( i = faces.size(); i--; ) {
    if ( dynamic_cast<CompositeSurface*>(faces.step_and_get()->get_surface_ptr()) )
      CompositeTool::instance()->uncomposite(faces.get());
  }
  
  faces.clean_out();
  vol_ptr->ref_faces( faces );
  for ( i = faces.size(); i--; ) {
    if ( !dynamic_cast<PartitionSurface*>(faces.step_and_get()->get_surface_ptr()))
      faces.change_to(0);
  }
  faces.remove_all_with_value(0);
  if ( faces.size() )
    PartitionTool::instance()->unpartitionAll(faces, part_results);
  
  if ( !all_children )  
    return;
  
  faces.clean_out();
  vol_ptr->ref_faces( faces );
  for ( i = faces.size(); i--; )
    remove_virtual_geometry( faces.step_and_get() );
}

void VirtualQueryEngine::remove_virtual_geometry( RefFace* face_ptr )
{
  int i;
  DLIList<RefEdge*> edges, part_results;
  face_ptr->ref_edges( edges );
  for ( i = edges.size(); i--; ) {
    if ( dynamic_cast<CompositeCurve*>(edges.step_and_get()->get_curve_ptr()) )
      CompositeTool::instance()->uncomposite(edges.get());
  }
  
  edges.clean_out();
  face_ptr->ref_edges( edges );
  for ( i = edges.size(); i--; ) {
    if ( !dynamic_cast<PartitionCurve*>(edges.step_and_get()->get_curve_ptr()))
      edges.change_to(0);
  }
  edges.remove_all_with_value(0);
  if ( edges.size() )
    PartitionTool::instance()->unpartitionAll(edges, part_results);
}
  
   
  
       
CubitStatus VirtualQueryEngine::export_solid_model( DLIList<TopologyBridge*>& ,
                                        const char* ,
                                        const char* ,
                                        const CubitString &,
                                                       const char* )  
{
  return CUBIT_FAILURE;
}

CubitStatus VirtualQueryEngine::save_temp_geom_file( DLIList<TopologyBridge*>& ,
                                        const char* ,
                                        const CubitString &,
                                        CubitString &, 
                                        CubitString &) 
{
  return CUBIT_FAILURE;
}
  
CubitStatus VirtualQueryEngine::import_temp_geom_file(FILE* , 
                                 const char* ,
                                 const char* ,
                                 DLIList<TopologyBridge*> &)
{ return CUBIT_FAILURE; }

CubitStatus VirtualQueryEngine::import_solid_model(
                                                      const char* ,
                                                      const char* ,
                                                      DLIList<TopologyBridge*>&,
                                                      CubitBoolean ,
                                                      const char* ,
                                                      CubitBoolean ,
                                                      CubitBoolean ,
                                                      CubitBoolean ,
                                                      CubitBoolean ,
                                                      CubitBoolean ,
                                                      CubitBoolean  )  
{
  PRINT_INFO("VirtualQueryEngine::import_solid_model\n");
  
  default_error_message( "import_solid_model()");
  return CUBIT_FAILURE;
}

void VirtualQueryEngine::delete_solid_model_entities(DLIList<BodySM*>& list) const
{
  for ( int i = list.size(); i--; )
    delete_solid_model_entities( list.get_and_step() );
}
      


CubitStatus VirtualQueryEngine::delete_solid_model_entities(
  BodySM* bodysm ) const
{ 
  CompositeBody* compbod = dynamic_cast<CompositeBody*>(bodysm);
  if (compbod)
  {
    while (compbod->num_bodies())
    {
      BodySM* dead_body = compbod->get_body( compbod->num_bodies() - 1 );
      dead_body->get_geometry_query_engine()->delete_solid_model_entities( dead_body );
    }
  
    CompositeEngine::instance().clean_out_deactivated_geometry();
    return CUBIT_SUCCESS;
  }
  
  
  PartitionBody* bod = dynamic_cast<PartitionBody*>(bodysm);
  if ( !bod )
    return CUBIT_FAILURE;
  
  PartitionEngine::instance().delete_solid_model_entities( bod, bodysm );
  if ( bodysm ) {
    bodysm->get_geometry_query_engine()->delete_solid_model_entities(bodysm);
  }
  return CUBIT_SUCCESS;
}

CubitStatus VirtualQueryEngine::delete_solid_model_entities(
  Surface* surf_ptr ) const
{
  PartitionSurface* ps_ptr = dynamic_cast<PartitionSurface*>(surf_ptr);
  CompositeSurface* cs_ptr = dynamic_cast<CompositeSurface*>(surf_ptr);
  
  if (cs_ptr)
  {
    while(cs_ptr->num_surfs())
    {
      Surface* dead_surf = cs_ptr->get_surface(cs_ptr->num_surfs()-1);
      dead_surf->get_geometry_query_engine()->delete_solid_model_entities(dead_surf);
    }
    
    CompositeEngine::instance().clean_out_deactivated_geometry();
  }
  
  else if( ps_ptr )
  {
    PartitionEngine::instance().delete_solid_model_entities( ps_ptr, surf_ptr );
    if ( surf_ptr )
      surf_ptr->get_geometry_query_engine()->delete_solid_model_entities(surf_ptr);
  }
  
  return CUBIT_SUCCESS;
}

CubitStatus VirtualQueryEngine::delete_solid_model_entities(
  Curve* curve_ptr ) const
{
  CompositeCurve* cc_ptr = dynamic_cast<CompositeCurve*>(curve_ptr);
  PartitionCurve* pc_ptr = dynamic_cast<PartitionCurve*>(curve_ptr);

  if( cc_ptr )
  {
    while(cc_ptr->num_curves())
    {
      Curve* dead_curv = cc_ptr->get_curve(cc_ptr->num_curves()-1);
      dead_curv->get_geometry_query_engine()->delete_solid_model_entities(dead_curv);
    }
    
    CompositeEngine::instance().clean_out_deactivated_geometry();
  }
  else if( pc_ptr )
  {
    PartitionEngine::instance().delete_solid_model_entities( pc_ptr, curve_ptr );
    if ( curve_ptr )
      curve_ptr->get_geometry_query_engine()->delete_solid_model_entities(curve_ptr);
  }
 
  return CUBIT_SUCCESS;
}

CubitStatus VirtualQueryEngine::delete_solid_model_entities( Point* ) const
{
  return CUBIT_FAILURE;
}

CubitStatus VirtualQueryEngine::fire_ray(BodySM *,
                                            const CubitVector &,
                                            const CubitVector &,
                                            DLIList<double>& ,
                                            DLIList<GeometryEntity*> *) const
{
  PRINT_INFO("VirtualQueryEngine::fire_ray\n");
  
  default_error_message( "fire_ray()");
  return CUBIT_FAILURE;
}

double VirtualQueryEngine::get_sme_resabs_tolerance() const
{
  return GEOMETRY_RESABS;
}

double VirtualQueryEngine::set_sme_resabs_tolerance( double )
{
  PRINT_INFO("VirtualQueryEngine::set_sme_resabs_tolerance\n");
  
  default_error_message( "set_sme_resabs_tolerance()");
  return 0;
}

CubitStatus VirtualQueryEngine::set_int_option( const char* , int )
{
  PRINT_INFO("VirtualQueryEngine::set_int_option\n");
  
  default_error_message( "set_int_option()");
  return CUBIT_FAILURE;
}

CubitStatus VirtualQueryEngine::set_dbl_option( const char* , double )
{
  PRINT_INFO("VirtualQueryEngine::set_dbl_option\n");
  
  default_error_message( "set_dbl_option()");
  return CUBIT_FAILURE;
}

CubitStatus VirtualQueryEngine::set_str_option( const char* , const char* )
{
  PRINT_INFO("VirtualQueryEngine::set_str_option\n");
  
  default_error_message( "set_str_option()");
  return CUBIT_FAILURE;
}

void VirtualQueryEngine::register_attributes()
{
  CubitStatus result;
  result = CGMApp::instance()->attrib_manager()->register_attrib_type(
                                               CA_PARTITION_VG, "partition vg", "PARTITION_VG",
                                               CAPartitionVG_creator, CUBIT_FALSE,
                                               CUBIT_FALSE, CUBIT_TRUE, CUBIT_TRUE,
                                               CUBIT_FALSE, CUBIT_FALSE);
  assert (CUBIT_SUCCESS == result);

  result = CGMApp::instance()->attrib_manager()->register_attrib_type(
                                               CA_COMPOSITE_VG, "composite vg", "COMPOSITE_VG",
                                               CACompositeVG_creator, CUBIT_FALSE,
                                               CUBIT_FALSE, CUBIT_TRUE, CUBIT_TRUE,
                                               CUBIT_FALSE, CUBIT_FALSE);
  assert (CUBIT_SUCCESS == result);

  result = CGMApp::instance()->attrib_manager()->register_attrib_type(
                                               CA_VIRTUAL_VG, "virtual vg", "VIRTUAL_VG",
                                               CAVirtualVG_creator, CUBIT_FALSE,
                                               CUBIT_FALSE, CUBIT_TRUE, CUBIT_TRUE,
                                               CUBIT_FALSE, CUBIT_FALSE);
  assert (CUBIT_SUCCESS == result);

}



CubitStatus VirtualQueryEngine::translate( BodySM* body, 
                                              const CubitVector& offset )
{
  if (PartitionBody* part = dynamic_cast<PartitionBody*>(body))
    return PartitionEngine::instance().translate( part, offset );
  else if(CompositeBody* comp = dynamic_cast<CompositeBody*>(body))
    return CompositeEngine::instance().translate( comp, offset );

  PRINT_ERROR("Non-virtual entity in VGE::translate.\n");
  return CUBIT_FAILURE;
}

CubitStatus VirtualQueryEngine::rotate   ( BodySM* body, 
                                              const CubitVector& axis, 
                                              double angle )
{
  if (PartitionBody* part = dynamic_cast<PartitionBody*>(body))
    return PartitionEngine::instance().rotate( part, axis, angle );
  else if(CompositeBody* comp = dynamic_cast<CompositeBody*>(body))
    return CompositeEngine::instance().rotate( comp, axis, angle );

  PRINT_ERROR("Non-virtual entity in VGE::rotate.\n");
  return CUBIT_FAILURE;
}

CubitStatus VirtualQueryEngine::scale    ( BodySM* body, 
                                              double factor )
{
  return scale( body, CubitVector(factor,factor,factor) );
}

CubitStatus VirtualQueryEngine::scale    ( BodySM* body, 
                                              const CubitVector& factors )
{
  if (PartitionBody* part = dynamic_cast<PartitionBody*>(body))
    return PartitionEngine::instance().scale( part, factors );
  else if(CompositeBody* comp = dynamic_cast<CompositeBody*>(body))
    return CompositeEngine::instance().scale( comp, factors );

  PRINT_ERROR("Non-virtual entity in VGE::scale.\n");
  return CUBIT_FAILURE;
}

CubitStatus VirtualQueryEngine::reflect  ( BodySM* body, 
                                              const CubitVector& axis )
{
  if (PartitionBody* part = dynamic_cast<PartitionBody*>(body))
    return PartitionEngine::instance().reflect( part, axis );
  else if(CompositeBody* comp = dynamic_cast<CompositeBody*>(body))
    return CompositeEngine::instance().reflect( comp, axis );

  PRINT_ERROR("Non-virtual entity in VGE::reflect.\n");
  return CUBIT_FAILURE;
}

CubitStatus VirtualQueryEngine::restore_transform( BodySM* body )
{
  if (PartitionBody* part = dynamic_cast<PartitionBody*>(body))
    return PartitionEngine::instance().restore_transform( part );
  else if(CompositeBody* comp = dynamic_cast<CompositeBody*>(body))
    return CompositeEngine::instance().restore_transform( comp );

  PRINT_ERROR("Non-virtual entity in VGE::restore_transform.\n");
  return CUBIT_FAILURE;
}

CubitStatus VirtualQueryEngine::translate( GeometryEntity* ent, 
                                              const CubitVector& offset )
{
  if (PartitionEntity* part = dynamic_cast<PartitionEntity*>(ent))
    return PartitionEngine::instance().translate( part, offset );
  else if(CompositeSurface* csurf = dynamic_cast<CompositeSurface*>(ent))
    return CompositeEngine::instance().translate( csurf, offset );
  else if(CompositeCurve* ccurve = dynamic_cast<CompositeCurve*>(ent))
    return CompositeEngine::instance().translate( ccurve, offset );

  PRINT_ERROR("Non-virtual entity in VGE::translate.\n");
  return CUBIT_FAILURE;
}

CubitStatus VirtualQueryEngine::rotate   ( GeometryEntity* ent, 
                                              const CubitVector& axis, 
                                              double degrees )
{
  if (PartitionEntity* part = dynamic_cast<PartitionEntity*>(ent))
    return PartitionEngine::instance().rotate( part, axis, degrees );
  else if(CompositeSurface* csurf = dynamic_cast<CompositeSurface*>(ent))
    return CompositeEngine::instance().rotate( csurf, axis, degrees );
  else if(CompositeCurve* ccurve = dynamic_cast<CompositeCurve*>(ent))
    return CompositeEngine::instance().rotate( ccurve, axis, degrees );

  PRINT_ERROR("Non-virtual entity in VGE::rotate.\n");
  return CUBIT_FAILURE;
}

CubitStatus VirtualQueryEngine::scale    ( GeometryEntity* ent, 
                                              double factor )
{
  return scale( ent, CubitVector( factor, factor, factor ) );
}

CubitStatus VirtualQueryEngine::scale    ( GeometryEntity* ent, 
                                              const CubitVector& factors )
{
  if (PartitionEntity* part = dynamic_cast<PartitionEntity*>(ent))
    return PartitionEngine::instance().scale( part, factors );
  else if(CompositeSurface* csurf = dynamic_cast<CompositeSurface*>(ent))
    return CompositeEngine::instance().scale( csurf, factors );
  else if(CompositeCurve* ccurve = dynamic_cast<CompositeCurve*>(ent))
    return CompositeEngine::instance().scale( ccurve, factors );

  PRINT_ERROR("Non-virtual entity in VGE::scale.\n");
  return CUBIT_FAILURE;
}

CubitStatus VirtualQueryEngine::reflect  ( GeometryEntity* ent, 
                                              const CubitVector& axis )
{
  if (PartitionEntity* part = dynamic_cast<PartitionEntity*>(ent))
    return PartitionEngine::instance().reflect( part, axis );
  else if(CompositeSurface* csurf = dynamic_cast<CompositeSurface*>(ent))
    return CompositeEngine::instance().reflect( csurf, axis );
  else if(CompositeCurve* ccurve = dynamic_cast<CompositeCurve*>(ent))
    return CompositeEngine::instance().reflect( ccurve, axis );

  PRINT_ERROR("Non-virtual entity in VGE::reflect.\n");
  return CUBIT_FAILURE;
}

//Function should only be called on volumes of multi-volume bodies
CubitBoolean VirtualQueryEngine::volumes_overlap (Lump *lump1, Lump *lump2 ) const
{
  PartitionLump* partition_lump1 = dynamic_cast<PartitionLump*>(lump1);
  PartitionLump* partition_lump2 = dynamic_cast<PartitionLump*>(lump2);

  Lump *tmp_lump1 = lump1; 
  if( partition_lump1 )
    tmp_lump1 = partition_lump1->real_lump();

  Lump *tmp_lump2 = lump2; 
  if( partition_lump2 )
    tmp_lump2 = partition_lump2->real_lump();

  GeometryQueryEngine *gqe = tmp_lump1->get_geometry_query_engine();
  if( gqe != tmp_lump2->get_geometry_query_engine() )
  {
    PRINT_ERROR("Volumes must be of the same type (ACIS, SolidWorks, etc) to\n"
                "find if they overlap.\n");
    return CUBIT_FALSE;
  }

  return gqe->volumes_overlap( tmp_lump1, tmp_lump2 );
}


CubitBoolean VirtualQueryEngine::bodies_overlap (BodySM *body_ptr_1, BodySM *body_ptr_2 ) const
{
  PartitionBody* partition_body1 = dynamic_cast<PartitionBody*>(body_ptr_1);
  PartitionBody* partition_body2 = dynamic_cast<PartitionBody*>(body_ptr_2);

  BodySM *body1 = body_ptr_1;
  if( partition_body1 )
    body1 = partition_body1->real_body();

  BodySM *body2 = body_ptr_2;
  if( partition_body2 )
    body2 = partition_body2->real_body();
  
  GeometryQueryEngine *gqe = body1->get_geometry_query_engine();
  if( gqe != body2->get_geometry_query_engine() )
  {
    PRINT_ERROR("Volumes must be of the same type (ACIS, SolidWorks, etc) to\n"
                "find if they overlap.\n");
    return CUBIT_FALSE;
  }

  return gqe->bodies_overlap( body1, body2 );

}


