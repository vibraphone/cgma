//-------------------------------------------------------------------------
// Filename      : OCCShell.cpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : David White
//
// Creation Date : 7/18/2000
//
//-------------------------------------------------------------------------

// ********** BEGIN STANDARD INCLUDES      **********
#include <stddef.h>
// ********** END STANDARD INCLUDES        **********
#include "config.h"
// ********** BEGIN CUBIT INCLUDES         **********
#include "CastTo.hpp"
#include "CubitUtil.hpp"

#include "OCCQueryEngine.hpp"
#include "OCCShell.hpp"
#include "ShellSM.hpp"
#include "Lump.hpp"
#include "Surface.hpp"

#include "OCCBody.hpp"
#include "OCCLump.hpp"
#include "OCCSurface.hpp"
#include "OCCLoop.hpp"
#include "OCCCoEdge.hpp"
#include "OCCCurve.hpp"
#include "OCCPoint.hpp"
#include "FacetEvalTool.hpp"
#include "CubitPoint.hpp"
#include "CubitFacetEdge.hpp"
#include "GfxDebug.hpp"

#include <TopExp.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
// ********** END CUBIT INCLUDES           **********

// ********** BEGIN STATIC DECLARATIONS    **********
// ********** END STATIC DECLARATIONS      **********

// ********** BEGIN PUBLIC FUNCTIONS       **********

//-------------------------------------------------------------------------
// Purpose       : A constructor with a pointer to a ACIS SHELL.
//
// Special Notes :
//
//-------------------------------------------------------------------------
OCCShell::OCCShell(TopoDS_Shell *theShell)
{
  myLump = NULL;
  myTopoDSShell = theShell;
}
OCCShell::OCCShell(Lump* my_lump,
                       DLIList<Surface*> &my_surfs )
{
  myLump = my_lump;
  mySurfs += my_surfs;
}

//-------------------------------------------------------------------------
// Purpose       : A constructor with list of surfaces.
//
// Special Notes : For use with save/restore
//
//-------------------------------------------------------------------------
OCCShell::OCCShell( DLIList<Surface*> &my_surfs )
{
  myLump = NULL; 
  mySurfs += my_surfs;
}
//-------------------------------------------------------------------------
// Purpose       : The destructor.
//
// Special Notes :
//
//-------------------------------------------------------------------------
OCCShell::~OCCShell()
{}


//-------------------------------------------------------------------------
// Purpose       : Get geometry modeling engine: OCCQueryEngine
//
// Special Notes :
//
//-------------------------------------------------------------------------
GeometryQueryEngine* 
                 OCCShell::get_geometry_query_engine() const
{
   return OCCQueryEngine::instance();
}                 

void OCCShell::append_simple_attribute_virt(CubitSimpleAttrib*)
{
}
void OCCShell::remove_simple_attribute_virt(CubitSimpleAttrib* )
{
}
void OCCShell::remove_all_simple_attribute_virt()
{
}
CubitStatus OCCShell::get_simple_attribute(DLIList<CubitSimpleAttrib*>&)
{
  return CUBIT_FAILURE;
}
CubitStatus OCCShell::get_simple_attribute(const CubitString&,
                                              DLIList<CubitSimpleAttrib*>&)
  { return CUBIT_FAILURE; }

/*
void OCCShell::bodysms(DLIList<BodySM*> &bodies) 
{
  if (myLump)
    myLump->bodysms(bodies);
}

void OCCShell::lumps(DLIList<Lump*> &lumps)
{
    lumps.append_unique(myLump);
}

void OCCShell::shellsms(DLIList<ShellSM*> &shellsms)
{
  shellsms.append_unique(this);
}

void OCCShell::surfaces(DLIList<Surface*> &surfaces)
{
  int ii;
  for ( ii = mySurfs.size(); ii > 0; ii-- )
  {
    surfaces.append_unique(mySurfs.get_and_step());
  }
}

void OCCShell::loopsms(DLIList<LoopSM*> &loopsms)
{
  int ii;
  for ( ii = mySurfs.size(); ii > 0; ii-- )
  {
    mySurfs.get_and_step()->loopsms(loopsms);
  }
}

void OCCShell::curves(DLIList<Curve*> &curves)
{
  int ii;
  for ( ii = mySurfs.size(); ii > 0; ii-- )
  {
    mySurfs.get_and_step()->curves(curves);
  }
}

void OCCShell::coedgesms(DLIList<CoEdgeSM*> &coedgesms)
{
  int ii;
  for ( ii = mySurfs.size(); ii > 0; ii-- )
  {
    mySurfs.get_and_step()->coedgesms(coedgesms);
  }
}

void OCCShell::points(DLIList<Point*> &points)
{
  int ii;
  for ( ii = mySurfs.size(); ii > 0; ii-- )
  {
    mySurfs.get_and_step()->points(points);
  }
}
*/

void OCCShell::add_lump(Lump* lump_ptr)
{
  assert(NULL == myLump);
  myLump = lump_ptr;
}


//-------------------------------------------------------------------------
// Purpose       : Query solid modeler topology
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 
//-------------------------------------------------------------------------
void OCCShell::get_parents_virt( DLIList<TopologyBridge*>& parents ) 
  { parents.append(myLump); }
void OCCShell::get_children_virt( DLIList<TopologyBridge*>& children )
{
  TopTools_IndexedMapOfShape M;
  TopExp::MapShapes(*myTopoDSShell, TopAbs_FACE, M);
  int ii;
  for (ii=1; ii<=M.Extent(); ii++) {
	  TopologyBridge *surface = OCCQueryEngine::occ_to_cgm(M(ii));
	  children.append_unique(surface);
  }
}


void OCCShell::get_lumps( DLIList<OCCLump*>& result_list )
{
  OCCLump* lump = dynamic_cast<OCCLump*>(myLump);
  if (lump)
    result_list.append(lump);
}

void OCCShell::get_surfaces( DLIList<OCCSurface*>& result_list )
{
  TopTools_IndexedMapOfShape M;
  TopExp::MapShapes(*myTopoDSShell, TopAbs_FACE, M);
  int ii;
  for (ii=1; ii<=M.Extent(); ii++) {
	  TopologyBridge *surface = OCCQueryEngine::occ_to_cgm(M(ii));
	  result_list.append_unique(dynamic_cast<OCCSurface*>(surface));
  }
}


void OCCShell::get_coedges( DLIList<OCCCoEdge*>& result_list )
{
  DLIList<OCCSurface*> surface_list;
  get_surfaces( surface_list );
  surface_list.reset();
  for ( int i = 0; i < surface_list.size(); i++ )
    surface_list.next(i)->get_coedges( result_list );
}

void OCCShell::get_curves( DLIList<OCCCurve*>& result_list )
{
  DLIList<OCCCoEdge*> coedge_list;
  get_coedges( coedge_list );
  coedge_list.reset();
  for ( int i = coedge_list.size(); i--; )
  {
    OCCCoEdge* coedge = coedge_list.get_and_step();
    OCCCurve* curve = dynamic_cast<OCCCurve*>(coedge->curve());
    if (curve)
      result_list.append_unique(curve);
  }
}

//-------------------------------------------------------------------------
// Purpose       : Disconnect input surfaces from "this" shell 
//
// Special Notes : 
//
// Creator       : Corey Ernst 
//
// Creation Date : 08/31/04
//-------------------------------------------------------------------------
void OCCShell::disconnect_surfaces( DLIList<OCCSurface*> &surfs_to_disconnect )
{
  for (int i = surfs_to_disconnect.size(); i--; )
  {
    OCCSurface* surface = surfs_to_disconnect.get_and_step();
    if( mySurfs.move_to( dynamic_cast<Surface*>(surface) ) )
      mySurfs.change_to(NULL);

    if (surface)
      surface->remove_shell(this);
  }
  mySurfs.remove_all_with_value( NULL );
}


//-------------------------------------------------------------------------
// Purpose       : Tear down topology
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/29/03
//-------------------------------------------------------------------------
void OCCShell::disconnect_all_surfaces()
{
  mySurfs.reset();
  for (int i = mySurfs.size(); i--; )
  {
    Surface* sm_ptr = mySurfs.get_and_step();
    OCCSurface* surface = dynamic_cast<OCCSurface*>(sm_ptr);
    if (surface)
      surface->remove_shell(this);
  }
  mySurfs.clean_out();
}

//-------------------------------------------------------------------------
// Purpose       : Flip surface-use sense
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/26/04
//-------------------------------------------------------------------------
void OCCShell::reverse()
{
  for (int i = mySurfs.size(); i--; )
  {
    OCCSurface* surf = dynamic_cast<OCCSurface*>(mySurfs.next(i));
    CubitSense sense = surf->get_shell_sense( this );
    assert( CUBIT_UNKNOWN != sense );
    surf->set_shell_sense( this, CubitUtil::opposite_sense( sense ) );
  }
}

//-------------------------------------------------------------------------
// Purpose       : Actually flip the underlying surfaces
//
// Special Notes : 
//
// Creator       : Michael Brewer
//
// Creation Date : 03/22/05
//-------------------------------------------------------------------------
void OCCShell::reverse_surfaces()
{
  for (int i = mySurfs.size(); i--; )
  {
    OCCSurface* surf = dynamic_cast<OCCSurface*>(mySurfs.next(i));
    surf->reverse_sense();
  }
}
//-------------------------------------------------------------------------
// Purpose       : Determines if point is contained within shell 
//
// Special Notes : If the shell is a void, it contains all points it doesn't enclose 
//
// Creator       : Corey Ernst 
//
// Creation Date : 09/01/04
//-------------------------------------------------------------------------
CubitPointContainment OCCShell::point_containment( const CubitVector &input_point )
{
    //just to avoid crashes, make sure we have surfaces in this list.
  if(mySurfs.size() < 1){
    return CUBIT_PNT_OUTSIDE;
  }
  CubitVector closest_location, point = input_point;
  OCCSurface *closest_surf = dynamic_cast<OCCSurface*>(mySurfs.get_and_step());
  OCCSurface *tmp_surf;
  closest_surf->closest_point_trimmed( point, closest_location ); 
  double shortest_dist = point.distance_between( closest_location );

  int i;
  //Look through all this shell's surfaces to find the closest
  //surface to input_point

  for(i=mySurfs.size()-1; i--;)
  {
    CubitVector tmp_closest_location; 
    tmp_surf = dynamic_cast<OCCSurface*>(mySurfs.get_and_step());
    tmp_surf->closest_point_trimmed( point, tmp_closest_location ); 

    double tmp_shortest_dist = point.distance_between( tmp_closest_location );
    if( tmp_shortest_dist < shortest_dist ) 
    {
      closest_location = tmp_closest_location;
      closest_surf = tmp_surf;
      shortest_dist = tmp_shortest_dist;
    }
  }

  //determine if it's on the surface, inside or outside
  if( shortest_dist <= GEOMETRY_RESABS )
    return CUBIT_PNT_BOUNDARY;

  //find the closest facet on that surface to input_point
  CubitFacet *closest_facet ; 
  //= closest_surf->get_eval_tool()->closest_facet( closest_location );

  //get the coordinates of the closest point 
  CubitPoint *pt1, *pt2;
  closest_facet->closest_point_trimmed( point, closest_location, pt1, pt2 );

  CubitVector normal;
  CubitPoint *on_point = NULL;
  //case 1: point is closest to an edge of the facet
  if( pt1 || pt2 )
  {
      //only returned one, so the it is closest to a node
    if(!pt2){
      on_point = pt1;
    }
    else if(!pt1){
      on_point = pt2;
    }
    else{//double-check that we are not closest to a single node
      
      if((pt1->coordinates().distance_between(closest_location)) <= GEOMETRY_RESABS ) 
        on_point = pt1;
      else if((pt2->coordinates().distance_between(closest_location)) <= GEOMETRY_RESABS ) 
        on_point = pt2;
    }
    
    if( on_point )
    {
      //get all facets that share this point
      DLIList<CubitFacet*> facets_sharing_point;
      on_point->facets( facets_sharing_point );

      for(i=facets_sharing_point.size(); i--; )
        normal += facets_sharing_point.get_and_step()->normal();

      //average the normals of all these facets 
      normal /= facets_sharing_point.size();
    }
    else
    {
      //case 2: get the 2 normals of the 2 neighboring facets and average them
      int index;
      CubitFacetEdge *shared_facet_edge; 
      shared_facet_edge = closest_facet->edge_from_pts( pt1, pt2, index);
      CubitFacet *other_facet = shared_facet_edge->other_facet( closest_facet );
      if(!other_facet){
        PRINT_ERROR("Edge is not connected to two facets.\n");
        normal = closest_facet->normal();
      }
      else{
        normal = ( closest_facet->normal() + other_facet->normal() ) / 2;
      }
    }
  }
  else
  {
    //case 3: just get the normal of this facet
    normal = closest_facet->normal();
  }

//   if ( closest_surf->get_relative_surface_sense() == CUBIT_REVERSED )
//   {
//     PRINT_WARNING("mbrewer:  This shouldn't happen anymore.\n");
//     normal = -1.0*( normal );
//   }
  
  ShellSM *shell_sm = static_cast<ShellSM*>(this);
  if( closest_surf->get_shell_sense( shell_sm ) == CUBIT_REVERSED )
  {
    normal = -1.0*( normal );
  }
  if( (point-closest_location)%(normal) > 0 )
    return CUBIT_PNT_OUTSIDE;
  else
    return CUBIT_PNT_INSIDE;
}

//Determine whether this shell is a sheet (ie, closed or not)
CubitBoolean OCCShell::is_sheet()
{
    //get a list of all the facets in the sheet
  DLIList<CubitFacet*> facet_list;
  int i;
  for (i = mySurfs.size(); i--; )
  {
    OCCSurface* surf = dynamic_cast<OCCSurface*>(mySurfs.next(i));
    //surf->tris(facet_list);
  }
    //should be unique... doesn't hurt anything if it isn't
  //facet_list.uniquify_ordered();
  CubitFacet* this_facet;
  CubitPoint* node_1;
  CubitPoint* node_2;
  CubitPoint* node_3;
    //if any of the facets don't have another facet across a given
    // edge... this is a sheet shell because it isn't closed.
  for (i = facet_list.size(); i--; ){
    this_facet=facet_list.get_and_step();
    if(!this_facet){
      PRINT_ERROR("Unexpected NULL pointer.");
      return CUBIT_TRUE;
    }
    this_facet->tri_nodes(node_1,node_2,node_3);
    if(!this_facet->shared_facet( node_1, node_2 ) ||
       !this_facet->shared_facet( node_1, node_3 ) ||
       !this_facet->shared_facet( node_2, node_3 ) ){
      return CUBIT_TRUE;
    }
  }
    //otherwise we made it through without finding a gap so it is close...
    // thus not a sheet.
  return CUBIT_FALSE;
}

       

// ********** END PUBLIC FUNCTIONS         **********

// ********** BEGIN PROTECTED FUNCTIONS    **********
// ********** END PROTECTED FUNCTIONS      **********

// ********** BEGIN PRIVATE FUNCTIONS      **********
// ********** END PRIVATE FUNCTIONS        **********

// ********** BEGIN HELPER CLASSES         **********
// ********** END HELPER CLASSES           **********

// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********

// ********** BEGIN STATIC FUNCTIONS       **********
// ********** END STATIC FUNCTIONS         **********

