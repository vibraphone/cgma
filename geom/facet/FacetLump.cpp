//-------------------------------------------------------------------------
// Filename      : FacetLump.cpp
//
// Purpose       : 
//
// Creator       : David White
//
// Creation Date : 7/18/2000
//
//-------------------------------------------------------------------------

// ********** BEGIN STANDARD INCLUDES      **********
#include <assert.h>
// ********** END STANDARD INCLUDES        **********


// ********** BEGIN CUBIT INCLUDES         **********
#include "FacetQueryEngine.hpp"
#include "FacetLump.hpp"
#include "CastTo.hpp"
#include "Surface.hpp"
#include "DLIList.hpp"
#include "CubitFacet.hpp"
#include "CubitPoint.hpp"
#include "CubitVector.hpp"
#include "CubitString.hpp"
#include "ShellSM.hpp"
#include "BodySM.hpp"
#include "Body.hpp"

#include "FacetBody.hpp"
#include "FacetShell.hpp"
#include "FacetSurface.hpp"
#include "FacetLoop.hpp"
#include "FacetCoEdge.hpp"
#include "FacetCurve.hpp"
#include "FacetPoint.hpp"

// ********** END CUBIT INCLUDES           **********

// ********** BEGIN FORWARD DECLARATIONS   **********
class RefVolume;
// ********** END FORWARD DECLARATIONS     **********

// ********** BEGIN STATIC DECLARATIONS    **********
// ********** END STATIC DECLARATIONS      **********

// ********** BEGIN PUBLIC FUNCTIONS       **********

//-------------------------------------------------------------------------
// Purpose       : The constructor with a pointer to the owning shells and body.
//
// Special Notes :
//
//-------------------------------------------------------------------------
FacetLump::FacetLump(DLIList<ShellSM*> &shells,
                     BodySM *body_ptr )
{
  myBodyPtr = body_ptr;
  myShells += shells;
}


FacetLump::~FacetLump()
{}

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to append a
//                 attribute to the GE. The name is attached to the 
//                 underlying solid model entity this one points to.
//
//
// Special Notes : 
//
// Creator       : Malcolm J. Panthaki
//
// Creation Date : 11/21/96
//-------------------------------------------------------------------------
void FacetLump::append_simple_attribute_virt(CubitSimpleAttrib *csa)
  { attribSet.append_attribute(csa); }

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to remove a simple 
//                 attribute attached to this geometry entity. The name is 
//                 removed from the underlying BODY this points to.
//
// Special Notes : 
//
// Creator       : David R. White
//
// Creation Date : 03/18/97
//-------------------------------------------------------------------------
void FacetLump::remove_simple_attribute_virt(CubitSimpleAttrib *csa )
  { attribSet.remove_attribute(csa); }

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to remove all simple 
//                 attributes attached to this geometry entity.  Also
//                 removes lingering GTC attributes.
//
//
// Special Notes : 
//
// Creator       : Greg Nielson
//
// Creation Date : 07/10/98
//-------------------------------------------------------------------------
void FacetLump::remove_all_simple_attribute_virt()
{ attribSet.remove_all_attributes(); }

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to get the  
//                 attributes attached to this geometry entity. The name is 
//                 attached to the underlying BODY this points to.
//
// Special Notes : 
//
//-------------------------------------------------------------------------
CubitStatus FacetLump::get_simple_attribute(DLIList<CubitSimpleAttrib*>& csa_list)
  { return attribSet.get_attributes( csa_list ); }

CubitStatus FacetLump::get_simple_attribute( const CubitString& name,
                                        DLIList<CubitSimpleAttrib*>& csa_list )
  { return attribSet.get_attributes( name, csa_list ); }

CubitStatus FacetLump::save_attribs( FILE *file_ptr )
  { return attribSet.save_attributes( file_ptr ); }

CubitStatus FacetLump::restore_attribs( FILE *file_ptr, unsigned int endian )
  { return attribSet.restore_attributes( file_ptr, endian ); }




//-------------------------------------------------------------------------
// Purpose       : Get the bounding box of the object.
//
// Special Notes :
//
//-------------------------------------------------------------------------
CubitBox FacetLump::bounding_box() const 
{
  CubitBox my_box, temp_box;
  DLIList<FacetSurface*> surfaces;
  int ii;
  const_cast<FacetLump*>(this)->get_surfaces(surfaces);
  if (surfaces.size() > 0)
  {
    Surface* surface = surfaces.get_and_step();
    my_box = surface->bounding_box();
    for ( ii = surfaces.size(); ii > 1; ii-- )
    {
      surface = surfaces.get_and_step();
      temp_box = surface->bounding_box();
        //unite the boxes..
      my_box |= temp_box;
    }
  }
  return my_box;
}

//-------------------------------------------------------------------------
// Purpose       : Get geometry modeling engine: AcisGeometryEngine
//
// Special Notes :
//
//-------------------------------------------------------------------------
GeometryQueryEngine* 
                 FacetLump::get_geometry_query_engine() const
{
   return FacetQueryEngine::instance();
}                 

//-------------------------------------------------------------------------
// Purpose       : Returns the volume of the Lump
//
// Special Notes :
//
// Creator       : 
//
// Creation Date : 
//-------------------------------------------------------------------------
double FacetLump::measure()
{
  DLIList<CubitFacet*> bounding_facets;
  DLIList<CubitPoint*> bounding_points;
  DLIList<FacetSurface*> surfaces;
  Surface *curr_surface;
  FacetSurface *facet_surface;
    //if this is a sheet body... return 0.0
  
    //Body *tmp_body = CAST_TO(myBodyPtr->topology_entity(), Body);
  if( is_sheet() ) 
    return 0.0;
  
  int ii;
  get_surfaces(surfaces);
  if (surfaces.size() > 0)
  { 
    for ( ii = surfaces.size(); ii > 0; ii-- )
    {
      curr_surface = surfaces.get_and_step();
      facet_surface = CAST_TO(curr_surface, FacetSurface);
      if ( facet_surface == NULL )
      {
        PRINT_ERROR("Facet lump has surfaces that aren't facets?");
        return 1;
      }
      facet_surface->get_my_facets(bounding_facets, bounding_points);
    }
  }
  double volume, curr_facet_area, summation = 0.0;
  CubitFacet *curr_facet;
  CubitVector normal_of_curr_facet, vector_of_point;
  CubitPoint *point_1, *point_2, *point_3;

  for( int jj = bounding_facets.size(); jj > 0; jj-- )
  {
    curr_facet = bounding_facets.get_and_step();
    curr_facet_area = curr_facet->area();  // Current facet's area
    
    normal_of_curr_facet = curr_facet->normal(); // Current facet's normal

    curr_facet->points(point_1, point_2, point_3); // Current facet's points

    vector_of_point = point_1->coordinates(); // One point's vector

    summation += ( double(vector_of_point % normal_of_curr_facet) * curr_facet_area);
  }

  volume = summation / 3;
  
  return volume;
}
void FacetLump::get_parents_virt(DLIList<TopologyBridge*> &bodies) 
{
  if (myBodyPtr != NULL )
    bodies.append_unique(myBodyPtr);
}

void FacetLump::get_children_virt(DLIList<TopologyBridge*> &shellsms)
{
  int ii;
  for ( ii = myShells.size(); ii > 0; ii-- )
  {
    shellsms.append_unique(myShells.get_and_step());
  }
}

void FacetLump::get_bodies( DLIList<FacetBody*>& result_list )
{
  FacetBody* body = dynamic_cast<FacetBody*>(myBodyPtr);
  if (body)
    result_list.append(body);
}

void FacetLump::get_shells( DLIList<FacetShell*>& result_list )
{
  myShells.reset();
  for ( int i = 0; i < myShells.size(); i++ )
    if ( FacetShell* shell = dynamic_cast<FacetShell*>(myShells.next(i)) )
      result_list.append(shell);
}

void FacetLump::get_surfaces( DLIList<FacetSurface*>& result_list )
{
  DLIList<FacetShell*> shell_list;
  DLIList<FacetSurface*> tmp_list;
  get_shells(shell_list);
  shell_list.reset();
  for ( int i = 0; i < shell_list.size(); i++ )
  {
    tmp_list.clean_out();
    shell_list.next(i)->get_surfaces( tmp_list );
    result_list.merge_unique( tmp_list );
  }
}


void FacetLump::get_coedges( DLIList<FacetCoEdge*>& result_list )
{
  DLIList<FacetSurface*> surface_list;
  get_surfaces( surface_list );
  surface_list.reset();
  for ( int i = 0; i < surface_list.size(); i++ )
    surface_list.next(i)->get_coedges( result_list );
}

void FacetLump::get_curves( DLIList<FacetCurve*>& result_list )
{
  DLIList<FacetCoEdge*> coedge_list;
  get_coedges( coedge_list );
  coedge_list.reset();
  for ( int i = coedge_list.size(); i--; )
  {
    FacetCoEdge* coedge = coedge_list.get_and_step();
    FacetCurve* curve = dynamic_cast<FacetCurve*>(coedge->curve());
    if (curve)
      result_list.append_unique(curve);
  }
}


void FacetLump::add_shell( FacetShell *shell )
{
    ShellSM* sm_ptr; 
    sm_ptr = dynamic_cast<ShellSM*>(shell);
    
    if( sm_ptr )
      myShells.append( sm_ptr );

    shell->add_lump( this );
}

void FacetLump::remove_shell( FacetShell *shell )
{
    ShellSM* sm_ptr; 
    sm_ptr = dynamic_cast<ShellSM*>(shell);
    
    if( sm_ptr )
      myShells.remove( sm_ptr );

    shell->remove_lump();
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
void FacetLump::disconnect_all_shells()
{
  myShells.reset();
  for (int i = myShells.size(); i--; )
  {
    ShellSM* sm_ptr = myShells.get_and_step();
    FacetShell* shell = dynamic_cast<FacetShell*>(sm_ptr);
    if (shell)
    {
      assert(shell->get_lump() == this);
      shell->remove_lump();
    }
  }
  myShells.clean_out();
}

//-------------------------------------------------------------------------
// Purpose       : Calculate centroid
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/12/04
//-------------------------------------------------------------------------
#include "GfxDebug.hpp"
CubitStatus FacetLump::mass_properties( CubitVector& centroid, double& volume )
{
  int i;
  
  DLIList<FacetShell*> shells( myShells.size() );
  CAST_LIST( myShells, shells, FacetShell );
  assert( myShells.size() == shells.size() );
  
  DLIList<FacetSurface*> surfaces;
  DLIList<FacetShell*> surf_shells;
  get_surfaces( surfaces );
  
  DLIList<CubitFacet*> facets, surf_facets;
  DLIList<CubitPoint*> junk;
  DLIList<CubitSense> senses;
  for (i = surfaces.size(); i--; )
  {
    FacetSurface* surf = surfaces.step_and_get();
    surf_shells.clean_out();
    surf->get_shells( surf_shells );
    surf_shells.intersect( shells );
    assert( surf_shells.size() );
    CubitSense sense = surf->get_shell_sense( surf_shells.get() );
    if (surf_shells.size() == 1 && CUBIT_UNKNOWN != sense)
    {
      surf_facets.clean_out();
      junk.clean_out();
      surf->get_my_facets( surf_facets, junk );
      facets += surf_facets;
      
      for (int j = surf_facets.size(); j--; )
        senses.append(sense);
    }
  }
  
  const CubitVector p0 = bounding_box().center();
  CubitVector p1, p2, p3, normal;
  centroid.set( 0.0, 0.0, 0.0 );
  volume = 0.0;
  
  facets.reset();
  senses.reset();
  for (i = facets.size(); i--; )
  {
    CubitFacet* facet = facets.get_and_step();
    CubitSense sense = senses.get_and_step();
    p1 = facet->point(0)->coordinates();
    p2 = facet->point(1)->coordinates();
    p3 = facet->point(2)->coordinates();
    normal = (p3 - p1) * (p2 - p1);

    double two_area = normal.length();
    if (two_area > CUBIT_RESABS )
    {
      if (CUBIT_REVERSED == sense)
        normal = -normal;

      normal /= two_area;

      double height = normal % (p0 - p1);
      double vol = two_area * height;

      volume += vol;
      centroid += vol * (p0 + p1 + p2 + p3);
    }
  }
  
  if (volume > CUBIT_RESABS)
    centroid /= 4.0 * volume;
  volume /= 6.0;
  return CUBIT_SUCCESS;
}

CubitPointContainment FacetLump::point_containment( const CubitVector &point )
{
  CubitPointContainment pc_value; 
  FacetShell *facet_shell;

  int i;
  for(i=myShells.size(); i--;)
  {
    facet_shell = dynamic_cast<FacetShell*>(myShells.get_and_step()); 
    pc_value = facet_shell->point_containment( point );
    if( pc_value == CUBIT_PNT_OUTSIDE )
      return CUBIT_PNT_OUTSIDE;
    else if( pc_value == CUBIT_PNT_BOUNDARY )
      return CUBIT_PNT_BOUNDARY;
  }

  return CUBIT_PNT_INSIDE;
  
}

//Determine whether this lump is really a sheet (that is, sheet-body).
CubitBoolean FacetLump::is_sheet( )
{
  FacetShell *facet_shell;
  int i;
    //if any of the shells are sheets, the body is assume to be a sheet
    // for our purposes...
  for(i=myShells.size(); i--;)
  {
    facet_shell = dynamic_cast<FacetShell*>(myShells.get_and_step()); 
    if(facet_shell->is_sheet()){
      return CUBIT_TRUE;
    }
  }
  return CUBIT_FALSE;
}

