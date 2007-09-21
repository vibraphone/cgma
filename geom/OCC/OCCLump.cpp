//-------------------------------------------------------------------------
// Filename      : OCCLump.cpp
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
#include "OCCQueryEngine.hpp"
#include "OCCLump.hpp"
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

#include "OCCBody.hpp"
#include "OCCShell.hpp"
#include "OCCSurface.hpp"
#include "OCCLoop.hpp"
#include "OCCCoEdge.hpp"
#include "OCCCurve.hpp"
#include "OCCPoint.hpp"

#include <TopExp.hxx>
#include <TopTools_IndexedMapOfShape.hxx>

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
OCCLump::OCCLump(TopoDS_Solid *theSolid)
{
  myTopoDSSolid = theSolid;
}

OCCLump::OCCLump(DLIList<ShellSM*> &shells,
                     BodySM *body_ptr )
{
  myBodyPtr = body_ptr;
  myShells += shells;
}


OCCLump::~OCCLump()
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
void OCCLump::append_simple_attribute_virt(CubitSimpleAttrib *csa)
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
void OCCLump::remove_simple_attribute_virt(CubitSimpleAttrib *csa )
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
void OCCLump::remove_all_simple_attribute_virt()
{ attribSet.remove_all_attributes(); }

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to get the  
//                 attributes attached to this geometry entity. The name is 
//                 attached to the underlying BODY this points to.
//
// Special Notes : 
//
//-------------------------------------------------------------------------
CubitStatus OCCLump::get_simple_attribute(DLIList<CubitSimpleAttrib*>& csa_list)
  { return attribSet.get_attributes( csa_list ); }

CubitStatus OCCLump::get_simple_attribute( const CubitString& name,
                                        DLIList<CubitSimpleAttrib*>& csa_list )
  { return attribSet.get_attributes( name, csa_list ); }

CubitStatus OCCLump::save_attribs( FILE *file_ptr )
  { return attribSet.save_attributes( file_ptr ); }

CubitStatus OCCLump::restore_attribs( FILE *file_ptr, unsigned int endian )
  { return attribSet.restore_attributes( file_ptr, endian ); }




//-------------------------------------------------------------------------
// Purpose       : Get the bounding box of the object.
//
// Special Notes :
//
//-------------------------------------------------------------------------
CubitBox OCCLump::bounding_box() const 
{
  CubitBox my_box, temp_box;
  DLIList<OCCSurface*> surfaces;
  int ii;
  const_cast<OCCLump*>(this)->get_surfaces(surfaces);
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
                 OCCLump::get_geometry_query_engine() const
{
   return OCCQueryEngine::instance();
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
double OCCLump::measure()
{
  DLIList<CubitFacet*> bounding_facets;
  DLIList<CubitPoint*> bounding_points;
  DLIList<OCCSurface*> surfaces;
  Surface *curr_surface;
  OCCSurface *facet_surface;
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
      facet_surface = CAST_TO(curr_surface, OCCSurface);
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
void OCCLump::get_parents_virt(DLIList<TopologyBridge*> &bodies) 
{
  if (myBodyPtr != NULL )
    bodies.append_unique(myBodyPtr);
}

void OCCLump::get_children_virt(DLIList<TopologyBridge*> &shellsms)
{
  TopTools_IndexedMapOfShape M;
  TopExp::MapShapes(*myTopoDSSolid, TopAbs_SHELL, M);
  int ii;
  for (ii=1; ii<=M.Extent(); ii++) {
	  TopologyBridge *shell = OCCQueryEngine::occ_to_cgm(M(ii));
	  shellsms.append_unique(shell);
  }
}

void OCCLump::get_bodies( DLIList<OCCBody*>& result_list )
{
  OCCBody* body = dynamic_cast<OCCBody*>(myBodyPtr);
  if (body)
    result_list.append(body);
}

void OCCLump::get_shells( DLIList<OCCShell*>& result_list )
{
  TopTools_IndexedMapOfShape M;
  TopExp::MapShapes(*myTopoDSSolid, TopAbs_SHELL, M);
  int ii;
  for (ii=1; ii<=M.Extent(); ii++) {
	  TopologyBridge *shell = OCCQueryEngine::occ_to_cgm(M(ii));
	  result_list.append_unique(dynamic_cast<OCCShell*>(shell));
  }
}

void OCCLump::get_surfaces( DLIList<OCCSurface*>& result_list )
{
  DLIList<OCCShell*> shell_list;
  DLIList<OCCSurface*> tmp_list;
  get_shells(shell_list);
  shell_list.reset();
  for ( int i = 0; i < shell_list.size(); i++ )
  {
    tmp_list.clean_out();
    shell_list.next(i)->get_surfaces( tmp_list );
    result_list.merge_unique( tmp_list );
  }
}


void OCCLump::get_coedges( DLIList<OCCCoEdge*>& result_list )
{
  DLIList<OCCSurface*> surface_list;
  get_surfaces( surface_list );
  surface_list.reset();
  for ( int i = 0; i < surface_list.size(); i++ )
    surface_list.next(i)->get_coedges( result_list );
}

void OCCLump::get_curves( DLIList<OCCCurve*>& result_list )
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


void OCCLump::add_shell( OCCShell *shell )
{
    ShellSM* sm_ptr; 
    sm_ptr = dynamic_cast<ShellSM*>(shell);
    
    if( sm_ptr )
      myShells.append( sm_ptr );

    shell->add_lump( this );
}

void OCCLump::remove_shell( OCCShell *shell )
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
void OCCLump::disconnect_all_shells()
{
  myShells.reset();
  for (int i = myShells.size(); i--; )
  {
    ShellSM* sm_ptr = myShells.get_and_step();
    OCCShell* shell = dynamic_cast<OCCShell*>(sm_ptr);
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
CubitStatus OCCLump::mass_properties( CubitVector& centroid, double& volume )
{
  int i;
  
  DLIList<OCCShell*> shells( myShells.size() );
  CAST_LIST( myShells, shells, OCCShell );
  assert( myShells.size() == shells.size() );
  
  DLIList<OCCSurface*> surfaces;
  DLIList<OCCShell*> surf_shells;
  get_surfaces( surfaces );
  
  DLIList<CubitFacet*> facets, surf_facets;
  DLIList<CubitPoint*> junk;
  DLIList<CubitSense> senses;
  for (i = surfaces.size(); i--; )
  {
    OCCSurface* surf = surfaces.step_and_get();
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

CubitPointContainment OCCLump::point_containment( const CubitVector &point )
{
  CubitPointContainment pc_value; 
  OCCShell *facet_shell;

  int i;
  for(i=myShells.size(); i--;)
  {
    facet_shell = dynamic_cast<OCCShell*>(myShells.get_and_step()); 
    pc_value = facet_shell->point_containment( point );
    if( pc_value == CUBIT_PNT_OUTSIDE )
      return CUBIT_PNT_OUTSIDE;
    else if( pc_value == CUBIT_PNT_BOUNDARY )
      return CUBIT_PNT_BOUNDARY;
  }

  return CUBIT_PNT_INSIDE;
  
}

//Determine whether this lump is really a sheet (that is, sheet-body).
CubitBoolean OCCLump::is_sheet( )
{
  OCCShell *facet_shell;
  int i;
    //if any of the shells are sheets, the body is assume to be a sheet
    // for our purposes...
  for(i=myShells.size(); i--;)
  {
    facet_shell = dynamic_cast<OCCShell*>(myShells.get_and_step()); 
    if(facet_shell->is_sheet()){
      return CUBIT_TRUE;
    }
  }
  return CUBIT_FALSE;
}

