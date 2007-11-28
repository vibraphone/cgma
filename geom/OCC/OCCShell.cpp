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
#include "TopTools_IndexedDataMapOfShapeListOfShape.hxx"
#include "TopoDS.hxx"
#include "TopTools_ListIteratorOfListOfShape.hxx"
#include "TopTools_DataMapOfShapeInteger.hxx"

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
  myTopoDSShell = theShell;
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

//-------------------------------------------------------------------------
// Purpose       : Query solid modeler topology
//
// Special Notes : 
//
// Author        : Jane Hu 
//
// Creation Date : 11/28/07
//-------------------------------------------------------------------------
void OCCShell::get_parents_virt( DLIList<TopologyBridge*>& parents ) 
{ 
  OCCBody * body = NULL;
  DLIList <OCCBody* > *bodies = OCCQueryEngine::BodyList;
  TopTools_IndexedDataMapOfShapeListOfShape M;
  for(int i = 0; i <  bodies->size(); i++)
  {
     body = bodies->get_and_step();
     TopExp::MapShapesAndAncestors(*(body->get_TopoDS_Shape()),
				   TopAbs_SHELL, TopAbs_SOLID, M);
     const TopTools_ListOfShape& ListOfShapes = 
				M.FindFromKey(*(get_TopoDS_Shell()));
     if (!ListOfShapes.IsEmpty()) 
     {
         TopTools_ListIteratorOfListOfShape it(ListOfShapes) ;
         for (;it.More(); it.Next())
         {
	   TopoDS_Solid Solid = TopoDS::Solid(it.Value());
           int k = OCCQueryEngine::OCCMap->Find(Solid);
	   parents.append((OCCLump*)(OCCQueryEngine::OccToCGM->find(k))->second);
	 }
     } 
  }
}


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

