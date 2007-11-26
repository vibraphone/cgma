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
#include "config.h"
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
#include "Bnd_Box.hxx"
#include "BRepBndLib.hxx"
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
  Bnd_Box box;
  const TopoDS_Shape shape=*myTopoDSSolid;
  //calculate the bounding box
  BRepBndLib::Add(shape, box);
  double min[3], max[3];

  //get values
  box.Get(min[0], min[1], min[2], max[0], max[1], max[2]);

  //update boundingbox.
  CubitBox cBox(min, max);
  return cBox;
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
      //facet_surface->get_my_facets(bounding_facets, bounding_points);
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


