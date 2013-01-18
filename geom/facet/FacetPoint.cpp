//-------------------------------------------------------------------------
// Filename      : FacetPoint.cpp
//
// Purpose       : 
//
// Special Notes :
//
// Creator       : Steven J. Owen
//
// Creation Date : 07/15/00
//
// Owner         : Steven J. Owen
//-------------------------------------------------------------------------

// ********** BEGIN STANDARD INCLUDES      **********
#include <assert.h>
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN CUBIT INCLUDES         **********
#include "CubitMessage.hpp"
#include "FacetPoint.hpp"
#include "FacetAttrib.hpp"
#include "CubitPoint.hpp"
#include "CubitPointData.hpp"
#include "RefVertex.hpp"
#include "FacetQueryEngine.hpp"
#include "CastTo.hpp"
#include "Curve.hpp"


#include "FacetBody.hpp"
#include "FacetLump.hpp"
#include "FacetShell.hpp"
#include "FacetSurface.hpp"
#include "FacetLoop.hpp"
#include "FacetCoEdge.hpp"
#include "FacetCurve.hpp"

// ********** END CUBIT INCLUDES           **********

// ********** BEGIN STATIC DECLARATIONS    **********
// ********** END STATIC DECLARATIONS      **********

// ********** BEGIN PUBLIC FUNCTIONS       **********

//-------------------------------------------------------------------------
// Purpose       : The constructor with a pointer to the associated VERTEX . 
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 07/16/00
//-------------------------------------------------------------------------
FacetPoint::FacetPoint( const CubitVector &location, DLIList<Curve*> &curves )
{
  myPoint = (CubitPoint *) new CubitPointData( location );
  myCurves += curves;
  iCreated = CUBIT_TRUE;
}

//-------------------------------------------------------------------------
// Purpose       : The constructor with a pointer to the associated CubitPoint . 
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 12/28/00
//-------------------------------------------------------------------------
FacetPoint::FacetPoint( CubitPoint *thePoint, DLIList<Curve*> &curves )
{
  myPoint = thePoint;
  myCurves += curves;
  iCreated = CUBIT_FALSE;
}

//-------------------------------------------------------------------------
// Purpose       : The constructor with a pointer to the associated CubitPoint . 
//
// Special Notes : Was make especially for save/restore
//
// Creator       : Corey Ernst 
//
// Creation Date : 02/03/03
//-------------------------------------------------------------------------
FacetPoint::FacetPoint( CubitPoint *thePoint ) 
{
  myPoint = thePoint;
  iCreated = CUBIT_FALSE;
}

//-------------------------------------------------------------------------
// Purpose       : The destructor. 
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 07/16/00
//-------------------------------------------------------------------------
FacetPoint::~FacetPoint() 
{
  if (iCreated && myPoint != NULL)
  {
    delete myPoint;
  }
  myPoint = NULL;
}

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to append a
//                 attribute to the GE. The name is attached to the 
//                 underlying solid model entity this one points to.
//
//
// Special Notes : 
//
// Creator       : Steve Owen
//
// Creation Date : 07/16/00
//-------------------------------------------------------------------------
void FacetPoint::append_simple_attribute_virt(const CubitSimpleAttrib &csa)
  { attribSet.append_attribute(csa); }

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to remove a simple 
//                 attribute attached to this geometry entity. The name is 
//                 removed from the underlying BODY this points to.
//
// Special Notes : 
//
// Creator       : Steve Owen
//
// Creation Date : 07/16/00
//-------------------------------------------------------------------------
void FacetPoint::remove_simple_attribute_virt(const CubitSimpleAttrib &csa)
  { attribSet.remove_attribute(csa); }

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to remove all simple 
//                 attributes attached to this geometry entity.  Also
//                 removes lingering GTC attributes.
//
//
// Special Notes : 
//
// Creator       : Steve Owen
//
// Creation Date : 07/16/00
//-------------------------------------------------------------------------
void FacetPoint::remove_all_simple_attribute_virt()
  { attribSet.remove_all_attributes(); }

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to get the  
//                 attributes attached to this geometry entity. The name is 
//                 attached to the underlying BODY this points to.
//
// Special Notes : 
//
// Creator       : Steve Owen
//
// Creation Date : 07/16/00
//-------------------------------------------------------------------------
CubitStatus FacetPoint::get_simple_attribute(DLIList<CubitSimpleAttrib>&
                                               csa_list)
  { return attribSet.get_attributes(csa_list); }
CubitStatus FacetPoint::get_simple_attribute(const CubitString& name,
                                     DLIList<CubitSimpleAttrib>& csa_list )
  { return attribSet.get_attributes( name, csa_list ); }

CubitStatus FacetPoint::save_attribs( FILE *file_ptr )
  { return attribSet.save_attributes(file_ptr); }

CubitStatus FacetPoint::restore_attribs( FILE *file_ptr, unsigned int endian )
  { return attribSet.restore_attributes(file_ptr, endian); }


//-------------------------------------------------------------------------
// Purpose       : Returns the coordinates of this Point. 
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 07/16/00
//-------------------------------------------------------------------------
CubitVector FacetPoint::coordinates() const
{
  return myPoint->coordinates();
}

//-------------------------------------------------------------------------
// Purpose       : Get geometry modeling engine: AcisGeometryEngine
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 07/16/00
//-------------------------------------------------------------------------
GeometryQueryEngine* FacetPoint::get_geometry_query_engine() const
{
  return FacetQueryEngine::instance();
}                 

//-------------------------------------------------------------------------
// Purpose       : Get the bounding box of the object.
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 07/16/00
//-------------------------------------------------------------------------
CubitBox FacetPoint::bounding_box() const 
{
  CubitVector temp_vector = this->coordinates();
  CubitBox temp_box(temp_vector);
  return temp_box;
}


/*
void FacetPoint::bodysms(DLIList<BodySM*> &bodies)
{
  int ii;
  for ( ii = myCurves.size(); ii > 0; ii-- )
  {
    myCurves.get_and_step()->bodysms(bodies);
  }
}

void FacetPoint::lumps(DLIList<Lump*> &lumps)
{
  int ii;
  for ( ii = myCurves.size(); ii > 0; ii-- )
  {
    myCurves.get_and_step()->lumps(lumps);
  }
}

void FacetPoint::shellsms(DLIList<ShellSM*> &shellsms)
{
  int ii;
  for ( ii = myCurves.size(); ii > 0; ii-- )
  {
    myCurves.get_and_step()->shellsms(shellsms);
  }
}

void FacetPoint::surfaces(DLIList<Surface*> &surfaces)
{
  int ii;
  for ( ii = myCurves.size(); ii > 0; ii-- )
  {
    myCurves.get_and_step()->surfaces(surfaces);
  }
}

void FacetPoint::loopsms(DLIList<LoopSM*> &loopsms)
{
  int ii;
  for ( ii = myCurves.size(); ii > 0; ii-- )
  {
    myCurves.get_and_step()->loopsms(loopsms);
  }
}

void FacetPoint::coedgesms(DLIList<CoEdgeSM*> &coedgesms)
{
  int ii;
  for ( ii = myCurves.size(); ii > 0; ii-- )
  {
     myCurves.get_and_step()->coedgesms(coedgesms);
  }
}


void FacetPoint::curves(DLIList<Curve*> &curves)
{
  int ii;
  for ( ii = myCurves.size(); ii > 0; ii-- )
  {
    curves.append_unique( myCurves.get_and_step() );
  }
}

void FacetPoint::points(DLIList<Point*> &points)
{
  points.append_unique( this );
}
*/


void FacetPoint::get_parents_virt( DLIList<TopologyBridge*>& parents ) 
  { CAST_LIST_TO_PARENT( myCurves, parents ); }
void FacetPoint::get_children_virt( DLIList<TopologyBridge*>& ) 
  {  }


void FacetPoint::get_lumps( DLIList<FacetLump*>& result_list )
{
  DLIList<FacetShell*> shell_list;
  get_shells( shell_list );
  shell_list.reset();
  for ( int i = shell_list.size(); i--; )
  {
    FacetShell* shell = shell_list.get_and_step();
    shell->get_lumps( result_list );
    FacetLump* lump = dynamic_cast<FacetLump*>(shell->get_lump());
    if (lump)
      result_list.append_unique(lump);
  }
}

void FacetPoint::get_shells( DLIList<FacetShell*>& result_list )
{
  DLIList<FacetSurface*> surface_list;
  DLIList<FacetShell*> temp_list;
  get_surfaces( surface_list );
  surface_list.reset();
  for ( int i = surface_list.size(); i--; )
  {
    FacetSurface* surface = surface_list.get_and_step();
    temp_list.clean_out();
    surface->get_shells( temp_list );
    result_list.merge_unique( temp_list );
  }
}

void FacetPoint::get_surfaces( DLIList<FacetSurface*>& result_list )
{
  DLIList<FacetLoop*> loop_list;
  get_loops( loop_list );
  loop_list.reset();
  for ( int i = loop_list.size(); i--; )
  {
    FacetLoop* loop = loop_list.get_and_step();
    FacetSurface* surface = dynamic_cast<FacetSurface*>(loop->get_surface());
    if (surface)
      result_list.append_unique(surface);
  }
}

void FacetPoint::get_loops( DLIList<FacetLoop*>& result_list )
{
  DLIList<FacetCoEdge*> coedge_list;
  get_coedges( coedge_list );
  coedge_list.reset();
  for ( int i = coedge_list.size(); i--; )
  {
    FacetCoEdge* coedge = coedge_list.get_and_step();
    FacetLoop* loop = dynamic_cast<FacetLoop*>(coedge->get_loop());
    if (loop)
      result_list.append_unique(loop);
  }
}

void FacetPoint::get_coedges( DLIList<FacetCoEdge*>& result_list )
{
  DLIList<FacetCurve*> curve_list;
  get_curves( curve_list );
  curve_list.reset();
  for ( int i = curve_list.size(); i--; )
    curve_list.get_and_step()->get_coedges( result_list );
}

void FacetPoint::get_curves( DLIList<FacetCurve*>& result_list )
{
  myCurves.reset();
  for ( int i = 0; i < myCurves.size(); i++ )
    if ( FacetCurve* curve = dynamic_cast<FacetCurve*>(myCurves.next(i)) )
      result_list.append(curve);
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
CubitStatus FacetPoint::disconnect_curve (FacetCurve* curve)
{
  if (!myCurves.move_to(curve))
    return CUBIT_FAILURE;
  myCurves.remove();
  
  if (curve->start_point() == this)
    curve->remove_start_point();
  
  if (curve->end_point() == this)
    curve->remove_end_point();
  
  return CUBIT_SUCCESS;
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
