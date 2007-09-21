//-------------------------------------------------------------------------
// Filename      : OCCPoint.cpp
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
#include "OCCPoint.hpp"
#include "OCCAttrib.hpp"
#include "CubitPoint.hpp"
#include "CubitPointData.hpp"
#include "RefVertex.hpp"
#include "OCCQueryEngine.hpp"
#include "CastTo.hpp"
#include "Curve.hpp"


#include "OCCBody.hpp"
#include "OCCLump.hpp"
#include "OCCShell.hpp"
#include "OCCSurface.hpp"
#include "OCCLoop.hpp"
#include "OCCCoEdge.hpp"
#include "OCCCurve.hpp"

#include <BRep_Tool.hxx>

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
OCCPoint::OCCPoint( CubitVector &location, DLIList<Curve*> &curves )
{
  assert(0);
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
OCCPoint::OCCPoint( CubitPoint *thePoint, DLIList<Curve*> &curves )
{
  assert(0);
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
OCCPoint::OCCPoint( CubitPoint *thePoint ) 
{
  assert(0);
  myPoint = thePoint;
  iCreated = CUBIT_FALSE;
}

//-------------------------------------------------------------------------
// Purpose       : The constructor with a pointer to the associated TopoDS_Vertex . 
//
// Special Notes :
//
// Creator       : Steve Owen
//
// Creation Date : 12/28/00
//-------------------------------------------------------------------------
OCCPoint::OCCPoint( TopoDS_Vertex *thePoint, DLIList<Curve*> &curves )
{
  assert(0);
  myTopoDSVertex = thePoint;
  myCurves += curves;
  iCreated = CUBIT_FALSE;
}

//-------------------------------------------------------------------------
// Purpose       : The constructor with a pointer to the associated TopoDS_Vertex . 
//
// Special Notes : Was make especially for save/restore
//
// Creator       : Corey Ernst 
//
// Creation Date : 02/03/03
//-------------------------------------------------------------------------
OCCPoint::OCCPoint( TopoDS_Vertex *thePoint ) 
{
  myTopoDSVertex = thePoint;
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
OCCPoint::~OCCPoint() 
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
void OCCPoint::append_simple_attribute_virt(CubitSimpleAttrib *csa)
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
void OCCPoint::remove_simple_attribute_virt(CubitSimpleAttrib *csa)
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
void OCCPoint::remove_all_simple_attribute_virt()
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
CubitStatus OCCPoint::get_simple_attribute(DLIList<CubitSimpleAttrib*>&
                                               csa_list)
  { return attribSet.get_attributes(csa_list); }
CubitStatus OCCPoint::get_simple_attribute(const CubitString& name,
                                     DLIList<CubitSimpleAttrib*>& csa_list )
  { return attribSet.get_attributes( name, csa_list ); }

CubitStatus OCCPoint::save_attribs( FILE *file_ptr )
  { return attribSet.save_attributes(file_ptr); }

CubitStatus OCCPoint::restore_attribs( FILE *file_ptr, unsigned int endian )
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
CubitVector OCCPoint::coordinates() const
{
  gp_Pnt Point = BRep_Tool::Pnt(*myTopoDSVertex);
  CubitVector p(Point.X(), Point.Y(), Point.Z());
  return p;
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
GeometryQueryEngine* OCCPoint::get_geometry_query_engine() const
{
  return OCCQueryEngine::instance();
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
CubitBox OCCPoint::bounding_box() const 
{
  CubitVector temp_vector = this->coordinates();
  CubitBox temp_box(temp_vector);
  return temp_box;
}


/*
void OCCPoint::bodysms(DLIList<BodySM*> &bodies)
{
  int ii;
  for ( ii = myCurves.size(); ii > 0; ii-- )
  {
    myCurves.get_and_step()->bodysms(bodies);
  }
}

void OCCPoint::lumps(DLIList<Lump*> &lumps)
{
  int ii;
  for ( ii = myCurves.size(); ii > 0; ii-- )
  {
    myCurves.get_and_step()->lumps(lumps);
  }
}

void OCCPoint::shellsms(DLIList<ShellSM*> &shellsms)
{
  int ii;
  for ( ii = myCurves.size(); ii > 0; ii-- )
  {
    myCurves.get_and_step()->shellsms(shellsms);
  }
}

void OCCPoint::surfaces(DLIList<Surface*> &surfaces)
{
  int ii;
  for ( ii = myCurves.size(); ii > 0; ii-- )
  {
    myCurves.get_and_step()->surfaces(surfaces);
  }
}

void OCCPoint::loopsms(DLIList<LoopSM*> &loopsms)
{
  int ii;
  for ( ii = myCurves.size(); ii > 0; ii-- )
  {
    myCurves.get_and_step()->loopsms(loopsms);
  }
}

void OCCPoint::coedgesms(DLIList<CoEdgeSM*> &coedgesms)
{
  int ii;
  for ( ii = myCurves.size(); ii > 0; ii-- )
  {
     myCurves.get_and_step()->coedgesms(coedgesms);
  }
}


void OCCPoint::curves(DLIList<Curve*> &curves)
{
  int ii;
  for ( ii = myCurves.size(); ii > 0; ii-- )
  {
    curves.append_unique( myCurves.get_and_step() );
  }
}

void OCCPoint::points(DLIList<Point*> &points)
{
  points.append_unique( this );
}
*/


void OCCPoint::get_parents_virt( DLIList<TopologyBridge*>& parents ) 
  { CAST_LIST_TO_PARENT( myCurves, parents ); }
void OCCPoint::get_children_virt( DLIList<TopologyBridge*>& ) 
  {  }


void OCCPoint::get_lumps( DLIList<OCCLump*>& result_list )
{
  DLIList<OCCShell*> shell_list;
  get_shells( shell_list );
  shell_list.reset();
  for ( int i = shell_list.size(); i--; )
  {
    OCCShell* shell = shell_list.get_and_step();
    shell->get_lumps( result_list );
    OCCLump* lump = dynamic_cast<OCCLump*>(shell->get_lump());
    if (lump)
      result_list.append_unique(lump);
  }
}

void OCCPoint::get_shells( DLIList<OCCShell*>& result_list )
{
  DLIList<OCCSurface*> surface_list;
  DLIList<OCCShell*> temp_list;
  get_surfaces( surface_list );
  surface_list.reset();
  for ( int i = surface_list.size(); i--; )
  {
    OCCSurface* surface = surface_list.get_and_step();
    temp_list.clean_out();
    surface->get_shells( temp_list );
    result_list.merge_unique( temp_list );
  }
}

void OCCPoint::get_surfaces( DLIList<OCCSurface*>& result_list )
{
  DLIList<OCCLoop*> loop_list;
  get_loops( loop_list );
  loop_list.reset();
  for ( int i = loop_list.size(); i--; )
  {
    OCCLoop* loop = loop_list.get_and_step();
    OCCSurface* surface = dynamic_cast<OCCSurface*>(loop->get_surface());
    if (surface)
      result_list.append_unique(surface);
  }
}

void OCCPoint::get_loops( DLIList<OCCLoop*>& result_list )
{
  DLIList<OCCCoEdge*> coedge_list;
  get_coedges( coedge_list );
  coedge_list.reset();
  for ( int i = coedge_list.size(); i--; )
  {
    OCCCoEdge* coedge = coedge_list.get_and_step();
    OCCLoop* loop = dynamic_cast<OCCLoop*>(coedge->get_loop());
    if (loop)
      result_list.append_unique(loop);
  }
}

void OCCPoint::get_coedges( DLIList<OCCCoEdge*>& result_list )
{
  DLIList<OCCCurve*> curve_list;
  get_curves( curve_list );
  curve_list.reset();
  for ( int i = curve_list.size(); i--; )
    curve_list.get_and_step()->get_coedges( result_list );
}

void OCCPoint::get_curves( DLIList<OCCCurve*>& result_list )
{
  myCurves.reset();
  for ( int i = 0; i < myCurves.size(); i++ )
    if ( OCCCurve* curve = dynamic_cast<OCCCurve*>(myCurves.next(i)) )
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
CubitStatus OCCPoint::disconnect_curve (OCCCurve* curve)
{
  assert(0);
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
