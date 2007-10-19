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
#include "config.h"
#include "OCCPoint.hpp"
#include "OCCQueryEngine.hpp"
#include "CastTo.hpp"
#include "OCCAttribSet.hpp"
#include "CubitSimpleAttrib.hpp"
#include "BRep_Tool.hxx"

// ********** END CUBIT INCLUDES           **********

// ********** BEGIN STATIC DECLARATIONS    **********
// ********** END STATIC DECLARATIONS      **********

// ********** BEGIN PUBLIC FUNCTIONS       **********

//-------------------------------------------------------------------------
// Purpose       : The constructor with a pointer to the location   
//                 
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 10/08/07
//-------------------------------------------------------------------------
OCCPoint::OCCPoint( const CubitVector &location )
{
  gp_Pnt pt = gp_Pnt( location.x(), location.y(), location.z());
  TopoDS_Vertex theVertex = BRepBuilderAPI_MakeVertex(pt);
  myTopoDSVertex = new TopoDS_Vertex(theVertex);
}

//-------------------------------------------------------------------------
// Purpose       : The constructor with a gp_Pnt point
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 10/19/07
//-------------------------------------------------------------------------
OCCPoint::OCCPoint( gp_Pnt& thePoint )
{
  TopoDS_Vertex theVertex = BRepBuilderAPI_MakeVertex(thePoint);
  myTopoDSVertex = new TopoDS_Vertex(theVertex);
}

//-------------------------------------------------------------------------
// Purpose       : The destructor. 
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 10/08/07
//-------------------------------------------------------------------------
OCCPoint::~OCCPoint() 
{
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
  { return attribSet.get_attributes(csa_list);
  }

CubitStatus OCCPoint::get_simple_attribute(const CubitString& name,
                                     DLIList<CubitSimpleAttrib*>& csa_list )
  { return attribSet.get_attributes( name, csa_list );
  }

CubitStatus OCCPoint::save_attribs( FILE *file_ptr )
  { return attribSet.save_attributes(file_ptr); 
  }

CubitStatus OCCPoint::restore_attribs( FILE *file_ptr, unsigned int endian )
  { return attribSet.restore_attributes(file_ptr, endian); 
  }

//-------------------------------------------------------------------------
// Purpose       : Returns the coordinates of this Point. 
//
// Special Notes :
//
// Creator       : Jane Hu
//
// Creation Date : 10/08/07
//-------------------------------------------------------------------------
CubitVector OCCPoint::coordinates() const
{
  gp_Pnt pt = BRep_Tool::Pnt(*myTopoDSVertex);
  CubitVector p(pt.X(), pt.Y(), pt.Z());
  return p;
}

CubitBoolean OCCPoint::is_equal(OCCPoint & other, double Tol)
{
  gp_Pnt pt = BRep_Tool::Pnt(*myTopoDSVertex);
  TopoDS_Vertex otherVertex = *(other.get_TopoDS_Vertex());
  const gp_Pnt otherPnt = BRep_Tool::Pnt(otherVertex);
  return pt.IsEqual(otherPnt, Tol);
}   
  
double OCCPoint::distance(OCCPoint & other)
{
  gp_Pnt pt = BRep_Tool::Pnt(*myTopoDSVertex);
  return pt.Distance(BRep_Tool::Pnt(*(other.get_TopoDS_Vertex())));
}

double OCCPoint::SquareDistance (OCCPoint & other)
{
  gp_Pnt pt = BRep_Tool::Pnt(*myTopoDSVertex);
  return pt.SquareDistance(BRep_Tool::Pnt(*(other.get_TopoDS_Vertex())));
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


void OCCPoint::get_parents_virt( DLIList<TopologyBridge*>& parents ) 
  {  }
void OCCPoint::get_children_virt( DLIList<TopologyBridge*>& ) 
  {  }

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
