#include "CompositeLoop.hpp"
#include "DLIList.hpp"

#include "CompositeSurface.hpp"
#include "CompositeCoEdge.hpp"
#include "VirtualQueryEngine.hpp"

//-------------------------------------------------------------------------
// Purpose       : Constructor
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/11/02
//-------------------------------------------------------------------------
CompositeLoop::CompositeLoop()
  : mySurface(0), myCoedge(0), loopNext(0), numCoedges(0)
{
}

//-------------------------------------------------------------------------
// Purpose       : Destructor
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/11/02
//-------------------------------------------------------------------------
CompositeLoop::~CompositeLoop()
{
  if( mySurface )
  {
    mySurface->remove( this );
    mySurface = 0;
  }
    
  CubitStatus s = remove_all_coedges();
  assert( s );
}

//-------------------------------------------------------------------------
// Purpose       : insert a child coedge
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/23/02
//-------------------------------------------------------------------------
CubitStatus CompositeLoop::insert_after( CompositeCoEdge* insert,
                                         CompositeCoEdge* after )
{
  if( insert->myLoop != 0 )
    return CUBIT_FAILURE;
  
  if( numCoedges == 0 )
  {
    if( after )
      return CUBIT_FAILURE;
    
    insert->nextCoedge = insert;
    insert->prevCoedge = insert;
    myCoedge = insert;
  }
  else
  {
    if( !after || after->myLoop != this )
      return CUBIT_FAILURE;
    
    insert->prevCoedge = after;
    insert->nextCoedge = after->nextCoedge;
    insert->nextCoedge->prevCoedge = insert;
    insert->prevCoedge->nextCoedge = insert;
  }
  
  insert->myLoop = this;
  numCoedges++;
  return CUBIT_SUCCESS;
}
CubitStatus CompositeLoop::insert_before( CompositeCoEdge* insert,
                                          CompositeCoEdge* before )
{
  if( insert->myLoop != 0 )
    return CUBIT_FAILURE;
  
  if( numCoedges == 0 )
  {
    if( before ) 
      return CUBIT_FAILURE;
    
    insert->nextCoedge = insert;
    insert->prevCoedge = insert;
    myCoedge = insert;
  }
  else
  {
    if( !before || before->myLoop != this )
      return CUBIT_FAILURE;
    
    insert->nextCoedge = before;
    insert->prevCoedge = before->prevCoedge;
    insert->nextCoedge->prevCoedge = insert;
    insert->prevCoedge->nextCoedge = insert;
  }
  
  insert->myLoop = this;
  numCoedges++;
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : remove a child coedge
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/23/02
//-------------------------------------------------------------------------
CubitStatus CompositeLoop::remove( CompositeCoEdge* coedge )
{
  if( coedge->myLoop != this )
    return CUBIT_FAILURE;
  
  if( numCoedges == 1 )
  {
    assert( coedge->nextCoedge == coedge &&
            coedge->prevCoedge == coedge );
    myCoedge = 0;
  }
  else
  {
    coedge->nextCoedge->prevCoedge = coedge->prevCoedge;
    coedge->prevCoedge->nextCoedge = coedge->nextCoedge;
    if( myCoedge == coedge )
      myCoedge = coedge->nextCoedge;
  }
  
  numCoedges--;
  coedge->myLoop = 0;
  coedge->nextCoedge = 0;
  coedge->prevCoedge = 0;
  
  return CUBIT_SUCCESS;
}




//-------------------------------------------------------------------------
// Purpose       : remove (and return) all coedges
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/24/02
//-------------------------------------------------------------------------
CubitStatus CompositeLoop::remove_all_coedges( DLIList<CompositeCoEdge*>* list )
{
  while( myCoedge )
  {
    if( list )
      list->append( myCoedge );
    remove( myCoedge );
  }
  
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : get parent bridges
//
// Special Notes : pure virtual in TopologyBridge
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/11/02
//-------------------------------------------------------------------------
void CompositeLoop::get_parents_virt( DLIList<TopologyBridge*>& parents )
  { if( mySurface ) parents.append( mySurface ); }


//-------------------------------------------------------------------------
// Purpose       : get child bridges
//
// Special Notes : pure virtual in TopologyBridge
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/11/02
//-------------------------------------------------------------------------
void CompositeLoop::get_children_virt( DLIList<TopologyBridge*>& children )
{
  CompositeCoEdge* coedge = myCoedge;
  
  if( coedge ) do
  {
    children.append( coedge );
    coedge = next_coedge( coedge );
  } while( coedge != myCoedge );
}


//-------------------------------------------------------------------------
// Purpose       : Get CompositeEngine
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/11/02
//-------------------------------------------------------------------------
GeometryQueryEngine* CompositeLoop::get_geometry_query_engine() const
  { return VirtualQueryEngine::instance(); }

//-------------------------------------------------------------------------
// Purpose       : Attribute functions
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/11/02
//-------------------------------------------------------------------------
void CompositeLoop::append_simple_attribute_virt( CubitSimpleAttrib* )
{ }
void CompositeLoop::remove_simple_attribute_virt( CubitSimpleAttrib* )
{ }
void CompositeLoop::remove_all_simple_attribute_virt()
{ }
CubitStatus CompositeLoop::get_simple_attribute( DLIList<CubitSimpleAttrib*>& )
{ return CUBIT_FAILURE; }
CubitStatus CompositeLoop::get_simple_attribute(
					const CubitString& , DLIList<CubitSimpleAttrib*>& )
{ return CUBIT_FAILURE; }

  
  
//-------------------------------------------------------------------------
// Purpose       : Reverse order and sense of coedges
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/17/02
//-------------------------------------------------------------------------
void CompositeLoop::reverse(bool b_reverse_coedges)
{
  CompositeCoEdge* coedge = myCoedge;
  if( coedge ) do
  {
    CompositeCoEdge* temp = coedge->nextCoedge;
    coedge->nextCoedge = coedge->prevCoedge;
    coedge->prevCoedge = temp;
    if (b_reverse_coedges) coedge->reverse();
    coedge = temp;
  } while( coedge != myCoedge );
}


void CompositeLoop::print_debug_info( const char* line_prefix )
{
  if( line_prefix == 0 )
    line_prefix = "";
    
  char* new_prefix = new char[strlen(line_prefix)+3];
  strcpy( new_prefix, line_prefix );
  strcat( new_prefix, "  " );
  
  PRINT_INFO("%sCompositeLoop @ %p : \n", line_prefix, this );
  CompositeCoEdge* coedge = first_coedge();
  if( !coedge )
    PRINT_INFO("%s  No CoEdges!!\n", line_prefix);
  else do
  {
    coedge->print_debug_info( new_prefix );
    coedge = next_coedge( coedge );
  } while( coedge != first_coedge() );
  
  delete [] new_prefix;
}

