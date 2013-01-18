#include "DLIList.hpp"

#include "PartitionCoEdge.hpp"
#include "PartitionLoop.hpp"
#include "PartitionSurface.hpp"
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
PartitionLoop::PartitionLoop(  )
  : mySurface(0), firstCoedge( 0 ), nextInSurface( 0 ), numCoedges( 0 )
{
}

//-------------------------------------------------------------------------
// Purpose       : The purpose of this function is to see if a loop is an external
//                  or internal loop of a surface.
//
// Special Notes : 
//
// Creator       : Jonathan Bugman
//
// Creation Date : 9/9/2008
//-------------------------------------------------------------------------
CubitBoolean PartitionLoop::is_external()
{
		  PRINT_ERROR( "This command is not supported with this engine.\n");
          return CUBIT_FAILURE;
}

LoopType PartitionLoop::loop_type()
{
  return LOOP_TYPE_UNKNOWN;
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
PartitionLoop::~PartitionLoop()
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
CubitStatus PartitionLoop::insert_after( PartitionCoEdge* insert,
                                         PartitionCoEdge* after )
{
  if( insert->myLoop != 0 )
    return CUBIT_FAILURE;
  
  if( numCoedges == 0 )
  {
    if( after )
      return CUBIT_FAILURE;
    
    insert->loopNext = insert;
    insert->loopPrev = insert;
    firstCoedge = insert;
  }
  else
  {
    if( !after || after->myLoop != this )
      return CUBIT_FAILURE;
    
    insert->loopPrev = after;
    insert->loopNext = after->loopNext;
    insert->loopNext->loopPrev = insert;
    insert->loopPrev->loopNext = insert;
  }
  
  insert->myLoop = this;
  numCoedges++;
  return CUBIT_SUCCESS;
}
CubitStatus PartitionLoop::insert_before( PartitionCoEdge* insert,
                                          PartitionCoEdge* before )
{
  if( insert->myLoop != 0 )
    return CUBIT_FAILURE;
  
  if( numCoedges == 0 )
  {
    if( before ) 
      return CUBIT_FAILURE;
    
    insert->loopNext = insert;
    insert->loopPrev = insert;
    firstCoedge = insert;
  }
  else
  {
    if( !before || before->myLoop != this )
      return CUBIT_FAILURE;
    
    insert->loopNext = before;
    insert->loopPrev = before->loopPrev;
    insert->loopNext->loopPrev = insert;
    insert->loopPrev->loopNext = insert;
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
CubitStatus PartitionLoop::remove( PartitionCoEdge* coedge )
{
  if( coedge->myLoop != this )
    return CUBIT_FAILURE;
  
  if( numCoedges == 1 )
  {
    assert( coedge->loopNext == coedge &&
            coedge->loopPrev == coedge );
    firstCoedge = 0;
  }
  else
  {
    coedge->loopNext->loopPrev = coedge->loopPrev;
    coedge->loopPrev->loopNext = coedge->loopNext;
    if( firstCoedge == coedge )
      firstCoedge = coedge->loopNext;
  }
  
  numCoedges--;
  coedge->myLoop = 0;
  coedge->loopNext = 0;
  coedge->loopPrev = 0;
  
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
CubitStatus PartitionLoop::remove_all_coedges( DLIList<PartitionCoEdge*>* list )
{
  while( firstCoedge )
  {
    if( list )
      list->append( firstCoedge );
    remove( firstCoedge );
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
void PartitionLoop::get_parents_virt( DLIList<TopologyBridge*>& parents )
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
void PartitionLoop::get_children_virt( DLIList<TopologyBridge*>& children )
{
  PartitionCoEdge* coedge = firstCoedge;
  
  if( coedge ) do
  {
    children.append( coedge );
    coedge = next_coedge( coedge );
  } while( coedge != firstCoedge );
}


//-------------------------------------------------------------------------
// Purpose       : get tb layer number
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/03/03
//-------------------------------------------------------------------------
int PartitionLoop::layer() const
  { return get_surface()->layer(); }

//-------------------------------------------------------------------------
// Purpose       : Get PartitionEngine
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/11/02
//-------------------------------------------------------------------------
GeometryQueryEngine* PartitionLoop::get_geometry_query_engine() const
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
void PartitionLoop::append_simple_attribute_virt( const CubitSimpleAttrib& )
{ }
void PartitionLoop::remove_simple_attribute_virt( const CubitSimpleAttrib& )
{ }
void PartitionLoop::remove_all_simple_attribute_virt()
{ }
CubitStatus PartitionLoop::get_simple_attribute( DLIList<CubitSimpleAttrib>& )
{ return CUBIT_FAILURE; }
CubitStatus PartitionLoop::get_simple_attribute( const CubitString&,
                                                 DLIList<CubitSimpleAttrib>& )
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
void PartitionLoop::reverse()
{
  PartitionCoEdge* coedge = firstCoedge;
  if( coedge ) do
  {
    PartitionCoEdge* temp = coedge->loopNext;
    coedge->loopNext = coedge->loopPrev;
    coedge->loopPrev = temp;
    coedge->reverse_sense();
    coedge = temp;
  } while( coedge != firstCoedge );
}


void PartitionLoop::print_debug_info( const char* /*line_prefix*/ )
{
/*
  if( line_prefix == 0 )
    line_prefix = "";
    
  char* new_prefix = new char[strlen(line_prefix)+3];
  strcpy( new_prefix, line_prefix );
  strcat( new_prefix, "  " );
  
  PRINT_INFO("%sPartitionLoop @ %p : \n", line_prefix, this );
  PartitionCoEdge* coedge = first_coedge();
  if( !coedge )
    PRINT_INFO("%s  No CoEdges!!\n");
  else do
  {
    coedge->print_debug_info( new_prefix );
    coedge = next_coedge( coedge );
  } while( coedge != first_coedge() );
  
  delete [] new_prefix;
*/
}

