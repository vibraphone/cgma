//-------------------------------------------------------------------------
// Filename      : CompositeShell.cpp
//
// Purpose       : ShellSM used in composite TopologyBridge graph
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/11/02
//-------------------------------------------------------------------------

#include "CompositeShell.hpp"
#include "CompositeCoSurf.hpp"
#include "CompositeSurface.hpp"
#include "VirtualQueryEngine.hpp"
#include "CompositeLump.hpp"

CompositeShell::CompositeShell()
  : myLump(0), lumpNext(0), firstCoSurf(0)
{
}

CompositeShell::~CompositeShell()
{
//  if( myLump ) 
//    myLump->remove( this );
  
  while( firstCoSurf )
  {
    CompositeCoSurf* cosurf = firstCoSurf;
    remove( cosurf );
    if( cosurf->get_surface() )
      cosurf->get_surface()->remove( cosurf );
    delete cosurf;
  }
  
  assert( !myLump );
}

CompositeCoSurf* CompositeShell::next_co_surf( CompositeCoSurf* prev ) const
{ 
  if (0 == prev)
    return firstCoSurf;
  else if (this == prev->myShell)
    return prev->shellNext;
  else
    return 0;
}

CubitStatus CompositeShell::add( CompositeCoSurf* cosurf )
{
  if( cosurf->myShell )
    return CUBIT_FAILURE;
  
  cosurf->shellNext = firstCoSurf;
  firstCoSurf = cosurf;
  cosurf->myShell = this;

  return CUBIT_SUCCESS;
}

CubitStatus CompositeShell::remove( CompositeCoSurf* cosurf )
{
  if( cosurf->myShell != this )
    return CUBIT_FAILURE;
  
  if( firstCoSurf == cosurf )
    firstCoSurf = cosurf->shellNext;
  else 
  {
    CompositeCoSurf* prev = firstCoSurf;
    while( prev && prev->shellNext != cosurf )
      prev = prev->shellNext;
    assert( prev != 0 );
    prev->shellNext = cosurf->shellNext;
  }
  
  cosurf->myShell = 0;
  cosurf->shellNext = 0;
  return CUBIT_SUCCESS;
}


CompositeCoSurf* CompositeShell::add( CompositeSurface* surface,
                                      CubitSense sense )
{
  CompositeCoSurf* cosurf = new CompositeCoSurf( sense );
  if( !surface->add( cosurf ) || !add( cosurf ) )
  {
    delete cosurf;
    cosurf = 0;
  }
  return cosurf;
}

CompositeCoSurf* CompositeShell::find_first( const CompositeSurface* surf ) const
{
  CompositeCoSurf* cosurf = firstCoSurf;
  while( cosurf && cosurf->get_surface() != surf )
    cosurf = cosurf->shellNext;
  return cosurf;
}

CompositeCoSurf* CompositeShell::find_next( const CompositeCoSurf* prev ) const
{
  CompositeCoSurf* cosurf = prev->shellNext;
  while( cosurf && cosurf->get_surface() != prev->get_surface() )
    cosurf = cosurf->shellNext;
  return cosurf;
}

CubitSense CompositeShell::find_sense( const CompositeSurface* surf ) const
{
  CompositeCoSurf* first = find_first( surf );
  if( !first )
    return CUBIT_UNKNOWN;
  
  CompositeCoSurf* next = first;
  while( (next = find_next( next )) )
    if( next->sense() != first->sense() )
      return CUBIT_UNKNOWN;
  
  return first->sense();
}


  
  


void CompositeShell::get_parents_virt( DLIList<TopologyBridge*>& parents )
  { if( myLump ) parents.append( static_cast<TopologyBridge*>(myLump) ); }

void CompositeShell::get_children_virt( DLIList<TopologyBridge*>& children )
{ 
  for( CompositeCoSurf* c = firstCoSurf; c; c = c->next_in_shell() )
    children.append( c->get_surface() );
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
GeometryQueryEngine* CompositeShell::get_geometry_query_engine() const
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
void CompositeShell::append_simple_attribute_virt( CubitSimpleAttrib* )
{ }
void CompositeShell::remove_simple_attribute_virt( CubitSimpleAttrib* )
{ }
void CompositeShell::remove_all_simple_attribute_virt()
{ }
CubitStatus CompositeShell::get_simple_attribute( DLIList<CubitSimpleAttrib*>& )
{ return CUBIT_FAILURE; }
CubitStatus CompositeShell::get_simple_attribute(
					const CubitString& , DLIList<CubitSimpleAttrib*>&  )
{ return CUBIT_FAILURE; }


//-------------------------------------------------------------------------
// Purpose       : Debug output
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/16/04
//-------------------------------------------------------------------------
void CompositeShell::print_debug_info( const char* prefix )
{
  if( prefix == 0 ) prefix = "";
  PRINT_INFO("%sCompositeShell @ %p : \n", prefix, this );
    
  char* new_prefix = new char[strlen(prefix)+3];
  strcpy( new_prefix, prefix );
  strcat( new_prefix, "  " );
  
  CompositeCoSurf* cosurf = first_co_surf();
  if( !cosurf )
    PRINT_INFO("%sNo Surfaces!!\n", prefix);
  else for (; cosurf; cosurf = next_co_surf( cosurf ))
    cosurf->print_debug_info( new_prefix );
  
  delete [] new_prefix;
}
