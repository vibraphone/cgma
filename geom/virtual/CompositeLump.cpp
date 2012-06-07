//-------------------------------------------------------------------------
// Filename      : CompositeLump.cpp
//
// Purpose       : Combine Lumps
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/11/02
//-------------------------------------------------------------------------

#include "CompositeLump.hpp"
#include "CompositeShell.hpp"
#include "CompositeBody.hpp"
#include "VirtualQueryEngine.hpp"
#include "CompositeEngine.hpp"

//-------------------------------------------------------------------------
// Purpose       : Public constructor - replace real Lump with Composite
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/20/02
//-------------------------------------------------------------------------
CompositeLump::CompositeLump( Lump* lump )
  : myBody(0), nextLump(0), firstShell(0), hiddenSet(0)
{ 
  assert( lump != NULL );
  compGeom = new CompositeGeom(1);
  compGeom->append( lump, CUBIT_FORWARD ); 
  if( lump->owner() )
    lump->owner()->swap_bridge( lump, this, false );
  lump->owner( this );
}

//-------------------------------------------------------------------------
// Purpose       : Internal constructor for splitting Composites
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/20/02
//-------------------------------------------------------------------------
CompositeLump::CompositeLump( CompositeGeom* geometry )
  : myBody(0),
    nextLump(0),
    compGeom( geometry ),
    firstShell(0),
    hiddenSet(0)
{
  assert( geometry != NULL );
  for( int i = 0; i < compGeom->num_entities(); i++ )
  {
    GeometryEntity* entity = compGeom->entity(i);
    assert( !entity->owner() );
    entity->owner(this);
  }
}

//-------------------------------------------------------------------------
// Purpose       : Destructor
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 07/20/02
//-------------------------------------------------------------------------
CompositeLump::~CompositeLump()
{
  int i;
  assert( !myBody );
  
  while( firstShell )
    remove( firstShell );
    
  for( i = 0; i < num_lumps(); i++ )
    if( get_lump(i)->owner() == this )
      get_lump(i)->owner(0);
      
  delete hiddenSet;
  delete compGeom;
}

CubitStatus CompositeLump::add( CompositeShell* shell )
{
  if( shell->myLump )
    return CUBIT_FAILURE;
  
  shell->lumpNext = firstShell;
  firstShell = shell;
  shell->myLump = this;
  return CUBIT_SUCCESS;
}

CubitStatus CompositeLump::remove( CompositeShell* shell )
{
  if( shell->myLump != this )
    return CUBIT_FAILURE;
  
  if( firstShell == shell )
  {
    firstShell = firstShell->lumpNext;
  }
  else
  {
    CompositeShell *prev = firstShell, *next = firstShell->lumpNext;
    while( next != shell )
    {
      assert( next != NULL );
      prev = next;
      next = next->lumpNext;
    }
    
    prev->lumpNext = shell->lumpNext;
  }
  
  shell->lumpNext = 0;
  shell->myLump = 0;
  
  return CUBIT_SUCCESS;
}


CubitBox CompositeLump::bounding_box() const
  { return compGeom->bounding_box(); }
  
double CompositeLump::measure()
  { return compGeom->measure( compGeom->num_entities() -1 ); }

void CompositeLump::get_parents_virt( DLIList<TopologyBridge*>& parents )
  { if( myBody ) parents.append( myBody ); }

void CompositeLump::get_children_virt( DLIList<TopologyBridge*>& children )
{
  for( CompositeShell* shell = firstShell; shell; shell = shell->lumpNext )
    children.append( shell );
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
GeometryQueryEngine* CompositeLump::get_geometry_query_engine() const
  { return VirtualQueryEngine::instance(); }


CubitStatus CompositeLump::add( Lump* lump )
{ 
  if( !lump->owner() && compGeom->append(lump,CUBIT_FORWARD) )
  {
    lump->owner(this);
    return CUBIT_SUCCESS;
  }
  return CUBIT_FAILURE;
}    

CubitStatus CompositeLump::remove( Lump* lump )
  { return remove_lump( index_of( lump ) ); }

CubitStatus CompositeLump::remove_lump( int index )
{
  if( index < 0 || index >= num_lumps() )
    return CUBIT_FAILURE;
  
  TopologyBridge* lump = compGeom->entity(index);
  if( compGeom->remove( index, false ) )
  {
    if( lump->owner() == this )
      lump->owner(0);
    return CUBIT_SUCCESS;
  }
  return CUBIT_FAILURE;
}
  
CubitStatus CompositeLump::remove_bridge( TopologyBridge* bridge )
{
  int i = compGeom->index_of(bridge);
  if (i < 0)
    return CUBIT_FAILURE;
  
  assert (bridge->owner() == this);
  bridge->owner(0);
  if (!compGeom->remove(i, true))
    return CUBIT_FAILURE;
  
  if (compGeom->num_entities() == 0)
    CompositeEngine::instance().notify_deactivated(this);
  
  return CUBIT_SUCCESS;
}
  
  
CubitStatus CompositeLump::swap_bridge( TopologyBridge* old_tb, 
                                        TopologyBridge* new_tb,
                                        bool )
{
  if( new_tb->owner() )
    return CUBIT_FAILURE;
  
  int index = compGeom->index_of( old_tb );
  GeometryEntity* ge = dynamic_cast<GeometryEntity*>(new_tb);
  if( index >= 0 && ge != 0 )
  {
    if( old_tb->owner() == this )
      old_tb->owner(0);
    new_tb->owner(this);
    return compGeom->swap( index, ge );
  }
  
  return CUBIT_FAILURE;
}
  
CubitBoolean CompositeLump::contains_bridge( TopologyBridge* bridge ) const
{
  return compGeom->index_of(bridge) < 0 ? CUBIT_FALSE : CUBIT_TRUE;
}

void CompositeLump::notify_reversed( TopologyBridge* )
  { assert(0); }


//-------------------------------------------------------------------------
// Purpose       : Attribute functions
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/11/02
//-------------------------------------------------------------------------
void CompositeLump::append_simple_attribute_virt( CubitSimpleAttrib* csa )
{ compGeom->add_attribute(csa); }
void CompositeLump::remove_simple_attribute_virt( CubitSimpleAttrib* csa )
{ compGeom->rem_attribute(csa); }
void CompositeLump::remove_all_simple_attribute_virt()
{ compGeom->rem_all_attributes(); }
CubitStatus CompositeLump::get_simple_attribute( DLIList<CubitSimpleAttrib*>& list )
{ compGeom->get_attributes( list ); return CUBIT_SUCCESS; }
CubitStatus CompositeLump::get_simple_attribute(
					const CubitString& name, DLIList<CubitSimpleAttrib*>& attrib_list )
{
  compGeom->get_attributes( name.c_str(), attrib_list );
  return CUBIT_SUCCESS;
}


//-------------------------------------------------------------------------
// Purpose       : Split this CompositeLump into two.
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/11/04
//-------------------------------------------------------------------------
CompositeLump* CompositeLump::split( VGArray<int>& indices_to_move )
{
  int i;
  
  for( i = 0; i < indices_to_move.size(); i++ )
    if( indices_to_move[i] < 0 || indices_to_move[i] >= num_lumps() )
      return 0;
  
  CompositeGeom* new_geom = compGeom->split( indices_to_move );
  if( !new_geom )
    return 0;
  
  for( i = 0; i < new_geom->num_entities(); i++ )
    new_geom->entity(i)->owner( 0 );
    
  return new CompositeLump( new_geom );
}

//-------------------------------------------------------------------------
// Purpose       : Combine composite volumes
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/11/04
//-------------------------------------------------------------------------
CubitStatus CompositeLump::combine( CompositeLump* dead_lump )
{
  int old_size = compGeom->num_entities();
  compGeom->merge( *(dead_lump->compGeom) );
  if( dead_lump->hiddenSet != 0 )
    hidden_entities().merge( dead_lump->hiddenSet );
  for( int i = old_size; i < compGeom->num_entities(); i++ )
  {
    TopologyBridge* bridge = compGeom->entity(i);
    assert( bridge->owner() == dead_lump );
    bridge->owner( this );
  }
  
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Debug output
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 09/16/04
//-------------------------------------------------------------------------
void CompositeLump::print_debug_info( const char* prefix, bool brief )
{
  if( prefix == 0 ) prefix = "";
  CompositeShell* shell = 0;
  if (brief)
  {
    int count = 0;
    while ((shell = next_shell(shell))  != NULL ) ++count;
    PRINT_INFO( "%sCompositeLump %p : %d shells, %d lumps.\n",
                prefix, static_cast<void*>(this), count, num_lumps() );
    return;
  }
  
  PRINT_INFO("%sCompositeLump %p:\n", prefix, static_cast<void*>(this) );
  
  char* new_prefix = new char[strlen(prefix)+3];
  strcpy( new_prefix, prefix );
  strcat( new_prefix, "  " );
  if (hiddenSet) hiddenSet->print_debug_info( new_prefix );
  else PRINT_INFO("%sNo Hidden Entities.\n", new_prefix );
  while ((shell = next_shell( shell )) != NULL )
    shell->print_debug_info( new_prefix );
  delete [] new_prefix;
}

  
void CompositeLump::get_hidden_surfaces( DLIList<Surface*>& surfs )
{
  if (hiddenSet)
    hiddenSet->hidden_surfaces( surfs );
}

CubitStatus CompositeLump::mass_properties( CubitVector &centroid, double &volume )
{
  PRINT_ERROR("CompositeLump::mass_properties is not implemented\n");
  centroid.set(0,0,0);
  volume = 0;
  return CUBIT_FAILURE;
}
