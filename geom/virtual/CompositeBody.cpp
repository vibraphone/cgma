//-------------------------------------------------------------------------
// Filename      : CompositeBody.cpp
//
// Purpose       : Composite of BodySMs
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/11/02
//-------------------------------------------------------------------------

#include "DLIList.hpp"
#include "CompositeBody.hpp"
#include "CompositeLump.hpp"
#include "VirtualQueryEngine.hpp"
#include "CompositeEngine.hpp"

CompositeBody::CompositeBody()
  : firstLump(0)
{
}

CompositeBody::~CompositeBody()
{
  int i;
  
  for( i = 0; i < realBodies.size(); i++ )
    if( realBodies[i]->owner() == this )
      realBodies[i]->owner(0);
  
  while( firstLump )
    remove( firstLump );
}


CompositeLump* CompositeBody::next_lump( CompositeLump* prev ) const
{
  return prev ? prev->nextLump : firstLump;
}

CubitStatus CompositeBody::add( CompositeLump* lump )
{
  if( lump->myBody != 0 )
    return CUBIT_FAILURE;
    
  lump->myBody = this;
  lump->nextLump = firstLump;
  firstLump = lump;
  return CUBIT_SUCCESS;
}

CubitStatus CompositeBody::remove( CompositeLump* lump )
{
  if( lump->myBody != this )
    return CUBIT_FAILURE;
  
  if( firstLump == lump )
  {
    firstLump = lump->nextLump;
  }
  else
  {
    CompositeLump* prev = firstLump; 
    while( prev && prev->nextLump != lump )
      prev = prev->nextLump;
    assert( prev != NULL );
    
    prev->nextLump = lump->nextLump;
  }
  
  lump->myBody = 0;
  lump->nextLump = 0;
  return CUBIT_SUCCESS;
}

CubitStatus CompositeBody::add( BodySM* body )
{
  if( index_of( body ) >= 0  )
    return CUBIT_FAILURE;
    
  if( body->owner() )
    body->owner()->swap_bridge( body, this, false );
  body->owner(this);
  
  realBodies.push( body );
  return CUBIT_SUCCESS;
}

CubitStatus CompositeBody::remove( BodySM* body )
  { return remove_body( index_of( body ) ); }

CubitStatus CompositeBody::remove_body( int index )
{
  if( index < 0 ) return CUBIT_FAILURE;
  
  if( realBodies[index]->owner() == this )
    realBodies[index]->owner(0);
  realBodies.remove(index);
  return CUBIT_SUCCESS;
}
/*
CubitStatus CompositeBody::move( const CubitVector& offset )
{
  int i;
  for (i = 0; i < realBodies.size(); i++)
    if (CUBIT_SUCCESS != realBodies[i]->move( offset ))
      break;
  
  if (i == realBodies.size())
    return CUBIT_SUCCESS;
  
  for (int j = 0; j < i; j++)
    realBodies[j]->move( -offset );
  return CUBIT_FAILURE;
}


CubitStatus CompositeBody::rotate( const CubitVector& axis, double angle )
{
  int i;
  for (i = 0; i < realBodies.size(); i++)
    if (CUBIT_SUCCESS != realBodies[i]->rotate( axis, angle ))
      break;
  
  if (i == realBodies.size())
    return CUBIT_SUCCESS;
  
  for (int j = 0; j < i; j++)
    realBodies[j]->rotate( axis, -angle );
  return CUBIT_FAILURE;
}

CubitStatus CompositeBody::scale( double factor )
{
  int i;
  for (i = 0; i < realBodies.size(); i++)
    if (CUBIT_SUCCESS != realBodies[i]->scale( factor ))
      break;
  
  if (i == realBodies.size())
    return CUBIT_SUCCESS;
  
  for (int j = 0; j < i; j++)
    realBodies[j]->scale( 1.0/factor );
  return CUBIT_FAILURE;
}

CubitStatus CompositeBody::scale( const CubitVector& factors )
{
  int i;
  for (i = 0; i < realBodies.size(); i++)
    if (CUBIT_SUCCESS != realBodies[i]->scale( factors ))
      break;
  
  if (i == realBodies.size())
    return CUBIT_SUCCESS;
  
  const CubitVector unscale( 1.0/factors.x(), 1.0/factors.y(), 1.0/factors.z() );
  for (int j = 0; j < i; j++)
    realBodies[j]->scale( unscale );
  return CUBIT_FAILURE;
}

  
CubitStatus CompositeBody::reflect( const CubitVector& axis )
{
  int i;
  for (i = 0; i < realBodies.size(); i++)
    if (CUBIT_SUCCESS != realBodies[i]->reflect( axis ))
      break;
  
  if (i == realBodies.size())
    return CUBIT_SUCCESS;
  
  for (int j = 0; j < i; j++)
    realBodies[j]->reflect( axis );
  return CUBIT_FAILURE;
}
    

CubitStatus CompositeBody::restore()
  { return CUBIT_FAILURE; }

CubitStatus CompositeBody::reverse()
  { return CUBIT_FAILURE; }
*/
CubitStatus CompositeBody::get_transforms( CubitTransformMatrix& )
  { return CUBIT_FAILURE; }

void CompositeBody::get_parents_virt( DLIList<TopologyBridge*>& )
  { }

void CompositeBody::get_children_virt( DLIList<TopologyBridge*>& children )
{
  for( CompositeLump* lump = firstLump; lump; lump = lump->nextLump )
    children.append( lump );
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
GeometryQueryEngine* CompositeBody::get_geometry_query_engine() const
  { return VirtualQueryEngine::instance(); }

  
CubitStatus CompositeBody::remove_bridge( TopologyBridge* bridge )
{
  int i;
  for (i = realBodies.size() - 1; i >= 0 && realBodies[i] != bridge; --i);
  if (i < 0)
    return CUBIT_FAILURE;
  
  assert( bridge->owner() == this );
  bridge->owner( 0 );
  realBodies.remove( i );
  
  if (realBodies.size() == 0)
    CompositeEngine::instance().notify_deactivated( this );
  
  return CUBIT_SUCCESS;
}
  
  
CubitStatus CompositeBody::swap_bridge( TopologyBridge* old_tb, 
                                        TopologyBridge* new_tb,
                                        bool )
{
  if( new_tb->owner() )
    return CUBIT_FAILURE;
  
  BodySM* new_body = dynamic_cast<BodySM*>(new_tb);
  BodySM* old_body = dynamic_cast<BodySM*>(old_tb);
  int index = realBodies.find( old_body );
  if( index >= 0 && new_body != 0 && realBodies.find(new_body) < 0 )
  {
    if( old_body->owner() == this )
      old_body->owner(0);
    new_body->owner(this);
    realBodies[index] = new_body;
    return CUBIT_SUCCESS;
  }
  
  return CUBIT_FAILURE;
}
  
CubitBoolean CompositeBody::contains_bridge( TopologyBridge* bridge ) const
{
  return index_of(dynamic_cast<BodySM*>(bridge)) < 0 ? CUBIT_FALSE : CUBIT_TRUE;
}

void CompositeBody::notify_reversed( TopologyBridge* )
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
void CompositeBody::append_simple_attribute_virt( const CubitSimpleAttrib& )
{ }
void CompositeBody::remove_simple_attribute_virt( const CubitSimpleAttrib& )
{ }
void CompositeBody::remove_all_simple_attribute_virt()
{ }
CubitStatus CompositeBody::get_simple_attribute( DLIList<CubitSimpleAttrib>& )
{ return CUBIT_FAILURE; }
CubitStatus CompositeBody::get_simple_attribute(
          const CubitString& , DLIList<CubitSimpleAttrib>& )
{ return CUBIT_FAILURE; }

//-------------------------------------------------------------------------
// Purpose       : 
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/10/04
//-------------------------------------------------------------------------
CubitPointContainment CompositeBody::point_containment( const CubitVector& pos )
{
  int inside = 0;
  int boundary = 0;
  
  for (int i = 0; i < realBodies.size(); ++i)
  {
    switch( realBodies[i]->point_containment( pos ) )
    {
      case CUBIT_PNT_BOUNDARY:
        boundary++;
        break;
      case CUBIT_PNT_INSIDE:
        inside++;
        break;
      case CUBIT_PNT_OUTSIDE:
        break;
      default:
        return CUBIT_PNT_UNKNOWN;
    }
  }
  
  if (inside)
    return CUBIT_PNT_INSIDE;
  else if (boundary > 1)
    return CUBIT_PNT_INSIDE;
  else if (boundary)
    return CUBIT_PNT_BOUNDARY;
  else
    return CUBIT_PNT_OUTSIDE;
}

//-------------------------------------------------------------------------
// Purpose       : 
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 05/10/04
//-------------------------------------------------------------------------
CubitStatus CompositeBody::mass_properties( CubitVector& result,
                                            double& volume )
{
  double vol;
  CubitVector centroid;
  result.set( 0.0, 0.0, 0.0 );
  volume = 0;
  
  for (int i = 0; i < realBodies.size(); ++i)
  {
    if (CUBIT_FAILURE == realBodies[i]->mass_properties( centroid, vol ))
      return CUBIT_FAILURE;
    
    result += vol * centroid;
    volume += vol;
  }
  
  if (volume > CUBIT_RESABS)
    result /= volume;
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : Combine
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 06/11/04
//-------------------------------------------------------------------------
void CompositeBody::combine( CompositeBody* other )
{
  int oldsize = realBodies.size();
  realBodies.size( oldsize + other->realBodies.size() );
  for (int i = 0; i < other->realBodies.size(); i++)
  {
    BodySM* bod = other->realBodies[i];
    realBodies[i+oldsize] = bod;
    bod->owner(this);
  }
  other->realBodies.size(0);
}
  
