//-------------------------------------------------------------------------
// Filename      : CompositeCoEdge.cpp
//
// Purpose       : Combined set of CoEdgeSMs
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/11/02
//-------------------------------------------------------------------------

#include "CompositeCoEdge.hpp"
#include "CompositeCurve.hpp"
#include "CompositeLoop.hpp"
#include "CompositeEngine.hpp"
#include "VirtualQueryEngine.hpp"

// for printing debug info
#include "RefEdge.hpp"

//-------------------------------------------------------------------------
// Purpose       : Constructor
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/11/02
//-------------------------------------------------------------------------
CompositeCoEdge::CompositeCoEdge( CoEdgeSM* coedge )
  : myLoop(0), 
    nextCoedge(0), 
    prevCoedge(0), 
    myCurve(0), 
    nextOnCurve(0)
{
  mySense = coedge->sense();
  coedgeSet.push( coedge );
  if( coedge->owner() )
    coedge->owner()->swap_bridge( coedge, this, false );
  coedge->owner(this);
}

CompositeCoEdge::CompositeCoEdge( CompositeCurve* point_curve )
  : mySense(CUBIT_FORWARD),
    myLoop(0), 
    nextCoedge(0), 
    prevCoedge(0), 
    myCurve(0), 
    nextOnCurve(0)
{
  assert(point_curve->num_curves() == 0);
  CubitStatus stat = point_curve->add(this);
  assert(stat);
}

CompositeCoEdge::CompositeCoEdge()
  : mySense(CUBIT_UNKNOWN),
    myLoop(0), 
    nextCoedge( 0 ),
    prevCoedge( 0 ), 
    myCurve(0),
    nextOnCurve(0)
  {}

//-------------------------------------------------------------------------
// Purpose       : Destructor
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/11/02
//-------------------------------------------------------------------------
CompositeCoEdge::~CompositeCoEdge()
{
  if( myLoop )
    myLoop->remove(this);
  if( myCurve )
    myCurve->remove(this);

  assert( !prevCoedge && !nextCoedge && !nextOnCurve );

  for( int i = 0; i < coedgeSet.size(); i++ )
    if( coedgeSet[i]->owner() == this )
      coedgeSet[i]->owner(0);
}
/*
//-------------------------------------------------------------------------
// Purpose       : Add an underlying coedge
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/11/02
//-------------------------------------------------------------------------
CubitStatus CompositeCoEdge::append( CoEdgeSM* coedge_ptr )
{
  if( index_of( coedge_ptr ) >= 0 || coedge_ptr->owner() ) 
    return CUBIT_FAILURE;
  
  coedge_ptr->owner( this );
  coedgeSet.push( coedge_ptr );
  
  return CUBIT_SUCCESS;
}
*/

//-------------------------------------------------------------------------
// Purpose       : find the index of the coedge that owns the passed curve
//
// Special Notes : returns -1 if not found
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/11/02
//-------------------------------------------------------------------------
int CompositeCoEdge::index_of( Curve* ptr ) const
{
  int i;
  DLIList<TopologyBridge*> curve_list;
  TopologyBridge* curve = ptr;
  for( i = coedgeSet.size() - 1; i > 0; i-- )
  {
    curve_list.clean_out();
    coedgeSet[i]->get_children( curve_list );
    if( curve_list.is_in_list( curve ) )
      break;
  }
  return i;
}

//-------------------------------------------------------------------------
// Purpose       : Split this CompositeCoEdge into two at the specified index
//
// Special Notes : new/other gets CoEdge at the passed index
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/11/02
//-------------------------------------------------------------------------
CompositeCoEdge* CompositeCoEdge::split( int index )
{
  if( index < 0 || index >= coedgeSet.size() )
    return 0;
  
  ++index;
  CompositeCoEdge* new_cce = new CompositeCoEdge();
  new_cce->mySense = mySense;
  
  int new_cce_count = coedgeSet.size() - index;
  new_cce->coedgeSet.size( new_cce_count );
  
  for( int i = 0; i < new_cce_count; i++ )
  {
    new_cce->coedgeSet[i] = coedgeSet[i+index];
    new_cce->coedgeSet[i]->owner( new_cce );
  }
  coedgeSet.size( index );
  
  if( myLoop )
  {
    CubitStatus s;
    if( mySense == CUBIT_FORWARD )
      s = myLoop->insert_after( new_cce, this );
    else
      s = myLoop->insert_before( new_cce, this );
    assert( s );
  }
  
  new_cce->mySense = mySense;
  
  return new_cce;
}

CubitStatus CompositeCoEdge::combine( CompositeCoEdge* dead, bool prepend )
{
  int insert;
  if ( prepend )
  {
    insert = 0;
    coedgeSet.size_end( coedgeSet.size() + dead->coedgeSet.size() );
  }
  else
  {
    insert = coedgeSet.size();
    coedgeSet.size( coedgeSet.size() + dead->coedgeSet.size() );
  }
  
  for( int i = 0; i < dead->coedgeSet.size(); i++ )
  {
    CoEdgeSM* coedge = dead->coedgeSet[i];
    assert( coedge->owner() == dead );
    coedge->owner(this);
    coedgeSet[insert++] = coedge;
  }
  dead->coedgeSet.size(0);
  return CUBIT_SUCCESS;
}

/*
//-------------------------------------------------------------------------
// Purpose       : dequeue underlying coedge
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/11/02
//-------------------------------------------------------------------------
CoEdgeSM* CompositeCoEdge::remove_first()
{
  CoEdgeSM* result = 0;
  if( coedgeSet.size() > 0 )
  {
    result = coedgeSet[0];
    coedgeSet.remove(0);
  }
  return result;
}

//-------------------------------------------------------------------------
// Purpose       : pop underlying coedge
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/11/02
//-------------------------------------------------------------------------
CoEdgeSM* CompositeCoEdge::remove_last()
{
  CoEdgeSM* result = 0;
  if( coedgeSet.size() > 0 )
  {
    result = coedgeSet.pop()
    result->owner(0);
  }
  return result;
}
*/
//-------------------------------------------------------------------------
// Purpose       : get parents (pure virtual in TopologyBridge)
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/11/02
//-------------------------------------------------------------------------
void CompositeCoEdge::get_parents_virt( DLIList<TopologyBridge*>& parents )
{ 
  LoopSM* result = get_parent_loop();
  if( result ) 
  {
    parents.append( result ); 
  }
}
  
//-------------------------------------------------------------------------
// Purpose       : get children (pure virtual in TopologyBridge)
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/11/02
//-------------------------------------------------------------------------
void CompositeCoEdge::get_children_virt( DLIList<TopologyBridge*>& children )
{ 
  if( myCurve ) 
  {
    children.append( myCurve->primary_stitched_curve() ); 
  }
/*
  else if( num_coedges() )
  {
    DLIList<TopologyBridge*> coedge_children;
    coedge(0)->get_children( coedge_children );
    assert( coedge_children.size() == 1 );
    children.append( coedge_children.get() );
  }
*/
}
    
/*
//-------------------------------------------------------------------------
// Purpose       : Set child curve
//
// Special Notes : If curve is composite, update curve link to this
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/11/02
//-------------------------------------------------------------------------
void CompositeCoEdge::curve( CompositeCurve* curve_ptr )
{
  if( myCurve )
  {
    myCurve->remove(this);
    assert( !myCurve );
  }  
  
  if( curve_ptr )
  {
    curve_ptr->add(this);
    assert( myCurve == curve_ptr );
  }
}

//-------------------------------------------------------------------------
// Purpose       : set loop pointer
//
// Special Notes : if loop is composite, add this to loop instead
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/11/02
//-------------------------------------------------------------------------
void CompositeCoEdge::loop( CompositeLoop* loop_ptr )
{
  if( myLoop )
    myLoop->remove( this );
  
  if( loop_ptr )
    loop_ptr->insert_before( loop_ptr->first_coedge(), this );
  
  myLoop = loop_ptr;
}
*/
//-------------------------------------------------------------------------
// Purpose       : Get parent loop no higher than composite level
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/05/02
//-------------------------------------------------------------------------
LoopSM* CompositeCoEdge::get_parent_loop()
{
  LoopSM* result = get_loop();
  
  if( !result && num_coedges() )
  {
    DLIList<TopologyBridge*> parents(1);
    get_coedge(0)->get_parents_virt( parents );
    assert( parents.size() == 1 );
    result = dynamic_cast<LoopSM*>(parents.get());
  }
  
  return result;
}

//-------------------------------------------------------------------------
// Purpose       : remove an underlying bridge
//
// Special Notes : pure virtual in TBOwner
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/11/02
//-------------------------------------------------------------------------
CubitStatus CompositeCoEdge::remove_bridge( TopologyBridge* bridge )
{
  int index;
  for( index = coedgeSet.size() - 1; index >= 0; index-- )
    if( coedgeSet[index] == bridge )
      break;
  if( index < 0 )
    return CUBIT_FAILURE;
  
  coedgeSet.remove( index );
  bridge->owner(0);
/*  
  if ( coedgeSet.size() > 0 )
    return CUBIT_SUCCESS;
  
  if ( get_curve() )
    get_curve()->remove(this);

  if ( get_loop() )
  {
    CompositeLoop* loop = get_loop();
    loop->remove(this);
    if ( loop->first_coedge() == 0 )
      delete loop;
  }
  
  delete this;  
*/    
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : exchange one underlying coedge for another
//
// Special Notes : pure virtual in TBOwner
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/11/02
//-------------------------------------------------------------------------
CubitStatus CompositeCoEdge::swap_bridge( TopologyBridge* old_tb,
                                          TopologyBridge* new_tb,
                                          bool )
{
  CoEdgeSM* old_coedge = dynamic_cast<CoEdgeSM*>(old_tb);
  CoEdgeSM* new_coedge = dynamic_cast<CoEdgeSM*>(new_tb);
  
  int index = index_of( old_coedge );
  if( index < 0 || !new_coedge || index_of(new_coedge) >= 0 )
    return CUBIT_FAILURE;
  
  coedgeSet[index] = new_coedge;

  old_tb->owner(0);
  if( new_tb->owner() )
    new_tb->owner()->remove_bridge( new_tb );
  new_tb->owner(this);
  
  return CUBIT_SUCCESS;
}

//-------------------------------------------------------------------------
// Purpose       : see if we are the owner of the passed TB
//
// Special Notes : pure virtual in TBOwner
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/11/02
//-------------------------------------------------------------------------
CubitBoolean CompositeCoEdge::contains_bridge( TopologyBridge* bridge ) const
{
  CompositeCoEdge* coedge = dynamic_cast<CompositeCoEdge*>(bridge);
  return (index_of(coedge) < 0) ? CUBIT_FALSE : CUBIT_TRUE;
}

void CompositeCoEdge::notify_reversed( TopologyBridge* )
  {}
  

//-------------------------------------------------------------------------
// Purpose       : Get CompositeEngine
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/11/02
//-------------------------------------------------------------------------
GeometryQueryEngine* CompositeCoEdge::get_geometry_query_engine() const
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
void CompositeCoEdge::append_simple_attribute_virt( const CubitSimpleAttrib& )
{ }
void CompositeCoEdge::remove_simple_attribute_virt( const CubitSimpleAttrib& )
{ }
void CompositeCoEdge::remove_all_simple_attribute_virt()
{ }
CubitStatus CompositeCoEdge::get_simple_attribute( DLIList<CubitSimpleAttrib>& )
{ return CUBIT_FAILURE; }
CubitStatus CompositeCoEdge::get_simple_attribute( const CubitString& ,
                                                DLIList<CubitSimpleAttrib>& )
{ return CUBIT_FAILURE; }




void CompositeCoEdge::reverse()
{
  switch( mySense ) {
    case CUBIT_FORWARD:
      mySense = CUBIT_REVERSED;
      break;
    case CUBIT_REVERSED:
      mySense = CUBIT_FORWARD;
      break;
    default:
      mySense = CUBIT_UNKNOWN;
  }
  
  int half = coedgeSet.size() / 2;
  for( int i = 0; i < half; i++ )
  {
    int j = coedgeSet.size() - i - 1;
    CoEdgeSM* tmp = coedgeSet[i];
    coedgeSet[i] = coedgeSet[j];
    coedgeSet[j] = tmp;
  }
    
}

//-------------------------------------------------------------------------
// Purpose       : Get start and end points reversed if sense is reversed
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 03/17/02
//-------------------------------------------------------------------------
CompositePoint* CompositeCoEdge::start_point()
{
  return mySense == CUBIT_FORWARD 
       ? myCurve->start_point() 
       : myCurve->end_point();
}
CompositePoint* CompositeCoEdge::end_point()
{
  return mySense == CUBIT_FORWARD 
       ? myCurve->end_point() 
       : myCurve->start_point();
}


CubitStatus CompositeCoEdge::remove_coedge( int index )
{
  if( index < 0 || index >= coedgeSet.size() )
    return CUBIT_FAILURE;
  
  coedgeSet[index]->owner(0);
  coedgeSet.remove( index );

  return CUBIT_SUCCESS;
}


CubitStatus CompositeCoEdge::insert_coedge( int index, CoEdgeSM* coedge )
{
  if( index < 0 || index > coedgeSet.size() )
    return CUBIT_FAILURE;
  
  coedgeSet.insert(coedge, index);
  coedge->owner(this);
  return CUBIT_SUCCESS;
}

void CompositeCoEdge::print_debug_info( const char* prefix, bool brief )
{
  if( prefix == 0 ) prefix = "";
  
  const char* sense = mySense == CUBIT_FORWARD ? "Forward" :
                      mySense == CUBIT_REVERSED ? "Reverse" : "UNKNOWN";
                      
  PRINT_INFO("%sCompCoEdge %p %s ", prefix, this, sense );
  if ( num_coedges() == 1 )
    PRINT_INFO("%s %p ", fix_type_name(typeid(*get_coedge(0)).name()),
      get_coedge(0));
  else
    PRINT_INFO("%d coedges ", num_coedges() );
    
  if( !myCurve )
    PRINT_INFO("NULL CURVE\n");
  else if( brief )
#ifdef TOPOLOGY_BRIDGE_IDS
    PRINT_INFO("curve %d\n", myCurve->get_id() );
#else
    PRINT_INFO("curve %p\n", myCurve );
#endif
  else
    { PRINT_INFO("\n  ");  myCurve->print_debug_info(prefix, true); }

/*  
  if( coedgeSet.size() == 0 )
    PRINT_INFO(" No CoEdgeSMs!\n");
  else if( coedgeSet.size() == 1 )
    PRINT_INFO(" CoEdgeSM=%p\n", coedgeSet[0] );
  else
  {
    PRINT_INFO("\n");
    for( int i = 0; i < coedgeSet.size(); i++ )
      PRINT_INFO("%s  CoEdgeSM[%d] = %p\n", prefix, i, coedgeSet[i] );
  }
*/
}
