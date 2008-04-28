//-------------------------------------------------------------------------
// Filename      : CompositeBody.hpp
//
// Purpose       : Composite of BodySMs
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/11/02
//-------------------------------------------------------------------------

#ifndef COMPOSITE_BODY_HPP
#define COMPOSITE_BODY_HPP

#include "VGDefines.h"
#include "BodySM.hpp"
#include "TBOwner.hpp"
#include "VGArray.hpp"

class CompositeLump;

class CompositeBody: public BodySM, public TBOwner
{
  public:
  
    CompositeBody();
    ~CompositeBody();
    
#ifdef BOYD15
    int num_lumps() const;
#endif
    CompositeLump* next_lump( CompositeLump* after_this = 0 ) const;
    
    CubitStatus add( CompositeLump* lump );
    CubitStatus remove( CompositeLump* lump );
    
    int num_bodies() const;
    BodySM* get_body( int index ) const;
    int index_of( BodySM* body ) const;
    
    CubitStatus add( BodySM* body );
    CubitStatus remove( BodySM* body );
    CubitStatus remove_body( int index );
/*    
    CubitStatus move( const CubitVector& offset );
    CubitStatus rotate( const CubitVector& axis, double degrees );
    CubitStatus scale( double factor );
    CubitStatus scale( const CubitVector& factors );
    CubitStatus reflect( const CubitVector& axis );
    CubitStatus restore();
    CubitStatus reverse();
*/
    CubitStatus get_transforms( CubitTransformMatrix& tfm );
    
    void get_parents_virt( DLIList<TopologyBridge*>& parents );
    void get_children_virt( DLIList<TopologyBridge*>& children );
    int layer() const { return COMPOSITE_LAYER; }
    GeometryQueryEngine* get_geometry_query_engine() const;

    virtual void append_simple_attribute_virt( CubitSimpleAttrib* simple_attrib_ptr );
    virtual void remove_simple_attribute_virt( CubitSimpleAttrib* simple_attrib_ptr );
    virtual void remove_all_simple_attribute_virt();
    virtual CubitStatus get_simple_attribute( DLIList<CubitSimpleAttrib*>& attrib_list );
    virtual CubitStatus get_simple_attribute( const CubitString& name,
                                    DLIList<CubitSimpleAttrib*>& attrib_list );
 
    CubitStatus remove_bridge( TopologyBridge* bridge );
    CubitStatus swap_bridge( TopologyBridge* old_tb, 
                             TopologyBridge* new_tb, 
                             bool reversed );
    CubitBoolean contains_bridge( TopologyBridge* bridge ) const;
    void notify_reversed( TopologyBridge* bridge );
    
    CubitPointContainment point_containment( const CubitVector& pos );
    
    CubitStatus mass_properties( CubitVector& centroid, double& volume );
    
    void combine( CompositeBody* other );
    
    CubitBoolean is_sheet_body(){return CUBIT_FALSE;}
  private:
  
    CompositeLump* firstLump;
    VGArray<BodySM*> realBodies;
    
};

inline int CompositeBody::num_bodies() const
  { return realBodies.size(); }

inline BodySM* CompositeBody::get_body( int index ) const
  { return realBodies[index]; }

inline int CompositeBody::index_of( BodySM* body ) const
  { return realBodies.find( body ); }



#endif
