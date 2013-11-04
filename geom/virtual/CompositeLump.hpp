//-------------------------------------------------------------------------
// Filename      : CompositeLump.hpp
//
// Purpose       : Combine Lumps
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/11/02
//-------------------------------------------------------------------------

#ifndef COMPOSITE_LUMP_HPP
#define COMPOSITE_LUMP_HPP

#include "VGDefines.h"
#include "Lump.hpp"
#include "CompositeGeom.hpp"
#include "TBOwner.hpp"
#include "HiddenEntitySet.hpp"
#include "CompositeShell.hpp"

class CompositeBody;
class HiddenEntitySet;

class CompositeLump : public Lump, public TBOwner
{
friend class CompositeBody;
public:
  
  CompositeLump( Lump* real_lump );
  CompositeLump( CompositeGeom* geom );
  virtual ~CompositeLump();
  
  int num_lumps() const;
  Lump* get_lump( int index ) const;
  int index_of( Lump* ) const;
  void update();
 
  CubitStatus add( Lump* lump );
  CubitStatus remove( Lump* lump );
  CubitStatus remove_lump( int index );
  
  HiddenEntitySet& hidden_entities();
  bool has_hidden_entities() const;
  void get_hidden_surfaces( DLIList<Surface*>& surfaces );
  
  CompositeShell* first_shell() const;
  CompositeShell* next_shell( CompositeShell* after_this ) const;
  
  CubitStatus add( CompositeShell* shell );
  CubitStatus remove( CompositeShell* shell );
  
  CompositeBody* get_body() const;
  
  CubitBox bounding_box() const;
  double measure();
  
  void get_parents_virt( DLIList<TopologyBridge*>& parents );
  void get_children_virt( DLIList<TopologyBridge*>& children );
  int layer() const { return COMPOSITE_LAYER; }
  GeometryQueryEngine* get_geometry_query_engine() const;

  void append_simple_attribute_virt( const CubitSimpleAttrib& simple_attrib_ptr );
  void remove_simple_attribute_virt( const CubitSimpleAttrib& simple_attrib_ptr );
  void remove_all_simple_attribute_virt();
  CubitStatus get_simple_attribute( DLIList<CubitSimpleAttrib>& attrib_list );
  CubitStatus get_simple_attribute( const CubitString& name,
                                    DLIList<CubitSimpleAttrib>& attrib_list );
  
  CubitStatus remove_bridge( TopologyBridge* bridge );
  CubitStatus swap_bridge( TopologyBridge* old_tb, TopologyBridge* new_tb, bool );
  CubitBoolean contains_bridge( TopologyBridge* bridge ) const;
  void notify_reversed( TopologyBridge* bridge );
  
  CompositeLump* split( VGArray<int>& indices_to_move );
  CubitStatus combine( CompositeLump* dead_vol );
  
  void print_debug_info( const char* line_prefix = 0, bool brief = false );

  virtual CubitStatus mass_properties( CubitVector &centroid, double &volume );

private:

  CompositeBody* myBody;
  CompositeLump* nextLump;

  CompositeGeom* compGeom;
  
  CompositeShell* firstShell;
  
  HiddenEntitySet* hiddenSet;
};

inline CompositeShell* CompositeLump::first_shell() const
  { return firstShell; }
  
inline CompositeShell* CompositeLump::next_shell( CompositeShell* after ) const
  { return !after ? firstShell : after->myLump == this ? after->lumpNext : 0; }

inline int CompositeLump::num_lumps() const
  { return compGeom->num_entities(); }

inline Lump* CompositeLump::get_lump( int index ) const
  { return dynamic_cast<Lump*>(compGeom->entity( index )); }

inline int CompositeLump::index_of( Lump* lump ) const
  { return compGeom->index_of( lump ); }

inline void CompositeLump::update()
  { compGeom->update_cached_data(); }

inline CompositeBody* CompositeLump::get_body() const
  { return myBody; }

inline HiddenEntitySet& CompositeLump::hidden_entities()
{ 
  if( !hiddenSet )
    hiddenSet = new HiddenEntitySet(this);
  return *hiddenSet;
}


#endif
