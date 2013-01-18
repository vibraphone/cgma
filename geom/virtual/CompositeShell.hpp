//-------------------------------------------------------------------------
// Filename      : CompositeShell.hpp
//
// Purpose       : ShellSM used in composite TopologyBridge graph
//
// Special Notes : 
//
// Creator       : Jason Kraftcheck
//
// Creation Date : 01/11/02
//-------------------------------------------------------------------------

#ifndef COMPOSITE_SHELL_HPP
#define COMPOSITE_SHELL_HPP

#include "VGDefines.h"
#include "ShellSM.hpp"

class CompositeSurface;
class CompositeCoSurf;
class CompositeLump;

class CompositeShell : public ShellSM
{
friend class CompositeLump;
public:

  CompositeShell();
  virtual ~CompositeShell();
  
  CompositeLump* get_lump() const;
  
  CompositeShell* next_shell() const;
  
  CompositeCoSurf* first_co_surf( ) const;
  CompositeCoSurf* next_co_surf( CompositeCoSurf* prev ) const; 
  
  CubitStatus add( CompositeCoSurf* cosurf );
  CubitStatus remove( CompositeCoSurf* cosurf );
  
  CompositeCoSurf* add( CompositeSurface* surface, CubitSense sense );
    // create a CoSurf
  CompositeCoSurf* find_first( const CompositeSurface* surface ) const;
    // find first CoSurf with the passed surface
  CompositeCoSurf* find_next( const CompositeCoSurf* prev ) const;
    // find next CoSurf with the same surface
  CubitSense find_sense( const CompositeSurface* surface ) const;
    // returns CUBIT_UNKNOWN if multiple CoSurfs

  void append_simple_attribute_virt( const CubitSimpleAttrib& simple_attrib_ptr );
  void remove_simple_attribute_virt( const CubitSimpleAttrib& simple_attrib_ptr );
  void remove_all_simple_attribute_virt();
  CubitStatus get_simple_attribute( DLIList<CubitSimpleAttrib>& attrib_list );
  CubitStatus get_simple_attribute( const CubitString& name,
                                    DLIList<CubitSimpleAttrib>& attrib_list );
  
  void get_parents_virt( DLIList<TopologyBridge*>& parents );
  void get_children_virt( DLIList<TopologyBridge*>& children );
  GeometryQueryEngine* get_geometry_query_engine() const;
  int layer() const { return COMPOSITE_LAYER; }

  void print_debug_info( const char* line_prefix = 0 );

private:

  CompositeLump* myLump;
  CompositeShell* lumpNext;
  
  CompositeCoSurf* firstCoSurf;
};

inline CompositeLump* CompositeShell::get_lump() const
  { return myLump; }

inline CompositeCoSurf* CompositeShell::first_co_surf() const
  { return firstCoSurf; }


#endif
