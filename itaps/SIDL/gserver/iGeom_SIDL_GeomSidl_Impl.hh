// 
// File:          iGeom_SIDL_GeomSidl_Impl.hh
// Symbol:        iGeom_SIDL.GeomSidl-v0.1
// Symbol Type:   class
// Babel Version: 0.10.10
// sidl Created:  20090126 14:50:19 CST
// Generated:     20090126 14:50:21 CST
// Description:   Server-side implementation for iGeom_SIDL.GeomSidl
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.10.10
// source-line   = 5
// source-url    = file:/home/jason/meshkit/cgm/itaps/SIDL/iGeom_SIDL.sidl
// 

#ifndef included_iGeom_SIDL_GeomSidl_Impl_hh
#define included_iGeom_SIDL_GeomSidl_Impl_hh

#ifndef included_sidl_cxx_hh
#include "sidl_cxx.hh"
#endif
#ifndef included_iGeom_SIDL_GeomSidl_IOR_h
#include "iGeom_SIDL_GeomSidl_IOR.h"
#endif
// 
// Includes for all method dependencies.
// 
#ifndef included_iBase_EntityType_hh
#include "iBase_EntityType.hh"
#endif
#ifndef included_iBase_Error_hh
#include "iBase_Error.hh"
#endif
#ifndef included_iBase_StorageOrder_hh
#include "iBase_StorageOrder.hh"
#endif
#ifndef included_iBase_TagValueType_hh
#include "iBase_TagValueType.hh"
#endif
#ifndef included_iGeom_SIDL_GeomSidl_hh
#include "iGeom_SIDL_GeomSidl.hh"
#endif
#ifndef included_sidl_BaseInterface_hh
#include "sidl_BaseInterface.hh"
#endif
#ifndef included_sidl_ClassInfo_hh
#include "sidl_ClassInfo.hh"
#endif


// DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl._includes)
// Insert-Code-Here {iGeom_SIDL.GeomSidl._includes} (includes or arbitrary code)
#include "iGeom.h"
// DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl._includes)

namespace iGeom_SIDL { 

  /**
   * Symbol "iGeom_SIDL.GeomSidl" (version 0.1)
   */
  class GeomSidl_impl
  // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl._inherits)
  // Insert-Code-Here {iGeom_SIDL.GeomSidl._inherits} (optional inheritance here)
    void processError() throw(::iBase::Error);
  // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl._inherits)
  {

  private:
    // Pointer back to IOR.
    // Use this to dispatch back through IOR vtable.
    GeomSidl self;

    // DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl._implementation)
    // Insert-Code-Here {iGeom_SIDL.GeomSidl._implementation} (additional details)
    iGeom_Instance igeomInstance;
    // DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl._implementation)

  private:
    // private default constructor (required)
    GeomSidl_impl() 
    {} 

  public:
    // sidl constructor (required)
    // Note: alternate Skel constructor doesn't call addref()
    // (fixes bug #275)
    GeomSidl_impl( struct iGeom_SIDL_GeomSidl__object * s ) : self(s,
      true) { _ctor(); }

    // user defined construction
    void _ctor();

    // virtual destructor (required)
    virtual ~GeomSidl_impl() { _dtor(); }

    // user defined destruction
    void _dtor();

    // static class initializer
    static void _load();

  public:

    /**
     * user defined non-static method.
     */
    void
    createTag (
      /* in */ const ::std::string& tag_name,
      /* in */ int32_t number_of_values,
      /* in */ ::iBase::TagValueType tag_type,
      /* out */ void*& tag_handle
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    destroyTag (
      /* in */ void* tag_handle,
      /* in */ bool forced
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getTagName (
      /* in */ void* tag_handle,
      /* out */ ::std::string& tag_name
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getTagSizeValues (
      /* in */ void* tag_handle,
      /* out */ int32_t& size_values
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getTagSizeBytes (
      /* in */ void* tag_handle,
      /* out */ int32_t& size_bytes
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getTagHandle (
      /* in */ const ::std::string& tag_name,
      /* out */ void*& tag_handle
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getTagType (
      /* in */ void* tag_handle,
      /* out */ ::iBase::TagValueType& tag_data_type
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getData (
      /* in */ void* entity_handle,
      /* in */ void* tag_handle,
      /* inout */ ::sidl::array<char>& tag_value,
      /* out */ int32_t& tag_value_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getIntData (
      /* in */ void* entity_handle,
      /* in */ void* tag_handle,
      /* out */ int32_t& int_data
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getDblData (
      /* in */ void* entity_handle,
      /* in */ void* tag_handle,
      /* out */ double& dbl_data
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getEHData (
      /* in */ void* entity_handle,
      /* in */ void* tag_handle,
      /* out */ void*& eh_data
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    setData (
      /* in */ void* entity_handle,
      /* in */ void* tag_handle,
      /* in */ ::sidl::array<char> tag_value,
      /* in */ int32_t tag_value_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    setIntData (
      /* in */ void* entity_handle,
      /* in */ void* tag_handle,
      /* in */ int32_t tag_value
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    setDblData (
      /* in */ void* entity_handle,
      /* in */ void* tag_handle,
      /* in */ double tag_value
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    setEHData (
      /* in */ void* entity_handle,
      /* in */ void* tag_handle,
      /* in */ void* tag_value
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getAllTags (
      /* in */ void* entity_handle,
      /* inout */ ::sidl::array<void*>& tag_handles,
      /* out */ int32_t& tag_handles_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    rmvTag (
      /* in */ void* entity_handle,
      /* in */ void* tag_handle
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getArrData (
      /* in */ ::sidl::array<void*> entity_handles,
      /* in */ int32_t entity_handles_size,
      /* in */ void* tag_handle,
      /* inout */ ::sidl::array<char>& value_array,
      /* out */ int32_t& value_array_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getIntArrData (
      /* in */ ::sidl::array<void*> entity_handles,
      /* in */ int32_t entity_handles_size,
      /* in */ void* tag_handle,
      /* inout */ ::sidl::array<int32_t>& value_array,
      /* out */ int32_t& value_array_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getDblArrData (
      /* in */ ::sidl::array<void*> entity_handles,
      /* in */ int32_t entity_handles_size,
      /* in */ void* tag_handle,
      /* inout */ ::sidl::array<double>& value_array,
      /* out */ int32_t& value_array_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getEHArrData (
      /* in */ ::sidl::array<void*> entity_handles,
      /* in */ int32_t entity_handles_size,
      /* in */ void* tag_handle,
      /* inout */ ::sidl::array<void*>& value_array,
      /* out */ int32_t& value_array_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    setArrData (
      /* in */ ::sidl::array<void*> entity_handles,
      /* in */ int32_t entity_handles_size,
      /* in */ void* tag_handle,
      /* in */ ::sidl::array<char> value_array,
      /* in */ int32_t value_array_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    setIntArrData (
      /* in */ ::sidl::array<void*> entity_handles,
      /* in */ int32_t entity_handles_size,
      /* in */ void* tag_handle,
      /* in */ ::sidl::array<int32_t> value_array,
      /* in */ int32_t value_array_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    setDblArrData (
      /* in */ ::sidl::array<void*> entity_handles,
      /* in */ int32_t entity_handles_size,
      /* in */ void* tag_handle,
      /* in */ ::sidl::array<double> value_array,
      /* in */ int32_t value_array_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    setEHArrData (
      /* in */ ::sidl::array<void*> entity_handles,
      /* in */ int32_t entity_handles_size,
      /* in */ void* tag_handle,
      /* in */ ::sidl::array<void*> value_array,
      /* in */ int32_t value_array_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    rmvArrTag (
      /* in */ ::sidl::array<void*> entity_handles,
      /* in */ int32_t entity_handles_size,
      /* in */ void* tag_handle
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    setEntSetData (
      /* in */ void* entity_set,
      /* in */ void* tag_handle,
      /* inout */ ::sidl::array<char>& tag_value,
      /* in */ int32_t tag_value_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    setEntSetIntData (
      /* in */ void* entity_set,
      /* in */ void* tag_handle,
      /* in */ int32_t tag_value
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    setEntSetDblData (
      /* in */ void* entity_set,
      /* in */ void* tag_handle,
      /* in */ double tag_value
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    setEntSetEHData (
      /* in */ void* entity_set,
      /* in */ void* tag_handle,
      /* in */ void* tag_value
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getEntSetData (
      /* in */ void* entity_set,
      /* in */ void* tag_handle,
      /* inout */ ::sidl::array<char>& tag_value,
      /* out */ int32_t& tag_value_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getEntSetIntData (
      /* in */ void* entity_set,
      /* in */ void* tag_handle,
      /* out */ int32_t& int_data
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getEntSetDblData (
      /* in */ void* entity_set,
      /* in */ void* tag_handle,
      /* out */ double& dbl_data
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getEntSetEHData (
      /* in */ void* entity_set,
      /* in */ void* tag_handle,
      /* out */ void*& eh_data
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getAllEntSetTags (
      /* in */ void* entity_set,
      /* inout */ ::sidl::array<void*>& tag_handles,
      /* out */ int32_t& tag_handles_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    rmvEntSetTag (
      /* in */ void* entity_set,
      /* in */ void* tag_handle
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    createEntSet (
      /* in */ bool isList,
      /* out */ void*& entity_set
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    destroyEntSet (
      /* in */ void* entity_set
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    isList (
      /* in */ void* entity_set,
      /* out */ int32_t& is_list
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getNumEntSets (
      /* in */ void* entity_set,
      /* in */ int32_t num_hops,
      /* out */ int32_t& num_sets
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getEntSets (
      /* in */ void* entity_set,
      /* in */ int32_t num_hops,
      /* inout */ ::sidl::array<void*>& contained_entset_handles,
      /* out */ int32_t& contained_entset_handles_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    addEntToSet (
      /* in */ void* entity_handle,
      /* in */ void* entity_set
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    rmvEntFromSet (
      /* in */ void* entity_handle,
      /* in */ void* entity_set
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    addEntArrToSet (
      /* in */ ::sidl::array<void*> entity_handles,
      /* in */ int32_t entity_handles_size,
      /* in */ void* entity_set
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    rmvEntArrFromSet (
      /* in */ ::sidl::array<void*> entity_handles,
      /* in */ int32_t entity_handles_size,
      /* in */ void* entity_set
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    addEntSet (
      /* in */ void* entity_set_to_add,
      /* in */ void* entity_set_handle
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    rmvEntSet (
      /* in */ void* entity_set_to_remove,
      /* in */ void* entity_set_handle
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    isEntContained (
      /* in */ void* containing_entity_set,
      /* in */ void* entity_handle,
      /* out */ int32_t& is_contained
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    isEntArrContained (
      /* in */ void* containing_set,
      /* in */ ::sidl::array<void*> entity_handles,
      /* in */ int32_t entity_handles_size,
      /* inout */ ::sidl::array<int32_t>& is_contained,
      /* out */ int32_t& is_contained_size
    )
    throw () 
    ;

    /**
     * user defined non-static method.
     */
    void
    isEntSetContained (
      /* in */ void* containing_entity_set,
      /* in */ void* contained_entity_set,
      /* out */ int32_t& is_contained
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    addPrntChld (
      /* in */ void* parent_entity_set,
      /* in */ void* child_entity_set
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    rmvPrntChld (
      /* in */ void* parent_entity_set,
      /* in */ void* child_entity_set
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    isChildOf (
      /* in */ void* parent_entity_set,
      /* in */ void* child_entity_set,
      /* out */ int32_t& is_child
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getNumChld (
      /* in */ void* entity_set,
      /* in */ int32_t num_hops,
      /* out */ int32_t& num_child
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getNumPrnt (
      /* in */ void* entity_set,
      /* in */ int32_t num_hops,
      /* out */ int32_t& num_parent
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getChldn (
      /* in */ void* from_entity_set,
      /* in */ int32_t num_hops,
      /* inout */ ::sidl::array<void*>& child_handles,
      /* out */ int32_t& child_handles_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getPrnts (
      /* in */ void* from_entity_set,
      /* in */ int32_t num_hops,
      /* inout */ ::sidl::array<void*>& parent_handles,
      /* out */ int32_t& parent_handles_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    subtract (
      /* in */ void* entity_set_1,
      /* in */ void* entity_set_2,
      /* out */ void*& result_entity_set
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    intersect (
      /* in */ void* entity_set_1,
      /* in */ void* entity_set_2,
      /* out */ void*& result_entity_set
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    unite (
      /* in */ void* entity_set_1,
      /* in */ void* entity_set_2,
      /* out */ void*& result_entity_set
    )
    throw ( 
      ::iBase::Error
    );


    /**
     * Load a model specified by name. Which formats are supported and the
     * specific meaning of this name string (e.g. file name, model name,
     * etc.) are implementation-dependent.  Options are also implementation-
     * dependent.
     * @param name Name of the model
     * @param options String options 
     */
    void
    load (
      /* in */ const ::std::string& name,
      /* in */ const ::std::string& options
    )
    throw ( 
      ::iBase::Error
    );


    /**
     * Save a model to a file specified by name. Which formats are supported and the
     * specific meaning of this name string (e.g. file name, model name,
     * etc.) are implementation-dependent.  Options are also implementation-
     * dependent.
     * @param name Name of the file to save to
     * @param options String options 
     */
    void
    save (
      /* in */ const ::std::string& name,
      /* in */ const ::std::string& options
    )
    throw ( 
      ::iBase::Error
    );


    /**
     * Return gentities of specified dimension in this set, or in whole model.
     * @param set_handle Entity set being queried (if 0, whole model)
     * @param gentity_dimension Dimension of entities being queried
     * @param gentity_handles Gentity handles
     */
    void
    getEntities (
      /* in */ void* set_handle,
      /* in */ ::iBase::EntityType gentity_type,
      /* inout */ ::sidl::array<void*>& gentity_handles,
      /* inout */ int32_t& gentity_handles_size
    )
    throw ( 
      ::iBase::Error
    );


    /**
     * Return number of gentities of specified dimension in this set, or in
     * whole model.
     * @param set_handle Entity set being queried (if 0, whole model)
     * @param gentity_dimension Dimension of entities being queried
     * @return Number of entities
     */
    void
    getNumOfType (
      /* in */ void* set_handle,
      /* in */ ::iBase::EntityType gentity_type,
      /* out */ int32_t& num_type
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    ::iBase::EntityType
    getEntType (
      /* in */ void* handle
    )
    throw ( 
      ::iBase::Error
    );


    /**
     *    Returns an integer array of topological dimensions for an input
     *    array of entity handles.
     */
    void
    getArrType (
      /* in */ ::sidl::array<void*> gentity_handles,
      /* in */ int32_t gentity_handles_size,
      /* inout */ ::sidl::array< ::iBase::EntityType>& gtype,
      /* inout */ int32_t& gtype_size
    )
    throw ( 
      ::iBase::Error
    );


    /**
     * Get the adjacent entities of a given dimension.
     * @param gentity_handle Entity for which adjacencies are requested
     * @param to_dimension Target dimension of adjacent entities
     * @param adj_gentities List returned with adjacent entities
     */
    void
    getEntAdj (
      /* in */ void* entity_handle,
      /* in */ ::iBase::EntityType to_dimension,
      /* inout */ ::sidl::array<void*>& adj_gentities,
      /* inout */ int32_t& adj_gentities_size
    )
    throw ( 
      ::iBase::Error
    );


    /**
     * Get the adjacent entities of a given dimension.
     * @param gentity_handle Entity for which adjacencies are requested
     * @param to_dimension Target dimension of adjacent entities
     * @param adj_gentities List returned with adjacent entities
     */
    void
    getArrAdj (
      /* in */ ::sidl::array<void*> entity_handles,
      /* in */ int32_t entity_handles_size,
      /* in */ ::iBase::EntityType requested_entity_type,
      /* inout */ ::sidl::array<void*>& adj_entity_handles,
      /* inout */ int32_t& adj_entity_handles_size,
      /* inout */ ::sidl::array<int32_t>& offset,
      /* inout */ int32_t& offset_size
    )
    throw ( 
      ::iBase::Error
    );


    /**
     * Get the "2nd order" adjacent entities, through a specified "bridge"
     * dimension, of a target dimension.  For example, given a region, return
     * the regions (to_dimension=3) sharing an edge (bridge_dimension=1)
     * with that region.  bridge_dimension must be less than dimension of 
     * gentity_handle, and to_dimension must be greater than bridge dimension.
     * 
     * @param gentity_handle Entity for which 2nd order adjacencies are requested
     * @param to_dimension Target dimension of 2nd order adjacent entities
     * @param bridge_dimension Dimension of "bridge" entities
     * @param adj_gentities List returned with 2nd order adjacent entities
     */
    void
    getEnt2ndAdj (
      /* in */ void* gentity_handle,
      /* in */ ::iBase::EntityType bridge_dimension,
      /* in */ ::iBase::EntityType to_dimension,
      /* inout */ ::sidl::array<void*>& adjacent_gentities,
      /* inout */ int32_t& adjacent_gentities_size
    )
    throw ( 
      ::iBase::Error
    );


    /**
     * Get the "2nd order" adjacent entities, through a specified "bridge"
     * dimension, of a target dimension.  For example, given a region, return
     * the regions (to_dimension=3) sharing an edge (bridge_dimension=1)
     * with that region.  bridge_dimension must be less than dimension of 
     * gentity_handle, and to_dimension must be greater than bridge dimension.
     * 
     * @param gentity_handle Entity for which 2nd order adjacencies are requested
     * @param order_adjacent_key Dimension of "bridge" entities
     * @param requested_entity_type Target dimension of 2nd order adjacent entities
     * @param adj_entity_handles List returned with 2nd order adjacent entities
     */
    void
    getArr2ndAdj (
      /* in */ ::sidl::array<void*> entity_handles,
      /* in */ int32_t entity_handles_size,
      /* in */ ::iBase::EntityType order_adjacent_key,
      /* in */ ::iBase::EntityType requested_entity_type,
      /* inout */ ::sidl::array<void*>& adj_entity_handles,
      /* inout */ int32_t& adj_entity_handles_size,
      /* inout */ ::sidl::array<int32_t>& offset,
      /* inout */ int32_t& offset_size
    )
    throw ( 
      ::iBase::Error
    );


    /**
     * Return whether or not entities are adjacent.
     * @param gentity_handle1 1st entity
     * @param gentity_handle2 2nd entity
     * @param are_adjacent If true, entities are adjacent
     */
    int32_t
    isEntAdj (
      /* in */ void* gentity_handle1,
      /* in */ void* gentity_handle2,
      /* out */ bool& are_adjacent
    )
    throw ( 
      ::iBase::Error
    );


    /**
     * Return whether or not arrays of entities are adjacent.
     * @param gentity_handle1 1st entity
     * @param gentity_handle2 2nd entity
     * @param are_adjacent If true, entities are adjacent
     */
    void
    isArrAdj (
      /* in */ ::sidl::array<void*> entity_handles_1,
      /* in */ int32_t entity_handles_1_size,
      /* in */ ::sidl::array<void*> entity_handles_2,
      /* in */ int32_t entity_handles_2_size,
      /* inout */ ::sidl::array<int32_t>& is_adjacent_info,
      /* inout */ int32_t& is_adjacent_info_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    int32_t
    getTopoLevel (
      /* in */ const ::std::string& model_name
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getEntClosestPt (
      /* in */ void* entity_handle,
      /* in */ double near_x,
      /* in */ double near_y,
      /* in */ double near_z,
      /* out */ double& on_x,
      /* out */ double& on_y,
      /* out */ double& on_z
    )
    throw ( 
      ::iBase::Error
    );


    /**
     * Return a points on specified entities closest to specified points
     * in space.  Input coordinates and output points are interleaved in 
     * the arrays.
     * @param gentity_handles The gentities being queried
     * @param near_coordinates Input coordinates
     * @param on_coordinates Closest point on gentity
     */
    void
    getArrClosestPt (
      /* in */ ::sidl::array<void*> gentity_handles,
      /* in */ int32_t gentity_handles_size,
      /* in */ ::iBase::StorageOrder storage_order,
      /* in */ ::sidl::array<double> near_coordinates,
      /* in */ int32_t near_coordinates_size,
      /* inout */ ::sidl::array<double>& on_coordinates,
      /* inout */ int32_t& on_coordinates_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getEntNrmlXYZ (
      /* in */ void* entity_handle,
      /* in */ double x,
      /* in */ double y,
      /* in */ double z,
      /* out */ double& nrml_i,
      /* out */ double& nrml_j,
      /* out */ double& nrml_k
    )
    throw ( 
      ::iBase::Error
    );


    /**
     * Return the normals at point on specified entities.  Returns error
     * if any input entity is not a gface.  Input coordinates and normals
     * are interleaved in the arrays.
     * @param gentity_handles The gentities being queried
     * @param coordinates Input coordinates, interleaved
     * @param normals The normals at the specified points, interleaved
     */
    void
    getArrNrmlXYZ (
      /* in */ ::sidl::array<void*> gentity_handles,
      /* in */ int32_t gentity_handles_size,
      /* in */ ::iBase::StorageOrder storage_order,
      /* in */ ::sidl::array<double> coordinates,
      /* in */ int32_t coordinates_size,
      /* inout */ ::sidl::array<double>& normals,
      /* inout */ int32_t& normals_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getEntNrmlPlXYZ (
      /* in */ void* entity_handle,
      /* in */ double x,
      /* in */ double y,
      /* in */ double z,
      /* out */ double& pt_x,
      /* out */ double& pt_y,
      /* out */ double& pt_z,
      /* out */ double& nrml_i,
      /* out */ double& nrml_j,
      /* out */ double& nrml_k
    )
    throw ( 
      ::iBase::Error
    );


    /**
     * Return points and normals on specified entities closest to specified points
     * in space.  Input coordinates and output points are interleaved in 
     * the arrays.
     * @param gentity_handles The gentities being queried
     * @param near_coordinates Input coordinates
     * @param on_coordinates Closest point on gentity
     * @param normals Normals on gentity
     */
    void
    getArrNrmlPlXYZ (
      /* in */ ::sidl::array<void*> gentity_handles,
      /* in */ int32_t gentity_handles_size,
      /* in */ ::iBase::StorageOrder storage_order,
      /* in */ ::sidl::array<double> near_coordinates,
      /* in */ int32_t near_coordinates_size,
      /* inout */ ::sidl::array<double>& on_coordinates,
      /* inout */ int32_t& on_coordinates_size,
      /* inout */ ::sidl::array<double>& normals,
      /* inout */ int32_t& normals_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getEntTgntXYZ (
      /* in */ void* entity_handle,
      /* in */ double x,
      /* in */ double y,
      /* in */ double z,
      /* out */ double& tgnt_i,
      /* out */ double& tgnt_j,
      /* out */ double& tgnt_k
    )
    throw ( 
      ::iBase::Error
    );


    /**
     * Return the tangent at point on specified entities.  Returns error
     * if any input entity is not a gedge.  Input coordinates and tangents
     * are interleaved in the arrays.
     * @param gentity_handles The gentities being queried
     * @param coordinates Input coordinates, interleaved
     * @param tangents The tangents at the specified points, interleaved
     */
    void
    getArrTgntXYZ (
      /* in */ ::sidl::array<void*> gentity_handles,
      /* in */ int32_t gentity_handles_size,
      /* in */ ::iBase::StorageOrder storage_order,
      /* in */ ::sidl::array<double> coordinates,
      /* in */ int32_t coordinates_size,
      /* inout */ ::sidl::array<double>& tangents,
      /* inout */ int32_t& tangents_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getFcCvtrXYZ (
      /* in */ void* face_handle,
      /* in */ double x,
      /* in */ double y,
      /* in */ double z,
      /* out */ double& cvtr1_i,
      /* out */ double& cvtr1_j,
      /* out */ double& cvtr1_k,
      /* out */ double& cvtr2_i,
      /* out */ double& cvtr2_j,
      /* out */ double& cvtr2_k
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getEgCvtrXYZ (
      /* in */ void* edge_handle,
      /* in */ double x,
      /* in */ double y,
      /* in */ double z,
      /* out */ double& cvtr_i,
      /* out */ double& cvtr_j,
      /* out */ double& cvtr_k
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getEntArrCvtrXYZ (
      /* in */ ::sidl::array<void*> entity_handles,
      /* in */ int32_t entity_handles_size,
      /* in */ ::iBase::StorageOrder storage_order,
      /* in */ ::sidl::array<double> coords,
      /* in */ int32_t coords_size,
      /* inout */ ::sidl::array<double>& cvtr_1,
      /* inout */ int32_t& cvtr_1_size,
      /* inout */ ::sidl::array<double>& cvtr_2,
      /* inout */ int32_t& cvtr_2_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getEgEvalXYZ (
      /* in */ void* edge_handle,
      /* in */ double x,
      /* in */ double y,
      /* in */ double z,
      /* out */ double& on_x,
      /* out */ double& on_y,
      /* out */ double& on_z,
      /* out */ double& tgnt_i,
      /* out */ double& tgnt_j,
      /* out */ double& tgnt_k,
      /* out */ double& cvtr_i,
      /* out */ double& cvtr_j,
      /* out */ double& cvtr_k
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getFcEvalXYZ (
      /* in */ void* face_handle,
      /* in */ double x,
      /* in */ double y,
      /* in */ double z,
      /* out */ double& on_x,
      /* out */ double& on_y,
      /* out */ double& on_z,
      /* out */ double& nrml_i,
      /* out */ double& nrml_j,
      /* out */ double& nrml_k,
      /* out */ double& cvtr1_i,
      /* out */ double& cvtr1_j,
      /* out */ double& cvtr1_k,
      /* out */ double& cvtr2_i,
      /* out */ double& cvtr2_j,
      /* out */ double& cvtr2_k
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getArrEgEvalXYZ (
      /* in */ ::sidl::array<void*> edge_handles,
      /* in */ int32_t edge_handles_size,
      /* in */ ::iBase::StorageOrder storage_order,
      /* in */ ::sidl::array<double> coords,
      /* in */ int32_t coords_size,
      /* inout */ ::sidl::array<double>& on_coords,
      /* inout */ int32_t& on_coords_size,
      /* inout */ ::sidl::array<double>& tangent,
      /* inout */ int32_t& tangent_size,
      /* inout */ ::sidl::array<double>& cvtr,
      /* inout */ int32_t& cvtr_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getArrFcEvalXYZ (
      /* in */ ::sidl::array<void*> face_handles,
      /* in */ int32_t face_handles_size,
      /* in */ ::iBase::StorageOrder storage_order,
      /* in */ ::sidl::array<double> coords,
      /* in */ int32_t coords_size,
      /* inout */ ::sidl::array<double>& on_coords,
      /* inout */ int32_t& on_coords_size,
      /* inout */ ::sidl::array<double>& normal,
      /* inout */ int32_t& normal_size,
      /* inout */ ::sidl::array<double>& cvtr_1,
      /* inout */ int32_t& cvtr_1_size,
      /* inout */ ::sidl::array<double>& cvtr_2,
      /* inout */ int32_t& cvtr_2_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getEntBoundBox (
      /* in */ void* entity_handle,
      /* out */ double& min_x,
      /* out */ double& min_y,
      /* out */ double& min_z,
      /* out */ double& max_x,
      /* out */ double& max_y,
      /* out */ double& max_z
    )
    throw ( 
      ::iBase::Error
    );


    /**
     * Return the bounding boxex of given entities; coordinates returned
     * interleaved.
     * @param gentity_handles The gentities being queried
     * @param min_corners Minimum corner coordinates of the boxes, interleaved
     * @param max_corners Maximum corner coordinates of the boxes, interleaved
     */
    void
    getArrBoundBox (
      /* in */ ::sidl::array<void*> gentity_handles,
      /* in */ int32_t gentity_handles_size,
      /* inout */ ::iBase::StorageOrder& storage_order,
      /* inout */ ::sidl::array<double>& min_corner,
      /* inout */ int32_t& min_corner_size,
      /* inout */ ::sidl::array<double>& max_corner,
      /* inout */ int32_t& max_corner_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getVtxCoord (
      /* in */ void* vertex_handle,
      /* out */ double& x,
      /* out */ double& y,
      /* out */ double& z
    )
    throw ( 
      ::iBase::Error
    );


    /**
     * Return the coordinates of the specified vertices; returns error if any
     * of the entities are not gvertices.  Coordinates returned interleaved.
     * @param gentity_handles The gentities being queried
     * @param coordinates The coordinates of the gvertices, interleaved.
     */
    void
    getVtxArrCoords (
      /* in */ ::sidl::array<void*> gentity_handles,
      /* in */ int32_t gentity_handles_size,
      /* inout */ ::iBase::StorageOrder& storage_order,
      /* inout */ ::sidl::array<double>& coordinates,
      /* inout */ int32_t& coordinates_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getPntIntsct (
      /* in */ double x,
      /* in */ double y,
      /* in */ double z,
      /* in */ double dir_x,
      /* in */ double dir_y,
      /* in */ double dir_z,
      /* inout */ ::sidl::array<void*>& intersect_entity_handles,
      /* inout */ int32_t& intersect_entity_handles_size,
      /* inout */ ::iBase::StorageOrder& storage_order,
      /* inout */ ::sidl::array<double>& intersect_coords,
      /* inout */ int32_t& intersect_coords_size,
      /* inout */ ::sidl::array<double>& param_coords,
      /* inout */ int32_t& param_coords_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getPntArrRayIntsct (
      /* in */ ::iBase::StorageOrder storage_order,
      /* in */ ::sidl::array<double> coords,
      /* in */ int32_t coords_size,
      /* in */ ::sidl::array<double> directions,
      /* in */ int32_t directions_size,
      /* inout */ ::sidl::array<void*>& intersect_entity_handles,
      /* inout */ int32_t& intersect_entity_handles_size,
      /* inout */ ::sidl::array<int32_t>& offset,
      /* inout */ int32_t& offset_size,
      /* inout */ ::sidl::array<double>& intersect_coords,
      /* inout */ int32_t& intersect_coords_size,
      /* inout */ ::sidl::array<double>& param_coords,
      /* inout */ int32_t& param_coords_size
    )
    throw ( 
      ::iBase::Error
    );


    /**
     * Return the sense of a gface with respect to a gregion.  Sense is either
     * forward (=1), reverse (=-1), both (=2), or unknown (=0).  Error is returned
     * if first entity is not a gface or second entity is not a gregion.
     * @param gface Gface whose sense is being queried.
     * @param gregion Gregion gface is being queried with respect to
     */
    int32_t
    getEntNrmlSense (
      /* in */ void* gface,
      /* in */ void* gregion
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getArrNrmlSense (
      /* in */ ::sidl::array<void*> face_handles,
      /* in */ int32_t face_handles_size,
      /* in */ ::sidl::array<void*> region_handles,
      /* in */ int32_t region_handles_size,
      /* inout */ ::sidl::array<int32_t>& sense,
      /* inout */ int32_t& sense_size
    )
    throw ( 
      ::iBase::Error
    );


    /**
     * Return the sense of a gedge with respect to a gface.  Sense is either
     * forward (=1), reverse (=-1), both (=2), or unknown (=0).  Error is returned
     * if first entity is not a gedge or second entity is not a gface.
     * @param gedge Gedge whose sense is being queried.
     * @param gface Gface gedge is being queried with respect to
     */
    int32_t
    getEgFcSense (
      /* in */ void* gedge,
      /* in */ void* gface
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getEgFcArrSense (
      /* in */ ::sidl::array<void*> edge_handles,
      /* in */ int32_t edge_handles_size,
      /* in */ ::sidl::array<void*> face_handles,
      /* in */ int32_t face_handles_size,
      /* inout */ ::sidl::array<int32_t>& sense,
      /* inout */ int32_t& sense_size
    )
    throw ( 
      ::iBase::Error
    );


    /**
     * Return the sense of a gedge with respect to a specified order of
     * vertices bounding the gedge.  Sense is either forward (=1), reverse (=-1), 
     * or unknown (=0).  Error is returned if any gentities are not the expected
     * type or if the gedge is bounded by only one gvertex (in this case, use
     * getGtangentSense).
     * @param gedge Gedge whose sense is being queried.
     * @param gvertex1 First gvertex
     * @param gvertex2 Second gvertex
     */
    int32_t
    getEgVtxSense (
      /* in */ void* gedge,
      /* in */ void* gvertex1,
      /* in */ void* gvertex2
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getEgVtxArrSense (
      /* in */ ::sidl::array<void*> edge_handles,
      /* in */ int32_t edge_handles_size,
      /* in */ ::sidl::array<void*> vertex_handles_1,
      /* in */ int32_t vertex_handles_1_size,
      /* in */ ::sidl::array<void*> vertex_handles_2,
      /* in */ int32_t vertex_handles_2_size,
      /* inout */ ::sidl::array<int32_t>& sense,
      /* inout */ int32_t& sense_size
    )
    throw ( 
      ::iBase::Error
    );


    /**
     * Return the arc length / area / volume of the entities
     * @param gentity_handles Entities for which measure is requested
     * @param gentity_handles_size Number of gentities
     * @param measures Arc length / area / volume of the entities
     * @param measures_length Number of entries in measures
     */
    void
    measure (
      /* in */ ::sidl::array<void*> gentity_handles,
      /* in */ int32_t gentity_handles_size,
      /* inout */ ::sidl::array<double>& measures,
      /* inout */ int32_t& measures_size
    )
    throw ( 
      ::iBase::Error
    );


    /**
     * Return the type of surface as a string; if not a surface, an error is returned
     * @param face_handle Face for which the type is requested
     * @param face_type Type of face, returned as a string
     */
    void
    getFaceType (
      /* in */ void* gface_handle,
      /* inout */ ::std::string& face_type
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    int32_t
    getParametric() throw ( 
      ::iBase::Error
    );

    /**
     * Return whether a given gentity is parametric or not.  If a gentity
     * is not parametric, all of the following functions will return an error
     * when called on that entity.
     * @param gentity_handle Gentity being queried.
     */
    int32_t
    isEntParametric (
      /* in */ void* gentity_handle
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    isArrParametric (
      /* in */ ::sidl::array<void*> entity_handles,
      /* in */ int32_t entity_handles_size,
      /* inout */ ::sidl::array<int32_t>& is_parametric,
      /* inout */ int32_t& is_parametric_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getEntUVtoXYZ (
      /* in */ void* entity_handle,
      /* in */ double u,
      /* in */ double v,
      /* out */ double& x,
      /* out */ double& y,
      /* out */ double& z
    )
    throw ( 
      ::iBase::Error
    );


    /**
     * Given sets of parametric coordinates, return the corresponding real
     * space coordinates on the gentities.  Input and output coordinates are
     * interleaved.
     * @param gentity_handles Gentities being queried.
     * @param uv Input parametric coordinates
     * @param xyz Output real space coordinates
     */
    void
    getArrUVtoXYZ (
      /* in */ ::sidl::array<void*> gentity_handles,
      /* in */ int32_t gentity_handles_size,
      /* in */ ::iBase::StorageOrder storage_order,
      /* in */ ::sidl::array<double> uv,
      /* in */ int32_t uv_size,
      /* inout */ ::sidl::array<double>& coordinates,
      /* inout */ int32_t& coordinates_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getEntUtoXYZ (
      /* in */ void* entity_handle,
      /* in */ double u,
      /* out */ double& x,
      /* out */ double& y,
      /* out */ double& z
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getArrUtoXYZ (
      /* in */ ::sidl::array<void*> entity_handles,
      /* in */ int32_t entity_handles_size,
      /* in */ ::sidl::array<double> u,
      /* in */ int32_t u_size,
      /* inout */ ::iBase::StorageOrder& storage_order,
      /* inout */ ::sidl::array<double>& on_coords,
      /* inout */ int32_t& on_coords_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getEntXYZtoUV (
      /* in */ void* entity_handle,
      /* in */ double x,
      /* in */ double y,
      /* in */ double z,
      /* out */ double& u,
      /* out */ double& v
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getEntXYZtoU (
      /* in */ void* entity_handle,
      /* in */ double x,
      /* in */ double y,
      /* in */ double z,
      /* out */ double& u
    )
    throw ( 
      ::iBase::Error
    );


    /**
     * Given sets of real space coordinates, return the corresponding 
     * parametric coordinates on the gentities.  Input and output coordinates 
     * are interleaved.
     * @param gentity_handles Gentities being queried.
     * @param xyz Input real space coordinates
     * @param uv Output parametric coordinates
     */
    void
    getArrXYZtoUV (
      /* in */ ::sidl::array<void*> gentity_handles,
      /* in */ int32_t gentity_handles_size,
      /* in */ ::iBase::StorageOrder storage_order,
      /* in */ ::sidl::array<double> coordinates,
      /* in */ int32_t coordinates_size,
      /* inout */ ::sidl::array<double>& uv,
      /* inout */ int32_t& uv_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getArrXYZtoU (
      /* in */ ::sidl::array<void*> gentity_handles,
      /* in */ int32_t gentity_handles_size,
      /* in */ ::iBase::StorageOrder storage_order,
      /* in */ ::sidl::array<double> coordinates,
      /* in */ int32_t coordinates_size,
      /* inout */ ::sidl::array<double>& u,
      /* inout */ int32_t& u_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getEntXYZtoUVHint (
      /* in */ void* entity_handle,
      /* in */ double x,
      /* in */ double y,
      /* in */ double z,
      /* inout */ double& u,
      /* inout */ double& v
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getArrXYZtoUVHint (
      /* in */ ::sidl::array<void*> entity_handles,
      /* in */ int32_t entity_handles_size,
      /* in */ ::iBase::StorageOrder storage_order,
      /* in */ ::sidl::array<double> coords,
      /* in */ int32_t coords_size,
      /* inout */ ::sidl::array<double>& uv,
      /* inout */ int32_t& uv_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getEntUVRange (
      /* in */ void* entity_handle,
      /* out */ double& u_min,
      /* out */ double& v_min,
      /* out */ double& u_max,
      /* out */ double& v_max
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getEntURange (
      /* in */ void* entity_handle,
      /* out */ double& u_min,
      /* out */ double& u_max
    )
    throw ( 
      ::iBase::Error
    );


    /**
     * Return the uv range of the specified gentities.  Parameters are interleaved.
     * @param gentity_handles Gentities being queried.
     * @param uv_min Minimum parameters of gentities, interleaved
     * @param uv_max Maximum parameters of gentities, interleaved
     */
    void
    getArrUVRange (
      /* in */ ::sidl::array<void*> gentity_handles,
      /* in */ int32_t gentity_handles_size,
      /* inout */ ::iBase::StorageOrder& storage_order,
      /* inout */ ::sidl::array<double>& uv_min,
      /* inout */ int32_t& uv_min_size,
      /* inout */ ::sidl::array<double>& uv_max,
      /* inout */ int32_t& uv_max_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getArrURange (
      /* in */ ::sidl::array<void*> entity_handles,
      /* in */ int32_t entity_handles_size,
      /* inout */ ::sidl::array<double>& u_min,
      /* inout */ int32_t& u_min_size,
      /* inout */ ::sidl::array<double>& u_max,
      /* inout */ int32_t& u_max_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getEntUtoUV (
      /* in */ void* edge_handle,
      /* in */ void* face_handle,
      /* in */ double in_u,
      /* out */ double& u,
      /* out */ double& v
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getVtxToUV (
      /* in */ void* vertex_handle,
      /* in */ void* face_handle,
      /* out */ double& u,
      /* out */ double& v
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getVtxToU (
      /* in */ void* vertex_handle,
      /* in */ void* edge_handle,
      /* out */ double& u
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getArrUtoUV (
      /* in */ ::sidl::array<void*> edge_handles,
      /* in */ int32_t edge_handles_size,
      /* in */ ::sidl::array<void*> face_handles,
      /* in */ int32_t face_handles_size,
      /* in */ ::sidl::array<double> u_in,
      /* in */ int32_t u_in_size,
      /* inout */ ::iBase::StorageOrder& storage_order,
      /* inout */ ::sidl::array<double>& uv,
      /* inout */ int32_t& uv_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getVtxArrToUV (
      /* in */ ::sidl::array<void*> vertex_handles,
      /* in */ int32_t vertex_handles_size,
      /* in */ ::sidl::array<void*> face_handles,
      /* in */ int32_t face_handles_size,
      /* inout */ ::iBase::StorageOrder& storage_order,
      /* inout */ ::sidl::array<double>& uv,
      /* inout */ int32_t& uv_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getVtxArrToU (
      /* in */ ::sidl::array<void*> vertex_handles,
      /* in */ int32_t vertex_handles_size,
      /* in */ ::sidl::array<void*> edge_handles,
      /* in */ int32_t edge_handles_size,
      /* inout */ ::sidl::array<double>& u,
      /* inout */ int32_t& u_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getEntNrmlUV (
      /* in */ void* entity_handle,
      /* in */ double u,
      /* in */ double v,
      /* out */ double& nrml_i,
      /* out */ double& nrml_j,
      /* out */ double& nrml_k
    )
    throw ( 
      ::iBase::Error
    );


    /**
     * Return the normals at specified uv positions on gfaces.  If any
     * gentity input is not a face, returns error.  Input parameters and 
     * output normals are interleaved.
     * @param gface_handles The entities being queried
     * @param parameters The uv parameters of points being queried, interleaved
     * @param normals Normals at specified points, interleaved
     */
    void
    getArrNrmlUV (
      /* in */ ::sidl::array<void*> gface_handles,
      /* in */ int32_t gface_handles_size,
      /* in */ ::iBase::StorageOrder storage_order,
      /* in */ ::sidl::array<double> parameters,
      /* in */ int32_t parameters_size,
      /* inout */ ::sidl::array<double>& normals,
      /* in */ int32_t normals_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getEntTgntU (
      /* in */ void* entity_handle,
      /* in */ double param_coord,
      /* out */ double& tngt_i,
      /* out */ double& tngt_j,
      /* out */ double& tngt_k
    )
    throw ( 
      ::iBase::Error
    );


    /**
     * Return the tangents at specified u positions on gedges.  If any
     * gentity input is not a face, returns error.  Output normals are 
     * interleaved.
     * @param gentity_handles The gedges being queried
     * @param parameters The u parameters of points being queried
     * @param tangents Tangents at specified points, interleaved
     */
    void
    getArrTgntU (
      /* in */ ::sidl::array<void*> gedge_handles,
      /* in */ int32_t gedge_handles_size,
      /* in */ ::iBase::StorageOrder storage_order,
      /* in */ ::sidl::array<double> parameters,
      /* in */ int32_t parameters_size,
      /* inout */ ::sidl::array<double>& tangents,
      /* inout */ int32_t& tangents_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getEnt1stDrvt (
      /* in */ void* entity_handle,
      /* in */ double u,
      /* in */ double v,
      /* inout */ ::sidl::array<double>& drvt_u,
      /* inout */ int32_t& drvt_u_size,
      /* inout */ ::sidl::array<double>& drvt_v,
      /* inout */ int32_t& drvt_v_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getArr1stDrvt (
      /* in */ ::sidl::array<void*> entity_handles,
      /* in */ int32_t entity_handles_size,
      /* in */ ::iBase::StorageOrder storage_order,
      /* in */ ::sidl::array<double> uv,
      /* in */ int32_t uv_size,
      /* inout */ ::sidl::array<double>& drvt_u,
      /* inout */ int32_t& drvt_u_size,
      /* inout */ ::sidl::array<int32_t>& u_offset,
      /* inout */ int32_t& u_offset_size,
      /* inout */ ::sidl::array<double>& drvt_v,
      /* inout */ int32_t& drvt_v_size,
      /* inout */ ::sidl::array<int32_t>& v_offset,
      /* inout */ int32_t& v_offset_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getEnt2ndDrvt (
      /* in */ void* entity_handle,
      /* in */ double u,
      /* in */ double v,
      /* inout */ ::sidl::array<double>& drvt_uu,
      /* inout */ int32_t& drvt_uu_size,
      /* inout */ ::sidl::array<double>& drvt_vv,
      /* inout */ int32_t& drvt_vv_size,
      /* inout */ ::sidl::array<double>& drvt_uv,
      /* inout */ int32_t& drvt_uv_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getArr2ndDrvt (
      /* in */ ::sidl::array<void*> entity_handles,
      /* in */ int32_t entity_handles_size,
      /* in */ ::iBase::StorageOrder storage_order,
      /* in */ ::sidl::array<double> uv,
      /* in */ int32_t uv_size,
      /* inout */ ::sidl::array<double>& drvt_uu,
      /* inout */ int32_t& drvt_uu_size,
      /* inout */ ::sidl::array<int32_t>& uu_offset,
      /* inout */ int32_t& uu_offset_size,
      /* inout */ ::sidl::array<double>& drvt_vv,
      /* inout */ int32_t& drvt_vv_size,
      /* inout */ ::sidl::array<int32_t>& vv_offset,
      /* inout */ int32_t& vv_offset_size,
      /* inout */ ::sidl::array<double>& drvt_uv,
      /* inout */ int32_t& drvt_uv_size,
      /* inout */ ::sidl::array<int32_t>& uv_offset,
      /* inout */ int32_t& uv_offset_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getFcCvtrUV (
      /* in */ void* entity_handle,
      /* in */ double u,
      /* in */ double v,
      /* out */ double& cvtr1_i,
      /* out */ double& cvtr1_j,
      /* out */ double& cvtr1_k,
      /* out */ double& cvtr2_i,
      /* out */ double& cvtr2_j,
      /* out */ double& cvtr2_k
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getFcArrCvtrUV (
      /* in */ ::sidl::array<void*> face_handles,
      /* in */ int32_t face_handles_size,
      /* in */ ::iBase::StorageOrder storage_order,
      /* in */ ::sidl::array<double> uv,
      /* in */ int32_t uv_size,
      /* inout */ ::sidl::array<double>& cvtr_1,
      /* inout */ int32_t& cvtr_1_size,
      /* inout */ ::sidl::array<double>& cvtr_2,
      /* inout */ int32_t& cvtr_2_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    isEntPeriodic (
      /* in */ void* entity_handle,
      /* out */ int32_t& in_u,
      /* out */ int32_t& in_v
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    isArrPeriodic (
      /* in */ ::sidl::array<void*> entity_handles,
      /* in */ int32_t entity_handles_size,
      /* inout */ ::sidl::array<int32_t>& in_uv,
      /* inout */ int32_t& in_uv_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    int32_t
    isFcDegenerate (
      /* in */ void* face_handle
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    isFcArrDegenerate (
      /* in */ ::sidl::array<void*> face_handles,
      /* in */ int32_t face_handles_size,
      /* inout */ ::sidl::array<int32_t>& degenerate,
      /* inout */ int32_t& degenerate_size
    )
    throw ( 
      ::iBase::Error
    );


    /**
     * Return the relative and absolute tolerances at the modeler level.  If
     * model does not have a modeler-wide tolerance, zero is returned for both
     * values.
     * @param relative_tolerance Relative tolerance for model as a whole
     * @param absolute_tolerance Absolute tolerance for model as a whole
     */
    void
    getTolerance (
      /* out */ int32_t& type,
      /* out */ double& tolerance
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    double
    getEntTolerance (
      /* in */ void* entity_handle
    )
    throw ( 
      ::iBase::Error
    );


    /**
     * Return the relative and absolute tolerances for specified gentities.  If
     * a gentity does not have a specific tolerance, zero is returned for both
     * values.
     * @param gentity_handles Gentities being queried
     * @param relative_tolerances Relative tolerances
     * @param absolute_tolerances Absolute tolerances
     */
    void
    getArrTolerance (
      /* in */ ::sidl::array<void*> entity_handles,
      /* in */ int32_t entity_handles_size,
      /* inout */ ::sidl::array<double>& tolerances,
      /* inout */ int32_t& tolerances_size
    )
    throw ( 
      ::iBase::Error
    );


    /**
     * Initialize an iterator over gentities of a specified dimension.
     * @param gentity_dimension Dimension of gentities to be iterated over
     * @param gentity_iterator Iterator initialized by this function
     */
    void
    initEntIter (
      /* in */ int32_t gentity_dimension,
      /* out */ void*& gentity_iterator
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    bool
    initEntArrIter (
      /* in */ void* entity_set_handle,
      /* in */ ::iBase::EntityType requested_entity_type,
      /* in */ int32_t requested_array_size,
      /* out */ void*& entArr_iterator
    )
    throw ( 
      ::iBase::Error
    );


    /**
     * Get the next entity for this iterator.
     * @param gentity_iterator Iterator being iterated over
     * @param gentity_handle Next gentity
     * @return If true, there are more gentities, if false, this is the last one
     */
    bool
    getNextEntIter (
      /* inout */ void*& gentity_iterator,
      /* out */ void*& gentity_handle
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    bool
    getNextEntArrIter (
      /* in */ void* entArr_iterator,
      /* inout */ ::sidl::array<void*>& entity_handles,
      /* inout */ int32_t& entity_handles_size
    )
    throw ( 
      ::iBase::Error
    );


    /**
     * Reset an iterator back to the first gentity
     * @param gentity_iterator Iterator reset by this function
     */
    void
    resetEntIter (
      /* in */ void* gentity_iterator
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    resetEntArrIter (
      /* in */ void* entArr_iterator
    )
    throw ( 
      ::iBase::Error
    );


    /**
     * Delete an iterator
     * @param gentity_iterator Iterator deleted by this function
     */
    void
    endEntIter (
      /* in */ void* Gentity_dim_iterator
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    endEntArrIter (
      /* in */ void* entArr_iterator
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    Copy (
      /* in */ void* geom_entity,
      /* out */ void*& geom_entity2
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    SweepAboutAxis (
      /* in */ void* geom_entity,
      /* in */ double angle,
      /* in */ double axis_normal_x,
      /* in */ double axis_normal_y,
      /* in */ double axis_normal_z,
      /* out */ void*& geom_entity2
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    Delete (
      /* in */ void* geom_entity
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    Sphere (
      /* in */ double radius,
      /* out */ void*& geom_entity
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    Brick (
      /* in */ double x,
      /* in */ double y,
      /* in */ double z,
      /* out */ void*& geom_entity
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    Cylinder (
      /* in */ double height,
      /* in */ double major_rad,
      /* in */ double minor_rad,
      /* out */ void*& geom_entity
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    Torus (
      /* in */ double major_rad,
      /* in */ double minor_rad,
      /* out */ void*& geom_entity
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    Move (
      /* inout */ void*& geom_entity,
      /* in */ double x,
      /* in */ double y,
      /* in */ double z
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    Rotate (
      /* inout */ void*& geom_entity,
      /* in */ double angle,
      /* in */ double axis_normal_x,
      /* in */ double axis_normal_y,
      /* in */ double axis_normal_z
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    Reflect (
      /* inout */ void*& geom_entity,
      /* in */ double plane_normal_x,
      /* in */ double plane_normal_y,
      /* in */ double plane_normal_z
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    Scale (
      /* inout */ void*& geom_entity,
      /* in */ double scale_x,
      /* in */ double scale_y,
      /* in */ double scale_z
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    Unite (
      /* in */ ::sidl::array<void*> geom_entities,
      /* in */ int32_t geom_entities_size,
      /* out */ void*& geom_entity
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    Subtract (
      /* in */ void* blank,
      /* in */ void* tool,
      /* out */ void*& geom_entity
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    Section (
      /* inout */ void*& geom_entity,
      /* in */ double plane_normal_x,
      /* in */ double plane_normal_y,
      /* in */ double plane_normal_z,
      /* in */ double offset,
      /* in */ bool reverse,
      /* out */ void*& geom_entity2
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    Imprint (
      /* in */ ::sidl::array<void*> geom_entities,
      /* in */ int32_t geom_entities_size
    )
    throw ( 
      ::iBase::Error
    );

    /**
     * user defined non-static method.
     */
    void
    Merge (
      /* in */ ::sidl::array<void*> geom_entities,
      /* in */ int32_t geom_entities_size,
      /* in */ double tolerance
    )
    throw ( 
      ::iBase::Error
    );

  };  // end class GeomSidl_impl

} // end namespace iGeom_SIDL

// DO-NOT-DELETE splicer.begin(iGeom_SIDL.GeomSidl._misc)
// Insert-Code-Here {iGeom_SIDL.GeomSidl._misc} (miscellaneous things)
// DO-NOT-DELETE splicer.end(iGeom_SIDL.GeomSidl._misc)

#endif
