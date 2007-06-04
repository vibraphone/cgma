// 
// File:          TSTTG_CGM_CgmGeom_Impl.hh
// Symbol:        TSTTG_CGM.CgmGeom-v0.1
// Symbol Type:   class
// Babel Version: 0.9.8
// sidl Created:  20051011 16:26:00 GMT-06:00
// Generated:     20051011 16:26:03 GMT-06:00
// Description:   Server-side implementation for TSTTG_CGM.CgmGeom
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.9.8
// source-line   = 5
// source-url    = file:/home/tjtautg/tstt/cvs/TSTTG/TSTTG_CGM.sidl
// 

#ifndef included_TSTTG_CGM_CgmGeom_Impl_hh
#define included_TSTTG_CGM_CgmGeom_Impl_hh

#ifndef included_sidl_cxx_hh
#include "sidl_cxx.hh"
#endif
#ifndef included_TSTTG_CGM_CgmGeom_IOR_h
#include "TSTTG_CGM_CgmGeom_IOR.h"
#endif
// 
// Includes for all method dependencies.
// 
#ifndef included_TSTTB_Error_hh
#include "TSTTB_Error.hh"
#endif
#ifndef included_TSTTB_TagValueType_hh
#include "TSTTB_TagValueType.hh"
#endif
#ifndef included_TSTTG_GentityType_hh
#include "TSTTG_GentityType.hh"
#endif
#ifndef included_TSTTG_CGM_CgmGeom_hh
#include "TSTTG_CGM_CgmGeom.hh"
#endif
#ifndef included_sidl_BaseInterface_hh
#include "sidl_BaseInterface.hh"
#endif
#ifndef included_sidl_ClassInfo_hh
#include "sidl_ClassInfo.hh"
#endif


// DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom._includes)
// Put additional includes or other arbitrary code here...
// DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom._includes)

namespace TSTTG_CGM { 

  /**
   * Symbol "TSTTG_CGM.CgmGeom" (version 0.1)
   */
  class CgmGeom_impl
  // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom._inherits)
  // Put additional inheritance here...
  // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom._inherits)
  {

  private:
    // Pointer back to IOR.
    // Use this to dispatch back through IOR vtable.
    CgmGeom self;

    // DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom._implementation)
    // Put additional implementation details here...
    void processError() throw(::TSTTB::Error);
    // DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom._implementation)

  private:
    // private default constructor (required)
    CgmGeom_impl() {} 

  public:
    // sidl constructor (required)
    // Note: alternate Skel constructor doesn't call addref()
    // (fixes bug #275)
    CgmGeom_impl( struct TSTTG_CGM_CgmGeom__object * s ) : self(s,
      true) { _ctor(); }

    // user defined construction
    void _ctor();

    // virtual destructor (required)
    virtual ~CgmGeom_impl() { _dtor(); }

    // user defined destruction
    void _dtor();

  public:


    /**
     * Return gentities of specified dimension in this set, or in whole model.
     * @param set_handle Entity set being queried (if 0, whole model)
     * @param gentity_dimension Dimension of entities being queried
     * @param gentity_handles Gentity handles
     */
    void
    gentitysetGetGentitiesOfType (
      /*in*/ void* set_handle,
      /*in*/ ::TSTTG::GentityType gentity_type,
      /*inout*/ ::sidl::array<void*>& gentity_handles,
      /*out*/ int32_t& gentity_handles_size
    )
    throw ( 
      ::TSTTB::Error
    );


    /**
     * Return number of gentities of specified dimension in this set, or in
     * whole model.
     * @param set_handle Entity set being queried (if 0, whole model)
     * @param gentity_dimension Dimension of entities being queried
     * @return Number of entities
     */
    int32_t
    gentitysetGetNumberGentitiesOfType (
      /*in*/ void* set_handle,
      /*in*/ ::TSTTG::GentityType gentity_type
    )
    throw ( 
      ::TSTTB::Error
    );


    /**
     *    Returns an integer array of topological dimensions for an input
     *    array of entity handles.
     */
    void
    gentityGetType (
      /*in*/ ::sidl::array<void*> gentity_handles,
      /*in*/ int32_t gentity_handles_size,
      /*inout*/ ::sidl::array< ::TSTTG::GentityType>& gtype,
      /*inout*/ int32_t& gtype_size
    )
    throw ( 
      ::TSTTB::Error
    );


    /**
     * Get the adjacent entities of a given dimension.
     * @param gentity_handle Entity for which adjacencies are requested
     * @param to_dimension Target dimension of adjacent entities
     * @param adj_gentities List returned with adjacent entities
     */
    void
    gentityGetAdjacencies (
      /*in*/ void* gentity_handle,
      /*in*/ int32_t to_dimension,
      /*inout*/ ::sidl::array<void*>& adj_gentities,
      /*inout*/ int32_t& adj_gentities_size
    )
    throw ( 
      ::TSTTB::Error
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
    gentityGet2OAdjacencies (
      /*in*/ void* gentity_handle,
      /*in*/ int32_t bridge_dimension,
      /*in*/ int32_t to_dimension,
      /*inout*/ ::sidl::array<void*>& adjacent_gentities,
      /*out*/ int32_t& adjacent_gentities_size
    )
    throw ( 
      ::TSTTB::Error
    );


    /**
     * Return whether or not entities are adjacent.
     * @param gentity_handle1 1st entity
     * @param gentity_handle2 2nd entity
     * @param are_adjacent If true, entities are adjacent
     */
    void
    gentityIsAdjacent (
      /*in*/ void* gentity_handle1,
      /*in*/ void* gentity_handle2,
      /*out*/ bool& are_adjacent
    )
    throw ( 
      ::TSTTB::Error
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
    gLoad (
      /*in*/ const ::std::string& name,
      /*in*/ ::sidl::array< ::std::string> options,
      /*in*/ int32_t options_size
    )
    throw ( 
      ::TSTTB::Error
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
    gSave (
      /*in*/ const ::std::string& name,
      /*in*/ ::sidl::array< ::std::string> options,
      /*in*/ int32_t options_size
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    createTag (
      /*in*/ const ::std::string& tag_name,
      /*in*/ int32_t number_of_values,
      /*in*/ ::TSTTB::TagValueType tag_type,
      /*out*/ void*& tag_handle
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    destroyTag (
      /*in*/ void* tag_handle,
      /*in*/ bool forced
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    ::std::string
    getTagName (
      /*in*/ void* tag_handle
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    int32_t
    getTagSizeValues (
      /*in*/ void* tag_handle
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    int32_t
    getTagSizeBytes (
      /*in*/ void* tag_handle
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void*
    getTagHandle (
      /*in*/ const ::std::string& tag_name
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    ::TSTTB::TagValueType
    getTagType (
      /*in*/ void* tag_handle
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getArrData (
      /*in*/ ::sidl::array<void*> entity_handles,
      /*in*/ int32_t entity_handles_size,
      /*in*/ void* tag_handle,
      /*inout*/ ::sidl::array<char>& value_array,
      /*out*/ int32_t& value_array_size
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getIntArrData (
      /*in*/ ::sidl::array<void*> entity_handles,
      /*in*/ int32_t entity_handles_size,
      /*in*/ void* tag_handle,
      /*inout*/ ::sidl::array<int32_t>& value_array,
      /*out*/ int32_t& value_array_size
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getDblArrData (
      /*in*/ ::sidl::array<void*> entity_handles,
      /*in*/ int32_t entity_handles_size,
      /*in*/ void* tag_handle,
      /*inout*/ ::sidl::array<double>& value_array,
      /*out*/ int32_t& value_array_size
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getEHArrData (
      /*in*/ ::sidl::array<void*> entity_handles,
      /*in*/ int32_t entity_handles_size,
      /*in*/ void* tag_handle,
      /*inout*/ ::sidl::array<void*>& value_array,
      /*out*/ int32_t& value_array_size
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    setArrData (
      /*in*/ ::sidl::array<void*> entity_handles,
      /*in*/ int32_t entity_handles_size,
      /*in*/ void* tag_handle,
      /*in*/ ::sidl::array<char> value_array,
      /*in*/ int32_t value_array_size
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    setIntArrData (
      /*in*/ ::sidl::array<void*> entity_handles,
      /*in*/ int32_t entity_handles_size,
      /*in*/ void* tag_handle,
      /*in*/ ::sidl::array<int32_t> value_array,
      /*in*/ int32_t value_array_size
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    setDblArrData (
      /*in*/ ::sidl::array<void*> entity_handles,
      /*in*/ int32_t entity_handles_size,
      /*in*/ void* tag_handle,
      /*in*/ ::sidl::array<double> value_array,
      /*in*/ int32_t value_array_size
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    setEHArrData (
      /*in*/ ::sidl::array<void*> entity_handles,
      /*in*/ int32_t entity_handles_size,
      /*in*/ void* tag_handle,
      /*in*/ ::sidl::array<void*> value_array,
      /*in*/ int32_t value_array_size
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    rmvArrTag (
      /*in*/ ::sidl::array<void*> entity_handles,
      /*in*/ int32_t entity_handles_size,
      /*in*/ void* tag_handle
    )
    throw ( 
      ::TSTTB::Error
    );


    /**
     * Initialize an iterator over gentities of a specified dimension.
     * @param gentity_dimension Dimension of gentities to be iterated over
     * @param gentity_iterator Iterator initialized by this function
     */
    void
    gentityIteratorInit (
      /*in*/ int32_t gentity_dimension,
      /*out*/ void*& gentity_iterator
    )
    throw ( 
      ::TSTTB::Error
    );


    /**
     * Get the next entity for this iterator.
     * @param gentity_iterator Iterator being iterated over
     * @param gentity_handle Next gentity
     * @return If true, there are more gentities, if false, this is the last one
     */
    bool
    gentityIteratorNext (
      /*inout*/ void*& gentity_iterator,
      /*out*/ void*& gentity_handle
    )
    throw ( 
      ::TSTTB::Error
    );


    /**
     * Reset an iterator back to the first gentity
     * @param gentity_iterator Iterator reset by this function
     */
    void
    gentityIteratorReset (
      /*inout*/ void*& gentity_iterator
    )
    throw ( 
      ::TSTTB::Error
    );


    /**
     * Delete an iterator
     * @param gentity_iterator Iterator deleted by this function
     */
    void
    gentityIteratorDelete (
      /*in*/ void* Gentity_dim_iterator
    )
    throw ( 
      ::TSTTB::Error
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
    gentityClosestPoint (
      /*in*/ ::sidl::array<void*> gentity_handles,
      /*in*/ int32_t gentity_handles_size,
      /*in*/ ::sidl::array<double> near_coordinates,
      /*in*/ int32_t near_coordinates_size,
      /*inout*/ ::sidl::array<double>& on_coordinates,
      /*out*/ int32_t& on_coordinates_size
    )
    throw ( 
      ::TSTTB::Error
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
    gentityNormal (
      /*in*/ ::sidl::array<void*> gentity_handles,
      /*in*/ int32_t gentity_handles_size,
      /*in*/ ::sidl::array<double> coordinates,
      /*in*/ int32_t coordinates_size,
      /*inout*/ ::sidl::array<double>& normals,
      /*out*/ int32_t& normals_size
    )
    throw ( 
      ::TSTTB::Error
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
    gentityClosestPointAndNormal (
      /*in*/ ::sidl::array<void*> gentity_handles,
      /*in*/ int32_t gentity_handles_size,
      /*in*/ ::sidl::array<double> near_coordinates,
      /*in*/ int32_t near_coordinates_size,
      /*inout*/ ::sidl::array<double>& on_coordinates,
      /*out*/ int32_t& on_coordinates_size,
      /*inout*/ ::sidl::array<double>& normals,
      /*out*/ int32_t& normals_size
    )
    throw ( 
      ::TSTTB::Error
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
    gentityTangent (
      /*in*/ ::sidl::array<void*> gentity_handles,
      /*in*/ int32_t gentity_handles_size,
      /*in*/ ::sidl::array<double> coordinates,
      /*in*/ int32_t coordinates_size,
      /*inout*/ ::sidl::array<double>& tangents,
      /*out*/ int32_t& tangents_size
    )
    throw ( 
      ::TSTTB::Error
    );


    /**
     * Return the bounding boxex of given entities; coordinates returned
     * interleaved.
     * @param gentity_handles The gentities being queried
     * @param min_corners Minimum corner coordinates of the boxes, interleaved
     * @param max_corners Maximum corner coordinates of the boxes, interleaved
     */
    void
    gentityBoundingBox (
      /*in*/ ::sidl::array<void*> gentity_handles,
      /*in*/ int32_t gentity_handles_size,
      /*inout*/ ::sidl::array<double>& min_corner,
      /*out*/ int32_t& min_corner_size,
      /*inout*/ ::sidl::array<double>& max_corner,
      /*out*/ int32_t& max_corner_size
    )
    throw ( 
      ::TSTTB::Error
    );


    /**
     * Return the coordinates of the specified vertices; returns error if any
     * of the entities are not gvertices.  Coordinates returned interleaved.
     * @param gentity_handles The gentities being queried
     * @param coordinates The coordinates of the gvertices, interleaved.
     */
    void
    getGvertexCoordinates (
      /*in*/ ::sidl::array<void*> gentity_handles,
      /*in*/ int32_t gentity_handles_size,
      /*inout*/ ::sidl::array<double>& coordinates,
      /*out*/ int32_t& coordinates_size
    )
    throw ( 
      ::TSTTB::Error
    );


    /**
     * Return the sense of a gface with respect to a gregion.  Sense is either
     * forward (=1), reverse (=-1), both (=2), or unknown (=0).  Error is returned
     * if first entity is not a gface or second entity is not a gregion.
     * @param gface Gface whose sense is being queried.
     * @param gregion Gregion gface is being queried with respect to
     */
    int32_t
    getGnormalSense (
      /*in*/ void* gface,
      /*in*/ void* gregion
    )
    throw ( 
      ::TSTTB::Error
    );


    /**
     * Return the sense of a gedge with respect to a gface.  Sense is either
     * forward (=1), reverse (=-1), both (=2), or unknown (=0).  Error is returned
     * if first entity is not a gedge or second entity is not a gface.
     * @param gedge Gedge whose sense is being queried.
     * @param gface Gface gedge is being queried with respect to
     */
    int32_t
    getGtangentSense (
      /*in*/ void* gedge,
      /*in*/ void* gface
    )
    throw ( 
      ::TSTTB::Error
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
    getGvertexTangentSense (
      /*in*/ void* gedge,
      /*in*/ void* gvertex1,
      /*in*/ void* gvertex2
    )
    throw ( 
      ::TSTTB::Error
    );


    /**
     * Return whether a given gentity is parametric or not.  If a gentity
     * is not parametric, all of the following functions will return an error
     * when called on that entity.
     * @param gentity_handle Gentity being queried.
     */
    int32_t
    gentityIsParametric (
      /*in*/ void* gentity_handle
    )
    throw ( 
      ::TSTTB::Error
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
    gentityUvToXyz (
      /*in*/ ::sidl::array<void*> gentity_handles,
      /*in*/ int32_t gentity_handles_size,
      /*in*/ ::sidl::array<double> uv,
      /*in*/ int32_t uv_size,
      /*inout*/ ::sidl::array<double>& coordinates,
      /*out*/ int32_t& coordinates_size
    )
    throw ( 
      ::TSTTB::Error
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
    gentityXyzToUv (
      /*in*/ ::sidl::array<void*> gentity_handles,
      /*in*/ int32_t gentity_handles_size,
      /*in*/ ::sidl::array<double> coordinates,
      /*in*/ int32_t coordinates_size,
      /*inout*/ ::sidl::array<double>& uv,
      /*out*/ int32_t& uv_size
    )
    throw ( 
      ::TSTTB::Error
    );


    /**
     * Return the uv range of the specified gentities.  Parameters are interleaved.
     * @param gentity_handles Gentities being queried.
     * @param uv_min Minimum parameters of gentities, interleaved
     * @param uv_max Maximum parameters of gentities, interleaved
     */
    void
    gentityUvRange (
      /*in*/ ::sidl::array<void*> gentity_handles,
      /*in*/ int32_t gentity_handles_size,
      /*inout*/ ::sidl::array<double>& uv_min,
      /*out*/ int32_t& uv_min_size,
      /*inout*/ ::sidl::array<double>& uv_max,
      /*out*/ int32_t& uv_max_size
    )
    throw ( 
      ::TSTTB::Error
    );


    /**
     * Given source gentities, parametric positions on those gentities, and 
     * bounding gentities, return the parametric positions on the bounding
     * gentities.  If a source gentity is a gvertex, parametric positions for
     * that entry are ignored.  In cases where multiple source entities are
     * input, two input and output parameters per gentity is assumed, even if
     * input consists only of gedges or gvertices.
     * @param src_gentity_handles Source gentities
     * @param src_uv Uv positions on source gentities
     * @param trg_gentity_handles Target gentities
     * @param trg_uv Uv positions on target gentities
     */
    void
    Greparam_edge_face (
      /*in*/ ::sidl::array<void*> src_gentity_handles,
      /*in*/ int32_t src_gentity_handles_size,
      /*in*/ ::sidl::array<double> src_uv,
      /*in*/ int32_t src_uv_size,
      /*in*/ ::sidl::array<void*> trg_gentity_handles,
      /*in*/ int32_t trg_gentity_handles_size,
      /*in*/ ::sidl::array<double> trg_uv,
      /*in*/ int32_t trg_uv_size
    )
    throw ( 
      ::TSTTB::Error
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
    gentityNormalUv (
      /*in*/ ::sidl::array<void*> gface_handles,
      /*in*/ int32_t gface_handles_size,
      /*in*/ ::sidl::array<double> parameters,
      /*in*/ int32_t parameters_size,
      /*inout*/ ::sidl::array<double>& normals,
      /*in*/ int32_t normals_size
    )
    throw ( 
      ::TSTTB::Error
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
    gentityTangentU (
      /*in*/ ::sidl::array<void*> gedge_handles,
      /*in*/ int32_t gedge_handles_size,
      /*in*/ ::sidl::array<double> parameters,
      /*in*/ int32_t parameters_size,
      /*inout*/ ::sidl::array<double>& tangents,
      /*out*/ int32_t& tangents_size
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getData (
      /*in*/ void* entity_handle,
      /*in*/ void* tag_handle,
      /*inout*/ ::sidl::array<char>& tag_value,
      /*out*/ int32_t& tag_value_size
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    int32_t
    getIntData (
      /*in*/ void* entity_handle,
      /*in*/ void* tag_handle
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    double
    getDblData (
      /*in*/ void* entity_handle,
      /*in*/ void* tag_handle
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void*
    getEHData (
      /*in*/ void* entity_handle,
      /*in*/ void* tag_handle
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    setData (
      /*in*/ void* entity_handle,
      /*in*/ void* tag_handle,
      /*in*/ ::sidl::array<char> tag_value,
      /*in*/ int32_t tag_value_size
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    setIntData (
      /*in*/ void* entity_handle,
      /*in*/ void* tag_handle,
      /*in*/ int32_t tag_value
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    setDblData (
      /*in*/ void* entity_handle,
      /*in*/ void* tag_handle,
      /*in*/ double tag_value
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    setEHData (
      /*in*/ void* entity_handle,
      /*in*/ void* tag_handle,
      /*in*/ void* tag_value
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getAllTags (
      /*in*/ void* entity_handle,
      /*inout*/ ::sidl::array<void*>& tag_handles,
      /*out*/ int32_t& tag_handles_size
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    rmvTag (
      /*in*/ void* entity_handle,
      /*in*/ void* tag_handle
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    createEntSet (
      /*in*/ bool isList,
      /*out*/ void*& entity_set
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    destroyEntSet (
      /*in*/ void* entity_set
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    bool
    isList (
      /*in*/ void* entity_set
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    int32_t
    getNumEntSets (
      /*in*/ void* entity_set,
      /*in*/ int32_t num_hops
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getEntSets (
      /*in*/ void* entity_set,
      /*in*/ int32_t num_hops,
      /*inout*/ ::sidl::array<void*>& contained_entset_handles,
      /*out*/ int32_t& contained_entset_handles_size
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    addEntToSet (
      /*in*/ void* entity_handle,
      /*inout*/ void*& entity_set
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    rmvEntFromSet (
      /*in*/ void* entity_handle,
      /*inout*/ void*& entity_set
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    addEntArrToSet (
      /*in*/ ::sidl::array<void*> entity_handles,
      /*in*/ int32_t entity_handles_size,
      /*inout*/ void*& entity_set
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    rmvEntArrFromSet (
      /*in*/ ::sidl::array<void*> entity_handles,
      /*in*/ int32_t entity_handles_size,
      /*inout*/ void*& entity_set
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    addEntSet (
      /*in*/ void* entity_set_to_add,
      /*inout*/ void*& entity_set_handle
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    rmvEntSet (
      /*in*/ void* entity_set_to_remove,
      /*inout*/ void*& entity_set_handle
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    bool
    isEntContained (
      /*in*/ void* containing_entity_set,
      /*in*/ void* entity_handle
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    bool
    isEntSetContained (
      /*in*/ void* containing_entity_set,
      /*in*/ void* contained_entity_set
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    subtract (
      /*in*/ void* entity_set_1,
      /*in*/ void* entity_set_2,
      /*out*/ void*& result_entity_set
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    intersect (
      /*in*/ void* entity_set_1,
      /*in*/ void* entity_set_2,
      /*out*/ void*& result_entity_set
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    unite (
      /*in*/ void* entity_set_1,
      /*in*/ void* entity_set_2,
      /*out*/ void*& result_entity_set
    )
    throw ( 
      ::TSTTB::Error
    );


    /**
     * Return the relative and absolute tolerances at the modeler level.  If
     * model does not have a modeler-wide tolerance, zero is returned for both
     * values.
     * @param relative_tolerance Relative tolerance for model as a whole
     * @param absolute_tolerance Absolute tolerance for model as a whole
     */
    void
    getGtolerance (
      /*out*/ double& relative_tolerance,
      /*out*/ double& absolute_tolerance
    )
    throw ( 
      ::TSTTB::Error
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
    getGentityTolerance (
      /*in*/ ::sidl::array<void*> gentity_handles,
      /*in*/ int32_t gentity_handles_size,
      /*inout*/ ::sidl::array<double>& relative_tolerances,
      /*out*/ int32_t& relative_tolerances_size,
      /*inout*/ ::sidl::array<double>& absolute_tolerances,
      /*out*/ int32_t& absolute_tolerances_size
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    setEntSetData (
      /*in*/ void* entity_set,
      /*in*/ void* tag_handle,
      /*inout*/ ::sidl::array<char>& tag_value,
      /*in*/ int32_t tag_value_size
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    setEntSetIntData (
      /*in*/ void* entity_set,
      /*in*/ void* tag_handle,
      /*in*/ int32_t tag_value
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    setEntSetDblData (
      /*in*/ void* entity_set,
      /*in*/ void* tag_handle,
      /*in*/ double tag_value
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    setEntSetEHData (
      /*in*/ void* entity_set,
      /*in*/ void* tag_handle,
      /*in*/ void* tag_value
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getEntSetData (
      /*in*/ void* entity_set,
      /*in*/ void* tag_handle,
      /*inout*/ ::sidl::array<char>& tag_value,
      /*out*/ int32_t& tag_value_size
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    int32_t
    getEntSetIntData (
      /*in*/ void* entity_set,
      /*in*/ void* tag_handle
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    double
    getEntSetDblData (
      /*in*/ void* entity_set,
      /*in*/ void* tag_handle
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void*
    getEntSetEHData (
      /*in*/ void* entity_set,
      /*in*/ void* tag_handle
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getAllEntSetTags (
      /*in*/ void* entity_set,
      /*inout*/ ::sidl::array<void*>& tag_handles,
      /*out*/ int32_t& tag_handles_size
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    rmvEntSetTag (
      /*in*/ void* entity_set,
      /*in*/ void* tag_handle
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    Brick (
      /*in*/ double x,
      /*in*/ double y,
      /*in*/ double z,
      /*out*/ void*& geom_entity
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    Cylinder (
      /*in*/ double height,
      /*in*/ double major_rad,
      /*in*/ double minor_rad,
      /*out*/ void*& geom_entity
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    Torus (
      /*in*/ double major_rad,
      /*in*/ double minor_rad,
      /*out*/ void*& geom_entity
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    addPrntChld (
      /*inout*/ void*& parent_entity_set,
      /*inout*/ void*& child_entity_set
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    rmvPrntChld (
      /*inout*/ void*& parent_entity_set,
      /*inout*/ void*& child_entity_set
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    bool
    isChildOf (
      /*in*/ void* parent_entity_set,
      /*in*/ void* child_entity_set
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    int32_t
    getNumChld (
      /*in*/ void* entity_set,
      /*in*/ int32_t num_hops
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    int32_t
    getNumPrnt (
      /*in*/ void* entity_set,
      /*in*/ int32_t num_hops
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getChldn (
      /*in*/ void* from_entity_set,
      /*in*/ int32_t num_hops,
      /*inout*/ ::sidl::array<void*>& child_handles,
      /*out*/ int32_t& child_handles_size
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    getPrnts (
      /*in*/ void* from_entity_set,
      /*in*/ int32_t num_hops,
      /*inout*/ ::sidl::array<void*>& parent_handles,
      /*out*/ int32_t& parent_handles_size
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    Move (
      /*inout*/ void*& geom_entity,
      /*in*/ double x,
      /*in*/ double y,
      /*in*/ double z
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    Rotate (
      /*inout*/ void*& geom_entity,
      /*in*/ double angle,
      /*in*/ double axis_normal_x,
      /*in*/ double axis_normal_y,
      /*in*/ double axis_normal_z
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    Reflect (
      /*inout*/ void*& geom_entity,
      /*in*/ double plane_normal_x,
      /*in*/ double plane_normal_y,
      /*in*/ double plane_normal_z
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    Scale (
      /*inout*/ void*& geom_entity,
      /*in*/ double scale_x,
      /*in*/ double scale_y,
      /*in*/ double scale_z
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    Copy (
      /*in*/ void* geom_entity,
      /*out*/ void*& geom_entity2
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    SweepAboutAxis (
      /*in*/ void* geom_entity,
      /*in*/ double angle,
      /*in*/ double axis_normal_x,
      /*in*/ double axis_normal_y,
      /*in*/ double axis_normal_z,
      /*out*/ void*& geom_entity2
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    Delete (
      /*in*/ void* geom_entity
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    Unite (
      /*in*/ ::sidl::array<void*> geom_entities,
      /*in*/ int32_t geom_entities_size,
      /*out*/ void*& geom_entity
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    Subtract (
      /*in*/ void* blank,
      /*in*/ void* tool,
      /*out*/ void*& geom_entity
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    Section (
      /*inout*/ void*& geom_entity,
      /*in*/ double plane_normal_x,
      /*in*/ double plane_normal_y,
      /*in*/ double plane_normal_z,
      /*in*/ double offset,
      /*in*/ bool reverse,
      /*out*/ void*& geom_entity2
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    Imprint (
      /*in*/ ::sidl::array<void*> geom_entities,
      /*in*/ int32_t geom_entities_size
    )
    throw ( 
      ::TSTTB::Error
    );

    /**
     * user defined non-static method.
     */
    void
    Merge (
      /*in*/ ::sidl::array<void*> geom_entities,
      /*in*/ int32_t geom_entities_size,
      /*in*/ double tolerance
    )
    throw ( 
      ::TSTTB::Error
    );

  };  // end class CgmGeom_impl

} // end namespace TSTTG_CGM

// DO-NOT-DELETE splicer.begin(TSTTG_CGM.CgmGeom._misc)
// Put miscellaneous things here...
// DO-NOT-DELETE splicer.end(TSTTG_CGM.CgmGeom._misc)

#endif
