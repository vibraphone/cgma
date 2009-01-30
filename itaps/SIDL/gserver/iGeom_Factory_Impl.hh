// 
// File:          iGeom_Factory_Impl.hh
// Symbol:        iGeom.Factory-v0.6.99
// Symbol Type:   class
// Babel Version: 0.10.10
// sidl Created:  20090126 14:50:17 CST
// Generated:     20090126 14:50:21 CST
// Description:   Server-side implementation for iGeom.Factory
// 
// WARNING: Automatically generated; only changes within splicers preserved
// 
// babel-version = 0.10.10
// source-line   = 1115
// source-url    = file:/home/jason/meshkit/cgm/itaps/SIDL/iGeom.sidl
// xml-url       = /home/jason/meshkit/cgm/itaps/SIDL/repo/iGeom.Factory-v0.6.99.xml
// 

#ifndef included_iGeom_Factory_Impl_hh
#define included_iGeom_Factory_Impl_hh

#ifndef included_sidl_cxx_hh
#include "sidl_cxx.hh"
#endif
#ifndef included_iGeom_Factory_IOR_h
#include "iGeom_Factory_IOR.h"
#endif
// 
// Includes for all method dependencies.
// 
#ifndef included_iBase_Error_hh
#include "iBase_Error.hh"
#endif
#ifndef included_iGeom_Factory_hh
#include "iGeom_Factory.hh"
#endif
#ifndef included_iGeom_Geometry_hh
#include "iGeom_Geometry.hh"
#endif
#ifndef included_sidl_BaseInterface_hh
#include "sidl_BaseInterface.hh"
#endif
#ifndef included_sidl_ClassInfo_hh
#include "sidl_ClassInfo.hh"
#endif


// DO-NOT-DELETE splicer.begin(iGeom.Factory._includes)
// Insert-Code-Here {iGeom.Factory._includes} (includes or arbitrary code)
// DO-NOT-DELETE splicer.end(iGeom.Factory._includes)

namespace iGeom { 

  /**
   * Symbol "iGeom.Factory" (version 0.6.99)
   */
  class Factory_impl
  // DO-NOT-DELETE splicer.begin(iGeom.Factory._inherits)
  // Insert-Code-Here {iGeom.Factory._inherits} (optional inheritance here)
  // DO-NOT-DELETE splicer.end(iGeom.Factory._inherits)
  {

  private:
    // Pointer back to IOR.
    // Use this to dispatch back through IOR vtable.
    Factory self;

    // DO-NOT-DELETE splicer.begin(iGeom.Factory._implementation)
    // Insert-Code-Here {iGeom.Factory._implementation} (additional details)
    // DO-NOT-DELETE splicer.end(iGeom.Factory._implementation)

  private:
    // private default constructor (required)
    Factory_impl() 
    {} 

  public:
    // sidl constructor (required)
    // Note: alternate Skel constructor doesn't call addref()
    // (fixes bug #275)
    Factory_impl( struct iGeom_Factory__object * s ) : self(s,true) { _ctor(); }

    // user defined construction
    void _ctor();

    // virtual destructor (required)
    virtual ~Factory_impl() { _dtor(); }

    // user defined destruction
    void _dtor();

    // static class initializer
    static void _load();

  public:
    /**
     * user defined static method.
     */
    static ::iGeom::Geometry
    newGeom (
      /* in */ const ::std::string& options
    )
    throw ( 
      ::iBase::Error
    );


  };  // end class Factory_impl

} // end namespace iGeom

// DO-NOT-DELETE splicer.begin(iGeom.Factory._misc)
// Insert-Code-Here {iGeom.Factory._misc} (miscellaneous things)
// DO-NOT-DELETE splicer.end(iGeom.Factory._misc)

#endif
