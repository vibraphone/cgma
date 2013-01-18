//- Class: CubitObserver
//- Description: Observer class - *very* rudimentary implementation
//- Owner:  Tim Tautges

#ifndef CUBITOBSERVER_HPP
#define CUBITOBSERVER_HPP

#include "CubitDefines.h"
#include "CubitEntity.hpp"
#include "CubitEvent.hpp"
#include "DLIList.hpp"
#include "CubitUtilConfigure.h"

class CubitObservable;

class CUBIT_UTIL_EXPORT CubitObserver
{
public:

  CubitObserver();
  
  virtual ~CubitObserver();

  CubitStatus watch_observable(CubitObservable *observable)
    { return register_observable(observable); }
    //- Cause this observer to watch the specified observable.
  
    //- Request that this observer stop watching the specified observable
  CubitStatus stop_watching_observable(CubitObservable *observable)
    { return unregister_observable(observable); }
  
   
  virtual CubitStatus notify_observer(
     CubitObservable *observable, const CubitEvent &observer_event) = 0;
    // This should be called through the
    // CubitObservable::notify_observers generic notify function.
    // NOTE: VERY IMPORTANT: when implementing this function, leaf
    // classes should call register_observable() or
    // unregister_observable() if the event is one of construction or
    // destruction, resp.; this makes sure the observable knows we're
    // observing it, or unhooks us from the observable.  A count of
    // observables being observed is also kept, to prevent stale
    // pointers on observables

  CubitStatus register_observable(CubitObservable *observable);
    //- Cause this observer to watch the specified observable.
  
  CubitStatus unregister_observable(CubitObservable *observable);
    //- remove this observer from the observable's list, and if successful,
    //- decrement the observable count  

  //! register observer with the global CubitEventDispatcher
  static CubitStatus register_observer(CubitObserver* obs);
  //! unregister observer with the global CubitEventDispatcher
  static CubitStatus unregister_observer(CubitObserver* obs);

private:

  int observableCount;
    //- number of observables to which this observer is registered
}; // 

#endif

