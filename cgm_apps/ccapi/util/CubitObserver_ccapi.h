#ifndef CUBITOBSERVER_CCAPI_H
#define CUBITOBSERVER_CCAPI_H

#include "CubitDefines.h"
#include "EntityType.h"

#ifdef __cplusplus
extern "C" {
#endif

  void *CubitObserver_get_address(void *this_observer, enum EntityType);
    /* - used for type casting */
        
  enum CubitStatus CubitObserver_notify_observer(void *this_observer,
                                                   /* CubitObservable * */ void *observable,
                                                   /* const CubitEvent & */ void *observer_event,
                                                 enum CubitBoolean from_observable);
    /* - generic notify function. */
    /* - NOTE: VERY IMPORTANT: when implementing this function, leaf classes should */
    /* - call register_observable() or unregister_observable() if the event is one */
    /* - of construction or destruction, resp.; this makes sure the observable */
    /* - knows we're observing it, or unhooks us from the observable.  A count of */
    /* - observables being observed is also kept, to prevent stale pointers on */
    /* - observables */

  enum CubitStatus CubitObserver_notify_static_observers(void *this_observer,
                                                                  /* CubitObservable * */ void *observable,
                                                                  /* const CubitEvent & */ void *observer_event,
                                                                enum CubitBoolean from_observable);
    /* - notify function for static observers; calls notify_observer for all static observers */

  void CubitObserver_init_static_observers(void *this_observer);
  void CubitObserver_term_static_observers(void *this_observer);
    /* - Access function to create and destroy the list of observed. */

#ifdef __cplusplus
}
#endif

#endif
