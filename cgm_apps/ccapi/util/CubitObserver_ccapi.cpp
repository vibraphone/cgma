#include "CubitObserver_ccapi.h"

#include "CubitObserver.hpp"
#include "CubitObservable.hpp"

#include "CubitEvent.hpp"

void *CubitObserver_get_address(void *this_observer, enum EntityType entity_type)
{
  CubitObserver *temp_observer = (CubitObserver *) this_observer;
  return temp_observer->get_address(entity_type);
}
        
enum CubitStatus CubitObserver_notify_observer(void *this_observer,
                                                 /* CubitObservable * */ void *observable,
                                                 /* const CubitEvent & */ void *observer_event,
                                               enum CubitBoolean from_observable)
{
  CubitObserver *temp_observer = (CubitObserver *) this_observer;
  CubitObservable *temp_observable = (CubitObservable *) observable;
  CubitEvent *temp_observer_event = (CubitEvent *) observer_event;
  
  return temp_observer->notify_observer(temp_observable, *temp_observer_event, from_observable);
}

enum CubitStatus CubitObserver_notify_static_observers(void *this_observer,
                                                                /* CubitObservable * */ void *observable,
                                                                /* const CubitEvent & */ void *observer_event,
                                                              enum CubitBoolean from_observable)
{
  CubitObserver *temp_observer = (CubitObserver *) this_observer;
  CubitObservable *temp_observable = (CubitObservable *) observable;
  CubitEvent *temp_observer_event = (CubitEvent *) observer_event;

  return temp_observer->notify_static_observers(temp_observable, *temp_observer_event, from_observable);
}

void CubitObserver_init_static_observers(void *)
{
  CubitObserver::init_static_observers();
}

void CubitObserver_term_static_observers(void *)
{
  CubitObserver::term_static_observers();
}

