#include "CubitObserver.hpp"
#include "CubitObservable.hpp"
#include "CubitEvent.hpp"
#include "AppUtil.hpp"

CubitObserver::CubitObserver()
{
  observableCount = 0;
}

CubitObserver::~CubitObserver()
{
  assert(observableCount == 0);
}

CubitStatus CubitObserver::register_observable(CubitObservable *observable)
{
  if (observable == NULL)
    return CUBIT_FAILURE;

    //- add this observer to the observable's list, and if successful,
    //- increment the observable count
  CubitStatus success = observable->add_observer(this);
  
  if (success == CUBIT_SUCCESS)
    observableCount++;
  
  return success;
}

CubitStatus CubitObserver::unregister_observable(CubitObservable *observable)
{
  if (observable == NULL)
    return CUBIT_FAILURE;

    //- remove this observer from the observable's list, and if successful,
    //- decrement the observable count
  CubitStatus success = CUBIT_SUCCESS;

    //- only call the remove function on the observable if we're not being called
    //- from the observable
  success = observable->remove_observer(this);

  if (success == CUBIT_SUCCESS)
    observableCount--;
  
  return success;
}

CubitStatus CubitObserver::register_observer(CubitObserver* obs)
{
    AppUtil::instance()->event_dispatcher().add_observer(obs);
    return CUBIT_SUCCESS;
}

CubitStatus CubitObserver::unregister_observer(CubitObserver* obs)
{
    AppUtil::instance()->event_dispatcher().remove_observer(obs);
    return CUBIT_SUCCESS;
}
