
#include "CubitEventDispatcher.hpp"
#include "CubitObserver.hpp"
#include "CubitObservable.hpp"
#include <algorithm>

CubitEventDispatcher::CubitEventDispatcher()
{
}

CubitEventDispatcher::~CubitEventDispatcher()
{
}

void CubitEventDispatcher::send_event(CubitObservable* observable, const CubitEvent& event)
{
  if(observable)
  {
      observable->notify_observers(event);
  }

  std::vector<CubitObserver*>::iterator iter;
  for(iter = mObservers.begin(); iter != mObservers.end(); ++iter)
  {
    (*iter)->notify_observer(observable, event);
  }
}

void CubitEventDispatcher::add_observer(CubitObserver* obs)
{
  std::vector<CubitObserver*>::iterator iter;
  iter = std::find(mObservers.begin(), mObservers.end(), obs);
  if(iter == mObservers.end())
  {
    mObservers.push_back(obs);
  }
}

void CubitEventDispatcher::remove_observer(CubitObserver* obs)
{
  mObservers.erase(
      std::remove(mObservers.begin(), mObservers.end(), obs),
      mObservers.end()
      );
}

