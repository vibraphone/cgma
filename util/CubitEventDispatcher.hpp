
#ifndef CubitEventDispatcher_hpp
#define CubitEventDispatcher_hpp

class CubitObserver;
class CubitEvent;
class CubitObservable;
#include <vector>
#include "CubitUtilConfigure.h"

//! An event dispatcher class to handle sending events to observers.
//! A CubitObserver may register with this class to receive events.
class CUBIT_UTIL_EXPORT CubitEventDispatcher
{
  public:
    CubitEventDispatcher();
    virtual ~CubitEventDispatcher();

    //! Send an event to observers.
    //! The event is sent two ways:
    //!  1. If an observable is given, the event is sent to the observers directly observing that entity.
    //!  2. The event is sent to all observers added to this dispatcher.
    //! Ownership of the event is taken by this object.
    void send_event(CubitObservable* observable, const CubitEvent& event);

    //! Add an observer.
    void add_observer(CubitObserver* obs);
    //! Remove an observer.
    void remove_observer(CubitObserver* obs);

  protected:

    std::vector<CubitObserver*> mObservers;

};

#endif
