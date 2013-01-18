//! \file CubitConcurrentApi.cpp

#include "CubitConcurrentApi.h"

#ifdef WIN32
#include <windows.h>
#else
#include <pthread.h>
#endif

#ifdef WIN32
  // struct for working with mutexes
  class RealMutex : public CubitConcurrent::Mutex
  {
  public:
    RealMutex()
    {
      InitializeCriticalSection(&mCriticalSection);
    }
    ~RealMutex()
    {
      DeleteCriticalSection(&mCriticalSection);
    }
    virtual void lock()
    {
      EnterCriticalSection(&mCriticalSection);
    }
    virtual void unlock()
    {
      LeaveCriticalSection(&mCriticalSection);
    }
  protected:
    CRITICAL_SECTION mCriticalSection;
  };

#else
  // struct for working with mutexes
  class RealMutex : public CubitConcurrent::Mutex
  {
  public:
    RealMutex()
    {
      pthread_mutex_init(&mMutex, NULL);
    }
    ~RealMutex()
    {
      pthread_mutex_destroy(&mMutex);
    }
    virtual void lock()
    {
      pthread_mutex_lock(&mMutex);
    }
    virtual void unlock()
    {
      pthread_mutex_unlock(&mMutex);
    }
  protected:
    pthread_mutex_t mMutex;
  };
#endif

CubitConcurrent *CubitConcurrent::mInstance = 0;


CubitConcurrent::CubitConcurrent()
{
}

CubitConcurrent::~CubitConcurrent()
{
}

const char* CubitConcurrent::get_base_type() const
{ 
    return "ConcurrentApi";
}

CubitConcurrent::Mutex* CubitConcurrent::create_mutex()
{
  return new RealMutex;
}

void CubitConcurrent::destroy_mutex(CubitConcurrent::Mutex* m)
{
  delete m;
}
