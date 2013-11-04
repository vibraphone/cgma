//! \file CubitConcurrentApi.h
/*! \brief Api for concurrency
 */ 

#ifndef _CUBIT_CONCURRENT_API_H_
#define _CUBIT_CONCURRENT_API_H_


#include "CubitUtilConfigure.h"
#include <vector>

// class to provide a way to run tasks concurrently
class  CUBIT_UTIL_EXPORT CubitConcurrent 
{
public:
  CubitConcurrent();
  virtual ~CubitConcurrent();

    //!   Gets the global Concurrent instance.
    //! \return
    //!   A Pointer to the global Concurrent instance.
    static CubitConcurrent *instance() {return mInstance;}


  // struct for working with mutexes
  class Mutex
  {
  public:
    virtual ~Mutex() {}
    virtual void lock() = 0;
    virtual void unlock() = 0;
  };

  // convenience class to lock/unlock mutex within a scope
  // use this whenever possible
  // it also provides a safe lock()/unlock()
  struct MutexLocker
  {
    MutexLocker(Mutex* mutex) : mMutex(mutex), mLocked(false)
    {
      lock();
    }
    ~MutexLocker()
    {
      unlock();
    }
    void lock()
    {
      if(!mLocked)
      {
        mMutex->lock();
        mLocked = true;
      }
    }
    void unlock()
    {
      if(mLocked)
      {
        mMutex->unlock();
        mLocked = false;
      }
    }
    Mutex* mMutex;
    bool mLocked;
  };

  virtual Mutex* create_mutex();
  virtual void destroy_mutex(Mutex* m);

  // struct for working with thread local storage
  class ThreadLocalStorageInterface
  {
  public:
    virtual void* local_data() = 0;
    virtual void set_local_data(void*) = 0;
  protected:
    virtual ~ThreadLocalStorageInterface() {}
  };

  virtual ThreadLocalStorageInterface* create_local_storage(void (*cleanup_function)(void*)) = 0;
  virtual void destroy_local_storage(ThreadLocalStorageInterface* s) = 0;

  template <class T>
  class ThreadLocalStorage
  {
  public:
    ThreadLocalStorage()
    {
      mTLS = CubitConcurrent::instance()->create_local_storage(ThreadLocalStorage::cleanup);
    }

    ~ThreadLocalStorage()
    {
      CubitConcurrent::instance()->destroy_local_storage(mTLS);
    }

    T* local_data()
    {
      return reinterpret_cast<T*>(mTLS->local_data());
    }

    void set_local_data(T* t)
    {
      mTLS->set_local_data(t);
    }

  private:

    static void cleanup(void* p)
    {
      delete reinterpret_cast<T*>(p);
    }

    ThreadLocalStorageInterface* mTLS;
  };


  // struct to encapsulate a task to execute
  struct Task
  {
      virtual ~Task() {}
      virtual void execute() = 0;
  };

  // struct to represent a group of tasks (several tasks started with one api call)
  // this sometimes gives a simpler setup of tasks and also provides cancellation ability
  struct TaskGroup
  {
    std::vector<Task*> tasks;
  };


  // create a schedule a task for a member function with no arguments
  // for example:
  /*
    class Foo
    {
      void foo()
      {
        Task* t = c->create_and_schedule(*this, &Foo::thread_func);
        c->wait(t);
        delete t;
      }
      
      void thread_func()
      {
      }
    };
  */
  template <typename X>
  Task* create_and_schedule(X& x, void (X::*fun)())
  {
    Task* t = new ClassFunctionTask<X>(x, fun);
    this->schedule(t);
    return t;
  };

  // create and schedule a task for a member function with one argument
  // for example:
  /*
    class Foo
    {
      void foo()
      {
        Task* t = c->create_and_schedule(*this, &Foo::thread_func, 5);
        c->wait(t);
        delete t;
      }
      
      void thread_func(int v)
      {
      }
    };
  */
  template <typename X, typename Param1, typename Arg1>
  Task* create_and_schedule(X& x, void (X::*fun)(Param1), const Arg1& arg1)
  {
    return create_task1<X,Param1,const Arg1&>(x,fun,arg1);
  };
  
  // same as above but to handle references passed through thread function
  template <typename X, typename Param1, typename Arg1>
  Task* create_and_schedule(X& x, void (X::*fun)(Param1), Arg1& arg1)
  {
    return create_task1<X,Param1,Arg1&>(x,fun,arg1);
  };
  
  // create a schedule a task for a member function with two arguments
  // for example:
  /*
    class Foo
    {
      void foo()
      {
        int result;
        Task* t = c->create_and_schedule(*this, &Foo::square, 5, result);
        c->wait(t);
        delete t;
      }
      
      void square(int v, int& result)
      {
        result = v*v;
      }
    };
  */
  template <typename X, typename Param1, typename Param2, typename Arg1, typename Arg2>
  Task* create_and_schedule(X& x, void (X::*fun)(Param1, Param2), const Arg1& arg1, const Arg2& arg2)
  {
    return create_task2<X,Param1,Param2, const Arg1&,const Arg2&>(x,fun,arg1, arg2);
  };
  
  template <typename X, typename Param1, typename Param2, typename Arg1, typename Arg2>
  Task* create_and_schedule(X& x, void (X::*fun)(Param1, Param2), Arg1& arg1, const Arg2& arg2)
  {
    return create_task2<X,Param1,Param2,Arg1&, const Arg2&>(x,fun,arg1, arg2);
  };
  
  template <typename X, typename Param1, typename Param2, typename Arg1, typename Arg2>
  Task* create_and_schedule(X& x, void (X::*fun)(Param1, Param2), const Arg1& arg1, Arg2& arg2)
  {
    return create_task2<X,Param1,Param2,const Arg1&,Arg2&>(x,fun,arg1, arg2);
  };
  
  template <typename X, typename Param1, typename Param2, typename Arg1, typename Arg2>
  Task* create_and_schedule(X& x, void (X::*fun)(Param1, Param2), Arg1& arg1, Arg2& arg2)
  {
    return create_task2<X,Param1,Param2,Arg1&,Arg2&>(x,fun,arg1, arg2);
  };

  // create a schedule a task group for a member function with one argument from a sequence
  // for example:
  /*
    class Foo
    {
      void foo()
      {
        std::vector<int> input(5, 3);
        TaskGroup* t = c->create_and_schedule_group(*this, &Foo::thread_func, input);
        c->wait(t);
        c->delete_group(t);
      }
      
      void thread_func(int v)
      {
      }
    };
  */
  template <typename X, typename Param, typename Sequence>
  TaskGroup* create_and_schedule_group(X& x, void (X::*fun)(Param), Sequence& seq)
  {
    return create_taskgroup1<X, Param, Sequence, typename Sequence::iterator>(x, fun, seq);
  };
  
  template <typename X, typename Param, typename Sequence>
  TaskGroup* create_and_schedule_group(X& x, void (X::*fun)(Param), const Sequence& seq)
  {
    return create_taskgroup1<X, Param, const Sequence, typename Sequence::const_iterator>(x, fun, seq);
  };

  // create a schedule a task group for a member function with two arguments from two sequences
  // for example:
  /*
    class Foo
    {
      void foo()
      {
        std::vector<int> input(5, 3);
        std::vector<int> output(5);
        TaskGroup* t = c->create_and_schedule_group(*this, &Foo::square, input, output);
        c->wait(t);
        c->delete_group(t);
      }
      
      void square(int v, int& result)
      {
        result = v*v;
      }
    };
  */
  template <typename X, typename Param1, typename Param2, typename Sequence1, typename Sequence2>
  TaskGroup* create_and_schedule_group(X& x, void (X::*fun)(Param1, Param2), Sequence1& seq1, Sequence2& seq2)
  {
    return create_taskgroup2<X, Param1, Param2, Sequence1, Sequence2, typename Sequence1::iterator, typename Sequence2::iterator>(x, fun, seq1, seq2);
  };
  
  template <typename X, typename Param1, typename Param2, typename Sequence1, typename Sequence2>
  TaskGroup* create_and_schedule_group(X& x, void (X::*fun)(Param1, Param2), const Sequence1& seq1, Sequence2& seq2)
  {
    return create_taskgroup2<X, Param1, Param2, const Sequence1, Sequence2, typename Sequence1::const_iterator, typename Sequence2::iterator>(x, fun, seq1, seq2);
  };
  
  template <typename X, typename Param1, typename Param2, typename Sequence1, typename Sequence2>
  TaskGroup* create_and_schedule_group(X& x, void (X::*fun)(Param1, Param2), Sequence1& seq1, const Sequence2& seq2)
  {
    return create_taskgroup2<X, Param1, Param2, Sequence1, const Sequence2, typename Sequence1::iterator, typename Sequence2::const_iterator>(x, fun, seq1, seq2);
  };
  
  template <typename X, typename Param1, typename Param2, typename Sequence1, typename Sequence2>
  TaskGroup* create_and_schedule_group(X& x, void (X::*fun)(Param1, Param2), const Sequence1& seq1, const Sequence2& seq2)
  {
    return create_taskgroup2<X, Param1, Param2, const Sequence1, const Sequence2, typename Sequence1::const_iterator, typename Sequence2::const_iterator>(x, fun, seq1, seq2);
  };

  // delete a task group created by create_and_schedule_group
  void delete_group(TaskGroup* tg)
    {
    for(size_t i=0; i<tg->tasks.size(); i++)
      {
      delete tg->tasks[i];
      }
    delete tg;
    }
 


  // wait for a task to finish
  virtual void wait(Task* task) = 0;

  // wait for a task to finish
  // implementation may decide to wait by running a local event loop, instead of pausing this thread or taking on the unstarted task
  // warning: do not use this with recursive tasks
  virtual void idle_wait(Task* task)
  {
    wait(task);
  }
  
  // wait for a set of tasks to finish
  virtual void wait(const std::vector<Task*>& task) = 0;

  //wait for any of a set of tasks to finish
  virtual void wait_for_any(const std::vector<Task*>& tasks,std::vector<Task*>& finished_tasks) = 0;

  // return whether a task is complete
  virtual bool is_completed(Task* task) = 0;
  
  // return whether a task is currently running (as opposed to waiting in the queue)
  virtual bool is_running(Task* task) = 0;
  
  // wait for a task group to complete
  virtual void wait(TaskGroup* task_group) = 0;
  
  // return whether a task group is complete
  virtual bool is_completed(TaskGroup* task_group) = 0;
  
  // return whether a task group is currently running (as opposed to waiting in the queue)
  virtual bool is_running(TaskGroup* task_group) = 0;
 
  // cancel a task group's execution
  // this only un-queues tasks that haven't started, and running tasks will run to completion
  // after canceling a task group, one still needs to call wait() for completion.
  virtual void cancel(TaskGroup* task_group) = 0;
  

protected: 
    static CubitConcurrent *mInstance; //!< Stores the global instance.

  // schedule a task for execution
  virtual void schedule(Task* task) = 0;
  
  // schedule a group of tasks for execution
  virtual void schedule(TaskGroup* task_group) = 0;
  
  
  
  
  // implementation helpers
  template <class X>
  struct ClassFunctionTask : public Task
  {
    public:
      ClassFunctionTask(X& _ptr, void (X::*_MemFun)()) : ptr(_ptr), MemFun(_MemFun) {}
      void execute()
      {
        (ptr.*MemFun)();
      }
    protected:
      X& ptr;
      void (X::*MemFun)();
        private:
          const ClassFunctionTask& operator=(const ClassFunctionTask&);
              ClassFunctionTask(const ClassFunctionTask&);
  };

  template <typename X, typename Param1>
  struct ClassFunctionTaskArg1 : public Task
  {
    public:
      ClassFunctionTaskArg1(X& _ptr, void (X::*_MemFun)(Param1), Param1 _arg1) : ptr(_ptr), MemFun(_MemFun), arg1(_arg1) {}
      void execute()
      {
        (ptr.*MemFun)(arg1);
      }
    protected:
      X& ptr;
      void (X::*MemFun)(Param1);
      Param1 arg1;
        private:
          const ClassFunctionTaskArg1& operator=(const ClassFunctionTaskArg1&);
              ClassFunctionTaskArg1(const ClassFunctionTaskArg1&);
  };

  template <typename X, typename Param1, typename Param2>
  struct ClassFunctionTaskArg2 : public Task
  {
    public:
      ClassFunctionTaskArg2(X& _ptr, void (X::*_MemFun)(Param1, Param2), Param1 _arg1, Param2 _arg2)
        : ptr(_ptr), MemFun(_MemFun), arg1(_arg1), arg2(_arg2) {}
      void execute()
      {
        (ptr.*MemFun)(arg1, arg2);
      }
    protected:
      X& ptr;
      void (X::*MemFun)(Param1, Param2);
      Param1 arg1;
      Param2 arg2;
        private:
          const ClassFunctionTaskArg2& operator=(const ClassFunctionTaskArg2&);
              ClassFunctionTaskArg2(const ClassFunctionTaskArg2&);
  };

  template <typename X, typename Param1, typename Arg1>
  inline Task* create_task1(X& x, void (X::*fun)(Param1), Arg1 arg1)
  {
    Task* t = new ClassFunctionTaskArg1<X, Param1>(x, fun, arg1);
    this->schedule(t);
    return t;
  }
  template <typename X, typename Param1, typename Param2, typename Arg1, typename Arg2>
  inline Task* create_task2(X& x, void (X::*fun)(Param1, Param2), Arg1 arg1, Arg2 arg2)
  {
    Task* t = new ClassFunctionTaskArg2<X, Param1, Param2>(x, fun, arg1, arg2);
    this->schedule(t);
    return t;
  }

  template <typename X, typename Param, typename Sequence, typename Iterator>
  inline TaskGroup* create_taskgroup1(X& x, void (X::*fun)(Param), Sequence& seq)
  {
    TaskGroup* tg = new TaskGroup;
    Iterator iter;
    for(iter = seq.begin(); iter != seq.end(); ++iter)
      {
      Task* t = new ClassFunctionTaskArg1<X, Param>(x, fun, *iter);
      tg->tasks.push_back(t);
      }
    this->schedule(tg);
    return tg;
  }
  
  template <typename X, typename Param1, typename Param2, typename Sequence1, typename Sequence2, typename Iterator1, typename Iterator2>
  inline TaskGroup* create_taskgroup2(X& x, void (X::*fun)(Param1, Param2), Sequence1& seq1, Sequence2& seq2)
  {
    if(seq1.size() != seq2.size())
      return NULL;

    TaskGroup* tg = new TaskGroup;
    Iterator1 iter1;
    Iterator2 iter2;
    for(iter1 = seq1.begin(), iter2 = seq2.begin(); iter1 != seq1.end(); ++iter1, ++iter2)
      {
      Task* t = new ClassFunctionTaskArg2<X, Param1, Param2>(x, fun, *iter1, *iter2);
      tg->tasks.push_back(t);
      }
    this->schedule(tg);
    return tg;
  }


public:

  const char* get_base_type() const;

};


#endif // CUBITCONCURRENT_API_H_

