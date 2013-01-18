//! \file CubitQtConcurrentApi.h
/*! \brief Api for concurrency based on Qt
 */ 

#ifndef CUBIT_QT_CONCURRENT_API_H_
#define CUBIT_QT_CONCURRENT_API_H_

#include "CubitConcurrentApi.h"
#include "CubitUtilConfigure.h"
#include <QMutex>
#include <QMap>
#include <QFuture>


// class to provide a way to run tasks concurrently
class CUBIT_UTIL_EXPORT CubitQtConcurrent : public CubitConcurrent
{
public:
  CubitQtConcurrent();
  virtual ~CubitQtConcurrent();
  
  const std::string& get_name() const;
  const char* get_type() const;

  //Mutex* create_mutex();
  //void destroy_mutex(Mutex* m);

  ThreadLocalStorageInterface* create_local_storage(void (*cleanup_function)(void*));
  void destroy_local_storage(ThreadLocalStorageInterface* i);

  // wait for a task to finish
  virtual void wait(Task* task);

  // wait for a task to finish, but do not take on the task if the thread pool hasn't taken it yet
  virtual void idle_wait(Task* task);

  // wait for a set of tasks to finish
  virtual void wait(const std::vector<Task*>& task);

  //wait for any of a set of tasks to finish
  void wait_for_any(const std::vector<Task*>& tasks,std::vector<Task*>& finished_tasks);

  // return whether a task is complete
  virtual bool is_completed(Task* task);
  
  // return whether a task is currently running (as opposed to waiting in the queue)
  virtual bool is_running(Task* task);
  
  // wait for a task group to complete
  virtual void wait(TaskGroup* task_group);
  
  // return whether a task group is complete
  virtual bool is_completed(TaskGroup* task_group);
  
  // return whether a task group is currently running (as opposed to waiting in the queue)
  virtual bool is_running(TaskGroup* task_group);
 
  // cancel a task group's execution
  // this only un-queues tasks that haven't started, and running tasks will run to completion
  // after canceling a task group, one still needs to call wait() for completion.
  virtual void cancel(TaskGroup* task_group);
  

protected: 
  // schedule a task for execution
  virtual void schedule(Task* task);
  
  // schedule a group of tasks for execution
  virtual void schedule(TaskGroup* task_group);
  
  static void execute(Task* t);
  
  QMap<TaskGroup*, QFuture<void> > taskgroupmap;
  QMutex m2;

  QMap<Task*, QFuture<void> > taskmap;
  QMutex m;
  
  std::string _name;
};


#endif // CUBIT_QT_CONCURRENT_API_H_

