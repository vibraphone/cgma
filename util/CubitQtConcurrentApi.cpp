//! \file CubitQtConcurrentApi.cpp

#include "CubitQtConcurrentApi.h"
#include <QtConcurrentRun>
#include <QtConcurrentMap>
#include <QFutureSynchronizer>
#include <QFutureWatcher>
#include <QEventLoop>
#include <QThreadStorage>

#include <QCoreApplication>

namespace {

  class QtMutex : public QMutex, public CubitConcurrent::Mutex
  {
  public:
    //QtMutex() : QMutex(QMutex::Recursive) {}
    QtMutex() : QMutex() {}
    ~QtMutex() {}
    void lock()
    {
      QMutex::lock();
    }
    void unlock()
    {
      QMutex::unlock();
    }
  };


  struct TLSWrapper
  {
    TLSWrapper(void* p, void (*cleanup)(void*)) : mP(p) {}
    ~TLSWrapper()
    {
      (*mCleanup)(mP);
    }

    void* mP;
    void (*mCleanup)(void*);
  };


  struct QtTLS : public QThreadStorage<TLSWrapper*>,  public CubitConcurrent::ThreadLocalStorageInterface
  {
    QtTLS(void (*cleanup)(void*)) : mCleanup(cleanup)
    {
    }

    void* local_data()
    {
      return localData()->mP;
    }

    void set_local_data(void* p)
    {
      TLSWrapper* w = new TLSWrapper(p, mCleanup);
      setLocalData(w);
    }

    void (*mCleanup)(void*);

  };
}
CubitQtConcurrent::CubitQtConcurrent()
{
    _name = "CubitQtConcurrent";

    // If there is no global instance, set this object as the instance.
  if(!CubitConcurrent::mInstance)
    CubitConcurrent::mInstance = this;
}

CubitQtConcurrent::~CubitQtConcurrent()
{
  // If this is the global instance, clear the pointer.
  if(this == CubitConcurrent::mInstance)
    CubitConcurrent::mInstance = 0;


}

const std::string& CubitQtConcurrent::get_name() const
{
    return _name;
}

const char* CubitQtConcurrent::get_type() const
{
    return _name.c_str();
}

/*
CubitConcurrent::Mutex* CubitQtConcurrent::create_mutex()
{
  return new QtMutex;
}

void CubitQtConcurrent::destroy_mutex(CubitConcurrent::Mutex* m)
{
  delete static_cast<QtMutex*>(m);
}
*/

CubitConcurrent::ThreadLocalStorageInterface* CubitQtConcurrent::create_local_storage(void (*cleanup_function)(void*))
{
  return new QtTLS(cleanup_function);
}

void CubitQtConcurrent::destroy_local_storage(ThreadLocalStorageInterface* i)
{
  delete static_cast<QtTLS*>(i);
}

void CubitQtConcurrent::schedule(CubitConcurrent::Task* task)
{
    QFuture<void> f = ::QtConcurrent::run(task, &Task::execute);    
    m.lock();
    taskmap[task] = f;
    m.unlock();
}

void CubitQtConcurrent::wait(CubitConcurrent::Task* task)
{
    m.lock();
    QMap<CubitConcurrent::Task*, QFuture<void> >::iterator iter = taskmap.find(task);
    m.unlock();
    iter->waitForFinished();
    m.lock();
    taskmap.erase(iter);
    m.unlock();
}

void CubitQtConcurrent::idle_wait(CubitConcurrent::Task* task)
{
  m.lock();
  QMap<CubitConcurrent::Task*, QFuture<void> >::iterator iter = taskmap.find(task);
  m.unlock();
  if(!iter->isFinished())
  {
    QEventLoop loop;
    QFutureWatcher<void> watcher;
    watcher.setFuture(*iter);
    QObject::connect(&watcher, SIGNAL(finished()), &loop, SLOT(quit()));
    loop.exec();
  }
  m.lock();
  taskmap.erase(iter);
  m.unlock();
}
void CubitQtConcurrent::wait_for_any(const std::vector<CubitConcurrent::Task*>& tasks,std::vector<CubitConcurrent::Task*>& finished_tasks)
{
    m.lock();
    for(size_t i=0; i<tasks.size(); i++)
    {
      QMap<CubitConcurrent::Task*, QFuture<void> >::iterator iter = taskmap.find(tasks[i]);
      if(iter->isFinished())
      {
        finished_tasks.push_back(tasks[i]);
        taskmap.erase(iter);
      }
    }
    m.unlock();

    if(!finished_tasks.empty())
      return;


    if(!QCoreApplication::instance())
    {
      int arg=0;
      new QCoreApplication(arg,NULL);
    }

    QEventLoop evLoop;

    m.lock();
    for(size_t i=0; i<tasks.size(); i++)
    {
      QFutureWatcher<void> *f= new QFutureWatcher<void>(&evLoop);
      QMap<CubitConcurrent::Task*, QFuture<void> >::iterator iter = taskmap.find(tasks[i]);
      f->setFuture(*iter);
      QObject::connect(f,SIGNAL(finished()),&evLoop,SLOT(quit()));
    }
    m.unlock();

    evLoop.exec();

    m.lock();
    for(size_t i=0; i<tasks.size(); i++)
    {
      QMap<CubitConcurrent::Task*, QFuture<void> >::iterator iter = taskmap.find(tasks[i]);
      if(iter->isFinished())
      {
        finished_tasks.push_back(tasks[i]);
        taskmap.erase(iter);
      }
    }
    m.unlock();

}
void CubitQtConcurrent::wait(const std::vector<CubitConcurrent::Task*>& task)
{
    m.lock();
    std::vector<QMap<CubitConcurrent::Task*, QFuture<void> >::iterator> iters;
    QFutureSynchronizer<void> f;
    for(size_t i=0; i<task.size(); i++)
    {
        QMap<CubitConcurrent::Task*, QFuture<void> >::iterator iter = taskmap.find(task[i]);
        iters.push_back(iter);
        f.addFuture(*iter);
    }
    m.unlock();

    f.waitForFinished();

    m.lock();
    for(size_t i=0; i<iters.size(); i++)
    {
        taskmap.erase(iters[i]);
    }
    m.unlock();
}

bool CubitQtConcurrent::is_completed(CubitConcurrent::Task* task)
{
    m.lock();
    QMap<CubitConcurrent::Task*, QFuture<void> >::iterator iter = taskmap.find(task);
    m.unlock();

    return iter->isFinished();
}

bool CubitQtConcurrent::is_running(CubitConcurrent::Task* task)
{
    m.lock();
    QMap<CubitConcurrent::Task*, QFuture<void> >::iterator iter = taskmap.find(task);
    m.unlock();

    return iter->isRunning();
}

// internal helper to start tasks within a task group
void CubitQtConcurrent::execute(CubitConcurrent::Task* t)
{
    t->execute();
}

void CubitQtConcurrent::schedule(CubitConcurrent::TaskGroup* tg)
{
    QFuture<void> f = ::QtConcurrent::map(tg->tasks, CubitQtConcurrent::execute);
    m2.lock();
    taskgroupmap[tg] = f;
    m2.unlock();
}

void CubitQtConcurrent::wait(CubitConcurrent::TaskGroup* tg)
{
    m2.lock();
    QMap<CubitConcurrent::TaskGroup*, QFuture<void> >::iterator iter = taskgroupmap.find(tg);
    m2.unlock();
    iter->waitForFinished();
    m2.lock();
    taskgroupmap.erase(iter);
    m2.unlock();
}

bool CubitQtConcurrent::is_completed(CubitConcurrent::TaskGroup* tg)
{
    m2.lock();
    QMap<CubitConcurrent::TaskGroup*, QFuture<void> >::iterator iter = taskgroupmap.find(tg);
    m2.unlock();

    return iter->isFinished();
}

bool CubitQtConcurrent::is_running(CubitConcurrent::TaskGroup* tg)
{
    m2.lock();
    QMap<CubitConcurrent::TaskGroup*, QFuture<void> >::iterator iter = taskgroupmap.find(tg);
    m2.unlock();

    return iter->isRunning();
}

void CubitQtConcurrent::cancel(CubitConcurrent::TaskGroup* tg)
{
    m2.lock();
    QMap<CubitConcurrent::TaskGroup*, QFuture<void> >::iterator iter = taskgroupmap.find(tg);
    m2.unlock();

    iter->cancel();
}

