//-       Class: MemoryManager
//-       Owner: Jim Hipp
//- Description: MemoryManager provides object information and stack
//-              storage pointer management used in association with
//-              MemoryBlock and MemoryAllocation classes.
//- Checked By: 
//-    Version:
//  Modified 5/7/96 R W Ostensen:  Fixed operator_new function as per recommendation
//  from J Hipp.  Eliminated assignment of headOfFreeList to wild pointer loaction when
//  memAllocatnSize=1.


#include "MemoryBlock.hpp"
#include "MemoryManager.hpp"
#include "CubitMessage.hpp"
#include "ArrayBasedContainer.hpp"
#include "AppUtil.hpp"
#include <cassert>
#include <cstring>
#ifndef WIN32
#include <unistd.h>
#endif


MemoryManager* MemoryManager::memoryManagerListHead = NULL;

// constructors
MemoryManager::MemoryManager() : useNew(false)
{
  // do not allow default constructor (must assign objectSize)
  assert(0);
}

MemoryManager::MemoryManager(const MemoryManager&) : useNew(false)
{
  // do not allow copy
  assert(0);
}

MemoryManager::MemoryManager(const char* name, size_t size, int mem_size,
                             int static_flag)
  : useNew(true)
{
    // allocate space for name and copy into buffer
  if (name) 
  {
    objectName = new char[strlen(name) + 1];
    strcpy(objectName, name);
  }
  else
    objectName = NULL;
  
    // initialize manager parameters
  staticManager   = static_flag;
  memBlockStack   = NULL;
  headOfFreeList  = NULL;
  
  assert(size >= sizeof(char*));
  objectSize = size;
  
  set_memory_allocation_increment(mem_size);
  
    // attach new manager to static stack
  next = memoryManagerListHead;
  memoryManagerListHead = this;
}

// destructor
MemoryManager::~MemoryManager()
{
  // delete objectName space and destroy block memory
  if (objectName)
    delete [] objectName;
  destroy_memory(staticManager);
  
  // find 'prev' manager pointer
  MemoryManager* prev        = NULL;
  MemoryManager* mem_manager = memoryManagerListHead;
  while (mem_manager != this)
  {
    prev        = mem_manager;
    mem_manager = mem_manager->next;
  }

  // remove memory manager from static list

  if (!next)
  {
    if (!prev)
    {
      // no memory managers left

      memoryManagerListHead = NULL;
    }
    else
    {
      // remove tail node

      prev->next = NULL;
    }
  }
  else
  {
    if (!prev)
    {
      // remove head node

      memoryManagerListHead = next;
    }
    else
    {
      // remove intermediate node

      prev->next = next;
    }
  }
}

// set memory block allocation increment
void MemoryManager::set_memory_allocation_increment(int mem_size)
{
  // set memory allocation increment

  if (mem_size <= 0)
  {
     memAllocatnSize = DEFAULT_MEMORY_ALLOC_SIZE;
  }
  else
  {
    memAllocatnSize = mem_size;
  }
}

// destroy allocated block memory
void MemoryManager::destroy_memory(int static_flag)
{
  // assert if this is not a static memory manager and some of it's objects
  // are still in use
  if (!static_flag) {
    assert (!get_used_objects());
  }

  // else delete block memory if any was allocated
  if (memBlockStack) delete memBlockStack;

  // reset stack pointers to empty
  headOfFreeList = NULL;
  memBlockStack  = NULL;
}

// return number of objects allocated
int MemoryManager::get_allocated_objects()
{
  // return total number of objects allocated for this memory manager

  if (memBlockStack)
  {
    return (memBlockStack->get_memory_allocation() / objectSize);
  }
  else
  {
    return 0;
  }
}

// return free objects
int MemoryManager::get_free_objects()
{
  // return total number of free objects (allocated but not used) for this
  // memory manager

  if (headOfFreeList)
  {
    int i = 0;
    char* list_ptr = headOfFreeList;
    while (list_ptr)
    {
      i++;
      list_ptr = *((char**) list_ptr);
    }

    return i;
  }
  else
  {
    return 0;
  }
}

// return used objects
int MemoryManager::get_used_objects()
{
  // return total number of used objects for this memory manager

  return (get_allocated_objects() - get_free_objects());
}

// print allocation information to the command line
void MemoryManager::show_object_memory(const char* name)
{
  int instance = 0;

  // find all instances of memory managers for object: name

  MemoryManager* mem_manager = memoryManagerListHead;
  while (mem_manager)
  {
    if (!strcmp(mem_manager->objectName, name))
    {
      // if found then print pertinent memory allocation information

      instance++;
      if (instance > 1)
      {
        PRINT_INFO("\nObject Name(%d): %s\n\n", instance, name);
      }
      else
      {
        PRINT_INFO("\nObject Name: %s\n\n", name);
      }

      int a_obj = mem_manager->get_allocated_objects();
      int f_obj = mem_manager->get_free_objects();
      int u_obj = a_obj - f_obj;

      PRINT_INFO("  Object Size: %zu     Allocation Increment: %d\n\n",
		 mem_manager->objectSize, mem_manager->memAllocatnSize);
      PRINT_INFO("  Allocated Objects: %d  (bytes) %d\n", a_obj,
		 a_obj * (int)(mem_manager->objectSize));
      if (a_obj)
      {
        PRINT_INFO("       Free Objects: %d  (bytes) %d (%d%%)\n", f_obj,
	  	   f_obj * (int)(mem_manager->objectSize),
		   (100*f_obj)/a_obj);
        PRINT_INFO("       Used Objects: %d  (bytes) %d (%d%%)\n", u_obj,
		   u_obj * (int)(mem_manager->objectSize),
		   (100*u_obj)/a_obj);
      }
    }

    // get next memory manager

    mem_manager = mem_manager->next;
  }

  // if none were found then announce

  if (!instance)
  {
    PRINT_INFO("\nObject: %s was not found ...\n",  name);
  }
}


// show allocation for all memory manager objects
void MemoryManager::show_all_object_memory()
{

#if defined(MACOSX)
	
  pid_t PID = getpid();
  char command1[60];
  unsigned long rss, vm;

  FILE *pipe;
  char buf[1024];
  
  //get size of real memory
  sprintf(command1,"ps -o rss -p %d | grep -v RSS",PID);
  pipe = popen(command1, "r");
  if (pipe) 
  {
    fgets(buf, 1024, pipe);
    rss = strtoul(buf, NULL, 0);
    pclose(pipe);
  }

  //get size of virtual memory
  sprintf(command1,"ps -o vsz -p %d | grep -v VSZ",PID);
  pipe = popen(command1, "r");
  if (pipe) 
  {
    fgets(buf, 1024, pipe);
    vm = strtoul(buf, NULL, 0);
    pclose(pipe);
  }

  PRINT_INFO("Total memory = %u\n", vm );
  PRINT_INFO("Resident memory = %u\n", rss );


/*
  struct rusage my_rusage;
  int ret_val = getrusage( RUSAGE_CHILDREN, &my_rusage ); 

  if( ret_val == 0 )
  {
    PRINT_INFO("It was a success\n");
    PRINT_INFO("Memory size = %d\n", my_rusage.ru_maxrss );
    PRINT_INFO("Unshared data size = %d\n", my_rusage.ru_idrss);
    PRINT_INFO("Integeral unshared data size = %d\n", my_rusage.ru_isrss);
    PRINT_INFO("more values: %d %d %d %d %d %d %d %d %d %d %d \n",
		my_rusage.ru_ixrss, my_rusage.ru_minflt, my_rusage.ru_majflt, my_rusage.ru_nswap, 
                my_rusage.ru_inblock, my_rusage.ru_oublock, my_rusage.ru_msgsnd, my_rusage.ru_msgrcv, 
                my_rusage.ru_nsignals, my_rusage.ru_nvcsw, my_rusage.ru_nivcsw ); 
  }
  else
    PRINT_INFO("It was a failure\n");
 */  

 
/* 
  int i, mib[4];
  size_t len;
  struct kinfo_proc kp;

  len = 4;
  sysctlnametomib("kern.proc.pid", mib, &len);
  len = sizeof(kp);
  int pid = getpid();
  mib[3] = pid;
  if (sysctl(mib, 4, &kp, &len, NULL, 0) == -1)
  {
    perror("sysctl");
    PRINT_INFO("Got problems\n");
  }
  else if (len > 0)
  {
    PRINT_INFO("The call was successful!!!\n");
  }
*/

/*
  int i, mib[4];
  size_t len;
  struct kinfo_proc kp;
 
  len = 4; 
  sysctlnametomib("kern.proc.pid", mib, &len);

  for (i = 0; i < 100; i++) 
  {         
    mib[3] = i;         
    len = sizeof(kp);         
    if (sysctl(mib, 4, &kp, &len, NULL, 0) == -1)
      perror("sysctl");         
    else if (len > 0)   
      PRINT_INFO("Call was successful!\n"); 
  }
*/ 

#endif


#if defined(CUBIT_LINUX)
  unsigned long vm, rss;
  process_mem_usage( vm, rss );

  unsigned long read, write;
  process_file_io( read, write );

  PRINT_INFO("Total memory = %lu\n", vm );
  PRINT_INFO("Resident memory = %lu\n", rss );
  PRINT_INFO("Bytes read = %lu\n", read );
  PRINT_INFO("Bytes written = %lu\n", write );
#endif

  /*
  long int a_obj_byte_total = 0;
  long int f_obj_byte_total = 0;
  long int u_obj_byte_total = 0;

  // loop over all memory managers

  PRINT_INFO("\nDynamic Memory Allocation per Object\n\n");

  MemoryManager* mem_manager = memoryManagerListHead;
  while (mem_manager) {
    long int a_obj = mem_manager->get_allocated_objects();
    long int f_obj = mem_manager->get_free_objects();
    long int u_obj = a_obj - f_obj;
    if (a_obj) {
      // sum total allocated memory parameters (bytes)
      a_obj_byte_total += a_obj * mem_manager->objectSize;
      f_obj_byte_total += f_obj * mem_manager->objectSize;
      u_obj_byte_total += u_obj * mem_manager->objectSize;
    }
    mem_manager = mem_manager->next;
  }

  mem_manager = memoryManagerListHead;
  while (mem_manager)
  {
    // print pertinent memory allocation information 

    long int a_obj = mem_manager->get_allocated_objects();
    long int f_obj = mem_manager->get_free_objects();
    long int u_obj = a_obj - f_obj;

    if (a_obj)
    {
      PRINT_INFO("\nObject Name: %s\n\n", mem_manager->objectName);

      PRINT_INFO("  Object Size: %d     Allocation Increment: %d\n\n",
	         mem_manager->objectSize, mem_manager->memAllocatnSize);
      if (a_obj_byte_total != 0)
      {
         PRINT_INFO("  Allocated Objects: %ld  (bytes) %ld (%d%% of Total)\n",
                    a_obj, a_obj * mem_manager->objectSize,
                    int((100.0*a_obj * mem_manager->objectSize)/a_obj_byte_total));
         PRINT_INFO("       Free Objects: %ld  (bytes) %ld (%d%%)\n", f_obj,
                    f_obj * mem_manager->objectSize,
                    int((100.0*f_obj)/a_obj));
         PRINT_INFO("       Used Objects: %ld  (bytes) %ld (%d%%)\n", u_obj,
                    u_obj * mem_manager->objectSize,
                    int((100.0*u_obj)/a_obj));
      }
      else
      {
         PRINT_INFO("  Allocated Objects: %ld  (bytes) %ld (100%% of Total)\n",
                    a_obj, a_obj * mem_manager->objectSize);
         PRINT_INFO("       Free Objects: %ld  (bytes) %ld (100%%)\n", f_obj,
                    f_obj * mem_manager->objectSize);
         PRINT_INFO("       Used Objects: %ld  (bytes) %ld (100%%)\n", u_obj,
                    u_obj * mem_manager->objectSize);
      }
    }

    mem_manager = mem_manager->next;
  }

  // print total memory allocation information

  char sizechar;

  PRINT_INFO("\nTotal Memory Allocation Information\n\n");
  int divisor;
  if (a_obj_byte_total > 10000000) {
    sizechar = 'M';
    divisor = 1000000;
  }
  else {
    sizechar = 'K';
    divisor = 1000;
  }
  
  PRINT_INFO("  Allocated Memory: %ld%c (%ld bytes)\n", a_obj_byte_total/divisor, 
             sizechar, a_obj_byte_total);

  if (a_obj_byte_total)
  {
    PRINT_INFO("       Free Memory: %ld%c (%d%%)\n", f_obj_byte_total/divisor, sizechar,
	       int((100.0*f_obj_byte_total)/a_obj_byte_total));
    PRINT_INFO("       Used Memory: %ld%c (%d%%)\n", u_obj_byte_total/divisor, 
               sizechar,
	       int((100.0*u_obj_byte_total)/a_obj_byte_total));
  }

#ifndef JANUS
#ifndef WIN32
  struct rusage r_usage;
  AppUtil::instance()->apputil_getrusage(r_usage);
  PRINT_INFO("       (System reports %ld%c used, incl. executable)\n", 
              r_usage.ru_maxrss*getpagesize()/divisor, sizechar);
#else
  PRINT_INFO("\n");
#endif // WIN32
#endif // JANUS

  // print DLList non-Pool allocation information
  PRINT_INFO("\nTotal non-pool ArrayBasedContainer memory allocation  = %u%c\n"
	       "Maximum non-pool ArrayBasedContainer memory allocated = %u%c\n",
             ArrayBasedContainer::current_allocated_memory()/divisor, sizechar, 
             ArrayBasedContainer::maximum_allocated_memory()/divisor, sizechar);
*/



#if 0
  // print HOOPS memory usage
  long allocated = 0;
  long in_use = 0;
    //  DrawingTool::instance()->show_memory(allocated, in_use);
  if (allocated != 0)
     PRINT_INFO("\nGraphics subsystem memory: Allocated = %u%c (%d bytes)\n"
                "                           In-Use    = %u%c (%d%%)\n",
                allocated/divisor, sizechar, allocated, 
                in_use/divisor, sizechar, int((100.0*in_use)/allocated));
  else
     PRINT_INFO("\nGraphics subsystem memory: Allocated = %u%c (%d bytes)\n"
                "                           In-Use    = %u%c (100%%)\n",
                allocated/divisor, sizechar, allocated, 
                in_use/divisor, sizechar);
#endif
}

// compress memory for the requested object
int MemoryManager::compress_object_memory(const char* name)
{
  // find all instances of memory managers for object: name

  int found        = 0;
  int saved_memory = 0;
  MemoryManager* mem_manager = memoryManagerListHead;
  while (mem_manager)
  {
    if (!strcmp(mem_manager->objectName, name))
    {
      // if found then compress memory

      saved_memory += mem_manager->compress_memory();
      found         = 1;
    }

    mem_manager = mem_manager->next;
  }

  return found;
}

// compress all object memory
int MemoryManager::compress_all_object_memory()
{
  // find all instances of memory managers
  int saved_memory = 0;
  MemoryManager* mem_manager = memoryManagerListHead;
  MemoryManager* block_manager = NULL;
  while (mem_manager)
  {
    if (!strcmp(mem_manager->objectName, "MemoryBlock"))
    {
      // save block_manager until end

      block_manager = mem_manager;
    }
    else
    {
      // compress memory

      saved_memory += mem_manager->compress_memory();
    }

    mem_manager = mem_manager->next;
  }

  if (block_manager)
  {
    saved_memory += block_manager->compress_memory();
  }

  return saved_memory;
}
      
// generic operator new call
void* MemoryManager::operator_new(size_t size)
{

  if (useNew) return malloc(size);

  // send requests of "wrong" size to ::new
  
  try
  { 
    if (size != objectSize) return ::new char[size];
  }
  catch(...) 
  {
    return (void*) NULL;
  }
  // get new element from head of free list

  char* p = headOfFreeList;

  try
  {
    if(!p)
        {
        // allocate new block

        int block_size  = memAllocatnSize * size;
        char* new_block = ::new char[block_size];
        if (!new_block) return (void*) NULL;

        // link new elements to form the free list

        int fill_limit = (memAllocatnSize - 1) * size;
        for (int j = 0; j < fill_limit; j += size)
        {
          *((char**) &new_block[j]) = &new_block[j + size];
        }
        *((char**) &new_block[fill_limit]) = (char*) NULL;

        // assign new element

        p = new_block;

        // save new block to memory block stack

        memBlockStack = new MemoryBlock(memBlockStack, new_block, block_size);
      }
  }
  catch(...)
  {
    return (void*) NULL;
  }
  //assign head of free list and return p

  headOfFreeList = *((char**) p);
  return (void*) p;
}

// generic operator delete call
void MemoryManager::operator_delete(void *deadObject, size_t size)
{
  if (useNew) 
  {
    free(deadObject);
    return;
  }

  // requests of "wrong" size to ::delete

  if (size != objectSize)
  {
    ::delete [] ((char*) deadObject);
    return;
  }

  // attach dead element to head of free list

  char* delete_object = (char*) deadObject;
  *((char**) delete_object) = headOfFreeList;
  headOfFreeList = delete_object;
}

// compress memory blocks
int MemoryManager::compress_memory()
{
  // if free objects exist then begin compression algorithm

  if (headOfFreeList)
  {
    // find total number of memory blocks attached to stack
 
    int n_blocks = 0;
    MemoryBlock* mem_block = memBlockStack;
    while (mem_block)
    {
      n_blocks++;
      mem_block = mem_block->next_block();
    }

    if (n_blocks == 0)
    {
       // no available memory to free ... return 0
       // this is here for safety ... n_blocks should never be zero if
       // headOfFreeList is not Null

       return 0;
    }
    else
    {
      // first determine if all objects are free ... if so then perform
      // the easy compression routine

      if (!get_used_objects())
      {
        // all objects are free ... delete all blocks

        int n_bytes = memBlockStack->get_memory_allocation();
        destroy_memory(staticManager);

        // return freed memory

        return n_bytes;
      }

      // else perform the complex routine to remove those memory blocks that
      // have all free elements

      // begin by constructing an integer array for each memory block to
      // tally the number of free objects that each block contains

        // if there are a lot of blocks, we can save a huge amount of
        // time by looking in the last few blocks that contained an
        // element.
      const int use_cache = n_blocks > 8;
      int i, j, k;
      i = j = k = 0;
      const int cache_size = 4;
      MemoryBlock* mem_block_sav[cache_size];
      int i_sav[cache_size];
      for ( i = cache_size; i--; )
      {
        mem_block_sav[i] = NULL;
        i_sav[i] = 0;
      }
      int found = 0;
      
      mem_block = NULL;
      char* list_ptr = NULL;
      
      unsigned int* free_tally = new unsigned int [n_blocks];
      for (i = 0; i < n_blocks; i++) free_tally[i] = 0;

      // loop through free list tallying free elements

      list_ptr = headOfFreeList;
      while (list_ptr)
      {
        // find memory block that owns this element
        
          // check last few blocks for speed
        found = CUBIT_FALSE;
        if ( use_cache )
        {
          for ( i = 0; i < cache_size; i++ )
          {
            mem_block = mem_block_sav[i];
            if ( mem_block &&
                 list_ptr >= mem_block->get_block() &&
                 list_ptr < (mem_block->get_block() + 
                             mem_block->block_size()) )
            {
              k = i_sav[i];
              free_tally[k]++;
              found = CUBIT_TRUE;
              break;
            }
          }
        }
        if ( !found )
        {
            // search through all blocks
          mem_block = memBlockStack;
          for (i = 0; i < n_blocks; i++)
          {
            if ((list_ptr >= mem_block->get_block()) &&
                (list_ptr < (mem_block->get_block() + mem_block->block_size())))
            {
                // increment tally and exit
              
              free_tally[i]++;
                // save
              if ( use_cache && mem_block_sav[j] != mem_block )
              {
                mem_block_sav[j] = mem_block;
                i_sav[j] = i;
                if ( ++j >= cache_size )
                  j = 0;
              }
              break;
            }
            
              //  get next memory block
            mem_block = mem_block->next_block();
          }
        }

        // get next element
        list_ptr = *((char**) list_ptr);
      }

      // zero tally for memory blocks that cannot be removed ... those that
      // have some used elements

      int all_blocks = 0;
      mem_block = memBlockStack;
      for (i = 0; i < n_blocks; i++)
      {
        if (free_tally[i] != (mem_block->block_size() / objectSize))
	{
          free_tally[i] = 0;
          all_blocks++;
        }

        mem_block = mem_block->next_block();
      }

      if (all_blocks == n_blocks)
      {
        // no memory can be saved ... all blocks have some used elements
        // return 0

        delete [] free_tally;
        return 0;
      }

      // adjust free list pointers to remove those that belong to
      // memory blocks that can be deleted
      char* prev_ptr = NULL;
      list_ptr = headOfFreeList;
      while (list_ptr)
      {
        // find memory block that owns this element
          // check last few blocks for speed
        found = CUBIT_FALSE;
        if ( use_cache )
        {
          for ( i = 0; i < cache_size; i++ )
          {
            mem_block = mem_block_sav[i];
            if ( mem_block &&
                 list_ptr >= mem_block->get_block() &&
                 list_ptr < (mem_block->get_block() + 
                             mem_block->block_size()) )
            {
              k = i_sav[i];
              found = CUBIT_TRUE;
              break;
            }
          }
        }
        if ( !found )
        {
          mem_block = memBlockStack;
          for (i = 0; i < n_blocks; i++)
          {
            if ((list_ptr >= mem_block->get_block()) &&
                (list_ptr < (mem_block->get_block() + mem_block->block_size())))
            {
              k = i;
                // save
              if ( use_cache && mem_block_sav[j] != mem_block )
              {
                mem_block_sav[j] = mem_block;
                i_sav[j] = i;
                if ( ++j >= cache_size )
                  j = 0;
              }
              break;
            }
              // get next memory block
            mem_block = mem_block->next_block();
          }
        }
        
        if (free_tally[k])
        {
            // remove element
          
          if (prev_ptr)
          {
            *((char**) prev_ptr) = *((char**) list_ptr);
          }
          else
          {
            headOfFreeList = *((char**) list_ptr);
          }
        }
        else
        {
            // advance prev_ptr  
          prev_ptr = list_ptr;
        }

        // get next element
        list_ptr = *((char**) list_ptr);
      }

      // delete all memory blocks that have free_tally[i] > 0

      i = 0;
      int save_bytes = 0;
      MemoryBlock* prev_block = NULL;
      mem_block               = memBlockStack;
      while (mem_block)
      {
        if (free_tally[i])
        {
          // set previous MemoryBlocks next pointer to skip this block

          if (prev_block)
          {
            prev_block->next_block(mem_block->next_block());
          }
          else
          {
            memBlockStack = mem_block->next_block();
          }

          // set MemoryBlock next pointer to NULL to avoid recusive delete
          // update saved memory and delete mem_block

          mem_block->next_block((MemoryBlock*) NULL);
          save_bytes += mem_block->block_size();
          delete mem_block;

          // update mem_block to point to new current MemoryBlock

          if (prev_block)
          {
            mem_block = prev_block->next_block();
          }
          else
          {
            mem_block = memBlockStack;
          }
        }
        else
        {
          // if block wasn't removed then update previous and current blocks

          prev_block = mem_block;
          mem_block = mem_block->next_block();
        }

        // increment to next block (used by free_tally array)

        ++i;
      }

      // return freed memory (bytes)

      delete [] free_tally;
      return save_bytes;
    }
  }
  else
  {
    // no memory allocated ... return 0

    return 0;
  }
}


void MemoryManager::process_mem_usage(unsigned long &vm_usage, unsigned long &resident_set)
{

#if defined(CUBIT_LINUX)
  using std::ios_base;
  using std::ifstream;
  using std::string;

  vm_usage     = 0;
  resident_set = 0;

  // 'file' stat seems to give the most reliable results
  //
  ifstream stat_stream("/proc/self/stat",ios_base::in);
  
  // dummy vars for leading entries in stat that we don't care about
  //
  string pid, comm, state, ppid, pgrp, session, tty_nr;
  string tpgid, flags, minflt, cminflt, majflt, cmajflt;
  string utime, stime, cutime, cstime, priority, nice;
  string O, itrealvalue, starttime;
  
  // the two fields we want...virtual memory size and resident memory size
  //
  unsigned long vsize;
  long rss;
  
  stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
              >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
              >> utime >> stime >> cutime >> cstime >> priority >> nice
              >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest
              
  long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
  vm_usage     = vsize / 1024;
  resident_set = rss * page_size_kb;
#endif
}
  
 
void MemoryManager::process_file_io(unsigned long &read, unsigned long &write )
{
#if defined(CUBIT_LINUX)
  using std::ios_base;
  using std::ifstream;
  using std::string;

  read = 0;
  write = 0;

  // 'file' stat seems to give the most reliable results
  //
  ifstream stat_stream("/proc/self/io",ios_base::in);
  
  // dummy vars for leading entries in stat that we don't care about
  //
  string char1, char2;

  //----------------------Getting two numbers out of this file 
  // I/O counter: chars read
  //The number of bytes which this task has caused to be read from storage. This
  //is simply the sum of bytes which this process passed to read() and pread().
  //It includes things like tty IO and it is unaffected by whether or not actual
  //physical disk IO was required (the read might have been satisfied from
  //pagecache)
  
  // I/O counter: chars written
  //The number of bytes which this task has caused, or shall cause to be written
  //to disk. Similar caveats apply here as with rchar.
  
  unsigned long tmp_read, tmp_write;
  
  stat_stream >> char1 >> tmp_read >> char2 >> tmp_write; //don't care about the rest 
  
  read = tmp_read;
  write = tmp_write;
#endif
}



