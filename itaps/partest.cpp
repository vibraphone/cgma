// Test parallel iGeom
// July 09
// Author : Hong-Jun Kim
// read OCC geometry file and distribute it to remote processors
// with various methods

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <sys/resource.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdlib.h>

#include "CGMmpi.h"
#include "iGeom.h"

#define IGEOM_ASSERT(ierr) if (ierr!=0) printf("igeom assert\n");
#define IGEOM_NULL 0

void util_getrusage(struct rusage &r_usage);

int main(int argc, char* argv[]){
  int rank, size;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  printf("Hello\n");

  iGeom_Instance igeom;
  int ierr;
  igeom = IGEOM_NULL;
  iGeom_newGeom("PARALLEL", &igeom, &ierr, 8);
  IGEOM_ASSERT(ierr);
 
  // check command line arg
  const char* filename = 0;
  int nMethod = 1;
  std::string options;
  
  if (argc == 3) {
    filename = argv[1];
    nMethod = atoi(argv[2]);
  }
  else {
    if (rank == 0) {
      printf("Usage: %s [<filename>] [<send_methond>]\n", argv[0]);
      printf("No file specified.  Defaulting to: %s\n", "Moto.brep");
      printf("No send method specified.  Defaulting to: read_delete\n");
    }
    filename = "../test/Moto.brep";
    nMethod = 1;
  }
  
  // set geometry send method
  if (nMethod == 0) options += "PARALLEL=READ;";
  else if (nMethod == 1) options += "PARALLEL=READ_DELETE;";
  else if (nMethod == 2) options += "PARALLEL=BCAST;";
  else if (nMethod == 3) options += "PARALLEL=BCAST_DELETE;";
  else if (nMethod == 4) options += "PARALLEL=SCATTER;";
  else if (nMethod == 5) options += "PARALLEL=SCATTER_DELETE;";
  else if (nMethod == 6) options += "PARALLEL=READ_PARALLEL;";
  else {
    printf("Send method %d is not supported. Defaulting to: read_delete\n", nMethod);
    options = "READ_DELETE;";
  }

  // do body partitioning with round robin distribution
  options += "PARTITION=GEOM_DIMENSION;PARTITION_VAL=3;PARTITION_DISTRIBUTE;";

  struct rusage r_usage;
  long int rss_before, rss_after, size_before, size_after;
  static size_t pagesize = getpagesize();

  if (rank == 0) {
    util_getrusage(r_usage);
    rss_before = r_usage.ru_maxrss*pagesize;
    size_before = r_usage.ru_idrss*pagesize;
  }
    
  double tStart = MPI_Wtime();
  iGeom_load(igeom, filename, options.c_str(), &ierr, strlen(filename), options.length());
	    
  MPI_Barrier(MPI_COMM_WORLD); 
  double tEnd = MPI_Wtime();

  if (rank == 0) {
    util_getrusage(r_usage);
    rss_after = r_usage.ru_maxrss*pagesize;
    size_after = r_usage.ru_idrss*pagesize;

    printf("proc=%d, rss_before=%ld, size_before=%ld\n", rank, rss_before, size_before);
    printf("proc=%d, rss_after=%ld, size_after=%ld\n", rank, rss_after, size_after);
  }

  IGEOM_ASSERT(ierr);
  
  printf("Done. Geometry load time is %f\n", tEnd - tStart);
  MPI_Finalize();
}

void util_getrusage(struct rusage &r_usage)
{
  getrusage(RUSAGE_SELF, &r_usage);
  
  // this machine doesn't return rss - try going to /proc
  // print the file name to open
  if (r_usage.ru_maxrss == 0) {
    char file_str[4096], dum_str[4096];
    int file_ptr = -1, file_len;
    file_ptr = open("/proc/self/stat", O_RDONLY);
    file_len = read(file_ptr, file_str, sizeof(file_str)-1);
    if (file_len == 0) return;
    close(file_ptr);
    file_str[file_len] = '\0';
    
    // read the preceeding fields and the ones we really want...
    int dum_int;
    unsigned int dum_uint, vm_size, rss;
    static int page_size = getpagesize();
    int num_fields = sscanf(file_str, 
			    "%d " // pid
			    "%s " // comm
			    "%c " // state
			    "%d %d %d %d %d " // ppid, pgrp, session, tty, tpgid
			    "%u %u %u %u %u " // flags, minflt, cminflt, majflt, cmajflt
			    "%d %d %d %d %d %d " // utime, stime, cutime, cstime, counter, priority
			    "%u %u " // timeout, itrealvalue
			    "%d " // starttime
			    "%u %u", // vsize, rss
			    &dum_int, 
			    dum_str, 
			    dum_str, 
			    &dum_int, &dum_int, &dum_int, &dum_int, &dum_int, 
			    &dum_uint, &dum_uint, &dum_uint, &dum_uint, &dum_uint,
			    &dum_int, &dum_int, &dum_int, &dum_int, &dum_int, &dum_int, 
			    &dum_uint, &dum_uint, 
			    &dum_int,
			    &vm_size, &rss);
    if (num_fields == 24) {
      r_usage.ru_maxrss = rss/page_size;
      r_usage.ru_idrss = vm_size/page_size;
    }
  }
}
