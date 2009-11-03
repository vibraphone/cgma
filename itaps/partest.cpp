// Test parallel iGeom

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "CGMmpi.h"
#include "iGeom.h"

#define IGEOM_ASSERT(ierr) if (ierr!=0) printf("igeom assert\n");
#define IGEOM_NULL 0

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
    if (rank == 0) printf("Usage: %s [<filename>] [<send_methond>]\n", argv[0]);
    if (argc != 1 || argc != 2) return 1;
    if (rank == 0) {
      printf("  No file specified.  Defaulting to: %s\n", "Moto.brep");
      printf("  No send method specified.  Defaulting to: read_delete\n");
    }
    filename = "../../../test_files/Moto.brep";
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

  double tStart = MPI_Wtime();
  iGeom_load(igeom, filename, options.c_str(), &ierr, strlen(filename), options.length());
	    
  MPI_Barrier(MPI_COMM_WORLD); 
  double tEnd = MPI_Wtime();

  IGEOM_ASSERT(ierr);
  
  printf("Done. Geometry load time is %f\n", tEnd - tStart);
  MPI_Finalize();
}
