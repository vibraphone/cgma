//
// Parallel meshing application
//
// This application implements parallel meshing by reading the solid
// model onto the master, then loading it onto the processors
//

#include "iostream.h"
#include "AppUtil.hpp"
#include "CGMApp.hpp"
#include "CubitString.hpp"
#include "CubitObserver.hpp"
#include "GeometryQueryTool.hpp"
#include "AcisQueryEngine.hpp"
#include "ProcData.hpp"
#include "AcisMemFile.hpp"
#include "DLIList.hpp"
#include <sys/resource.h>
#include <unistd.h>

void testsab_startup(int &argc, char **argv, int &send_method);
void read_files(const int argc, char **argv);

extern "C" void gl_cleanup(void) {}

int memFragments = 0;
int memTotal = 0;


int main(int argc, char **argv) 
{
  DEBUG_FLAG(100, CUBIT_TRUE);
  DEBUG_FLAG(137, CUBIT_TRUE);

    // Initialize MPI
  ProcData::instance()->initialize(argc, argv);

  PRINT_DEBUG_100("Processor %d initialized.\n", ProcData::instance()->myRank);

  if (ProcData::instance()->is_master() && argc == 1) {
    PRINT_ERROR("Usage: test_sab <sabt_file_name> [<sabt_file_name> ...][bcast | bcast_del | scatter]\n");
    return 0;
  }

    // Initialize the GeometryTool
  int send_method;
  testsab_startup(argc, argv, send_method);

  AcisQueryEngine *aqe = AcisQueryEngine::instance();
  GeometryQueryTool *gti = GeometryQueryTool::instance();
  assert(gti);

    // open the file(s), read the contents (if I'm the master)
  struct rusage r_usage;
  long int rss_bare, rss_geom, size_bare, size_geom,
    rss_mess, size_mess;
  static size_t pagesize = getpagesize();
  if (ProcData::instance()->is_master() &&
      CUBIT_TRUE == DEBUG_FLAG(137)) {
    AppUtil::instance()->apputil_getrusage(r_usage);
    rss_bare = r_usage.ru_maxrss*pagesize;
    size_bare = r_usage.ru_idrss*pagesize;
    
  }
  
  if (ProcData::instance()->is_master()) {
    read_files(argc, argv);
    if (CUBIT_TRUE == DEBUG_FLAG(137)) {
      AppUtil::instance()->apputil_getrusage(r_usage);
      rss_geom = r_usage.ru_maxrss*pagesize;
      size_geom = r_usage.ru_idrss*pagesize;
    }
  }
  
    // pass the model to the processors
  DLIList<RefEntity*> ref_entity_list;
  if (ProcData::instance()->is_master()) {
    gti->ref_entity_list("Body",
                         ref_entity_list,
                         CUBIT_FALSE);
    PRINT_DEBUG_100("Master: read file.\n");
  }
  else {
    PRINT_DEBUG_100("Slave %d: entering bcast.\n", ProcData::instance()->myRank);
  }

  AcisMemFile mem_file(aqe);

  // select a method to distribute geometry
  if (1 == send_method) {
    mem_file.bcast_and_delete_entity_list(ref_entity_list, argv[1]);
  }
  else if (2 == send_method) {
    mem_file.scatter_entity_list(ref_entity_list, argv[1]);
  }
  else {
    mem_file.bcast_entity_list(ref_entity_list, argv[1]);
  }
  
  PRINT_DEBUG_100("Proc = %d; number of entities = %d\n", 
             ProcData::instance()->myRank,
             ref_entity_list.size());

  if (ProcData::instance()->is_master()) {
    if (CUBIT_TRUE == DEBUG_FLAG(137)) {
      AppUtil::instance()->apputil_getrusage(r_usage);
      rss_mess = r_usage.ru_maxrss*pagesize;
      size_mess = r_usage.ru_idrss*pagesize;
    }
    
    PRINT_DEBUG_137("Max RSS before, after read & after messages (bytes) = %d, %d, %d\n",
                    rss_bare, rss_geom, rss_mess);
    PRINT_DEBUG_137("Max size before, after read & after messages (bytes) = %d, %d, %d\n",
                    size_bare, size_geom, size_mess);
  }

  delete ProcData::instance();

  return 0;
}

void read_files(const int argc, char **argv) 
{
  int i;
  CubitString file_type;
  for (i = 1; i < argc; i++) {
    if (strstr(argv[i], ".sab")) file_type = "ACIS_SAB";
    else if (strstr(argv[i], ".sat")) file_type = "ACIS_SAT";
    else if (strstr(argv[i], ".stp")) file_type = "STEP";
    else if (strstr(argv[i], ".igs")) file_type = "IGES";
    else {
      PRINT_ERROR("File type not known for file %s; skipping.\n",
                  argv[i]);
      continue;
    }
    
    FILE *fp = fopen(argv[i], "r");
    DLIList<TopologyBridge*> imported_bridges;
    GeometryQueryTool::instance()->import_solid_model(fp, argv[i],
                                                      file_type.c_str());
    fclose(fp);
  }
}

void testsab_startup(int &argc, char **argv, int &send_method)
{
    // initialize util
    AppUtil::instance()->startup(argc, argv);

   AcisQueryEngine::instance();

   // initialize CGM
   CGMApp::instance()->startup(argc, argv);

  if (!strcmp(argv[argc-1], "bcast_del")) {
    send_method = 1;
    argc--;
  }
  else if (!strcmp(argv[argc-1], "scatter")) {
    send_method = 2;
    argc--;
  }
  else if (!strcmp(argv[argc-1], "bcast")) {
    send_method = 0;
    argc--;
  }
  else
    send_method = 0;
}

  
