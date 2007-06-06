#include "ParallelGeomTool.hpp"
#include "AcisMemFile.hpp"
#include "AcisQueryEngine.hpp"
#include "TSTTG_CGM.h"
#include "ProcData.hpp"

ParallelGeomTool *ParallelGeomTool::instance_ = NULL;

int ParallelGeomTool::load_parallel(TSTTG_Instance geom,
                                    const char *name, 
                                    int par_load_option) 
{
    // declare an acismemfile object
  AcisMemFile amf(AcisQueryEngine::instance());

  DLIList<RefEntity*> ref_ents, ref_ents_master;

    // if I'm the master, load the file and get the entities
  int result;
  if (ProcData::instance()->is_master()) {
    result = TSTTG_load(geom, name, NULL, 0);
    if (TSTTB_SUCCESS != result) return result;

      // get all the volumes
    TSTTG_EntityHandle *ents = NULL;
    int ents_alloc = 0, ents_size;
    result = TSTTG_getEntities(geom, 0, TSTTG_REGION, 
                               &ents, &ents_alloc, &ents_size);
    if (TSTTB_SUCCESS != result) return result;

      // make this the storage for the ent list
    ref_ents.copy_from((RefEntity**)ents, ents_size);
  }

    // now communicate the entities
  if (BCAST == par_load_option) {
      // broadcast
    result = amf.bcast_entity_list(ref_ents);
  }
  else if (BCAST_AND_DELETE == par_load_option) {
      // bcast and delete
    result = amf.bcast_and_delete_entity_list(ref_ents, ref_ents_master);
  }
  else if (SCATTER == par_load_option) {
      // scatter
    result = amf.scatter_entity_list(ref_ents);
  }
  
  return result;
}

