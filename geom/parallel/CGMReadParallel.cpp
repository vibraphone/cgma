#include <cstdio>

#include "CubitString.hpp"
#include "CubitMessage.hpp"
#include "DLList.hpp"
#include "RefEntity.hpp"
//#include "RefVolume.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include "CubitEntity.hpp"
#include "Body.hpp"
#include "CastTo.hpp"
#include "CubitUtil.hpp"
#include "CubitAttrib.hpp"
#include "CADefines.hpp"
#include "CABodies.hpp"
#include "TDParallel.hpp"
#include "CAMergePartner.hpp"

#include "TopologyBridge.hpp"
#include "GeometryQueryTool.hpp"
#include "CGMReadParallel.hpp"
#include "CGMParallelConventions.h"
#include "CGMParallelComm.hpp"

const bool debug = false;

enum CGMParallelActions {PA_READ=0, PA_BROADCAST, PA_DELETE_NONLOCAL,
			 PA_SCATTER, PA_SCATTER_DELETE, PA_BALANCE,
			 PA_CHECK_GIDS_SERIAL, PA_RESOLVE_SHARED_ENTS,
			 PA_EXCHANGE_GHOSTS};

enum CGMPartitionActions {PT_GEOM_DIM=0, PT_PAR_PART};

const char *CGMParallelActionsNames[] = {
  "PARALLEL READ",
  "PARALLEL BROADCAST", 
  "PARALLEL DELETE NONLOCAL",
  "PARALLEL SCATTER",
  "PARALLEL CHECK_GIDS_SERIAL",
  "PARALLEL GET_FILESET_ENTS",
  "PARALLEL RESOLVE_SHARED_ENTS",
  "PARALLEL EXCHANGE_GHOSTS"
};

const char* CGMReadParallel::CGMparallelOptsNames[] = { "NONE", "READ", "READ_DELETE", "BCAST", 
							"BCAST_DELETE", "SCATTER", "SCATTER_DELETE",
							"READ_PARALLEL", "FORMAT", "", 0 };

const char* CGMReadParallel::CGMpartitionOptsNames[] = { "NONE", "GEOM_DIMENSION",
							 "PARARELL_PARTITION", "", 0 };

CGMReadParallel::CGMReadParallel(GeometryQueryTool* gqt, CGMParallelComm *pc)
  : m_gqt(gqt), m_pcomm(pc) 
{
  if (!m_pcomm) {
    m_pcomm = CGMParallelComm::get_pcomm(0);
    if (NULL == m_pcomm) m_pcomm = new CGMParallelComm();
  }

  m_scatter = false;
  m_rank = m_pcomm->proc_config().proc_rank();
  m_proc_size = m_pcomm->proc_config().proc_size();
}

CubitStatus CGMReadParallel::load_file(const char *file_name,
					const char *options,
					const char* set_tag_name,
					const int* set_tag_values,
					int num_set_tag_values) 
{
  CGMFileOptions opts(options);

  // Get parallel settings
  int parallel_mode;
  CGMFOErrorCode result = opts.match_option("PARALLEL", CGMparallelOptsNames, 
					    parallel_mode);
  if (FO_FAILURE == result) {
    PRINT_ERROR( "Unexpected value for 'PARALLEL' option\n" );
    return CUBIT_FAILURE;
  }
  else if (FO_ENTITY_NOT_FOUND == result) {
    parallel_mode = 0;
  }

  bool surf_partition = false;
  std::string partition_tag_name;
  std::vector<int> partition_tag_vals;
  int geom_dim;

  // Get partition tag value(s), if any, and whether they're to be
  result = opts.get_ints_option("PARTITION_VAL", partition_tag_vals);

  // Get partition setting
  result = opts.get_option("PARTITION", partition_tag_name);

  if (FO_ENTITY_NOT_FOUND == result || partition_tag_name.empty()) {
    partition_tag_name = "GEOM_DIMENSION";
    m_round_robin = true;
    geom_dim = 3;
  }
  else {
    // use geom dimension for partition
    if (partition_tag_name == "GEOM_DIMENSION") {
      int geom_dim = 0;
      for (std::vector<int>::iterator pit = partition_tag_vals.begin(); 
	   pit != partition_tag_vals.end(); pit++) {
	geom_dim = *pit;
	if (geom_dim == 2) surf_partition = true; // body & surface partition
	else if (geom_dim == 3) surf_partition = false; // body partition only
	else {
	  PRINT_ERROR("Geometry dimension %d is not supported.\n", geom_dim);
	  return CUBIT_FAILURE;
	}
      }
    }
    // static partition, use chaco
    else if (partition_tag_name == "PAR_PARTITION_STATIC") {

    }
    // dynamic partition, use zoltan
    else if (partition_tag_name == "PAR_PARTITION_DYNAMIC") {

    }

    // round-robin
    result = opts.get_null_option("PARTITION_DISTRIBUTE");
    if (FO_SUCCESS == result) m_round_robin = true;
  }

  // get MPI IO processor rank
  int reader_rank;
  result = opts.get_int_option("MPI_IO_RANK", reader_rank);
  if (FO_ENTITY_NOT_FOUND == result) reader_rank = 0;
  else if (FO_SUCCESS != result) {
    PRINT_ERROR( "Unexpected value for 'MPI_IO_RANK' option\n" );
    return CUBIT_FAILURE;
  }
  m_pcomm->proc_config().set_master(reader_rank); // set master processor
  bool reader = (reader_rank == (int) m_rank);

  // now that we've parsed all the parallel options, make an instruction
  // queue
  std::vector<int> pa_vec;
  switch (parallel_mode) {

  case POPT_READ:
    pa_vec.push_back(PA_READ);
    pa_vec.push_back(PA_BALANCE);
    break;
  
  case POPT_DEFAULT:
  case POPT_READ_DELETE:
    pa_vec.push_back(PA_READ);
    pa_vec.push_back(PA_BALANCE);
    pa_vec.push_back(PA_DELETE_NONLOCAL);
    break;

  case POPT_BCAST:
    if (reader) {
      pa_vec.push_back(PA_READ);
      pa_vec.push_back(PA_BALANCE);
    }
    pa_vec.push_back(PA_BROADCAST);
    break;
    
  case POPT_BCAST_DELETE:
    if (reader) {
      pa_vec.push_back(PA_READ);
      pa_vec.push_back(PA_BALANCE);
    }
    pa_vec.push_back(PA_BROADCAST);
    pa_vec.push_back(PA_DELETE_NONLOCAL);
    break;
    
  case PORT_SCATTER:
    if (reader) {
      pa_vec.push_back(PA_READ);
      pa_vec.push_back(PA_BALANCE);
    }
    pa_vec.push_back(PA_SCATTER);
    m_scatter = true;
    break;
  
  case POPT_FORMAT:
    PRINT_ERROR( "Access to format-specific parallel read not implemented.\n");
    return CUBIT_FAILURE;

  case POPT_READ_PARALLEL:
    PRINT_ERROR( "Partitioning for PARALLEL=READ_PARALLEL not supported yet.\n");
    return CUBIT_FAILURE;

  default:
    return CUBIT_FAILURE;
  }

  return load_file(file_name, parallel_mode, 
                   partition_tag_name,
                   partition_tag_vals, pa_vec, opts,
                   set_tag_name, set_tag_values,
		   num_set_tag_values,
                   reader_rank, surf_partition
		   );
}

CubitStatus CGMReadParallel::load_file(const char *file_name,
				       int parallel_mode, 
				       std::string &partition_tag_name, 
				       std::vector<int> &partition_tag_vals, 
				       std::vector<int> &pa_vec,
				       const CGMFileOptions &opts,
				       const char* set_tag_name,
				       const int* set_tag_values,
				       const int num_set_tag_values,
				       const int reader_rank,
				       const bool surf_partition
				       )
{
  

  // actuate CA_BODIES and turn on auto flag for other attributes
  CGMApp::instance()->attrib_manager()->register_attrib_type(CA_BODIES, "bodies", "BODIES",
							     CABodies_creator, CUBIT_TRUE,
							     CUBIT_TRUE, CUBIT_TRUE, CUBIT_TRUE,
							     CUBIT_TRUE, CUBIT_FALSE);
  CGMApp::instance()->attrib_manager()->auto_flag(CUBIT_TRUE);
  
  if (debug) {
    DEBUG_FLAG(90, CUBIT_TRUE);
    DEBUG_FLAG(138, CUBIT_TRUE);
  }
  
  // do the work by options
  bool i_read = false;
  std::vector<int>::iterator vit;
  int i;
  //DLIList<RefEntity*> surf_entity_list, body_entity_list;
  for (i = 1, vit = pa_vec.begin(); vit != pa_vec.end(); vit++, i++) {
    CubitStatus result = CUBIT_SUCCESS;
    switch (*vit) {
//==================
    case PA_READ:
      i_read = true;
      double tStart, tEnd;

      if (debug) {
	std::cout << "Reading file " << file_name << std::endl;
	tStart = MPI_Wtime();
      }
      
      result = read_entities(file_name);

      if (CUBIT_SUCCESS != result) {
	PRINT_ERROR("Reading file %s failed.\n", file_name);
	return CUBIT_FAILURE;
      }
      else if (debug) {
	tEnd = MPI_Wtime();
	PRINT_INFO("Read time in proc %d is %f.\n", m_rank,
		   tEnd - tStart);
	PRINT_INFO("Read done.\n");
      }
      
      break;

//==================
    case PA_BALANCE:
      if (debug) std::cout << "Balancing entities." << std::endl;
      
      result = balance();
      if (CUBIT_SUCCESS != result) return result;

      if (debug) PRINT_INFO("Balancing entities done.\n");

      break;
      
//==================     
    case PA_DELETE_NONLOCAL:
      if (debug) {
	PRINT_INFO("Deleting nonlocal entities.\n");
	tStart = MPI_Wtime();
      }
     
      result = delete_nonlocal_entities(reader_rank,
					partition_tag_name, 
					partition_tag_vals);
					//surf_entity_list,
					//body_entity_list,
					//round_robin);
      
      if (CUBIT_SUCCESS != result) {
	PRINT_ERROR("Delete failed.\n");
	return CUBIT_FAILURE;
      }
      else if (debug) {
	tEnd = MPI_Wtime();
	PRINT_INFO("Delete done.\n");
	PRINT_INFO("Delete time in proc %d is %f.\n", m_rank,
		   tEnd - tStart);
      }
      break;

//==================      
    case PA_BROADCAST:
      // do the actual broadcast; if single-processor, ignore error
      if (m_proc_size > 1) {
	//if (body_partition) {
	  if (debug) {
	    PRINT_INFO("Broadcasting Body entities.\n");
	    tStart = MPI_Wtime();
	  }

	  result = m_pcomm->broadcast_entities(reader_rank,
					       m_pcomm->partition_body_list());
	  
	  if (CUBIT_SUCCESS != result) {
	    PRINT_ERROR("Broadcasting Body entities failed.\n");
	    return CUBIT_FAILURE;
	  }
	  else if (debug) {
	    tEnd = MPI_Wtime();
	    PRINT_INFO("Bcast bodies done.\n");
	    PRINT_INFO("Broadcast bodies time in proc %d is %f.\n", m_proc_size,
		       tEnd - tStart);
	  }
      }
      
      break;

//==================      
    case PA_SCATTER:
      // do the actual scatter
      if (m_proc_size > 1) {
	  if (debug) {
	    PRINT_INFO("Scattering body entities.\n");
	    tStart = MPI_Wtime();
	  }
	  result = m_pcomm->scatter_entities(reader_rank,
					     m_pcomm->partition_body_list());
	  
	  if (CUBIT_SUCCESS != result) {
	    PRINT_ERROR("Scattering body entities failed.\n");
	    return CUBIT_FAILURE;
	  }
	  else if (debug) {
	    tEnd = MPI_Wtime();
	    PRINT_INFO("Scatter bodies done.\n");
	    PRINT_INFO("Scatter bodies time in proc %d is %f.\n", m_proc_size,
		       tEnd - tStart);
	  }
      }
      if (debug) PRINT_INFO("Scatter done.\n");
      
      break;

//==================    
    default:
      return CUBIT_FAILURE;
    }
  }

  return CUBIT_SUCCESS;
}

CubitStatus CGMReadParallel::read_entities(const char* file_name)
{
  // check file type
  CubitString file_type;
  //if (strstr(file_name, ".sab")) file_type = "ACIS_SAB";
  //else if (strstr(file_name, ".sat")) file_type = "ACIS_SAT";
  if (strstr(file_name, ".stp")) file_type = "STEP";
  else if (strstr(file_name, ".igs")) file_type = "IGES";
  else if (strstr(file_name, ".occ") ||
	   strstr(file_name, ".OCC") ||
	   strstr(file_name, ".brep") ||
	   strstr(file_name, ".BREP")) file_type = "OCC";
  else {
    PRINT_ERROR("File type not known for file %s; skipping.\n", file_name);
    return CUBIT_FAILURE;
  }

  // import solid model
  CubitStatus result = m_gqt->import_solid_model(file_name, file_type.c_str());
  if (CUBIT_SUCCESS != result) {
    PRINT_ERROR("Reading file %s failed.\n", file_name);
    return CUBIT_FAILURE;
  }
  
  // get body entities
  DLIList<RefEntity*>& body_entity_list = m_pcomm->partition_body_list();
  body_entity_list.clean_out();
  result = m_gqt->ref_entity_list("body", body_entity_list, CUBIT_FALSE);
  if (CUBIT_SUCCESS != result) {
    PRINT_ERROR("Getting Body entities failed.\n");
    return result;
  }

  return result;
}

CubitStatus CGMReadParallel::balance()
{
  // get bodies
  int i, j;
  DLIList<RefEntity*>& body_entity_list = m_pcomm->partition_body_list();
  int n_proc = m_proc_size;
  double* loads = new double[n_proc]; // estimated loads for each processor
  for (i = 0; i < n_proc; i++) loads[i] = 0.0;

  if (m_round_robin) { // round-robin case
    int n_entity = body_entity_list.size();
    int n_entity_proc = n_entity/n_proc; // # of entities per processor
    int i_entity_proc = n_entity_proc; // entity index limit for each processor
    int proc = 0;
    RefEntity* entity;

    // assign processors to bodies
    body_entity_list.reset();
    for (i = 0; i < n_entity; i++) {
      if (i == i_entity_proc) {
	proc++;
	if (proc < n_proc) i_entity_proc += n_entity_proc;
	else {
	  proc %= n_proc;
	  i_entity_proc++;
	}
      }

      // assign to bodies
      entity = body_entity_list.get_and_step();
      DLIList<int> shared_procs;
      shared_procs.append(proc);
      TDParallel *td_par = (TDParallel *) entity->get_TD(&TDParallel::is_parallel);
      if (td_par == NULL) td_par = new TDParallel(entity, NULL, &shared_procs);
      //loads[proc] += entity->measure();

      // assign to volumes
      DLIList<RefVolume*> volumes;
      (dynamic_cast<TopologyEntity*> (entity))->ref_volumes(volumes);
      int n_vol = volumes.size();
      volumes.reset();
      for (j = 0; j < n_vol; j++) {
	RefEntity *vol = volumes.get_and_step();
	td_par = (TDParallel *) vol->get_TD(&TDParallel::is_parallel);
	if (td_par == NULL) td_par = new TDParallel(vol, NULL, &shared_procs);
	loads[proc] += vol->measure();
      }
    }

    // Get all child entities
    DLIList<RefEntity*> child_list;
    RefEntity::get_all_child_ref_entities(body_entity_list, child_list);
    int n_child = child_list.size();

    // assign processors to interface entities
    for (i = 0; i < n_child; i++) {
      entity = child_list.get_and_step();
      CubitAttrib* att = entity->get_cubit_attrib(CA_MERGE_PARTNER,
						  CUBIT_FALSE);
      
      if (att != NULL) { // if it is shared entity
	DLIList<Body*> parent_bodies;
	//DLIList<RefVolume*> parent_volumes;
	DLIList<int> shared_procs;
	(dynamic_cast<TopologyEntity*> (entity))->bodies(parent_bodies);
	//(dynamic_cast<TopologyEntity*> (entity))->ref_volumes(parent_volumes);
	int n_parent = parent_bodies.size();
	//int n_parent = parent_volumes.size();
	
	for (j = 0; j < n_parent; j++) {
	  RefEntity *parent_body = parent_bodies.get_and_step();
	  //RefVolume *parent_vol = parent_volumes.get_and_step();
	  //RefEntity *parent_vol = CAST_TO(parent_volumes.get_and_step(), RefEntity);
	  TDParallel *parent_td = (TDParallel *) parent_body->get_TD(&TDParallel::is_parallel);
	  //TDParallel *parent_td = (TDParallel *) parent_vol->get_TD(&TDParallel::is_parallel);
	  
	  if (parent_td == NULL) {
	    PRINT_ERROR("parent Volume has to be partitioned.");
	    return CUBIT_FAILURE;
	  }
	  shared_procs.append_unique(parent_td->get_charge_proc());
	}

	if (shared_procs.size() > 1) { // if it is interface
	  TDParallel *td_par = (TDParallel *) entity->get_TD(&TDParallel::is_parallel);
	  if (td_par == NULL) {
	    CAMergePartner *camp_ptr = CAST_TO(att, CAMergePartner);
	    int merge_id = camp_ptr->merge_id();
	    if (entity->entity_type_info() == typeid(RefFace)) { // face
	      if (shared_procs.size() != 2) {
		PRINT_ERROR("Error: # of shared processors of interface surface should be 2.");
		return CUBIT_FAILURE;
	      }
	      // make the first shared processor is charging mesh
	      if (loads[shared_procs[0]] > loads[shared_procs[1]]) {
		int temp_proc = shared_procs.pop();
		shared_procs.append(temp_proc);
	      }
	      td_par = new TDParallel(entity, NULL, &shared_procs, merge_id, 1);
	    }
	    else if (entity->entity_type_info() == typeid(RefEdge) ||
		     entity->entity_type_info() == typeid(RefVertex)) {
	      td_par = new TDParallel(entity, NULL, &shared_procs, merge_id, 1);
	    }
	  }
	}
      }
    }
  }

  return CUBIT_SUCCESS;
}

CubitStatus CGMReadParallel::delete_nonlocal_entities(int reader,
						      std::string &ptag_name,
						      std::vector<int> &ptag_vals)
{
  // find volumes deleted
  int i;
  CubitStatus result;
  DLIList<RefEntity*>& body_entity_list = m_pcomm->partition_body_list();
  //DLIList<RefEntity*>& vol_entity_list = m_pcomm->partition_vol_list();
  DLIList<RefEntity*> partition_list, delete_body_list;
  int nEntity = body_entity_list.size();
  body_entity_list.reset();
  //int nEntity = vol_entity_list.size();
  //vol_entity_list.reset();
  for (i = 0; i < nEntity; i++) {
    RefEntity* entity = body_entity_list.get_and_step();
    TDParallel *td_par = (TDParallel *) entity->get_TD(&TDParallel::is_parallel);

    if (td_par == NULL) {
      //RefEntity* entity = vol_entity_list.get_and_step();
      DLIList<RefEntity*> volumes;
      entity->get_child_ref_entities(volumes);
      
      // check if the first Volume is partitioned here
      volumes.reset();
      RefEntity *vol = volumes.get();
      if (vol == NULL || vol->entity_type_info() != typeid(RefVolume)) {
	PRINT_ERROR("Partitioned Body should have at least one Volume.");
	return CUBIT_FAILURE;
      }
      td_par = (TDParallel *) vol->get_TD(&TDParallel::is_parallel);
      
      if (td_par == NULL) {
	PRINT_ERROR("Partitioned Volume should have TDParallel data.");
	return CUBIT_FAILURE;
      }
    }

    //if (td_par->get_charge_proc() != m_rank) delete_body_list.append(entity);
    if (td_par->get_charge_proc() != m_rank) delete_body_list.append(entity);
    else partition_list.append(entity);
  }
  
  // delete bodies
  if (m_rank != reader) {
    nEntity = delete_body_list.size();
    //nEntity = delete_vol_list.size();
    for (i = 0; i < nEntity; i++) {
      GeometryQueryTool::instance()->delete_RefEntity(delete_body_list[i]);
      //GeometryQueryTool::instance()->delete_RefEntity(delete_vol_list[i]);
    }
  }

  // update volume list in ParallelComm
  body_entity_list.clean_out();
  body_entity_list += partition_list;
  //vol_entity_list.clean_out();
  //vol_entity_list += partition_list;
  
  // print info
  char pre_body[100];
  DLIList<CubitEntity*> tmp_body_list;
  if (debug) {
    if (m_rank != reader) {
      CAST_LIST_TO_PARENT(delete_body_list, tmp_body_list);
      sprintf(pre_body, "Deleted %d Bodies: ", tmp_body_list.size());
      CubitUtil::list_entity_ids(pre_body, tmp_body_list );
    }
    std::cerr << "Partitioned Body list size after delete: "
	      << partition_list.size() << std::endl;
  }

  return CUBIT_SUCCESS;
}
