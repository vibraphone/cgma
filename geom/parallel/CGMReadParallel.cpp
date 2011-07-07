#include <cstdio>

#include "CubitString.hpp"
#include "CubitMessage.hpp"
#include "DLList.hpp"
#include "RefEntity.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"
#include "CubitEntity.hpp"
#include "Body.hpp"
#include "CastTo.hpp"
#include "CubitUtil.hpp"
#include "CADefines.hpp"
#include "CABodies.hpp"
#include "TDParallel.hpp"
#include "CAMergePartner.hpp"
#include "TDUniqueId.hpp"

#include "TopologyBridge.hpp"
#include "GeometryQueryTool.hpp"
#include "CGMReadParallel.hpp"
#include "CGMParallelConventions.h"
#include "CGMParallelComm.hpp"
#include "CubitCompat.hpp"

const bool CGM_read_parallel_debug = false;

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

  m_round_robin = false;
  m_partition_static = false;
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
      m_partition_static = true;
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
  
  if (CGM_read_parallel_debug) {
    DEBUG_FLAG(90, CUBIT_TRUE);
    DEBUG_FLAG(138, CUBIT_TRUE);
  }
  
  // do the work by options
  bool i_read = false;
  std::vector<int>::iterator vit;
  int i;

  for (i = 1, vit = pa_vec.begin(); vit != pa_vec.end(); vit++, i++) {
    CubitStatus result = CUBIT_SUCCESS;
    switch (*vit) {
//==================
    case PA_READ:
      i_read = true;
      double tStart, tEnd;

      if (CGM_read_parallel_debug) {
	std::cout << "Reading file " << file_name << std::endl;
	tStart = MPI_Wtime();
      }
      
      result = read_entities(file_name);

      if (CUBIT_SUCCESS != result) {
	PRINT_ERROR("Reading file %s failed.\n", file_name);
	return CUBIT_FAILURE;
      }
      else if (CGM_read_parallel_debug) {
	tEnd = MPI_Wtime();
	PRINT_INFO("Read time in proc %d is %f.\n", m_rank,
		   tEnd - tStart);
	PRINT_INFO("Read done.\n");
      }
      
      break;

//==================
    case PA_BALANCE:
      if (CGM_read_parallel_debug) std::cout << "Balancing entities." << std::endl;
      if (m_round_robin) result = balance_round_robin();
      if (CUBIT_SUCCESS != result) return result;

      if (CGM_read_parallel_debug) PRINT_INFO("Balancing entities done.\n");

      break;
      
//==================     
    case PA_DELETE_NONLOCAL:
      if (CGM_read_parallel_debug) {
	PRINT_INFO("Deleting nonlocal entities.\n");
	tStart = MPI_Wtime();
      }
     
      result = delete_nonlocal_entities(reader_rank,
					partition_tag_name, 
					partition_tag_vals);
      
      if (CUBIT_SUCCESS != result) {
	PRINT_ERROR("Delete failed.\n");
	return CUBIT_FAILURE;
      }
      else if (CGM_read_parallel_debug) {
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
        if (CGM_read_parallel_debug) {
          PRINT_INFO("Broadcasting Body entities.\n");
          tStart = MPI_Wtime();
        }
        
        result = m_pcomm->broadcast_entities(reader_rank,
                                             m_pcomm->partition_body_list());
        
        if (CUBIT_SUCCESS != result) {
          PRINT_ERROR("Broadcasting Body entities failed.\n");
          return CUBIT_FAILURE;
        }
        else if (CGM_read_parallel_debug) {
          tEnd = MPI_Wtime();
          PRINT_INFO("Bcast bodies done.\n");
          PRINT_INFO("Broadcast bodies time in proc %d is %f.\n", m_proc_size,
                     tEnd - tStart);
        }

        if (!check_partition_info()) {
          PRINT_ERROR("Check partition info failed.\n");
          return CUBIT_FAILURE;
        }
      }
      
      break;

//==================      
    case PA_SCATTER:
      // do the actual scatter
      if (m_proc_size > 1) {
        if (CGM_read_parallel_debug) {
          PRINT_INFO("Scattering body entities.\n");
          tStart = MPI_Wtime();
        }
        result = m_pcomm->scatter_entities(reader_rank,
                                           m_pcomm->partition_body_list());
        
        if (CUBIT_SUCCESS != result) {
          PRINT_ERROR("Scattering body entities failed.\n");
          return CUBIT_FAILURE;
        }
        else if (CGM_read_parallel_debug) {
          tEnd = MPI_Wtime();
          PRINT_INFO("Scatter bodies done.\n");
          PRINT_INFO("Scatter bodies time in proc %d is %f.\n", m_proc_size,
                     tEnd - tStart);
        }

        if (!check_partition_info()) {
          PRINT_ERROR("Check partition info failed.\n");
          return CUBIT_FAILURE;
        }
      }
      if (CGM_read_parallel_debug) PRINT_INFO("Scatter done.\n");
      
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
  if (strstr(file_name, ".sab")) file_type = "ACIS_SAB";
  else if (strstr(file_name, ".sat")) file_type = "ACIS_SAT";
  else if (strstr(file_name, ".stp")) file_type = "STEP";
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
  CubitStatus result = CubitCompat_import_solid_model(file_name,
                                                      file_type.c_str());
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

CubitStatus CGMReadParallel::balance_round_robin()
{
  // get bodies
  int i, j;
  DLIList<RefEntity*>& body_entity_list = m_pcomm->partition_body_list();
  int n_proc = m_proc_size;
  double* loads = new double[n_proc]; // estimated loads for each processor
  double* ve_loads = new double[n_proc]; // estimated loads for each processor
  for (i = 0; i < n_proc; i++) {
    loads[i] = 0.0;
    ve_loads[i] = 0.0;
  }

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
      loads[proc] += entity->measure();

      // assign to volumes, it should be removed in future
      DLIList<RefVolume*> volumes;
      (dynamic_cast<TopologyEntity*> (entity))->ref_volumes(volumes);
      int n_vol = volumes.size();
      volumes.reset();
      for (j = 0; j < n_vol; j++) {
	RefEntity *vol = volumes.get_and_step();
	td_par = (TDParallel *) vol->get_TD(&TDParallel::is_parallel);
	if (td_par == NULL) td_par = new TDParallel(vol, NULL, &shared_procs);
      }

      // add local surface load
      DLIList<RefFace*> faces;
      (dynamic_cast<TopologyEntity*> (entity))->ref_faces(faces);
      int n_face = faces.size();
      faces.reset();
      for (j = 0; j < n_face; j++) {
        RefFace* face = faces.get_and_step();
        TopologyEntity *te = CAST_TO(face, TopologyEntity);
        if (te->bridge_manager()->number_of_bridges() < 2) {
          loads[proc] = loads[proc] + face->measure();
        }
      }
    }

    // Get all child entities
    DLIList<RefEntity*> child_list;
    RefEntity::get_all_child_ref_entities(body_entity_list, child_list);
    int n_child = child_list.size();

    // assign processors to interface entities
    child_list.reset();
    for (i = 0; i < n_child; i++) {
      entity = child_list.get_and_step();
      TopologyEntity *te = CAST_TO(entity, TopologyEntity);
      
      if (te->bridge_manager()->number_of_bridges() > 1) {
        DLIList<Body*> parent_bodies;
	DLIList<int> shared_procs;
	(dynamic_cast<TopologyEntity*> (entity))->bodies(parent_bodies);
	int n_parent = parent_bodies.size();
	
	for (j = 0; j < n_parent; j++) {
	  RefEntity *parent_vol = CAST_TO(parent_bodies.get_and_step(), RefEntity);
	  TDParallel *parent_td = (TDParallel *) parent_vol->get_TD(&TDParallel::is_parallel);
	  
	  if (parent_td == NULL) {
	    PRINT_ERROR("parent Volume has to be partitioned.");
	    return CUBIT_FAILURE;
	  }
	  shared_procs.append_unique(parent_td->get_charge_proc());
	}

	if (shared_procs.size() > 1) { // if it is interface
	  TDParallel *td_par = (TDParallel *) entity->get_TD(&TDParallel::is_parallel);
	  if (td_par == NULL) {
            int merge_id = TDUniqueId::get_unique_id(entity);
	    if (entity->entity_type_info() == typeid(RefFace)) { // face
	      if (shared_procs.size() != 2) {
		PRINT_ERROR("Error: # of shared processors of interface surface should be 2.");
		return CUBIT_FAILURE;
	      }

	      // balance interface surface loads
              if (loads[shared_procs[0]] > loads[shared_procs[1]]) {
                shared_procs.reverse();
	      }
              loads[shared_procs[0]] = loads[shared_procs[0]] + entity->measure();
	      td_par = new TDParallel(entity, NULL, &shared_procs, merge_id, 1);
	    }
	    else if (entity->entity_type_info() == typeid(RefEdge) ||
		     entity->entity_type_info() == typeid(RefVertex)) {
              // balance interface surface loads
              int min_p = shared_procs[0];
              int n_shared_proc = shared_procs.size();
              for (int i = 1; i < n_shared_proc; i++) {
                if (ve_loads[shared_procs[i]] < ve_loads[min_p]) {
                  min_p = shared_procs[i];
                }
              }
              ve_loads[min_p] = ve_loads[min_p] + entity->measure();
              shared_procs.remove(min_p);
              shared_procs.insert_first(min_p);
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
  // find bodies deleted
  int i;
  CubitStatus result;
  DLIList<RefEntity*>& body_entity_list = m_pcomm->partition_body_list();
  DLIList<RefEntity*> partition_list, delete_body_list;
  int nEntity = body_entity_list.size();
  body_entity_list.reset();

  for (i = 0; i < nEntity; i++) {
    RefEntity* entity = body_entity_list.get_and_step();
    TDParallel *td_par = (TDParallel *) entity->get_TD(&TDParallel::is_parallel);
    if (td_par == NULL) {
      PRINT_ERROR("Partitioned Volume should have TDParallel data.");
      return CUBIT_FAILURE;
    }

    if (td_par->get_charge_proc() != m_rank) { // candidate to be deleted
      // check child surfaces
      DLIList<RefFace*> face_list;
      (dynamic_cast<TopologyEntity*> (entity))->ref_faces(face_list);
      bool b_partitioned_surf = false;
      int n_face = face_list.size();
      face_list.reset();
      for (int j = 0; j < n_face; j++) {
        RefEntity* face = face_list.get_and_step();
        TDParallel *td_par_face = (TDParallel *) face->get_TD(&TDParallel::is_parallel);
        if (td_par_face != NULL) {
          DLIList<int>* shared_procs = td_par_face->get_shared_proc_list();
          int n_shared = shared_procs->size();
          shared_procs->reset();
          for (int k = 0; k < n_shared; k++) {
            if (shared_procs->get_and_step() == m_rank) {
              b_partitioned_surf = true;
              break;
            }
          }
        }
      }
      if (b_partitioned_surf) partition_list.append(entity);
      else delete_body_list.append(entity);
    }
    else partition_list.append(entity);
  }
  
  // print info
  char pre_body[100];
  DLIList<CubitEntity*> tmp_body_list;
  if (CGM_read_parallel_debug) {
    if (m_rank != reader) {
      CAST_LIST_TO_PARENT(delete_body_list, tmp_body_list);
      sprintf(pre_body, "Will delete %d Bodies: ", tmp_body_list.size());
      CubitUtil::list_entity_ids(pre_body, tmp_body_list );
    }
    std::cout << "Partitioned Body list size after delete: "
	      << partition_list.size() << std::endl;
  }

  // delete bodies
  if (m_rank != reader) {
    nEntity = delete_body_list.size();
    delete_body_list.reset();
    for (i = 0; i < nEntity; i++) {
      GeometryQueryTool::instance()->delete_RefEntity(delete_body_list.get_and_step());
    }
  }

  // update Body list in ParallelComm
  body_entity_list.clean_out();
  body_entity_list += partition_list;

  return CUBIT_SUCCESS;
}

CubitStatus CGMReadParallel::check_partition_info()
{
  int i, j;
  DLIList<RefEntity*>& body_entity_list = m_pcomm->partition_body_list();
  int nEntity = body_entity_list.size();
  body_entity_list.reset();

  for (i = 0; i < nEntity; i++) {
    RefEntity* entity = body_entity_list.get_and_step();
    TDParallel *td_par = (TDParallel *) entity->get_TD(&TDParallel::is_parallel);
    if (td_par == NULL) { // if body is not partitioned
      DLIList<RefEntity*> volumes;
      entity->get_child_ref_entities(volumes);

      // check if the first Volume is partitioned here, should be removed in future
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

      DLIList<int> s_procs;
      s_procs.append(td_par->get_charge_proc());
      td_par = new TDParallel(entity, NULL, &s_procs);
    }
  }

  return CUBIT_SUCCESS;
}
