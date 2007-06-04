//#include <vector>
#include <set>
#include "string.h"
#include "CubitDefines.h"
#include "CGMMemFile.hpp"
#include "CubitMessage.hpp"
#include "GeometryQueryTool.hpp"
#include "TDUniqueId.hpp"
#include "Body.hpp"
#include "ProcData.hpp"
#include "TDParallel.hpp"

//std::vector<int> CGMMemFile::bodyToProc;
DLIList <int> CGMMemFile::bodyToProc;
//int *CGMMemFile::bodyToProc = NULL;
//int CGMMemFile::nBodyToProc = 0

/*
int larger_vol(RefEntity* &x, RefEntity* &y)
{
  double x_vol = x->measure();
  double y_vol = y->measure();
  if (x_vol > y_vol)
    return -1;
  else if (x_vol < y_vol)
    return 1;
  else
    return 0;
}
*/
CGMMemFile::CGMMemFile(void* buffer,
                       size_t bufSize)
{
  PRINT_DEBUG_100("Initializing memfile in proc %d, size %d.\n", 
                  ProcData::instance()->myRank, bufSize);
  m_pBuffer = (unsigned char *) buffer;
  m_sizeBuffer = bufSize;
  m_currentPosition = 0;
  //balanceSurf = 0;
  //totalLoads = NULL;
  //balancedLists = NULL;
}

CGMMemFile::~CGMMemFile() 
{
  // delete the buffer
  delete m_pBuffer;
}

CubitStatus CGMMemFile::bcast_buffer() 
{
    //- broadcasts the buffer contained in this object
    // get a ProcData object, then simply broadcast
  ProcData *proc_data = ProcData::instance();
  if (proc_data->is_master()) {
    PRINT_DEBUG_100("Broadcasting buffer size from master.\n");
#ifdef USE_MPI
    MPI_Bcast(&m_sizeBuffer, 1, MPI_INT, proc_data->masterRank, 
              MPI_COMM_WORLD);
#endif
    PRINT_DEBUG_100("Broadcasting buffer from master, %d bytes.\n", m_sizeBuffer);
#ifdef USE_MPI
    MPI_Bcast(m_pBuffer, m_sizeBuffer, MPI_BYTE, proc_data->masterRank, 
              MPI_COMM_WORLD);
#endif
  }
  else {
    int this_size;
    PRINT_DEBUG_100("Broadcasting buffer size from proc %d.\n", ProcData::instance()->myRank);
#ifdef USE_MPI
    MPI_Bcast(&this_size, 1, MPI_INT, proc_data->masterRank, 
              MPI_COMM_WORLD);
#endif
    PRINT_DEBUG_100("Processor %d: received size of %d.\n", proc_data->myRank, this_size);
    check_size(this_size);
    PRINT_DEBUG_100("Broadcasting buffer from proc %d, %d bytes.\n", 
                    ProcData::instance()->myRank, this_size);
#ifdef USE_MPI
    MPI_Bcast(m_pBuffer, this_size, MPI_BYTE, proc_data->masterRank, 
              MPI_COMM_WORLD);
#endif
  }
  
  return CUBIT_SUCCESS;
}

CubitStatus CGMMemFile::bcast_entity_list(DLIList<RefEntity*> &ref_entity_list) 
{
    //- broadcasts the entity list; root takes list as input, 
    //- others write as output
  CubitStatus status;
  double twrite, tbcast;
  if (ProcData::instance()->is_master()) {
    if (DEBUG_FLAG(137)) {
#ifdef USE_MPI
      twrite = MPI_Wtime();
#endif
      tbcast = twrite;
    }
    write_refentity_list(ref_entity_list);
    if (CUBIT_TRUE == DEBUG_FLAG(137))
#ifdef USE_MPI
      twrite = MPI_Wtime() - twrite;
#endif
  }
  
  status = bcast_buffer();
  
  if (!ProcData::instance()->is_master()) {
    read_refentity_list(ref_entity_list);
  }

/*
  if (CUBIT_TRUE == DEBUG_FLAG(137)) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (ProcData::instance()->is_master()) {
      tbcast = MPI_Wtime() - tbcast;
      GeometryQueryTool *gti = GeometryQueryTool::instance();
      
      PRINT_DEBUG_137("buffer size (bytes), time to write buffer, total time to bcast (sec) = %d, %f, %f\n",
                      m_sizeBuffer, twrite, tbcast);
      PRINT_DEBUG_137("Complexity: vertices, edges, faces, bodies = %d, %d, %d, %d\n",
                      gti->num_ref_vertices(), gti->num_ref_edges(), gti->num_ref_faces(), 
                      gti->num_bodies());
      PRINT_DEBUG_137("Stats: %s, %d, %d, %f, %f, %d, %d, %d, %d\n",
                      filename, ProcData::instance()->numProcs, m_sizeBuffer, twrite, tbcast, 
                      gti->num_ref_vertices(), gti->num_ref_edges(), gti->num_ref_faces(), 
                      gti->num_bodies());
    }
  }
*/
  return status;
}

CubitStatus CGMMemFile::bcast_and_delete_entity_list(DLIList<RefEntity*> &ref_entity_list,
                                                     DLIList<RefEntity*> &ref_entity_list_master) 
{
  // reset the attribImporteds flags to facilitate attribute reporting
  //CubitAttrib::clear_attrib_importeds();

  //- broadcasts the entity list and delete the rest which this process doesn't need
  CubitStatus status;
  int nEntity, totalRefSize, endEntity, index, nBarEntity, restEntity;
  double tbcastDelete;
  ProcData *proc_data = ProcData::instance();
  int nProcs = proc_data->numProcs;
  int myid = proc_data->myRank;

  if (proc_data->is_master()) {
    totalRefSize = get_refentity_list_size(ref_entity_list);

    if (DEBUG_FLAG(137)) {
#ifdef USE_MPI
      tbcastDelete = MPI_Wtime();
#endif
    }
    write_refentity_list(ref_entity_list);
  }

  status = bcast_buffer();

  if (!proc_data->is_master()) {
    read_refentity_list(ref_entity_list);
    
    nEntity = ref_entity_list.size();
    nBarEntity = nEntity/nProcs;
    restEntity = nEntity%nProcs;
    
    // distribute with Round-Robin
    // tjt - don't need the 'else' part of the next if test,
    // since if restEntity == 0 you'll get the same results (i.e. the 'else' clause
    // of the next (nested) if will be the one evaluated
    // H.J.K - I thought another 'else' part can reduce some flops.
    //if (restEntity) {
    if (myid < restEntity) {
      index = myid*(nBarEntity+1) + 1;
      endEntity = (myid+1)*(nBarEntity+1) + 1;
    }
    else {
      index = restEntity*(nBarEntity+1) + (myid-restEntity)*nBarEntity + 1;
      endEntity = restEntity*(nBarEntity+1) + (myid-restEntity+1)*nBarEntity + 1;
    }
    //}
    /* else {
       index = myid*nBarEntity + 1;
       endEntity = (myid+1)*nBarEntity + 1;
       }*/
    
    // tjt - huh?  Where do you delete the entities?
    // H.J.K - I changed code and deleted ref_entity_list directly without temp_list
    ref_entity_list.reset();
    for (int i = 1; i < index; i++) {
      ref_entity_list.remove();
    }
    
    ref_entity_list.step(endEntity-index);
    
    for (int i = endEntity; i < nEntity+1; i++) {
      ref_entity_list.remove();
    }

    // report attribs imported
    //CubitAttrib::report_attrib_importeds();
  }
  else {
    nEntity = ref_entity_list.size();
    nBarEntity = nEntity/nProcs;
    restEntity = nEntity%nProcs;

    ref_entity_list_master.clean_out();
    ref_entity_list.reset();

    for (int j = 0; j < nBarEntity; j++)
      ref_entity_list_master.append(ref_entity_list.get_and_step());

    if (restEntity > 0)
      ref_entity_list_master.append(ref_entity_list.get_and_step());

    //cout << "ref_entity_list_master.size()="<< ref_entity_list_master.size() << endl;
  }

/*  
  if (CUBIT_TRUE == DEBUG_FLAG(137)) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (ProcData::instance()->is_master()) {
      tbcastDelete = MPI_Wtime() - tbcastDelete;
      GeometryQueryTool *gti = GeometryQueryTool::instance();
      
      PRINT_DEBUG_137("total time to bcast and delete(sec) = %f\n", tbcastDelete);
      PRINT_DEBUG_137("total reference entity list size (bytes) = %d\n", totalRefSize);
      PRINT_DEBUG_137("Complexity: vertices, edges, faces, bodies = %d, %d, %d, %d\n",
                      gti->num_ref_vertices(), gti->num_ref_edges(), gti->num_ref_faces(), 
                      gti->num_bodies());
      double twrite = 0.0;
      PRINT_DEBUG_137("Stats: %s, %d, %d, %f, %f, %d, %d, %d, %d\n",
                      filename, ProcData::instance()->numProcs, m_sizeBuffer, twrite, tbcastDelete, 
                      gti->num_ref_vertices(), gti->num_ref_edges(), gti->num_ref_faces(),                       gti->num_bodies());
    }
  }
*/
  return status;
}
/*
CubitStatus CGMMemFile::send_body_to_procs(DLIList<RefEntity*> &ref_entity_list) 
{
  int nEntity;
  DLIList<RefEntity*> temp_list, temp_ref_list;

  ProcData *proc_data = ProcData::instance();

  int nProcs = proc_data->numProcs;

  if (proc_data->is_master()) {
    nEntity = ref_entity_list.size();
    int nBarEntity = nEntity/nProcs;
    int restEntity = nEntity%nProcs;
    int nEndEntity = nBarEntity + restEntity;
    
    ref_entity_list.reset();
    bodyToProc.clean_out();

    // make temporary lists to contain geometry information for each processors
    for (int i = 0; i < nProcs; i++) {
      for ( int j = 0; j < nBarEntity; j++) {
	RefEntity* body_ptr = ref_entity_list.get_and_step();
	temp_list.append(body_ptr);
	bodyToProc.append(TDUniqueId::get_unique_id(body_ptr));
	bodyToProc.append(i);
      }
      
      if (i < restEntity) {
	RefEntity* body_ptr = ref_entity_list.get_and_step();
	temp_list.append(body_ptr);
	bodyToProc.append(TDUniqueId::get_unique_id(body_ptr));
	bodyToProc.append(i);
      }

      temp_list.clean_out();
    }
  }

  // broadcast number of bodies
  MPI_Bcast(&nEntity, 1, MPI_INT, proc_data->masterRank, MPI_COMM_WORLD);

  int body_to_proc[2*nEntity];

  if (proc_data->is_master()) // make temporary array
    bodyToProc.copy_to(body_to_proc);
  
  MPI_Bcast(body_to_proc, 2*nEntity, MPI_INT, proc_data->masterRank, MPI_COMM_WORLD);

  if (!proc_data->is_master()) // copy from temporary array
    bodyToProc.copy_from(body_to_proc, 2*nEntity);
  
  return CUBIT_SUCCESS;
}
*/

CubitStatus CGMMemFile::scatter_entity_list(DLIList<RefEntity*> &ref_entity_list) 
{
  // scatter exact amount of geometry information to each processors
  CubitStatus status;
  double tScatter;
  int mySendCount;
  int nEntity;
  DLIList<RefEntity*> temp_list, temp_ref_list;

  ProcData *proc_data = ProcData::instance();
  if (DEBUG_FLAG(137) && proc_data->is_master())
#ifdef USE_MPI
    tScatter = MPI_Wtime();
#endif

  int nProcs = proc_data->numProcs;
  int *sendCounts = new int[nProcs];
  int *displacements = new int[nProcs];
  int nBarEntity;
  int restEntity;
  int nEndEntity;
  displacements[0] = 0;

  if (proc_data->is_master()) {
    nEntity = ref_entity_list.size();
    nBarEntity = nEntity/nProcs;
    restEntity = nEntity%nProcs;
    nEndEntity = nBarEntity + restEntity;
    
    ref_entity_list.reset();
    int sum = 0;

    // make temporary lists to contain geometry information for each processors
    for (int i = 0; i < nProcs; i++) {
      for ( int j = 0; j < nBarEntity; j++) {
	RefEntity* body_ptr = ref_entity_list.get_and_step();
	temp_list.append(body_ptr);
      }
      
      if (i < restEntity) {
	RefEntity* body_ptr = ref_entity_list.get_and_step();
	temp_list.append(body_ptr);
      }

      sendCounts[i] = get_refentity_list_size(temp_list);
      sum += sendCounts[i];
      temp_list.clean_out();
    }

      // check the size of the buffer and resize if necessary
    check_size(sum);
    
      // now actually append the information
    ref_entity_list.reset();
    for (int i = 0; i < nProcs; i++) {
      for ( int j = 0; j < nBarEntity; j++) {
	temp_list.append(ref_entity_list.get_and_step());
      }
      
      if (i < restEntity) {
	temp_list.append(ref_entity_list.get_and_step());
      }

      append_refentity_list(temp_list, sendCounts[i]);
      temp_list.clean_out();
    }
  }

  // broadcast buffer size array
  PRINT_DEBUG_100("Broadcasting buffer size array from master.\n");
#ifdef USE_MPI
  MPI_Bcast(sendCounts, nProcs, MPI_INT, proc_data->masterRank, MPI_COMM_WORLD);
#endif
  
  for (int i = 1; i < nProcs; i++) {
    displacements[i] = displacements[i-1] + sendCounts[i-1];
  }
  
  mySendCount = sendCounts[proc_data->myRank];

  if (!proc_data->is_master()) check_size(mySendCount);

  PRINT_DEBUG_100("Scattering buffer from master.\n");

  // scatter geometry
#ifdef USE_MPI
  MPI_Scatterv(m_pBuffer, sendCounts, displacements, MPI_BYTE, m_pBuffer, 
	       mySendCount, MPI_BYTE, proc_data->masterRank, MPI_COMM_WORLD);
#endif

  if (!proc_data->is_master()) {
    status = read_refentity_list(ref_entity_list);
  }
  /*
  ref_entity_list.reset();
  for (int i = 0; i < ref_entity_list.size(); i++) {
    DLIList<RefEntity*> entities;
    RefEntity *entity = ref_entity_list.get_and_step();
    entity->get_all_child_ref_entities(entities);
    //cout << "my_id=" << ProcData::instance()->myRank 
    // << " num of children=" << entities.size() << endl;
    //entities.reset();
    //for (int j = 0; j < entities.size(); j++) {
    //cout << "my_id=" << ProcData::instance()->myRank 
    //   << " measure=" << entities.get_and_step()->measure() << endl;
    //}
  }
  */
  /*
  // delete unneeded part in master node
  if (proc_data->is_master()) {
    int index, endEntity;
    if (0 < restEntity) {
      index = 1;
      endEntity = nBarEntity + 2;
    }
    else {
      index = restEntity*(nBarEntity+1) + 1;
      endEntity = nBarEntity + 1;
    }
    ref_entity_list.reset();
    for (int i = 1; i < index; i++) {
      ref_entity_list.remove();
    }
    
    ref_entity_list.step(endEntity-index);
    
    for (int j = endEntity; j < nEntity+1; j++) {
      ref_entity_list.remove();
    }
  }
  */
  /*
  int body_to_proc[2*nEntity];

  if (proc_data->is_master()) // make temporary array
    bodyToProc.copy_to(body_to_proc);
  
  MPI_Bcast(body_to_proc, 2*nEntity, MPI_INT, proc_data->masterRank, MPI_COMM_WORLD);

  if (!proc_data->is_master()) // copy from temporary array
    bodyToProc.copy_from(body_to_proc, 2*nEntity);

  cout << "myid=" << proc_data->myRank << " bodyToProc.size() = " << bodyToProc.size() << endl;

  DLIList <int> *body_to_proc_ptr = get_body_to_proc();

  cout << "myid=" << proc_data->myRank << " body_to_proc_ptr->size() = " << body_to_proc_ptr->size() << endl;

  if (!proc_data->is_master())
    bodyToProc.resize(2*nEntity);

  MPI_Bcast(bodyToProc, nBodyToProc, MPI_INT, proc_data->masterRank, MPI_COMM_WORLD);

  cout << "myid=" << proc_data->myRank << " bodyToProc.size() = " << bodyToProc.size() << endl;
  
  std::vector<int> *body_to_proc_ptr = get_body_to_proc();
  
  cout << "myid=" << proc_data->myRank << " body_to_proc_ptr->size() = " << body_to_proc_ptr->size() << endl;
  cout << "myid=" << proc_data->myRank << " nEntity pointer = " << &nEntity << endl;
  cout << "myid=" << proc_data->myRank << " body_to_proc_ptr = " << body_to_proc_ptr << endl;
  cout << "myid=" << proc_data->myRank << " bodyToProcPtr = " << &bodyToProc << endl;
  */  
  /*
  // broadcast body to processors map
  MPI_Bcast(bodyToProc, nBodyToProc, MPI_INT, proc_data->masterRank, MPI_COMM_WORLD);
  
  if (!proc_data->is_master()) {
    for (int k = 0; k < nBodyToProc; k++) {
      cout << "bodyToProc[" << k << "]=" << bodyToProc[k] << endl;
    }
  }
  */
/*
  if (CUBIT_TRUE == DEBUG_FLAG(137)) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (proc_data->is_master()) {
      tScatter = MPI_Wtime() - tScatter;
      GeometryQueryTool *gti = GeometryQueryTool::instance();
      
      PRINT_DEBUG_137("total time to scatter (sec) = %f\n", tScatter);
      PRINT_DEBUG_137("Complexity: vertices, edges, faces, bodies = %d, %d, %d, %d\n",
		      gti->num_ref_vertices(), gti->num_ref_edges(), gti->num_ref_faces(), 
		      gti->num_bodies());
      PRINT_DEBUG_137("Stats: %s, %d, %d, %f, %f, %d, %d, %d, %d\n",
                      filename, ProcData::instance()->numProcs, m_sizeBuffer, 0.0, tScatter, 
                      gti->num_ref_vertices(), gti->num_ref_edges(), gti->num_ref_faces(), 
                      gti->num_bodies());
    }
  }
*/
  return status;
}
/*
CubitStatus CGMMemFile::scatter_entity_list(DLIList<RefEntity*> &ref_entity_list,
                                            const char *filename) 
{
  // scatter exact amount of geometry information to each processors
  CubitStatus status;
  double tScatter;
  int mySendCount;
  int nEntity;
  DLIList<RefEntity*> temp_list, temp_ref_list;

  ProcData *proc_data = ProcData::instance();
  if (DEBUG_FLAG(137) && proc_data->is_master())
    tScatter = MPI_Wtime();

  int nProcs = proc_data->numProcs;
  int *sendCounts = new int[nProcs];
  int *displacements = new int[nProcs];
  int nBarEntity;
  int restEntity;
  int nEndEntity;
  displacements[0] = 0;

  if (proc_data->is_master()) {
    int curPosition = 0;
    nEntity = ref_entity_list.size();
    nBarEntity = nEntity/nProcs;
    restEntity = nEntity%nProcs;
    nEndEntity = nBarEntity + restEntity;
    
    ref_entity_list.reset();
    int sum = 0;

    // make temporary lists to contain geometry information for each processors
    for (int i = 0; i < nProcs; i++) {
      for ( int j = 0; j < nBarEntity; j++) {
	RefEntity* body_ptr = ref_entity_list.get_and_step();
	temp_list.clean_out();
	temp_list.append(body_ptr);
	sendCounts[i] += get_refentity_list_size(temp_list);
      }
      
      if (i < restEntity) {
	RefEntity* body_ptr = ref_entity_list.get_and_step();
	temp_list.clean_out();
	temp_list.append(body_ptr);
	sendCounts[i] += get_refentity_list_size(temp_list);
      }

      sum += sendCounts[i];
    }

      // check the size of the buffer and resize if necessary
    check_size(sum);
    
      // now actually append the information
    ref_entity_list.reset();
    for (int i = 0; i < nProcs; i++) {
      for ( int j = 0; j < nBarEntity; j++) {
	temp_list.clean_out();
	temp_list.append(ref_entity_list.get_and_step());
	int count = get_refentity_list_size(temp_list);
	append_refentity_list(temp_list, count);
      }
      
      if (i < restEntity) {
	temp_list.clean_out();
	temp_list.append(ref_entity_list.get_and_step());
	int count = get_refentity_list_size(temp_list);
	append_refentity_list(temp_list, count);
      }

      //append_refentity_list(temp_list, sendCounts[i]);
      //temp_list.clean_out();
    }
  }

  // broadcast buffer size array
  PRINT_DEBUG_100("Broadcasting buffer size array from master.\n");
  MPI_Bcast(sendCounts, nProcs, MPI_INT, proc_data->masterRank, MPI_COMM_WORLD);
  
  for (int i = 1; i < nProcs; i++) {
    displacements[i] = displacements[i-1] + sendCounts[i-1];
  }
  
  mySendCount = sendCounts[proc_data->myRank];

  if (!proc_data->is_master()) check_size(mySendCount);

  PRINT_DEBUG_100("Scattering buffer from master.\n");

  // scatter geometry
  MPI_Scatterv(m_pBuffer, sendCounts, displacements, MPI_BYTE, m_pBuffer, 
	       mySendCount, MPI_BYTE, proc_data->masterRank, MPI_COMM_WORLD);

  if (!proc_data->is_master()) {
    status = read_refentity_list(ref_entity_list);
  }

  ref_entity_list.reset();
  for (int i = 0; i < ref_entity_list.size(); i++) {
    DLIList<RefEntity*> entities;
    RefEntity *entity = ref_entity_list.get_and_step();
    entity->get_all_child_ref_entities(entities);
    cout << "my_id=" << ProcData::instance()->myRank 
	 << " num of children=" << entities.size() << endl;
    entities.reset();
    for (int j = 0; j < entities.size(); j++) {
      cout << "my_id=" << ProcData::instance()->myRank 
	   << " measure=" << entities.get_and_step()->measure() << endl;
    }
  }

  if (CUBIT_TRUE == DEBUG_FLAG(137)) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (proc_data->is_master()) {
      tScatter = MPI_Wtime() - tScatter;
      GeometryQueryTool *gti = GeometryQueryTool::instance();
      
      PRINT_DEBUG_137("total time to scatter (sec) = %f\n", tScatter);
      PRINT_DEBUG_137("Complexity: vertices, edges, faces, bodies = %d, %d, %d, %d\n",
		      gti->num_ref_vertices(), gti->num_ref_edges(), gti->num_ref_faces(), 
		      gti->num_bodies());
      PRINT_DEBUG_137("Stats: %s, %d, %d, %f, %f, %d, %d, %d, %d\n",
                      filename, ProcData::instance()->numProcs, m_sizeBuffer, 0.0, tScatter, 
                      gti->num_ref_vertices(), gti->num_ref_edges(), gti->num_ref_faces(), 
                      gti->num_bodies());
    }
  }

  return status;
}
*/
/*
CubitStatus CGMMemFile::vol_balace_send_body_processors(DLIList<RefEntity*> &ref_entity_list) 
{
  ProcData *proc_data = ProcData::instance();
  bodyToProc.clean_out();
  int nEntity = 0;

  if (proc_data->is_master()) {
    int nProcs = proc_data->numProcs;
    nEntity = ref_entity_list.size();
    balancedLists = new DLIList<RefEntity*>*[nProcs];
    
    for (int k = 0; k < nProcs; k++)
      balancedLists[k] = new DLIList<RefEntity*>;
    
    totalLoads = new double[nProcs];
    ref_entity_list.sort(larger_vol);
    ref_entity_list.reset();
    
    // fill first bodies to processors
    for (int l = 0; l < nProcs; l++) {
      RefEntity *body_ptr = ref_entity_list.get_and_step();
      // just calculate volumes for load
      totalLoads[l] = body_ptr->measure();
      bodyToProc.append(TDUniqueId::get_unique_id(body_ptr));
      bodyToProc.append(l);
      balancedLists[l]->append(body_ptr);
    }
    
    // fill the rest
    for (int i = nProcs; i < nEntity; i++) {
      RefEntity *body_ptr = ref_entity_list.get_and_step();
      double min_vol = totalLoads[0];
      int min_i = 0;
      
      for (int j = 1; j < nProcs; j++) {
	if (min_vol > totalLoads[j]) {
	  min_vol = totalLoads[j];
	  min_i = j;
	}
      }
      
      totalLoads[min_i] += body_ptr->measure();
      bodyToProc.append(TDUniqueId::get_unique_id(body_ptr));
      bodyToProc.append(min_i);
      balancedLists[min_i]->append(body_ptr);
    }
  }

  // broadcast number of bodies
  MPI_Bcast(&nEntity, 1, MPI_INT, proc_data->masterRank, MPI_COMM_WORLD);
  int body_to_proc[2*nEntity];

  if (proc_data->is_master()) // make temporary array
    bodyToProc.copy_to(body_to_proc);

  MPI_Bcast(body_to_proc, 2*nEntity, MPI_INT, proc_data->masterRank, MPI_COMM_WORLD);

  if (!proc_data->is_master()) // copy from temporary array
    bodyToProc.copy_from(body_to_proc, 2*nEntity);

  return CUBIT_SUCCESS;
}

CubitStatus CGMMemFile::send_body_to_procs_balanced(DLIList<RefEntity*> &ref_entity_list) 
{
  ProcData *proc_data = ProcData::instance();
  int nProcs = proc_data->numProcs;
  int nEntity;

  // load balance bodies by volume
  if (proc_data->is_master()) {
    double *total_volumes = new double[nProcs];
    nEntity = ref_entity_list.size();
    ref_entity_list.sort(larger_vol);
    bodyToProc.clean_out();
    ref_entity_list.reset();

    // fill first bodies to processors
    for (int l = 0; l < nProcs; l++) {
      RefEntity *body_ptr = ref_entity_list.get_and_step();
      total_volumes[l] = body_ptr->measure();
      bodyToProc.append(TDUniqueId::get_unique_id(body_ptr));
      bodyToProc.append(l);
    }

    // fill the rest
    for (int i = nProcs; i < nEntity; i++) {
      RefEntity *body_ptr = ref_entity_list.get_and_step();
      double min_vol = total_volumes[0];
      int min_i = 0;

      for (int j = 1; j < nProcs; j++) {
	if (min_vol > total_volumes[j]) {
	  min_vol = total_volumes[j];
	  min_i = j;
	}
      }
      
      total_volumes[min_i] += body_ptr->measure();
      bodyToProc.append(TDUniqueId::get_unique_id(body_ptr));
      bodyToProc.append(min_i);
    }

    delete total_volumes;
  }

  // broadcast number of bodies
  MPI_Bcast(&nEntity, 1, MPI_INT, proc_data->masterRank, MPI_COMM_WORLD);
  int body_to_proc[2*nEntity];

  if (proc_data->is_master()) // make temporary array
    bodyToProc.copy_to(body_to_proc);

  MPI_Bcast(body_to_proc, 2*nEntity, MPI_INT, proc_data->masterRank, MPI_COMM_WORLD);

  if (!proc_data->is_master()) // copy from temporary array
    bodyToProc.copy_from(body_to_proc, 2*nEntity);

  return CUBIT_SUCCESS;
}
*/

/*
CubitStatus CGMMemFile::scatter_vol_balanced_entity_list(DLIList<RefEntity*> &ref_entity_list,
							 DLIList<RefEntity*> &ref_entity_list_master,
							 const char *filename) 
{
  // scatter exact amount of geometry information to each processors
  CubitStatus status;
  double tScatter;
  int mySendCount;
  int nEntity;

  ProcData *proc_data = ProcData::instance();
  if (DEBUG_FLAG(137) && proc_data->is_master())
    tScatter = MPI_Wtime();

  int nProcs = proc_data->numProcs;
  int *sendCounts = new int[nProcs];
  int *displacements = new int[nProcs];
  int nBarEntity;
  int restEntity;
  int nEndEntity;
  displacements[0] = 0;
  ref_entity_list_master;

  if (proc_data->is_master()) {
    int curPosition = 0;
    nEntity = ref_entity_list.size();
    nBarEntity = nEntity/nProcs;
    restEntity = nEntity%nProcs;
    nEndEntity = nBarEntity + restEntity;

    double *body_volumes = new double[nEntity];
    double *total_volumes = new double[nProcs];
    DLIList<RefEntity*> **load_balace_array = new DLIList<RefEntity*>*[nProcs];

    for (int k = 0; k < nProcs; k++)
      load_balace_array[k] = new DLIList<RefEntity*>;

    ref_entity_list.reset();
    for (int l = 0; l < nProcs; l++) {
      RefEntity *body_ptr = ref_entity_list.get_and_step();
      load_balace_array[l]->append(body_ptr);
      total_volumes[l] = body_ptr->measure();
    }

    for (int i = nProcs; i < nEntity; i++) {
      RefEntity *body_ptr = ref_entity_list.get_and_step();
      double min_vol = total_volumes[0];
      int min_i = 0;

      for (int j = 1; j < nProcs; j++) {
	if (min_vol > total_volumes[j]) {
	  min_vol = total_volumes[j];
	  min_i = j;
	}
      }
      
      load_balace_array[min_i]->append(body_ptr);
      total_volumes[min_i] += body_ptr->measure();
    }

    int sum = 0;
    // make temporary lists to contain geometry information for each processors
    for (int i = 0; i < nProcs; i++) {
      PRINT_DEBUG_100("load_balace_array[i]=%d\n", load_balace_array[i]);
      PRINT_DEBUG_100("*load_balace_array[i]=%d\n", *load_balace_array[i]);
      PRINT_DEBUG_100("sendCounts[i]=%d\n", sendCounts[i]);
      sendCounts[i] = get_refentity_list_size(*(load_balace_array[i]));
      sum += sendCounts[i];
    }

    // check the size of the buffer and resize if necessary
    check_size(sum);
    
    // now actually append the information
    ref_entity_list.clean_out();
    for (int i = 0; i < nProcs; i++) {
      append_refentity_list(*load_balace_array[i], sendCounts[i]);
      ref_entity_list += *load_balace_array[i]; // serialize ref_entity_list
    }
    
    // check and save ref_entity_list information for master
    ref_entity_list_master.clean_out();
    ref_entity_list_master += *load_balace_array[0];
    //ref_entity_list_master = *load_balace_array[0];
  }

  // broadcast buffer size array
  PRINT_DEBUG_100("Broadcasting buffer size array from master.\n");
  MPI_Bcast(sendCounts, nProcs, MPI_INT, proc_data->masterRank, MPI_COMM_WORLD);

  for (int i = 1; i < nProcs; i++) {
    displacements[i] = displacements[i-1] + sendCounts[i-1];
  }
  
  mySendCount = sendCounts[proc_data->myRank];

  if (!proc_data->is_master()) check_size(mySendCount);

  PRINT_DEBUG_100("Scattering buffer from master.\n");

  // scatter geometry
  MPI_Scatterv(m_pBuffer, sendCounts, displacements, MPI_BYTE, m_pBuffer, 
	       mySendCount, MPI_BYTE, proc_data->masterRank, MPI_COMM_WORLD);

  if (!proc_data->is_master())
    status = read_refentity_list(ref_entity_list);

  if (CUBIT_TRUE == DEBUG_FLAG(137)) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (proc_data->is_master()) {
      tScatter = MPI_Wtime() - tScatter;
      GeometryQueryTool *gti = GeometryQueryTool::instance();
      
      PRINT_DEBUG_137("total time to scatter (sec) = %f\n", tScatter);
      PRINT_DEBUG_137("Complexity: vertices, edges, faces, bodies = %d, %d, %d, %d\n",
		      gti->num_ref_vertices(), gti->num_ref_edges(), gti->num_ref_faces(), 
		      gti->num_bodies());
      PRINT_DEBUG_137("Stats: %s, %d, %d, %f, %f, %d, %d, %d, %d\n",
                      filename, ProcData::instance()->numProcs, m_sizeBuffer, 0.0, tScatter, 
                      gti->num_ref_vertices(), gti->num_ref_edges(), gti->num_ref_faces(), 
                      gti->num_bodies());
    }
  }

  return status;
}
*/
/*
// before partition input
CubitStatus CGMMemFile::scatter_balanced_entity_list(DLIList<RefEntity*> &ref_entity_list,
							 DLIList<RefEntity*> &ref_entity_list_master,
							 const char *filename,
							 DLIList<RefEntity*> **balanced_lists) 
{
  //cout << "my_id=" << ProcData::instance()->myRank << "CGMMemFile::scatter_balanced_entity_list1" << endl;
  // scatter exact amount of geometry information to each processors
  CubitStatus status;
  double tScatter;
  int mySendCount;
  ProcData *proc_data = ProcData::instance();

  if (DEBUG_FLAG(137) && proc_data->is_master())
    tScatter = MPI_Wtime();

  int nProcs = proc_data->numProcs;
  int *sendCounts = new int[nProcs];
  int *displacements = new int[nProcs];
  int nBarEntity;
  int restEntity;
  int nEndEntity;
  displacements[0] = 0;
  //ref_entity_list_master;

  //cout << "my_id=" << ProcData::instance()->myRank << "CGMMemFile::scatter_balanced_entity_list2" << endl;
  if (proc_data->is_master()) {
    //cout << "my_id=" << ProcData::instance()->myRank << "CGMMemFile::scatter_balanced_entity_list3" << endl;
    int sum = 0;

    // make temporary lists to contain geometry information for each processors
    for (int i = 0; i < nProcs; i++) {
      //cout << "CGMMemFile::scatter_balanced_entity_list3.11" << endl;
      sendCounts[i] = get_refentity_list_size(*(balanced_lists[i]));
      //cout << "sendCounts[" << i << "]=" << sendCounts[i] << ",balanced_lists[i]->size()="
      //<< balanced_lists[i]->size() << endl;
      sum += sendCounts[i];
    }

    //cout << "CGMMemFile::scatter_balanced_entity_list3.1" << endl;
    // check the size of the buffer and resize if necessary
    check_size(sum);
    //cout << "CGMMemFile::scatter_balanced_entity_list3.2" << endl;

    // now actually append the information
    ref_entity_list.clean_out();
    for (int i = 0; i < nProcs; i++) {
      append_refentity_list(*balanced_lists[i], sendCounts[i]);
      ref_entity_list += *balanced_lists[i]; // serialize ref_entity_list
    }
    
    //cout << "CGMMemFile::scatter_balanced_entity_list3.3" << endl;
    // check and save ref_entity_list information for master
    ref_entity_list_master.clean_out();
    ref_entity_list_master += *balanced_lists[0];
    //cout << "CGMMemFile::scatter_balanced_entity_list3.4" << endl;
  }

  //cout << "my_id=" << ProcData::instance()->myRank << "CGMMemFile::scatter_balanced_entity_list4" << endl;
  // broadcast buffer size array
  PRINT_DEBUG_100("Broadcast
ing buffer size array from master.\n");
  
  //for (int n = 0; n < nProcs; n++)
  //cout << "sendCounts[" << n << "]=" << sendCounts[n] << endl;

  MPI_Bcast(sendCounts, nProcs, MPI_INT, proc_data->masterRank, MPI_COMM_WORLD);
  //cout << "my_id=" << ProcData::instance()->myRank << "CGMMemFile::scatter_balanced_entity_list4.1" << endl;

  for (int i = 1; i < nProcs; i++) {
    displacements[i] = displacements[i-1] + sendCounts[i-1];
  }
  //cout << "my_id=" << ProcData::instance()->myRank << "CGMMemFile::scatter_balanced_entity_list4.2" << endl;
  
  mySendCount = sendCounts[proc_data->myRank];
  //cout << "my_id=" << ProcData::instance()->myRank << "CGMMemFile::scatter_balanced_entity_list4.3" << endl;

  if (!proc_data->is_master()) check_size(mySendCount);
  //cout << "my_id=" << ProcData::instance()->myRank << "CGMMemFile::scatter_balanced_entity_list4.4" << endl;

  PRINT_DEBUG_100("Scattering buffer from master.\n");

  //cout << "CGMMemFile::scatter_balanced_entity_list5" << endl;
  // scatter geometry
  MPI_Scatterv(m_pBuffer, sendCounts, displacements, MPI_BYTE, m_pBuffer, 
	       mySendCount, MPI_BYTE, proc_data->masterRank, MPI_COMM_WORLD);

  if (!proc_data->is_master())
    status = read_refentity_list(ref_entity_list);

  //cout << "CGMMemFile::scatter_balanced_entity_list6" << endl;
  if (CUBIT_TRUE == DEBUG_FLAG(137)) {
    //MPI_Barrier(MPI_COMM_WORLD);
    if (proc_data->is_master()) {
      tScatter = MPI_Wtime() - tScatter;
      GeometryQueryTool *gti = GeometryQueryTool::instance();
      
      PRINT_DEBUG_137("total time to scatter (sec) = %f\n", tScatter);
      PRINT_DEBUG_137("Complexity: vertices, edges, faces, bodies = %d, %d, %d, %d\n",
		      gti->num_ref_vertices(), gti->num_ref_edges(), gti->num_ref_faces(), 
		      gti->num_bodies());
      PRINT_DEBUG_137("Stats: %s, %d, %d, %f, %f, %d, %d, %d, %d\n",
                      filename, ProcData::instance()->numProcs, m_sizeBuffer, 0.0, tScatter, 
                      gti->num_ref_vertices(), gti->num_ref_edges(), gti->num_ref_faces(), 
                      gti->num_bodies());
    }
  }
  //cout << "CGMMemFile::scatter_balanced_entity_list7" << endl;

  return status;
}
*/
/*
// partition input without balancedLists
CubitStatus CGMMemFile::scatter_balanced_entity_list(DLIList<RefEntity*> &ref_entity_list,
						     DLIList<RefEntity*> &ref_entity_list_master,
						     const char *filename)
{
  //cout << "CGMMemFile::scatter_balanced_entity_list1"
  //   << ",balanced_lists=" << balanced_lists << endl;
  // scatter exact amount of geometry information to each processors
  CubitStatus status;
  double tScatter;
  int mySendCount;
  ProcData *proc_data = ProcData::instance();

  if (DEBUG_FLAG(137) && proc_data->is_master())
    tScatter = MPI_Wtime();

  int nProcs = proc_data->numProcs;
  int *sendCounts = new int[nProcs];
  int *displacements = new int[nProcs];
  //int nBarEntity;
  //int restEntity;
  //int nEndEntity;
  displacements[0] = 0;

  cout << "CGMMemFile::scatter_balanced_entity_list2" << endl;

  if (proc_data->is_master()) {
    cout << "CGMMemFile::scatter_balanced_entity_list3.1" << endl;

    DLIList<RefEntity*> **balance_lists = new DLIList<RefEntity*>*[nProcs];
    
    for (int i = 0; i < nProcs; i++)
      balance_lists[i] = new DLIList<RefEntity*>;

    cout << "CGMMemFile::scatter_balanced_entity_list3.2" << endl;

    int n_refs = ref_entity_list.size();
    ref_entity_list.reset();

    for (; n_refs > 0; n_refs--) {
      cout << "CGMMemFile::scatter_balanced_entity_list3.3" << endl;
      RefEntity *entity = ref_entity_list.get_and_step();
      TDParallel *td_par = (TDParallel *) entity->get_TD(&TDParallel::is_parallel);

      if (td_par == NULL) {
	td_par = new TDParallel(entity);
      }

      //cout << "td_par->get_proc_id()=" << td_par->get_proc_id() << endl;
      //cout << "balance_lists=" << balance_lists << endl;
      //cout << "balance_lists[td_par->get_proc_id()]="
      //   << balance_lists[td_par->get_proc_id()] << endl;

      int temp_id = td_par->get_proc_id();

      if (temp_id > -1)
	balance_lists[temp_id]->append(entity);
    }

    cout << "CGMMemFile::scatter_balanced_entity_list3.4" << endl;
    int sum = 0;
    for (int j = 0; j < nProcs; j++) {
      cout << "CGMMemFile::scatter_balanced_entity_list3.41" << endl;

      sendCounts[j] = get_refentity_list_size(*(balance_lists[j]));
      cout << "CGMMemFile::scatter_balanced_entity_list3.42" << endl;
      sum += sendCounts[j];
      cout << "CGMMemFile::scatter_balanced_entity_list3.43" << endl;
    }

    cout << "CGMMemFile::scatter_balanced_entity_list3.5" << endl;
    // check the size of the buffer and resize if necessary
    check_size(sum);

    cout << "CGMMemFile::scatter_balanced_entity_list3.6" << endl;
    // now actually append the information
    ref_entity_list.clean_out();
    for (int k = 0; k < nProcs; k++) {
      append_refentity_list(*balance_lists[k], sendCounts[k]);
      ref_entity_list += *balance_lists[k]; // serialize ref_entity_list
    }
    cout << "CGMMemFile::scatter_balanced_entity_list3.7" << endl;    
    // check and save ref_entity_list information for master
    ref_entity_list_master.clean_out();
    ref_entity_list_master += *balance_lists[0];

    for (int l = 0; l < nProcs; l++)
      delete balance_lists[l];
    cout << "CGMMemFile::scatter_balanced_entity_list3.8" << endl;
    delete [] balance_lists;

  }

  cout << "my_id=" << ProcData::instance()->myRank << "CGMMemFile::scatter_balanced_entity_list4" << endl;
  // broadcast buffer size array
  PRINT_DEBUG_100("Broadcast
ing buffer size array from master.\n");
  
  //for (int n = 0; n < nProcs; n++)
  //cout << "sendCounts[" << n << "]=" << sendCounts[n] << endl;

  MPI_Bcast(sendCounts, nProcs, MPI_INT, proc_data->masterRank, MPI_COMM_WORLD);
  cout << "my_id=" << ProcData::instance()->myRank << "CGMMemFile::scatter_balanced_entity_list4.1" << endl;

  for (int i = 1; i < nProcs; i++) {
    displacements[i] = displacements[i-1] + sendCounts[i-1];
  }
  cout << "my_id=" << ProcData::instance()->myRank << "CGMMemFile::scatter_balanced_entity_list4.2" << endl;
  
  mySendCount = sendCounts[proc_data->myRank];
  cout << "my_id=" << ProcData::instance()->myRank << "CGMMemFile::scatter_balanced_entity_list4.3" << endl;

  if (!proc_data->is_master()) check_size(mySendCount);
  cout << "my_id=" << ProcData::instance()->myRank << "CGMMemFile::scatter_balanced_entity_list4.4" << endl;

  PRINT_DEBUG_100("Scattering buffer from master.\n");

  cout << "CGMMemFile::scatter_balanced_entity_list5" << endl;
  // scatter geometry
  MPI_Scatterv(m_pBuffer, sendCounts, displacements, MPI_BYTE, m_pBuffer, 
	       mySendCount, MPI_BYTE, proc_data->masterRank, MPI_COMM_WORLD);

  cout << "CGMMemFile::scatter_balanced_entity_list5.1" << endl;
  if (!proc_data->is_master())
    status = read_refentity_list(ref_entity_list);

  cout << "CGMMemFile::scatter_balanced_entity_list6" << endl;
  if (CUBIT_TRUE == DEBUG_FLAG(137)) {
    //MPI_Barrier(MPI_COMM_WORLD);
    if (proc_data->is_master()) {
      tScatter = MPI_Wtime() - tScatter;
      GeometryQueryTool *gti = GeometryQueryTool::instance();
      
      PRINT_DEBUG_137("total time to scatter (sec) = %f\n", tScatter);
      PRINT_DEBUG_137("Complexity: vertices, edges, faces, bodies = %d, %d, %d, %d\n",
		      gti->num_ref_vertices(), gti->num_ref_edges(), gti->num_ref_faces(), 
		      gti->num_bodies());
      PRINT_DEBUG_137("Stats: %s, %d, %d, %f, %f, %d, %d, %d, %d\n",
                      filename, ProcData::instance()->numProcs, m_sizeBuffer, 0.0, tScatter, 
                      gti->num_ref_vertices(), gti->num_ref_edges(), gti->num_ref_faces(), 
                      gti->num_bodies());
    }
  }
  cout << "CGMMemFile::scatter_balanced_entity_list7" << endl;

  //for (int l = 0; l < nProcs; l++)
  //delete balance_lists[l];



  return status;
}
*/
// partition input
CubitStatus CGMMemFile::scatter_balanced_entity_list(DLIList<RefEntity*> &ref_entity_list,
                                                     DLIList<RefEntity*> &ref_entity_list_master,
                                                     DLIList<RefEntity*> **balanced_lists) 
{
  //cout << "CGMMemFile::scatter_balanced_entity_list1"
  //   << ",balanced_lists=" << balanced_lists << endl;
  // scatter exact amount of geometry information to each processors
  CubitStatus status;
  double tScatter;
  int mySendCount;
  ProcData *proc_data = ProcData::instance();

  if (DEBUG_FLAG(137) && proc_data->is_master())
#ifdef USE_MPI
    tScatter = MPI_Wtime();
#endif

  int nProcs = proc_data->numProcs;
  int *sendCounts = new int[nProcs];
  int *displacements = new int[nProcs];
  //int nBarEntity;
  //int restEntity;
  //int nEndEntity;
  displacements[0] = 0;

  //cout << "CGMMemFile::scatter_balanced_entity_list2" << endl;

  if (proc_data->is_master()) {
    //cout << "CGMMemFile::scatter_balanced_entity_list3.1" << endl;

    //DLIList<RefEntity*> **balance_lists = new DLIList<RefEntity*>*[nProcs];
    
    //for (int i = 0; i < nProcs; i++)
    //balance_lists[i] = new DLIList<RefEntity*>;

    //cout << "CGMMemFile::scatter_balanced_entity_list3.4" << endl;
    int sum = 0;
    for (int j = 0; j < nProcs; j++) {
      //cout << "CGMMemFile::scatter_balanced_entity_list3.41" << endl;
      //cout << "balanced_lists[j]->size()=" << balanced_lists[j]->size()
      //   << endl;

      sendCounts[j] = get_refentity_list_size(*(balanced_lists[j]));
      //cout << "CGMMemFile::scatter_balanced_entity_list3.42" << endl;
      sum += sendCounts[j];
      //cout << "CGMMemFile::scatter_balanced_entity_list3.43" << endl;
    }

    //cout << "CGMMemFile::scatter_balanced_entity_list3.5" << endl;
    // check the size of the buffer and resize if necessary
    check_size(sum);

    //cout << "CGMMemFile::scatter_balanced_entity_list3.6" << endl;
    // now actually append the information
    ref_entity_list.clean_out();
    for (int k = 0; k < nProcs; k++) {
      append_refentity_list(*balanced_lists[k], sendCounts[k]);
      ref_entity_list += *balanced_lists[k]; // serialize ref_entity_list
    }
    //cout << "CGMMemFile::scatter_balanced_entity_list3.7" << endl;    
    // check and save ref_entity_list information for master
    ref_entity_list_master.clean_out();
    //cout << "CGMMemFile::scatter_balanced_entity_list3.71" << endl;
    //cout << "balanced_lists[0]=" << balanced_lists[0] << endl;
    ref_entity_list_master += *balanced_lists[0];
    //cout << "CGMMemFile::scatter_balanced_entity_list3.72" << endl;
  }

  //cout << "my_id=" << ProcData::instance()->myRank << "CGMMemFile::scatter_balanced_entity_list4" << endl;
  // broadcast buffer size array
  PRINT_DEBUG_100("Broadcasting buffer size array from master.\n");
  
  //for (int n = 0; n < nProcs; n++)
  //cout << "sendCounts[" << n << "]=" << sendCounts[n] << endl;

#ifdef USE_MPI
  MPI_Bcast(sendCounts, nProcs, MPI_INT, proc_data->masterRank, MPI_COMM_WORLD);
#endif
  //cout << "my_id=" << ProcData::instance()->myRank << "CGMMemFile::scatter_balanced_entity_list4.1" << endl;

  for (int i = 1; i < nProcs; i++) {
    displacements[i] = displacements[i-1] + sendCounts[i-1];
  }
  //cout << "my_id=" << ProcData::instance()->myRank << "CGMMemFile::scatter_balanced_entity_list4.2" << endl;
  
  mySendCount = sendCounts[proc_data->myRank];
  //cout << "my_id=" << ProcData::instance()->myRank << "CGMMemFile::scatter_balanced_entity_list4.3" << endl;

  if (!proc_data->is_master()) check_size(mySendCount);
  //cout << "my_id=" << ProcData::instance()->myRank << "CGMMemFile::scatter_balanced_entity_list4.4" << endl;

  PRINT_DEBUG_100("Scattering buffer from master.\n");

  //cout << "CGMMemFile::scatter_balanced_entity_list5" << endl;
  // scatter geometry
#ifdef USE_MPI
  MPI_Scatterv(m_pBuffer, sendCounts, displacements, MPI_BYTE, m_pBuffer, 
	       mySendCount, MPI_BYTE, proc_data->masterRank, MPI_COMM_WORLD);
#endif

  //cout << "CGMMemFile::scatter_balanced_entity_list5.1" << endl;
  if (!proc_data->is_master())
    status = read_refentity_list(ref_entity_list);

  //cout << "CGMMemFile::scatter_balanced_entity_list6" << endl;
//  if (CUBIT_TRUE == DEBUG_FLAG(138)) {
    //MPI_Barrier(MPI_COMM_WORLD);
//    if (proc_data->is_master()) {
//      tScatter = MPI_Wtime() - tScatter;
//      GeometryQueryTool *gti = GeometryQueryTool::instance();
      
//      PRINT_DEBUG_138("total time to scatter (sec) = %f\n", tScatter);
//      PRINT_DEBUG_138("Complexity: vertices, edges, faces, bodies = %d, %d, %d, %d\n",
//		      gti->num_ref_vertices(), gti->num_ref_edges(), gti->num_ref_faces(), 
//		      gti->num_bodies());
//      PRINT_DEBUG_138("Stats: %s, %d, %d, %f, %f, %d, %d, %d, %d\n",
//                      filename, ProcData::instance()->numProcs, m_sizeBuffer, 0.0, tScatter, 
//                      gti->num_ref_vertices(), gti->num_ref_edges(), gti->num_ref_faces(), 
//                      gti->num_bodies());
//    }
//  }
  //cout << "CGMMemFile::scatter_balanced_entity_list7" << endl;

  //for (int l = 0; l < nProcs; l++)
  //delete balance_lists[l];



  return status;
}

CubitStatus CGMMemFile::check_size(int &target_size, CubitBoolean keep) 
{
  PRINT_DEBUG_100("Checking buffer size on proc %d, target size %d.\n", 
                  ProcData::instance()->myRank, target_size);
  if ((int)m_sizeBuffer < target_size) {
    PRINT_DEBUG_100("Increasing buffer size on proc %d.\n", ProcData::instance()->myRank);
    void *temp_buffer = malloc(target_size);
    if (keep && 0 != m_currentPosition) memcpy(temp_buffer, m_pBuffer, m_currentPosition);
    delete m_pBuffer;
    m_pBuffer = (unsigned char *) temp_buffer;
    m_sizeBuffer = target_size;
  }

  return CUBIT_SUCCESS;
}
DLIList <int> *CGMMemFile::get_body_to_proc() {
  return &bodyToProc;
}
/*
std::vector<int> *CGMMemFile::get_body_to_proc() {
  return &bodyToProc;
}
*/
/*
int *CGMMemFile::get_body_to_proc() {
  return bodyToProc;
}

int CGMMemFile::get_num_body_to_proc() {
  return nBodyToProc;
}
*/
