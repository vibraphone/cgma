//- Class:       UndoCommands
//- Description: Provides undo capability for Cubit system
//- Owner:       Ray J. Meyers
//- Checked by:
//- Version: $Id: 


#ifndef CUBITUNDO_HPP
#define CUBITUNDO_HPP

#include <stack>
#include "DLIList.hpp"

class UndoCommand;
class CubitString;
class RefEntity;
class Body;
class RefVolume;
class RefFace;
class RefEdge;
class RefVertex;
class MRefEntity;
class MRefVolume;
class MeshEntity;

class CUBIT_UTIL_EXPORT CubitUndo 
{
  // the one instance of the Cubit Undo object
  // creating a second instance will delete and replace the first instance.
  static CubitUndo* mInstance;

public:
  
  CubitUndo();
  virtual ~CubitUndo();

  static void start_undo_group();
  static void end_undo_group();
  static void set_undo_enabled(bool value);
  static bool get_undo_enabled();
  static void clear_undo_groups();
  static int  execute_undo();
  static void set_undo_by_restoring_state();
  static void set_undo_by_command(CubitString command);

  static void save_state_with_cubit_file( DLIList<Body*> &bodies_to_save,
                                          DLIList<RefEntity*> *free_ents_to_save = NULL );
  static void save_state_with_cubit_file( DLIList<RefVolume*> &volumes_to_save );
  static void save_state_with_cubit_file( DLIList<RefFace*> &faces_to_save );
  static void save_state_with_cubit_file( RefFace* face_to_save );
  static void save_state_with_cubit_file( RefVolume* volume_to_save );
  static void save_state_with_cubit_file( DLIList<MRefVolume*> &volumes_to_save );
  static void save_state_with_cubit_file( DLIList<RefEdge*> &edges_to_save );
  static void save_state_with_cubit_file( DLIList<RefVertex*> &vertices_to_save, 
                                          bool save_if_free = false );

  static void save_state(); 
  static  int number_undo_groups();
  static void note_result_bodies( DLIList<Body*> &bodies_to_delete );
  static void note_result_body( Body* body_to_delete );
  static void note_result_entities( DLIList<RefEntity*> &entities_to_delete );
  static void note_result_entities( DLIList<MRefVolume*> &vols_to_delete );
  static void note_result_entity( RefEntity *entity_to_delete );
  static void remove_last_undo();
  static void remove_last_undo_group();
  static void create_dummy_undo_object();
  static bool inside_undo_group();

  //meshing data
  static void save_mesh_size_data( MRefEntity *entity );
  static void save_mesh_scheme_data( MRefEntity *entity );
  static void save_smooth_scheme_data( MRefEntity *entity );
  static void save_all_mesh_data( MRefEntity *entity );

  //merging
  static void save_merge_data( DLIList<RefEntity*> &ents, bool merging = true );

  //smoothing
  static void save_mesh_state_before_smoothing( MRefEntity *entity );
  static void save_mesh_state_before_smoothing( DLIList<MRefEntity*> &ents );
  static void save_mesh_state_before_smoothing( DLIList<MeshEntity*> &ents );

protected:
  
  virtual void p_start_undo_group() = 0;
  virtual void p_end_undo_group() = 0;

  virtual void p_set_undo_enabled(bool value)=0;
  virtual bool p_get_undo_enabled()=0;
  virtual void p_clear_undo_groups()=0;
  virtual int p_execute_undo()=0;
  virtual void p_set_undo_by_restoring_state()=0;
  virtual void p_set_undo_by_command(CubitString command)=0;

  virtual void p_save_state_with_cubit_file( DLIList<Body*> &bodies_to_save, 
                                             DLIList<RefEntity*> *free_ref_ents = NULL )=0;
  virtual void p_save_state_with_cubit_file( DLIList<RefVolume*> &volumes_to_save )=0;
  virtual void p_save_state_with_cubit_file( DLIList<RefFace*> &faces_to_save )=0;
  virtual void p_save_state_with_cubit_file( RefFace *face_to_save )=0;
  virtual void p_save_state_with_cubit_file( RefVolume* volume_to_save )=0;
  virtual void p_save_state_with_cubit_file( DLIList<MRefVolume*> &volumes_to_save )=0;
  virtual void p_save_state_with_cubit_file( DLIList<RefEdge*> &edges_to_save )=0;
  virtual void p_save_state_with_cubit_file( DLIList<RefVertex*> &vertices_to_save,
                                             bool save_only_if_free = false )=0;

  virtual void p_save_state()=0; 
  virtual  int p_number_undo_groups()=0;
  virtual void p_note_result_bodies( DLIList<Body*> &bodies_to_delete )=0;
  virtual void p_note_result_body( Body* body_to_delete )=0;
  virtual void p_note_result_entities( DLIList<RefEntity*> &entities_to_delete )=0;
  virtual void p_note_result_entities( DLIList<MRefVolume*> &vols_to_delete )=0;
  virtual void p_note_result_entity( RefEntity *entity_to_delete )=0;
  virtual void p_remove_last_undo_group()=0;
  virtual void p_remove_last_undo()=0;
  virtual void p_create_dummy_undo_object()=0;
  virtual bool p_inside_undo_group()=0;

  //meshing data
  virtual void p_save_mesh_size_data( MRefEntity *entity )=0;
  virtual void p_save_mesh_scheme_data( MRefEntity *entity )=0;
  virtual void p_save_smooth_scheme_data( MRefEntity *entity )=0;
  virtual void p_save_all_mesh_data( MRefEntity *entity )=0;
  
  //merging
  virtual void p_save_merge_data( DLIList<RefEntity*> &ents, bool merging = true)=0;

  //smoothing
 virtual void p_save_mesh_state_before_smoothing( DLIList<MRefEntity*> &ents_to_smooth )=0;
 virtual void p_save_mesh_state_before_smoothing( DLIList<MeshEntity*> &ents_to_smooth )=0;

}; 

#endif




