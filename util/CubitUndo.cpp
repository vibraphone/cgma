//- Class:       UndoCommands
//- Description: Provides main undo capability for Cubit system
//- Owner:       Ray J. Meyers
//- Checked by:
//- Version: $Id: 

#include "CubitUndo.hpp"
#include "CubitString.hpp"

CubitUndo* CubitUndo::mInstance = 0;

CubitUndo::CubitUndo()
{
  if(mInstance)
  {
    assert(0);
    // if you want a new instance in place of a previous one,
    // delete the previous one first
  }
  mInstance = this;
}

CubitUndo::~CubitUndo()
{
  mInstance = 0;
}

void CubitUndo::start_undo_group()
{
  if(!mInstance ) return;
  mInstance->p_start_undo_group(); 
}

void CubitUndo::end_undo_group()
{
  if(!mInstance ) return;
  mInstance->p_end_undo_group(); 
}

void CubitUndo::set_undo_enabled(bool value)
{
  if(!mInstance ) return;
  mInstance->p_set_undo_enabled( value );
}

bool CubitUndo::get_undo_enabled()
{
  if(!mInstance ) return false;
  return mInstance->p_get_undo_enabled();
}

void CubitUndo::clear_undo_groups()
{
  if(!mInstance ) return;
  mInstance->p_clear_undo_groups();
}

int CubitUndo::execute_undo()
{
  if(!mInstance ) 
    return 0;

  return mInstance->p_execute_undo();
}

void CubitUndo::remove_last_undo()
{
  if(!mInstance ) return;
  mInstance->p_remove_last_undo();
}


void CubitUndo::remove_last_undo_group()
{
  if(!mInstance ) return;
  mInstance->p_remove_last_undo_group();
}

void CubitUndo::save_state()
{
  if(!mInstance ) return;
  mInstance->p_save_state();
}

void CubitUndo::save_state_with_cubit_file(DLIList<RefVolume*> &volumes_to_save )
{
  if(!mInstance ) return;
  mInstance->p_save_state_with_cubit_file( volumes_to_save );
}

void CubitUndo::save_state_with_cubit_file(RefFace *face_to_save )
{
  if(!mInstance ) return;
  mInstance->p_save_state_with_cubit_file( face_to_save );
}

void CubitUndo::save_state_with_cubit_file(RefVolume *volume_to_save )
{
  if(!mInstance ) return;
  mInstance->p_save_state_with_cubit_file( volume_to_save );
}

void CubitUndo::save_state_with_cubit_file( DLIList<MRefVolume*> &volumes_to_save )
{
  if(!mInstance ) return;
  mInstance->p_save_state_with_cubit_file( volumes_to_save );
}
void CubitUndo::save_state_with_cubit_file(DLIList<RefFace*> &faces_to_save )
{
  if(!mInstance ) return;
  mInstance->p_save_state_with_cubit_file( faces_to_save );
}

void CubitUndo::save_state_with_cubit_file(DLIList<RefEdge*> &edges_to_save )
{
  if(!mInstance ) return;
  mInstance->p_save_state_with_cubit_file( edges_to_save );
}

void CubitUndo::save_state_with_cubit_file(DLIList<RefVertex*> &verts_to_save,
                                           bool save_only_if_free )
{
  if(!mInstance ) return;
  mInstance->p_save_state_with_cubit_file( verts_to_save, save_only_if_free );
}

void CubitUndo::save_state_with_cubit_file( DLIList<Body*> &bodies_to_save,
                                            DLIList<RefEntity*> *free_ents_to_save )
{
  if(!mInstance ) return;
  mInstance->p_save_state_with_cubit_file( bodies_to_save, free_ents_to_save );
}

void CubitUndo::note_result_body( Body *body_to_delete )
{
  if(!mInstance ) return;
  mInstance->p_note_result_body( body_to_delete );
}

void CubitUndo::note_result_bodies( DLIList<Body*> &bodies_to_delete )
{
  if(!mInstance ) return;
  mInstance->p_note_result_bodies( bodies_to_delete );
}
void CubitUndo::note_result_boundary_conditions( RefEntity *owner_entity )
{
    if(!mInstance ) return;
  mInstance->p_note_result_boundary_conditions( owner_entity );
}

void CubitUndo::note_result_entity( RefEntity *entity_to_delete ) 
{
  if(!mInstance ) return;
  mInstance->p_note_result_entity( entity_to_delete );
}

void CubitUndo::note_result_entities( DLIList<MRefVolume*> &vols_to_delete )
{
  if(!mInstance ) return;
  mInstance->p_note_result_entities( vols_to_delete );
}

void CubitUndo::note_result_entities( DLIList<RefEntity*> &entities_to_delete )
{
  if(!mInstance ) return;
  mInstance->p_note_result_entities( entities_to_delete );
}

void CubitUndo::set_undo_by_restoring_state()
{
  if(!mInstance ) return;
  mInstance->p_set_undo_by_restoring_state();
}

void CubitUndo::set_undo_by_command(CubitString undo_command_string)
{
  if(!mInstance ) return;
  mInstance->p_set_undo_by_command( undo_command_string );
}

int CubitUndo::number_undo_groups()
{
  if(!mInstance ) return 0;
  return mInstance->p_number_undo_groups();
}

void CubitUndo::create_dummy_undo_object()
{
  if(!mInstance ) return;
  mInstance->p_create_dummy_undo_object();
}

bool CubitUndo::inside_undo_group()
{
  if(!mInstance ) return false;
  return mInstance->p_inside_undo_group();
}

void CubitUndo::save_mesh_size_data( MRefEntity *entity )
{
  if(!mInstance ) return;
  mInstance->p_save_mesh_size_data( entity );
}

void CubitUndo::save_mesh_scheme_data( MRefEntity *entity )
{
  if(!mInstance ) return;
  mInstance->p_save_mesh_scheme_data( entity );
}

void CubitUndo::save_smooth_scheme_data( MRefEntity *entity )
{
  if(!mInstance ) return;
  mInstance->p_save_smooth_scheme_data( entity );
}

void CubitUndo::save_all_mesh_data( MRefEntity *entity )
{
  if(!mInstance ) return;
  mInstance->p_save_all_mesh_data( entity );
}

void CubitUndo::save_merge_data( DLIList<RefEntity*> &ents, bool merging )
{
  if(!mInstance ) return;
  mInstance->p_save_merge_data( ents, merging );
}

void CubitUndo::save_mesh_state_before_smoothing( MRefEntity *entity )
{
  if(!mInstance ) return;
  DLIList<MRefEntity*> ents_to_smooth(1);
  ents_to_smooth.append( entity );
  mInstance->p_save_mesh_state_before_smoothing( ents_to_smooth ); 
}

void CubitUndo::save_mesh_state_before_smoothing( DLIList<MRefEntity*> &ents_to_smooth )
{
  if(!mInstance ) return;
  mInstance->p_save_mesh_state_before_smoothing( ents_to_smooth ); 
}

void CubitUndo::save_mesh_state_before_smoothing( DLIList<MeshEntity*> &ents_to_smooth )
{
  if(!mInstance ) return;
  mInstance->p_save_mesh_state_before_smoothing( ents_to_smooth ); 
}

