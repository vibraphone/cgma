#ifndef __CUBIT_EVENT_DEFINES__
#define __CUBIT_EVENT_DEFINES__

/* Event types. These are used in the Subject-Observer patterns in the
 * code to specify what event occurred.  This allows the Observers to
 *          update themselves accordingly, if required. */
enum CubitEventType { INVALID_EVENT_TYPE = -1,
                 VGI_BODY_CONSTRUCTED,
                   /* A ModelEntity was created */
                 MODEL_ENTITY_CONSTRUCTED,
                   /* A ModelEntity was created */
                 MODEL_ENTITY_MODIFIED,
                   /* A ModelEntity was changed */
                 MODEL_ENTITY_DESTRUCTED,
                   /* A ModelEntity was deleted */
                 COMPARISON_FOUND,
                   /* A partner entity was found (they compare successfully) */
                 DEVELOPER_COMMAND_FLAG_MODIFIED,
                  /* The 'set developer command on/off' was issued */
                 ENTITY_SURVIVED_MERGE,
                   /* An entity was involved in a merge and survived */
                 MERGE_COMPLETED,
                   /* A merge operation is completed successfully */
                 MERGE_ABORTED,
                   /* A merge operation was aborted */
                 MESH_DELETED,
                   /* Mesh associated with 1 or more RefEntities was deleted */
                 MESH_CREATED,
                   /* 1 or more RefEntities was meshed */
                 MESH_MODIFIED,
                   /* 1 or more RefEntities had its exiting mesh modified */
                 MESH_TRANSFORMED,
                   /* Mesh nodes moved*/
                 GEOMETRY_TOPOLOGY_MODIFIED,
                   /* Both geometry and topology was
                    modified. e.g. partitioned.
                   */
                 TOPOLOGY_MODIFIED,
                   /* The topology of a RefVolume or one of its sub-entities
                    * was changed.  Used for virtual topology. */
                 GEOMETRY_MODIFIED,
                   /* The geometry of a RefEntity was altered. */
                 HEALER_COMPLETED,
                   /* Notifies that the healer completed */
                 COMPOSITE_CREATION_COMPLETED,
                   /* Notifies that a composite was created */
                 SPLIT_SURFACE_COMPLETED,
                   /* Notifies that a split face operation is done */
                 SEPARATE_OPERATION_COMPLETED,
                   /* Notifies that a split face operation is done */
                 COLLAPSE_CURVE_COMPLETED,
                   /* Notifies that a collapse curve operation is done */
                 AUTO_CLEAN_COMPLETED,
                   /* Notifies that an auto_clean operation is done */
                 REMOVE_SURFACE_COMPLETED,
                   /* Notifies that a remove surface operation is done */
                 REMOVE_TOPOLOGY_COMPLETED,
                   /* Notifies that a remove_topology operation is done */
                 REGULARIZE_ENTITY_COMPLETED,
                   /* Notifies that a regularize entity operation is done */
                 NEW_ENTITY_UNMERGED,
                   /*A surface, curve, or vertex was unmerged. */
                 DAG_NODE_DESTRUCTED,
                   /* A DAG Node was destructed
                    * CHILD_TO_BE_SWAPPED
                    * A child is about to be swapped */
                 FREE_REF_ENTITY_GENERATED,
                   /* A vertex that is not a part of a curve,
                    * or a curve that is not a part of a surface,
                    * was just generated */
                 ENTITY_CONSTRUCTED,
                 ENTITY_DESTRUCTED,
                   /* An generic entity is constructed or destructed; used to
                    * notify observers of destruction of observables */
                 TOP_LEVEL_ENTITY_DESTRUCTED,
                   /* A body or free entity and all it's children
                    * are about to be destructed.
                    * The pointer is still valid at this point for all
                    * calls except geometry queries. 
                    * The pointer can be cast to a higher level entity
                    * if desired. */

                 ENTITY_NAME_CHANGED,
                   /* the name for an entity changed */
                 
                 ENTITY_GEOMETRY_COLOR_CHANGED,
                 ENTITY_MESH_COLOR_CHANGED,
                   /* the color for an entity changed */
                 ENTITY_VISIBILITY_CHANGED,
                  /* the visibility of an entity changed */

                 ENTITIES_MERGED,
                   /* two entities are merged together. See MergeEvent. */
                 ID_SET,
                   /* an entity had its id set (changed). See IdSetEvent. */
                 IDS_COMPRESSED,
                   /* ids have been compressed. See model.cpp */
                 MODEL_ENTITY_HIDDEN,
                   /* A sub-type of MODEL_ENTITY_DESTRUCTED to indicate
                    * that the entity in question is being hidden (VG)
                    * rather than actually being destructed.           */
                 MODEL_ENTITY_RESTORED,
                   /* The reverse of MODEL_ENTITY_HIDDEN, and a sub-type
                      of MODEL_ENTITY_CONSTRUCTED */
                 MODEL_RESET,
                   /* Indicates a 'reset' command has been issued and commpleted. */
                 GROUP_MODIFIED,
                   // Group was modified.  Sent only once after all
                   // modifications are complete.
                 GRAPHICS_MODE_MODIFIED,
                   // The overall graphics mode has been changed in some way.
                   // This could result from the perspective being changed, the
                   // rendering mode changed, the pick mode changed, etc.
                 
                 GENESIS_ENTITY_CREATED,
                   // Block, nodeset or sideset was created.
                 GENESIS_ENTITY_DELETED,
                   // Block, nodeset or sideset was deleted.
                 GENESIS_ENTITY_MODIFIED,
                   // Block, nodeset or sideset was modified.
                 SUSPEND_GENESIS_PROCESSING,
                   // Suspend processing (in CubitInterface) for Genesis events
                 RESUME_GENESIS_PROCESSING,
                   // Resume processing (in CubitInterface) for Genesis events
                 UPDATE_GENESIS_DISPLAY,
                   // Force an update of all genesis entity display
                 SUSPEND_TREE_PROCESSING,
                   // suspend processing of events sent to the GUI tree
                 RESUME_TREE_PROCESSING,
                   // resume processing of events sent to the GUI tree
                 MATERIAL_CREATED,
                   // Material created
                 MATERIAL_MODIFIED,
                   // Material modified
                 MATERIAL_DELETED,
                   // Material deleted
                 MEDIA_CREATED,
                   // Media created
                 MEDIA_MODIFIED,
                   // Media modified
                 MEDIA_DELETED,
                   // Media deleted
                 BC_ENTITY_CREATED,
                   // BCS
                 BC_ENTITY_DELETED,
                   // BCS
                 BC_ENTITY_MODIFIED,
                   // BCS
                 IC_ENTITY_CREATED,
                   // ICS
                 IC_ENTITY_DELETED,
                   // ICS
                 IC_ENTITY_MODIFIED,
                   // ICS
                 BC_CONTAINER_ENTITY_CREATED,
                   // BC Set, Contact Pair
                 BC_CONTAINER_ENTITY_DELETED,
                 // BC Set, Contact Pair
                 BC_CONTAINER_ENTITY_MODIFIED,
                 // BC Set, Contact Pair
                 CONSTRAINT_ENTITY_CREATED,
                 // Constraints
                 CONSTRAINT_ENTITY_DELETED,
                 // Constraints
                 CONSTRAINT_ENTITY_MODIFIED,
                 // Constraints
                 COORDINATE_SYSTEM_CREATED,
                   // Coordinate system
                 COORDINATE_SYSTEM_DELETED,
                   // Coordinate system
                 COORDINATE_SYSTEM_MODIFIED,
                   // Coordinate system
                 UNDO_STATE_CHANGED,
                   // Undo has been changed in some way
                 UNDO_COMPLETE,
                   // Undo has completed its processing

                   // ******** Assembly Events  *********
                   //   All assembly events pass an AssemblyEvent
                   //   object.  The functions used to get valid
                   //   data for a particular event type is listed
                   //   right after the event type name.
                 WEBCUT_COMPLETED,
                   // Inform regarding the completion of a webcut operation
                 
                 ASSEMBLY_ADD_CHILD,
                   // (get_assembly(), get_node())
                   // Sent just after a new child is added to an assembly
                 
                 ASSEMBLY_REMOVE_CHILD,
                   // (get_assembly(), get_node())
                   // Sent just after a child is removed from an assembly
                 
                 ASSEMBLY_NODE_NAME_CHANGE,
                   // (get_node())
                   // Sent just after an assembly's name is changed. Note
                   // that the instance number may change along with the
                   // name, but this won't trigger a separate
                   // ASSEMBLY_NODE_INSTANCE_CHANGE event.
                 
                 ASSEMBLY_NODE_INSTANCE_CHANGE,
                   // (get_node())
                   // Sent just after an assembly or part's instance
                   // number is changed
                 
                 ASSEMBLY_NODE_STRING_PROPERTY_CHANGE,
                   // (get_node(), get_property())
                   // Sent just after an assembly or node has one of its
                   // properties changed, and that property is a string
                   // value
                 
                 PART_ADD_VOLUME,
                   // (get_part(), get_volume())
                   // Sent just after a volume is associated with a part
                 
                 PART_REMOVE_VOLUME,
                   // (get_part(), get_volume())
                   // Sent just after a volume is removed from a part
                 
                 ASSEMBLY_TREE_STRING_PROPERTY_CHANGE,
                   // (get_property())
                   // Sent just after a property of the tree changes.
                   // These are properties that apply to the tree as a whole,
                   // not to any individual AssemblyNode.

                   // 
                   // ********* End of Assembly Events **********

                 ACIS_FILE_IMPORTED,
                   // Import an ACIS file
                   
                 GRANITE_FILE_IMPORTED,
                   // Import a granite file

                 GEOMETRY_ENGINE_CHANGED,
                   // The geometry engine changed

                 WORKING_DIRECTORY_CHANGED,
                   // The working directory has changed

                 APREPRO_MODIFIED,
                   // The working directory has changed

                 FATAL_ERROR_ENCOUNTERED
                   // If this event occurs, you're done

};

#endif // __CUBIT_EVENT_DEFINES__
