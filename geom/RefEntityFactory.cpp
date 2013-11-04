#include "RefEntityFactory.hpp"
#include "RefVertex.hpp"
#include "RefEdge.hpp"
#include "RefFace.hpp"
#include "RefVolume.hpp"
#include "Body.hpp"
#include "RefGroup.hpp"
#include "GeometryQueryTool.hpp"
#include "DLIList.hpp"
#include "AppUtil.hpp"

static int sort_by_ascending_ids(CubitEntity*& a, CubitEntity*& b)
{
  if( a->id() < b->id() )
    return -1;
  else
    return 1;
} 

static int sort_by_ascending_ids(RefVertex*& a, RefVertex*& b)
{
  if( a->id() < b->id() ) 
    return -1;
  else
    return 1;
} 

static int sort_by_ascending_ids(RefVolume*& a, RefVolume*& b)
{
  if( a->id() < b->id() )
    return -1;
  else
    return 1;
} 

static int sort_by_ascending_ids(RefEdge*& a, RefEdge*& b)
{
  if( a->id() < b->id() )
    return -1;
  else
    return 1;
} 

static int sort_by_ascending_ids(RefFace*& a, RefFace*& b)
{
  if( a->id() < b->id() )
    return -1;
  else
    return 1;
} 

static int sort_by_ascending_ids(Body*& a, Body*& b)
{
  if( a->id() < b->id() )
    return -1;
  else
    return 1;
} 

static int sort_by_ascending_ids(RefGroup*& a, RefGroup*& b)
{
  if( a->id() < b->id() )
    return -1;
  else
    return 1;
} 

RefEntityFactory *RefEntityFactory::instance_ = NULL;

RefEntityFactory *RefEntityFactory::instance()
{
  if (instance_ == NULL) {
    instance_ = new RefEntityFactory();
  }

  return instance_;
}

RefEntityFactory::RefEntityFactory(RefEntityFactory *factory) 
{
  instance_ = factory;

  AppUtil::instance()->event_dispatcher().add_observer(this);


  if (factory == NULL) {
    refVertexList = new DLIList<RefVertex*>();
    refEdgeList = new DLIList<RefEdge*>();
    refFaceList = new DLIList<RefFace*>();
    refGroupList = new DLIList<RefGroup*>();
    refVolumeList = new DLIList<RefVolume*>();
    bodyList = new DLIList<Body*>();
    register_observer(this);
  }
  else {
    refVertexList = NULL;
    refEdgeList = NULL;
    refFaceList = NULL;
    refGroupList = NULL;
    refVolumeList = NULL;
    bodyList = NULL;
  }

  refVertexListIsSorted = true;
  refEdgeListIsSorted = true;
  refFaceListIsSorted = true;
  refVolumeListIsSorted = true;
  bodyListIsSorted = true;
  refGroupListIsSorted = true;
  ManageListSorting = true;
  
  reset_ids();
}

void RefEntityFactory::delete_instance()
{
  if( NULL != instance_ )
  {
    delete instance_;
    instance_ = NULL;
  }
}

RefEntityFactory::~RefEntityFactory()
{
  if (refVertexList != NULL) delete refVertexList;
  if (refEdgeList != NULL) delete refEdgeList;
  if (refFaceList != NULL) delete refFaceList;
  if (refGroupList != NULL) delete refGroupList;
  if (refVolumeList != NULL) delete refVolumeList;
  if (bodyList != NULL) delete bodyList;
  unregister_observer(this);
  instance_ = NULL;
}

RefVertex *RefEntityFactory::construct_RefVertex(TBPoint *point)
{
  RefVertex *temp = new RefVertex(point);
  AppUtil::instance()->send_event( temp, MODEL_ENTITY_CONSTRUCTED );
  return temp;
}

RefEdge *RefEntityFactory::construct_RefEdge(Curve *curve)
{
  RefEdge *temp = new RefEdge(curve);
  AppUtil::instance()->send_event( temp, MODEL_ENTITY_CONSTRUCTED );
  return temp;
}

RefFace *RefEntityFactory::construct_RefFace(Surface *surface)
{
  RefFace *temp = new RefFace(surface);
  AppUtil::instance()->send_event( temp, MODEL_ENTITY_CONSTRUCTED );
  return temp;
}

RefVolume *RefEntityFactory::construct_RefVolume(Lump *lump)
{
  RefVolume *temp = new RefVolume(lump);
  AppUtil::instance()->send_event( temp, MODEL_ENTITY_CONSTRUCTED );
  return temp;
}

Body *RefEntityFactory::construct_Body(BodySM *body_sm)
{
  Body *temp = new Body(body_sm);
  AppUtil::instance()->send_event( temp, MODEL_ENTITY_CONSTRUCTED );
  return temp;
}

RefGroup *RefEntityFactory::construct_RefGroup(const char* name)
{
  RefGroup *temp = new RefGroup(name);
  AppUtil::instance()->send_event( temp, MODEL_ENTITY_CONSTRUCTED );
  return temp;
}

RefGroup *RefEntityFactory::construct_RefGroup (DLIList<RefEntity*>& entity_list)
{
  RefGroup *temp = new RefGroup(entity_list);
  AppUtil::instance()->send_event( temp, MODEL_ENTITY_CONSTRUCTED );
  return temp;
}

CubitStatus RefEntityFactory::ref_entity_list( char const* keyword,
                                               DLIList<RefEntity*> &entity_list,
                                               const CubitBoolean print_errors )
{
   if ( !strcmp( keyword, "body" ) || !strcmp( keyword, "Body" )) {
     CAST_LIST_TO_PARENT((*bodyList), entity_list) ;
   }
   else if ( !strcmp( keyword, "curve" ) || !strcmp( keyword, "Curve" ) ) {
     CAST_LIST_TO_PARENT((*refEdgeList), entity_list) ;
   }
   else if ( !strcmp( keyword, "volume" ) || !strcmp( keyword, "Volume" ) ) {
     CAST_LIST_TO_PARENT((*refVolumeList), entity_list) ; 
   }
   else if ( !strcmp( keyword, "vertex" ) || !strcmp( keyword, "Vertex" ) ) {
     CAST_LIST_TO_PARENT((*refVertexList), entity_list) ;
   }
   else if ( !strcmp( keyword, "surface" ) || !strcmp( keyword, "Surface" ) ) {
     CAST_LIST_TO_PARENT((*refFaceList), entity_list) ; 
   }
   else if ( !strcmp( keyword, "group" ) || !strcmp( keyword, "Group" ) ) {
     CAST_LIST_TO_PARENT((*refGroupList), entity_list) ;
   }
   else 
   {
     if (print_errors)
       PRINT_ERROR("Invalid list type for the ref_entity_list "
                   "function: %s\n", keyword);
     return CUBIT_FAILURE;
   }

   return CUBIT_SUCCESS;
}

void RefEntityFactory::add(RefEntity* ref_entity)
{
  assert(ref_entity != 0);
  
  RefGroup *group;
  RefVolume *volume;
  RefFace *face;
  RefEdge *edge;
  RefVertex *vertex;
  Body *body;

  if ( (group = CAST_TO(ref_entity, RefGroup)) != NULL ) {
    add(group);
  }
  else if ( (volume = CAST_TO(ref_entity, RefVolume)) != NULL ) {
    add(volume);
  }
  else if ( (face = CAST_TO(ref_entity, RefFace)) != NULL ) {
    add(face);
  }
  else if ( (edge = CAST_TO(ref_entity, RefEdge)) != NULL ) {
    add(edge);
  }
  else if ( (vertex = CAST_TO(ref_entity, RefVertex)) != NULL ) {
    add(vertex);
  }
  else if ( (body = CAST_TO(ref_entity, Body)) != NULL ) {
    add(body);
  }
}

void RefEntityFactory::add(RefGroup* refGPtr)
{
   assert(!refGroupList->is_in_list(refGPtr));
   if(refGroupList->size() > 0)
   {
    refGroupList->last();
    if(ManageListSorting && refGPtr->id() < refGroupList->get()->id())
      refGroupListIsSorted = false;
   }
   refGroupList->append(refGPtr);
   if (refGPtr->id() > maxRefGroupId) maxRefGroupId = refGPtr->id();
}

void RefEntityFactory::add(Body* bodyPtr)
{
   assert(!bodyList->is_in_list(bodyPtr));
   if(bodyList->size() > 0)
   {
    bodyList->last();
    if(ManageListSorting && bodyPtr->id() < bodyList->get()->id())
      bodyListIsSorted = false;
   }
   bodyList->append(bodyPtr);
   if (bodyPtr->id() > maxBodyId) maxBodyId = bodyPtr->id();
}

void RefEntityFactory::add(RefVolume* refVPtr)
{
   assert(!refVolumeList->is_in_list(refVPtr));
   if(refVolumeList->size() > 0)
   {
    refVolumeList->last();
    if(ManageListSorting && refVPtr->id() < refVolumeList->get()->id())
      refVolumeListIsSorted = false;
   }
   refVolumeList->append(refVPtr);
   if (refVPtr->id() > maxRefVolumeId)
       maxRefVolumeId = refVPtr->id();
}

void RefEntityFactory::add(RefFace* refFPtr)
{
   assert(!refFaceList->is_in_list(refFPtr));
   if(refFaceList->size() > 0)
   {
    refFaceList->last();
    if(ManageListSorting && refFPtr->id() < refFaceList->get()->id())
      refFaceListIsSorted = false;
   }
   refFaceList->append(refFPtr);
   if (refFPtr->entityId > maxRefFaceId)
       maxRefFaceId = refFPtr->entityId;
}

void RefEntityFactory::add(RefEdge* refEPtr)
{
   assert(!refEdgeList->is_in_list(refEPtr));
   if(refEdgeList->size() > 0)
   {
    refEdgeList->last();
    if(ManageListSorting && refEPtr->id() < refEdgeList->get()->id())
      refEdgeListIsSorted = false;
   }
   refEdgeList->append(refEPtr);
   if (refEPtr->id() > maxRefEdgeId)
       maxRefEdgeId = refEPtr->id();
}

void RefEntityFactory::add(RefVertex* refVPtr)
{
   assert(!refVertexList->is_in_list(refVPtr));
   if(refVertexList->size() > 0)
   {
    refVertexList->last();
    if(ManageListSorting && refVPtr->id() < refVertexList->get()->id())
      refVertexListIsSorted = false;
   }
   refVertexList->append(refVPtr);
   if (refVPtr->id() > maxRefVertexId)
       maxRefVertexId = refVPtr->id();
}


void RefEntityFactory::remove(RefEntity* ref_entity)
{
  assert(ref_entity != 0);
  
  RefGroup *group;
  RefVolume *volume;
  RefFace *face;
  RefEdge *edge;
  RefVertex *vertex;
  Body *body;

  if ( (group = CAST_TO(ref_entity, RefGroup)) != NULL ) {
    remove(group);
  }
  else if ( (volume = CAST_TO(ref_entity, RefVolume)) != NULL ) {
    remove(volume);
  }
  else if ( (face = CAST_TO(ref_entity, RefFace)) != NULL ) {
    remove(face);
  }
  else if ( (edge = CAST_TO(ref_entity, RefEdge)) != NULL ) {
    remove(edge);
  }
  else if ( (vertex = CAST_TO(ref_entity, RefVertex)) != NULL ) {
    remove(vertex);
  }
  else if ( (body = CAST_TO(ref_entity, Body)) != NULL ) {
    remove(body);
  }
}

void RefEntityFactory::remove (RefVertex *ref_vertex_ptr)
{
  if (!refVertexList) return;
  
   if (ref_vertex_ptr != NULL)
   {
      if (refVertexList->move_to ( ref_vertex_ptr )) 
          refVertexList->remove ();
      else assert(0);
   }
}

void RefEntityFactory::remove (RefEdge *ref_edge_ptr)
{
  if (!refEdgeList) return;

   if (ref_edge_ptr != NULL)
   {
      if (refEdgeList->move_to ( ref_edge_ptr ))
          refEdgeList->remove ();
      else assert(0);
   }
}

void RefEntityFactory::remove (RefFace *ref_face_ptr)
{
  if (!refFaceList) return;

   if (ref_face_ptr != NULL)
   {
      if (refFaceList->move_to ( ref_face_ptr ))
          refFaceList->remove ();
      else assert(0);
   }
}

void RefEntityFactory::remove (RefVolume *ref_volume_ptr)
{
  if (!refVolumeList) return;

  if (ref_volume_ptr != NULL)
   {
      if (refVolumeList->move_to ( ref_volume_ptr ))
          refVolumeList->remove ();
      else assert(0);
   }
}

void RefEntityFactory::remove (RefGroup *ref_group_ptr)
{
  if (!refGroupList) return;

  if (ref_group_ptr != NULL)
  {
    if (refGroupList->move_to ( ref_group_ptr ))
      refGroupList->remove ();
    else assert(0);
  }
}

void RefEntityFactory::remove (Body *body_ptr)
{
  if (!bodyList) return;

  if (body_ptr != NULL)
  {
    if (bodyList->move_to ( body_ptr ))
      bodyList->remove ();
    else assert(0);
  }
}

void RefEntityFactory::bodies(DLIList<Body*> &bodies) 
{
  bodies = *bodyList;
}

void RefEntityFactory::ref_volumes(DLIList<RefVolume*> &ref_volumes)
{
  ref_volumes = *refVolumeList;
}

void RefEntityFactory::ref_groups(DLIList<RefGroup*> &ref_groups)
{
  ref_groups = *refGroupList;
}

void RefEntityFactory::ref_faces(DLIList<RefFace*> &ref_faces)
{
  ref_faces = *refFaceList;
}

void RefEntityFactory::ref_edges(DLIList<RefEdge*> &ref_edges)
{
  ref_edges = *refEdgeList;
}

void RefEntityFactory::ref_vertices(DLIList<RefVertex*> &ref_vertices)
{
  ref_vertices = *refVertexList;
}

#ifdef PROE
void RefEntityFactory::ref_parts(DLIList<RefPart*> &ref_parts)
{
	ref_parts = *refPartList;
}

void RefEntityFactory::ref_assemblies(DLIList<RefAssembly*> &ref_assemblies)
{
	ref_assemblies = *refAssemblyList;
}
#endif

int RefEntityFactory::num_bodies() const {return bodyList->size();}
int RefEntityFactory::num_ref_volumes() const {return refVolumeList->size();}
int RefEntityFactory::num_ref_groups() const {return refGroupList->size();}
int RefEntityFactory::num_ref_faces() const {return refFaceList->size();}
int RefEntityFactory::num_ref_edges() const {return refEdgeList->size();}
int RefEntityFactory::num_ref_vertices() const {return refVertexList->size();}

RefEntity* RefEntityFactory::get_ref_entity(const char *type_string, int id)
{
     // Call the right function, cast it.
  switch( toupper(*type_string) )
  {
    case 'B':
      if( !strcmp(type_string+1,"ody") )
        return get_body(id);
      break;
    case 'V':
      if( !strcmp(type_string+1,"olume" ) )
        return get_ref_volume(id);
      else if(!strcmp(type_string+1,"ertex") )
        return get_ref_vertex(id);
      break;
    case 'S':
      if( !strcmp(type_string+1,"urface") )
        return get_ref_face(id);
      break;
    case 'C':
      if( !strcmp(type_string+1,"urve") )
        return get_ref_edge(id);
      break;
    case 'G':
      if( !strcmp(type_string+1,"roup") )
        return get_ref_group(id);
      break;
  }   

  PRINT_ERROR("Invalid list type for the get_ref_entity "
              "function: %s\n", type_string);
  return NULL;
}
//_________  Add Code by DZ of Cat,  3/17/99 11:24:22 AM  _________
RefEntity* RefEntityFactory::get_ref_entity (const type_info& type, int id)
{
   // Call the right function, cast it.
  if( type == typeid(RefGroup) )
     return CAST_TO(get_ref_group(id), RefEntity);
  else if( type == typeid(Body) )
     return CAST_TO(get_body(id), RefEntity);
  else if( type == typeid(RefVolume) )
     return CAST_TO(get_ref_volume(id), RefEntity);
  else if( type == typeid(RefFace) )
     return CAST_TO(get_ref_face(id), RefEntity);
  else if( type == typeid(RefEdge) )
     return CAST_TO(get_ref_edge(id), RefEntity);
  else if( type == typeid(RefVertex) )
     return CAST_TO(get_ref_vertex(id), RefEntity);
      
   return NULL;
}
//_________  Code End by DZ of Cat,  3/17/99 11:24:22 AM  _________

Body* RefEntityFactory::get_body (int id)
{
   if (!bodyList->size())
       return NULL;
   
   // Make sure the list is sorted for the binary search.
   if(ManageListSorting && !bodyListIsSorted)
   {
     bodyList->sort(sort_by_ascending_ids);
     bodyListIsSorted = true;
   }

   // Do a binary search on the sorted list.  We are making the assumption
   // here that there are no NULL entries in the list.
   bool found = false;
   int left_index = 0;
   int right_index = bodyList->size()-1;
   int mid_index = (left_index + right_index)/2;
   int mid_id = ((*bodyList)[mid_index])->id();

   while(!found && (right_index-left_index) > 1)
   {
     if(mid_id == id)
       found = true;
     else
     {
       if(mid_id > id)
         right_index = mid_index;
       else
         left_index = mid_index;
       mid_index = (left_index + right_index)/2;
       mid_id = ((*bodyList)[mid_index])->id();
     }
   }

   if(!found)
   {
     if(((*bodyList)[left_index])->id() == id)
       return ((*bodyList)[left_index]);
     else if(((*bodyList)[right_index])->id() == id)
       return ((*bodyList)[right_index]);
   }
   else
     return ((*bodyList)[mid_index]);
   
   return NULL ;
}

RefGroup* RefEntityFactory::get_ref_group   (int id)
{
   if (!refGroupList->size())
       return NULL;
   
   // Make sure the list is sorted for the binary search.
   if(ManageListSorting && !refGroupListIsSorted)
   {
     refGroupList->sort(sort_by_ascending_ids);
     refGroupListIsSorted = true;
   }

   // Do a binary search on the sorted list.  We are making the assumption
   // here that there are no NULL entries in the list.
   bool found = false;
   int left_index = 0;
   int right_index = refGroupList->size()-1;
   int mid_index = (left_index + right_index)/2;
   int mid_id = ((*refGroupList)[mid_index])->id();

   while(!found && (right_index-left_index) > 1)
   {
     if(mid_id == id)
       found = true;
     else
     {
       if(mid_id > id)
         right_index = mid_index;
       else
         left_index = mid_index;
       mid_index = (left_index + right_index)/2;
       mid_id = ((*refGroupList)[mid_index])->id();
     }
   }

   if(!found)
   {
     if(((*refGroupList)[left_index])->id() == id)
       return ((*refGroupList)[left_index]);
     else if(((*refGroupList)[right_index])->id() == id)
       return ((*refGroupList)[right_index]);
   }
   else
     return ((*refGroupList)[mid_index]);

   return NULL ;
}

RefVolume* RefEntityFactory::get_ref_volume (int id) 
{
   if (!refVolumeList->size())
       return NULL;
   
   // Make sure the list is sorted for the binary search.
   if(ManageListSorting && !refVolumeListIsSorted)
   {
     refVolumeList->sort(sort_by_ascending_ids);
     refVolumeListIsSorted = true;
   }

   // Do a binary search on the sorted list.  We are making the assumption
   // here that there are no NULL entries in the list.
   bool found = false;
   int left_index = 0;
   int right_index = refVolumeList->size()-1;
   int mid_index = (left_index + right_index)/2;
   int mid_id = ((*refVolumeList)[mid_index])->id();

   while(!found && (right_index-left_index) > 1)
   {
     if(mid_id == id)
       found = true;
     else
     {
       if(mid_id > id)
         right_index = mid_index;
       else
         left_index = mid_index;
       mid_index = (left_index + right_index)/2;
       mid_id = ((*refVolumeList)[mid_index])->id();
     }
   }

   if(!found)
   {
     if(((*refVolumeList)[left_index])->id() == id)
       return ((*refVolumeList)[left_index]);
     else if(((*refVolumeList)[right_index])->id() == id)
       return ((*refVolumeList)[right_index]);
   }
   else
     return ((*refVolumeList)[mid_index]);

   return NULL ;
}

RefFace* RefEntityFactory::get_ref_face (int id) 
{
   if (!refFaceList->size())
       return NULL;
   
   // Make sure the list is sorted for the binary search.
   if(ManageListSorting && !refFaceListIsSorted)
   {
     refFaceList->sort(sort_by_ascending_ids);
     refFaceListIsSorted = true;
   }

   // Do a binary search on the sorted list.  We are making the assumption
   // here that there are no NULL entries in the list.
   bool found = false;
   int left_index = 0;
   int right_index = refFaceList->size()-1;
   int mid_index = (left_index + right_index)/2;
   int mid_id = ((*refFaceList)[mid_index])->id();

   while(!found && (right_index-left_index) > 1)
   {
     if(mid_id == id)
       found = true;
     else
     {
       if(mid_id > id)
         right_index = mid_index;
       else
         left_index = mid_index;
       mid_index = (left_index + right_index)/2;
       mid_id = ((*refFaceList)[mid_index])->id();
     }
   }

   if(!found)
   {
     if(((*refFaceList)[left_index])->id() == id)
       return ((*refFaceList)[left_index]);
     else if(((*refFaceList)[right_index])->id() == id)
       return ((*refFaceList)[right_index]);
   }
   else
     return ((*refFaceList)[mid_index]);

   return NULL ;
}

RefEdge* RefEntityFactory::get_ref_edge (int id) 
{
   if (!refEdgeList->size())
       return NULL;
   
   // Make sure the list is sorted for the binary search.
   if(ManageListSorting && !refEdgeListIsSorted)
   {
     refEdgeList->sort(sort_by_ascending_ids);
     refEdgeListIsSorted = true;
   }

   // Do a binary search on the sorted list.  We are making the assumption
   // here that there are no NULL entries in the list.
   bool found = false;
   int left_index = 0;
   int right_index = refEdgeList->size()-1;
   int mid_index = (left_index + right_index)/2;
   int mid_id = ((*refEdgeList)[mid_index])->id();

   while(!found && (right_index-left_index) > 1)
   {
     if(mid_id == id)
       found = true;
     else
     {
       if(mid_id > id)
         right_index = mid_index;
       else
         left_index = mid_index;
       mid_index = (left_index + right_index)/2;
       mid_id = ((*refEdgeList)[mid_index])->id();
     }
   }

   if(!found)
   {
     if(((*refEdgeList)[left_index])->id() == id)
       return ((*refEdgeList)[left_index]);
     else if(((*refEdgeList)[right_index])->id() == id)
       return ((*refEdgeList)[right_index]);
   }
   else
     return ((*refEdgeList)[mid_index]);

   return NULL ;
}

RefVertex* RefEntityFactory::get_ref_vertex (int id) 
{
   if (!refVertexList->size())
       return NULL;
   
   // Make sure the list is sorted for the binary search.
   if(ManageListSorting && !refVertexListIsSorted)
   {
     refVertexList->sort(sort_by_ascending_ids);
     refVertexListIsSorted = true;
   }

   // Do a binary search on the sorted list.  We are making the assumption
   // here that there are no NULL entries in the list.
   bool found = false;
   int left_index = 0;
   int right_index = refVertexList->size()-1;
   int mid_index = (left_index + right_index)/2;
   int mid_id = ((*refVertexList)[mid_index])->id();

   while(!found && (right_index-left_index) > 1)
   {
     if(mid_id == id)
       found = true;
     else
     {
       if(mid_id > id)
         right_index = mid_index;
       else
         left_index = mid_index;
       mid_index = (left_index + right_index)/2;
       mid_id = ((*refVertexList)[mid_index])->id();
     }
   }

   if(!found)
   {
     if(((*refVertexList)[left_index])->id() == id)
       return ((*refVertexList)[left_index]);
     else if(((*refVertexList)[right_index])->id() == id)
       return ((*refVertexList)[right_index]);
   }
   else
     return ((*refVertexList)[mid_index]);

   return NULL ;
}

//* Methods: next*Id
//** Increments the maximum object ID and returns this value.
//==
int RefEntityFactory::next_body_id () 
{
   return ++maxBodyId;
}

int RefEntityFactory::next_ref_group_id () 
{
   return ++maxRefGroupId;
}

int RefEntityFactory::next_ref_volume_id    () 
{
   return ++maxRefVolumeId;
}

int RefEntityFactory::next_ref_face_id () 
{
   return ++maxRefFaceId;
}

int RefEntityFactory::next_ref_edge_id () 
{
   return ++maxRefEdgeId;
}

int RefEntityFactory::next_ref_vertex_id () 
{
   return ++maxRefVertexId;
}

#ifdef PROE
int RefEntityFactory::next_ref_assembly_id ()
{
	return ++maxRefAssemblyId;
}

int RefEntityFactory::next_ref_part_id ()
{
	return ++maxRefPartId;
}
#endif

Body *RefEntityFactory::get_first_body()
{
	bodyList->reset(); 
	return bodyList->size() ? bodyList->get() : 0;
}

RefVolume *RefEntityFactory::get_first_ref_volume()
{
	refVolumeList->reset(); 
	return refVolumeList->size() ? refVolumeList->get() : 0;
}

RefGroup *RefEntityFactory::get_first_ref_group()
{
	refGroupList->reset(); 
	return refGroupList->size() ? refGroupList->get() : 0;
}

RefFace *RefEntityFactory::get_first_ref_face()
{
	refFaceList->reset(); 
	return refFaceList->size() ? refFaceList->get() : 0;
}

RefEdge *RefEntityFactory::get_first_ref_edge()
{
	refEdgeList->reset(); 
	return refEdgeList->size() ? refEdgeList->get() : 0;
}

RefVertex *RefEntityFactory::get_first_ref_vertex()
{
	refVertexList->reset(); 
	return refVertexList->size() ? refVertexList->get() : 0;
}

        
Body *RefEntityFactory::get_next_body()
{return bodyList->size() ? bodyList->step_and_get() : 0;}

RefVolume *RefEntityFactory::get_next_ref_volume()
{return refVolumeList->size() ? refVolumeList->step_and_get() : 0;}

RefGroup *RefEntityFactory::get_next_ref_group()
{return refGroupList->size() ? refGroupList->step_and_get() : 0;}

RefFace *RefEntityFactory::get_next_ref_face()
{return refFaceList->size() ? refFaceList->step_and_get() : 0;}

RefEdge *RefEntityFactory::get_next_ref_edge()
{return refEdgeList->size() ? refEdgeList->step_and_get() : 0;}

RefVertex *RefEntityFactory::get_next_ref_vertex()
{return refVertexList->size() ? refVertexList->step_and_get() : 0;}

Body *RefEntityFactory::get_last_body()
{
	bodyList->last(); 
	return bodyList->size() ? bodyList->get() : 0;
}

RefVolume *RefEntityFactory::get_last_ref_volume()
{
	refVolumeList->last(); 
	return refVolumeList->size() ? refVolumeList->get() : 0;
}

RefGroup *RefEntityFactory::get_last_ref_group()
{
	refGroupList->last(); 
	return refGroupList->size() ? refGroupList->get() : 0;
}

RefFace *RefEntityFactory::get_last_ref_face()
{
	refFaceList->last(); 
	return refFaceList->size() ? refFaceList->get() : 0;
}

RefEdge *RefEntityFactory::get_last_ref_edge()
{
	refEdgeList->last(); 
	return refEdgeList->size() ? refEdgeList->get() : 0;
}

RefVertex *RefEntityFactory::get_last_ref_vertex()
{
	refVertexList->last(); 
	return refVertexList->size() ? refVertexList->get() : 0;
}

void RefEntityFactory::incorporate_id (RefEntity *ref_ent)
{
  if (ref_ent)
  {
//________  Change Code by DZ of Cat,  3/17/99 10:47:39 AM  ________
     int max_ent_id = 0;

     max_ent_id = ref_ent->id();
        // Do nothing with groups
     if( CAST_TO( ref_ent, RefGroup ) ) {
       if (maxRefGroupId < max_ent_id)
         maxRefGroupId = max_ent_id;
     }
     else if( CAST_TO( ref_ent, Body ) )
     {
       if (maxBodyId < max_ent_id)
          maxBodyId = max_ent_id;
     }
     else if( CAST_TO( ref_ent, RefVolume ) )
     {
       if (maxRefVolumeId < max_ent_id)
          maxRefVolumeId = max_ent_id;
     }
     else if( CAST_TO( ref_ent, RefFace ) )
     {
       if (maxRefFaceId < max_ent_id)
          maxRefFaceId = max_ent_id;
     }
     else if( CAST_TO( ref_ent, RefEdge ) )
     {
       if (maxRefEdgeId < max_ent_id)
          maxRefEdgeId = max_ent_id;
     }
     else if( CAST_TO( ref_ent, RefVertex ) )
     {
       if (maxRefVertexId < max_ent_id)
          maxRefVertexId = max_ent_id;
     }
//________  Change End by DZ of Cat,  3/17/99 10:47:39 AM  ________
  }
}

int RefEntityFactory::maximum_id(const char *entity_type)
{
   if (strcmp("body", entity_type) == 0)
       return maxBodyId;
   else if (strcmp("curve", entity_type) == 0)
       return maxRefEdgeId;
   else if (strcmp("group", entity_type) == 0)
       return maxRefGroupId;
   else if (strcmp("volume", entity_type) == 0)
       return maxRefVolumeId;
   else if (strcmp("vertex", entity_type) == 0)
       return maxRefVertexId;
   else if (strcmp("surface", entity_type) == 0)
       return maxRefFaceId;
   else {
      PRINT_ERROR("Unrecognized entity_type: '%s'\n", entity_type);
      return 0;
   }
}

void RefEntityFactory::maximum_id (const type_info& type, int max_id)
{
  if( type == typeid(RefGroup) )
     maxRefGroupId = max_id;
  else if( type == typeid(Body) )
     maxBodyId = max_id;
  else if( type == typeid(RefVolume) )
     maxRefVolumeId = max_id;
  else if( type == typeid(RefFace) )
     maxRefFaceId = max_id;
  else if( type == typeid(RefEdge) )
     maxRefEdgeId = max_id;
  else if( type == typeid(RefVertex) )
     maxRefVertexId = max_id;
}
      
int RefEntityFactory::maximum_id (RefEntity *ref_ent)
{
   if (!ref_ent)
      return 0;
   else if( CAST_TO( ref_ent, RefGroup ) )
      return maxRefGroupId;
   else if( CAST_TO( ref_ent, Body ) )
      return maxBodyId;
   else if( CAST_TO( ref_ent, RefVolume ) )
      return maxRefVolumeId;
   else if( CAST_TO( ref_ent, RefFace ) )
      return maxRefFaceId;
   else if( CAST_TO( ref_ent, RefEdge ) )
      return maxRefEdgeId;
   else if( CAST_TO( ref_ent, RefVertex ) )
      return maxRefVertexId;
   else
      return 0;
}

void RefEntityFactory::compress_ref_ids(const char *entity_type, int retain_max_id)
{
  const type_info& type = RefEntity::get_entity_type_info(entity_type);
  assert(type != typeid(InvalidEntity));
  
  /*
#ifdef VIRTUAL_GEOMETRY_ENGINE_HPP
  VirtualGeometryEngine::instance()->
    hidden_entity_mngr.compress_hidden_ids(type);
#endif
  */
/*  
  DLIList<GeometryQueryEngine*> gqeList;
  GeometryQueryTool::instance()->get_gqe_list(gqeList);//Get the gqeList from the GQT 
  gqeList.reset();
  int i;
  for (i = 0; i < gqeList.size(); i++)//Step through the list and call compress_ids.
  {                                   //The VGE is the only engine that will do work. 
    gqeList.get_and_step()->compress_ids(type);
  }
*/

  DLIList<CubitEntity*> list;
  DLIList<RefEntity*> temp_list;
  CubitStatus result = ref_entity_list(entity_type, temp_list, CUBIT_FALSE);
  if (  result == CUBIT_SUCCESS ) {
    CAST_LIST_TO_PARENT(temp_list, list);
    compress_ids(list);
  }
  
    // set the maximum entity id to the new max values.
  if (!retain_max_id)
  {
        // must use the num_xxx() functions here, since we may have a derived factory!
    if( type == typeid(RefGroup) )
       maxRefGroupId = num_ref_groups();
    else if( type == typeid(Body) )
       maxBodyId = num_bodies();
    else if( type == typeid(RefVolume) )
       maxRefVolumeId = num_ref_volumes();
    else if( type == typeid(RefFace) )
       maxRefFaceId = num_ref_faces();
    else if( type == typeid(RefEdge) )
       maxRefEdgeId = num_ref_edges();
    else if( type == typeid(RefVertex) )
       maxRefVertexId = num_ref_vertices();
  }
}

void RefEntityFactory::compress_ids(DLIList<CubitEntity*> &list)
{
  int id = 1;
  CubitEntity* entity;
  
  if (list.size())
  {
    list.reset();
      // if these are ref volumes, recompute color
    CubitBoolean set_color = CUBIT_FALSE;
    if (CAST_TO(list.get(), RefVolume)) set_color = CUBIT_TRUE;

    // sort list so ids of entities are ascending order
    list.sort( sort_by_ascending_ids); 

    for (int i=list.size(); i > 0; i--)
    {
      entity = list.get_and_step();
      if (entity->id() != id){
        entity->set_id(id);
        if (set_color) {
          entity->color(CUBIT_DEFAULT_COLOR);
        }
      }
      id++;
    }
  }
}

void RefEntityFactory::reset_ids()
{
  maxBodyId           = 0;
  maxRefVolumeId      = 0;
  maxRefGroupId       = 0;
  maxRefFaceId        = 0;
  maxRefEdgeId        = 0;
  maxRefVertexId      = 0;
  maxSurfSubDomainId  = 0;
  maxCurveSubDomainId = 0;
  maxRefCoordSysId    = 0;
#ifdef PROE
  maxRefAssemblyId    = 0;
  maxRefPartId        = 0;
#endif
}

#if 0
static int re_factory_sort_volumes(RefVolume*& /*a*/, RefVolume*& /*b*/)
{
    // This is currently not implemented.  When it is added to CGM,
    // look at MRefEntityFactory to see what to do.
  return 0;
}
#endif

void RefEntityFactory::renumber_geometry_by_properties(CubitBoolean retain_max)
{
    // See MRefEntityFactory::renumber_geometry_by_properties.
    // If you want this implemented, copy what was done there.
  PRINT_ERROR("'Sort' option not supported.\n"
              "       Compressing IDs without sorting.");
  compress_ref_ids("group", retain_max);
  compress_ref_ids("body", retain_max);
  compress_ref_ids("volume", retain_max);
  compress_ref_ids("surface", retain_max);
  compress_ref_ids("curve", retain_max);
  compress_ref_ids("vertex", retain_max);
}

CubitStatus RefEntityFactory::notify_observer(CubitObservable *observable,
                                              const CubitEvent &observer_event)
{
  RefEntity *entity = CAST_TO(observable, RefEntity);
  if (!entity) return CUBIT_FAILURE;
  
    //- handle MODEL_ENTITY_DESTRUCTED/MODEL_ENTITY_CONSTRUCTED events
  int event = observer_event.get_event_type();
  if (event == MODEL_ENTITY_CONSTRUCTED) 
    add(entity);
  else if (event == MODEL_ENTITY_DESTRUCTED) 
    remove(entity);
  else if (event == ID_SET) 
  {
    if(CAST_TO(entity, RefEdge) && ManageListSorting)
      refEdgeListIsSorted = false;
    else if(CAST_TO(entity, RefFace) && ManageListSorting)
      refFaceListIsSorted = false;
    else if(CAST_TO(entity, RefVertex) && ManageListSorting)
      refVertexListIsSorted = false;
    else if(CAST_TO(entity, RefVolume) && ManageListSorting)
      refVolumeListIsSorted = false;
    else if(CAST_TO(entity, RefGroup) && ManageListSorting)
      refGroupListIsSorted = false;
    else if(CAST_TO(entity, Body) && ManageListSorting)
      bodyListIsSorted = false;
  }

  return CUBIT_SUCCESS;
}

