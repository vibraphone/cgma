// Implementation of the mark attribute

#include "attrib_mark.h"


#include "api.hxx"
#include "alltop.hxx"
#include "datamsc.hxx"

//Define macros for this attribute and its parent, to simplify later
// stuff and make it suitable for making into macros.

#define THIS() ATTRIB_MARK
#define PARENT() ATTRIB_GTC
#define PARENT_LIB GTCATTRIB
#define THIS_LIB GTCATTRIB


// Identifier used externally to identify a particular entity type.
// This is used in the save/restore system for translating to/from
// external file format.

#define ATTRIB_MARK_NAME "mark"

// Implement the standard attribute functions.

ATTCOPY_DEF( "mark_attribute" )

LOSE_DEF

DTOR_DEF

DEBUG_DEF

SAVE_DEF
   write_int(mark_data); //Save specific data

RESTORE_DEF
   {
      int buffer;
      buffer = read_int();
      set_mark(buffer);
   }

COPY_DEF
   set_mark(from->mark() );

SCAN_DEF
   // (no specific pointer data)

FIX_POINTER_DEF
   // (no specific pointer data)

TERMINATE_DEF

// make a mark attribute

ATTRIB_MARK::ATTRIB_MARK(
   ENTITY* owner,
   int mark_
) : ATTRIB_GTC(owner)
{
   //Initialize members
   set_mark(mark_);
}

// Set the member data.

void
ATTRIB_MARK::set_mark(
   int mark_
)
{
   backup();
   mark_data = mark_;
}

// Virtual function called when an owner entity is being split.

void
ATTRIB_MARK::split_owner(
   ENTITY* new_ent
)
{
   // Duplicate the attribute on the new entity
   new ATTRIB_MARK(new_ent, mark_data);
}

// Virtual function called when two entities are to be merged.

void
ATTRIB_MARK::merge_owner(
   ENTITY* other_ent,
   logical delete_owner
)
{
   // If the owner of this attribute is goint to be deleted, and there is
   // no attribute attached to the other entity, the we transfer ourself
   // to that other entity.

   if (delete_owner) {
      ATTRIB* other_att = find_attrib(
         other_ent,
         ATTRIB_GTC_TYPE,
         ATTRIB_MARK_TYPE
      );
      if (other_att == NULL ) {
         // No mark on other entity, so transfer ourself
         move(other_ent);
      }
   }
}


void set_mark(ENTITY* entity, int mark)
{
   // Search for an existing mark attribute and set it.

   ATTRIB_MARK *mark_att = (ATTRIB_MARK *) find_attrib(entity, ATTRIB_GTC_TYPE, ATTRIB_MARK_TYPE);

   if (mark_att != NULL)
   {
      //cout << "already had a mark\n" << endl;
      mark_att->backup();
      mark_att->set_mark(mark);
   }
   else
   {
      //cout << "adding a mark\n" << endl;
      // Create a new attribute if none found
      new ATTRIB_MARK(entity, mark);
   }
}

void set_mark(ENTITY* entity)
// mark entity with something, don't care what
{
   set_mark(entity, 1);
}

int get_mark(ENTITY* entity)
{
  int get_mark_result = 0 ;
  ATTRIB_MARK *mark_att = (ATTRIB_MARK *) find_attrib(entity, ATTRIB_GTC_TYPE,
						      ATTRIB_MARK_TYPE);
  if (mark_att) get_mark_result = mark_att->mark();
  return get_mark_result;
}

logical marked(ENTITY* entity)
{
   ATTRIB_MARK *mark_att = (ATTRIB_MARK *) find_attrib(entity, ATTRIB_GTC_TYPE, ATTRIB_MARK_TYPE);
   if (!mark_att) return FALSE;
   return TRUE;
}


void mark_edges(BODY* body, int mark_)
//mark all the edges of body with the value of mark
{
   for (LUMP* lump = body->lump(); lump; lump = lump->next())
      for (SHELL* shell = lump->shell(); shell; shell = shell->next())
         for (FACE* face = shell->first_face(); face; face = face->next_face())
            for (LOOP* loop = face->loop(); loop; loop = loop->next())
            {
               //ALA: I don't like to split up the control for loops as I've
               //done here.  Somewhere I've written a loop for all coedges
               //that places all the control in a for statement, so you can't
               //be bitten so easilty when you copy code.  I'll find it
               //and replace the following control structure.
               COEDGE* start = loop->start();
               COEDGE* coedge = start;
               while (coedge)
               {
                  EDGE* edge = coedge->edge();
                  set_mark(edge, mark_);
                  coedge = coedge->next();
                  if (coedge == start) coedge = NULL;
               }
            }
   if (body && body->lump())
      if (body->lump()->shell())
         for (WIRE* wire = body->lump()->shell()->wire(); 
              wire; wire = wire->next())
         {
            COEDGE* start = wire->coedge();
            COEDGE* coedge = start;
            while (coedge)
            {
               EDGE* edge = coedge->edge();
               set_mark(edge, mark_);
               coedge = coedge->next();
               if (coedge == start) coedge = NULL;
            }
         }
}

void mark_edges(BODY* body)
// mark all the edges of body.  Here we don't care about the value of the mark,
// only that the mark exists.
{
   mark_edges(body, 1);
}

void mark_coedges(BODY* body, int mark_)
//mark all the coedges of body with the value of mark
{
   for (LUMP* lump = body->lump(); lump; lump = lump->next())
      for (SHELL* shell = lump->shell(); shell; shell = shell->next())
         for (FACE* face = shell->first_face(); face; face = face->next_face())
            for (LOOP* loop = face->loop(); loop; loop = loop->next())
            {
               //ALA: I don't like to split up the control for loops as I've
               //done here.  Somewhere I've written a loop for all coedges
               //that places all the control in a for statement, so you can't
               //be bitten so easilty when you copy code.  I'll find it
               //and replace the following control structure.
               COEDGE* start = loop->start();
               COEDGE* coedge = start;
               while (coedge)
               {
                  set_mark(coedge, mark_);
                  coedge = coedge->next();
                  if (coedge == start) coedge = NULL;
               }
            }
   if (body && body->lump())
      if (body->lump()->shell())
         for (WIRE* wire = body->lump()->shell()->wire(); 
              wire; wire = wire->next())
         {
            COEDGE* start = wire->coedge();
            COEDGE* coedge = start;
            while (coedge)
            {
               set_mark(coedge, mark_);
               coedge = coedge->next();
               if (coedge == start) coedge = NULL;
            }
         }
}


void mark_coedges(BODY* body)
// mark all the coedges of body.  Here we don't care about the value of the mark,
// only that the mark exists.
{
   mark_coedges(body, 1);
}


void unmark_(ENTITY* entity)
{
API_BEGIN

   ATTRIB_MARK * attribute =
      (ATTRIB_MARK *) find_attrib(entity,
                                  ATTRIB_GTC_TYPE,
                                  ATTRIB_MARK_TYPE);
   if (attribute)
   {
      attribute->unhook();
      attribute->lose();
   }

API_END
}


void unmark(VERTEX* vertex)
{
   unmark_((ENTITY*) vertex);
}

void unmark(EDGE* edge)
{
   unmark_((ENTITY*) edge);
   VERTEX* start_vertex = edge->start();
   VERTEX* end_vertex = edge->end();
   unmark(start_vertex);
   unmark(end_vertex);
}

void unmark(COEDGE* coedge)
{
   unmark_((ENTITY*) coedge);
   EDGE* edge = coedge->edge();
   unmark(edge);
}

void unmark(LOOP* loop)
{
   unmark_((ENTITY*) loop);
   COEDGE* coedge = NULL;
   COEDGE* first_coedge = loop->start();
   for (coedge = loop->start(); coedge; coedge = coedge->next(), coedge = (coedge == first_coedge)?NULL:coedge)
      unmark(coedge);
}

void unmark(FACE* face)
{
   unmark_((ENTITY*) face);
   for (LOOP* loop = face->loop(); loop; loop = loop->next())
      unmark(loop);
}

void unmark(SHELL* shell)
{
   unmark_((ENTITY*) shell);
   for (FACE* face = shell->first_face(); face; face = face->next_face())
      unmark(face);
}

void unmark(LUMP* lump)
{
   unmark_((ENTITY*) lump);
   for (SHELL * shell = lump->shell(); shell; shell = shell->next())
      unmark(shell);
}

void unmark(WIRE* wire)
{
   ENTITY_LIST wire_coedges;
   wire_coedges.add(wire->coedge());

   int i=0;
   while (wire_coedges[i] != NULL)
   {
      COEDGE* coedge = (COEDGE *) wire_coedges[i];
      unmark(coedge);
      i++;
      wire_coedges.add(coedge->next());
      wire_coedges.add(coedge->previous());
   }
}

void unmark(BODY* body)
{
   unmark_((ENTITY*) body);
   for (LUMP* lump = body->lump(); lump; lump = lump->next())
      unmark(lump);
   if (body && body->lump())
      if (body->lump()->shell())
         for (WIRE* wire = body->lump()->shell()->wire(); 
              wire; wire = wire->next())
            unmark(wire);
}

