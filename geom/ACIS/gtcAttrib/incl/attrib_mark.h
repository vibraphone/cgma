#if !defined( ATTRIB_MARK_CLASS )
#define ATTRIB_MARK_CLASS

#include "attrib_gtc.h"
#include "attrib_gtc_export.h"

#include "logical.h"

class EDGE ;
class COEDGE ;
class FACE ;
class LOOP ;
class LUMP ;
class SHELL ;
class VERTEX ;
class WIRE ;


extern DECL_GTCATTRIB int ATTRIB_MARK_TYPE;
#define ATTRIB_MARK_LEVEL (ATTRIB_GTC_LEVEL + 1)

class DECL_GTCATTRIB ATTRIB_MARK: public ATTRIB_GTC {
   int mark_data;

public:

ATTRIB_MARK( ENTITY* = NULL, int = 0);

int mark() const { return mark_data; }
void set_mark(int);

// Functions called to aid attribute migration during modeling
// operations.

virtual void split_owner( ENTITY * );

// The owner of this attribute is about to be split in two - the
// argument is the new piece.  The owner of this attribute is about to
// be merged with the given entity.  The logical argument is TRUE if the
// owner is to be deleted in the merge.

   virtual void merge_owner( ENTITY *, // "other entity"
                             logical   // deleting_owner
                           );

   ATTRIB_FUNCTIONS( ATTRIB_MARK, GTCATTRIB )
};

DECL_GTCATTRIB void set_mark(ENTITY*, int);
DECL_GTCATTRIB void set_mark(ENTITY*);
DECL_GTCATTRIB int get_mark(ENTITY*);
DECL_GTCATTRIB logical marked(ENTITY*);

DECL_GTCATTRIB void mark_edges(ENTITY*, int);
DECL_GTCATTRIB void mark_edges(BODY*);
DECL_GTCATTRIB void mark_coedges(BODY*);

// functions to unmark an entity and everything underneath it

DECL_GTCATTRIB void unmark(VERTEX*);
DECL_GTCATTRIB void unmark(EDGE*);
DECL_GTCATTRIB void unmark(COEDGE*);
DECL_GTCATTRIB void unmark(LOOP*);
DECL_GTCATTRIB void unmark(FACE*);
DECL_GTCATTRIB void unmark(SHELL*);
DECL_GTCATTRIB void unmark(LUMP*);
DECL_GTCATTRIB void unmark(WIRE*);
DECL_GTCATTRIB void unmark(BODY*);

// function to remove a mark from a single entity
DECL_GTCATTRIB void unmark_(ENTITY*);
#endif

