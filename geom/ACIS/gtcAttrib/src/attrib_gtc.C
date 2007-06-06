#include <stdio.h>

#include "attrib_gtc.h"

#include "datamsc.hxx"

// Define macros for theis attribute and its parent, to provide the
// information to the definition macro.

#define THIS() ATTRIB_GTC
#define PARENT() ATTRIB

#define PARENT_LIB KERN
#define THIS_LIB GTCATTRIB


// Identifier used externally to identify a particular entity type.
// This is only used within the save/restore system for translating
// to/from external file format, but must be unique amongst
// attributes derived directly from ATTRIB, across all application
// developers

#define ATTRIB_GTC_NAME "gtc"

MASTER_ATTRIB_DEFN( "gtc master attribute" );

