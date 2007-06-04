/**********************************************************************/
/*    ACIS (R) is a registered trademark of Spatial Technology Inc.   */
/*                                                                    */
/*      This program is the sole property  of Spatial Technology      */
/*      Inc. and Three-Space Ltd.  and is protected by copyright      */
/*      under the  laws of  the United  States.  This program is      */
/*      confidential, proprietary, and a trade secret, not to be      */
/*      disclosed  without  written authorization  from  Spatial      */
/*      Technology   Inc.  or   Three-Space    Ltd.    Any  use,      */
/*      duplication, or disclosure  of  this  program  by  other      */
/*      than  Spatial Technology  Inc. or  Three-Space Ltd., and      */
/*      their  assigned  licensees  and  customers  is  strictly      */
/*      forbidden by law.                                             */
/*                                                                    */
/*       Copyright (c) 1987-1993 by                                   */
/*           Spatial Technology Inc. and Three-Space Ltd.             */
/*                       All rights reserved                          */
/**********************************************************************/
// @(#)attrib_snl.cxx	1.1	11/10/89

// Implementation of container attribute for a specific application
// developer.

// This one is for Three-Space Ltd internal use.

#include <stdio.h>

// ACIS Includes
#if CUBIT_ACIS_VERSION < 1100
#include "kernel/dcl_kern.h"
#include "kernel/kerndata/attrib/at_tsl.hxx"
#include "kernel/kerndata/data/datamsc.hxx"
#else
#include "dcl_kern.h"
#include "at_tsl.hxx"
#include "datamsc.hxx"
#endif

#include "attrib_snl.hpp"

// Define macros for this attribute and its parent, to provide the
// information to the definition macro.

#define THIS() ATTRIB_SNL
#define THIS_LIB NONE
#define PARENT() ATTRIB
#define PARENT_LIB KERN

 
// Identifier used externally to identify an particular entity type.
// This is only used within the save/restore system for translating
// to/from external file format, but must be unique amongst attributes
// derived directly from ATTRIB, across all application developers.

#define ATTRIB_SNL_NAME "snl"

MASTER_ATTRIB_DEFN( "snl master attribute" )

