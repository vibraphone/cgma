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
// @(#)attrib_snl.hpp   1.1     11/10/89

// Attribute declaration for a private container attribute. Each
// application developer receives one of these customised for himself.
// All attributes specific to the application developer are then made
// derived classes of this attribute, ensuring that different
// developers can assign identifiers independently without mutual
// interference.

// This one is for Sandia Labs

#if !defined( ATTRIB_SNL_CLASS )
#define ATTRIB_SNL_CLASS

// CUBIT Includes
#include "decl_none.h"
// ACIS Includes
#if CUBIT_ACIS_VERSION < 1100
#include "kernel/kerndata/attrib/attrib.hxx"
#else
#include "attrib.hxx"
#endif

// This attribute type is a derived class of ATTRIB.

extern DECL_NONE int ATTRIB_SNL_TYPE;

#define ATTRIB_SNL_LEVEL (ATTRIB_LEVEL + 1)

MASTER_ATTRIB_DECL( ATTRIB_SNL, NONE)

#endif
