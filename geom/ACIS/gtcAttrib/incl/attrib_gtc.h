// @(#)attrib_gtc.hxx	1.1	11/10/89

// Attribute declaration for a private container attribute. Each
// application developer receives one of these customised for himself.
// All attributes specific to the application developer are then made
// derived classes of this attribute, ensuring that different
// developers can assign identifiers independently without mutual
// interference.

// This one is for Goodyear Tire and Rubber Company
// Author: Arlo L. Ames, Sandia National Laboratories

#if !defined( ATTRIB_GTC_CLASS )
#define ATTRIB_GTC_CLASS

#include "attrib.hxx"
#include "attrib_gtc_export.h"
// This attribute type is a derived class of ATTRIB.

extern DECL_GTCATTRIB int ATTRIB_GTC_TYPE;

#define ATTRIB_GTC_LEVEL (ATTRIB_LEVEL + 1) 
MASTER_ATTRIB_DECL(ATTRIB_GTC, GTCATTRIB)

#endif

