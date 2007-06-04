/************************************************************************\
 * Copyright (c) 2002 Sandia Corporation.                               *
 *					         			*
 * Under the terms of Contract DE-AC04-94AL85000, there is a            *
 * non-exclusive license for use of this work by or on behalf of the    *
 * U.S. Government. Export of this program may require a license from   *
 * the United States Government.                                        *
 *					         			*
 * This software is the property of Sandia Corporation and discloses    *
 * material protectable under copyright laws of the United States.      *
 * Use, Duplication, or Disclosure is prohibited, unless allowed        *
 * subject to the terms of a separate license agreement.                *
 *					         			*
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES          *
 * DEPARTMENT OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR       *
 * EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY    *
 * LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS,    *
 * OR USEFULNESS OF ANY INFORMATION, APPARATUS OR PROCESS DISCLOSED,    *
 * OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.*
 *					         			*
\************************************************************************/

/***************************************************************************
   Filename      : mm_cgm_typedefs.h
  
   Purpose       : This file contains the typedefs associated with the 
                 : CGM C Interface functions.
  
   Special Notes :
  
   Creator       : Robert A. Kerr
  
   Creation Date : 09/26/2001
  
   Current Owner : Robert A. Kerr
****************************************************************************/

#ifndef _MM_CGM_TYPEDEFS_H
#define _MM_CGM_TYPEDEFS_H

/******************** BEGIN STANDARD INCLUDES   ****************************/
/******************** END STANDARD INCLUDES     ****************************/

/******************** BEGIN CGM INCLUDES   ****************************/
#include "CubitDefines.h"
#include "GeometryDefines.h"
/******************** END CGM INCLUDES     ****************************/

/******************** BEGIN STRUCT DECLARATIONS ****************************/

/* Handles to various CGM Topology Entities.
 * These Handles are passed through the CGM C Interface. 
*/

/*
 * Topology Entities
 */
typedef struct BodyHandle_    BodyHandle;
typedef struct VolumeHandle_  VolumeHandle;
typedef struct FaceHandle_    FaceHandle;
typedef struct EdgeHandle_    EdgeHandle;
typedef struct VertexHandle_  VertexHandle;

/*
 * Geometry Entities
 */
typedef struct PlaneHandle_  PlaneHandle;

/******************** END STRUCT DECLARATIONS   ****************************/

/******************** BEGIN EXTERN FUNCTIONS    ****************************/
/******************** END EXTERN FUNCTIONS    ****************************/

#endif

