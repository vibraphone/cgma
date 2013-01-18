/* Class:       GeometryDefines 
 * Description: GeometryDefines - all global geometric definitions for CUBIT.
 * Owner:       David White
 * Checked by:
 * Version:
 */

#ifndef GEOMETRY_DEFINES_HPP
#define GEOMETRY_DEFINES_HPP

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include "CubitString.hpp"

/* CUBIT_RESABS - Values less than CUBIT_RESABS are considered to be zero. */
#ifdef __cplusplus
const double CUBIT_RESABS = 1.0E-12;
#else
#define CUBIT_RESABS 1.0E-12
#endif

/* GEOMETRY_RESABS - If the distance between two points is less */
/* than GEOMETRY_RESABS the points are considered to be identical. */
#ifdef __cplusplus
const double GEOMETRY_RESABS = 1.0E-6;
const double DEFAULT_GEOM_FACTOR = 500.0;
#else
#define GEOMETRY_RESABS 1.0E-6
#define DEFAULT_GEOM_FACTOR 500.0
#endif

/* Types of solid modeler engines. */
#ifdef ACIS
	#undef ACIS
#endif
enum SolidModelerType 
{ 
   NOT_A_SOLID_MODELER,
   OCC,
   ACIS,
   PROE_GEOM,                  /* Normal Pro/E model */
   PROE_PART,
   PROE_ASSEMBLY,
   PROE_FEM_MESH_SOLID,        /* Pro/Mesh models... */
   PROE_FEM_MESH_SHELL,
   PROE_FEM_MESH_MIXED,
   PROE_FEM_MESH_BOUNDARY,
   PROE_FEM_MESH_QUILT
} ;

/* Types of underlying geometric representations. */
enum GeometricRepresentationType 
{
   GEOMTYPE_NONE,
   SOLID_MODEL,
   FACETTED,
   MESH_BASED,
   COMPOSITE
};

/* Types of geometry for the various GeometryEntities. */
enum GeometryType
{
    /* Point types */
  UNDEFINED_POINT_TYPE,
  
    /* Curve types */
  ARC_CURVE_TYPE,
  ELLIPSE_CURVE_TYPE,
  PARABOLA_CURVE_TYPE,
  SEGMENTED_CURVE_TYPE, /*curve defined by a chain of linear segments. */
  SPLINE_CURVE_TYPE,
  STRAIGHT_CURVE_TYPE,
  POINT_CURVE_TYPE,
  HELIX_CURVE_TYPE,
  HYPERBOLA_CURVE_TYPE,
  UNDEFINED_CURVE_TYPE,
  
    /* Surface types */
  CONE_SURFACE_TYPE,
  PLANE_SURFACE_TYPE,
  SPHERE_SURFACE_TYPE,
  SPLINE_SURFACE_TYPE,
  TORUS_SURFACE_TYPE,
  BEST_FIT_SURFACE_TYPE,
  FACET_SURFACE_TYPE,
  CYLINDER_SURFACE_TYPE, // only currently defined in ACIS Engine
  REVOLUTION_SURFACE_TYPE,  //OCC surface type
  EXTRUSION_SURFACE_TYPE,   //OCC surface type
  OFFSET_SURFACE_TYPE,      //OCC surface type
  UNDEFINED_SURFACE_TYPE,
  
    /* Lump types */
  UNDEFINED_LUMP_TYPE
};

enum ImprintType
{
  NO_IMPRINT=0,
  ONLY_INVOLVED_BODIES,
  INCLUDE_NEIGHBORS,
  TOL_IMPRINT,
  TOL_IMPRINT_INCLUDE_NEIGHBORS
};

/* loops may be the following types */
enum LoopType
{
  /* Unknown loop type */
  LOOP_TYPE_UNKNOWN,
  
  /* The external loop of a surface */
  LOOP_TYPE_EXTERNAL,
  
  /* The loop is a hole */
  LOOP_TYPE_HOLE,
  
  /* The loop is a u or v periodic loop (only applies to periodic surfaces)
   * An example of this is a cylindrical surface with only 2 loops with each
   * loop defining a cap of the cylinder. 
   * If its a u periodic loop, the face is periodic in u */
  LOOP_TYPE_U_PERIODIC,
  LOOP_TYPE_V_PERIODIC
};

struct ModelExportOptions
{
  //this option only relevant to IGES export
  unsigned short export_as_solid : 1;

  //Granite, Acis option
  CubitString logfile_name;
};



#endif

