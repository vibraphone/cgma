/*
 *
 *
 * Copyright (C) 2004 Sandia Corporation.  Under the terms of Contract DE-AC04-94AL85000
 * with Sandia Corporation, the U.S. Government retains certain rights in this software.
 *
 * This file is part of facetbool--contact via cubit@sandia.gov
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *
 *
 */

#ifndef _FBDEFINES
#define _FBDEFINES

const double EPSILON = 1.e-7;
const double EPSILON2 = 1.e-14;
const int NUMHASHBINS = 101;
const int INTERIOR_VERT = 7;
const int UNKNOWN_VERT = 8;
const int VERTEX_0 = 1;
const int VERTEX_1 = 2;
const int VERTEX_2 = 3;
const int EDGE_0 = 4;
const int EDGE_1 = 5;
const int EDGE_2 = 6;
const int BDRY_EDGE = 9;
const int INTERIOR_EDGE = -999999999;
const int INTERSECTION_EDGE = -55555555;
const int UNKNOWN = 999999999;
const int LEFT = 1;
const int RIGHT = 2;
const int BOTH = 3;
const int UNSET = -1;
const int BDRY_VERT = 777777777;
const int FB_INSIDE = 1;
const int FB_OUTSIDE = 2;
const int FB_SAME = 4;
const int FB_OPPOSITE = 8;
const int FB_ORIENTATION_INSIDE = 1;
const int FB_ORIENTATION_ON = 0;
const int FB_ORIENTATION_OUTSIDE = 2;
const int FB_ORIENTATION_UNDEFINED = -2;
const int FB_ORIENTATION_SAME = 4;
const int FB_ORIENTATION_OPPOSITE = 8;
const int FB_INVALID_BODY = -3;

#endif
