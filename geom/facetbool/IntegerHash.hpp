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

#ifndef _HASHBOOL
#define _HASHBOOL
#include <memory.h>

class IntegerHash
{

public:

	IntegerHash(int nBins = 101, int binSizeIncr = 100);
	~IntegerHash();
	void getNumberofBins(int *numBins) const;
	int *getHashBin(int hashValue, int *binSize);
	void addtoHashList(int hashValue, int value);
private:
	int numberofBins;
	int binSizeIncrement;
	int **binArray;
	int *binSize, *maxBinSize, hashIndex;
	void allocateMoreHash(int index);
 
};

#endif
