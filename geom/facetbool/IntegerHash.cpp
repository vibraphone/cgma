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

#include "IntegerHash.hpp"

IntegerHash::IntegerHash(int numBins, int binSizeIncr)
{
int i;

	numberofBins = numBins;
	binSizeIncrement = binSizeIncr;
	binSize = new int[numberofBins];
	maxBinSize = new int[numberofBins];
	binArray = new int *[numberofBins];

	for ( i = 0; i < numberofBins; i++ ) {
		maxBinSize[i] = binSizeIncrement;
		binSize[i] = 0;
		binArray[i] = new int[maxBinSize[i]];
	}
}

IntegerHash::~IntegerHash()
{
int i;

	if ( binSize) delete[] binSize;
	if ( maxBinSize) delete[] maxBinSize;
	for ( i = 0; i < numberofBins; i++ ) {
		if ( binArray[i] ) delete[] binArray[i];
	}
	if ( binArray ) delete [] binArray;	
}

void IntegerHash::getNumberofBins(int *numBins) const
{
  *numBins = numberofBins;
}

int *IntegerHash::getHashBin(int hashValue, int *binnSize)
{
	hashIndex = hashValue%numberofBins;
	*binnSize = binSize[hashIndex];
	return binArray[hashIndex];
}

void IntegerHash::addtoHashList(int hashValue, int value)
{
int i;
	hashIndex = hashValue%numberofBins;
// Is it already there?
	for ( i = 0; i < binSize[hashIndex]; i++ ) {
		if ( binArray[hashIndex][i] == value ) return;
	}
	if ( binSize[hashIndex] > maxBinSize[hashIndex] - 1 ) {
//  Add more memory to this bin.
		allocateMoreHash(hashIndex);
	}

	binArray[hashIndex][binSize[hashIndex]] = value;
	binSize[hashIndex] += 1;

}

void IntegerHash::allocateMoreHash(int index)
{
int *itemp;

	maxBinSize[hashIndex] += binSizeIncrement;
	itemp = new int[maxBinSize[hashIndex]];
	memcpy(itemp,binArray[index],binSize[hashIndex]*sizeof(int));
	delete [] binArray[index];
	binArray[index] = itemp;
}

