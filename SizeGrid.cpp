/*
 *  SizeGrid.cpp
 *  Created by Jérémy Riousset on 10/25/07.
 */

#include "SizeGrid.h"

/**************************************************************************************/
SizeGrid::SizeGrid(){r = 1; z = 1;};

SizeGrid::SizeGrid(int rr=1, int zz=1, int tt=1)
{SizeGrid::init(rr,zz,tt);};

double SizeGrid::max()
{return (r >= z)*r + (r < z)*z;};

bool SizeGrid::init(int rr, int zz, int tt)
{
	r = rr; z = zz; t = tt;
	return true;
};

SizeGrid::~SizeGrid(){};
/**************************************************************************************/

