/*
 *  ResGrid.cpp
 *  Created by Jérémy Riousset on 10/25/07.
 */

#include "ResGrid.h"

/**************************************************************************************/
ResGrid::ResGrid() {r = 1; z = 1; t=1;};			// Default Constructor

ResGrid::ResGrid(SizeDomain L, SizeGrid N)
{ResGrid::init(L,N);};								// Surcharged constructor

bool ResGrid::init(SizeDomain L, SizeGrid N)
{
	r	= L.r/(N.r-1);
	z	= L.z/(N.z-1);
	t	= L.t/(N.t-1);
	rz	= sqrt(pow(r,2) + pow(z,2));
	return true;
}												// Initiation after declaration

ResGrid::~ResGrid(){};
/**************************************************************************************/
