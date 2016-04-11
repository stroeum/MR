/*
 *  ResGrid.h
 *  Created by Jérémy Riousset on 10/25/07.
 */

#ifndef RESGRID_H
#define RESGRID_H

#include <cmath>
#include "SizeDomain.h"
#include "SizeGrid.h"

/**************************************************************************************/
class ResGrid
{
public:
	double r,z,rz, t;
	ResGrid();									// Default Constructor
	ResGrid(SizeDomain L, SizeGrid N);			// Surcharged constructor
	bool init(SizeDomain L, SizeGrid N);			// Initiation after declaration
	~ResGrid();									// Destructor
};
/**************************************************************************************/

#endif // RESGRID_H
