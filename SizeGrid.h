/*
 *  SizeGrid.h
 *  Created by Jérémy Riousset on 10/25/07.
 */

#ifndef SIZEGRID_H
#define SIZEGRID_H
#include <iostream>
using namespace std;

/**************************************************************************************/
class SizeGrid
{
public:
	int r,z,t;
	SizeGrid();
	SizeGrid(int rr, int zz, int tt);	
	double max();
	bool init(int rr, int zz, int tt);
	~SizeGrid();
}; // Number of steps in the box
/**************************************************************************************/

#endif // SIZEGRID_H
