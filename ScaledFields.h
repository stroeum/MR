/*
 *  ScaledFields.h
 *  Created by Jérémy Riousset on 10/25/07.
 */

#ifndef SCALED_FIELDS_H
#define SCALED_FIELDS_H

#include "ResGrid.h"
#include "Matrix.h"
#include "AtmScaling.h"

/**************************************************************************************/
class ScaledFields
{
public:
	static CMatrix1D init(double EE_gnd, double zz_gnd, ResGrid dd, SizeGrid NN);		// Initiate after declaration
};
/**************************************************************************************/

#endif SCALED_FIELDS_H