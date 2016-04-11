/*
 *  BoundaryConditions.h
 *  Created by Jérémy Riousset on 5/21/08.
 */

#ifndef BOUNDARYCONDITIONS_H
#define BOUNDARYCONDITIONS_H

#include "MpiInput.h"
#include "Input.h"
#include "MathFunctions.h"

class BC
{
public:
	static bool Update(string BBCType, CMatrix2D& llocal_rhos, CMatrix2D& llocal_rho, CMatrix2D& llocal_phi);
};

#endif // BOUNDARYCONDITIONS_H
