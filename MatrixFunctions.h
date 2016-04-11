/*
 *  MatrixFunctions.h
 *  Created by Jérémy Riousset on 4/23/08.
 */

#ifndef MATRIXFUNCTIONS_H
#define MATRIXFUNCTIONS_H
#include "Matrix.h"
#include "MpiInput.h"
#include "ResGrid.h"
#include "SizeGrid.h"

class M_foo
{
public:
	static	CMatrix2D	Average(CMatrix2D MM, SizeGrid NN);
	static	CMatrix2D	Dr(CMatrix2D MM, ResGrid dd, SizeGrid NN, int ddirection);
	static	CMatrix2D	Dz(CMatrix2D MM, ResGrid dd, SizeGrid NN, int ddirection);
	static	CMatrix1D	DzAxis(CMatrix2D MM, ResGrid dd, SizeGrid NN, int ddirection); // z-derivative at r=0
};

#endif // MATRIXFUNCTIONS_H
