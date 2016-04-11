/*
 *  SorSolution.h
 *  Created by Jérémy Riousset on 11/19/07.
 */

#ifndef SORSOLUTION_H
#define SORSOLUTION_H

#include "Sources.h"
#include "MpiInput.h"

/**************************************************************************************/
enum SourceType {ChargeDistribution, PotentialDistribution};
enum PointsColor {black, red};
// Allowed type of sources to use SOR solver
/**************************************************************************************/

/**************************************************************************************/
class SorSolution
{
private:
	SourceType	type;																	// Type of source
	double		epsilon;																// Maximum tolerable error
	int			MaxStep;																// Maximum number of iterations
	CMatrix1D	a,b;																	// Laplace's Equation coefficients (cont.)
	double		c[2], d[2], e[2];														// Laplace's Equation coefficients (cont.)
	CMatrix2D	f;																		// Laplace's Equation coefficients (cont.)
	double		ErrDen;																	// Normalisation of the error
//	int			is, ie, ir;

public:	
	SorSolution(){};																	// Default constructor
	SorSolution(double eepsilon,int MMaxStep, ResGrid dd,SizeGrid NN, Charge& CC);	
	// Constructor surcharge
	~SorSolution(){};																	// Destructor
	void init(double eepsilon,int MMaxStep, ResGrid dd,SizeGrid NN, Charge& CC);
	void Update_GhostRows(CMatrix2D& local_MM, SizeGrid local_NN,PointsColor CC);		// Exchange rows
	void Solve(ResGrid dd, SizeGrid NN, CMatrix2D& pphi);
};

#endif // SORSOLUTION_H
