/*
 *  Input.h
 *  Created by Jérémy Riousset on 11/19/07.
 */

#ifndef INPUT_H
#define INPUT_H
#include "Cloud.h"
#include "ConductivityProfile.h"
#include "SorSolution.h"
#include "ScaledFields.h"

/*GLOBAL VARIABLES DECLARATION*****************************************************/
class Var
{
public:
	/*GLOBAL VARIABLES*/
	static	SizeGrid				N;												// Number of discretization points
	static	SizeDomain				L;												// Dimensions of the simulation domain
	static	ResGrid					d;												// Lengths of the discretization-grid
	
	static	const double			epsilon;										// SOR precision
	static	const int				MaxStep;										// Allowed maximum number of point per SOR iteration

	static	double					z_gnd;											// _m Altitude of the ground plane
	static	double					Q, Zq, Rq1,Rq2;									// Charge center parameters:
																					// * Charge (Q), 
																					// * Z-coordinates of a charge center (Zq),
																					// * 1st and 2nd geometrical parameters (Rq1, Rq2)
	static	Charge					C;												// Charge center	
	static	Cloud					cloud;											// Cloud
	static	CMatrix1D				Ec;												// _V/_m	Scaled E-field initiation threshold
	
	static	CMatrix2D				phi;											// _V		Total electric potential
	static	CMatrix2D				Er;												// _V/_m	r-component of the electric field
	static	CMatrix2D				Ez;												// _V/_m	z-component of the electric field
	static	CMatrix2D				rho;											// _C/_m3	Induced charge density
	static	CMatrix2D				rhos;											// _C/_m3	Source charge density
	static	ConductivityProfile		sigma;											// _S/_m, _S/_m2, _S/_m2	atmospheric, r-derivative, z-derivative conductivity profile
	static int						flash_type;										//  2-IC, 1-CG, 3-BJ, 0---, -1-?? Type of lightning defined by altitude of initiation
	static CMatrix1D				store_Ez1d;										// _V/_m	1-D matrix for storage of z-component of the electric field at each step
	static CMatrix1D				store_rho1d;									// _nC/_m3	1-D matrix for storage of charge density at each step
	static CMatrix1D				store_rhos1d;									// _nC/_m3	1-D matrix for storage of charge density at each step
	static CMatrix2D				store_rho2d;									// _nC/_m3	2-D matrix for storage of charge density at each step
	static CMatrix2D				store_rhos2d;									// _nC/_m3	2-D matrix for storage of charge density at each step
	static CMatrix2D				store_Q;										// _C		2-D matrix for storage of net charge content of each layer
	
	/*SLICE 2-D MATRICES*/
	static	SizeGrid				local_N;
	static	CMatrix2D				local_phi;
	static	CMatrix2D				local_Er;
	static	CMatrix2D				local_Ez;
	static	CMatrix2D				local_rho;
	static	CMatrix2D				local_rhos;
	static	ConductivityProfile		local_sigma;
	static	Charge					local_C;
	
	/*PSOR ALGORITHM*/
	static	SorSolution				SOR;
};
/**********************************************************************************/

#endif // INPUT_H
