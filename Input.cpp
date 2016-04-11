/*
 *  Input.cpp
 *  Created by Jérémy Riousset on 11/19/07.
 */

#include "Input.h"

/*Variables Declaration************************************************************/
SizeGrid				Var::N;														// Number of discretization points
SizeDomain				Var::L;														// Dimensions of the simulation domain
ResGrid					Var::d;														// Lengths of the discretization-grid

const double			Var::epsilon		= eps;									// SOR precision
const int				Var::MaxStep		= mu;									// Allowed maximum number of point per SOR iteration

double					Var::z_gnd;													// _m Altitude of the ground plane
double					Var::Q, Var::Zq, Var::Rq1,Var::Rq2;							// Charge center parameters:
																					// * Charge (Q), 
																					// * Z-coordinates of a charge center (Zq),
																					// * 1st and 2nd geometrical parameters (Rq1, Rq2) 		
Charge					Var::C;														// Charge center	
Cloud					Var::cloud;													// Cloud
CMatrix1D				Var::Ec;													// _V/_m	Scaled E-field initiation threshold

CMatrix2D				Var::phi;													// _V		Total electric potential
CMatrix2D				Var::Er;													// _V/_m	r-component of the electric field
CMatrix2D				Var::Ez;													// _V/_m	z-component of the electric field
CMatrix2D				Var::rhos;													// _C/m3	Source charge density
CMatrix2D				Var::rho;													// _C/m3	Induced charge density
ConductivityProfile		Var::sigma;													// _S/_m, _S/_m2, _S/_m2	atmospheric, r-derivative, z-derivative conductivity profile

/*SLICE 2-D MATRICES*/
SizeGrid				Var::local_N;												// Size of local matrices
CMatrix2D				Var::local_phi;												// Scattered potential matrix
CMatrix2D				Var::local_Er;												// Scattered r-component of the electric field
CMatrix2D				Var::local_Ez;												// Scattered z-component of the electric field
CMatrix2D				Var::local_rho;												// Scattered charge density matrix
CMatrix2D				Var::local_rhos;											// Scattered charge density matrix
ConductivityProfile		Var::local_sigma;											// _S/_m, _S/_m2, _S/_m2	atmospheric, r-derivative, z-derivative conductivity profile		
Charge					Var::local_C;												// Scattered charge


int						Var::flash_type;											// IC, CG, BJ, --, ?? Type of lightning defined by altitude of initiation
CMatrix1D				Var::store_Ez1d;											// _V/_m	1-D matrix for storage of z-component of the electric field at each step
CMatrix1D				Var::store_rho1d;											// _nC/_m3	1-D matrix for storage of charge density at each step 
CMatrix1D				Var::store_rhos1d;											// _nC/_m3	1-D matrix for storage of charge density at each step 
CMatrix2D				Var::store_rho2d;											// _nC/_m3	2-D matrix for storage of charge density at each step 
CMatrix2D				Var::store_rhos2d;											// _nC/_m3	2-D matrix for storage of charge density at each step 
CMatrix2D				Var::store_Q;												// _C		2-D matrix for storage of net charge content of each layer

/*PSOR ALGORITHM*/
SorSolution				Var::SOR;													// Coefficient of SOR solver
/**********************************************************************************/

