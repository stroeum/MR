/*
 *  ConductivityProfile.h
 *  Created by Jérémy Riousset on 4/23/08.
 */

#ifndef CONDUCTIVITYPROFILE_H
#define CONDUCTIVITYPROFILE_H

#include "Matrix.h"
#include "MatrixFunctions.h"
#include <string>

/**************************************************************************************/
class ConductivityProfile
{
private:
	string type;																		// Type of conductivity profile
	double za, zb;																		// Altitudes of the conductivity regions 
	double ra, rb;																		// Radii	 of the conductivity regions
	double sigma_a, sigma_b;															// Value of the propagation negative threshold at z = z_gnd
	double alpha;																		// Slope of tanh

public:
	CMatrix2D atm;																		// Atmospheric conductivity
	CMatrix2D dr;																		// r-derivative or the atmospheric conductivity
	CMatrix2D dz;																		// z-derivative or the atmospheric conductivity
	
	ConductivityProfile(){};															// Default constructor
	ConductivityProfile(string ttype, double zza,double rra,double ssigma_a, double zzb,double rrb,double ssigma_b, double alpha, ResGrid dd,SizeGrid NN);
																						// Constructor surcharge
	ConductivityProfile(string ttype, double zz_gnd, double zza,double rra,double ssigma_a, double zzb, double aalpha, ResGrid dd,SizeGrid NN);
																						// Constructor surcharge
	CMatrix1D getParams();																// Retrieve private parameters
	bool init(string ttype, double zza,double rra,double ssigma_a, double zzb,double rrb,double ssigma_b, double alpha, ResGrid dd,SizeGrid NN);
																						// Initiate after declaration
	bool init(string ttype, double zz_gnd, double zza,double rra,double ssigma_a, double zzb, double aalpha, ResGrid dd,SizeGrid NN);
																						// Initiate after declaration
	~ConductivityProfile(){};															// Destructor
};
/**************************************************************************************/

#endif //CONDUCTIVITYPROFILE_H