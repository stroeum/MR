/*
 *  Cloud.h
 *  Created by Jérémy Riousset on 5/11/08.
 */

#ifndef CLOUD_H
#define CLOUD_H
#include "Sources.h"
#include "MatrixFunctions.h"
#include "MpiInput.h"
#include <vector>

class Cloud
{
private:
	double			Emax;																// _V/_m	Max E-field in the cloud
	double			zmax;																// _m		Altitude of Emax
	double			z_screen;															// _m		Upper limit of the screening charge
	double			rhos_ref;															// _C/_m3	Reference charge density of the screening charge
	
public:
	string			type;																// Cloud type
	vector<Charge>	layer;																// Description	of each layer
	vector<double>	I;																	// Loading current densities

	Cloud(){};																			// Default Constructor
	~Cloud(){};																			// Destructor
	
	bool	tripole();																	// Tripolar Cloud
	int		flash(CMatrix2D& pphi, CMatrix1D EEc, ResGrid dd, SizeGrid NN);				// Define which type of lightning can develop in a given cloud, if any. 
	bool	update(int flash, CMatrix2D& local_rrhos,  CMatrix2D& local_rrho, ResGrid dd, SizeGrid local_NN);
																						// Evolution of each layer
	CMatrix2D	LWupdate(ResGrid dd, SizeGrid local_NN);								// Evolution of each layer during Lax-Wendroff algorithm
	void	NetCharge(CMatrix2D& local_rrhos, CMatrix2D& local_rrho, ResGrid dd, SizeGrid local_NN);
																						// Net charge in a given charge volume
	double	LayersCharge(CMatrix2D& local_rrho, int nn, ResGrid dd, SizeGrid local_NN);	// Net charge in the volume of layer nn
	double	ScreenCharge(CMatrix2D& local_rrhos, CMatrix2D& local_rrho, ResGrid dd, SizeGrid local_NN);
																						// Net charge in the screening charge region
	double	ScreenChargeBJ(CMatrix2D& local_rrhos, CMatrix2D& local_rrho, ResGrid dd, SizeGrid local_NN);
																						// Net charge in the screening charge region
	void	LayersReDistribute(CMatrix2D& llocal_rho, int nn, double QQ, ResGrid dd, SizeGrid local_NN);
	void	ScreenReDistribute(CMatrix2D& local_rrho, ResGrid dd, SizeGrid local_NN);
	void	ScreenReDistributeBJ(CMatrix2D& local_rrho, ResGrid dd, SizeGrid local_NN);
};

#endif // CLOUD_H
