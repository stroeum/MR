/*
 *  Lax.h
 *  Created by Jérémy Riousset on 4/23/08.
 */

#ifndef LAX_H
#define LAX_H

#include "BoundaryConditions.h"
#include "MpiFunctions.h"

class LAX
{
private:
	static CMatrix2D	phi_dz;
	static CMatrix2D	phi_dr;	
	static CMatrix2D	rho_avg;	
	static CMatrix2D	tmp_rho;
	static CMatrix2D	tmp_rhos;
	static CMatrix2D	tmp_phi;
	static CMatrix2D	tmp_phi_dz;
	static CMatrix2D	tmp_phi_dr;	
	static SorSolution	tmp_SOR;
	static Charge		tmp_C;
	
public:
	static void Solve(CMatrix2D& pphi, CMatrix2D& rrhos, CMatrix2D& rrho, 
					  ConductivityProfile& ssigma, 
					  ResGrid dd,SizeGrid NN, bool ccompensation);	
	static void Wendroff(Cloud& ccloud, CMatrix2D& pphi, CMatrix2D& rrhos, CMatrix2D& rrho, 
					  ConductivityProfile& ssigma, 
					  ResGrid dd,SizeGrid NN, bool ccompensation);
};

#endif // LAX_H