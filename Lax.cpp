/*
 *  Lax.cpp
 *  Created by Jérémy Riousset on 4/23/08.
 */

#include "Lax.h"

CMatrix2D	LAX::phi_dr;
CMatrix2D	LAX::phi_dz;
CMatrix2D	LAX::rho_avg;
CMatrix2D	LAX::tmp_rhos;
CMatrix2D	LAX::tmp_rho;
CMatrix2D	LAX::tmp_phi;
CMatrix2D	LAX::tmp_phi_dr;
CMatrix2D	LAX::tmp_phi_dz;
SorSolution	LAX::tmp_SOR;
Charge		LAX::tmp_C;

void LAX::Solve(CMatrix2D& pphi, CMatrix2D& rrhos, CMatrix2D& rrho, ConductivityProfile& ssigma, ResGrid dd,SizeGrid NN, bool ccompensation)
{
	// Lax W/ and W/O COMPENSATION
	// W/  if COMP == 1
	// W/O if COMP == 0
	// Prefix tmp refers to quantities derived at mid-time step (n+1/2)

	MPI_foo::UpdateInterfaceColumns(pphi, NN);
	phi_dr = M_foo::Dr(pphi,dd,NN,1);
	phi_dz = M_foo::Dz(pphi,dd,NN,1);

	if(ccompensation == 0)
	{
		for(int ii=MPI_Var::is ; ii<=MPI_Var::ie ; ii++) for(int kk=0 ; kk<=NN.z-1 ; kk++)
			rrho(ii,kk) = (ssigma.dr(ii,kk)*phi_dr(ii,kk) + ssigma.dz(ii,kk)*phi_dz(ii,kk))*dd.t - ssigma.atm(ii,kk)/eps0*dd.t*(rrhos(ii,kk) + rrho(ii,kk))	+ rrho(ii,kk);
	}
	else if(ccompensation == 1)
	{
		for(int ii=0 ; ii<=NN.r-1 ; ii++) for(int kk=0 ; kk<=NN.z-1 ; kk++)
			rrho(ii,kk) = (ssigma.dr(ii,kk)*phi_dr(ii,kk) + ssigma.dz(ii,kk)*phi_dz(ii,kk))*dd.t - ssigma.atm(ii,kk)/eps0*dd.t*rrho(ii,kk)					+ rrho(ii,kk);
	}
}

void LAX::Wendroff(Cloud& ccloud, CMatrix2D& pphi, CMatrix2D& rrhos, CMatrix2D& rrho, ConductivityProfile& ssigma, ResGrid dd,SizeGrid NN, bool ccompensation)
{
	// Lax-Wendroff W/ and W/O COMPENSATION
	// W/  if COMP == 1
	// W/O if COMP == 0
	// Prefix tmp refers to quantities derived at mid-time step (n+1/2)
	
	MPI_foo::UpdateInterfaceColumns(pphi, NN);
	phi_dr = M_foo::Dr(pphi,dd,NN,1);
	phi_dz = M_foo::Dz(pphi,dd,NN,1);
	rho_avg = M_foo::Average(rrho,NN);
	
	tmp_rho.init(NN.r,NN.z);
	tmp_rhos.init(NN.r,NN.z);
	tmp_phi.init(NN.r,NN.z);
	
	if(ccompensation == 0)
	{
		for(int ii=MPI_Var::is ; ii<=MPI_Var::ie ; ii++) for(int kk=0 ; kk<=NN.z-1 ; kk++)
			tmp_rho(ii,kk) = (ssigma.dr(ii,kk)*phi_dr(ii,kk) + ssigma.dz(ii,kk)*phi_dz(ii,kk))*dd.t/2 - ssigma.atm(ii,kk)/eps0*dd.t/2*(rrhos(ii,kk) + rrho(ii,kk))	+ rho_avg(ii,kk);
	}
	else if(ccompensation == 1)
	{
		for(int ii=0 ; ii<=NN.r-1 ; ii++) for(int kk=0 ; kk<=NN.z-1 ; kk++)
			tmp_rho(ii,kk) = (ssigma.dr(ii,kk)*phi_dr(ii,kk) + ssigma.dz(ii,kk)*phi_dz(ii,kk))*dd.t/2 - ssigma.atm(ii,kk)/eps0*dd.t/2*rrho(ii,kk)					+ rho_avg(ii,kk);
	}

	if(ccloud.type.empty() == 1) 
		tmp_rhos = rrhos;
	else
		tmp_rhos = ccloud.LWupdate(dd, NN);
	
	tmp_C.distribution(tmp_rhos+tmp_rho, NN);
	tmp_SOR.init(eps, mu, dd,NN, tmp_C);
	BC::Update("capacitor", tmp_rhos, tmp_rho, tmp_phi);
	tmp_SOR.Solve(dd, NN, tmp_phi);
	
	MPI_foo::UpdateInterfaceColumns(tmp_phi, NN);
	tmp_phi_dr = M_foo::Dr(tmp_phi,dd,NN,1);
	tmp_phi_dz = M_foo::Dz(tmp_phi,dd,NN,1);

	if(ccompensation == 0)
	{
		for(int ii=MPI_Var::is ; ii<=MPI_Var::ie ; ii++) for(int kk=0 ; kk<=NN.z-1 ; kk++)
			rrho(ii,kk) = (ssigma.dr(ii,kk)*tmp_phi_dr(ii,kk) + ssigma.dz(ii,kk)*tmp_phi_dz(ii,kk))*dd.t - ssigma.atm(ii,kk)/eps0*dd.t*(rrhos(ii,kk) + rrho(ii,kk))	+ rrho(ii,kk);
	}
	else if(ccompensation == 1)
	{
		for(int ii=0 ; ii<=NN.r-1 ; ii++) for(int kk=0 ; kk<=NN.z-1 ; kk++)
			rrho(ii,kk) = (ssigma.dr(ii,kk)*tmp_phi_dr(ii,kk) + ssigma.dz(ii,kk)*tmp_phi_dz(ii,kk))*dd.t - ssigma.atm(ii,kk)/eps0*dd.t*rrho(ii,kk)					+ rrho(ii,kk);
	}
}
