/*
 *  ConductivityProfile.cpp
 *  Created by Jérémy Riousset on 4/23/08.
 */

#include "ConductivityProfile.h"

ConductivityProfile::ConductivityProfile(string ttype, double zza,double rra,double ssigma_a, double zzb,double rrb,double ssigma_b, double aalpha, ResGrid dd,SizeGrid NN)
{ConductivityProfile::init(ttype,zza,rra,ssigma_a, zzb,rrb,ssigma_b, aalpha, dd,NN);}

ConductivityProfile::ConductivityProfile(string ttype, double zz_gnd, double zza,double rra,double ssigma_a, double zzb, double aalpha, ResGrid dd,SizeGrid NN)
{ConductivityProfile::init(ttype,zz_gnd, zza,rra,ssigma_a, zzb, aalpha, dd,NN);}

CMatrix1D ConductivityProfile::getParams()
{
	CMatrix1D PParams(7);
	PParams[0] = za;
	PParams[1] = ra;
	PParams[2] = sigma_a;
	PParams[3] = zb;
	PParams[4] = rb;
	PParams[5] = sigma_b;
	PParams[6] = alpha;
	return PParams;
}

bool ConductivityProfile::init(string ttype, double zza,double rra,double ssigma_a, double zzb,double rrb,double ssigma_b, double aalpha, ResGrid dd,SizeGrid NN)
{
	type	= ttype;
	za		= zza;
	ra		= rra;
	sigma_a	= ssigma_a;
	zb		= zzb;
	rb		= rrb;
	sigma_b	= ssigma_b;
	alpha	= aalpha;
	
	int		is, ie;
	double	r;
	is = (MPI_Var::r_rank != 0);
	ie = (MPI_Var::r_rank == MPI_Var::n_processes-1)*(NN.r-1) +(MPI_Var::r_rank != MPI_Var::n_processes-1)*(NN.r-2);
	atm.init(NN.r,NN.z);
	
	for(int kk=0 ; kk<NN.z ; kk++) for(int ii=0 ; ii<NN.r ; ii++)
	{
		if(type.compare("tanh") == 0)
		{
			r = sqrt(pow(ii*dd.r,2) + pow(kk*dd.z - za,2));
			atm(ii,kk) = sigma_a+(sigma_b-sigma_a)*(1+tanh((r-rb)/alpha))/2;
		}
		if(type.compare("dirac") == 0)
		{
			if(pow(ii*dd.r,2) + pow(kk*dd.z - za,2) < pow(rrb,2))	atm(ii,kk) = sigma_a;
			else atm(ii,kk) = sigma_b;
		}
	}
	
	dr = M_foo::Dr(atm, dd, NN, 1); 
	dz = M_foo::Dz(atm, dd, NN, 1); 
	
	return true;
}

bool ConductivityProfile::init(string ttype, double zz_gnd, double zza,double rra,double ssigma_a, double zzb, double aalpha, ResGrid dd,SizeGrid NN)
{
	type	= ttype;
	za		= zza;
	ra		= rra;
	sigma_a	= ssigma_a;
	zb		= zzb;
	rb		= rra;
	alpha	= aalpha;
	
	int		is, ie;
	is = (MPI_Var::r_rank != 0);
	ie = (MPI_Var::r_rank == MPI_Var::n_processes-1)*(NN.r-1) +(MPI_Var::r_rank != MPI_Var::n_processes-1)*(NN.r-2);
	atm.init(NN.r,NN.z);
	
	for(int kk=0 ; kk<NN.z ; kk++) for(int ii=0 ; ii<NN.r ; ii++)
	{
		if(type.compare("cloud") ==0)
			atm(ii,kk) = (1*(1-(1-tanh((ii*dd.r-ra)/alpha))/2 * (1-tanh((kk*dd.z-za)/alpha))/2) + 0)*sigma_a*exp((kk*dd.z+zz_gnd)/zb);
	}
	
	dr = M_foo::Dr(atm, dd, NN, 1); 
	dz = M_foo::Dz(atm, dd, NN, 1); 
	
	return true;
}
