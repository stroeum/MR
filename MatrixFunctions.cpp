/*
 *  CMatrix2DUtils.cpp
 *  Created by Jérémy Riousset on 4/23/08.
 */

#include "MatrixFunctions.h"

CMatrix2D M_foo::Average(CMatrix2D MM, SizeGrid NN)
{
	CMatrix2D Avg(NN.r, NN.z);
	
	// First column //
	if(MPI_Var::r_rank == MPI_Var::root)
	{
		Avg(0,0)		= (2*MM(1,0)+MM(0,1))/3;
		Avg(0,NN.z-1)	= (2*MM(1,NN.z-1)+MM(0,NN.z-2))/3;
		
		for(int kk=1 ; kk<NN.z-1 ; kk++)
			Avg(0,kk)	= (2*MM(1,kk)+MM(0,kk-1)+MM(0,kk+1))/4;		
	}
	
	// Body of the domain //
	for(int ii=1 ; ii<NN.r-1 ; ii++)
	{
		Avg(ii,0)		= (MM(ii-1,0)+MM(ii+1,0)+MM(ii,1))/3;
		for(int kk=1 ; kk<NN.z-1 ; kk++)
			Avg(ii,kk)	= (MM(ii-1,kk)+MM(ii+1,kk)+MM(ii,kk-1)+MM(ii,kk+1))/4;
		Avg(ii,NN.z-1)	= (MM(ii-1,NN.z-1)+MM(ii+1,NN.z-1)+MM(ii,NN.z-2))/3;
	}
	
	// Last column //
	if(MPI_Var::r_rank == MPI_Var::n_processes-1)
	{
		Avg(NN.r-1,0)		= (MM(NN.r-2,0)		+MM(NN.r-1,1))/2;
		Avg(NN.r-1,NN.z-1)	= (MM(NN.r-2,NN.z-1)+MM(NN.r-1,NN.z-2))/2;
		
		for(int kk=1 ; kk<NN.z-1 ; kk++)
			Avg(NN.r-1,kk)	= (MM(NN.r-1,kk)+MM(NN.r-1,kk-1)+MM(NN.r-1,kk+1))/3;		
	}
	
	return Avg;
}

CMatrix2D	M_foo::Dr(CMatrix2D MM, ResGrid dd, SizeGrid NN, int ddirection)
{
	// direction must be either one or minus one //
	CMatrix2D DD(NN.r, NN.z);
	
	// First column //
//	if(MPI_Var::r_rank == MPI_Var::root)
		for(int kk=0 ; kk<NN.z ; kk++)	DD(0,kk)	= 0;		
	
	// Body of the domain //
	for(int ii=1 ; ii<NN.r-1 ; ii++)
		for(int kk=0 ; kk<NN.z ; kk++)
			DD(ii,kk)	= ddirection*(MM(ii+1,kk)-MM(ii-1,kk))/2/dd.r;
	
	// Last column //
//	if(MPI_Var::r_rank == MPI_Var::n_processes-1)
		for(int kk=0 ; kk<NN.z ; kk++)
			DD(NN.r-1,kk)	= ddirection*(MM(NN.r-1,kk)-MM(NN.r-2,kk))/dd.r;
	
	return DD;	
}

CMatrix2D	M_foo::Dz(CMatrix2D MM, ResGrid dd, SizeGrid NN, int ddirection)
{
	CMatrix2D DD(NN.r, NN.z);
	int is, ie;
	is = (MPI_Var::r_rank != 0);
	ie = (MPI_Var::r_rank == MPI_Var::n_processes-1)*(NN.r-1) +(MPI_Var::r_rank != MPI_Var::n_processes-1)*(NN.r-2);
	
	// First column //
	for(int ii=is ; ii<=ie ; ii++)
		DD(ii,0)		= ddirection*(MM(ii,1)-MM(ii,0))/dd.z;
	
	// Body of the domain //
	for(int kk=1 ; kk<NN.z-1 ; kk++) for(int ii=is ; ii<=ie ; ii++)
		DD(ii,kk)		= ddirection*(MM(ii,kk+1)-MM(ii,kk-1))/2/dd.z;
	
	// Last column //
	for(int ii=is ; ii<=ie ; ii++)
		DD(ii,NN.z-1)	= ddirection*(MM(ii,NN.z-1)-MM(ii,NN.z-2))/dd.z;
	
	return DD;	
}

CMatrix1D	M_foo::DzAxis(CMatrix2D MM, ResGrid dd, SizeGrid NN, int ddirection)
{
	CMatrix1D DD(NN.z);
	
	DD(0)		= ddirection*(MM(0,1)-MM(0,0))/dd.z;
	for(int kk=1 ; kk<NN.z-1 ; kk++) 
		DD(kk)	= ddirection*(MM(0,kk+1)-MM(0,kk-1))/2/dd.z;
	DD(NN.z-1)	= ddirection*(MM(0,NN.z-1)-MM(0,NN.z-2))/dd.z;
	
	return DD;	
}
