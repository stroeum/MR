/*
 *  MpiUtils.cpp
 *  Created by Jérémy Riousset on 11/19/07.
 */

#include "MpiFunctions.h"

void MPI_foo::CreateComm(void)
{
	MPI_Comm_rank(MPI_COMM_WORLD, &MPI_Var::world_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &MPI_Var::n_processes);
};

void MPI_foo::CreateGridComm(void)
{
	//	int processes_per_dimensions = (int)pow(MPI_Var::n_processes,1/3.);
	MPI_Var::dim_sizes		= MPI_Var::n_processes;										//processes_per_dimensions;
	MPI_Var::wrap_around	= 0;														// The domain is NOT periodic in z
	int tmp_grid_rank;
	
	MPI_Cart_create(MPI_COMM_WORLD, MPI_Var::n_dims, &MPI_Var::dim_sizes, &MPI_Var::wrap_around, MPI_Var::reorder, &MPI_Var::grid_comm);
	MPI_Comm_rank(MPI_Var::grid_comm, &tmp_grid_rank);
	MPI_Cart_coords(MPI_Var::grid_comm, tmp_grid_rank, MPI_Var::max_dims, &MPI_Var::coordinates);
	MPI_Cart_rank(MPI_Var::grid_comm, &MPI_Var::coordinates, &MPI_Var::grid_rank);
}

void MPI_foo::CreateCartComm(void)
{
	MPI_Var::free_coords = 1;
    MPI_Cart_sub(	MPI_Var::grid_comm, &MPI_Var::free_coords, &MPI_Var::r_comm);
	MPI_Comm_rank(	MPI_Var::r_comm,	&MPI_Var::r_rank);
	MPI_Cart_shift(	MPI_Var::r_comm,	MPI_Var::direction, MPI_Var::displacement, &MPI_Var::prev_r_rank, &MPI_Var::next_r_rank);
};

void MPI_foo::InitLocalDimensions(void)
{
	Var::local_N.r = (Var::N.r-2)/MPI_Var::dim_sizes+2;
	Var::local_N.z = Var::N.z;
}

void MPI_foo::FreeComm(void)
{
	MPI_Comm_free( &MPI_Var::grid_comm);
	MPI_Comm_free( &MPI_Var::r_comm);
};

CMatrix2D MPI_foo::Scatter( CMatrix2D MM)
{
	CMatrix2D local_MM;
	if( (Var::N.r-2)%MPI_Var::n_processes == 0 )
	{
		int			source;
		int			dest;
		int			tag		=0;
		MPI_Status	status;
		
		local_MM.init(Var::local_N.r,Var::local_N.z);
		
		//COPY 1st line of rho to 1st line of local_rho on root process//
		if(MPI_Var::r_rank == MPI_Var::root)
		{
			memcpy(&local_MM(0,0), &MM(0,0), Var::local_N.z*sizeof(double));
		};
		
		//SEND last line of rho to last line of local_rho on last process//
		if(MPI_Var::r_rank == MPI_Var::root)
		{
			dest = MPI_Var::dim_sizes-1;
			MPI_Send(	&MM(Var::N.r-1,0), Var::N.z,	MPI_DOUBLE, dest, tag, MPI_Var::r_comm);
		}
		if(MPI_Var::r_rank == MPI_Var::dim_sizes-1)
		{
			source = MPI_Var::root;
			MPI_Recv(	&local_MM(Var::local_N.r-1,0), Var::local_N.z,	MPI_DOUBLE, source, tag, MPI_Var::r_comm, &status);
		}
		
		//SCATTER all intermediate lines//
		MPI_Scatter(&MM(1,0), (Var::local_N.r-2)*Var::local_N.z, MPI_DOUBLE, &local_MM(1,0), (Var::local_N.r-2)*Var::local_N.z, MPI_DOUBLE, MPI_Var::root, MPI_Var::r_comm);
	}
	else
	{
		cout<<"\n***number of processes incompatible with discretization***"<<endl;
		exit(10);
	}
	return local_MM;
};

CMatrix2D	MPI_foo::Scatter(CMatrix2D	MM, bool UpdateGhostRows)
{
	CMatrix2D local_MM;
	local_MM = Scatter(MM);
	UpdateInterfaceColumns(local_MM, Var::local_N);
	return local_MM;
}

Charge MPI_foo::Scatter(Charge CC)
{
	// WE DO NOT CARRY UNECESSARY INFORMATION ABOUT THE CHARGE type TO SAVE COMMUNICATION (BCAST) TIME //
	Charge local_CC;
	local_CC.rho  = Scatter(CC.rho);
	return local_CC;
};

ConductivityProfile	MPI_foo::Scatter(ConductivityProfile	SSigma)
{
	// WE DO NOT CARRY UNECESSARY INFORMATION ABOUT THE CONDUCTIVITY PROFILE TO SAVE COMMUNICATION (BCAST) TIME //
	ConductivityProfile local_SSigma;
	local_SSigma.atm	= Scatter(SSigma.atm);
	local_SSigma.dr		= Scatter(SSigma.dr);
	local_SSigma.dz		= Scatter(SSigma.dz);
	return local_SSigma;
};

ConductivityProfile	MPI_foo::Scatter(ConductivityProfile	SSigma, bool UpdateGhostRows)
{
	// WE DO NOT CARRY UNECESSARY INFORMATION ABOUT THE CONDUCTIVITY PROFILE TO SAVE COMMUNICATION (BCAST) TIME //
	ConductivityProfile local_SSigma;
	local_SSigma.atm	= Scatter(SSigma.atm);
	UpdateInterfaceColumns(local_SSigma.atm, Var::local_N);
	local_SSigma.dr		= Scatter(SSigma.dr);
	UpdateInterfaceColumns(local_SSigma.dr, Var::local_N);
	local_SSigma.dz		= Scatter(SSigma.dz);
	UpdateInterfaceColumns(local_SSigma.dr, Var::local_N);
	return local_SSigma;
};

CMatrix2D MPI_foo::Gather( CMatrix2D local_MM)
{
	CMatrix2D MM;
	if( (Var::N.r-2)%MPI_Var::n_processes == 0 )
	{
		int			source;
		int			dest;
		int			tag		=0;
		MPI_Status	status;
		
		MM.init(Var::N.r,Var::N.z);
		
		//COPY 1st line of local_MM to 1st line of MM on root process//
		if(MPI_Var::r_rank == MPI_Var::root)
			memcpy(&MM(0,0), &local_MM(0,0), Var::local_N.z*sizeof(double));
		MPI_Bcast(&MM(0,0), Var::local_N.z, MPI_DOUBLE, MPI_Var::root, MPI_Var::r_comm);  
		
		//SEND last line of local_MM to last line of MM on last process//
		if(MPI_Var::r_rank == MPI_Var::dim_sizes-1)
		{
			dest = MPI_Var::root;
			MPI_Send(	&local_MM(Var::local_N.r-1,0), Var::N.z,	MPI_DOUBLE, dest, tag, MPI_Var::r_comm);
		}
		if(MPI_Var::r_rank == MPI_Var::root)
		{
			source = MPI_Var::dim_sizes-1;
			MPI_Recv(	&MM(Var::N.r-1,0),		Var::local_N.z,		MPI_DOUBLE, source, tag, MPI_Var::r_comm, &status);
		}
		MPI_Bcast(&MM(Var::N.r-1,0), Var::local_N.z, MPI_DOUBLE, MPI_Var::root, MPI_Var::r_comm);  
		
		//GATHER all intermediate lines//
		MPI_Allgather(&local_MM(1,0), (Var::local_N.r-2)*Var::local_N.z, MPI_DOUBLE, &MM(1,0), (Var::local_N.r-2)*Var::local_N.z, MPI_DOUBLE, MPI_Var::r_comm);
	}
	else
	{
		cout<<"\n***number of processes incompatible with discretization***"<<endl;
		exit(101);
	}
	return MM;
};

void MPI_foo::CreateGhostVector(void)
{
	// Var::local_N.z MUST BE EVEN //
	MPI_Type_vector(Var::local_N.z/2, 1, 2, MPI_DOUBLE, &MPI_Var::ghost_vector);
	MPI_Type_commit(&MPI_Var::ghost_vector);
};

void MPI_foo::UpdateInterfaceColumns(CMatrix2D& local_MM, SizeGrid local_NN)
{
	//SEND/RECV FRONT ROW//
	MPI_Send(&local_MM[1][0],				local_NN.z,	MPI_DOUBLE,	MPI_Var::prev_r_rank,	MPI_Var::tag,	MPI_Var::r_comm );		
	MPI_Recv(&local_MM[local_NN.r-1][0],	local_NN.z,	MPI_DOUBLE,	MPI_Var::next_r_rank,	MPI_Var::tag,	MPI_Var::r_comm, &MPI_Var::status);		
	
	//SEND/RECV BACK ROW//
	MPI_Send(&local_MM[local_NN.r-2][0],	local_NN.z,	MPI_DOUBLE,	MPI_Var::next_r_rank,	MPI_Var::tag,	MPI_Var::r_comm );		
	MPI_Recv(&local_MM[0][0],				local_NN.z,	MPI_DOUBLE,	MPI_Var::prev_r_rank,	MPI_Var::tag,	MPI_Var::r_comm, &MPI_Var::status);		
};
