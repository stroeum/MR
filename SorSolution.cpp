/*
 *  SorSolution.cpp
 *  Created by Jérémy Riousset on 11/19/07.
 */

#include "SorSolution.h"

SorSolution::SorSolution(double eepsilon,int MMaxStep, ResGrid dd,SizeGrid local_NN, Charge& local_CC)
{SorSolution::init(eepsilon,MMaxStep, dd, local_NN, local_CC);};

void SorSolution::init(double eepsilon, int MMaxStep, ResGrid dd, SizeGrid local_N, Charge& local_C)
{
	type	= ChargeDistribution;
	
	double local_ErrDen;

	CMatrix1D rr(local_N.r);
	for (int ii = MPI_Var::is ; ii <= MPI_Var::ie ; ii++)	rr(ii) = (MPI_Var::r_rank * (local_N.r-2) +ii)*dd.r;
	
	e[0] = -2*(1/(dd.z*dd.z) + 2/(dd.r*dd.r));
	e[1] = -2*(1/(dd.z*dd.z) + 1/(dd.r*dd.r));
	
	epsilon	= eepsilon;
	MaxStep	= MMaxStep;
	
	a.init(local_N.r);
	b.init(local_N.r);
	
	b.init(local_N.r);
	if(MPI_Var::r_rank == MPI_Var::root) 
	{
		a[0]	= 0;
		b[0]	= 4/(dd.r*dd.r) / e[0];
	}
	
	c[0]	= 1/(dd.z*dd.z) / e[0];
	d[0]	= c[0];
	
	c[1]	= 1 /(dd.z*dd.z) / e[1];
	d[1]	= c[1];
	
	for (int ii = 1 ; ii<=MPI_Var::ie ; ii++)
	{
		a[ii]	= (1/(dd.r*dd.r) - 1/(rr(ii) * 2*dd.r)) / e[1];
		b[ii]	= (1/(dd.r*dd.r) + 1/(rr(ii) * 2*dd.r)) / e[1];
	}
	
	ErrDen			= 0 ;
	f				= local_C.rho;
	for(int  kk=0 ; kk<local_N.z ; kk++) for(int  ii=MPI_Var::is ; ii<=MPI_Var::ie ; ii++)
	{
		if (ii==0)	f[ii][kk]	/= -eps0*e[0];
		else		f[ii][kk]	/= -eps0*e[1];
		local_ErrDen += f[ii][kk]*f[ii][kk];
	}
	MPI_Allreduce(&local_ErrDen, &ErrDen,	1, MPI_DOUBLE,	MPI_SUM, MPI_COMM_WORLD);
}; // Init SOR

void SorSolution::Update_GhostRows(CMatrix2D& local_MM, SizeGrid local_NN, PointsColor CC)			// Exchange rows
{	
	// local_N.r MUST BE EVEN //
	if(CC == black)
	{
		//SEND/RECV FRONT ROW//
		MPI_Send(&local_MM[local_NN.r-2][1],	1,	MPI_Var::ghost_vector,	MPI_Var::next_r_rank,	MPI_Var::tag,	MPI_Var::r_comm );
		MPI_Recv(&local_MM[0][1],				1,	MPI_Var::ghost_vector,	MPI_Var::prev_r_rank,	MPI_Var::tag,	MPI_Var::r_comm, &MPI_Var::status);
		
		//SEND/RECV BACK ROW//
		MPI_Send(&local_MM[1][0],				1,	MPI_Var::ghost_vector,	MPI_Var::prev_r_rank,	MPI_Var::tag,	MPI_Var::r_comm );
		MPI_Recv(&local_MM[local_NN.r-1][0],	1,	MPI_Var::ghost_vector,	MPI_Var::next_r_rank,	MPI_Var::tag,	MPI_Var::r_comm, &MPI_Var::status);
	}
	if(CC == red)
	{
		//SEND/RECV FRONT ROW//
		MPI_Send(&local_MM[local_NN.r-2][0],	1,	MPI_Var::ghost_vector,	MPI_Var::next_r_rank,	MPI_Var::tag,	MPI_Var::r_comm );
		MPI_Recv(&local_MM[0][0],				1,	MPI_Var::ghost_vector,	MPI_Var::prev_r_rank,	MPI_Var::tag,	MPI_Var::r_comm, &MPI_Var::status);
		
		//SEND/RECV BACK ROW//
		MPI_Send(&local_MM[1][1],				1,	MPI_Var::ghost_vector,	MPI_Var::prev_r_rank,	MPI_Var::tag,	MPI_Var::r_comm );
		MPI_Recv(&local_MM[local_NN.r-1][1],	1,	MPI_Var::ghost_vector,	MPI_Var::next_r_rank,	MPI_Var::tag,	MPI_Var::r_comm, &MPI_Var::status);
	}
	/*
	if(CC == black)
	{
		//SEND/RECV FRONT ROW//
		MPI_Isend(&local_MM[local_NN.r-2][1],	1,	MPI_Var::ghost_vector,	MPI_Var::next_r_rank,	MPI_Var::tag,	MPI_Var::r_comm, &MPI_Var::request);
		MPI_Recv(&local_MM[0][1],				1,	MPI_Var::ghost_vector,	MPI_Var::prev_r_rank,	MPI_Var::tag,	MPI_Var::r_comm, &MPI_Var::status);
		MPI_Barrier(MPI_Var::r_comm);
		
		//SEND/RECV BACK ROW//
		MPI_Isend(&local_MM[1][0],				1,	MPI_Var::ghost_vector,	MPI_Var::prev_r_rank,	MPI_Var::tag,	MPI_Var::r_comm, &MPI_Var::request);
		MPI_Recv(&local_MM[local_NN.r-1][0],	1,	MPI_Var::ghost_vector,	MPI_Var::next_r_rank,	MPI_Var::tag,	MPI_Var::r_comm, &MPI_Var::status);
		MPI_Barrier(MPI_Var::r_comm);
	}
	if(CC == red)
	{
		//SEND/RECV FRONT ROW//
		MPI_Isend(&local_MM[local_NN.r-2][0],	1,	MPI_Var::ghost_vector,	MPI_Var::next_r_rank,	MPI_Var::tag,	MPI_Var::r_comm, &MPI_Var::request);
		MPI_Recv(&local_MM[0][0],				1,	MPI_Var::ghost_vector,	MPI_Var::prev_r_rank,	MPI_Var::tag,	MPI_Var::r_comm, &MPI_Var::status);
		MPI_Barrier(MPI_Var::r_comm);
		
		//SEND/RECV BACK ROW//
		MPI_Isend(&local_MM[1][1],				1,	MPI_Var::ghost_vector,	MPI_Var::prev_r_rank,	MPI_Var::tag,	MPI_Var::r_comm, &MPI_Var::request);
		MPI_Recv(&local_MM[local_NN.r-1][1],	1,	MPI_Var::ghost_vector,	MPI_Var::next_r_rank,	MPI_Var::tag,	MPI_Var::r_comm, &MPI_Var::status);
		MPI_Barrier(MPI_Var::r_comm);
	}
	*/
}; // Update Ghost Rows

void SorSolution::Solve(ResGrid dd, SizeGrid local_NN, CMatrix2D& local_pphi)
{
	SizeGrid NN;
	NN.r = (local_NN.r-2)*MPI_Var::dim_sizes+2;
	NN.z = local_NN.z;

	double start, finish, runtime(0);
	MPI_Barrier(MPI_COMM_WORLD);
	start = MPI_Wtime();

	double		rres			= 0;
	double		rrb				= cos( M_PI / NN.max());
	double		wwb				= 2 / ( 1 + sqrt( 1 - pow(rrb,2) ) );
	double		local_ErrNum	= 0;						// Numerator of the Remaining Error
	double		Err				= 0;						// Remaining Error
	double		ErrNum			= 0;
	int			step			= 0;						// Step Counter

	ErrNum	= ErrDen;
	Err		= ErrNum/ErrDen;
	while(Err>epsilon && step<MaxStep)
	{
		local_ErrNum	= 0;
		
		//CALCULATE RED POINTS//
		Update_GhostRows(local_pphi, local_NN, black);
		for(int kk=1 ; kk<local_NN.z-1 ; kk++) for(int  ii=MPI_Var::is ; ii<local_NN.r-1 ; ii++)
			if( (ii+kk)%2 == 0 )
			{
				if(ii == 0) rres	=									b[0]  * local_pphi[ii+1][kk] + c[0] * local_pphi[ii][kk-1] + d[0] * local_pphi[ii][kk+1] + local_pphi[ii][kk] - f[ii][kk];
				else		rres	= a[ii] * local_pphi[ii-1][kk] +	b[ii] *	local_pphi[ii+1][kk] + c[1] * local_pphi[ii][kk-1] + d[1] * local_pphi[ii][kk+1] + local_pphi[ii][kk] - f[ii][kk];
				local_pphi[ii][kk]	-= wwb*rres;
				local_ErrNum		+= (rres*rres);
			};	 
		
		//CALCULATE BLACK POINTS//
		Update_GhostRows(local_pphi, local_NN, red);
		for(int kk=1 ; kk<local_NN.z-1 ; kk++) for(int  ii=MPI_Var::is ; ii<local_NN.r-1 ; ii++)
			if( (ii+kk)%2 != 0 )
			{
				if(ii == 0) rres	=									b[0]  * local_pphi[ii+1][kk] + c[0] * local_pphi[ii][kk-1] + d[0] * local_pphi[ii][kk+1] + local_pphi[ii][kk] - f[ii][kk];
				else		rres	= a[ii] * local_pphi[ii-1][kk] +	b[ii] *	local_pphi[ii+1][kk] + c[1] * local_pphi[ii][kk-1] + d[1] * local_pphi[ii][kk+1] + local_pphi[ii][kk] - f[ii][kk];
				local_pphi[ii][kk]	-= wwb*rres;
				local_ErrNum		+= rres*rres;
			};
		if((step+1)%25 == 0)
			MPI_Allreduce(&local_ErrNum, &ErrNum, 1,	MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		Err = ErrNum/ErrDen;
//		if(MPI_Var::r_rank == 0)	cout<<"@ step = "<<step<<" ErrNum = "<<ErrNum<<" ErrDen = "<< ErrDen<<" Err = "<<Err<<endl;

		step++;		
	};

	MPI_Barrier(MPI_COMM_WORLD);
	finish	 = MPI_Wtime();
	runtime	+= (finish - start);
	
	if(MPI_Var::r_rank == 0)
		if(step==MaxStep)
		{
			cout<<"*** Allowed computation time exceeded.\n";
			cout<<"*** Maximum allowed iteration step reached.\n";
			cout<<"*** Precision of the result : "<< Err<<".\n";
		}
	
//	if(MPI_Var::world_rank==MPI_Var::root)	
//	{
//		printf("Run time for SOR solver: %f s\n",runtime);
//		printf("steps in SOR solver    : %d\n",step);
//	}
}; // SOR Algorithm

//void SorSolution::Solve(ResGrid dd, SizeGrid local_N, const CMatrix2D& local_Un, CMatrix2D& local_phi)
//{
//	
//	double sstartTime	= MPI_Wtime();
//	double eendTime		= 0;
//	double rrunTime		= 0;
//	
//	double		rres	= 0;
//	double		rrb		= cos( M_PI / local_N.max());
//	double		wwb		= 2 / ( 1 + sqrt( 1 - pow(rrb,2) ) );
//	double		local_ErrNum;
//	double		ErrNum  = 0;						// Numerator of the Remaining Error
//	double		Err		= 0;						// Remaining Error
//	int			step	= 0;						// Step Counter
//	
//	ErrNum			= ErrDen;
//	Err				= ErrNum/ErrDen;
//	
//	while(Err>epsilon && step<MaxStep)
//	{
//		Err				= 0;
//		ErrNum			= 0;
//		local_ErrNum	= 0;
//		
//		//	for (int js = 1 ; js<local_N.y-1 ; js++) for (int is = 1 ; is <local_N.x-1 ; is++)
//		//	{
//		//		rres =   a * local_phi[is-1][js]	+  b * local_phi[is+1][js] + 
//		//				 c * local_phi[is][js-1]	+  d * local_phi[is][js+1] + 
//		//				 e * local_phi[is][js]	-  f[is][js];
//		//		local_phi[is][js]	-= wwb*rres;
//		//		local_ErrNum			+= (rres*rres);
//		//	};
//		//	ErrNum = local_ErrNum;
//		//	Err = ErrNum/ ErrDen;
//		
//		//UPDATE RED ROWS//
//		Update_GhostRows(local_phi,local_N,red);
//		
//		//CALCULATE BLACK POINTS ODD LINES//
//		for (int js = 1 ; js <local_N.y-1 ; js+=2)	for (int is = 1 ; is <local_N.x-1 ; is+=2)
//			if(local_Un[is][js]==0)
//			{	
//				rres =  a * local_phi[is-1][js]	+ b * local_phi[is+1][js] + 
//						c * local_phi[is][js-1]	+ d * local_phi[is][js+1] + 
//						e * local_phi[is][js]		-  f[is][js];
//				local_phi[is][js]	-= wwb*rres;
//				local_ErrNum		+= (rres*rres);
//			};
//		//CALCULATE BLACK POINTS EVEN LINES//
//		for (int js = 2 ; js <local_N.y-1 ; js+=2)	for (int is = 2 ; is <local_N.x-1 ; is+=2)
//			if(local_Un[is][js]==0)
//			{	
//				rres =	a * local_phi[is-1][js]	+ b * local_phi[is+1][js] + 
//						c * local_phi[is][js-1]	+ d * local_phi[is][js+1] + 
//						e * local_phi[is][js]		- f[is][js];
//				local_phi[is][js]	-= wwb*rres;
//				local_ErrNum		+= (rres*rres);
//			};
//		
//		//UPDATE BLACK ROWS//
//		Update_GhostRows(local_phi,local_N,black);
//		
//		//CALCULATE RED POINTS ODD LINES//
//		for (int js = 2 ; js <local_N.y-1 ; js+=2)	for (int is = 1 ; is <local_N.x-1 ; is+=2)
//			if(local_Un[is][js]==0)
//			{	
//				rres =	a * local_phi[is-1][js]	+ b * local_phi[is+1][js] + 
//						c * local_phi[is][js-1]	+ d * local_phi[is][js+1] + 
//						e * local_phi[is][js]		- f[is][js];
//				local_phi[is][js]	-= wwb*rres;
//				local_ErrNum		+= (rres*rres);
//			};
//		//CALCULATE RED POINTS EVEN LINES//
//		for (int js = 1 ; js <local_N.y-1 ; js+=2)	for (int is = 2 ; is <local_N.x-1 ; is+=2)
//			if(local_Un[is][js]==0)
//			{	
//				rres =	a * local_phi[is-1][js]	+ b * local_phi[is+1][js] + 
//						c * local_phi[is][js-1]	+ d * local_phi[is][js+1] + 
//						e * local_phi[is][js]		- f[is][js];
//				local_phi[is][js]	-= wwb*rres;
//				local_ErrNum		+= (rres*rres);
//			};
//		
//		//DERIVE GLOBAL ERROR//
//		MPI_Allreduce(&local_ErrNum, &ErrNum, 1,	MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//		Err		= ErrNum/ ErrDen;
//		if(MPI_Var::x_rank == 0)	cout<<"@ step = "<<step<<" local_ErrNum = "<<local_ErrNum<<" local_ErrDen = "<< ErrDen<<" Err = "<<Err<<endl;
//		
//		step++;
//	};
//	eendTime	= MPI_Wtime();
//	rrunTime	+= (eendTime - sstartTime);
//	if(MPI_Var::x_rank == 0)
//	{
//		if(step==MaxStep)
//		{
//			cout<<"*** Allowed computation time exceeded.\n";
//			cout<<"*** Maximum allowed iteration step reached.\n";
//			cout<<"*** Precision of the result : "<< Err<<".\n";
//		}
//		printf("steps in SOR solver    : %d\n",step);
//		printf("Run time for SOR solver: %fs\n",rrunTime);
//	};
//	
//	// cout<<"Max Steps = "<<MaxStep<<endl;
//	// printf("epsilon: %e\n",eErr);
//	// printf("steps in SOR solver    : %d\n",step);
//	// printf("Run time for SOR solver: %fs\n",(double)rrunTime/100);
//	
//	// UUn.fwrite("results/UUn.dat");	
//}
/**************************************************************************************/
