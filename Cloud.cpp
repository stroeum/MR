/*
 *  Cloud.cpp
 *  Created by Jérémy Riousset on 5/11/08.
 */

#include "Cloud.h"

bool Cloud::tripole()
{
	type = "tripole";
	I.resize(2);
	layer.resize(3);
	return true;
}

int Cloud::flash(CMatrix2D& pphi, CMatrix1D EEc, ResGrid dd, SizeGrid NN)
{
	int fflashtype;
	
	if(MPI_Var::r_rank == MPI_Var::root)
	{
		CMatrix1D EE = M_foo::DzAxis(pphi,dd,NN,-1);
		Emax = 0;
		for (int kk=0 ; kk<NN.z ; kk++)
			if (abs(EE[kk])-EEc[kk] >= Emax)
			{
				Emax = abs(EE[kk])-EEc[kk];
				zmax = kk*dd.z;
			};
		
		if(Emax>0)																			// Discharge can occur
		{	
			if(type.compare("tripole") == 0)												// Discharge can develop in a tripolar structure
			{
				if(layer[0].Zq <= zmax && zmax <= layer[1].Zq)	fflashtype = 1;				// Cloud-to-Ground
				if(layer[1].Zq <= zmax && zmax <= layer[2].Zq)	fflashtype = 2;				// Intracloud
				if(layer[2].Zq <= zmax && zmax <= 20e3		 )	fflashtype = 3;				// Blue Jet
				if(		  20e3 <= zmax || zmax <= layer[0].Zq)	fflashtype = -1;			// Unknown
				printf("Flash @ z = %f km\n",zmax*1e-3);
			}
		}
		else fflashtype = 0;																// None
	}
	
	MPI_Bcast(&fflashtype,	1,	MPI_INT,	MPI_Var::root,	MPI_Var::r_comm );				// Bcast lightning type to all processes
	return fflashtype;
} // return discharge type

bool	Cloud::update(int fflash, CMatrix2D& local_rrhos,  CMatrix2D& local_rrho, ResGrid dd, SizeGrid local_NN)
{
	double QQ1, QQ2, QQ3;
	
	if(type.compare("tripole") == 0)												// Discharge can develop in a tripolar structure
	{
		if(fflash == 1)																// Cloud-to-Ground
		{
			QQ1 = LayersCharge(local_rrhos, 0, dd, local_NN);
			LayersReDistribute(local_rrhos, 0, QQ1/2, dd, local_NN);
			layer[0].Q = QQ1/2;
			QQ2 = LayersCharge(local_rrhos, 1, dd, local_NN);
			LayersReDistribute(local_rrhos, 1, QQ2/2, dd, local_NN);
			layer[1].Q = QQ2/2;
		}
		if(fflash == 2)																// Intracloud/Gigantic Jet/Bolt-from-the-Blue
		{
			QQ1 = LayersCharge(local_rrhos, 1, dd, local_NN);
			QQ2 = LayersCharge(local_rrhos, 2, dd, local_NN);
			QQ3 = ScreenCharge(local_rrhos, local_rrho, dd, local_NN);
			
			//if(QQ3<2*QQ1+QQ2)
			if(abs(QQ2+QQ3)<=.5*abs(QQ1))											// Gigantic Jet/Bolt-from-the-Blue (with large excess of - charge)
			{
				LayersReDistribute(local_rrhos, 1, QQ1/2, dd, local_NN);
				layer[1].Q = QQ1/2;
				LayersReDistribute(local_rrhos, 2, QQ2/2, dd, local_NN);
				layer[2].Q = QQ2/2;
				ScreenReDistribute(local_rrho, dd, local_NN);
			}
			else if(.5*abs(QQ1)<=abs(QQ2+QQ3) && abs(QQ2+QQ3)<=abs(QQ1))			// Intracloud (with small excess of - charge)
			{
				LayersReDistribute(local_rrhos, 1, QQ1+QQ2/2+QQ3/2, dd, local_NN);
				layer[1].Q = QQ1+QQ2/2+QQ3/2;
				LayersReDistribute(local_rrhos, 2, QQ2-QQ2/2, dd, local_NN);
				layer[2].Q = QQ2-QQ2/2;
				ScreenReDistribute(local_rrho, dd, local_NN);
			}
			else if(abs(QQ1)<=abs(QQ2+QQ3))											// Intracloud (with excess of + charge)
			{
				LayersReDistribute(local_rrhos, 1, QQ1-QQ1/2, dd, local_NN);
				layer[1].Q = QQ1-QQ1/2;
				LayersReDistribute(local_rrhos, 2, QQ2+QQ1/2+QQ3/2, dd, local_NN);
				layer[2].Q = QQ2+QQ1/2+QQ3/2;
				ScreenReDistribute(local_rrho, dd, local_NN);
			}
		}
		if(fflash == 3)																// Blue Jet
		{
			QQ2 = LayersCharge(local_rrhos, 2, dd, local_NN);
			LayersReDistribute(local_rrhos, 2, QQ2-QQ2/2, dd, local_NN);
			layer[2].Q = QQ2-QQ2/2;
			//QQ2 = ScreenChargeBJ(local_rrhos, local_rrho, dd, local_NN);
			ScreenReDistributeBJ(local_rrho, dd, local_NN);
		}
		else //if(fflash == 0)														// No-discharge
		{
			layer[0].Q -= I[1]*dd.t;
			LayersReDistribute(local_rrhos, 0, layer[0].Q, dd, local_NN);
			
			layer[1].Q += (-I[0]+I[1])*dd.t;
			LayersReDistribute(local_rrhos, 1, layer[1].Q, dd, local_NN);

			layer[2].Q += I[0]*dd.t;
			LayersReDistribute(local_rrhos, 2, layer[2].Q, dd, local_NN);
		}			
	}
	return true;
}// Evolution of each layer

CMatrix2D	Cloud::LWupdate(ResGrid dd, SizeGrid local_NN)
{

	CMatrix2D local_rrhos(local_NN.r,local_NN.z);
	if(type.compare("tripole") == 0)												// Discharge can develop in a tripolar structure
	{
		double Q[3];
		Q[0] = layer[0].Q - I[1]*dd.t/2;
		LayersReDistribute(local_rrhos, 0, Q[0] , dd, local_NN);
		
		Q[1] = layer[1].Q + (-I[0]+I[1])*dd.t/2;
		LayersReDistribute(local_rrhos, 1, Q[1], dd, local_NN);
		
		Q[2] = layer[2].Q + I[0]*dd.t/2;
		LayersReDistribute(local_rrhos, 2, Q[2], dd, local_NN);
	}
	return local_rrhos;
}// Evolution of each layer

void	Cloud::NetCharge(CMatrix2D& local_rrhos, CMatrix2D& local_rrho, ResGrid dd, SizeGrid local_NN)
{
	double rr,zz,dQ;
	
	for (int nn = 0 ; nn<(int)layer.size() ; nn++)
	{
		dQ = 0;
		for(int  ii=MPI_Var::is ; ii<=MPI_Var::ie ; ii++)
		{
			rr = (MPI_Var::r_rank * (local_NN.r-2) + ii)*dd.r;
			for(int kk=0 ; kk<local_NN.z ; kk++) 
			{
				zz = kk*dd.z;
				
				if(layer[nn].type.compare("sphere")==0)
					if(sqrt(pow(rr,2)+pow(zz-layer[nn].Zq,2))<=layer[nn].Rq1)
						dQ += (local_rrhos[ii][kk]+local_rrho[ii][kk])*((ii==0) * M_PI*(dd.r/2)*(dd.r/2)*dd.z + (ii!=0)*2*M_PI*rr*dd.r*dd.z);
				
				if(layer[nn].type.compare("disk")==0)
					if ( rr<=layer[nn].Rq1 && abs(zz-layer[nn].Zq) <= layer[nn].Rq2/2 )
						dQ += (local_rrhos[ii][kk]+local_rrho[ii][kk])*((ii==0) * M_PI*(dd.r/2)*(dd.r/2)*dd.z + (ii!=0)*2*M_PI*rr*dd.r*dd.z);
			}

		}
		MPI_Allreduce(&dQ, &layer[nn].Q, 1, MPI_DOUBLE, MPI_SUM, MPI_Var::r_comm);
	}
};

double	Cloud::LayersCharge(CMatrix2D& local_rrho, int nn, ResGrid dd, SizeGrid local_NN)
{
	double rr,zz,dQ,NetQ;
	
	dQ = 0;
	for(int  ii=MPI_Var::is ; ii<=MPI_Var::ie ; ii++)
	{
		rr = (MPI_Var::r_rank * (local_NN.r-2) + ii)*dd.r;
		for(int kk=0 ; kk<local_NN.z ; kk++) 
		{
			zz = kk*dd.z;
			
			if(layer[nn].type.compare("sphere")==0)
				if(sqrt(pow(rr,2)+pow(zz-layer[nn].Zq,2))<=layer[nn].Rq1)
					dQ += local_rrho[ii][kk]*((ii==0) * M_PI*(dd.r/2)*(dd.r/2)*dd.z + (ii!=0)*2*M_PI*rr*dd.r*dd.z);
			
			if(layer[nn].type.compare("disk")==0)
				if ( rr<=layer[nn].Rq1 && abs(zz-layer[nn].Zq) <= layer[nn].Rq2/2 )
					dQ += local_rrho[ii][kk]*((ii==0) * M_PI*(dd.r/2)*(dd.r/2)*dd.z + (ii!=0)*2*M_PI*rr*dd.r*dd.z);
		}
	}
	MPI_Allreduce(&dQ, &NetQ, 1, MPI_DOUBLE, MPI_SUM, MPI_Var::r_comm);

	return NetQ;
};

double	Cloud::ScreenCharge(CMatrix2D& local_rrhos, CMatrix2D& local_rrho, ResGrid dd, SizeGrid local_NN)
{
	double rr,zz,dQ,NetQ;
	if(type.compare("tripole") == 0)												// Discharge can develop in a tripolar structure
	{
		if(MPI_Var::r_rank == MPI_Var::root)
		{
			int kkq = (int)round(layer[2].Zq/dd.z);
			rhos_ref = abs(local_rrhos[0][kkq]);
			for (int kk = kkq ; kk<local_NN.z ; kk++)
				if(abs(local_rrho[0][kk])>=.001*rhos_ref)
					z_screen = kk*dd.z;
		}
		MPI_Bcast(&z_screen, 1,	MPI_DOUBLE,	MPI_Var::root,	MPI_Var::r_comm );
//		MPI_Bcast(&z_screen, 1,	MPI_DOUBLE,	MPI_Var::root,	MPI_Var::r_comm );
		int kks = (int)((layer[2].Zq-layer[2].Rq2/2)/dd.z);
		int kke = (int)((layer[2].Zq+layer[2].Rq2/2)/dd.z);
		
		dQ = 0;
		for(int  ii=MPI_Var::is ; ii<=MPI_Var::ie ; ii++)
		{
			rr = (MPI_Var::r_rank * (local_NN.r-2) + ii)*dd.r;
			for(int kk=kks ; kk<=kke ; kk++)
			{
				zz = kk*dd.z;
				// Estimation of the screening charge to the entire space above the thunderstorm //
				// dQ += local_rrho[ii][kk]*((ii==0) * M_PI*(dd.r/2)*(dd.r/2)*dd.z + (ii!=0)*2*M_PI*rr*dd.r*dd.z);
				// Limit Estimation of the screening charge to the column above the P-layer //
				if(layer[2].type.compare("disk")==0)
					if ( rr<=layer[2].Rq1 && abs(zz-layer[2].Zq) <= layer[2].Rq2/2 )
						dQ += local_rrho[ii][kk]*((ii==0) * M_PI*(dd.r/2)*(dd.r/2)*dd.z + (ii!=0)*2*M_PI*rr*dd.r*dd.z);
			}
		}
		MPI_Allreduce(&dQ, &NetQ, 1, MPI_DOUBLE, MPI_SUM, MPI_Var::r_comm);
		
		return NetQ;
	}
	return 0;
};	// Net charge in the screening charge region

void	Cloud::LayersReDistribute(CMatrix2D& llocal_rho, int nn, double QQ, ResGrid dd, SizeGrid local_NN)
{
//	MPI_Barrier(MPI_Var::r_comm);
	double rr,zz, rrho0;
	rrho0 = QQ/layer[nn].V;
	for(int  ii=MPI_Var::is ; ii<=MPI_Var::ie ; ii++)
	{
		rr = (MPI_Var::r_rank * (local_NN.r-2) + ii)*dd.r;
		for(int kk=0 ; kk<local_NN.z ; kk++) 
		{
			zz = kk*dd.z;
			
			if(layer[nn].type.compare("sphere")==0)
				if(sqrt(pow(rr,2)+pow(zz-layer[nn].Zq,2))<=layer[nn].Rq1)
					llocal_rho[ii][kk] = rrho0;
			
			if(layer[nn].type.compare("disk")==0)
				if ( rr<=layer[nn].Rq1 && abs(zz-layer[nn].Zq) <= layer[nn].Rq2/2 )
					llocal_rho[ii][kk] = rrho0;
		}
	}
};

void	Cloud::ScreenReDistribute(CMatrix2D& local_rrho, ResGrid dd, SizeGrid local_NN)
{
	double rr,zz;
	if(type.compare("tripole") == 0)												// Discharge can develop in a tripolar structure
	{
		int kks = (int)((layer[2].Zq-layer[2].Rq2/2)/dd.z);
		int kke = (int)((layer[2].Zq+layer[2].Rq2/2)/dd.z);
		for(int  ii=MPI_Var::is ; ii<=MPI_Var::ie ; ii++)
		{
			rr = (MPI_Var::r_rank * (local_NN.r-2) + ii)*dd.r;
			for(int kk=kks ; kk<=kke ; kk++)
			{
				zz = kk*dd.z;
				if(layer[2].type.compare("disk")==0)
				{	
					if ( rr<=layer[2].Rq1 && abs(zz-layer[2].Zq) <= layer[2].Rq2/2 )
						local_rrho[ii][kk] /=2;
				}
			}
		}
	}
};

double	Cloud::ScreenChargeBJ(CMatrix2D& local_rrhos, CMatrix2D& local_rrho, ResGrid dd, SizeGrid local_NN)
{
	double rr,zz,dQ,NetQ;
	if(type.compare("tripole") == 0)												// Discharge can develop in a tripolar structure
	{
		if(MPI_Var::r_rank == MPI_Var::root)
		{
			int kkq = (int)round(layer[2].Zq/dd.z);
			rhos_ref = abs(local_rrhos[0][kkq]);
			for (int kk = kkq ; kk<local_NN.z ; kk++)
				if(abs(local_rrho[0][kk])>=.01*rhos_ref)
					z_screen = kk*dd.z;
		}
		MPI_Bcast(&z_screen, 1,	MPI_DOUBLE,	MPI_Var::root,	MPI_Var::r_comm );
		MPI_Bcast(&z_screen, 1,	MPI_DOUBLE,	MPI_Var::root,	MPI_Var::r_comm );
		int kks = (int)((layer[2].Zq+layer[2].Rq2/2)/dd.z);
		int kke = (int)(z_screen/dd.z);
		
		dQ = 0;
		for(int  ii=MPI_Var::is ; ii<=MPI_Var::ie ; ii++)
		{
			rr = (MPI_Var::r_rank * (local_NN.r-2) + ii)*dd.r;
			for(int kk=kks ; kk<=kke ; kk++)
			{
				zz = kk*dd.z;
				
				// Estimation of the screening charge to the entire space above the thunderstorm //
				dQ += local_rrho[ii][kk]*((ii==0) * M_PI*(dd.r/2)*(dd.r/2)*dd.z + (ii!=0)*2*M_PI*rr*dd.r*dd.z);
				/* // Limit Estimation of the screening charge to the column above the P-layer //
				 if(layer[2].type.compare("sphere")==0 || layer[2].type.compare("disk")==0)
				 if ( rr<=layer[2].Rq1)
				 dQ += local_rrho[ii][kk]*((ii==0) * M_PI*(dd.r/2)*(dd.r/2)*dd.z + (ii!=0)*2*M_PI*rr*dd.r*dd.z);
				 */
			}
		}
		MPI_Allreduce(&dQ, &NetQ, 1, MPI_DOUBLE, MPI_SUM, MPI_Var::r_comm);
		
		return NetQ;
	}
	return 0;
};	// Net charge in the screening charge region

void	Cloud::ScreenReDistributeBJ(CMatrix2D& local_rrho, ResGrid dd, SizeGrid local_NN)
{
	double rr,zz;
	if(type.compare("tripole") == 0)												// Discharge can develop in a tripolar structure
	{
		int kks = (int)((layer[2].Zq+layer[2].Rq2/2)/dd.z);
		for(int  ii=MPI_Var::is ; ii<=MPI_Var::ie ; ii++)
		{
			rr = (MPI_Var::r_rank * (local_NN.r-2) + ii)*dd.r;
			for(int kk=kks ; kk<local_NN.z ; kk++)
			{
				zz = kk*dd.z;
				if(layer[2].type.compare("sphere")==0 || layer[2].type.compare("disk")==0)
				{	
					// 	Limit to points which induced charge density exceeds 1% of the source charge of the P-layer //				
					/*
					 if(abs(local_rrho[ii][kk])>=.01*rhos_ref)
					 local_rrho[ii][kk] /=2;
					 */
					
					// 	Limit to points in the entire space above the thunderstorm //				
					local_rrho[ii][kk] /=2;
				}
			}
		}
	}
};

