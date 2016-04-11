/*
 *  main.cpp
 *  Created by Jérémy Riousset on 11/19/07.
 */

#include "Lax.h"
#include "IOFunctions.h"

void	Store(void);
double	TotalCharge(CMatrix2D& local_rrhos, CMatrix2D& local_rrho, ResGrid dd, SizeGrid local_NN);

int main(int argc, char** argv)
{
	/*DEFINE SIMULATION VARIABLES**************************************************/
	Var::N.r =  130;//98;//
	Var::N.z =  290;//130;//
	Var::N.t =  1000001;
	Var::L.init(64.5e3,72.25e3,(Var::N.t-1)*.0004);
//	Var::L.init(48.5e3,32.25e3,(Var::N.t-1)*.01);
	Var::d.init(Var::L,Var::N);
	
	/******************************************************************************/

	/*START MPI********************************************************************/
	MPI_Init(&argc, &argv);
	/******************************************************************************/
	
	//CREATE COMMUNICATORS AND TOPOLOGIES//
	MPI_foo::CreateComm();
	MPI_foo::CreateGridComm();
	MPI_foo::CreateCartComm();	
	
	MPI_foo::InitLocalDimensions();
	MPI_foo::CreateGhostVector();
	MPI_Var::is = (MPI_Var::r_rank != 0);
	MPI_Var::ie = (MPI_Var::r_rank == MPI_Var::n_processes-1)*(Var::local_N.r-1) + (MPI_Var::r_rank != MPI_Var::n_processes-1)*(Var::local_N.r-2);

	Var::z_gnd	= 0.00e3;
	Var::Ec		= ScaledFields::init(2.16e5, Var::z_gnd, Var::d, Var::N);

	// Cloud Initiation //
	Var::cloud.tripole();
	Var::cloud.I[0] = 3.0;//5.0;
	Var::cloud.I[1] = -.25;//-1.;
	
	Var::Q   = 0*  5; Var::Zq  = 4.50e3; Var::Rq1 = 2.75e3; Var::Rq2 = 1.50e3;
	Var::C.disk(Var::Q, Var::Zq, Var::Rq1,Var::Rq2, Var::d,Var::N);
	Var::cloud.layer[0] = Var::C;
	
	Var::Q   = 0*-55; Var::Zq  = 7.50e3; Var::Rq1 = 4.25e3; Var::Rq2 = 3.00e3;
	Var::C.disk(Var::Q, Var::Zq, Var::Rq1,Var::Rq2, Var::d,Var::N);
	Var::cloud.layer[1] = Var::C;
	
	Var::Q   = 0*55; Var::Zq  = 14.00e3; Var::Rq1 = 4.50e3; Var::Rq2 = 2.75e3;
	Var::C.disk(Var::Q, Var::Zq, Var::Rq1,Var::Rq2, Var::d,Var::N);
	Var::cloud.layer[2] = Var::C;
	
	Var::store_Q.init(2*(int)Var::cloud.layer.size() +3 , Var::N.t); // include time, LP, N, P, SC and total charges, respectively

	if (MPI_Var::world_rank == MPI_Var::root)
	{
		Var::phi.init(Var::N.r,Var::N.z);
		Var::rho.init(Var::N.r,Var::N.z);
		Var::rhos.init(Var::N.r,Var::N.z);
		
		Var::sigma.init("cloud", Var::z_gnd, 14.0e+3, 5.50e+3, 5e-14, 6e3, .75e3, Var::d,Var::N);
		Var::sigma.atm.fwrite("results/sigma.dat");

		Var::rhos = Var::cloud.layer[0].rho + Var::cloud.layer[1].rho + Var::cloud.layer[2].rho;
		
//		Var::Q   =-40; Var::Zq  = 11.00e3; Var::Rq1 = 4.00e3; Var::Rq2 = .50e3;
//		Var::C.disk(Var::Q, Var::Zq, Var::Rq1,Var::Rq2, Var::d,Var::N);
//		Var::rho = Var::C.rho;
	}
	//SCATTER MATRICES//
	Var::store_rhos1d.init(Var::local_N.z);
	Var::store_rho1d.init(Var::local_N.z);
	Var::store_rhos2d.init(Var::local_N.r,Var::local_N.z);
	Var::store_rho2d.init(Var::local_N.r,Var::local_N.z);
	
	Var::local_phi		= MPI_foo::Scatter(Var::phi,1);
	Var::local_rhos		= MPI_foo::Scatter(Var::rhos);
	Var::local_rho		= MPI_foo::Scatter(Var::rho);
	Var::local_sigma	= MPI_foo::Scatter(Var::sigma);

	for(int nn = 0 ; nn<Var::N.t ; nn++) //
	{
		// Calculate potential given the charge distribution //
		Var::local_C.distribution(Var::local_rhos+Var::local_rho, Var::local_N);
		Var::SOR.init(Var::epsilon, Var::MaxStep, Var::d,Var::local_N, Var::local_C);
//		BC::Update("closed", Var::local_rhos, Var::local_rho, Var::local_phi);
		Var::SOR.Solve(Var::d, Var::local_N, Var::local_phi);
		
		// Calculate the charge distribution given the potential //
		LAX::Solve(Var::local_phi, Var::local_rhos, Var::local_rho, Var::local_sigma, Var::d, Var::local_N, 0);
//		LAX::Wendroff(Var::cloud, Var::local_phi, Var::local_rhos, Var::local_rho, Var::local_sigma, Var::d, Var::local_N, 1);
		
		Var::flash_type = Var::cloud.flash(Var::local_phi, Var::Ec, Var::d, Var::local_N);
		
//		if (MPI_Var::world_rank == MPI_Var::root) 
//			cout<<nn<<" "<<Var::flash_type<<endl;
		
		Var::cloud.update(Var::flash_type, Var::local_rhos, Var::local_rho, Var::d, Var::local_N);
		
		if(Var::flash_type == -5) 
		{
			Var::N.t=nn; 
			break;
		}
		
		if(nn%1000 == 0)
		{
			Var::local_Er = M_foo::Dr(Var::local_phi, Var::d, Var::local_N, -1);
			Var::local_Ez = M_foo::Dz(Var::local_phi, Var::d, Var::local_N, -1);
			for(int kk=0 ; kk<=Var::N.z-1 ; kk++)
			{
				for(int ii=MPI_Var::is ; ii<=MPI_Var::ie ; ii++) 
				{
					Var::store_rho2d(ii,kk)		= 1e9*Var::local_rho(ii,kk);
					Var::store_rhos2d(ii,kk)	= 1e9*Var::local_rhos(ii,kk);
				}
			}
			
			//GATHER SCATTERED MATRICES//
			Var::phi	= MPI_foo::Gather(Var::local_phi);
			Var::Er		= MPI_foo::Gather(Var::local_Er);
			Var::Ez		= MPI_foo::Gather(Var::local_Ez);
			Var::rho	= MPI_foo::Gather(Var::store_rho2d);
			Var::rhos	= MPI_foo::Gather(Var::store_rhos2d);
			
			if (MPI_Var::world_rank == MPI_Var::root)
			{
				printf("step %4d.\n", nn);
				char		nname[50];
				
				sprintf(nname,"results/phi%dd%d.dat", 2,nn);
				IO::write(Var::phi,nname);
				
				sprintf(nname,"results/Er%dd%d.dat", 2,nn);
				IO::write(Var::Er,nname);
				
				sprintf(nname,"results/Ez%dd%d.dat", 2,nn);
				IO::write(Var::Ez,nname);
				
				sprintf(nname,"results/rhos%dd%d.dat", 2,nn);
				IO::write(Var::rhos,nname);
				
				sprintf(nname,"results/rho%dd%d.dat", 2,nn);
				IO::write(Var::rho,nname);
			}			
		}
		/*
		if(nn%50 == 0)
		{
			for(int kk=0 ; kk<=Var::N.z-1 ; kk++)
			{
				if (MPI_Var::world_rank == MPI_Var::root)
				{
					Var::store_Ez1d			= M_foo::DzAxis(Var::local_phi, Var::d, Var::local_N, -1);
					Var::store_rho1d(kk)	= 1e9*Var::local_rho(0,kk);
					Var::store_rhos1d(kk)	= 1e9*Var::local_rhos(0,kk);
				}
			}
			
			if (MPI_Var::world_rank == MPI_Var::root)
			{
				char		nname[50];
				sprintf(nname,"results/Ez%dd%d.dat", 1,nn);
				IO::write(Var::store_Ez1d,nname);
				sprintf(nname,"results/rhos%dd%d.dat", 1,nn);
				IO::write(Var::store_rhos1d,nname);
				sprintf(nname,"results/rho%dd%d.dat", 1,nn);
				IO::write(Var::store_rho1d,nname);
			}			
		}
		*/
		Var::store_Q[0][nn] = nn*Var::d.t;
		Var::store_Q[1][nn] = Var::cloud.LayersCharge(Var::local_rhos, 0,				Var::d, Var::local_N);
		Var::store_Q[2][nn] = Var::cloud.LayersCharge(Var::local_rho , 0,				Var::d, Var::local_N);
		Var::store_Q[3][nn] = Var::cloud.LayersCharge(Var::local_rhos, 1,				Var::d, Var::local_N);
		Var::store_Q[4][nn] = Var::cloud.LayersCharge(Var::local_rho , 1,				Var::d, Var::local_N);
		Var::store_Q[5][nn] = Var::cloud.LayersCharge(Var::local_rhos, 2,				Var::d, Var::local_N);
		Var::store_Q[6][nn] = Var::cloud.LayersCharge(Var::local_rho , 2,				Var::d, Var::local_N);
		Var::store_Q[7][nn] = Var::cloud.ScreenCharge(Var::local_rhos, Var::local_rho,	Var::d, Var::local_N);
		Var::store_Q[8][nn] = TotalCharge(Var::local_rhos, Var::local_rho, Var::d, Var::N);
	}
	
	/*STORE EXTRA PARAMETERS*******************************************************/
	if (MPI_Var::world_rank == MPI_Var::root) Store();
	/******************************************************************************/
	
	/*STOP MPI*********************************************************************/
	MPI_foo::FreeComm();
	MPI_Finalize();
	/******************************************************************************/
	
	return 0;
}

void Store(void)
{
	CMatrix1D Conductivity(Var::sigma.getParams());
	
	IO::write(Var::z_gnd,					"results/z_gnd.dat"			);
	IO::write(Var::d.r, Var::d.z, Var::d.t,	"results/d.dat"				);
	IO::write(Var::N.r,Var::N.z,Var::N.t,	"results/N.dat"				);
	IO::write(Conductivity,					"results/Conductivity.dat"	);
	IO::write(Var::Ec,						"results/Ec.dat"			);
	IO::write(Var::store_Q,					"results/Q.dat"				);
}

double	TotalCharge(CMatrix2D& local_rrhos, CMatrix2D& local_rrho, ResGrid dd, SizeGrid local_NN)
{
	double rr,zz,dQ,Q;
	
	dQ = 0;
	for(int  ii=MPI_Var::is ; ii<=MPI_Var::ie ; ii++)
	{
		rr = (MPI_Var::r_rank * (local_NN.r-2) + ii)*dd.r;
		for(int kk=0 ; kk<local_NN.z ; kk++) 
		{
			zz = kk*dd.z;
			dQ += (local_rrhos[ii][kk]+local_rrho[ii][kk])*((ii==0) * M_PI*(dd.r/2)*(dd.r/2)*dd.z + (ii!=0)*2*M_PI*rr*dd.r*dd.z);
		}
	}
	MPI_Reduce(&dQ, &Q, 1, MPI_DOUBLE, MPI_SUM, MPI_Var::root, MPI_Var::r_comm);
	return Q;
};
