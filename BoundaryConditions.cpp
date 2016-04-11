/*
 *  BoundaryConditions.cpp
 *  Created by Jérémy Riousset on 5/21/08.
 */

#include "BoundaryConditions.h"

bool BC::Update(string BBCType, CMatrix2D& llocal_rhos, CMatrix2D& llocal_rho, CMatrix2D& llocal_phi)
{
	CMatrix1D	top(Var::N.r-2);
	CMatrix1D	bot(Var::N.r-2);
	CMatrix1D	side(Var::N.z);
	double		A(0);
	double		B(0);
	double		dq, R1,R2,RA,RB, Z1,ZA,ZB,Zb,Zt, Zg_im,Zi_im, Rp,Zp;
	double		D1,D2,DA,DB,D1g,D2g,DAg,DBg,D1i, M1,M2,MA,MB,M1g,M2g,MAg,MBg,M1i;
	double		dqmax(0), local_dqmax(0);
	double		start, finish, runtime(0);
	int			M(25);																	// Account for M ground images and M ionospheric images
	int			k_Ion    = Var::N.z-1;		
	
	MPI_Barrier(MPI_Var::r_comm);
	start = MPI_Wtime();
	
	// Determine reference charge density //
	for(int ip=MPI_Var::is ; ip<=MPI_Var::ie ; ip++)
	{
		Rp	= (MPI_Var::r_rank * (Var::local_N.r-2) + ip)*Var::d.r;
		for(int kp=0 ; kp<Var::local_N.z ; kp++)
		{
			dq =	(Rp == 0) * M_PI * Var::d.r * Var::d.r * Var::d.z/4 * abs(llocal_rhos(ip,kp)+llocal_rho(ip,kp)) + 
			(Rp != 0) * M_PI * 2 * Rp   * Var::d.r * Var::d.z   * abs(llocal_rhos(ip,kp)+llocal_rho(ip,kp));
			if(dq>=local_dqmax)
				local_dqmax = dq;
		}
	}
	MPI_Allreduce(&local_dqmax, &dqmax, 1, MPI_DOUBLE, MPI_MAX, MPI_Var::r_comm);
	R1 = (Var::N.r-1)*Var::d.r;															//	B-------top-------------+
	RA = 0;																				//							|
	RB = 0;																				//							s
	ZA = 0;																				//							i
	ZB = (Var::local_N.z-1)*Var::d.z;													//							d
	Zt = ZB;																			//							e
	Zb = ZA;																			//							|
																						//	A-------bot-------------+
	for(int ip=MPI_Var::is ; ip<=MPI_Var::ie ; ip++)
	{
		Rp = (MPI_Var::r_rank * (Var::local_N.r-2) + ip)*Var::d.r;
		for(int kp=0 ; kp<Var::local_N.z ; kp++)
		{
			Zp	= kp*Var::d.z;
			dq	=	(Rp == 0) * M_PI * Var::d.r * Var::d.r * Var::d.z/4 * (llocal_rhos(ip,kp)+llocal_rho(ip,kp)) + 
			(Rp != 0) * M_PI * 2 * Rp   * Var::d.r * Var::d.z   * (llocal_rhos(ip,kp)+llocal_rho(ip,kp));
			if(abs(dq)>=0.01*dqmax)
			{
				if(BBCType.compare("open") == 0)										// Open boundary conditions
				{
					for(int kk = 0 ; kk<Var::local_N.z ; kk++)
					{
						Z1	= kk*Var::d.z;
						
						D1  = sqrt(pow(R1+Rp,2) + pow(Z1-Zp,2));
						DA	= sqrt(pow(RA+Rp,2) + pow(ZA-Zp,2));
						DB	= sqrt(pow(RB+Rp,2) + pow(ZB-Zp,2));
						
						if(D1 != 0)
						{
							M1 = 4*R1*Rp/(D1*D1);
							//if(MPI_Var::world_rank==MPI_Var::root) printf("@ %d, R1 = %f D1 = %f M1 = %f \n", MPI_Var::r_rank, R1,D1,M1);
							if(M1!=1) side(kk) +=	1/(4*M_PI*eps0)*dq/D1 /2/M_PI * 4 * SF::ellipke(M1);
						}
						if(DA != 0 && kk==0)
						{
							MA = 4*RA*Rp/(DA*DA);
							if(MA!=1) A += 1/(4*M_PI*eps0)*dq/DA /2/M_PI * 4 * SF::ellipke(MA);
						}
						if(DB != 0 && kk==Var::local_N.z-1)
						{
							MB = 4*RB*Rp/(DB*DB);
							if(MB!=1) B += 1/(4*M_PI*eps0)*dq/DB /2/M_PI * 4 * SF::ellipke(MB);
						}
					}
					for(int ii=1 ; ii<Var::N.r-1 ; ii++)
					{
						R2 = ii*Var::d.r;
						D2 = sqrt(pow(R2+Rp,2) + pow(Zb-Zp,2));
						if(D2 != 0)
						{
							M2 = 4*R2*Rp/(D2*D2);
							if(M2!=1) bot[ii-1]	+= 1/(4*M_PI*eps0)*dq/D2 /2/M_PI * 4 *SF::ellipke(M2);
						}
						D2 = sqrt(pow(R2+Rp,2) + pow(Zt-Zp,2));
						if(D2 != 0)
						{
							M2 = 4*R2*Rp/(D2*D2);
							if(M2!=1) top[ii-1]	+= 1/(4*M_PI*eps0)*dq/D2 /2/M_PI * 4 * SF::ellipke(M2);
						}
					}
				}
				
				if(BBCType.compare("ground") == 0)										// One ground image conditions
				{
					for(int kk = 0 ; kk<Var::local_N.z ; kk++)
					{
						Z1	= kk*Var::d.z;
						
						D1  = sqrt(pow(R1+Rp,2) + pow(Z1-Zp,2));
						DA	= sqrt(pow(RA+Rp,2) + pow(ZA-Zp,2));
						DB	= sqrt(pow(RB+Rp,2) + pow(ZB-Zp,2));
						D1g = sqrt(pow(R1+Rp,2) + pow(Z1+Zp,2));
						DAg	= sqrt(pow(RA+Rp,2) + pow(ZA+Zp,2));
						DBg	= sqrt(pow(RB+Rp,2) + pow(ZB+Zp,2));

						if(D1 != 0)
						{
							M1 = 4*R1*Rp/(D1*D1);
							if(M1!=1) side(kk) +=	1/(4*M_PI*eps0)*dq/D1 /2/M_PI * 4 * SF::ellipke(M1);
						}
						if(D1g != 0)
						{
							M1g = 4*R1*Rp/(D1g*D1g);
							if(M1g!=1) side(kk) -=	1/(4*M_PI*eps0)*dq/D1g /2/M_PI * 4 * SF::ellipke(M1g);
						}
						if(DA != 0 && kk==0)
						{
							MA = 4*RA*Rp/(DA*DA);
							if(MA!=1) A += 1/(4*M_PI*eps0)*dq/DA /2/M_PI * 4 * SF::ellipke(MA);
						}
						if(DAg != 0 && kk==0)
						{
							MAg = 4*RA*Rp/(DAg*DAg);
							if(MAg!=1) A -= 1/(4*M_PI*eps0)*dq/DAg /2/M_PI * 4 * SF::ellipke(MAg);
						}
						if(DB != 0 && kk==Var::local_N.z-1)
						{
							MB = 4*RB*Rp/(DB*DB);
							if(MB!=1) B += 1/(4*M_PI*eps0)*dq/DB /2/M_PI * 4 * SF::ellipke(MB);
						}
						if(DBg != 0 && kk==Var::local_N.z-1)
						{
							MBg = 4*RB*Rp/(DBg*DBg);
							if(MBg!=1) B -= 1/(4*M_PI*eps0)*dq/DBg /2/M_PI * 4 * SF::ellipke(MBg);
						}
					}
					for(int ii=1 ; ii<Var::N.r-1 ; ii++)
					{
						R2  = ii*Var::d.r;
						D2  = sqrt(pow(R2+Rp,2) + pow(Zb-Zp,2));
						D2g = sqrt(pow(R2+Rp,2) + pow(Zb+Zp,2));
						if(D2 != 0)
						{
							M2 = 4*R2*Rp/(D2*D2);
							if(M2!=1) bot[ii-1]	+= 1/(4*M_PI*eps0)*dq/D2 /2/M_PI * 4 *SF::ellipke(M2);
						}
						if(D2g != 0)
						{
							M2g = 4*R2*Rp/(D2g*D2g);
							if(M2g!=1) bot[ii-1] -= 1/(4*M_PI*eps0)*dq/D2g /2/M_PI * 4 *SF::ellipke(M2g);
						}
						D2 = sqrt(pow(R2+Rp,2) + pow(Zt-Zp,2));
						D2g = sqrt(pow(R2+Rp,2) + pow(Zt+Zp,2));
						if(D2 != 0)
						{
							M2 = 4*R2*Rp/(D2*D2);
							if(M2!=1) top[ii-1]	+= 1/(4*M_PI*eps0)*dq/D2 /2/M_PI * 4 * SF::ellipke(M2);
						}
						if(D2g != 0)
						{
							M2g = 4*R2*Rp/(D2g*D2g);
							if(M2g!=1) top[ii-1] -= 1/(4*M_PI*eps0)*dq/D2g /2/M_PI * 4 * SF::ellipke(M2g);
						}
					}
				}
				
				if(BBCType.compare("capacitor") == 0)									// Earth-Ionosphere cavity
				{
					for(int kk = 0 ; kk<Var::local_N.z ; kk++)
					{
						Z1	= kk*Var::d.z;
						
						D1  = sqrt(pow(R1+Rp,2) + pow(Z1-Zp,2));
						DA	= sqrt(pow(RA+Rp,2) + pow(ZA-Zp,2));
						DB	= sqrt(pow(RB+Rp,2) + pow(ZB-Zp,2));
						
						//Ambient Charge direct contribution //
						if(D1 != 0)
						{
							M1 = 4*R1*Rp/(D1*D1);
							if(M1!=1) side(kk) +=	1/(4*M_PI*eps0)*dq/D1 /2/M_PI * 4 * SF::ellipke(M1);
						}
												//Ambient Charge mirror contributions //

						Zg_im = kp*Var::d.z;													// Altitude of ground images				//
						Zi_im = kp*Var::d.z;													// Altitude of iono/electrosphere images	//
						
						for(int mm=1; mm<=M; mm++)
						{
							Zg_im -= ( mm%2*2*kp + (mm-1)%2*2*(k_Ion-kp) )*Var::d.z;
							Zi_im += ( (mm-1)%2*2*kp + mm%2*2*(k_Ion-kp) )*Var::d.z;
							
							D1g = sqrt(pow(R1+Rp,2) + pow(Z1-Zg_im,2));
							D1i = sqrt(pow(R1+Rp,2) + pow(Z1-Zi_im,2));							if(D1g != 0)
							{
								M1g = 4*R1*Rp/(D1g*D1g);
								if(M1g!=1) side(kk) +=	pow(-1.0,mm)*1/(4*M_PI*eps0)*dq/D1g /2/M_PI * 4 * SF::ellipke(M1g);
							}
							if(D1i != 0)
							{
								M1i = 4*R1*Rp/(D1i*D1i);
								if(M1i!=1) side(kk) +=	pow(-1.0,mm)*1/(4*M_PI*eps0)*dq/D1i /2/M_PI * 4 * SF::ellipke(M1i);
							}
						}
						side[0] = 0;
						side[Var::N.z-1] = 0;
					}
				}
			}	
		}
	}
	MPI_Reduce(&A, &llocal_phi[0][0],									1, MPI_DOUBLE, MPI_SUM, MPI_Var::root,			MPI_Var::r_comm);
	MPI_Reduce(&B, &llocal_phi[0][Var::local_N.z-1],					1, MPI_DOUBLE, MPI_SUM, MPI_Var::root,			MPI_Var::r_comm);
	MPI_Reduce(&side[0], &llocal_phi[Var::local_N.r-1][0], Var::local_N.z, MPI_DOUBLE, MPI_SUM, MPI_Var::n_processes-1, MPI_Var::r_comm);
	for(int nn = MPI_Var::root ; nn<MPI_Var::n_processes ; nn++) for(int ii=1 ; ii<Var::local_N.r-1 ; ii++)
	{		
		MPI_Reduce(&bot[nn * (Var::local_N.r-2) + ii-1], &llocal_phi[ii][0],				1, MPI_DOUBLE, MPI_SUM, nn, MPI_Var::r_comm);
		MPI_Reduce(&top[nn * (Var::local_N.r-2) + ii-1], &llocal_phi[ii][Var::local_N.z-1],	1, MPI_DOUBLE, MPI_SUM, nn, MPI_Var::r_comm);
	}

	MPI_Barrier(MPI_Var::r_comm);
	finish	 = MPI_Wtime();
	runtime	+= (finish - start);
	
	//if(MPI_Var::world_rank==MPI_Var::root)	printf("Run time for BC  solver: %f s\n",runtime);
	return true;
}
