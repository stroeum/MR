/* Sources.cpp */

#include "Sources.h"
/**************************************************************************************/
Charge::Charge()
{
	Charge::init(0, 0, 0,0);
}	// No initialization of rho and Un //

Charge::Charge(ResGrid dd, SizeGrid NN)
{Charge::reset(dd,NN);}

Charge::Charge(double QQ, double ZZq, double RRq1, double RRq2)
{
	Charge::init(QQ, ZZq, RRq1,RRq2);
}	// No initialization of rho and Un //

bool Charge::init(double QQ, double ZZq, double RRq1, double RRq2)
{
	type= "undefined";
	Q	= QQ;
	Zq	= ZZq;
	Rq1	= RRq1;	Rq2	= RRq2;
	// No initialization of rho and Un //
	return true;
}

bool Charge::reset(ResGrid dd, SizeGrid NN)
{
	Charge::init(0, 0, 0,0);
	rho.init(NN.r,NN.z);
	return true;
}
// Reset all charge attributes

bool Charge::distribution(CMatrix2D rrho, SizeGrid NN)
{
	Charge::init(0, 0, 0,0);
	rho = rrho;
	return true;
}

CMatrix1D Charge::getParams()
{
	CMatrix1D PParams(4);
	PParams[0] = Q;
	PParams[1] = Zq;
	PParams[2] = Rq1;
	PParams[3] = Rq2;
	return PParams;
}

string	Charge::getType()
{return type;}

void	Charge::volume(ResGrid dd, SizeGrid NN)
{
	V	= 0;
	for(int ii=0 ; ii<NN.r-1 ; ii++) for(int kk=0 ; kk<NN.z-1 ; kk++)
	{
		if(type.compare("sphere")==0)
		{
			if ( (ii*dd.r)*(ii*dd.r)+(kk*dd.z-Zq)*(kk*dd.z-Zq) <= Rq1*Rq1 )
				V +=  (ii==0) * M_PI*(dd.r/2)*(dd.r/2)*dd.z + (ii!=0)*2*M_PI*ii*dd.r*dd.r*dd.z;
			// V = 4/3 * M_PI * pow(Rq1,3)
		}
		if(type.compare("disk")==0)
		{
			if ( (ii*dd.r)<=Rq1 && abs(kk*dd.z-Zq) <= Rq2/2 )
				V +=  (ii==0) * M_PI*(dd.r/2)*(dd.r/2)*dd.z + (ii!=0)*2*M_PI*ii*dd.r*dd.r*dd.z;
			// V +=  (ii==1) * 2*M_PI*(dd.r/2)*(dd.r/2)*dd.z + (ii!=1)*2*M_PI*ii*dd.r*dd.r*dd.z;
			// V = M_PI * pow(Rq1,2)*Rq2;
		}
	}
}

void	Charge::density(ResGrid dd, SizeGrid NN)
{
	rho.init(NN.r,NN.z);
	
	if (V!=0) rho0 = Q / V;
	
	for(int ii=0 ; ii<NN.r ; ii++) for(int kk=0 ; kk<NN.z ; kk++)
	{
		if(type.compare("sphere")==0)
			if(sqrt(pow(ii*dd.r,2)+pow(kk*dd.z-Zq,2))<=Rq1)
				rho[ii][kk] = rho0;
		
		if(type.compare("disk")==0)
			if ( (ii*dd.r)<=Rq1 && abs(kk*dd.z-Zq) <= Rq2/2 )
				rho[ii][kk] = rho0;
	}
}

bool	Charge::sphere(double QQ, double ZZq, double RR, ResGrid dd, SizeGrid NN)
{
	Charge::init(QQ, ZZq, RR,RR);
	type	= "sphere";
	Charge::volume(dd,NN);
	Charge::density(dd,NN);
	
	return true;	
}

bool	Charge::disk(double QQ, double ZZq, double RRq1,double RRq2, ResGrid dd, SizeGrid NN)
{
	Charge::init(QQ, ZZq, RRq1,RRq2);
	type	= "disk";
	Charge::volume(dd,NN);
	Charge::density(dd,NN);
	
	return true;	
}

CMatrix1D	Charge::MonopoleAnalyticalSolution(	ResGrid dd, SizeGrid NN)
{
	CMatrix1D pphiAn(NN.z);
	// int kkq = (int)round(Zq/dd.z);
	
	for(int kk=0 ; kk<NN.z ; kk++)
	{
		if(fabs(kk*dd.z-Zq)<=Rq1)
			pphiAn[kk] = -Q/(4*eps0*M_PI)*(pow(kk*dd.z-Zq,2)/(2*pow(Rq1,3)) -3/(2*Rq1));
		else
			pphiAn[kk] = Q/(4*eps0*M_PI*fabs(kk*dd.z-Zq));
	}
	return pphiAn;	
}// Analytical solution: Monopole case

CMatrix1D	Charge::DipoleAnalyticalSolution(	ResGrid dd, SizeGrid NN)
{
	CMatrix1D pphiAn(NN.z);
	//	int kkq = (int)round(Zq/dd.z);
	
	for(int kk=0 ; kk<NN.z ; kk++)
	{
		if(fabs(kk*dd.z-Zq)<=Rq1)
			pphiAn[kk] = -Q/(4*eps0*M_PI)*(pow(kk*dd.z-Zq,2)/(2*pow(Rq1,3)) -3/(2*Rq1)) - Q/(4*eps0*M_PI*(kk*dd.z+Zq));
		else
			pphiAn[kk] = Q/(4*eps0*M_PI*fabs(kk*dd.z-Zq))								  - Q/(4*eps0*M_PI*(kk*dd.z+Zq));
	}
	return pphiAn;	
}// Analytical solution: Dipole case

CMatrix1D	Charge::MultipoleAnalyticalSolution(ResGrid dd, SizeGrid NN)
{
	CMatrix1D pphiAn(NN.z);
	/******************************************************************************/
	/* We neglect the ambient Laplacian field on Earth (100 V/m = 1e-3kV/cm)	  */
	/* The potential due to the ambient Laplacian field VL = 0.					  */
	/* To take into account this potential, VL = 100 * z (in V)					  */
	/******************************************************************************/
	double	Eambient = 0;
	double	VL;
	
	/******************************************************************************/
	/* Define number of images and store their positions						  */
	/******************************************************************************/
	int		M(1000);						// Account for M ground images and M ionospheric images
	double	z_GndImg = 0;				// Altitude of ground images
	double	z_IonImg = 0;				// Altitude of iono/electrosphere images
	double	z_Ion    = (NN.z-1)*dd.z;	// Altitude coordinate of the iono/electrosphere
	
	/******************************************************************************/
	/* Derive potential on the central axis										  */
	/******************************************************************************/
	for(int kk=0 ; kk<NN.z ; kk++)
	{
		VL = Eambient * kk *dd.z;
		pphiAn[kk] = VL;		
	}
	for(int kk=0 ; kk<NN.z ; kk++)
	{
		if(fabs(kk*dd.z-Zq)<=Rq1)
			pphiAn[kk] += -Q/(4*eps0*M_PI)*(pow(kk*dd.z-Zq,2)/(2*pow(Rq1,3)) -3/(2*Rq1));
		else
			pphiAn[kk] += Q/(4*eps0*M_PI*fabs(kk*dd.z-Zq));
		
		z_GndImg = Zq; // Altitude of ground images			    //
		z_IonImg = Zq; // Altitude of iono/electrosphere images //
					   //		cout<<"m = "<<setw(3)<<0<<"; z_Ion = "<<setw(8)<<z_Ion<<"; z_GndImg = "<<setw(8)<<z_GndImg<<"; z_IonImg = "<<setw(8)<<z_IonImg<<endl;
		for(int mm=1; mm<=M; mm++)
		{
			z_GndImg	= z_GndImg - ( mm%2*2*Zq + (mm-1)%2*2*(z_Ion-Zq) );
			z_IonImg	= z_IonImg + ( (mm-1)%2*2*Zq + mm%2*2*(z_Ion-Zq) );
			//			cout<<"m = "<<setw(3)<<mm<<"; z_Ion = "<<setw(8)<<z_Ion<<"; z_GndImg = "<<setw(8)<<z_GndImg<<"; z_IonImg = "<<setw(8)<<z_IonImg<<endl;
			
			pphiAn[kk] += pow(-1.0,mm)*
				(Q/(4*eps0*M_PI*fabs(kk*dd.z-z_GndImg)) +	// Ground Images
				 Q/(4*eps0*M_PI*fabs(kk*dd.z-z_IonImg)) ); // Ionosphere Images
		};
	}
	return pphiAn;
}// Analytical solution: Multipole case

Charge Charge::operator+=(const Charge& CC) 
{
	type = "undefined";
	Q	+= CC.Q;
	Zq	 = 0;
	rho += CC.rho;
	return *this;
} // operator +=

Charge& Charge::operator=(const Charge& CC)
{
	type	= CC.type;	
	Q		= CC.Q;
	V		= CC.V;
	rho0	= CC.rho0;
	Zq		= CC.Zq;
	Rq1		= CC.Rq1;	
	Rq2		= CC.Rq2;	
	rho		= CC.rho;
	return *this;
} // operator +=

Charge Charge::operator+(const Charge& CC) const
{
	Charge result(*this);
	result+=CC;
	return result;
} // operator =

void	Charge::fwrite(char * title)
{
	FILE * file = fopen (title, "a");	
	char * ChargeType = &type[0];
	if(file)
	{
		fprintf(file,"type: %s\n",ChargeType);
		fprintf(file," [Q]        = [%f]\n",Q);
		fprintf(file," [Zq]       = [%f]\n",Zq );		
		fprintf(file," [Rq1, Rq2] = [%f %f]\n",Rq1,Rq2);		
	}
	fclose(file);
} // fwrite


ostream & operator<< (ostream & os, const Charge & C)
{
	return os<<"type: "<<C.type<<"\n [Q]        = ["<<C.Q<<"]\n [Zq]       = ["<<C.Zq<<"]\n [Rq1, Rq2] = ["<<C.Rq1<<" "<<C.Rq2<<"]\n"/*<<"rho: "<<C.rho*/;
}

Charge::~Charge(){}
/**************************************************************************************/

/**************************************************************************************/
Potential::Potential()
{
	EquiPotential	= true;
	Vo				= 0;
	Zc				= 0;
	L				= 0;
	H				= 0;
}

Potential::~Potential(){}
/**************************************************************************************/

void write(Charge& CC,	char * fname)
{CC.fwrite(fname);}