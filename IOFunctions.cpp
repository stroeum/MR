/*
 *  IOFunctions.cpp
 *  Created by Jérémy Riousset on 10/26/07.
 */

#include "IOFunctions.h"

/**************************************************************************************/
/* Read a file with doubles, store in a list of double								  */
/**************************************************************************************/

list<double> IO::read(char * fname)
{
	list<double> LL;
	double xx;
	ifstream inFile;
	
	inFile.open(fname);
	if(!inFile)
		cerr<<"cannot open file";
	//nrerror("cannot open file"); 
	
	while(inFile>>xx)
		LL.push_back(xx);
	inFile.close();
	return LL;
}
/**************************************************************************************/

/**************************************************************************************/
/* write a double, 3 ints, 3 doubles, 4 doubles,									  */
/* a list of doubles, vectors, 1-D matrices											  */
/**************************************************************************************/
void IO::write(double dd, char * fname)
{
	FILE * file = fopen (fname, "w");
	
	if(file) fprintf(file,"%f\n", dd);
	fclose(file);
}

void IO::write(int nn, int mm, int pp,	char * fname)
{
	FILE * file = fopen (fname, "w");
	
	if(file) fprintf(file,"%d\n%d\n%d\n", nn, mm, pp);
	fclose(file);
}

void IO::write(double aa, double bb, double cc,	char * fname)
{
	FILE * file = fopen (fname, "w");
	
	if(file) fprintf(file,"%f\n%f\n%f\n", aa, bb, cc);
	fclose(file);
}

void IO::write(double aa, double bb, double cc, double dd,	char * fname)
{
	FILE * file = fopen (fname, "w");
	
	if(file) fprintf(file,"%f\n%f\n%f\n%f\n", aa, bb, cc, dd);
	fclose(file);
}

void IO::write(list<double>& LL, char * fname)
{
	FILE * file = fopen (fname, "w");
	list<double>::iterator it;
	
	if(file)
		for (it=LL.begin() ; it!=LL.end() ; it++)
			fprintf(file,"%f\n", *it);
	fclose(file);
}

void IO::write(list<CMatrix1D>& LL, char * fname)
{
	FILE * file = fopen (fname, "w");
	list<CMatrix1D>::iterator it;
	
	if(file)
		for (it=LL.begin() ; it!=LL.end() ; it++)
		{
			for (int ii=0 ; ii<it->getNbElem() ; ii++)
				fprintf(file,"%f ",it->getElem(ii));
			fprintf(file,"\n");
		};
	fclose(file);
}

void IO::write(CMatrix1D& MM, char * fname)
{MM.fwrite(fname);}

void IO::write(CMatrix2D& MM, char * fname)
{MM.fwrite(fname);}

void IO::write(CMatrix3D& MM, char * fname)
{MM.fwrite(fname);}
/**************************************************************************************/