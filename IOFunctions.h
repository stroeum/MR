/*
 *  IOfunctions.h
 *  Created by Jérémy Riousset on 10/26/07.
 */

#ifndef IOFUNCTIONS_H
#define IOFUNCTIONS_H
#include "Input.h"

class IO
{
public:
	static list<double> read(char *);
	static void write(double,							char *);
	static void write(int,int, int,						char *);
	static void write(double,double, double,			char *);
	static void write(double, double, double, double,	char *);
	static void write(list<double>&,					char *);
	static void write(list<CMatrix1D>&,					char *);
	static void write(CMatrix1D&,						char *);
	static void write(CMatrix2D&,						char *);
	static void write(CMatrix3D&,						char *);
};

#endif // IOFUNCTIONS_H
