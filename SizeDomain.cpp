/*
 *  SizeDomain.cpp
 *  Created by Jérémy Riousset on 10/25/07.
 */

#include "SizeDomain.h"

/**************************************************************************************/
SizeDomain::SizeDomain(){r = 1; z = 1; t = 1;};
SizeDomain::SizeDomain(double rr=1,double zz=1, double tt=1)
{SizeDomain::init(rr,zz,tt);};

bool SizeDomain::init(double rr, double zz, double tt)
{
	r = rr; z = zz; t = tt;
	return true;
};

SizeDomain::~SizeDomain(){};
/**************************************************************************************/
