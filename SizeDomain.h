/*
 *  SizeDomain.h
 *  Created by Jérémy Riousset on 10/25/07.
 */

#ifndef SIZEDOMAIN_H
#define SIZEDOMAIN_H
#include <iostream>
using namespace std;

/**************************************************************************************/
class SizeDomain
{
public:
	double r,z,t;
	SizeDomain();
	SizeDomain(double rr, double zz, double tt);
	bool init( double rr, double zz, double tt);
	~SizeDomain();
}; // Dimensions of the box
/**************************************************************************************/

#endif // SIZEDOMAIN_H
