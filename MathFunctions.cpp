/*
 *  MathFunctions.cpp
 *  MR
 *
 *  Created by Jérémy Riousset on 5/21/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "MathFunctions.h"

double SF::ellipke(double k)
{
	/*
	 * ellipke
	 * --------------- 
	 * Returns value of complete elliptic integral of the first kind for  
	 * specified parameter k. Uses Algorithm 55 by John Herndon (see Netlib
	 * for further details).
	 */

	double result, t;
	
	t=1.-k;
	
	result=(((0.032024666*t+0.054544409)*t + 0.097932891)*t + 1.3862944)
	- (((0.010944912*t +0.060118519)*t+0.12475074)*t+0.5)*log(t);
	
	return result;	
}
