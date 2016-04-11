/*
 *  ScaledFields.cpp
 *  Created by Jérémy Riousset on 10/25/07.
 */

#include "ScaledFields.h"

/**************************************************************************************/
CMatrix1D ScaledFields::init(double EE_gnd, double zz_gnd, ResGrid dd, SizeGrid NN)
{
	CMatrix1D MM(NN.z);
	int kk_gnd = (int)round(zz_gnd/dd.z);
	for (int kk=0 ; kk<NN.z ; kk++)	MM[kk] = EE_gnd*Scaling::StdAtmosphere((kk+kk_gnd)*dd.z);
	return MM;
}
/**************************************************************************************/