/*
 *  MPI_Utils.h
 *  Created by Jérémy Riousset on 11/19/07.
 */

#ifndef MPIFUNCTIONS_H
#define MPIFUNCTIONS_H
#include "Input.h"
#include "MpiInput.h"
class MPI_foo
{
public:	
	static	void				CreateComm(void);
	static	void				CreateGridComm(void);
	static	void				CreateCartComm(void);
	static	void				InitLocalDimensions(void);
	static	void				FreeComm(void);
	
	static	CMatrix2D			Scatter(CMatrix2D	MM);
	static	CMatrix2D			Scatter(CMatrix2D	MM, bool UpdateGhostRows);
	static	Charge				Scatter(Charge		CC);
	static	ConductivityProfile	Scatter(ConductivityProfile	SSigma);
	static	ConductivityProfile	Scatter(ConductivityProfile	SSigma, bool UpdateGhostRows);
	static	CMatrix2D			Gather(	CMatrix2D	local_MM);
	static	void				CreateGhostVector(void);
	static	void				UpdateInterfaceColumns(CMatrix2D& local_MM, SizeGrid local_NN);
};

#endif // MPIFUNCTIONS_H
