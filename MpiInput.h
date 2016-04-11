/*
 *  MPI_Input.h
 *  Created by Jérémy Riousset on 11/19/07.
 */

#ifndef MPI_INPUT_H
#define MPI_INPUT_H

#undef SEEK_SET
#undef SEEK_CUR
#undef SEEK_END
#include <mpi.h>
#include <iostream>

using namespace std;

/*MPI VARIABLES DECLARATION********************************************************/
class MPI_Var
{
public:
	/*CREATE COMM*/
	static	int						root;
	static	int						world_rank;
	static	int						n_processes;
	
	/*CREATE GRID COMMUNICATOR*/
	static	int						grid_rank;
	static	int						n_dims;
	static	int						dim_sizes;
	static	int						max_dims;
	static	int						coordinates;
	static	int						wrap_around;
	static	int						reorder;
	static	MPI_Comm				grid_comm;
	
	/*CREATE CARTESIAN COMMUNICATOR*/
	static	int						free_coords;
	static	int						direction;
	static	int						displacement;
	static	int						prev_r_rank;
	static	int						     r_rank;
	static	int						next_r_rank;
	static	MPI_Comm				r_comm;
	static	int						is,ie;
	
	/*CREATE LOCAL PLANES*/
	static	MPI_Datatype			ghost_vector;
	
	/*SEND/RECV PARAMETERS*/
	static	int						tag;
	static	MPI_Status				status;
	static	MPI_Request				request;
};
/**********************************************************************************/

#endif // MPI_INPUT_H
