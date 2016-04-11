/*
 *  MPI_Input.cpp
 *  Created by Jérémy Riousset on 11/19/07.
 */

#include "MpiInput.h"
/*MPI VARIABLES DECLARATION********************************************************/

/*CREATE COMM*/
int						MPI_Var::root;												// Root process rank
int						MPI_Var::world_rank;										// Rank in WORLD_COMM
int						MPI_Var::n_processes;										// Number of processes

/*CREATE GRID COMMUNICATOR*/
int						MPI_Var::grid_rank;											// Rank in grid_communicator
int						MPI_Var::n_dims			= 1;								// 3-D scattering
int						MPI_Var::max_dims		= 1;								// Number of dimensions in the communicator
int						MPI_Var::dim_sizes;											// Number of processes in the communicator
int						MPI_Var::coordinates;										// Coordinates of the process
int						MPI_Var::wrap_around;										// Periodicity of the domain
int						MPI_Var::reorder		= 1;								// Let hardware optimize the choice of processes
MPI_Comm				MPI_Var::grid_comm;											// Grid Communicator

/*CREATE CARTESIAN COMMUNICATOR*/
int						MPI_Var::free_coords;										// Allow column index to vary b/w processes of the column communicator
int						MPI_Var::direction		= 0;								// Direction (col or row) to find the neighbors (in 1-D this can only be 0)
int						MPI_Var::displacement	= 1;								// Look for neighbors in +axis (target is above current process, source is below)
int						MPI_Var::prev_r_rank;										// Rank of the source  in z_comm
int						MPI_Var::r_rank;											// Rank of the process in z_comm
int						MPI_Var::next_r_rank;										// Rank of the target  in z_comm
MPI_Comm				MPI_Var::r_comm;											// Communicator through z
int						MPI_Var::is,MPI_Var::ie;									// local r-limits of the domain

/*CREATE LOCAL PLANES*/
MPI_Datatype			MPI_Var::ghost_vector;										// Vector used for data transmission through ghost rows

/*SEND/RECV PARAMETERS*/
int						MPI_Var::tag;												// Tag for point to point communication
MPI_Status				MPI_Var::status;											// Status of point to point communication
MPI_Request				MPI_Var::request;											// Request for non-blocking point to point communication
/**********************************************************************************/
