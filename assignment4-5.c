/***************************************************************************/
/* Genevieve Grivas ********************************************************/
/* Adehla Habib ************************************************************/
/* Brandon Thorne **********************************************************/
/***************************************************************************/
/* Assignment04_05 *********************************************************/
/***************************************************************************/

/***************************************************************************/
/* Includes ****************************************************************/
/***************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<errno.h>
#include<math.h>

#include<clcg4.h>

#include<mpi.h>


/***************************************************************************/
/* Defines *****************************************************************/
/***************************************************************************/

#define ALIVE 1
#define DEAD  0

/***************************************************************************/
/* Global Vars *************************************************************/
/***************************************************************************/

int mpi_myrank; //rank id
int mpi_commsize; //number of total ranks
//MPI Send and Recv - not sure if we need
MPI_Status status; 
MPI_Request send_request, recv_request;
int ierr;
int *sendbuff, *recvbuff;
//input arguments
int num_pthreads; //number of pthreads
int threshold; //threshold value

/***************************************************************************/
/* Function Decs ***********************************************************/
/***************************************************************************/

// You define these


/***************************************************************************/
/* Function: Main **********************************************************/
/***************************************************************************/

int main(int argc, char *argv[])
{
//    int i = 0;
    int mpi_myrank;
    int mpi_commsize;
// Example MPI startup and using CLCG4 RNG
    MPI_Init( &argc, &argv);
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_commsize);
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_myrank);
    
// Init 16,384 RNG streams - each rank has an independent stream
    InitDefault();
    
// Note, used the mpi_myrank to select which RNG stream to use.
// You must replace mpi_myrank with the right row being used.
// This just show you how to call the RNG.    
    printf("Rank %d of %d has been started and a first Random Value of %lf\n", 
	   mpi_myrank, mpi_commsize, GenVal(mpi_myrank));
    
    MPI_Barrier( MPI_COMM_WORLD );
    
// Insert your code

//Get number of pthreads and threshold value
    num_pthread = argv[1]; //number of pthreads
    threshold = argv[2];

// END -Perform a barrier and then leave MPI
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize();
    return 0;
}

/***************************************************************************/
/* Other Functions - You write as part of the assignment********************/
/***************************************************************************/
