#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<errno.h>
#include<math.h>
#include"clcg4.h"
#include<mpi.h>
#include<pthread.h>

/***************************************************************************/
/* Defines *****************************************************************/
/***************************************************************************/
#define ALIVE 1
#define DEAD  0

#define u_size 20 
#define num_ticks 5 

/***************************************************************************/
/* Global Vars *************************************************************/
/***************************************************************************/

int mpi_myrank; //rank id
int mpi_commsize; //number of total ranks
//MPI Variables
/*MPI_Status status; 
MPI_Request send_request, recv_request;*/
//input arguments
int num_pthreads; //number of pthreads
double threshold; //threshold value
//Run Time
double start_time, end_time; //beginning and end time of program
//Program Variables
int rows_per_rank; //number of rows per rank
int rows_per_thread;
int **my_rows; //allocated number of rows per rank
int *top_ghost_row; //ghost row (top)
int *bottom_ghost_row; //ghost row (bottom)
int pthread_id;

pthread_barrier_t pbarrier;

/***************************************************************************/
/* Function Decs ***********************************************************/
/***************************************************************************/

void run_simulation(); //runs the pthreads simulation
void pthread_id0(); //function for pthread_id0

/***************************************************************************/
/* Function: Main **********************************************************/
/***************************************************************************/

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);
	MPI_Comm_size( MPI_COMM_WORLD, &mpi_commsize);
	MPI_Comm_rank( MPI_COMM_WORLD, &mpi_myrank);

	InitDefault();

	sscanf(argv[1], "%d", &num_pthreads);
	sscanf(argv[2], "%lf", &threshold);

	rows_per_rank = u_size / mpi_commsize;

	//allocate space for rows and columns on each rank
	my_rows = (int **)calloc(rows_per_rank, sizeof(int*));
	for(int i=0; i < rows_per_rank; i++)
	{
		my_rows[i] = (int *)calloc(u_size, sizeof(int));
	}

	//Randomly initialize universe
	for (int i=0; i < rows_per_rank; i++)
	{
		for(int j=0; j < u_size; j++)
		{
			if (GenVal(mpi_myrank*rows_per_rank+i) > threshold)
			{
				my_rows[i][j] = ALIVE;
     			}
			else
			{
				my_rows[i][j] = DEAD;
			}
		}
	}

	//Create pthreads
	pthread_t p_threads[num_pthreads-1];
	pthread_attr_t attr;
	pthread_attr_init (&attr);

	pthread_id=0;
	if(num_pthreads>0)
	  {
	    rows_per_thread = rows_per_rank/num_pthreads;

	    //allocate memory for ghost rows

	    /*
	    top_ghost_row = (int *)calloc(u_size, sizeof(int));
	    bottom_ghost_row = (int *)calloc(u_size, sizeof(int));

	    top_ghost_row = 0;
	    bottom_ghost_row = 0;

	    printf("top ghost row size = %lu, first and last entry: %d, %d\n", sizeof(top_ghost_row), top_ghost_row[0], top_ghost_row[u_size-1]);
	    */
	    
	    pthread_barrier_init(&pbarrier, NULL, num_pthreads);

	    run_simulation();
	  }

	MPI_Barrier( MPI_COMM_WORLD );
	MPI_Finalize();

	//Free allocated memory
	for(int i=0; i < rows_per_rank; i++)
	{
		free(my_rows[i]);
	}
	free(my_rows);

	//End
	
	return 0;
}

void run_simulation()
{
    //play the game
    for(int t = 0; t < num_ticks; t++)
    {
      printf("Rank %d playing the game, at tick %d. Number of ranks %d.\n", mpi_myrank, t, mpi_commsize);
      
      if (pthread_id == 0)
	  {
	    pthread_id0();
	  }
      MPI_Barrier(MPI_COMM_WORLD);

      printf("Rank %d. End of tick for loop. Number of ranks %d.\n", mpi_myrank, mpi_commsize);
    }
}

void pthread_id0()
{

  MPI_Status status; 
  MPI_Request send_request, send_request2, recv_request, recv_request2;

  top_ghost_row = (int *)malloc(sizeof(int)*u_size);
  bottom_ghost_row = (int *)malloc(sizeof(int)*u_size);

  for(int i=0; i<u_size; i++)
    {
      top_ghost_row[i] = 0;
      bottom_ghost_row[i] = 0;
    }

  //printf("I am Rank %d, first and last entry of top ghost row = %d, %d\n", mpi_myrank, top_ghost_row[0], top_ghost_row[u_size-1]);
  

  printf("Rank %d: Number of ranks %d\n", mpi_myrank, mpi_commsize);

  if (mpi_myrank == 0) //first rank
    {
      printf("Rank 0 before send/receive. Number of ranks %d.\n", mpi_commsize);
      MPI_Irecv(&top_ghost_row, u_size, MPI_INT, mpi_commsize-1, 
		1, MPI_COMM_WORLD, &recv_request);
      MPI_Irecv(&bottom_ghost_row, u_size, MPI_INT, mpi_myrank+1,
		2, MPI_COMM_WORLD, &recv_request2);
      printf("Rank 0 after receive requests posted. Number of ranks %d.\n", mpi_commsize);


      int arg = 0;
      //printf("Rank 0: SEND BUFFER:  first row  first entry: %d\n", &my_rows[&arg]);
      
      //send top row to last rank
      MPI_Isend(&my_rows[0], u_size, MPI_INT, mpi_commsize-1,
		2, MPI_COMM_WORLD, &send_request);
      MPI_Wait(&send_request, &status);
      //send last row to rank 1
      MPI_Isend(&my_rows[rows_per_rank-1], u_size, MPI_INT, mpi_myrank+1,
		1, MPI_COMM_WORLD, &send_request2);
      MPI_Wait(&send_request2, &status);

      printf("Rank 0 send requests satisfied. Number of ranks %d\n", mpi_commsize);
      
      //Receive Waits
      int ierr = MPI_Wait(&recv_request, &status);
      printf("TESTING. IERR: %d\n", ierr);
      printf("Rank 0, INBETWEEN receive requests. Number of ranks %d\n", mpi_commsize);
      MPI_Wait(&recv_request2, &status);

      printf("Rank 0 end of send/recv. Number of ranks %d,\n", mpi_commsize);
      printf("top_ghost_row = %d, bottom ghost row = %d\n", top_ghost_row[0], bottom_ghost_row[0]);
    }
  else if (mpi_myrank == mpi_commsize-1) //last rank
    {
      MPI_Irecv(&top_ghost_row, u_size, MPI_INT, mpi_myrank-1,
		1, MPI_COMM_WORLD, &recv_request);
      MPI_Irecv(&bottom_ghost_row, u_size, MPI_INT, 0,
		2, MPI_COMM_WORLD, &recv_request2);
      
      //send top row to previous rank
      MPI_Isend(&my_rows[0], u_size, MPI_INT, mpi_myrank-1,
		2, MPI_COMM_WORLD, &send_request);
      MPI_Wait(&send_request, &status);
      //send last row to first rank
      MPI_Isend(&my_rows[rows_per_rank-1], u_size, MPI_INT, mpi_myrank-1,
		1, MPI_COMM_WORLD, &send_request2);
      MPI_Wait(&send_request2, &status);

      //Receive Waits
      MPI_Wait(&recv_request, &status);
      MPI_Wait(&recv_request2, &status);
    }
  else //middle ranks
    {    
      MPI_Irecv(&top_ghost_row, u_size, MPI_INT, mpi_myrank-1, 
		1, MPI_COMM_WORLD, &recv_request);
      MPI_Irecv(&bottom_ghost_row, u_size, MPI_INT, mpi_myrank+1,
		2, MPI_COMM_WORLD, &recv_request2);
      
      //send top row to previous rank
      MPI_Isend(&my_rows[0], u_size, MPI_INT, mpi_myrank-1,
		2, MPI_COMM_WORLD, &send_request);
      MPI_Wait(&send_request, &status);
      //send bottom row to next rank
      MPI_Isend(&my_rows[rows_per_rank-1], u_size, MPI_INT, mpi_myrank+1,
		1, MPI_COMM_WORLD, &send_request2);
      MPI_Wait(&send_request2, &status);

      //Receive Waits
      MPI_Wait(&recv_request, &status);
      MPI_Wait(&recv_request2, &status);
    }

  printf("Rank %d. End of pthread0 function. Number of Ranks %d.\n", mpi_myrank, mpi_commsize);
}
