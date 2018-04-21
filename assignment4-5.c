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
#include"clcg4.h"
#include<mpi.h>
#include<pthread.h>

/***************************************************************************/
/* Defines *****************************************************************/
/***************************************************************************/
#define ALIVE 1
#define DEAD  0

#define u_size 1024 //world size 1024x1024
#define num_ticks 128 //128 ticks

/***************************************************************************/
/* Global Vars *************************************************************/
/***************************************************************************/

int mpi_myrank; //rank id
int mpi_commsize; //number of total ranks
//input arguments
int num_pthreads; //number of pthreads
double threshold; //threshold value
//Run Time
double start_time, end_time; //beginning and end time of program
//Program Variables
int rows_per_rank; //number of rows per rank
int rows_per_thread;
int **my_rows; //allocated number of rows per rank
int **ghost_rows; //parallel universe for each rank
int *top_ghost_row; //ghost row (top)
int *bottom_ghost_row; //ghost row (bottom)
pthread_barrier_t pbarrier;

/***************************************************************************/
/* Function Decs ***********************************************************/
/***************************************************************************/

void *run_simulation(void *); //runs the pthreads simulation
void sim_tick(int starting_row); //simulates game for single time "tick" for given rows
void pthread_id0(void); //function for pthread_id0

/***************************************************************************/
/* Function: Main **********************************************************/
/***************************************************************************/

int main(int argc, char *argv[])
{
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
    
    //Get number of pthreads and threshold value
    sscanf(argv[1], "%d", &num_pthreads); //number of pthreads
    sscanf(argv[2], "%lf", &threshold); //threshold value
    
    //Start time of program
    if (mpi_myrank == 0)
    {
        start_time = MPI_Wtime();
    }   
    //determine number of rows for each rank and each thread
    rows_per_rank = u_size/mpi_commsize; //number of rows per rank
    //allocate space for rows and columns on each rank
    my_rows = (int **)calloc(rows_per_rank, sizeof(int*));
	for (int i=0; i<rows_per_rank; i++)
	{
		my_rows[i] = (int *)calloc(u_size, sizeof(int));
	}
	//allocate ghost universe with in each rank

    ghost_rows = (int **)calloc(rows_per_rank, sizeof(int*));
	for (int i=0; i<rows_per_rank; i++)
	{
		ghost_rows[i] = (int *)calloc(u_size, sizeof(int));
	}
	
   //allocate space for ghost rows
   top_ghost_row = (int *)malloc(sizeof(int)*u_size);
   bottom_ghost_row = (int *)malloc(sizeof(int)*u_size);
    
   //Randomly initialize universe
    for (int i=0; i<rows_per_rank; i++)
    {
        for (int j=0; j < u_size; j++)
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

    if (mpi_myrank == 0){ printf("Done Making Universe.\n"); }
    
    //Create Pthreads
    pthread_t p_threads[num_pthreads-1];
    pthread_attr_t attr;
    pthread_attr_init (&attr);
    //2 or more pthreads
    if (num_pthreads > 1)
    {
		rows_per_thread = rows_per_rank/num_pthreads; //number of rows per thread
		printf("Number of rows per thread %d\n", rows_per_thread);
		//Create barrier to sync ticks in pthreads
		pthread_barrier_init(&pbarrier, NULL, num_pthreads);
		
		//pthreads 1 to num_pthreads create and run simulation function
		int i;
		for(i=1; i < num_pthreads; i++)
		{ 
			int *pthread_id_ptr = malloc(sizeof(int));
			*pthread_id_ptr = i;
			pthread_create(&p_threads[i-1], &attr, run_simulation, (void *)pthread_id_ptr);
		}
    }
    else
	{
		rows_per_thread = rows_per_rank;
	}
	
    int *pthread_id_ptr=malloc(sizeof(int));
	*pthread_id_ptr = 0;
	printf("Number of rows per thread %d\n", rows_per_thread);
	run_simulation(pthread_id_ptr);     

    //join threads to make sure all pthreads in rank are done
    for(int i = 0; i < num_pthreads - 1; i++)
	{
		pthread_join(p_threads[i], NULL);
	}
	//printf("Before join\n");
    MPI_Barrier(MPI_COMM_WORLD);

    //End time of program
    if (mpi_myrank == 0)
    {
        end_time = MPI_Wtime();
    }
    //End MPI
    MPI_Finalize(); 
 
    //Free allocated memory
    free(top_ghost_row);
    free(bottom_ghost_row);
    for(int i=0; i < rows_per_rank; i++)
    {
        free(my_rows[i]);
		free(ghost_rows[i]);
    }
    free(my_rows);
	free(ghost_rows);
    
    return 0;
}

/***************************************************************************/
/* Other Functions - You write as part of the assignment********************/
/***************************************************************************/

//Funtion ran in each pthread
void *run_simulation(void *void_pthread_id_ptr)
{
	int pthread_id;
	pthread_id = *((int *)void_pthread_id_ptr);
	free(void_pthread_id_ptr);
	//play the game
	for(int t = 0; t < num_ticks; t++)
	{		
		//(i.) Exchange row data with different MPI ranks for the ghost rows (the rows shared between MPI ranks).
		if (pthread_id == 0)
		{
			pthread_id0();
			printf("Done MPI receive and send\n");
		}
		//If there are 0 or 1 pthreads
		if ( (num_pthreads == 0) | (num_pthreads == 1) )
		{
			sim_tick(0);
			printf("Done with sim_tick\n");
			for (int i=0; i<rows_per_thread; i++)
			{
				for (int j=0; j<u_size; j++)
				{
					my_rows[i][j] = ghost_rows[i][j];
				}
			}
		}
      //multiple pthreads
		else 
		{
	  //(ii.) Read all the statuses of the neighbors depending upon the number of live/dead neighbors for the current tick
	  //and set the appropiate live/dead count.        
			int starting_row = pthread_id * rows_per_thread;
			printf("Rank %d, thread %d has starting row at %d \n", mpi_myrank, pthread_id, starting_row);
			sim_tick(starting_row);
			printf("sim tick finished successfully \n");
	  
			//(iii.) Finally, update the cell status
			for (int i=starting_row; i<rows_per_thread; i++)
			{
				for (int j=0; j<u_size; j++)
				{
					my_rows[i][j] = ghost_rows[i][j];
				}
			}
		}
		//Use pthread barrier to make sure each thread finishes a tick before moving on
		if(num_pthreads > 1) 
		{
			pthread_barrier_wait(&pbarrier);
		}
    }
    return NULL;
}
//runs a simulation for a tick for given rows
void sim_tick(int starting_row)
{
	//For each row in the table
	for(int i = starting_row; i < rows_per_thread; i++)
    {    
		for(int curr_col = 0; curr_col < u_size; curr_col++)
        {
			int num_alive_neighbors = 0; //keep count of live neighbors  
			//set direction
			int left = (curr_col - 1 + u_size) % u_size;
			int right = (curr_col + 1) % u_size;
			int up = i-1;
			int down = i+1;
			
			//count number of alive neighbors
			if (i == 0) //top thread row - use top ghost row
			{	
				num_alive_neighbors = my_rows[i][left] + my_rows[i][right] + top_ghost_row[curr_col] + my_rows[down][curr_col] +
										top_ghost_row[left] + top_ghost_row[right] + my_rows[down][left] + my_rows[down][right];
				
			}
			else if (i == starting_row + rows_per_thread-1) //bottom thread row - use bottom ghost row
			{
				num_alive_neighbors = my_rows[i][left] + my_rows[i][right] + my_rows[up][curr_col] + bottom_ghost_row[curr_col] +
										my_rows[up][left] + my_rows[up][right]+ bottom_ghost_row[left] + bottom_ghost_row[right];
						
			}
			else //middle thread row
			{
				num_alive_neighbors = my_rows[i][left] + my_rows[i][right] + my_rows[up][curr_col] + my_rows[down][curr_col] +
										my_rows[up][left] + my_rows[up][right] + my_rows[down][left] + my_rows[down][right];
				
			}
			//if random value is greater than threshold - play the game
			if (GenVal(mpi_myrank*i+curr_col) > threshold)
			{
				//1- Any live cell with fewer than two live neighbors dies, as if caused by under-population.
				if (num_alive_neighbors < 2)
				{
					ghost_rows[i][curr_col] = DEAD;
				}
				//2- Any live cell with two or three live neighbors lives on to the next generation
				//4- Any dead cell with exactly three live neighbors becomes a live cell, as if by reproduction
				else if(num_alive_neighbors < 4)
				{
					if((num_alive_neighbors == 3) & (my_rows[i][curr_col]==DEAD))
					{
						ghost_rows[i][curr_col] = ALIVE; 
					}
					else
					{
						ghost_rows[i][curr_col] = my_rows[i][curr_col];
					}
				}
				//3- Any live cell with more than three live neighbors dies, as if by over-population
				else if(num_alive_neighbors > 3)
				{
					ghost_rows[i][curr_col] = DEAD;
				}
			}
			//else, randomly select alive or dead
			else
			{
				if (GenVal(mpi_myrank*i+curr_col) > 0.5 )
				{ 
					ghost_rows[i][curr_col] = ALIVE;
				}
				else
				{
					ghost_rows[i][curr_col] = DEAD;
				}
			}
		}
	} 
} 

void pthread_id0()
	{
	MPI_Status status; 
	MPI_Request send_request, send_request2, recv_request, recv_request2;

	if (mpi_myrank == 0) //first rank
	{
		MPI_Irecv(top_ghost_row, u_size, MPI_INT, mpi_commsize-1, 
			1, MPI_COMM_WORLD, &recv_request);
		MPI_Irecv(bottom_ghost_row, u_size, MPI_INT, mpi_myrank+1,
			2, MPI_COMM_WORLD, &recv_request2);
		//send top row to last rank
		MPI_Isend(my_rows[0], u_size, MPI_INT, mpi_commsize-1,
			2, MPI_COMM_WORLD, &send_request);
		MPI_Wait(&send_request, &status);
		//send last row to rank 1
		MPI_Isend(my_rows[rows_per_rank-1], u_size, MPI_INT, mpi_myrank+1,
			1, MPI_COMM_WORLD, &send_request2);
		MPI_Wait(&send_request2, &status);
		//Receive Waits
		MPI_Wait(&recv_request, &status);
		MPI_Wait(&recv_request2, &status);
	}
	else if (mpi_myrank == mpi_commsize-1) //last rank
	{
		MPI_Irecv(top_ghost_row, u_size, MPI_INT, mpi_myrank-1,
			1, MPI_COMM_WORLD, &recv_request);
		MPI_Irecv(bottom_ghost_row, u_size, MPI_INT, 0,
			2, MPI_COMM_WORLD, &recv_request2);
		
		//send top row to previous rank
		MPI_Isend(my_rows[0], u_size, MPI_INT, mpi_myrank-1,
			2, MPI_COMM_WORLD, &send_request);
		MPI_Wait(&send_request, &status);
		//send last row to first rank
		MPI_Isend(my_rows[rows_per_rank-1], u_size, MPI_INT, 0,
			1, MPI_COMM_WORLD, &send_request2);
		MPI_Wait(&send_request2, &status);
		//Receive Waits
		MPI_Wait(&recv_request, &status);
		MPI_Wait(&recv_request2, &status);
	}
	else //middle ranks
	{    
		MPI_Irecv(top_ghost_row, u_size, MPI_INT, mpi_myrank-1, 
			1, MPI_COMM_WORLD, &recv_request);
		MPI_Irecv(bottom_ghost_row, u_size, MPI_INT, mpi_myrank+1,
			2, MPI_COMM_WORLD, &recv_request2);
		
		//send top row to previous rank
		MPI_Isend(my_rows[0], u_size, MPI_INT, mpi_myrank-1,
			2, MPI_COMM_WORLD, &send_request);
		MPI_Wait(&send_request, &status);
		//send bottom row to next rank
		MPI_Isend(my_rows[rows_per_rank-1], u_size, MPI_INT, mpi_myrank+1,
			1, MPI_COMM_WORLD, &send_request2);
		MPI_Wait(&send_request2, &status);
		//Receive Waits
		MPI_Wait(&recv_request, &status);
		MPI_Wait(&recv_request2, &status);
	}
}
