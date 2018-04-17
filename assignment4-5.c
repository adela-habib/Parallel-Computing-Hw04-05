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
//MPI Variables
MPI_Status status; 
MPI_Request send_request, recv_request;
MPI_Request recv_request2;
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

void *run_simulation(void *); //runs the pthreads simulation
int** sim_tick(int starting_row); //simulates game for single time "tick" for given rows
void pthread_id0(void); //function for pthread_id0

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
    sscanf(argv[1], "%d", &num_pthreads);
    //num_pthreads = argv[1]; //number of pthreads
    
    //Take in arg as double
    sscanf(argv[2],"%lf",&threshold);
    //threshold = argv[2];
    
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
    
//Randomly initialize universe
    for (int i=0; i<rows_per_rank; i++)
    {
        for (int j=0; j < u_size; j++)
        {
            //TODO - check this is the right way to initialize
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
    
//Create Pthreads
    pthread_t p_threads[num_pthreads];
    pthread_attr_t attr;
    pthread_attr_init (&attr);
  
    pthread_id = 0;
//1 or more pthreads
    if (num_pthreads > 0)
    {
      rows_per_thread = rows_per_rank/num_pthreads; //number of rows per thread

      //allocate memory for ghost rows
      top_ghost_row = (int *)calloc(u_size, sizeof(int));
      bottom_ghost_row = (int *)calloc(u_size, sizeof(int));
      
      //Create barrier to sync ticks in pthreads
      pthread_barrier_init(&pbarrier, NULL, num_pthreads);
      
      //pthreads 0 to num_pthreads run program function
      for(int i=0; i < num_pthreads; i++)
      {      
        //create pthreads
        pthread_create(&p_threads[i], 
        &attr, run_simulation, (void *) &pthread_id);
        
        pthread_id++;  
      }

      //join threads to make sure all pthreads in rank are done
      for(int i = 0; i < num_pthreads; i++){
              pthread_join(p_threads[i], NULL);
      }
    }
//No pthreads
    else
    {
      rows_per_thread = rows_per_rank;
      
      //allocate memory for ghost rows
      top_ghost_row = (int *)calloc(u_size, sizeof(int));
      bottom_ghost_row = (int *)calloc(u_size, sizeof(int));
      
      //TODO compiler does not like this
      run_simulation(0);
      
    }
      
    //Synchronize ranks
    MPI_Barrier(MPI_COMM_WORLD);  
    
//End time of program
    if (mpi_myrank == 0)
    {
        end_time = MPI_Wtime();
    }

//Free allocated memory
    free(top_ghost_row);
    free(bottom_ghost_row);
    for(int i=0; i < rows_per_rank; i++)
    {
        free(my_rows[i]);
    }
    free(my_rows);
    
// END -Perform a barrier and then leave MPI
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize();
    return 0;
}

/***************************************************************************/
/* Other Functions - You write as part of the assignment********************/
/***************************************************************************/

//runs a simulation for a tick for given rows
int** sim_tick(int starting_row)
{
  int curr_row = starting_row;
  //array to hold temporary table
  int ** temp_table = (int**)calloc(rows_per_thread, sizeof(int *));
  
  //For each row in the table
  for(int i = 0; i < rows_per_thread; i++)
  {    
    //array to hold temp results
    int * temp_cells = (int *)calloc(u_size, sizeof(int));
    
    //for each column in a row
    for(int curr_col = 0; curr_col < u_size; curr_col++)
    {
      int num_alive_neighbors = 0; //keep count of live neighbors  
      //set direction
      int left = (curr_col-1) % u_size;  //wrapping included by mod arithmatic
      int right = (curr_col + 1) % u_size; 
      int up = i+1;
      int down = i-1;
      
//count number of alive neighbors
      if (curr_row == 0) //top thread row - use top ghost row
      { 
        num_alive_neighbors = my_rows[curr_row][left] + my_rows[curr_row][right] + top_ghost_row[curr_col] + my_rows[down][curr_col] +
          top_ghost_row[left] + top_ghost_row[right] + my_rows[down][left] + my_rows[down][right];         
      }
      else if (curr_row == starting_row + rows_per_thread-1) //bottom thread row - use bottom ghost row
      {
        num_alive_neighbors = my_rows[curr_row][left] + my_rows[curr_row][right] + my_rows[up][curr_col] + bottom_ghost_row[curr_col] +
          my_rows[up][left] + my_rows[up][right]+ bottom_ghost_row[left] + bottom_ghost_row[right];
      }
      else //middle thread row
      {
        num_alive_neighbors = my_rows[curr_row][left] + my_rows[curr_row][right] + my_rows[up][curr_col] + my_rows[down][curr_col] +
          my_rows[up][left] + my_rows[up][right] + my_rows[down][left] + my_rows[down][right];
      }

//(ii.) Read all the statuses of the neighbors depending upon the number of live/dead neighbors for the current tick and set the appropiate live/dead count.
      //ADDITIONAL RANDOMNESS
      //if random value is greater than threshold - play the game
      if (GenVal(mpi_myrank*curr_row+curr_col) > threshold)
      {
          //1- Any live cell with fewer than two live neighbors dies, as if caused by under-population.
        if (num_alive_neighbors < 2)
        {
          temp_cells[curr_col] = DEAD; 
        }
        //2- Any live cell with two or three live neighbors lives on to the next generation
        //4- Any dead cell with exactly three live neighbors becomes a live cell, as if by reproduction
        else if(num_alive_neighbors < 4)
        {
          if((num_alive_neighbors == 3) & (my_rows[curr_row][curr_col]==DEAD))
          {
            temp_cells[curr_col] = ALIVE; // by reproduction
          }
          else
          {
            temp_cells[curr_col] = my_rows[curr_row][curr_col]; //stays the same
          }
              }
        //3- Any live cell with more than three live neighbors dies, as if by over-population
        else if(num_alive_neighbors > 3)
        {
          temp_cells[curr_col] = DEAD;
        }
      }
      //else, randomly select alive or dead
      else
      {
        //Genval gives a double from 0-1
        if (GenVal(mpi_myrank*curr_row+curr_col) > 0.5)
        { 
          temp_cells[curr_col] = ALIVE; //if random number if even, it is alive
        }
        else
        {
          temp_cells[curr_col] = DEAD; //random number is odd, it is dead
        }
      }
    } 
    
    //set the temporary table
    temp_table[i] = temp_cells;
    
    //increment the current row
    curr_row++;  
  } 
  return temp_table;
} 


//Funtion ran in each pthread
void *run_simulation(void *void_pthread_id_ptr)
{
    int *pthread_id_ptr = (int *)void_pthread_id_ptr;
    int pthread_id = *pthread_id_ptr;
   
    //play the game
    for(int t = 0; t < num_ticks; t++)
    {                              
        //(i.) Exchange row data with different MPI ranks for the ghost rows (the rows shared between MPI ranks). You exchange this (ghost) row data from thread 0 of each MPI rank. 
        //Be careful here on how to use mutexes so as to avoid deadlocks. Also note that this row exchange uses alive/dead status data from the previous tick.
        if (pthread_id == 0)
        {
          pthread_id0();
        }
        //If there are 0 or 1 pthreads
        if ((num_pthreads == 0) | (num_pthreads == 1))
        {
            int** temp_table = sim_tick(0);
          
            for (int i=0; i<rows_per_thread; i++)
            {
              for (int j=0; j<u_size; j++)
              {
                my_rows[i][j] = temp_table[i][j];
              }
            }
          
            //free temp table
            for(int i=0; i < rows_per_thread; i++)
            {
                  free(temp_table[i]);
            }
            free(temp_table);
        }
        //multiple pthreads
        else
        {
            //(ii.) Read all the statuses of the neighbors depending upon the number of live/dead neighbors for the current tick and set the appropiate live/dead count.        
            int starting_row = pthread_id * rows_per_thread;
            int** temp_table = sim_tick(starting_row);

            //(iii.) Finally, update the cell status
            for (int i=starting_row; i<rows_per_thread; i++)
            {
              for (int j=0; j<u_size; j++)
              {
                my_rows[i][j] = temp_table[i][j];
              }
            }

            //Free the temp table
            for(int i = 0; i < rows_per_thread; i++)
            {
              free(temp_table[i]);
            }
            free(temp_table);
        }
      
        //Use pthread barrier to make sure each thread finishes a tick before moving on
        if(num_pthreads > 0) 
        {
          pthread_barrier_wait(&pbarrier);

        }
    }   
}

void pthread_id0()
{
    if(mpi_myrank == 0) //first rank
    {
      //MPI_IRECV
      //receive top ghost row from rank 0
      MPI_Irecv(&top_ghost_row, u_size, MPI_INT, mpi_commsize-1, MPI_ANY_TAG, MPI_COMM_WORLD, &recv_request);
      //receive bottom ghost row from last rank
      MPI_Irecv(&bottom_ghost_row, u_size, MPI_INT, mpi_myrank+1, MPI_ANY_TAG, MPI_COMM_WORLD, &recv_request2);

      //MPI_ISEND
      //send top row to last rank
      MPI_Isend(&my_rows[0], u_size, MPI_INT, mpi_commsize-1, 1, MPI_COMM_WORLD, &send_request);
      MPI_Wait(&send_request, &status);
      //send last row to rank 1
      MPI_Isend(&my_rows[rows_per_rank-1], u_size, MPI_INT, mpi_myrank+1, 1, MPI_COMM_WORLD, &send_request);
      MPI_Wait(&send_request, &status);

      MPI_Wait(&recv_request, &status);
      MPI_Wait(&recv_request2, &status);

    }
    else if (mpi_myrank == mpi_commsize-1) //last rank
    {
      //MPI_IRECV
      //receive top ghost row from previous rank
      MPI_Irecv(&top_ghost_row, u_size, MPI_INT, mpi_myrank-1, MPI_ANY_TAG, MPI_COMM_WORLD, &recv_request);
      //receive bottom ghost row from first rank
      MPI_Irecv(&bottom_ghost_row, u_size, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &recv_request2);

      //MPI_ISEND
      //send top row to previous rank
      MPI_Isend(&my_rows[0], u_size, MPI_INT, mpi_myrank-1, 1, MPI_COMM_WORLD, &send_request);
      MPI_Wait(&send_request, &status);
      //send bottom row to first rank
      MPI_Isend(&my_rows[rows_per_rank-1], u_size, MPI_INT, 0, 1, MPI_COMM_WORLD, &send_request);
      MPI_Wait(&send_request, &status);

      MPI_Wait(&recv_request, &status);          
      MPI_Wait(&recv_request2, &status);

    }
    else //middle ranks
    {
      //MPI_IRECV
      //receive bottom ghost row from previous last rank
      MPI_Irecv(&bottom_ghost_row, u_size, MPI_INT, mpi_myrank-1, MPI_ANY_TAG, MPI_COMM_WORLD, &recv_request);
      //receives top ghost row from next rank
      MPI_Irecv(&top_ghost_row, u_size, MPI_INT, mpi_myrank+1, MPI_ANY_TAG, MPI_COMM_WORLD, &recv_request);

      //MPI_ISEND
      //send top row to previous rank
      MPI_Isend(&my_rows[0], u_size, MPI_INT, mpi_myrank-1, 1, MPI_COMM_WORLD, &send_request);
      MPI_Wait(&send_request, &status);
      //send bottom row to next rank
      MPI_Isend(&my_rows[rows_per_rank-1], u_size, MPI_INT, mpi_myrank+1, 1, MPI_COMM_WORLD, &send_request);
      MPI_Wait(&send_request, &status); 

      MPI_Wait(&recv_request, &status);          
      MPI_Wait(&recv_request2, &status);
    }
}
