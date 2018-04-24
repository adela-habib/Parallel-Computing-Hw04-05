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

#define u_size 16384 //world size
#define num_ticks 128 //128 ticks
#define s 16 //block size for heat map

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
int rows_per_thread; //number of rows per pthread
int **my_rows; //allocated number of rows per rank
int **ghost_rows; //parallel universe for each rank
int *top_ghost_row; //ghost row (top)
int *bottom_ghost_row; //ghost row (bottom)
pthread_barrier_t pbarrier; //pthread barrier

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
  //MPI Initialization
  MPI_Init( &argc, &argv);
  MPI_Comm_size( MPI_COMM_WORLD, &mpi_commsize);
  MPI_Comm_rank( MPI_COMM_WORLD, &mpi_myrank);
  
  // Init 16,384 RNG streams - each rank has an independent stream
  InitDefault();
  MPI_Barrier( MPI_COMM_WORLD );
    
  //Input vlues 
  sscanf(argv[1], "%d", &num_pthreads); //number of pthreads
  sscanf(argv[2], "%lf", &threshold); //threshold value
    
  //Start time of program
  if (mpi_myrank == 0)
    {
      start_time = MPI_Wtime();
    }   

  //determine number of rows for each rank and each thread
  rows_per_rank = u_size/mpi_commsize; 

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
	  if (GenVal(mpi_myrank*rows_per_rank+i) > 0.5)
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
  pthread_t p_threads[num_pthreads-1];
  pthread_attr_t attr;
  pthread_attr_init (&attr);

  //If there are 2 or more pthreads
  if (num_pthreads > 1)
    {
      //calculate the number of rows per thread
      rows_per_thread = rows_per_rank/num_pthreads; 
     
      //Create barrier to sync pthreads each tick
      pthread_barrier_init(&pbarrier, NULL, num_pthreads);
		
      //pthreads 1 to num_pthreads create and run simulation function
      int i;
      for(i=1; i < num_pthreads; i++)
	{ 
	  int *pthread_id_ptr = malloc(sizeof(int)); //initialize pthread id
	  *pthread_id_ptr = i; //set pthread id
	  
	  //create pthread, send to run_simulation function, pass pthread id
	  pthread_create(&p_threads[i-1], &attr, run_simulation, (void *)pthread_id_ptr);
	}
    }
  //Else there are less than 2 pthreads
  else
    {
      //the number of rows/thread = rows/rank
      rows_per_thread = rows_per_rank;
    }

  int *pthread_id_ptr=malloc(sizeof(int)); //initialize pthread id
  *pthread_id_ptr = 0; //set pthread id for thread 0

  //send thread 0 to run_simulation function
  run_simulation(pthread_id_ptr); 

  //join threads to make sure all pthreads in rank are done
  for(int i = 0; i < num_pthreads - 1; i++)
    {
      pthread_join(p_threads[i], NULL);
    }

  MPI_Barrier(MPI_COMM_WORLD);

  //End time of program
  if (mpi_myrank == 0)
    {
      end_time = MPI_Wtime();
    }

  /****************************************************************
   ***** OUTPUT CODE **********************************************
   ***** (1) 1,204x1,204 grid *************************************
   ***** (2) entire universe  *************************************
   *****************************************************************/

  //CODE FOR (1) - If statement used in order to bypass if needed
  int run_heatmap = 1;
  if(run_heatmap == 1){
    //Output Variables
    int rows_per_output = rows_per_rank/s; //number of rows to output (per rank)
    int output_size = u_size/s; //output grid size
    int **output_rows; //the calculated values to output
    MPI_File fh; //MPI file
    MPI_Status status; //MPI status
    char *output_char; //character array to output
    int char_size = rows_per_output*(4*output_size+1); //size of character array

    //Allocate space for char output array
    output_char = (char *)calloc(char_size, sizeof(char *));

    //Allocate space for rows & columns of outputted values
    output_rows = (int **)calloc(rows_per_output, sizeof(int*));
    for(int i=0; i<rows_per_output; i++)
      {
        output_rows[i] = (int *)calloc(output_size, sizeof(int));
      }

    //Calculate (16 x 16) values
    for(int i=0; i<rows_per_output; i++)
      {
        for(int j=0; j<output_size; j++)
	  {
	    //initialize value
	    int value = 0;
	    //sum over 16x16 squares
	    for(int k=0; k<s; k++)
	      {
		for(int l=0; l<s; l++)
		  {
		    value += my_rows[i*s+k][j*s+l];
		  }
	      }
	    //add it to output grid
	    output_rows[i][j] = value;
	  }
      }

    //Convert to Single Char Array
    int ind = 0;
    int x = 0;
    int y = 0;
    int z = 0;

    for(int i=0; i<rows_per_output; i++)
      {
        for(int j=0; j<output_size; j++)
	  {
	    //If the value is less than zero
	    if (output_rows[i][j]<10)
	      {
		output_char[ind] = '0'; //add a zero
		ind++;
		output_char[ind] = '0'; //add a zero
		ind++;
		output_char[ind] = '0' + output_rows[i][j]; //add the value
		ind++;
	      }
	    //If the value is less than 100
	    else if (output_rows[i][j]<100)
	      {
		z = output_rows[i][j]%10; //calculate the ones place
		y = (output_rows[i][j]-z)/10; //calculate the tens place
		output_char[ind] = '0'; //add a zero
		ind++;
		output_char[ind] = '0' + y; //add tens place
		ind++;
		output_char[ind] = '0' + z; //add ones place
		ind++;
	      }
	    //Else the value is greater than 100
	    else
	      {
		z = output_rows[i][j]%10; //calculate ones place
		y = ((output_rows[i][j]-z)/10)%10; //calculate tens place
		x = (output_rows[i][j]-z-(y*10))/100; //calculate hundreds place
		output_char[ind] = '0' + x; //add hundres place
		ind++;
		output_char[ind] = '0' + y; //add tens place
		ind++;
		output_char[ind] = '0' + z; //add ones place
		ind++;
	      }
	    output_char[ind] = ','; //add comma between each value
	    ind++;
	  }
        output_char[ind] = '\n'; //add newline between each outputted row
        ind++;
      }

    //MPI File operations
    int offset = mpi_myrank*char_size; //offset for each rank to write to file
    char *filename = "output.csv"; //outputted file name

    //MPI File operations
    MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    MPI_File_write_at(fh, offset, output_char, char_size, MPI_CHAR, &status);
    MPI_File_close(&fh);

    //Free allocated space
    for(int i=0; i<rows_per_output; i++)
    {
      free(output_rows[i]);
    }
    free(output_rows);
    free(output_char);  
  }

  //CODE FOR (2) - If statement used in order to bypass if needed
  int run_io = 1;
  if(run_io == 1)
    {
      //Total Universe Output Variables
      MPI_File file2; //MPI file 
      MPI_Status status2; //MPI status
      char *universe_char; //character array
      int uout_size = rows_per_rank*(2*u_size+1); //size of character array
      
      //Allocate space for the outputted universe array
      universe_char = (char *)calloc(uout_size, sizeof(char *));
      
      //Convert to single char array
      int ind2 = 0;
      for (int i=0; i<rows_per_rank; i++)
	{
	  for (int j=0; j<u_size; j++)
	    {
	      //if the value is 0
	      if (my_rows[i][j] == 0)
		{
		  universe_char[ind2] = '0'; //add a zero
		  ind2++;
		}
	      //else the value is 1
	      else
		{
		  universe_char[ind2] = '1'; //add a one
		  ind2++;
		}
	      universe_char[ind2] = ','; //add a comma between each value
	      ind2++;
	    }
	  universe_char[ind2] = '\n'; //add a newline between each outputted row
	  ind2++;
	}
      
      //MPI Second File Operations
      int offset2 = mpi_myrank*uout_size; //offset for each rank to write at
      char *filename2 = "Universe.csv"; //outputted file name

      //MPI File operations
      MPI_File_open(MPI_COMM_WORLD, filename2, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &file2);
      MPI_File_write_at(file2, offset2, universe_char, uout_size, MPI_CHAR, &status2);
      MPI_File_close(&file2);
      
      //Free Allocated Memory
      free(universe_char);
    }
  
  MPI_Barrier(MPI_COMM_WORLD);
  //Print compute time
  if(mpi_myrank == 0)
    {
      printf("\nFinished. Run Time: %f\n\n", end_time - start_time);
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
  //Set the pthread id (passed as a void pointer)
  int pthread_id; 
  pthread_id = *((int *)void_pthread_id_ptr); 
  free(void_pthread_id_ptr);

  //play the game
  for(int t = 0; t < num_ticks; t++)
    {
      //Print each tick
      if(mpi_myrank==0)
	{
	  printf("Tick %d out of %d\n", t+1, num_ticks);
	}
      //Exchange row data with different MPI ranks for the ghost rows
      if (pthread_id == 0)
	{
	  pthread_id0();
	}
      //If there are 0 or 1 pthreads
      if ( (num_pthreads == 0) | (num_pthreads == 1) )
	{
	  //Simulate a tick at starting row 0
	  sim_tick(0);

	  //update cell states
	  for (int i=0; i<rows_per_thread; i++)
	    {
	      for (int j=0; j<u_size; j++)
		{
		  my_rows[i][j] = ghost_rows[i][j];
		}
	    }
	}
      //Else there are multiple pthreads
      else 
	{
	  int starting_row = pthread_id * rows_per_thread; //set the starting row for each thread
      	  sim_tick(starting_row); //simulate one tick
	  
	  //Update the cell status
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
      //for each column in a row
      for(int curr_col = 0; curr_col < u_size; curr_col++)
        {
	  int num_alive_neighbors = 0; //keep count of live neighbors  

	  //set direction
	  int left = (curr_col - 1 + u_size) % u_size;
	  int right = (curr_col + 1) % u_size;
	  int up = i-1;
	  int down = i+1;
			
	  //Read all the statuses of the neighbors and count number of cells alive
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
	  else //middle thread rows
	    {
	      num_alive_neighbors = my_rows[i][left] + my_rows[i][right] + my_rows[up][curr_col] + my_rows[down][curr_col] +
		my_rows[up][left] + my_rows[up][right] + my_rows[down][left] + my_rows[down][right];
	    }

	  //if random value is greater than threshold - play the game
	  if (GenVal(mpi_myrank*rows_per_rank+i) > threshold)
	    {
	      //RULE 1- Any live cell with fewer than two live neighbors dies, as if caused by under-population.
	      if (num_alive_neighbors < 2)
		{
		  ghost_rows[i][curr_col] = DEAD;
		}
	      else if(num_alive_neighbors < 4)
		{
		  //RULE 4- Any dead cell with exactly three live neighbors becomes a live cell, as if by reproduction
		  if((num_alive_neighbors == 3) & (my_rows[i][curr_col]==DEAD))
		    {
		      ghost_rows[i][curr_col] = ALIVE; 
		    }
		  //RULE 2- Any live cell with two or three live neighbors lives on to the next generation
		  else
		    {
		      ghost_rows[i][curr_col] = my_rows[i][curr_col];
		    }
		}
	      //RULE 3- Any live cell with more than three live neighbors dies, as if by over-population
	      else if(num_alive_neighbors > 3)
		{
		  ghost_rows[i][curr_col] = DEAD;
		}
	    }
	  //else, randomly select alive or dead
	  else
	    {
	      if (GenVal(mpi_myrank*rows_per_rank+i) > 0.5 )
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

//(i.) Exchange row data with different MPI ranks for the ghost rows (the rows shared between MPI ranks).
void pthread_id0()
{
  //MPI Send/Recv Variables
  MPI_Status status; 
  MPI_Request send_request, send_request2, recv_request, recv_request2;
  
  if (mpi_myrank == 0) //first rank
    {
      //POST RECEIVES
      //receive from last rank into top ghost row
      MPI_Irecv(top_ghost_row, u_size, MPI_INT, mpi_commsize-1, 
		1, MPI_COMM_WORLD, &recv_request);
      //receive from rank 2 into bottom ghost row
      MPI_Irecv(bottom_ghost_row, u_size, MPI_INT, mpi_myrank+1,
		2, MPI_COMM_WORLD, &recv_request2);

      //SEND
      //send top row to last rank
      MPI_Isend(my_rows[0], u_size, MPI_INT, mpi_commsize-1,
		2, MPI_COMM_WORLD, &send_request);
      MPI_Wait(&send_request, &status);
      //send last row to rank 1
      MPI_Isend(my_rows[rows_per_rank-1], u_size, MPI_INT, mpi_myrank+1,
		1, MPI_COMM_WORLD, &send_request2);
      MPI_Wait(&send_request2, &status);
      
      //RECEIVE ghost rows
      MPI_Wait(&recv_request, &status);
      MPI_Wait(&recv_request2, &status);
    }
  else if (mpi_myrank == mpi_commsize-1) //last rank
    {
      //POST RECEIVES
      //receive from previous rank into top ghost row
      MPI_Irecv(top_ghost_row, u_size, MPI_INT, mpi_myrank-1,
		1, MPI_COMM_WORLD, &recv_request);
      //receive from rank 0 into bottom ghost row
      MPI_Irecv(bottom_ghost_row, u_size, MPI_INT, 0,
		2, MPI_COMM_WORLD, &recv_request2);

      //SEND
      //send top row to previous rank
      MPI_Isend(my_rows[0], u_size, MPI_INT, mpi_myrank-1,
		2, MPI_COMM_WORLD, &send_request);
      MPI_Wait(&send_request, &status);
      //send last row to first rank
      MPI_Isend(my_rows[rows_per_rank-1], u_size, MPI_INT, 0,
		1, MPI_COMM_WORLD, &send_request2);
      MPI_Wait(&send_request2, &status);

      //RECEIVE ghost rows
      MPI_Wait(&recv_request, &status);
      MPI_Wait(&recv_request2, &status);
    }
  else //middle ranks
    {
      //POST RECEIVES
      //receive from previous rank into top ghost row
      MPI_Irecv(top_ghost_row, u_size, MPI_INT, mpi_myrank-1, 
		1, MPI_COMM_WORLD, &recv_request);
      //receive from next rank into bottom ghost row
      MPI_Irecv(bottom_ghost_row, u_size, MPI_INT, mpi_myrank+1,
		2, MPI_COMM_WORLD, &recv_request2);

      //SENDS
      //send top row to previous rank
      MPI_Isend(my_rows[0], u_size, MPI_INT, mpi_myrank-1,
		2, MPI_COMM_WORLD, &send_request);
      MPI_Wait(&send_request, &status);
      //send bottom row to next rank
      MPI_Isend(my_rows[rows_per_rank-1], u_size, MPI_INT, mpi_myrank+1,
		1, MPI_COMM_WORLD, &send_request2);
      MPI_Wait(&send_request2, &status);

      //RECEIVE ghost rows
      MPI_Wait(&recv_request, &status);
      MPI_Wait(&recv_request2, &status);
    }
}
