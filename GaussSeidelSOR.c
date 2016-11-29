#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <mpi.h>
#include "utils.h"

int main(int argc, char ** argv) {
		int rank,size;
		int global[2],local[2]; //global matrix dimensions and local matrix dimensions (2D-domain, 2D-subdomain)
		int global_padded[2];   //padded global matrix dimensions (if padding is not needed, global_padded=global)
		int grid[2];            //processor grid dimensions
		int i,j,t;
		int global_converged=0,converged=0; //flags for convergence, global and per process
		MPI_Datatype dummy;     //dummy datatype used to align user-defined datatypes in memory
		double omega; 			//relaxation factor - useless for Jacobi

		struct timeval tts,ttf,tcs,tcf;   //Timers: total-tts,ttf, computation-tcs,tcf
		double ttotal=0,tcomp=0,total_time,comp_time;

		double ** U, ** u_current, ** u_previous, ** swap; //Global matrix, local current and previous matrices, pointer to swap between current and previous
		
		gettimeofday(&tts,NULL);

		MPI_Init(&argc,&argv);
		MPI_Comm_size(MPI_COMM_WORLD,&size);
		MPI_Comm_rank(MPI_COMM_WORLD,&rank);

		//----Read 2D-domain dimensions and process grid dimensions from stdin----//

		if (argc!=5) {
				fprintf(stderr,"Usage: mpirun .... ./exec X Y Px Py");
				exit(-1);
		}
		else {
				global[0]=atoi(argv[1]);
				global[1]=atoi(argv[2]);
				grid[0]=atoi(argv[3]);
				grid[1]=atoi(argv[4]);
		}

		//----Create 2D-cartesian communicator----//
		//----Usage of the cartesian communicator is optional----//

		MPI_Comm CART_COMM;         //CART_COMM: the new 2D-cartesian communicator
		int periods[2]={0,0};       //periods={0,0}: the 2D-grid is non-periodic
		int rank_grid[2];           //rank_grid: the position of each process on the new communicator

		MPI_Cart_create(MPI_COMM_WORLD,2,grid,periods,0,&CART_COMM);    //communicator creation
		MPI_Cart_coords(CART_COMM,rank,2,rank_grid);	                //rank mapping on the new communicator

		//----Compute local 2D-subdomain dimensions----//
		//----Test if the 2D-domain can be equally distributed to all processes----//
		//----If not, pad 2D-domain----//

		for (i=0;i<2;i++) {
				if (global[i]%grid[i]==0) {
						local[i]=global[i]/grid[i];
						global_padded[i]=global[i];
				}
				else {
						local[i]=(global[i]/grid[i])+1;
						global_padded[i]=local[i]*grid[i];
				}
		}

		//Initialization of omega
		omega=2.0/(1+sin(3.14/global[0]));

		//----Allocate global 2D-domain and initialize boundary values----//
		//----Rank 0 holds the global 2D-domain----//
		if (rank==0) {
				U=allocate2d(global_padded[0],global_padded[1]);   
				init2d(U,global[0],global[1]);
		}

		//----Allocate local 2D-subdomains u_current, u_previous----//
		//----Add a row/column on each size for ghost cells----//

		u_previous=allocate2d(local[0]+2,local[1]+2);
		u_current=allocate2d(local[0]+2,local[1]+2);   

		//----Distribute global 2D-domain from rank 0 to all processes----//

		//----Appropriate datatypes are defined here----//
		/*****The usage of datatypes is optional*****/

		//----Datatype definition for the 2D-subdomain on the global matrix----//

		MPI_Datatype global_block;
		MPI_Type_vector(local[0],local[1],global_padded[1],MPI_DOUBLE,&dummy);
		MPI_Type_create_resized(dummy,0,sizeof(double),&global_block);
		MPI_Type_commit(&global_block);

		//----Datatype definition for the 2D-subdomain on the local matrix----//

		MPI_Datatype local_block;
		MPI_Type_vector(local[0],local[1],local[1]+2,MPI_DOUBLE,&dummy);
		MPI_Type_create_resized(dummy,0,sizeof(double),&local_block);
		MPI_Type_commit(&local_block);

		//----Rank 0 defines positions and counts of local blocks (2D-subdomains) on global matrix----//
		int * scatteroffset, * scattercounts;
		if (rank==0) {
				scatteroffset=(int*)malloc(size*sizeof(int));
				scattercounts=(int*)malloc(size*sizeof(int));
				for (i=0;i<grid[0];i++)
						for (j=0;j<grid[1];j++) {
								scattercounts[i*grid[1]+j]=1;
								scatteroffset[i*grid[1]+j]=(local[0]*local[1]*grid[1]*i+local[1]*j);
						}
		}

		//----Rank 0 scatters the global matrix----//
		if(rank==0){

				//Copy local[0] rows instead of sending
				for(int i=0; i<local[0]; ++i)
						for(int j=0; j<local[1]; ++j)
								u_current[i+1][j+1]=U[i][j];

				//Send the corresponding block to each process
				for(int i=1; i<size; ++i) MPI_Send(&U[0][0]+scatteroffset[i],	1, global_block, i,	i, MPI_COMM_WORLD);
		}


		//Each process receives the local data
		if(rank!=0) MPI_Recv(&u_current[1][1],	1, local_block,	0, rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
																						//&req);

		/**** DEBUG SECTION ****/
		//if(rank==0) print2d(U, global_padded[0], global_padded[1]);
		/***********************/

		/*Make sure u_current and u_previous are
		  both initialized*/

		if (rank==0)
				free2d(U,global_padded[0],global_padded[1]);

		//----Find the 4 neighbors with which a process exchanges messages----//

		int north, south, east, west;
		
		north = south = east = west = -1;
		if(rank_grid[0]!=0 && rank_grid[1]!=0 && rank_grid[0]!=grid[0]-1 && rank_grid[1]!=grid[1]-1){
			north=rank-grid[1];
			south=rank+grid[1];
			east=rank+1;
			west=rank-1;
		}
		if(rank_grid[0]!=0) north = rank-grid[1];
		if(rank_grid[1]!=0) west = rank-1;
		if(rank_grid[0]!=grid[0]-1) south = rank+grid[1];
		if(rank_grid[1]!=grid[1]-1) east = rank+1;
		
		
		/**** DEBUG SECTION ****/
		//for(int i=0; i<size; ++i) { if(rank==i) printf("Process[%d]: N: %d S: %d W: %d E: %d \n", rank, north, south, west, east); }		
		/***********************/

		/*Fill your code here*/


		/*Make sure you handle non-existing
		  neighbors appropriately*/

		//---Define the iteration ranges per process-----//

		int i_min,i_max,j_min,j_max;
		
		//If the process is not in the right or the bottom border then:
		//Starting X point:		[1][1]
		//Ending X point:		[local[0]][1]
		//Starting Y point:		[1][1]
		//Ending Y point:		[local[0]][local[1]]
		
		i_min=1;
		i_max=local[0]+1;
		j_min=1;
		j_max=local[1]+1;
		
		//If the process is in the top border:
		//Starting X point: 	[2][1]
		//If the process is in the bottom border:
		//Ending X point:		[local[0]-x_padding_offset][1]
		//If the process is int the left border:
		//Starting Y point: 	[1][2]
		//If the process is in the right border:
		//Ending Y point:		[local[0]][local[1]-y_padding_offset][1]
		if(rank_grid[0]==0) i_min++;
		if(rank_grid[0]==grid[0]-1) i_max-=(global_padded[0]-global[0]+1);
		if(rank_grid[1]==0) j_min++;
		if(rank_grid[1]==grid[1]-1) j_max-=(global_padded[1]-global[1]+1);

		/**** DEBUG SECTION ****/
		//for(int i=0; i<size; ++i) { if(rank==i) printf("Process[%d]: Xmin: %d Xmax: %d Ymin: %d Ymax: %d \n", rank, i_min, i_max, j_min, j_max); }
		/***********************/
		
		/* Define appropriate datatype for table column */
		MPI_Datatype COL;
		MPI_Type_vector(i_max-i_min+1,1,local[1]+2,MPI_DOUBLE,&COL);
		MPI_Type_commit(&COL);
		

		/*Fill your code here*/


						if(south!=-1) MPI_Send(&u_current[i_max-1][1], j_max-j_min, MPI_DOUBLE, south, 17, MPI_COMM_WORLD);
						if(east!=-1) MPI_Send(&u_current[1][j_max-1], 1, COL, east, 17, MPI_COMM_WORLD);



		/*Three types of ranges:
		  -internal processes
		  -boundary processes
		  -boundary processes and padded global array
		 */





		//************************************//


		gettimeofday(&tcs, NULL);

		for(int i=0; i<local[0]+2; ++i)
			for(int j=0; j<local[1]+2; ++j)
				u_previous[i][j]=u_current[i][j];

		gettimeofday(&tcf, NULL);
		tcomp+=(tcf.tv_sec-tcs.tv_sec)+(tcf.tv_usec-tcs.tv_usec)*0.000001;

		//----Computational core----//   
#ifdef TEST_CONV
		for (t=0;t<T && !global_converged;t++) {
#endif
#ifndef TEST_CONV
#undef T
#define T 256
				for (t=0;t<T;t++) {
#endif
						/*Compute and Communicate*/

						/*Add appropriate timers for computation*/

						swap=u_previous;
						u_previous=u_current;
						u_current=swap;

						if(north!=-1) MPI_Recv(&u_current[0][1], j_max-j_min, MPI_DOUBLE, north, 17, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						if(west!=-1) MPI_Recv(&u_current[1][0], 1, COL, west, 17, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

						gettimeofday(&tcs, NULL);
						//printf("Process[%d]: Alles gut\n", rank); //DEBUG PRINTING
						for(int i=i_min; i<i_max; ++i){
							for(int j=j_min; j<j_max; ++j){
								u_current[i][j]=u_previous[i][j]+(u_current[i-1][j]+u_previous[i+1][j]+u_current[i][j-1]+u_previous[i][j+1]-4*u_previous[i][j])*omega/4.0;
							}
						}

						gettimeofday(&tcf, NULL);
						tcomp+=(tcf.tv_sec-tcs.tv_sec)+(tcf.tv_usec-tcs.tv_usec)*0.000001;

						if(north!=-1) MPI_Send(&u_current[1][1], j_max-j_min, MPI_DOUBLE, north, 17, MPI_COMM_WORLD);
						if(west!=-1) MPI_Send(&u_current[1][1], 1, COL, west, 17, MPI_COMM_WORLD);

						if(south!=-1) MPI_Recv(&u_current[i_max][1], j_max-j_min, MPI_DOUBLE, south, 17, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						if(east!=-1) MPI_Recv(&u_current[1][j_max], 1, COL, east, 17, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

						if(south!=-1) MPI_Send(&u_current[i_max-1][1], j_max-j_min, MPI_DOUBLE, south, 17, MPI_COMM_WORLD);
						if(east!=-1) MPI_Send(&u_current[1][j_max-1], 1, COL, east, 17, MPI_COMM_WORLD);

#ifdef TEST_CONV
						if (t%C==0) {
								/*Test convergence*/
								converged=converge(u_previous,u_current,local[0]+2,local[1]+2); //Check local convergence
								//MPI_Barrier(MPI_COMM_WORLD); //Wait for all processes to finish
								//printf("Process[%d]: %d\n", rank, converged); //DEBUG PRINTING
								MPI_Allreduce(&converged, &global_converged, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD); //global_converged=1 only if all have converged
						}		
#endif
				}
				
				gettimeofday(&ttf,NULL);

				ttotal=(ttf.tv_sec-tts.tv_sec)+(ttf.tv_usec-tts.tv_usec)*0.000001;

				MPI_Reduce(&ttotal,&total_time,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
				MPI_Reduce(&tcomp,&comp_time,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);



				//----Rank 0 gathers local matrices back to the global matrix----//

				if (rank==0) {
						U=allocate2d(global_padded[0],global_padded[1]);
				}

				//Each process sends the local data
				if(rank!=0) MPI_Send(&u_current[1][1],	1, local_block,	0, rank, MPI_COMM_WORLD);

				//----Rank 0 gathers the global matrix----//
				if(rank==0){

					//Copy local[0] rows instead of receiving
					for(int i=0; i<local[0]; ++i)
						for(int j=0; j<local[1]; ++j)
							U[i][j]=u_current[i+1][j+1];

					//Receive the corresponding block from each process
					for(int i=1; i<size; ++i) MPI_Recv(&U[0][0]+scatteroffset[i],	1, global_block, i,	i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

				}

				//----Printing results----//

				//**************TODO: Change "Jacobi" to "GaussSeidelSOR" or "RedBlackSOR" for appropriate printing****************//
				
				if (rank==0) {
						printf("Jacobi X %d Y %d Px %d Py %d Iter %d ComputationTime %lf TotalTime %lf midpoint %lf\n",global[0],global[1],grid[0],grid[1],t,comp_time,total_time,U[global[0]/2][global[1]/2]);	


#ifdef PRINT_RESULTS
						char * s=malloc(50*sizeof(char));
						sprintf(s,"resJacobiMPI_%dx%d_%dx%d",global[0],global[1],grid[0],grid[1]);
						fprint2d(s,U,global[0],global[1]);
						free(s);
#endif
				}
				MPI_Finalize();
				return 0;
		}
