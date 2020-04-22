#ifndef WAVE_SOLVE
#define WAVE_SOLVE

#include <string>   
#include <fstream>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <memory>
#include "problem_setup.hpp"
#include "Darray1.h"
#include "Darray2.h"
#include "Iarray2.h"
#include "twilight_2d.h"
#include "domain_decomposition.h"
using namespace std;


class Wave_Solve {
public:
	ProblemSetup setup;

	// MPI derived data types for communication
    MPI_Datatype sub_send_left;
    MPI_Datatype sub_recv_left;
    MPI_Datatype sub_send_right;
    MPI_Datatype sub_recv_right;
    MPI_Datatype sub_send_up;
    MPI_Datatype sub_recv_up;
    MPI_Datatype sub_send_down;
    MPI_Datatype sub_recv_down;

    // Requests for persistent communication
    MPI_Request send_req[4];
    MPI_Request recv_req[4];

	// Communication data
	int up_neigh, down_neigh, left_neigh, right_neigh;
	int size_lr, size_ud;
	double* IO_buf;

	Iarray2 mask, ghost_list_left, ghost_list_right, ghost_list_up, ghost_list_down, bdry_list;
	Darray1 x, y, dist_list_left, dist_list_right, dist_list_up, dist_list_down;
	double t;
	int n_ghost_right, n_ghost_left, n_ghost_up, n_ghost_down, n_bdry;

	int istart, iend, jstart, jend, N, M;
	int j_left_start, j_left_end, i_down_start, i_down_end;
	int j_right_start, j_right_end, i_up_start, i_up_end;



	Wave_Solve(Subdomain Local_Grid, int node_ID, MPI_Comm CART_COMM, MPI_Comm IO_Comm);
	virtual ~Wave_Solve() {if( IO_buf != 0 ) delete[] IO_buf;};
	void    Enforce_BC(Darray2 &v);
	void    Compute_Laplacian(Darray2& w, Darray2& lap);
	void    Compute_Laplacian_NB(Darray2& w, Darray2& lap, MPI_Request* recv_req);
	void    Set_Initial_Data(Darray2& wm, Darray2& w, Darray2& lap,Darray2& u);
	double  Compute_Laplacian_Error(const int flag, Darray2& lap, MPI_Comm CART_COMM);
	double  Compute_Solution_Error(const int flag, Darray2& w, MPI_Comm CART_COMM);
	void    Time_Step(Darray2& wm, Darray2& w, Darray2& wp, Darray2& lap);
	void    Taylor_Expand(Darray2& wm, Darray2& w, Darray2& lap);
	double  forcing(const double x, const double y, const double t);
	void    Solve_PDE(Darray2& wm, Darray2& w, Darray2& wp, Darray2& lap,Darray2& u, double* w_ptr, MPI_Comm CART_COMM);


// Define the mask and sweep arrays which will be used to enforce boundary conditions.
void Compute_Mask(){
     double level_set_val;
     double distance;
     double x_val;
     int ghost_ctr = 0;
     int num_int = 0;
     int ctr = 1;
     int n_left = 0;
     int n_right = 0;
     int n_up = 0;
     int n_down = 0;
     int dir_ctr[4];
     int n_phys = 0;

     bool bdry_check[4];
     for(int i = istart;i<=iend;i++){
         x_val = x(i);
         for(int j = jstart;j<=jend;j++){
             level_set_val = setup.Level_Set(x_val,y(j));

             // Interior point
             if(level_set_val < 0.0){
                 mask(i,j) = 1;
                 num_int++;
             }
             // Exterior point
             else if(level_set_val > 0.0){
                 mask(i,j) = 0;
             }
             // Boundary point
             else{
                 mask(i,j) = -2;
                 n_phys++;
             }
         }
     }
    // AAL need to identify which direction to use secant method along
    // Possibly use the componets of the normal direction??

     // Identify an interior ghost point candidate
    for(int i = istart+1;i<=iend-1;i++){
        for(int j = jstart+1;j<=jend-1;j++){
            if(mask(i,j) == 1){
                bdry_check[0] = ((mask(i-1,j) == 0) || (mask(i-1,j) == -3));
                bdry_check[1] = ((mask(i+1,j) == 0) || (mask(i+1,j) == -3));
                bdry_check[2] = ((mask(i,j-1) == 0) || (mask(i,j-1) == -3));
                bdry_check[3] = ((mask(i,j+1) == 0) || (mask(i,j+1) == -3));


	             // Update mask for exterior ghost points and related interior point
	            if(bdry_check[0]){
	             	mask(i,j) = -1;
	             	mask(i-1,j) = -3;
	             	ghost_ctr++;
	             	n_left++;
	             	num_int--;
	            }
	            if(bdry_check[1]){
	             	mask(i,j) = -1;
	             	mask(i+1,j) = -3;
	             	ghost_ctr++;
	             	n_right++;
	             	num_int--;
	            }
	            if(bdry_check[2]){
	             	mask(i,j) = -1;
	             	mask(i,j-1) = -3;
	             	ghost_ctr++;
	             	n_down++;
	             	num_int--;
	            }
	            if(bdry_check[3]){
	             	mask(i,j) = -1;
	             	mask(i,j+1) = -3;
	             	ghost_ctr++;
	             	n_up++;
	             	num_int--;
	            }
            }
        }
    }


    // Store list of indices, directions, and distances for each ghost point
    if(n_left > 0){
    	ghost_list_left.define(1,1,n_left,1,2);
    	dist_list_left.define(1,1,n_left);
    }
    if(n_right > 0){
    	ghost_list_right.define(1,1,n_right,1,2);
    	dist_list_right.define(1,1,n_right);
    }
    if(n_up > 0){
    	ghost_list_up.define(1,1,n_up,1,2);
    	dist_list_up.define(1,1,n_up);
    }
    if(n_down > 0){
    	ghost_list_down.define(1,1,n_down,1,2);
    	dist_list_down.define(1,1,n_down);
    }
    if(n_phys > 0){
    	bdry_list.define(1,1,n_phys,1,2);
    }

    n_ghost_right = n_right;
    n_ghost_left = n_left;
    n_ghost_up = n_up;
    n_ghost_down = n_down;
    n_bdry = n_phys;

    for(int i = 0;i<4;i++) dir_ctr[i] = 1;

    for(int i = istart+1;i<=iend-1;i++){
        for(int j = jstart+1;j<=jend-1;j++){
        	if(mask(i,j) == -1){

        		// Right
        		if(mask(i+1,j) == -3){
        			ghost_list_right(dir_ctr[0],1) = i;
        			ghost_list_right(dir_ctr[0],2) = j;
        			dist_list_right(dir_ctr[0]) = setup.Dist_to_bdry(x(i), x(i+1), y(j), y(j),1)/setup.hx;
        			dir_ctr[0] = dir_ctr[0]+1;
        		}

        		// Up
        		if(mask(i,j+1) == -3){
        			ghost_list_up(dir_ctr[1],1) = i;
        			ghost_list_up(dir_ctr[1],2) = j;
        			dist_list_up(dir_ctr[1]) = setup.Dist_to_bdry(x(i), x(i), y(j), y(j+1),2)/setup.hy;
        			dir_ctr[1] = dir_ctr[1]+1;
        		}

        		// Left
        		if(mask(i-1,j) == -3){
        			ghost_list_left(dir_ctr[2],1) = i;
        			ghost_list_left(dir_ctr[2],2) = j;
        			dist_list_left(dir_ctr[2]) = setup.Dist_to_bdry(x(i), x(i-1), y(j), y(j),1)/setup.hx;
        			dir_ctr[2] = dir_ctr[2]+1;
        		}

        		// Down
        		if(mask(i,j-1) == -3){
        			ghost_list_down(dir_ctr[3],1) = i;
        			ghost_list_down(dir_ctr[3],2) = j;
        			dist_list_down(dir_ctr[3]) = setup.Dist_to_bdry(x(i), x(i), y(j), y(j-1),2)/setup.hy;
        			dir_ctr[3] = dir_ctr[3]+1;
        		}
        	}
        } 
    }


    for(int i = istart;i<=iend;i++){
        for(int j = jstart;j<=jend;j++){
        	// If a grid point is on the physical boundary, track its indices
        	if(mask(i,j) == -2){
        		bdry_list(ctr,1) = i;
        		bdry_list(ctr,2) = j;
        		ctr++;
        	} 
        }
    }
};

// Compute discrete energy for this problem. Note: Only Dirich/Neumann problems
// with homogeneous boundary conditions will gives energy conservation.
double Compute_Energy(Darray2& wm, Darray2& w, Darray2& lap, MPI_Comm CART_COMM){
	double dt2 = setup.dt2;
	double energy = 0.0;
	double total_energy;
	double idt2 = setup.idt2;

	for(int i=istart+1;i<=iend-1;i++){
		for(int j=jstart+1;j<=jend-1;j++){
			if(mask(i,j) == 1){
				energy = energy + idt2*pow(wm(i,j) - w(i,j),2) - w(i,j)*lap(i,j);
			}
		}
	}
	MPI_Allreduce(&energy, &total_energy, 1, MPI_DOUBLE, MPI_SUM, CART_COMM);

	return 0.5*setup.hx*setup.hy*total_energy;
};

// Set up subarray type for MPI communication
void Setup_Subarrays(const int nolp){
	int ndims = 2;
	int M = size_ud - 2*nolp;
	int N = size_lr - 2*nolp;

	/************************************************************/
	/*	    		COMMUNICATORS IN X DIRECTION				*/
	/************************************************************/

	// Left send
    int array_of_subsizes[2] = {nolp,M};
    int array_of_sizes[2] = {size_lr,size_ud};

    int array_of_starts[2] = {nolp,nolp};
    int order = MPI_ORDER_C;

    MPI_Type_create_subarray(ndims,
                             array_of_sizes,
                             array_of_subsizes,
                             array_of_starts,
                             order, MPI_DOUBLE, &sub_send_left);
    MPI_Type_commit(&sub_send_left);
    
    // Right receive
    array_of_starts[0] = size_lr-1;
    array_of_starts[1] = nolp;
    MPI_Type_create_subarray(ndims,
                             array_of_sizes,
                             array_of_subsizes,
                             array_of_starts,
                             order, MPI_DOUBLE, &sub_recv_right);
    MPI_Type_commit(&sub_recv_right);

    // Right send
    array_of_starts[0] = size_lr - nolp - 1;
    array_of_starts[1] = nolp;    
    MPI_Type_create_subarray(ndims,
                             array_of_sizes,
                             array_of_subsizes,
                             array_of_starts,
                             order, MPI_DOUBLE, &sub_send_right);
    MPI_Type_commit(&sub_send_right);
    
    // Left receive
    array_of_starts[0] = 0;
    array_of_starts[1] = nolp;
    MPI_Type_create_subarray(ndims,
                             array_of_sizes,
                             array_of_subsizes,
                             array_of_starts,
                             order, MPI_DOUBLE, &sub_recv_left);
    MPI_Type_commit(&sub_recv_left);


	/************************************************************/
	/*	    		COMMUNICATORS IN Y DIRECTION				*/
	/************************************************************/
	
	array_of_subsizes[0] = N;
	array_of_subsizes[1] = nolp;

	// Up send
	array_of_starts[0] = nolp;
	array_of_starts[1] = size_ud-1-nolp;
    MPI_Type_create_subarray(ndims,
                             array_of_sizes,
                             array_of_subsizes,
                             array_of_starts,
                             order, MPI_DOUBLE, &sub_send_up);
    MPI_Type_commit(&sub_send_up);

    // Down receive
    array_of_starts[0] = nolp;
    array_of_starts[1] = 0;
    MPI_Type_create_subarray(ndims,
                             array_of_sizes,
                             array_of_subsizes,
                             array_of_starts,
                             order, MPI_DOUBLE, &sub_recv_down);
    MPI_Type_commit(&sub_recv_down);

	// Down send
	array_of_starts[0] = nolp;
	array_of_starts[1] = nolp;
    MPI_Type_create_subarray(ndims,
                             array_of_sizes,
                             array_of_subsizes,
                             array_of_starts,
                             order, MPI_DOUBLE, &sub_send_down);
    MPI_Type_commit(&sub_send_down);

    // Up receive
    array_of_starts[0] = nolp;
    array_of_starts[1] = size_ud - 1;
    MPI_Type_create_subarray(ndims,
                             array_of_sizes,
                             array_of_subsizes,
                             array_of_starts,
                             order, MPI_DOUBLE, &sub_recv_up);
    MPI_Type_commit(&sub_recv_up);
};

// Post all send and receive messages. Note that this returns a requests array, this is manually 
// checked later.
void Communicate_Solution(MPI_Comm CART_COMM, double* w_ptr, MPI_Request* send_req, MPI_Request* recv_req){
	int tag = 22;
	MPI_Status status;

	// Send left, Receive right    
    MPI_Isend(w_ptr, 1, sub_send_left, left_neigh, tag, CART_COMM, &send_req[0]);
    MPI_Irecv(w_ptr, 1, sub_recv_right, right_neigh, tag, CART_COMM, &recv_req[0]);

	// Send right, Receive left    
    tag = 23;
    MPI_Isend(w_ptr, 1, sub_send_right, right_neigh, tag, CART_COMM, &send_req[1]);
    MPI_Irecv(w_ptr, 1, sub_recv_left, left_neigh, tag, CART_COMM, &recv_req[1]);

	// Send up, Receive down 
	tag = 24;   
    MPI_Isend(w_ptr, 1, sub_send_up, up_neigh, tag, CART_COMM, &send_req[2]);
    MPI_Irecv(w_ptr, 1, sub_recv_down, down_neigh, tag, CART_COMM, &recv_req[2]);

	// Send down, Receive up    
    tag = 25;
    MPI_Isend(w_ptr, 1, sub_send_down, down_neigh, tag, CART_COMM, &send_req[3]);
    MPI_Irecv(w_ptr, 1, sub_recv_up, up_neigh, tag, CART_COMM, &recv_req[3]);

};

// Output solution on each node to FileName
void Print_Solution(char* FileName,double* IO_buf, Darray2& w, int* size_array, int proc_per_node, int node_rank, MPI_Comm IO_Comm){
    MPI_Status status;
    int tag = 22;
    int source;
    int counter;

    if(node_rank == 0){
	    FILE *extFile = fopen(FileName, "w");

	    // First write own data to file
		for(int i=istart+1;i<=iend-1;i++){
			for(int j=jstart+1;j<=jend-1;j++){
				fprintf(extFile, "%18.10e\n", w(i,j));
			}
		}

		// Receive data from each process on the node and write to file.
        for(int i = 1;i<proc_per_node;i++){
        	source = i;
            MPI_Recv(IO_buf,size_array[i], MPI_DOUBLE, source, tag, IO_Comm, &status);
		    
		    for (int k=0; k<size_array[i]; k++){
		    	fprintf(extFile, "%18.10e\n", IO_buf[k]);
		    }
        }
    	fclose(extFile);
    }
    else{
    	counter = 0;
		for(int i=istart+1;i<=iend-1;i++){
			for(int j=jstart+1;j<=jend-1;j++){
				IO_buf[counter] = w(i,j);
				counter++;
			}
		}
        MPI_Send(IO_buf, counter, MPI_DOUBLE, 0, tag, IO_Comm);
    }
    MPI_Barrier(MPI_COMM_WORLD);
};


// Free all MPI_Datatypes needed for communication
void Finalize(){
	MPI_Type_free( &sub_send_left  ); 
	MPI_Type_free( &sub_recv_left  );
	MPI_Type_free( &sub_send_right ); 
	MPI_Type_free( &sub_recv_right );
	MPI_Type_free( &sub_send_up    );
	MPI_Type_free( &sub_recv_up    );
	MPI_Type_free( &sub_send_down  ); 
	MPI_Type_free( &sub_recv_down  );	
};

};
#endif
