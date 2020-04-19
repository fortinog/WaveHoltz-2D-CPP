#include "wave_solve.hpp"

Wave_Solve::Wave_Solve(Subdomain Local_Grid, int node_ID, MPI_Comm CART_COMM, MPI_Comm IO_Comm){
	int nsteps, nolp;
	int ndims = 2;
	int x_off, y_off;
	int rank;
	int* size_array;
	int* mask_buf;
	int* full_sz;
	int node_rank;
	int buf_size;
	int proc_per_node;
	int loc_sz[4];
	double lap_err[3];
	double sol_err[3];
	double err;
	double dt, dt2, idt2; 
    char fName[100];   
	FILE *extFile;
    MPI_Status status;
	Darray2 w, lap, wm, wp;

	t = 0.0;
    MPI_Comm_rank( CART_COMM, &rank );

	N 	   = Local_Grid.nx_loc;
	M 	   = Local_Grid.ny_loc;
	nolp   = Local_Grid.nolp;

	// MPI rank of neighbor in each direction
	up_neigh    = Local_Grid.up_neigh;
	down_neigh  = Local_Grid.down_neigh;
	left_neigh  = Local_Grid.left_neigh;
	right_neigh = Local_Grid.right_neigh;

	/* Remove a point from a process if it is has part of a physical boundary. This is 
	   to avoid assigning a ghost point when a process already has a boundary point and 
	   save on storage. */
	N = (left_neigh  == MPI_PROC_NULL) ? N-1 : N;
	N = (right_neigh == MPI_PROC_NULL) ? N-1 : N;
	M = (down_neigh  == MPI_PROC_NULL) ? M-1 : M;
	M = (up_neigh    == MPI_PROC_NULL) ? M-1 : M;

	istart = 1-nolp;
	iend   = N+nolp;
	jstart = 1-nolp;
	jend   = M+nolp;

	/* Start and stop indices along a boundary in each direction for a process. 
	   Set up so that we avoid doing logical checks each time 
	   we enforce boundary conditions and loops are empty if no physical boundary
	   present */
	j_left_start  = 0;
	j_left_end    = -1;
	j_right_start = 0;
	j_right_end   = -1;
	i_up_start    = 0;
	i_up_end      = -1;
	i_down_start  = 0;
	i_down_end    = -1;

	x_off = nolp;
	y_off = nolp;

	// Only assign sensible boundary loop indices if next to a boundary.
	if(left_neigh  == MPI_PROC_NULL){
		j_left_start = jstart+1;
		j_left_end   = jend-1;
		x_off = 0;
	}
	if(right_neigh  == MPI_PROC_NULL){
		j_right_start = jstart+1;
		j_right_end   = jend-1;
	}
	if(up_neigh  == MPI_PROC_NULL){
		i_up_start = istart+1;
		i_up_end   = iend-1;
	}
	if(down_neigh  == MPI_PROC_NULL){
		y_off = 0;
		i_down_start = istart+1;
		i_down_end   = iend-1;
	}

	// Time stepping info
	nsteps = setup.nsteps;
	dt     = setup.dt;
	dt2    = pow(setup.dt,2);
	idt2   = 1.0/dt2;

	x.define(1,istart,iend);
	y.define(1,jstart,jend);


	// Setup local grid: Note i=0,N+1, j=0,M+1 correspond to either the boundary
	// or a ghost region for the given subdomain.
	for(int i=istart;i<=iend;i++) x(i) = setup.x_L + (i + Local_Grid.x_s - x_off)*setup.hx;
	for(int i=jstart;i<=jend;i++) y(i) = setup.y_L + (i + Local_Grid.y_s - y_off)*setup.hy;


	// Initialize working arrays
	w.define(1,istart,iend,jstart,jend);
	w.set_value(0.0);
	wp.define(1,istart,iend,jstart,jend);
	wp.set_value(0.0);
	wm.define(1,istart,iend,jstart,jend);
	wm.set_value(0.0);
	lap.define(1,istart,iend,jstart,jend);
	lap.set_value(0.0);
	mask.define(1,istart,iend,jstart,jend);
	mask.set_value(10);


	// Bring out data pointer for message passing
    double* w_ptr = w.c_ptr();

    // Length of data to be sent and received left/right or up/down
    size_lr = iend-istart+1;
    size_ud = jend-jstart+1;
    buf_size = (size_lr-2)*(size_ud-2);
    loc_sz[0] = (size_lr-2);
    loc_sz[1] = (size_ud-2);
    loc_sz[2] = istart + 1 + Local_Grid.x_s - x_off;
    loc_sz[3] = jstart + 1 + Local_Grid.y_s - y_off;

    // For local compute node, find my rank and number of processes
    MPI_Comm_rank(IO_Comm, &node_rank);
    MPI_Comm_size(IO_Comm, &proc_per_node);

    // Grab solution array sizes from each process on the node
    if(node_rank == 0){
    	size_array = new int[proc_per_node];
    	full_sz = new int[4*proc_per_node];
    }
    MPI_Gather(loc_sz, 4, MPI_INT, full_sz, 4, MPI_INT, 0, IO_Comm);

    // Print node data to file
    if(node_rank == 0){
    	sprintf(fName, "Node_Info_%4.4i.txt", node_ID);
    	extFile = fopen(fName, "w");
	    for(int i = 0;i<proc_per_node;i++){
	    	size_array[i] = full_sz[4*i]*full_sz[4*i+1];
			fprintf(extFile, "%1i %1i %1i %1i\n", full_sz[4*i], full_sz[4*i+1], full_sz[4*i+2], full_sz[4*i+3]);
	    }
    	fclose(extFile);
    }
    MPI_Barrier(IO_Comm);


    // If I am the root process on a node allocate a receive buffer for IO
    // large enough for all data, otherwise allocate just enough to copy
    // local data sans halo region.
    if(node_rank == 0){
    	int* max_sz = max_element(size_array, size_array+proc_per_node);
    	IO_buf = new double[*max_sz];
    	mask_buf = new int[*max_sz];
    } else{
    	IO_buf = new double[buf_size];
    	mask_buf = new int[buf_size];
    }


    /******************************************************************
	* 				    Compute and Print mask 						  *
    ******************************************************************/
    Compute_Mask();
    if(node_rank == 0){
    	sprintf(fName, "mask_%4.4i.txt", node_ID);

	    extFile = fopen(fName, "w");

	    // First write own data to file
		for(int i=istart+1;i<=iend-1;i++){
			for(int j=jstart+1;j<=jend-1;j++){
				fprintf(extFile, "%1i\n", mask(i,j));
			}
		}

		// Receive data from each process on the node and write to file.
        for(int i = 1;i<proc_per_node;i++){
        	int source = i;
            MPI_Recv(mask_buf,size_array[i], MPI_INT, source, 420, IO_Comm, &status);
		    
		    for (int k=0; k<size_array[i]; k++){
		    	fprintf(extFile, "%1i\n", mask_buf[k]);
		    }
        }
    	fclose(extFile);
    }
    else{
    	int counter = 0;
		for(int i=istart+1;i<=iend-1;i++){
			for(int j=jstart+1;j<=jend-1;j++){
				mask_buf[counter] = mask(i,j);
				counter++;
			}
		}
        MPI_Send(mask_buf, counter, MPI_INT, 0, 420, IO_Comm);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    delete[] mask_buf;
    

    // Set up all communication subarrays.
    Setup_Subarrays(nolp);

	// Fill in w, wm, lap so that we can start time stepping immediately
	Set_Initial_Data(wm, w, lap);

	// Solve PDE forward to final time specified in Problem_Setup
    Solve_PDE(wm, w, wp, lap, w_ptr, CART_COMM);

    // Compute solution error (1, 2, Inf norms)
    sol_err[0] = Compute_Solution_Error(1, w, CART_COMM);
    sol_err[1] = Compute_Solution_Error(2, w, CART_COMM);
    sol_err[2] = Compute_Solution_Error(3, w, CART_COMM);

    // Compute Laplacian error (1, 2, Inf norms)
    lap_err[0] = Compute_Laplacian_Error(1, lap, CART_COMM);
    lap_err[1] = Compute_Laplacian_Error(2, lap, CART_COMM);
    lap_err[2] = Compute_Laplacian_Error(3, lap, CART_COMM);

	MPI_Barrier(CART_COMM);
    if(rank == 0){
    	cout << "Error in computing the Laplacian : ";
		cout << scientific << right << setw(14)<< lap_err[0];
		cout << scientific << right << setw(14)<< lap_err[1];
		cout << scientific << right << setw(14)<< lap_err[2];
		cout << "\n";

    	cout << "Error in computing the solution : ";
		cout << scientific << right << setw(14)<< sol_err[0];
		cout << scientific << right << setw(14)<< sol_err[1];
		cout << scientific << right << setw(14)<< sol_err[2];
		cout << "\n";
    }
	MPI_Barrier(CART_COMM);

	// Output the solution to a text file.
	Print_Solution("out.txt", IO_buf, w, size_array, proc_per_node, node_rank, IO_Comm);
	Finalize();
}

// Compute the mask for all points in the domain via 
// the level set function (FMG fix)
void Wave_Solve::Compute_Mask(){

     double level_set_val;
     double distance;
     double x_val;
     int ghost_ctr = 0;
     int num_int = 0;
     int bdry_ctr = 0;
     int ctr =0;
     int n_left = 0;
     int n_right = 0;
     int n_up = 0;
     int n_down = 0;
     int dir_ctr[4];

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
                 bdry_ctr++;
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
    ghost_list_left.define(1,1,n_left,1,2);
    dist_list_left.define(1,1,n_left);
    ghost_list_right.define(1,1,n_right,1,2);
    dist_list_right.define(1,1,n_right);
    ghost_list_up.define(1,1,n_up,1,2);
    dist_list_up.define(1,1,n_up);
    ghost_list_down.define(1,1,n_down,1,2);
    dist_list_down.define(1,1,n_down);

    for(int i = 0;i<4;i++) dir_ctr[i] = 0;

    for(int i = istart+1;i<=iend-1;i++){
        for(int j = jstart+1;j<=jend-1;j++){
        	if(mask(i,j) == -1){

        		// Right
        		if(mask(i+1,j) == -3){
        			ghost_list_right(dir_ctr[0],1) = i;
        			ghost_list_right(dir_ctr[0],2) = j;
        			dist_list_right(dir_ctr[0]) = setup.Dist_to_bdry(x(i), x(i+1), y(j), y(j),1)/setup.hx;
        			dir_ctr[0]++;
        		}

        		// Up
        		if(mask(i,j+1) == -3){
        			ghost_list_up(dir_ctr[1],1) = i;
        			ghost_list_up(dir_ctr[1],2) = j;
        			dist_list_up(dir_ctr[1]) = setup.Dist_to_bdry(x(i), x(i), y(j), y(j+1),2)/setup.hy;
        			dir_ctr[1]++;
        		}

        		// Left
        		if(mask(i-1,j) == -3){
        			ghost_list_left(dir_ctr[2],1) = i;
        			ghost_list_left(dir_ctr[2],2) = j;
        			dist_list_left(dir_ctr[2]) = setup.Dist_to_bdry(x(i), x(i-1), y(j), y(j),1)/setup.hx;
        			dir_ctr[2]++;
        		}

        		// Down
        		if(mask(i,j-1) == -3){
        			ghost_list_down(dir_ctr[3],1) = i;
        			ghost_list_down(dir_ctr[3],2) = j;
        			dist_list_down(dir_ctr[3]) = setup.Dist_to_bdry(x(i), x(i), y(j), y(j-1),2)/setup.hy;
        			dir_ctr[3]++;
        		}
        	} 
        }
    }
    
}

// Fill in with initial data, compute the Laplacian, and do a 
// Taylor expansion so that we may start time stepping after this.
void Wave_Solve::Set_Initial_Data(Darray2& wm, Darray2& w, Darray2& lap){
	Twilight mms(setup);
	double x0;

	for(int i=istart+1;i<=iend-1;i++){
		x0 = x(i);
		for(int j=jstart+1;j<=jend-1;j++){
			w(i,j) = mms.trigTwilight(0,x0,0,y(j),0,0.0);
		}
	}
	Enforce_BC(w);
	Compute_Laplacian(w, lap);
	Taylor_Expand(wm,w,lap);
}

// Forcing (method of manufactured solutions)
double Wave_Solve::forcing(const double x0, const double y0, const double t0){
	Twilight mms(setup);
	double wxx, wyy, val;
	wxx = mms.trigTwilight(2,x0,0,y0,0,t0);
	wyy = mms.trigTwilight(0,x0,2,y0,0,t0);
	val = mms.trigTwilight(0,x0,0,y0,2,t0) - wxx - wyy;
	return val;
}


// Enforce the boundary conditions. In this simple case, fill in
// the boundary data.
void Wave_Solve::Enforce_BC(Darray2 &v){
	Twilight mms(setup);
	int i,j;
	double x0, y0;
	double t = (*this).t;

	// Left boundary
	i = istart;
	x0 = x(i);
	for(j=j_left_start;j<=j_left_end;j++){
		v(i,j) = mms.trigTwilight(0,x0,0,y(j),0,t);
	}

	// Right boundary
	i = iend;
	x0 = x(i);
	for(j=j_right_start;j<=j_right_end;j++){
		v(i,j) = mms.trigTwilight(0,x0,0,y(j),0,t);
	}

	// Bottom boundary
	j = jstart;
	y0 = y(j);
	for(i=i_down_start;i<=i_down_end;i++){
		v(i,j) = mms.trigTwilight(0,x(i),0,y0,0,t);
	}

	// Top boundary
	j = jend;
	y0 = y(j);
	for(i=i_up_start;i<=i_up_end;i++){
		v(i,j) = mms.trigTwilight(0,x(i),0,y0,0,t);
	}
}


void Wave_Solve::Compute_Laplacian(Darray2& w, Darray2& lap){
	double ihx2 = 1.0/pow(setup.hx,2);
	double ihy2 = 1.0/pow(setup.hy,2);
	double c1,c2,c3;
	double x0,xp,xm;
	double y0,yp,ym;

	// Compute the usual Laplacian (FULL)
	for(int i=istart+1;i<=iend-1;i++){
		xm = x(i-1);
		x0 = x(i);
		xp = x(i+1);
		for(int j=jstart+1;j<=jend-1;j++){
			ym = y(j-1);
			y0 = y(j);
			yp = y(j+1);

			// Compute d_{yy}
			c1 = setup.c2(x0,ym);
			c2 = setup.c2(x0,y0);
			c3 = setup.c2(x0,yp);
			lap(i,j) = 0.5*ihy2*((c2+c3)*w(i,j+1) - (c1+2.0*c2+c3)*w(i,j) + (c1+c2)*w(i,j-1));

			// Compute d_{xx}
			c1 = setup.c2(xm,y0);
			c2 = setup.c2(x0,y0);
			c3 = setup.c2(xp,y0);

			// Add all together
			lap(i,j) += 0.5*ihx2*((c2+c3)*w(i+1,j) - (c1+2.0*c2+c3)*w(i,j) + (c1+c2)*w(i-1,j));
		}
	}
   

}
// This version of the routine computes the Laplacian on the innermost grid points first 
// while communication is happening. Once all data has been communicated, compute the Laplacian
// on the edges of current computational subdomain.
void Wave_Solve::Compute_Laplacian_NB(Darray2& w, Darray2& lap, MPI_Request* recv_req){
	double ihx2 = 1.0/pow(setup.hx,2);
	double ihy2 = 1.0/pow(setup.hy,2);
	double c1,c2,c3;
	double x0,xp,xm;
	double y0,yp,ym;
	int i,j;

	// Compute the usual Laplacian (inner only)
	for(i=istart+2;i<=iend-2;i++){
		xm = x(i-1);
		x0 = x(i);
		xp = x(i+1);
		for(j=jstart+2;j<=jend-2;j++){
			ym = y(j-1);
			y0 = y(j);
			yp = y(j+1);

			// Compute d_{yy}
			c1 = setup.c2(x0,ym);
			c2 = setup.c2(x0,y0);
			c3 = setup.c2(x0,yp);
			lap(i,j) = 0.5*ihy2*((c2+c3)*w(i,j+1) - (c1+2.0*c2+c3)*w(i,j) + (c1+c2)*w(i,j-1));

			// Compute d_{xx}
			c1 = setup.c2(xm,y0);
			c2 = setup.c2(x0,y0);
			c3 = setup.c2(xp,y0);

			// Add all together
			lap(i,j) += 0.5*ihx2*((c2+c3)*w(i+1,j) - (c1+2.0*c2+c3)*w(i,j) + (c1+c2)*w(i-1,j));
		}
	}
    
    // Once all receives are processed, compute Laplacian on the 
    // edges of the current process subdomain
    MPI_Waitall( 4, recv_req, MPI_STATUSES_IGNORE );

    // bottom
	j = jstart+1;
	ym = y(j-1);
	y0 = y(j);
	yp = y(j+1);
	for(i=istart+1;i<=iend-1;i++){
		xm = x(i-1);
		x0 = x(i);
		xp = x(i+1);

		// Compute d_{yy}
		c1 = setup.c2(x0,ym);
		c2 = setup.c2(x0,y0);
		c3 = setup.c2(x0,yp);
		lap(i,j) = 0.5*ihy2*((c2+c3)*w(i,j+1) - (c1+2.0*c2+c3)*w(i,j) + (c1+c2)*w(i,j-1));

		// Compute d_{xx}
		c1 = setup.c2(xm,y0);
		c2 = setup.c2(x0,y0);
		c3 = setup.c2(xp,y0);

		// Add all together
		lap(i,j) += 0.5*ihx2*((c2+c3)*w(i+1,j) - (c1+2.0*c2+c3)*w(i,j) + (c1+c2)*w(i-1,j));
	}

	// top
	j = jend-1;
	ym = y(j-1);
	y0 = y(j);
	yp = y(j+1);
	for(i=istart+1;i<=iend-1;i++){
		xm = x(i-1);
		x0 = x(i);
		xp = x(i+1);

		// Compute d_{yy}
		c1 = setup.c2(x0,ym);
		c2 = setup.c2(x0,y0);
		c3 = setup.c2(x0,yp);
		lap(i,j) = 0.5*ihy2*((c2+c3)*w(i,j+1) - (c1+2.0*c2+c3)*w(i,j) + (c1+c2)*w(i,j-1));

		// Compute d_{xx}
		c1 = setup.c2(xm,y0);
		c2 = setup.c2(x0,y0);
		c3 = setup.c2(xp,y0);

		// Add all together
		lap(i,j) += 0.5*ihx2*((c2+c3)*w(i+1,j) - (c1+2.0*c2+c3)*w(i,j) + (c1+c2)*w(i-1,j));
	}


	// left
	i = istart + 1;
	xm = x(i-1);
	x0 = x(i);
	xp = x(i+1);
	for(j=jstart+1;j<=jend-1;j++){

		ym = y(j-1);
		y0 = y(j);
		yp = y(j+1);

		// Compute d_{yy}
		c1 = setup.c2(x0,ym);
		c2 = setup.c2(x0,y0);
		c3 = setup.c2(x0,yp);
		lap(i,j) = 0.5*ihy2*((c2+c3)*w(i,j+1) - (c1+2.0*c2+c3)*w(i,j) + (c1+c2)*w(i,j-1));

		// Compute d_{xx}
		c1 = setup.c2(xm,y0);
		c2 = setup.c2(x0,y0);
		c3 = setup.c2(xp,y0);

		// Add all together
		lap(i,j) += 0.5*ihx2*((c2+c3)*w(i+1,j) - (c1+2.0*c2+c3)*w(i,j) + (c1+c2)*w(i-1,j));
	}

	// right
	i = iend - 1;
	xm = x(i-1);
	x0 = x(i);
	xp = x(i+1);
	for(j=jstart+1;j<=jend-1;j++){

		ym = y(j-1);
		y0 = y(j);
		yp = y(j+1);

		// Compute d_{yy}
		c1 = setup.c2(x0,ym);
		c2 = setup.c2(x0,y0);
		c3 = setup.c2(x0,yp);
		lap(i,j) = 0.5*ihy2*((c2+c3)*w(i,j+1) - (c1+2.0*c2+c3)*w(i,j) + (c1+c2)*w(i,j-1));

		// Compute d_{xx}
		c1 = setup.c2(xm,y0);
		c2 = setup.c2(x0,y0);
		c3 = setup.c2(xp,y0);

		// Add all together
		lap(i,j) += 0.5*ihx2*((c2+c3)*w(i+1,j) - (c1+2.0*c2+c3)*w(i,j) + (c1+c2)*w(i-1,j));
	}

}

// Only in the constant coefficient case, flag=1 L1 norm, flag=2 L2 norm, flag=3 L_inf norm
double Wave_Solve::Compute_Laplacian_Error(const int flag, Darray2& lap, MPI_Comm CART_COMM){
	double wxx,true_lap;
	double err = 0.0;
	double global_err;
	double norm_sol= 0.0;
	double global_norm;
	double x0;
	double max_val = 0.0;
	double max_sol = 0.0;
	double t = (*this).t;
	Twilight mms(setup);

	switch(flag){
		case 1:
			for(int i=istart+1;i<=iend-1;i++){
				x0 = x(i);
				for(int j=jstart+1;j<=jend-1;j++){
					wxx = mms.trigTwilight(2,x0,0,y(j),0,t);
					true_lap = wxx + mms.trigTwilight(0,x0,2,y(j),0,t);
					err += abs(lap(i,j)-true_lap);
					norm_sol += abs(true_lap);
				}
			}
			MPI_Allreduce(&err, &global_err, 1, MPI_DOUBLE, MPI_SUM, CART_COMM);
			MPI_Allreduce(&norm_sol, &global_norm, 1, MPI_DOUBLE, MPI_SUM, CART_COMM);
			err = global_err/global_norm;
			break;
		case 2:
			for(int i=istart+1;i<=iend-1;i++){
				x0 = x(i);
				for(int j=jstart+1;j<=jend-1;j++){
					wxx = mms.trigTwilight(2,x0,0,y(j),0,t);
					true_lap = wxx + mms.trigTwilight(0,x0,2,y(j),0,t);
					err += pow(lap(i,j)-true_lap,2);
					norm_sol += pow(true_lap,2);
				}
			}

			MPI_Allreduce(&err, &global_err, 1, MPI_DOUBLE, MPI_SUM, CART_COMM);
			MPI_Allreduce(&norm_sol, &global_norm, 1, MPI_DOUBLE, MPI_SUM, CART_COMM);
			err = sqrt(global_err/global_norm);
			break;
		case 3:
			for(int i=istart+1;i<=iend-1;i++){
				x0 = x(i);
				for(int j=jstart+1;j<=jend-1;j++){
					wxx = mms.trigTwilight(2,x0,0,y(j),0,t);
					true_lap = wxx + mms.trigTwilight(0,x0,2,y(j),0,t);
					if(abs(lap(i,j)-true_lap) > max_val) max_val = abs(lap(i,j)-true_lap);
					if(abs(true_lap) > max_sol) max_sol = abs(true_lap);
				}
			}
			MPI_Allreduce(&max_val, &global_err, 1, MPI_DOUBLE, MPI_MAX, CART_COMM);
			MPI_Allreduce(&max_sol, &global_norm, 1, MPI_DOUBLE, MPI_MAX, CART_COMM);
			err = global_err/global_norm;
			break;
	}

	return err;
}

// Only in the constant coefficient case, flag=1 L1 norm, flag=2 L2 norm, flag=3 L_inf norm
double Wave_Solve::Compute_Solution_Error(const int flag, Darray2& w, MPI_Comm CART_COMM){
	double true_sol;
	double err = 0.0;
	double norm_sol= 0.0;
	double global_err;
	double global_norm;
	double x0;
	double max_val = 0.0;
	double max_sol = 0.0;
	double t = (*this).t;
	Twilight mms(setup);

	switch(flag){
		case 1:
			for(int i=istart+1;i<=iend-1;i++){
				x0 = x(i);
				for(int j=jstart+1;j<=jend-1;j++){
					true_sol = mms.trigTwilight(0,x0,0,y(j),0,t);
					err += abs(w(i,j)-true_sol);
					norm_sol += abs(true_sol);
				}
			}
			MPI_Allreduce(&err, &global_err, 1, MPI_DOUBLE, MPI_SUM, CART_COMM);
			MPI_Allreduce(&norm_sol, &global_norm, 1, MPI_DOUBLE, MPI_SUM, CART_COMM);
			err = global_err/global_norm;
			break;
		case 2:
			for(int i=istart+1;i<=iend-1;i++){
				x0 = x(i);
				for(int j=jstart+1;j<=jend-1;j++){
					true_sol = mms.trigTwilight(0,x0,0,y(j),0,t);
					err += pow(w(i,j)-true_sol,2);
					norm_sol += pow(true_sol,2);
				}
			}
			MPI_Allreduce(&err, &global_err, 1, MPI_DOUBLE, MPI_SUM, CART_COMM);
			MPI_Allreduce(&norm_sol, &global_norm, 1, MPI_DOUBLE, MPI_SUM, CART_COMM);
			err = sqrt(global_err/global_norm);
			break;
		case 3:
			for(int i=istart+1;i<=iend-1;i++){
				x0 = x(i);
				for(int j=jstart+1;j<=jend-1;j++){
					true_sol = mms.trigTwilight(0,x0,0,y(j),0,t);
					if(abs(w(i,j)-true_sol) > max_val) max_val = abs(w(i,j)-true_sol);
					if(abs(true_sol) > max_sol) max_sol = abs(true_sol);
				}
			}
			MPI_Allreduce(&max_val, &global_err, 1, MPI_DOUBLE, MPI_MAX, CART_COMM);
			MPI_Allreduce(&max_sol, &global_norm, 1, MPI_DOUBLE, MPI_MAX, CART_COMM);
			err = global_err/global_norm;
			break;
	}

	return err;
}

//	Do a Taylor expansion to fill in w_{-1} \approx w(-dt).
void Wave_Solve::Taylor_Expand(Darray2& wm, Darray2& w, Darray2& lap){
	Twilight mms(setup);
	double dt = setup.dt;
	double dt2 = setup.dt2;
	double x0;
	for(int i=istart+1;i<=iend-1;i++){
		x0 = x(i);
		for(int j=jstart+1;j<=jend-1;j++){
			wm(i,j) = w(i,j) - dt*mms.trigTwilight(0,x0,0,y(j),1,0.0) + 0.5*dt2*(lap(i,j) + forcing(x0,y(j),0.0));
			// wm(i,j) = w(i,j) + 0.5*dt2*(lap(i,j) + forcing(x0,y(j),0.0));
		}
	}
}

// Perform a single time step forward, centered difference in time
void Wave_Solve::Time_Step(Darray2& wm, Darray2& w, Darray2& wp, Darray2& lap){
	int N = setup.N;
	int M = setup.M;
	double dt2 = setup.dt2;
	double x0;
	double t = (*this).t;

	for(int i=istart+1;i<=iend-1;i++){
		x0 = x(i);
		for(int j=jstart+1;j<=jend-1;j++){
			wp(i,j) = 2.0*w(i,j) - wm(i,j) + dt2*(lap(i,j) + forcing(x0,y(j),t)); 
		}
	}
	Enforce_BC(wp);

	wm.copy(w);
	w.copy(wp);
	(*this).t += setup.dt;
}

// Compute the energy at a given time step. (Calculation from Patrick Joly's paper)
double Wave_Solve::Compute_Energy(Darray2& wm, Darray2& w, Darray2& lap, MPI_Comm CART_COMM){
	double dt2 = setup.dt2;
	double energy = 0.0;
	double total_energy;
	double idt2 = setup.idt2;

	for(int i=istart+1;i<=iend-1;i++){
		for(int j=jstart+1;j<=jend-1;j++){
			energy = energy + idt2*pow(wm(i,j) - w(i,j),2) - w(i,j)*lap(i,j);
		}
	}
	MPI_Allreduce(&energy, &total_energy, 1, MPI_DOUBLE, MPI_SUM, CART_COMM);

	return 0.5*setup.hx*setup.hy*total_energy;
}

// Do the full wave solve. Note: This currently prints out the 
// energy for every time step.
void Wave_Solve::Solve_PDE(Darray2& wm, Darray2& w, Darray2& wp, Darray2& lap, double* w_ptr, MPI_Comm CART_COMM){
	double energy_old = 1e10;
	double energy;
    MPI_Request send_req[4];
    MPI_Request recv_req[4];

	for(int i =0;i<setup.nsteps;i++){
	    Time_Step(wm,w,wp,lap);
	    energy = Compute_Energy(wm, w, lap,CART_COMM);
        cout << scientific << right << setw(14)<< abs(energy-energy_old)/abs(energy_old) << "\n";
	    energy_old = energy;
	    Communicate_Solution(CART_COMM,w_ptr,send_req,recv_req);
	    Compute_Laplacian_NB(w, lap, recv_req);
        MPI_Waitall( 4, send_req, MPI_STATUSES_IGNORE );
	}
}

// This routine sets up the communication subarrays to simplify MPI calls.
void Wave_Solve::Setup_Subarrays(const int nolp){
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

}

// Communicate all ghost points. Non-Blocking communication used.
void Wave_Solve::Communicate_Solution(MPI_Comm CART_COMM, double* w_ptr, MPI_Request* send_req, MPI_Request* recv_req){
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
}


// Communicate data to root process of each node and then write to a file.
void Wave_Solve::Print_Solution(char* FileName, double* IO_buf, Darray2& w, int* size_array, 
								int proc_per_node, 
								int node_rank, 
								MPI_Comm IO_Comm){
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
}

// Free MPI communication types
void Wave_Solve::Finalize(){
	MPI_Type_free( &sub_send_left  ); 
	MPI_Type_free( &sub_recv_left  );
	MPI_Type_free( &sub_send_right ); 
	MPI_Type_free( &sub_recv_right );
	MPI_Type_free( &sub_send_up    );
	MPI_Type_free( &sub_recv_up    );
	MPI_Type_free( &sub_send_down  ); 
	MPI_Type_free( &sub_recv_down  );
}