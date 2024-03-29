#include "wave_solve.hpp"
#include <algorithm>
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
	int WHiter= 0;
	double lap_err[3];
	double sol_err[3];
	double err;
	double dt, dt2, idt2; 
	double res, global_res, global_res_old;
	double alpha, beta, tmp;
	double res_delta;
    char fName[100];
    MPI_Request send_req[4];
    MPI_Request recv_req[4];
	Twilight mms(setup);      
	FILE *extFile;
    MPI_Status status;
	Darray2 w, lap, wm, wp, u, uold, b, r, u_sol;

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
	uold.define(1,istart,iend,jstart,jend);
    u.define(1,istart,iend,jstart,jend);
    u.set_value(0.0);
	mask.define(1,istart,iend,jstart,jend);
	mask.set_value(10);
	b.define(1,istart,iend,jstart,jend);
   	r.define(1,istart,iend,jstart,jend);
	u_sol.define(1,istart,iend,jstart,jend);
	u_sol.set_value(0.0);

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
    	int* max_sz = std::max_element(size_array, size_array+proc_per_node);
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

	// Get RHS of SPD WaveHoltz system (stored in b)
	Set_Initial_Data(wm, b, lap, u);
	Solve_PDE(wm, w, wp, lap, u, w_ptr, CART_COMM);
	b.copy(u);
	(*this).t = 0;

	// Fill in w, wm, lap so that we can start time stepping immediately
	Set_Initial_Data(wm, w, lap, u);
    Evolve_and_Project(b, uold, wm, w, wp, lap, u, w_ptr, CART_COMM);

	// Initialize residual and direction vector for Conjugate Gradient
    r.copy(b);
    res = Compute_Inner_Product(r, r, CART_COMM);

	// Perform standard Conjugate Gradient
    while((WHiter<setup.cg_iter) && (res > setup.cg_tol)){  
    	Evolve_and_Project(b, uold, wm, w, wp, lap, u, w_ptr, CART_COMM);
		alpha = res/Compute_Inner_Product(uold, u, CART_COMM);
		for(int i=istart+1;i<=iend-1;i++){
			for(int j=jstart+1;j<=jend-1;j++){
				r(i,j) -= alpha*u(i,j);
				u_sol(i,j) += alpha*uold(i,j);
			}
		}
		tmp = Compute_Inner_Product(r, r, CART_COMM);
		beta = tmp/res;
		res = tmp;
		for(int i=istart+1;i<=iend-1;i++){
			for(int j=jstart+1;j<=jend-1;j++){
				w(i,j) = r(i,j) + beta*uold(i,j);
			}
		}

		// Print convergence history to screen
		if(rank == 0){
			fprintf(stderr, "ITER : %3.3i\t RESIDUAL: %4.4e \n", WHiter, res);
		}
		WHiter++;

		if(WHiter%setup.print_freq == 0){
	    	sprintf(fName, "data/u_%4.4i_%4.4i.txt", node_ID,WHiter);
			Print_Solution(fName, IO_buf, u_sol, size_array, proc_per_node, node_rank, IO_Comm);
		}

		// Print solution to file
	}
	w.copy(u_sol);
	Enforce_BC(w);

    // Compute solution error (1, 2, Inf norms)
    sol_err[0] = Compute_Solution_Error(1, w, CART_COMM);
    sol_err[1] = Compute_Solution_Error(2, w, CART_COMM);
    sol_err[2] = Compute_Solution_Error(3, w, CART_COMM);

	MPI_Barrier(CART_COMM);
    if(rank == 0){
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


// Fill in with initial data, compute the Laplacian, and do a 
// Taylor expansion so that we may start time stepping after this.
void Wave_Solve::Set_Initial_Data(Darray2& wm, Darray2& w, Darray2& lap, Darray2& u){
	Twilight mms(setup);
	double x0;

	for(int i=istart;i<=iend;i++){
		x0 = x(i);
		for(int j=jstart;j<=jend;j++){
            w(i,j) = 0.0;
            u(i,j) = w(i,j)*0.5*setup.dt*0.75;
		}
	}
	Compute_Laplacian(w, lap);
	// Enforce_BC(w);
	Taylor_Expand(wm,w,lap);
}

// Forcing (method of manufactured solutions)
double Wave_Solve::forcing(const double x0, const double y0, const double t0){
	Twilight mms(setup);
	double wxx, wyy, val;
    double w;
	wxx = mms.trigTwilight(2,x0,0,y0,0,t0);
	wyy = mms.trigTwilight(0,x0,2,y0,0,t0);
    w = mms.trigTwilight(0,x0,0,y0,0,t0);
	val =  -(wxx + wyy+ pow(setup.omega,2)*w)*cos(setup.omega*t0);
	// val =  wxx + wyy - mms.trigTwilight(0,x0,0,y0,2,t0);
	return val;
}


// Enforce the boundary conditions. FMG: Dirichlet boundary conditions only.
void Wave_Solve::Enforce_BC(Darray2 &v){
	Twilight mms(setup);
	int i,j;
	double x0, y0;
	double alpha;
	double bc_val;


	// If a node point is a physical boundary, fill in value for 
	// Dirichlet boundary condition.
	for(int k = 1;k<=n_bdry;k++){
		i = bdry_list(k,1);
		j = bdry_list(k,2);
		v(i,j) =  cos(setup.omega*(*this).t)*(mms.trigTwilight(0,x(i),0,y(j),0,(*this).t));
	}

	// Boundary is right
	for(int k=1;k<=n_ghost_right;k++){
		alpha = dist_list_right(k);
		i = ghost_list_right(k,1);
		j = ghost_list_right(k,2);
		bc_val =  cos(setup.omega*(*this).t)*(mms.trigTwilight(0,x(i) + alpha*setup.hx,0,y(j),0,(*this).t));
		v(i,j) = (bc_val + alpha*v(i-1,j))/(1.0+alpha);
	}

	// Boundary is left
	for(int k=1;k<=n_ghost_left;k++){
		alpha = dist_list_left(k);
		i = ghost_list_left(k,1);
		j = ghost_list_left(k,2);
		bc_val =  cos(setup.omega*(*this).t)*(mms.trigTwilight(0,x(i) - alpha*setup.hx,0,y(j),0,(*this).t));
		v(i,j) = (bc_val + alpha*v(i+1,j))/(1.0+alpha);
	}

	// Boundary is above
	for(int k=1;k<=n_ghost_up;k++){
		alpha = dist_list_up(k);
		i = ghost_list_up(k,1);
		j = ghost_list_up(k,2);
		bc_val =  cos(setup.omega*(*this).t)*(mms.trigTwilight(0,x(i),0,y(j) + alpha*setup.hy,0,(*this).t));
		v(i,j) = (bc_val + alpha*v(i,j-1))/(1.0+alpha);
	}

	// Boundary is below
	for(int k=1;k<=n_ghost_down;k++){
		alpha = dist_list_down(k);
		i = ghost_list_down(k,1);
		j = ghost_list_down(k,2);
		bc_val =  cos(setup.omega*(*this).t)*(mms.trigTwilight(0,x(i),0,y(j) - alpha*setup.hy,0,(*this).t));
		v(i,j) = (bc_val + alpha*v(i,j+1))/(1.0+alpha);
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
	// double t = (*this).t;
	Twilight mms(setup);

	switch(flag){
		case 1:
			for(int i=istart+1;i<=iend-1;i++){
				x0 = x(i);
				for(int j=jstart+1;j<=jend-1;j++){
					if(mask(i,j) == 1){
						wxx = mms.trigTwilight(2,x0,0,y(j),0,(*this).t);
						true_lap = wxx + mms.trigTwilight(0,x0,2,y(j),0,(*this).t);
						err += abs(lap(i,j)-true_lap);
						norm_sol += abs(true_lap);
					}
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
					if(mask(i,j) == 1){
						wxx = mms.trigTwilight(2,x0,0,y(j),0,(*this).t);
						true_lap = wxx + mms.trigTwilight(0,x0,2,y(j),0,(*this).t);
						err += pow(lap(i,j)-true_lap,2);
						norm_sol += pow(true_lap,2);
					}
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
					if(mask(i,j) == 1){
						wxx = mms.trigTwilight(2,x0,0,y(j),0,(*this).t);
						true_lap = wxx + mms.trigTwilight(0,x0,2,y(j),0,(*this).t);
						if(abs(lap(i,j)-true_lap) > max_val) max_val = abs(lap(i,j)-true_lap);
						if(abs(true_lap) > max_sol) max_sol = abs(true_lap);		
					}
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
	Twilight mms(setup);

	switch(flag){
		// 1-norm
		case 1:
			for(int i=istart+1;i<=iend-1;i++){
				x0 = x(i);
				for(int j=jstart+1;j<=jend-1;j++){
					if((mask(i,j)) == 1){
						true_sol = mms.trigTwilight(0,x0,0,y(j),0,(*this).t);
						err += abs(w(i,j)-true_sol);
						norm_sol += abs(true_sol);
					}
				}
			}
			MPI_Allreduce(&err, &global_err, 1, MPI_DOUBLE, MPI_SUM, CART_COMM);
			MPI_Allreduce(&norm_sol, &global_norm, 1, MPI_DOUBLE, MPI_SUM, CART_COMM);
			err = global_err/global_norm;
			break;
		// 2-norm
		case 2:
			for(int i=istart+1;i<=iend-1;i++){
				x0 = x(i);
				for(int j=jstart+1;j<=jend-1;j++){
					if((mask(i,j)) == 1){
						true_sol = mms.trigTwilight(0,x0,0,y(j),0,(*this).t);
						err += pow(w(i,j)-true_sol,2);
						norm_sol += pow(true_sol,2);
					}
				}
			}
			MPI_Allreduce(&err, &global_err, 1, MPI_DOUBLE, MPI_SUM, CART_COMM);
			MPI_Allreduce(&norm_sol, &global_norm, 1, MPI_DOUBLE, MPI_SUM, CART_COMM);
			err = sqrt(global_err/global_norm);
			break;
		// Inf-norm
		case 3:
			for(int i=istart+1;i<=iend-1;i++){
				x0 = x(i);
				for(int j=jstart+1;j<=jend-1;j++){
					if((mask(i,j)) == 1){
						true_sol = mms.trigTwilight(0,x0,0,y(j),0,(*this).t);
						if(abs(w(i,j)-true_sol) > max_val) max_val = abs(w(i,j)-true_sol);
						if(abs(true_sol) > max_sol) max_sol = abs(true_sol);
					}
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
			wm(i,j) = w(i,j) + 0.5*dt2*(lap(i,j) + forcing(x0,y(j),0.0));
		}
	}
}

// Perform a single time step forward, centered difference in time
void Wave_Solve::Time_Step(Darray2& wm, Darray2& w, Darray2& wp, Darray2& lap){
	int N = setup.N;
	int M = setup.M;
	double dt2 = setup.dt2;
	double x0;
	for(int i=istart+1;i<=iend-1;i++){
		x0 = x(i);
		for(int j=jstart+1;j<=jend-1;j++){
			wp(i,j) = 2.0*w(i,j) - wm(i,j) + dt2*(lap(i,j) + forcing(x0,y(j),(*this).t)); 
		}
	}
	wm.copy(w);
	w.copy(wp);
	(*this).t += setup.dt;
	Enforce_BC(w);
}

// Do the full wave solve. Note: This currently prints out the 
// energy for every time step.
void Wave_Solve::Solve_PDE(Darray2& wm, Darray2& w, Darray2& wp, Darray2& lap, Darray2& u, double* w_ptr, MPI_Comm CART_COMM){
	double energy_old = 1e10;
	double energy;
    MPI_Request send_req[4];
    MPI_Request recv_req[4];

    Enforce_BC(w);

    // First step in Trap Rule
	for(int i=istart;i<=iend;i++){
        for (int j=jstart;j<=jend;j++){
            u(i,j) = 0.5*setup.dt*w(i,j)*0.75;
        }
    }

    // Compute Laplacian at t=0 and Taylor expand to get 
    // w(-dt,x,y)
	Communicate_Solution(CART_COMM,w_ptr,send_req,recv_req);
	Compute_Laplacian_NB(w, lap, recv_req);
	MPI_Waitall( 4, send_req, MPI_STATUSES_IGNORE );
	Taylor_Expand(wm,w,lap);


	for(int i =0;i<setup.nsteps;i++){
	    Time_Step(wm,w,wp,lap);
        for(int k=istart;k<=iend;k++){
            for(int l=jstart;l<=jend;l++){
                u(k,l) = u(k,l) + setup.dt*w(k,l)*(cos(setup.omega*(*this).t)-0.25);
            }
        }
        
		// energy = Compute_Energy(wm, w, lap,CART_COMM);
		// cout << scientific << right << setw(14)<< abs(energy-energy_old)/abs(energy_old) << "\n";
		// energy_old = energy;
	    Communicate_Solution(CART_COMM,w_ptr,send_req,recv_req);
	    Compute_Laplacian_NB(w, lap, recv_req);
        MPI_Waitall( 4, send_req, MPI_STATUSES_IGNORE );
	}

	// Correct trap rule at final integration point and scale by normalization constant
    for(int k=istart;k<=iend;k++){
        for(int l=jstart;l<=jend;l++){
            u(k,l) = u(k,l) - 0.5*setup.dt*w(k,l)*(cos(setup.omega*(*this).t)-0.25);
            u(k,l) = setup.omega*u(k,l)/M_PI;
        }
    }
}


// Perform one iteration of WaveHoltz.
void Wave_Solve::Evolve_and_Project(Darray2& b, Darray2& uold, Darray2& wm, Darray2& w, Darray2& wp, Darray2& lap,Darray2& u, double* w_ptr, MPI_Comm CART_COMM){
    MPI_Request send_req[4];
    MPI_Request recv_req[4];

    // Copy initial data for application of operator
    uold.copy(w);

	Solve_PDE(wm, w, wp, lap ,u, w_ptr, CART_COMM);
	(*this).t = 0;

	// Action of SPD WaveHoltz operator
	for(int i=istart;i<=iend;i++){
        for (int j=jstart;j<=jend;j++){
            u(i,j) = uold(i,j)-u(i,j)+b(i,j);
        }
    }
}
