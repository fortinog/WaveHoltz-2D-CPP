#include "wave_solve.hpp"

// Set up the grid
Wave_Solve::Wave_Solve(Subdomain Local_Grid, MPI_Comm CART_COMM){
	int nsteps, nolp;
	int ndims = 2;
	int x_off, y_off;
	int rank;
	double err;
	double dt, dt2, idt2; 
	Darray2 w, lap, wm, wp;

	t = 0.0;
    MPI_Comm_rank( CART_COMM, &rank );

	N 	   = Local_Grid.nx_loc;
	M 	   = Local_Grid.ny_loc;
	nolp   = Local_Grid.nolp;
	// istart = 1-nolp;
	// iend   = N+nolp;
	// jstart = 1-nolp;
	// jend   = M+nolp;

	up_neigh    = Local_Grid.up_neigh;
	down_neigh  = Local_Grid.down_neigh;
	left_neigh  = Local_Grid.left_neigh;
	right_neigh = Local_Grid.right_neigh;


	// istart = (left_neigh == MPI_PROC_NULL) ? istart+1 : istart;
	// iend   = (right_neigh == MPI_PROC_NULL) ? iend-1 : iend;
	// istart = (left_neigh == MPI_PROC_NULL) ? istart+1 : istart;
	// istart = (left_neigh == MPI_PROC_NULL) ? istart+1 : istart;

	N = (left_neigh  == MPI_PROC_NULL) ? N-1 : N;
	N = (right_neigh == MPI_PROC_NULL) ? N-1 : N;
	M = (down_neigh  == MPI_PROC_NULL) ? M-1 : M;
	M = (up_neigh    == MPI_PROC_NULL) ? M-1 : M;

	istart = 1-nolp;
	iend   = N+nolp;
	jstart = 1-nolp;
	jend   = M+nolp;

	// boundary loop indices
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

	nsteps = setup.nsteps;
	dt     = setup.dt;
	dt2    = pow(setup.dt,2);
	idt2   = 1.0/dt2;

	x.define(1,istart,iend);
	y.define(1,jstart,jend);


	// cout << "istart " << istart << " iend : " << iend << "\n";
	// cout << "jstart " << jstart << " jend : " << jend << "\n";

	// Setup local grid: Note i=0,N+1, j=0,M+1 correspond to either the boundary
	// or a ghost region for the given subdomain.

	// FMG: if I have a ghost, adjust either end
	for(int i=istart;i<=iend;i++) x(i) = setup.x_L + (i + Local_Grid.x_s - x_off)*setup.hx;
	for(int i=jstart;i<=jend;i++) y(i) = setup.y_L + (i + Local_Grid.y_s - y_off)*setup.hy;

	// for(int i =istart;i<=iend; i++) cout << x(i) << " ";
	// cout << "\n";
	// for(int i=jstart;i<=jend;i++) cout << y(i) << " ";
	// cout << "\n";

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
	mask.set_value(0);

    double* w_ptr = w.c_ptr();

    size_lr = iend-istart+1;
    size_ud = jend-jstart+1;

    // Set up all communication subarrays.
    // Setup_Subarrays(nolp, M, N);
    Setup_Subarrays(nolp);

	// Fill in w, wm, lap so that we can start time stepping immediately
	Set_Initial_Data(wm, w, lap);
	// Communicate_Solution(CART_COMM, w_ptr);
	// for(int i=0;i<3;i++){
	//     Time_Step(wm, w, wp, lap);
	//     Compute_Laplacian(w,lap);
	// 	Communicate_Solution(CART_COMM, w_ptr);
	// }
    Solve_PDE(wm, w, wp, lap, w_ptr, CART_COMM);

	cout << scientific << right << setw(14)<< Compute_Solution_Error(1, w, CART_COMM);
	cout << scientific << right << setw(14)<< Compute_Solution_Error(2, w, CART_COMM);
	cout << scientific << right << setw(14)<< Compute_Solution_Error(3, w, CART_COMM);
	cout << "\n";  
	MPI_Barrier(CART_COMM);

	cout << scientific << right << setw(14)<< Compute_Laplacian_Error(1, lap, CART_COMM);
	cout << scientific << right << setw(14)<< Compute_Laplacian_Error(2, lap, CART_COMM);
	cout << scientific << right << setw(14)<< Compute_Laplacian_Error(3, lap, CART_COMM);
	cout << "\n";  
	MPI_Barrier(CART_COMM);

	Finalize();
}

// Compute the mask for all points in the domain via 
// the level set function (FMG fix)
void Wave_Solve::Compute_Mask(){

	// double level_set_val;
	// double x_val;
	// int ghost_ctr = 0;
	// int num_int = 0;
	// int N = setup.N;
	// int M = setup.M;
	// bool bdry_check[4];

	// for(int i = 0;i<=M;i++){
	// 	x_val = x(i);
	// 	for(int j = 0;j<=N;j++){
	// 		level_set_val = setup.Level_Set(x_val,y(j));


	// 		// Interior point
	// 		if(level_set_val < 0.0){
	// 			mask(i,j) = 1;
	// 			num_int++;
	// 		}
	// 		// Exterior point
	// 		else if(level_set_val > 0.0){
	// 			mask(i,j) = 0;
	// 		} 
	// 		// Boundary point
	// 		else{
	// 			mask(i,j) = -1;
	// 			ghost_ctr++;
	// 		}
	// 	}
	// }

	// // Identify an interior ghost point candidate
	// for(int i = 0;i<=M;i++){
	// 	x_val = x(i);
	// 	for(int j = 0;j<=N;j++){
	// 		if(mask(i,j) == 1){
	// 			bdry_check[0] = mask(i-1,j) == 0;
	// 			bdry_check[1] = mask(i+1,j) == 0;
	// 			bdry_check[2] = mask(i,j-1) == 0;
	// 			bdry_check[3] = mask(i,j+1) == 0;
	// 		}

	// 		if(bdry_check[0] || bdry_check[1] || bdry_check[2] || bdry_check[3]){
	// 			mask(i,j) = -1;
	// 			ghost_ctr++;
	// 			num_int -= 1;
	// 		}
	// 	}
	// }
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
	Enforce_BC(wm);
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
// the boundary data. FMG:: Update to only update if I have a physical
// boundary
void Wave_Solve::Enforce_BC(Darray2 &v){
	Twilight mms(setup);
	int i,j;
	double x0, y0;
	double t = (*this).t;

	// // Left boundary
	// i = istart;
	// x0 = x(i);
	// for(j=jstart+1;j<=jend-1;j++){
	// 	v(i,j) = mms.trigTwilight(0,x0,0,y(j),0,t);
	// }

	// // Right boundary
	// i = iend;
	// x0 = x(i);
	// for(j=jstart+1;j<=jend-1;j++){
	// 	v(i,j) = mms.trigTwilight(0,x0,0,y(j),0,t);
	// }

	// // Bottom boundary
	// j = jstart;
	// y0 = y(j);
	// for(i=istart+1;i<=iend-1;i++){
	// 	v(i,j) = mms.trigTwilight(0,x(i),0,y0,0,t);
	// }

	// // Right boundary
	// j = jend;
	// y0 = y(j);
	// for(i=istart+1;i<=iend-1;i++){
	// 	v(i,j) = mms.trigTwilight(0,x(i),0,y0,0,t);
	// }
	// Left boundary
	i = istart;
	x0 = x(i);
	// cout << "x0 :" << x0 << "\n";
	for(j=j_left_start;j<=j_left_end;j++){
		v(i,j) = mms.trigTwilight(0,x0,0,y(j),0,t);
	}

	// Right boundary
	i = iend;
	x0 = x(i);
	// cout << "x0 :" << x0 << "\n";
	for(j=j_right_start;j<=j_right_end;j++){
		v(i,j) = mms.trigTwilight(0,x0,0,y(j),0,t);
	}

	// Bottom boundary
	j = jstart;
	y0 = y(j);
	// cout << "y0 :" << y0 << "\n";
	for(i=i_down_start;i<=i_down_end;i++){
		v(i,j) = mms.trigTwilight(0,x(i),0,y0,0,t);
	}

	// Top boundary
	j = jend;
	y0 = y(j);
	// cout << "y0 :" << y0 << "\n";
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

	// Compute the usual Laplacian (interior only)
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

// // Double the number of grid points (for the purposes of testing convergence)
// void Wave_Solve::Refine_Grid(){
// 	int N = 2*setup.N;
// 	int M = 2*setup.M;

// 	// Double grid points and recalculate grid size
// 	setup.N = N;
// 	setup.M = M;
// 	setup.hx = (setup.x_R-setup.x_L)/double(N+1);
// 	setup.hy = (setup.y_R-setup.y_L)/double(M+1);

// 	// Setup grid
// 	x.define(1,0,N+1);
// 	y.define(1,0,M+1);
// 	for(int i=0;i<=N+1;i++) x(i) = setup.x_L + i*setup.hx;
// 	for(int i=0;i<=M+1;i++) y(i) = setup.y_L + i*setup.hy;

// 	// Resize working arrays
// 	w.define(1,0,N+1,0,M+1);
// 	w.set_value(0.0);
// 	lap.define(1,0,N+1,0,M+1);
// 	lap.set_value(0.0);
// 	mask.define(1,0,N+1,0,M+1);
// 	mask.set_value(0);	
// }

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

// // Set everything to zero for convergence tests
// void Wave_Solve::Clear_Data(){
// 	(*this).t = 0.0;
// 	w.set_value(0.0);
// 	wp.set_value(0.0);
// 	wm.set_value(0.0);
// 	lap.set_value(0.0);
// }

// Do the full wave solve. Note: This currently prints out the 
// energy for every time step.
void Wave_Solve::Solve_PDE(Darray2& wm, Darray2& w, Darray2& wp, Darray2& lap, double* w_ptr, MPI_Comm CART_COMM){
	double energy_old = 1e10;
	double energy;

	for(int i =0;i<setup.nsteps;i++){
	    Time_Step(wm,w,wp,lap);
	    energy = Compute_Energy(wm, w, lap,CART_COMM);
        cout << scientific << right << setw(14)<< abs(energy-energy_old)/abs(energy_old) << "\n";
	    energy_old = energy;
	    Communicate_Solution(CART_COMM,w_ptr);
	    Compute_Laplacian(w, lap);
	}

}

// This routine sets up the communication subarrays.
void Wave_Solve::Setup_Subarrays(const int nolp){
	int ndims = 2;
	int M = size_ud - 2*nolp;
	int N = size_lr - 2*nolp;

	/************************************************************/
	/*	    		COMMUNICATORS IN X DIRECTION				*/
	/************************************************************/

	// Left send
    int array_of_subsizes[2] = {nolp,M};
    // int array_of_sizes[2] = {N+2*nolp,M+2*nolp};
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
    // array_of_starts[0] = N+nolp;
    // array_of_starts[1] = nolp;
    array_of_starts[0] = size_lr-1;
    array_of_starts[1] = nolp;
    MPI_Type_create_subarray(ndims,
                             array_of_sizes,
                             array_of_subsizes,
                             array_of_starts,
                             order, MPI_DOUBLE, &sub_recv_right);
    MPI_Type_commit(&sub_recv_right);

    // Right send
    // array_of_starts[0] = N;
    // array_of_starts[1] = nolp;
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
	// array_of_starts[0] = nolp;
	// array_of_starts[1] = M;
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
	// array_of_starts[0] = nolp;
	// array_of_starts[1] = nolp;
	array_of_starts[0] = nolp;
	array_of_starts[1] = nolp;
    MPI_Type_create_subarray(ndims,
                             array_of_sizes,
                             array_of_subsizes,
                             array_of_starts,
                             order, MPI_DOUBLE, &sub_send_down);
    MPI_Type_commit(&sub_send_down);

    // Up receive
    // array_of_starts[0] = nolp;
    // array_of_starts[1] = M+nolp;
    array_of_starts[0] = nolp;
    array_of_starts[1] = size_ud - 1;
    MPI_Type_create_subarray(ndims,
                             array_of_sizes,
                             array_of_subsizes,
                             array_of_starts,
                             order, MPI_DOUBLE, &sub_recv_up);
    MPI_Type_commit(&sub_recv_up);

}

// Communicate all ghost points. TODO: Non-blocking communication
// to hide latency.
void Wave_Solve::Communicate_Solution(MPI_Comm CART_COMM, double* w_ptr){
	int tag = 22;
	MPI_Status status;

	// Send left, Receive right    
    MPI_Send(w_ptr, 1, sub_send_left, left_neigh, tag, CART_COMM);
    MPI_Recv(w_ptr, 1, sub_recv_right, right_neigh, tag, CART_COMM, &status);

	// Send right, Receive left    
    tag = 23;
    MPI_Send(w_ptr, 1, sub_send_right, right_neigh, tag, CART_COMM);
    MPI_Recv(w_ptr, 1, sub_recv_left, left_neigh, tag, CART_COMM, &status);

	// Send up, Receive down 
	tag = 24;   
    MPI_Send(w_ptr, 1, sub_send_up, up_neigh, tag, CART_COMM);
    MPI_Recv(w_ptr, 1, sub_recv_down, down_neigh, tag, CART_COMM, &status);

	// Send down, Receive up    
    tag = 25;
    MPI_Send(w_ptr, 1, sub_send_down, down_neigh, tag, CART_COMM);
    MPI_Recv(w_ptr, 1, sub_recv_up, up_neigh, tag, CART_COMM, &status);
    MPI_Barrier(CART_COMM);
}

void Wave_Solve::Finalize(){
	MPI_Type_free( &sub_send_left ); 
	MPI_Type_free( &sub_recv_left);
	MPI_Type_free( &sub_send_right ); 
	MPI_Type_free( &sub_recv_right);
	MPI_Type_free( &sub_send_up );
	MPI_Type_free( &sub_recv_up);
	MPI_Type_free( &sub_send_down ); 
	MPI_Type_free( &sub_recv_down);
}