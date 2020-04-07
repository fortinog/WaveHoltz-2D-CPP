#include "wave_solve.hpp"

// Set up the grid
Wave_Solve::Wave_Solve(){
	int N, M, nsteps;
	double err;
	double dt, dt2, idt2; 
	t = 0.0;

	N 	   = setup.N;
	M 	   = setup.M;
	nsteps = setup.nsteps;
	dt     = setup.dt;
	dt2    = pow(setup.dt,2);
	idt2   = 1.0/dt2;

	x.define(1,0,N+1);
	y.define(1,0,M+1);

	// Setup grid: Note i=0,N+1, j=0,M+1 correspond to the boundary
	// if we consider solving in a box.
	for(int i=0;i<=N+1;i++) x(i) = setup.x_L + i*setup.hx;
	for(int i=0;i<=M+1;i++) y(i) = setup.y_L + i*setup.hy;

	w.define(1,0,N+1,0,M+1);
	w.set_value(0.0);
	wp.define(1,0,N+1,0,M+1);
	wp.set_value(0.0);
	wm.define(1,0,N+1,0,M+1);
	wm.set_value(0.0);
	lap.define(1,0,N+1,0,M+1);
	lap.set_value(0.0);
	mask.define(1,0,N+1,0,M+1);
	mask.set_value(0);

	// Fill in w, wm, lap so that we can start time stepping immediately
	Set_Initial_Data();
}

// Compute the mask for all points in the domain via 
// the level set function
void Wave_Solve::Compute_Mask(){

	double level_set_val;
	double x_val;
	int ghost_ctr = 0;
	int num_int = 0;
	int N = setup.N;
	int M = setup.M;
	bool bdry_check[4];

	for(int i = 0;i<=M;i++){
		x_val = x(i);
		for(int j = 0;j<=N;j++){
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
				mask(i,j) = -1;
				ghost_ctr++;
			}
		}
	}

	// Identify an interior ghost point candidate
	for(int i = 0;i<=M;i++){
		x_val = x(i);
		for(int j = 0;j<=N;j++){
			if(mask(i,j) == 1){
				bdry_check[0] = mask(i-1,j) == 0;
				bdry_check[1] = mask(i+1,j) == 0;
				bdry_check[2] = mask(i,j-1) == 0;
				bdry_check[3] = mask(i,j+1) == 0;
			}

			if(bdry_check[0] || bdry_check[1] || bdry_check[2] || bdry_check[3]){
				mask(i,j) = -1;
				ghost_ctr++;
				num_int -= 1;
			}
		}
	}
}

// Fill in with initial data, compute the Laplacian, and do a 
// Taylor expansion so that we may start time stepping after this.
void Wave_Solve::Set_Initial_Data(){
	Twilight mms(setup);
	double x0;

	for(int i=0;i<=setup.N+1;i++){
		x0 = x(i);
		for(int j=0;j<=setup.M+1;j++){
			w(i,j) = mms.trigTwilight(0,x0,0,y(j),0,0.0);
		}
	}
	Compute_Laplacian();
	Taylor_Expand();
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
	i = 0;
	x0 = x(i);
	for(j=1;j<=setup.M;j++){
		v(i,j) = mms.trigTwilight(0,x0,0,y(j),0,t);
	}

	// Right boundary
	i = setup.N+1;
	x0 = x(i);
	for(j=1;j<=setup.M;j++){
		v(i,j) = mms.trigTwilight(0,x0,0,y(j),0,t);
	}

	// Bottom boundary
	j = 0;
	y0 = y(j);
	for(i=1;i<=setup.N;i++){
		v(i,j) = mms.trigTwilight(0,x(i),0,y0,0,t);
	}

	// Right boundary
	j = setup.M+1;
	y0 = y(j);
	for(i=1;i<=setup.N;i++){
		v(i,j) = mms.trigTwilight(0,x(i),0,y0,0,t);
	}
}

void Wave_Solve::Compute_Laplacian(){
	int N = setup.N;
	int M = setup.M;

	double ihx2 = 1.0/pow(setup.hx,2);
	double ihy2 = 1.0/pow(setup.hy,2);
	double c1,c2,c3;
	double x0,xp,xm;
	double y0,yp,ym;


	// Compute the usual Laplacian (interior only)
	for(int i=1;i<=N;i++){
		xm = x(i-1);
		x0 = x(i);
		xp = x(i+1);
		for(int j=1;j<=M;j++){
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
double Wave_Solve::Compute_Laplacian_Error(const int flag){
	double wxx,true_lap;
	double err = 0.0;
	double norm_sol= 0.0;
	double x0;
	double max_val = 0.0;
	double max_sol = 0.0;
	double t = (*this).t;
	Twilight mms(setup);

	switch(flag){
		case 1:
			for(int i=1;i<=setup.N;i++){
				x0 = x(i);
				for(int j=1;j<=setup.M;j++){
					wxx = mms.trigTwilight(2,x0,0,y(j),0,t);
					true_lap = wxx + mms.trigTwilight(0,x0,2,y(j),0,t);
					err += abs(lap(i,j)-true_lap);
					norm_sol += abs(true_lap);
				}
			}
			err = err/norm_sol;
			break;
		case 2:
			for(int i=1;i<=setup.N;i++){
				x0 = x(i);
				for(int j=1;j<=setup.M;j++){
					wxx = mms.trigTwilight(2,x0,0,y(j),0,t);
					true_lap = wxx + mms.trigTwilight(0,x0,2,y(j),0,t);
					err += pow(lap(i,j)-true_lap,2);
					norm_sol += pow(true_lap,2);
				}
			}
			err = sqrt(err/norm_sol);
			break;
		case 3:
			for(int i=1;i<=setup.N;i++){
				x0 = x(i);
				for(int j=1;j<=setup.M;j++){
					wxx = mms.trigTwilight(2,x0,0,y(j),0,t);
					true_lap = wxx + mms.trigTwilight(0,x0,2,y(j),0,t);
					if(abs(lap(i,j)-true_lap) > max_val) max_val = abs(lap(i,j)-true_lap);
					if(abs(true_lap) > max_sol) max_sol = abs(true_lap);
				}
			}
			err = max_val/max_sol;
			break;
	}

	return err;
}

// Only in the constant coefficient case, flag=1 L1 norm, flag=2 L2 norm, flag=3 L_inf norm
double Wave_Solve::Compute_Solution_Error(const int flag){
	double true_sol;
	double err = 0.0;
	double norm_sol= 0.0;
	double x0;
	double max_val = 0.0;
	double max_sol = 0.0;
	double t = (*this).t;
	Twilight mms(setup);

	switch(flag){
		case 1:
			for(int i=1;i<=setup.N;i++){
				x0 = x(i);
				for(int j=1;j<=setup.M;j++){
					true_sol = mms.trigTwilight(0,x0,0,y(j),0,t);
					err += abs(w(i,j)-true_sol);
					norm_sol += abs(true_sol);
				}
			}
			err = err/norm_sol;
			break;
		case 2:
			for(int i=1;i<=setup.N;i++){
				x0 = x(i);
				for(int j=1;j<=setup.M;j++){
					true_sol = mms.trigTwilight(0,x0,0,y(j),0,t);
					err += pow(w(i,j)-true_sol,2);
					norm_sol += pow(true_sol,2);
				}
			}
			err = sqrt(err/norm_sol);
			break;
		case 3:
			for(int i=1;i<=setup.N;i++){
				x0 = x(i);
				for(int j=1;j<=setup.M;j++){
					true_sol = mms.trigTwilight(0,x0,0,y(j),0,t);
					if(abs(w(i,j)-true_sol) > max_val) max_val = abs(w(i,j)-true_sol);
					if(abs(true_sol) > max_sol) max_sol = abs(true_sol);
				}
			}
			err = max_val/max_sol;
			break;
	}

	return err;
}

// Double the number of grid points (for the purposes of testing convergence)
void Wave_Solve::Refine_Grid(){
	int N = 2*setup.N;
	int M = 2*setup.M;

	// Double grid points and recalculate grid size
	setup.N = N;
	setup.M = M;
	setup.hx = (setup.x_R-setup.x_L)/double(N+1);
	setup.hy = (setup.y_R-setup.y_L)/double(M+1);

	// Setup grid
	x.define(1,0,N+1);
	y.define(1,0,M+1);
	for(int i=0;i<=N+1;i++) x(i) = setup.x_L + i*setup.hx;
	for(int i=0;i<=M+1;i++) y(i) = setup.y_L + i*setup.hy;

	// Resize working arrays
	w.define(1,0,N+1,0,M+1);
	w.set_value(0.0);
	lap.define(1,0,N+1,0,M+1);
	lap.set_value(0.0);
	mask.define(1,0,N+1,0,M+1);
	mask.set_value(0);	
}

//	Do a Taylor expansion to fill in w_{-1} \approx w(-dt).
void Wave_Solve::Taylor_Expand(){
	Twilight mms(setup);
	double dt = setup.dt;
	double dt2 = setup.dt2;
	double x0;
	for(int i=1;i<=setup.N;i++){
		x0 = x(i);
		for(int j=1;j<=setup.M;j++){
			wm(i,j) = w(i,j) - dt*mms.trigTwilight(0,x0,0,y(j),1,0.0) + 0.5*dt2*(lap(i,j) + forcing(x0,y(j),0.0));
		}
	}
	Enforce_BC(wm);
}

// Perform a single time step forward, centered difference in time
void Wave_Solve::Time_Step(){
	int N = setup.N;
	int M = setup.M;
	double dt2 = setup.dt2;
	double x0;
	double t = (*this).t;

	for(int i=1;i<=N;i++){
		x0 = x(i);
		for(int j=1;j<=M;j++){
			wp(i,j) = 2.0*w(i,j) - wm(i,j) + dt2*(lap(i,j) + forcing(x0,y(j),t)); 
		}
	}
	Enforce_BC(wp);

	wm.copy(w);
	w.copy(wp);
	(*this).t += setup.dt;
}

// Compute the energy at a given time step. (Calculation from Patrick Joly's paper)
double Wave_Solve::Compute_Energy(){
	int N = setup.N;
	int M = setup.M;
	double dt2 = setup.dt2;
	double energy = 0.0;
	double idt2 = setup.idt2;

	for(int i=1;i<=N;i++){
		for(int j=1;j<=M;j++){
			energy = energy + idt2*pow(wm(i,j) - w(i,j),2) - w(i,j)*lap(i,j);
		}
	}

	return 0.5*setup.hx*setup.hy*energy;
}

// Set everything to zero for convergence tests
void Wave_Solve::Clear_Data(){
	(*this).t = 0.0;
	w.set_value(0.0);
	wp.set_value(0.0);
	wm.set_value(0.0);
	lap.set_value(0.0);
}

// Do the full wave solve. Note: This currently prints out the 
// energy for every time step.
void Wave_Solve::Solve_PDE(){
	double energy_old = 1e10;
	double energy;

	for(int i =0;i<setup.nsteps;i++){
	    Time_Step();
	    energy = Compute_Energy();
        cout << scientific << right << setw(14)<< abs(energy-energy_old)/abs(energy_old) << "\n";
	    energy_old = energy;
	    Compute_Laplacian();
	}

}