#include "wave_solve.hpp"

// Set up the grid
Wave_Solve::Wave_Solve(){
	int N, M;
	double err;
	N = setup.N;
	M = setup.M;

	x.define(1,-1,N);
	y.define(1,-1,M);

	// Setup grid
	for(int i=-1;i<=N;i++) x(i) = setup.x_L + i*setup.hx;
	for(int i=-1;i<=M;i++) y(i) = setup.y_L + i*setup.hy;

	w.define(1,-1,N,-1,M);
	w.set_value(0.0);
	lap.define(1,-1,N,-1,M);
	lap.set_value(0.0);
	mask.define(1,-1,N,-1,M);
	mask.set_value(0);

	// Set_Initial_Data();
	// Compute_Laplacian();
	// err = Compute_Laplacian_Error(3);
	// cout << "The error is :" << err << "\n";
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

void Wave_Solve::Set_Initial_Data(){
	Twilight mms(setup);

	for(int i=-1;i<=setup.N;i++){
		for(int j=-1;j<=setup.M;j++){
			w(i,j) = mms.trigTwilight(0,x(i),0,y(j),0,0.5);
		}
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


	// Compute the usual Laplacian
	// lap.set_value(0.0);
	for(int i=0;i<N;i++){
		xm = x(i-1);
		x0 = x(i);
		xp = x(i+1);
		for(int j=0;j<M;j++){
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
	Twilight mms(setup);

	switch(flag){
		case 1:
			for(int i=0;i<setup.N;i++){
				x0 = x(i);
				for(int j=0;j<setup.M;j++){
					wxx = mms.trigTwilight(2,x0,0,y(j),0,0.5);
					true_lap = wxx + mms.trigTwilight(0,x0,2,y(j),0,0.5);
					err += abs(lap(i,j)-true_lap);
					norm_sol += abs(true_lap);
				}
			}
			err = err/norm_sol;
			break;
		case 2:
			for(int i=0;i<setup.N;i++){
				x0 = x(i);
				for(int j=0;j<setup.M;j++){
					wxx = mms.trigTwilight(2,x0,0,y(j),0,0.5);
					true_lap = wxx + mms.trigTwilight(0,x0,2,y(j),0,0.5);
					err += pow(lap(i,j)-true_lap,2);
					norm_sol += pow(true_lap,2);
				}
			}
			err = sqrt(err/norm_sol);
			break;
		case 3:
			for(int i=0;i<setup.N;i++){
				x0 = x(i);
				for(int j=0;j<setup.M;j++){
					wxx = mms.trigTwilight(2,x0,0,y(j),0,0.5);
					true_lap = wxx + mms.trigTwilight(0,x0,2,y(j),0,0.5);
					if(abs(lap(i,j)-true_lap) > max_val) max_val = abs(lap(i,j)-true_lap);
					if(abs(true_lap) > max_sol) max_sol = abs(true_lap);
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
	setup.hx = (setup.x_R-setup.x_L)/double(N-1);
	setup.hy = (setup.y_R-setup.y_L)/double(M-1);

	// Setup grid
	x.define(1,-1,N);
	y.define(1,-1,M);
	for(int i=-1;i<=N;i++) x(i) = setup.x_L + i*setup.hx;
	for(int i=-1;i<=M;i++) y(i) = setup.y_L + i*setup.hy;

	// Resize working arrays
	w.define(1,-1,N,-1,M);
	w.set_value(0.0);
	lap.define(1,-1,N,-1,M);
	lap.set_value(0.0);
	mask.define(1,-1,N,-1,M);
	mask.set_value(0);	

}