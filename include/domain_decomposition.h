#ifndef DOMAIN_DECOMPOSITION
#define DOMAIN_DECOMPOSITION
#include <mpi.h>

typedef struct Subdomain{
	int nolp;
	int x_s;
	int x_e;
	int y_s;
	int y_e;
	int nx_loc;
	int ny_loc;
	int up_neigh;
	int down_neigh;
	int left_neigh;
	int right_neigh;
	int px;
	int py;
} Subdomain;

void decompose1d( int n, int m, int i, int* s, int* e );

MPI_Comm mpi_decompose_domain(int num_proc, int NX, int NY, int rerank, int* rank, int* dims,
							  Subdomain* orig_grid);

#endif