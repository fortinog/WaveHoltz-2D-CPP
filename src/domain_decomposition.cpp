#include "domain_decomposition.h"
#include <iostream>


// FMG double check this
void decompose1d( int n, int m, int i, int* s, int* e )
{
    const int length  = n / m;
    const int deficit = n % m;
    *s =  i * length + ( i < deficit ? i : deficit );
    *e = *s + length - ( i < deficit ? 0 : 1 );
    if ( ( *e >= n ) || ( i == m - 1 ) ) *e = n - 1;
}

// Set up Cartesian grid decomposition of processes.
MPI_Comm mpi_decompose_domain(int num_proc, int NX, int NY, int rerank, int* rank, int* dims,
							  Subdomain* orig_grid){
	MPI_Comm cart_comm;
	int coords[2];
	int periodic[2] = {0,0}; // assume domain is not periodic

	MPI_Dims_create(num_proc, 2, dims);
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periodic, rerank, &cart_comm);
	MPI_Cart_get(cart_comm, 2, dims, periodic, coords);
	MPI_Comm_rank(cart_comm, rank);

	// Determine all of my neighbor processes.
    MPI_Cart_shift( cart_comm, 0, 1, &(orig_grid->left_neigh),
                    &(orig_grid->right_neigh) );
    MPI_Cart_shift( cart_comm, 1, 1, &(orig_grid->down_neigh),
                    &(orig_grid->up_neigh) );

    // Figure out the extents and size of my portion of the grid.
    decompose1d( NX, dims[0], coords[0], &(orig_grid->x_s), &(orig_grid->x_e) );
    decompose1d( NY, dims[1], coords[1], &(orig_grid->y_s), &(orig_grid->y_e) );
    orig_grid->nx_loc = orig_grid->x_e - orig_grid->x_s + 1;
    orig_grid->ny_loc = orig_grid->y_e - orig_grid->y_s + 1;


    // Get my local coordinates
    orig_grid->px = coords[0];
    orig_grid->py = coords[1];

    // Adjust domain parameters to account for inter-domain halo
    // boundary data.  If we have a neighbor in a given direction
    // (rank of neighbor is non-negative) then we need to adjust the
    // starting or ending index.

    // ghost_grid->left_neigh  = orig_grid->left_neigh;
    // ghost_grid->right_neigh = orig_grid->right_neigh;
    // ghost_grid->up_neigh = orig_grid->up_neigh;
    // ghost_grid->down_neigh = orig_grid->down_neigh;

    // ghost_grid->x_s = orig_grid->x_s - (ghost_grid->left_neigh  >= 0 ? 1 : 0 );
    // ghost_grid->x_e = orig_grid->x_e + (ghost_grid->right_neigh >= 0 ? 1 : 0 );
    // ghost_grid->y_s = orig_grid->y_s - (ghost_grid->down_neigh >= 0 ? 1 : 0 );
    // ghost_grid->y_e = orig_grid->y_e + (ghost_grid->up_neigh >= 0 ? 1 : 0 );
    // ghost_grid->nx_loc = ghost_grid->x_e - ghost_grid->x_s + 1;
    // ghost_grid->ny_loc = ghost_grid->y_e - ghost_grid->y_s + 1;

    // ghost_grid->nolp = 1;
    orig_grid->nolp = 1;

	return cart_comm;
}