#ifdef POISSON_TEST
#include "grid3D.h"
#include "io.h"
#include "math.h"

void Grid3D::poissonErrorNorm(){

	Real l2norm;
	Real deltasq = 0.;

	#ifndef MPI_CHOLLA
	int nx_global = H.nx_real;
	int ny_global = H.ny_real;
	int nz_global = H.nz_real;
	#endif

	int apotidx, potidx;
	for ( int k = 0; k < H.nz_real; k++ ){
		for ( int j = 0; j < H.ny_real; j++ ){
			for ( int i = 0; i < H.nx_real; i++ ){
				apotidx = ( i + H.n_ghost ) + ( j + H.n_ghost ) * H.nx + ( k + H.n_ghost ) * H.nx * H.ny;
				potidx = (i+N_GHOST_POTENTIAL) + (j+N_GHOST_POTENTIAL)*(Grav.nx_local+2*N_GHOST_POTENTIAL) + (k+N_GHOST_POTENTIAL)*(Grav.nx_local+2*N_GHOST_POTENTIAL)*(Grav.ny_local+2*N_GHOST_POTENTIAL);

//				printf("apot pot %.10e %.10e\n", Grav.F.potential_h[potidx], C.analyticalPotential[apotidx]);
				if ( fabs(Grav.F.potential_h[potidx]) < 1.e-10 || fabs(C.analyticalPotential[apotidx]) < 1.e-10 ) chprintf("Potential is wrong...\n");

				deltasq += pow( C.analyticalPotential[apotidx] - Grav.F.potential_h[potidx], 2. );
			}
		}
	}


	#ifdef MPI_CHOLLA
	MPI_Allreduce(MPI_IN_PLACE, &deltasq, 1, MPI_CHREAL, MPI_SUM, world);
	#endif

	l2norm = sqrt( deltasq / nx_global / ny_global / nz_global );
	chprintf("L2 norm = %.20e\n", l2norm);

}

#endif
