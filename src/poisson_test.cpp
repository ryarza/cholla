#ifdef POISSON_TEST
#include "grid3D.h"
#include "io.h"
#include "math.h"

int Grid3D::poissonErrorNorm(){

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

  Real correctl2norm;
  if (nx_global == 64 ){
    correctl2norm = 0.00012863286755550373;
  }
  else if ( nx_global == 128 ){
    correctl2norm = 3.209325613406217e-5;
  }
  else if ( nx_global == 256 ){
    correctl2norm = 8.004489514884087e-6;
  }
  else if (nx_global == 512 ){
    correctl2norm = 1.9801531450853054e-6;
  }
  else{
    chprintf("Unsupported resolution!");
    exit(-1);
  }

  if ( fabs(l2norm / correctl2norm - 1.) < 1.e-10 ){
    return 1;
  }
  else{
    return 0;
  }


}

#endif
