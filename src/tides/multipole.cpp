#include "../global.h"
#include "../grid3D.h"
#include "../io.h"

#ifdef MPI_CHOLLA
#include "../mpi_routines.h"
#endif

void Grid3D::getMoments(){

	int id;
	Real r, phi, theta, rhosq, totrhosq, totrhosqtemp;
	Real pos[3], centertemp[3];
	std::complex<Real> Qtemp;

// Find the center of the multipole expansion according to Couch et al. 2013
	totrhosqtemp = 0.;
	for ( int ii = 0; ii < 3; ii++ ) centertemp[ii] = 0.;

	for (int k = H.n_ghost; k < H.nz - H.n_ghost; k++) {
		for (int j = H.n_ghost; j < H.ny - H.n_ghost; j++) {
			for (int i = H.n_ghost; i < H.nx - H.n_ghost; i++) {

				id = i + j*H.nx + k*H.nx*H.ny;
				Get_Position(i, j, k, &pos[0], &pos[1], &pos[2]);
				rhosq = C.density[id] * C.density[id];

				for ( int ii = 0; ii < 3; ii++ ) centertemp[ii] += rhosq * pos[ii];
				totrhosqtemp += rhosq;
				
			}
		}
	}

	#ifdef MPI_CHOLLA
	totrhosq = ReduceRealSum(totrhosqtemp);
	for ( int ii = 0; ii < 3; ii++ ) Grav.center[ii] = ReduceRealSum(centertemp[ii]);
	#else
	totrhosq = totrhosqtemp;
	for ( int ii = 0; ii < 3; ii++ ) Grav.center[ii] = centertemp[ii];
	#endif//MPI_CHOLLA

	for ( int ii = 0; ii < 3; ii++ ) Grav.center[ii] /= totrhosq;

/*
	Grav.center[0] = 0.;
	Grav.center[1] = 0.;
	Grav.center[2] = 0.;
*/
	chprintf("Center: %.5e, %.5e, %.5e\n", Grav.center[0], Grav.center[1], Grav.center[2]);

//	Get the multipole moments
	for ( int l = 0; l < Grav.lmaxBoundaries + 1; l++ ){
		for ( int m = - l; m < l + 1; m++ ){

			Qtemp = 0.;

			for (int k = H.n_ghost; k < H.nz - H.n_ghost; k++) {
				for (int j = H.n_ghost; j < H.ny - H.n_ghost; j++) {
					for (int i = H.n_ghost; i < H.nx - H.n_ghost; i++) {

//					Get the x, y, and z coordinates of this point with respect to the center of the expansion
						id = i + j*H.nx + k*H.nx*H.ny;
						Get_Position(i, j, k, &pos[0], &pos[1], &pos[2]);
						for ( int ii = 0; ii < 3; ii++) pos[ii] -= Grav.center[ii];

//					Turn into spherical coordinates
						r = sqrt( pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2] );
						phi = atan2(pos[1], pos[0]);
						theta = acos( pos[2] / r );

//					Integrate
						Qtemp += pow(r, l) * conj(Grav.Y(l, m, theta, phi)) * C.density[id];
					}
				}
			}

			#ifdef MPI_CHOLLA
				Grav.Q[l][m + l] = ReduceComplexSum(Qtemp);
			#else
				Grav.Q[l][m + l] = Qtemp;
			#endif//MPI_CHOLLA

			Grav.Q[l][m + l] *= H.dx * H.dy * H.dz;

			chprintf("Q[%i][%i]=%.5e+(%.5e)i\n", l, m, real(Grav.Q[l][m+l]), imag(Grav.Q[l][m+l]));

		}
	}

	Real testphi = 0.2532;
	Real testtheta = 0.9758;
	int testm = -1;
	int testl = 2;
	chprintf("Y(%i, %i, %.5e, %.5e)=%.5e+(%.5e)i\n", testl, testm, testtheta, testphi, real(Grav.Y(testl,testm,testtheta, testphi)), imag(Grav.Y(testl, testm, testtheta, testphi)));

}
