#ifdef STARS

#include "../global.h"
#include "../grid3D.h"
//#include "../io.h"
#include <math.h>

#ifdef MPI_CHOLLA
#include "../mpi_routines.h"
#endif

//void Grid3D::IntegrateOrbit(Real )

//Boneless integration method
	

// Given the density in the cells and the position of the black hole, this function returns the three components of the acceleration of the black hole
void Grid3D::AccBh(Real posBhx, Real posBhy, Real posBhz, Real *accBhx, Real *accBhy, Real *accBhz){

//These variables hold the positions of the cell centers
	Real posx, posy, posz, r;
//Will hold density
	Real rho;

//We need the volume of the cell because we compute the acceleration between two point masses: the BH and a mass rho * dV at the center of the cell
	Real dV = H.dx * H.dy * H.dz;

//Temp
	Real tmpAccBhx = 0.;
	Real tmpAccBhy = 0.;
	Real tmpAccBhz = 0.;

// To compute the acceleration the BH experiences, we sum over the acceleration caused by every cell in the star.

	int id;
//The disgraceful N^3
	for (int k=H.n_ghost; k<H.nz-H.n_ghost; k++) {
		for (int j=H.n_ghost; j<H.ny-H.n_ghost; j++) {
			for (int i=H.n_ghost; i<H.nx-H.n_ghost; i++) {
        id = i + j*H.nx + k*H.nx*H.ny;
				Get_Position(i, j, k, &posx, &posy, &posz);

				r = sqrt( ( posx - posBhx ) * ( posx - posBhx ) + ( posy - posBhy ) * ( posy - posBhy ) + ( posz - posBhz ) * ( posz - posBhz ) );
				rho = C.density[id];

				tmpAccBhx += rho * ( posx - posBhx ) / pow(r, 3.);
				tmpAccBhy += rho * ( posy - posBhy ) / pow(r, 3.);
				tmpAccBhz += rho * ( posz - posBhz ) / pow(r, 3.);
			}
		}
	}

	tmpAccBhx *= dV * G_CGS;
	tmpAccBhy *= dV * G_CGS;
	tmpAccBhz *= dV * G_CGS;

	#ifdef MPI_CHOLLA
	*accBhx = ReduceRealSum(tmpAccBhx);
	*accBhy = ReduceRealSum(tmpAccBhy);
	*accBhz = ReduceRealSum(tmpAccBhz);
	#else
	*accBhx = tmpAccBhx;
	*accBhy = tmpAccBhy;
	*accBhz = tmpAccBhz;
	#endif

}

/*

void Grid3D::IntegrateOrbit(Real posBhx, Real posBhy, Real posBhz, Real velBhx, Real velBhy, Real velBhz, Real accBhx, Real accBhy, Real accBhz, Real posStx, Real posSty, Real posStz, Real velStx, Real velSty, Real velStz, Real accStx, Real accSty, Real accStz){

//New positions and velocities for the star and bh
	Real newPosBhx, newPosBhy, newPosBhz, newPosStx, newPosSty, newPosStz;
	Real newVelBhx, newVelBhy, newVelBhz, newVelStx, newVelSty, newVelStz;

//Boneless integration
	newPosBhx = posBhx + H.dt * velBhx + 0.5 * H.dt * H.dt * accBhx;
	newPosBhy = posBhy + H.dt * velBhy + 0.5 * H.dt * H.dt * accBhy;
	newPosBhz = posBhz + H.dt * velBhz + 0.5 * H.dt * H.dt * accBhz;

	newPosStx = posStx + H.dt * velStx + 0.5 * H.dt * H.dt * accStx;
	newPosSty = posSty + H.dt * velSty + 0.5 * H.dt * H.dt * accSty;
	newPosStz = posStz + H.dt * velStz + 0.5 * H.dt * H.dt * accStz;
}

*/

#endif
