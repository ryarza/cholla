#ifdef TIDES

#include "../global.h"
#include "../grid3D.h"
#include "tides.h"
#include "../io.h"
#include <math.h>

#ifdef MPI_CHOLLA
#include "../mpi_routines.h"
#endif

void Grid3D::updateCOM(){

	Real posx, posy, posz;
	Real velxStTemp = 0.;
	Real velyStTemp = 0.;
	Real velzStTemp = 0.;
	Real posxStTemp = 0.;
	Real posyStTemp = 0.;
	Real poszStTemp = 0.;
	Real totrhoTemp = 0.;

	#ifdef MPI_CHOLLA
	Real totrho;
	#endif

	Real rho;
	int i, j, k, id;

	for (k=H.n_ghost; k<H.nz-H.n_ghost; k++) {
		for (j=H.n_ghost; j<H.ny-H.n_ghost; j++) {
			for (i=H.n_ghost; i<H.nx-H.n_ghost; i++) {
				id = i + j*H.nx + k*H.nx*H.ny;

//			Get the centered cell positions at (i,j,k)
        Get_Position(i, j, k, &posx, &posy, &posz);
				rho = C.density[id];

//			Position of the center of mass
				totrhoTemp += rho;
				posxStTemp += posx * rho;
				posyStTemp += posy * rho;
				poszStTemp += posz * rho;

//			Velocity of the center of mass
				velxStTemp += C.momentum_x[id];
				velyStTemp += C.momentum_y[id];
				velzStTemp += C.momentum_z[id];

			}
		}
	}

  #ifdef MPI_CHOLLA
	totrho = ReduceRealSum(totrhoTemp);
	S.Mbox = totrho * H.dx * H.dy * H.dz;

  S.posSt[0] = ReduceRealSum(posxStTemp) / totrho;
  S.posSt[1] = ReduceRealSum(posyStTemp) / totrho;
  S.posSt[2] = ReduceRealSum(poszStTemp) / totrho;

	S.velSt[0] = ReduceRealSum(velxStTemp) / totrho;
	S.velSt[1] = ReduceRealSum(velyStTemp) / totrho;
	S.velSt[2] = ReduceRealSum(velzStTemp) / totrho;

  #else

	S.Mbox = totrhoTemp * H.dx * H.dy * H.dz;

  S.posSt[0] = posxStTemp / totrho;
  S.posSt[1] = posyStTemp / totrho;
  S.posSt[2] = poszStTemp / totrho;

	S.velSt[0] /= velxStTemp / totrho;
	S.velSt[1] /= velyStTemp / totrho;
	S.velSt[2] /= velzStTemp / totrho;

	#endif//MPI_CHOLLA

}

void Star::update(Real t, Real dt){

//The time along the orbit is different from the hydro time because we relax the star, and because the time along the orbit is measured with t = 0 at periapsis. When the star is not relaxed, we hold the star at the initial position along the orbit (but we compute no tidal forces!). After it's relaxed, we compute the time along the orbit accounting for the offset from the initial conditions. Remember that t0 is negative.

//When relaxed == 0, the tidal tensors aren't used anyways.

	if ( relaxed == 0 ){
		tOrb = t0;
	}
	else{
		tOrb = t + t0;
	}

// Important: first do frames, then tidal tensors since they depend on the frames!
	eta = geteta(tOrb);
	updateFrameCoords (tOrb, dt);
	updateBhCoords    (tOrb, dt);
	updateTidalTensors(tOrb, dt);

}

Real Star::getEorb(){

	Real mu = Mbh * Mstar / ( Mbh + Mstar );

	return 0.5 * mu * ( pow(velFrame[0] - velBh[0], 2.) + pow(velFrame[1] - velBh[1], 2.) ) - G_CGS * Mbh * Mstar / sqrt( pow(posFrame[0] - posBh[0],2.) + pow(posFrame[1] - posBh[1],2.) );

}

Real Star::geteta(Real t){

	Real num;
	Real den;
	Real totdyn = t / tdynOrb;
	Real aux = 3. * totdyn + sqrt( 8. + 9. * totdyn * totdyn );

	num = - 2. + pow( aux , 2./3.);
	den = sqrt(2.) * pow(aux, 1./3.);

	return num / den;

}

Real Star::getdeta(Real t){

	Real totdyn = t / tdynOrb;
	Real aux = 3. * totdyn + sqrt( 8. + 9. * totdyn * totdyn );
	Real num = 2. + pow(aux, 2./3.);
	Real den = sqrt( 16. + 18. * totdyn * totdyn ) * pow(aux, 1./3.);

	return num / den / tdynOrb;
	
}

Real Star::getddeta(Real t){

	Real totdyn = t / tdynOrb;
	Real aux0 = 8. + 9. * totdyn * totdyn;
	Real aux1 = 3. * totdyn + sqrt( 8. + 9. * totdyn * totdyn);

	Real prefac = pow(aux1,1./3.);
	Real term1 = 2. * ( sqrt(aux0) - 2. * pow(aux1, 1./3.) );
	Real term2 = 3. * totdyn * ( -6. - ( - 3. * totdyn + sqrt(aux0) ) * pow(aux1, 1./3.));
	Real num = prefac * ( term1 + term2 );
	Real den = 2. * sqrt(2.) * pow(aux0, 3./2.);

	return num / den / tdynOrb / tdynOrb;

}

// Updates the coordinates of the frame for t and t + dt / 2
void Star::updateFrameCoords(Real t, Real dt){

	Real deta    = getdeta (t);
	Real ddeta   = getddeta(t);
	Real MstFrac = Mbh / ( Mstar + Mbh );

	Real etaExt   = geteta  (t + dt / 2.);
	Real detaExt  = getdeta (t + dt / 2.);
	Real ddetaExt = getddeta(t + dt / 2.);

	posFrame[0] = MstFrac * rp * ( 1. - eta * eta );
	posFrame[1] = MstFrac * rp * 2. * eta;
	posFrame[2] = 0.;

	velFrame[0] = MstFrac * (-2.) * rp * eta * deta;
	velFrame[1] = MstFrac *   2.  * rp * deta;
	velFrame[2] = 0.;

	accFrame[0] = MstFrac * (-2.) * rp * ( deta * deta + eta * ddeta );
	accFrame[1] = MstFrac *   2.  * rp * ddeta;
	accFrame[2] = 0.;

	posFrameExt[0] = MstFrac * rp * ( 1. - etaExt * etaExt );
	posFrameExt[1] = MstFrac * rp * 2. * etaExt;
	posFrameExt[2] = 0.;

	velFrameExt[0] = MstFrac * (-2.) * rp * etaExt * detaExt;
	velFrameExt[1] = MstFrac *   2.  * rp * detaExt;
	velFrameExt[2] = 0.;

	accFrameExt[0] = MstFrac * (-2.) * rp * ( detaExt * detaExt + etaExt * ddetaExt );
	accFrameExt[1] = MstFrac *   2.  * rp * ddetaExt;
	accFrameExt[2] = 0.;

}
/*
void Star::getFrameAndBhPos(Real t, Real *xFrame, Real *yFrame, Real *zFrame, Real *xBh, Real *yBh, Real *zBh){

	Real eta = geteta(t);
	Real MstFrac = Mbh   / ( Mstar + Mbh );
	Real MbhFrac = Mstar / ( Mstar + Mbh );

	*xFrame = MstFrac * rp * ( 1. - eta * eta );
	*yFrame = MstFrac * rp * 2. * eta;
	*zFrame = MstFrac * 0.;

	*xBh = MbhFrac * rp * ( 1. - eta * eta );
	*yBh = MbhFrac * rp * 2. * eta;
	*zBh = MbhFrac * 0.;

}

// Given the density in the cells and the position of the black hole, this function returns the three components of the acceleration of the black hole
void Grid3D::updateBhAcc(){

//These variables hold the positions of the cell centers and the distance between the cell center and the bh
	Real posx, posy, posz, r;
//Will hold density
	Real rho;

//We need the volume of the cell because we compute the acceleration between two point masses: the BH and a mass rho * dV at the center of the cell
	Real dV = H.dx * H.dy * H.dz;

// To compute the acceleration the BH experiences, we sum over the acceleration caused by every cell in the star.
	int id;
	Real accBhTemp[3];
	for (int k=H.n_ghost; k<H.nz-H.n_ghost; k++) {
		for (int j=H.n_ghost; j<H.ny-H.n_ghost; j++) {
			for (int i=H.n_ghost; i<H.nx-H.n_ghost; i++) {
        id = i + j*H.nx + k*H.nx*H.ny;
				Get_Position(i, j, k, &posx, &posy, &posz);

				r = sqrt( ( posx - S.posBh[0] ) * ( posx - S.posBh[0] ) + ( posy - S.posBh[1] ) * ( posy - S.posBh[1] ) + ( posz - S.posBh[2] ) * ( posz - S.posBh[2] ) );
				rho = C.density[id];

				accBhTemp[0] += rho * ( posx - posBh[0] ) / pow(r, 3.);
				accBhTemp[1] += rho * ( posy - posBh[1] ) / pow(r, 3.);
				accBhTemp[2] += rho * ( posz - posBh[2] ) / pow(r, 3.);
			}
		}
	}

	accBhTemp[0] *= dV * G_CGS;
	accBhTemp[1] *= dV * G_CGS;
	accBhTemp[2] *= dV * G_CGS;

	#ifdef MPI_CHOLLA
	S.accBh[0] = ReduceRealSum(accBhTemp[0]);
	S.accBh[1] = ReduceRealSum(accBhTemp[1];
	S.accBh[2] = ReduceRealSum(accBhTemp[2]);
	#else
	S.accBh[0] = accBhTemp[0];
	S.accBh[1] = accBhTemp[1];
	S.accBh[2] = accBhTemp[2];
	#endif

}

*/

void Star::updateBhCoords(Real t, Real dt){

// Case with hardcoded BH trajectory

	Real deta    = getdeta (t);
	Real ddeta   = getddeta(t);
	Real etaExt   = geteta  (t + dt / 2.);
	Real detaExt  = getdeta (t + dt / 2.);
	Real ddetaExt = getddeta(t + dt / 2.);
	Real MbhFrac = Mstar / ( Mstar + Mbh );

	posBh[0] = - MbhFrac * rp * ( 1. - eta * eta );
	posBh[1] = - MbhFrac * rp * 2. * eta;
	posBh[2] = 0.;

	velBh[0] = - MbhFrac * (-2.) * rp * eta * deta;
	velBh[1] = - MbhFrac *   2.  * rp * deta;
	velBh[2] = 0.;

	accBh[0] = - MbhFrac * (-2.) * rp * ( deta * deta + eta * ddeta ) ;
	accBh[1] = - MbhFrac *   2.  * rp * ddeta;
	accBh[2] = 0.;

	posBhExt[0] = - MbhFrac * rp * ( 1. - etaExt * etaExt );
	posBhExt[1] = - MbhFrac * rp * 2. * etaExt;
	posBhExt[2] = 0.;

	velBhExt[0] = - MbhFrac * (-2.) * rp * etaExt * detaExt;
	velBhExt[1] = - MbhFrac *   2.  * rp * detaExt;
	velBhExt[2] = 0.;

	accBhExt[0] = - MbhFrac * (-2.) * rp * ( detaExt * detaExt + etaExt * ddetaExt ) ;
	accBhExt[1] = - MbhFrac *   2.  * rp * ddetaExt;
	accBhExt[2] = 0.;

}

/*
void Grid3D::IntegrateOrbit(){

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
