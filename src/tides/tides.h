#ifdef TIDES

#ifndef TIDES_H
#define TIDES_H

#include <stdio.h>
#include <cmath>
#include "../global.h"
//#include "../grid3D.h"

class Star
{
public:


//Star
	Real Mstar;
	Real Rstar;
	Real polyN;
	Real tRelax;
	Real tdynStar;
	int relaxed;

//Total mass in the box. Used to compute mass loss for close encounters
	Real Mbox;

//Tides
	Real rp;
	Real r0;
	Real tdynOrb;
	Real Mbh;
	Real q;
	Real mu;
	Real rt;
	Real t0;
	Real tOrb;
	Real eta0;
	Real eta;
	Real E0orb;
	Real Eorb;
	Real E0star;
	Real DEoE;

//Coordinates of the center of the box
	Real posFrame[3];
	Real velFrame[3];
	Real accFrame[3];

//Extrapolated
	Real posFrameExt[3];
	Real velFrameExt[3];
	Real accFrameExt[3];

//	Coordinates of the bh
	Real posBh[3];
	Real velBh[3];
	Real accBh[3];

//Extrapolated
	Real posBhExt[3];
	Real velBhExt[3];
	Real accBhExt[3];

//	Coordinates of the star
	Real posSt[3];
	Real velSt[3];
	Real accSt[3];

	Real posStExt[3];
	Real velStExt[3];
	Real accStExt[3];

//Tidal tensors at the current time and at t + dt / 2

	Real extCij[3][3];
	Real extCijk[3][3][3];
	Real extCijkl[3][3][3][3];

	Real Cij[3][3];
	Real Cijk[3][3][3];
	Real Cijkl[3][3][3][3];

//Functions that change the state of S
	void initialize(struct parameters &P, Real t, Real dt);
	void update(Real t, Real dt);
	void updateFrameCoords(Real t, Real dt);
	void updateBhCoords(Real t, Real dt);
	void updateTidalTensors(Real t, Real dt);

//Value of eta (proxy for time) and its first two derivatives with respect to time. These are used to track the coordinates of the center of the frame at all times analytically
	Real geteta(Real t);
	Real getdeta(Real t);
	Real getddeta(Real t);

//Returns the tidal potential given a set of tidal tensors
	Real getTidalPotential(Real x, Real y, Real z, Real Cij[3][3], Real Cijk[3][3][3], Real Cijkl[3][3][3][3]);

//Returns the position of the center of the frame at time t
	void getFrameAndBhPos(Real t, Real *xFrame, Real *yFrame, Real *zFrame, Real *xBh, Real *yBh, Real *zBh);

	Real getEorb();

};

#endif
#endif
