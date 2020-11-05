#ifdef TIDES

#include "tides.h"
#include "../global.h"
#include "../io.h"
#include <string.h>

// Kronecker delta
int kronDelta(int i, int j){

  if ( i == j ) return 1;
  else{
    return 0;
  }

}

void Star::initialize(struct parameters &P, Real t, Real dt, int nx, int ny, int nz){

  Mstar   = P.Mstar;
  Mbh     = P.Mbh;
  Rstar   = P.Rstar;
  polyN   = P.polyN;

//We currently assume that simulations are restarted after the relaxation
  if (strcmp(P.init, "Read_Grid")==0){
    relaxed = 1;
  }
  else{
    relaxed = 0;
  }

//Mass ratio
  q = Mbh / Mstar;

//Reduced mass
  mu = Mbh * Mstar / ( Mbh + Mstar );

//Tidal radius
  rt = Rstar * pow(q, 1./3.);

//Pericenter distance
  rp = P.rprt * rt;

//Initial eta
  eta0 = - sqrt(P.r0rt / P.rprt - 1.);

//Initial distance
  r0 = P.r0rt * rt;

//Dynamical time of the star and of the orbit
  tdynStar = sqrt( pow(Rstar, 3.) / G_CGS / Mstar           );
  tdynOrb  = sqrt( pow(rp   , 3.) / G_CGS / ( Mstar + Mbh ) );

//Initial time
  t0 = sqrt(2.) * tdynOrb * eta0 * ( 1. + eta0 * eta0 / 3. );

// Total energy of the star
  E0star = ( G_CGS * Mstar * Mstar / Rstar ) * ( 3. / ( polyN - 5. ) + 1. / ( 5. - polyN ) / ( P.gamma - 1. ) );

//Relaxation time
  tRelax = P.tRelaxtDyn * tdynStar;
  relaxRate0 = P.relaxRate0;
  relaxRateBkgnd = P.relaxRateBkgnd;

  update(t, dt);

  comBlocks = ceil ( nx * ny * nz / COMTPB );
  bufferxstar = (Real *) malloc( sizeof(Real) * comBlocks * 3);
  buffervstar = (Real *) malloc( sizeof(Real) * comBlocks * 3);

  chprintf("  Star:\n");
  chprintf("    Mass  : %.10e g\n", Mstar);
  chprintf("    Radius: %.10e cm\n", Rstar);
  chprintf("    n_poly: %.10e\n", polyN);
  chprintf("    t_dyn : %.10e s\n", tdynStar);

  chprintf("  Orbit:\n");
  chprintf("    Mass ratio    : %.10e\n", q);
  chprintf("    t_dyn         : %.10e s\n", tdynOrb);
  chprintf("    Tidal radius  : %.10e cm\n", rt      );
  chprintf("    Initial dist  : %.10e cm\n", r0      );
  chprintf("    Periapsis dist: %.10e cm\n", rp      );
  chprintf("    Periapsis time: %.10e s\n", -t0     );

  if ( tRelax > 0 ) chprintf("  Relaxation enabled. Initial relax rate: %f. Background relax rate: %.f\n", relaxRate0, relaxRateBkgnd);

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



// Given a position and a set of tidal tensors, return the tidal potential
Real Star::getTidalPotential(Real x, Real y, Real z, Real Cij[3][3], Real Cijk[3][3][3], Real Cijkl[3][3][3][3]){

//Coordinates where the potential is requested
  Real coords[3];
  coords[0] = x;
  coords[1] = y;
  coords[2] = z;

  Real tidalPot;

// Using tidal tensors
  tidalPot = 0.;
  for ( int i = 0; i < 3; i++ ){
    for ( int j = 0; j < 3; j++){
      tidalPot += 0.5 * Cij[i][j] * coords[i] * coords[j];
      for ( int k = 0; k < 3; k++){
        tidalPot += (1./6.) + Cijk[i][j][k] * coords[i] * coords[j] * coords[k];
        for ( int l = 0; l < 3; l++){
          tidalPot += (1./24.) * Cijkl[i][j][k][l] * coords[i] * coords[j] * coords[k] * coords[l];
        }
      }
    }
  }

// Using the exact Newtonian potential
//  Real r0 = sqrt(  );
//  tidalPot = - G_CGS * Mbh / rOrb

  return tidalPot;
}

// Updates the tidal tensors, which only depend on the position of the center of the frame.
// Updates tensors for t and t + dt / 2, since the latter will be used in the extrapolated potential.
void Star::updateTidalTensors(Real t, Real dt){

//Coordinates
  Real r  = sqrt( pow( posFrame[0] - posBh[0], 2. ) + pow( posFrame[1] - posBh[1], 2. ) + pow( posFrame[2] - posBh[2], 2. ));
  Real r2 = r * r;
  Real r3 = r2 * r;
  Real r4 = r3 * r;
  Real r5 = r4 * r;

  Real rExt  = sqrt( pow( posFrameExt[0] - posBhExt[0], 2. ) + pow( posFrameExt[1] - posBhExt[1], 2. ) + pow( posFrameExt[2] - posBhExt[2], 2. ));
  Real r2Ext = rExt  * rExt;
  Real r3Ext = r2Ext * rExt;
  Real r4Ext = r3Ext * rExt;
  Real r5Ext = r4Ext * rExt;

  for ( int i = 0; i < 3; i++ ){
    for ( int j = 0; j < 3; j++ ){

//    Quadrupole tensor at t
      Cij[i][j] = kronDelta(i, j) - 3. * posFrame[i] * posFrame[j] / r2;
      Cij[i][j] *= G_CGS * Mbh / r3;

//    Quadrupole tensor at t + dt/2
      extCij[i][j] = kronDelta(i, j) - 3. * posFrameExt[i] * posFrameExt[j] / r2Ext;
      extCij[i][j] *= G_CGS * Mbh / r3Ext;

      for ( int k = 0; k < 3; k++){

//      Octupole tensor at t
        Cijk[i][j][k]  = 15. * posFrame[i] * posFrame[j] * posFrame[k] / r3
                         - 3.  * ( posFrame[i] * kronDelta(j, k) + posFrame[j] * kronDelta(i, k) + posFrame[k] * kronDelta(i, j) ) / r;
        Cijk[i][j][k] *= G_CGS * Mbh / r4;

//      Octupole tensor at t + dt / 2
        extCijk[i][j][k]  = 15. * posFrameExt[i] * posFrameExt[j] * posFrameExt[k] / r3Ext
                         - 3.  * ( posFrameExt[i] * kronDelta(j, k) + posFrameExt[j] * kronDelta(i, k) + posFrameExt[k] * kronDelta(i, j) ) / rExt;
        extCijk[i][j][k] *= G_CGS * Mbh / r4Ext;

        for ( int l = 0; l < 3; l++){

//        Hexadecapole tensor
          Cijkl[i][j][k][l]  = - 105. * posFrame[i] * posFrame[j] * posFrame[k] * posFrame[l] / r4
                               + 15. * (   kronDelta(i, l) * posFrame[j] * posFrame[k]
                                     + kronDelta(j, l) * posFrame[i] * posFrame[k]
                                     + kronDelta(k, l) * posFrame[i] * posFrame[j]
                                     + kronDelta(i, j) * posFrame[k] * posFrame[l]
                                     + kronDelta(j, k) * posFrame[i] * posFrame[l]
                                     + kronDelta(i, k) * posFrame[j] * posFrame[l]
                                       ) / r2
                                - 3. *  (   kronDelta(i, j) * kronDelta(k, l)
                                     + kronDelta(j, k) * kronDelta(i, l)
                                     + kronDelta(i, k) * kronDelta(j, l)
                                        );
          Cijkl[i][j][k][l] *= G_CGS * Mbh / r5;

//        Hexadecapole tensor at t + dt / 2
          extCijkl[i][j][k][l]  = - 105. * posFrameExt[i] * posFrameExt[j] * posFrameExt[k] * posFrameExt[l] / r4Ext
                               + 15. * (   kronDelta(i, l) * posFrameExt[j] * posFrameExt[k]
                                     + kronDelta(j, l) * posFrameExt[i] * posFrameExt[k]
                                     + kronDelta(k, l) * posFrameExt[i] * posFrameExt[j]
                                     + kronDelta(i, j) * posFrameExt[k] * posFrameExt[l]
                                     + kronDelta(j, k) * posFrameExt[i] * posFrameExt[l]
                                     + kronDelta(i, k) * posFrameExt[j] * posFrameExt[l]
                                       ) / r2Ext
                                - 3. *  (   kronDelta(i, j) * kronDelta(k, l)
                                     + kronDelta(j, k) * kronDelta(i, l)
                                     + kronDelta(i, k) * kronDelta(j, l)
                                        );
          extCijkl[i][j][k][l] *= G_CGS * Mbh / r5Ext;

        }
      }
    }
  }

}

#endif
