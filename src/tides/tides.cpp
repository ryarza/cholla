#ifdef TIDES

#include "tides.h"
#include "../global.h"
#include "../io.h"

// Kronecker delta
Real kronDelta(int i, int j){

  if ( i == j ) return 1.;
  else return 0.;

}

void Star::initialize(struct parameters P, Real t, Real dt, int nx, int ny, int nz){

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

  chprintf(" Star:\n");
  chprintf("  Mass  : %.10e g\n", Mstar);
  chprintf("  Radius: %.10e cm\n", Rstar);
  chprintf("  n_poly: %.10e\n", polyN);
  chprintf("  t_dyn : %.10e s\n", tdynStar);

  chprintf(" Orbit:\n");
  chprintf("  Mass ratio    : %.10e\n", q);
  chprintf("  t_dyn         : %.10e s\n", tdynOrb);
  chprintf("  Tidal radius  : %.10e cm\n", rt      );
  chprintf("  Initial dist  : %.10e cm\n", r0      );
  chprintf("  Periapsis dist: %.10e cm\n", rp      );
  chprintf("  Periapsis time: %.10e s\n", -t0     );

  if ( tRelax > 0 ) chprintf(" Relaxation enabled. Initial relax rate: %f. Background relax rate: %.f\n", relaxRate0, relaxRateBkgnd);

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
  updateTidalTensors();

}



// Given a position and a set of tidal tensors, return the tidal potential
Real Star::getTidalPotential(Real *x, Real argCij[3][3], Real argCijk[3][3][3], Real argCijkl[3][3][3][3]){

// Using tidal tensors
  Real tidalPot = 0.;
  for ( int i = 0; i < 3; i++ ){
    for ( int j = 0; j < 3; j++){
      tidalPot += 0.5 * argCij[i][j] * x[i] * x[j];
      for ( int k = 0; k < 3; k++){
        tidalPot += (1./6.) * argCijk[i][j][k] * x[i] * x[j] * x[k];
        for ( int l = 0; l < 3; l++){
          tidalPot += (1./24.) * argCijkl[i][j][k][l] * x[i] * x[j] * x[k] * x[l];
        }
      }
    }
  }

  return tidalPot;

}

// Updates the tidal tensors, which only depend on the position of the center of the frame.
// Updates tensors for t and t + dt / 2, since the latter will be used in the extrapolated potential.
void Star::updateTidalTensors(){

  Real r2 = pow( posFrame[0] - posBh[0], 2. ) + pow( posFrame[1] - posBh[1], 2. ) + pow( posFrame[2] - posBh[2], 2. );
  Real r  = sqrt(r2);
  Real r3 = r2 * r;
  Real r4 = r3 * r;
  Real r5 = r4 * r;

  Real r2Ext = pow( posFrameExt[0] - posBhExt[0], 2. ) + pow( posFrameExt[1] - posBhExt[1], 2. ) + pow( posFrameExt[2] - posBhExt[2], 2. );
  Real rExt  = sqrt( r2Ext );
  Real r3Ext = r2Ext * rExt;
  Real r4Ext = r3Ext * rExt;
  Real r5Ext = r4Ext * rExt;

  Real bigx[3], bigxExt[3];
  for ( int i = 0; i < 3; i++ ) bigx[i] = posFrame[i] - posBh[i];
  for ( int i = 0; i < 3; i++ ) bigxExt[i] = posFrameExt[i] - posBhExt[i];

  for ( int i = 0; i < 3; i++ ){
    for ( int j = 0; j < 3; j++ ){

//    Quadrupole tensor at t
      Cij[i][j] = kronDelta(i, j) - 3. * bigx[i] * bigx[j] / r2;
      Cij[i][j] *= G_CGS * Mbh / r3;

//    Quadrupole tensor at t + dt/2
      extCij[i][j] = kronDelta(i, j) - 3. * bigxExt[i] * bigxExt[j] / r2Ext;
      extCij[i][j] *= G_CGS * Mbh / r3Ext;

      for ( int k = 0; k < 3; k++){

//      Octupole tensor at t
        Cijk[i][j][k]  = 15. * bigx[i] * bigx[j] * bigx[k] / r3
                         - 3.  * ( bigx[i] * kronDelta(j, k) + bigx[j] * kronDelta(i, k) + bigx[k] * kronDelta(i, j) ) / r;
        Cijk[i][j][k] *= G_CGS * Mbh / r4;

//      Octupole tensor at t + dt / 2
        extCijk[i][j][k]  = 15. * bigxExt[i] * bigxExt[j] * bigxExt[k] / r3Ext
                         - 3.  * ( bigxExt[i] * kronDelta(j, k) + bigxExt[j] * kronDelta(i, k) + bigxExt[k] * kronDelta(i, j) ) / rExt;
        extCijk[i][j][k] *= G_CGS * Mbh / r4Ext;

        for ( int l = 0; l < 3; l++){

//        Hexadecapole tensor
          Cijkl[i][j][k][l]  = - 105. * bigx[i] * bigx[j] * bigx[k] * bigx[l] / r4
                               +  15. * (   kronDelta(i, l) * bigx[j] * bigx[k]
                                          + kronDelta(j, l) * bigx[i] * bigx[k]
                                          + kronDelta(k, l) * bigx[i] * bigx[j]
                                          + kronDelta(i, j) * bigx[k] * bigx[l]
                                          + kronDelta(j, k) * bigx[i] * bigx[l]
                                          + kronDelta(i, k) * bigx[j] * bigx[l]
                                        ) / r2
                                -  3. * (   kronDelta(i, j) * kronDelta(k, l)
                                          + kronDelta(j, k) * kronDelta(i, l)
                                          + kronDelta(i, k) * kronDelta(j, l)
                                       );
          Cijkl[i][j][k][l] *= G_CGS * Mbh / r5;

//        Hexadecapole tensor at t + dt / 2
          extCijkl[i][j][k][l]  = - 105. * bigxExt[i] * bigxExt[j] * bigxExt[k] * bigxExt[l] / r4Ext
                               + 15. * (   kronDelta(i, l) * bigxExt[j] * bigxExt[k]
                                     + kronDelta(j, l) * bigxExt[i] * bigxExt[k]
                                     + kronDelta(k, l) * bigxExt[i] * bigxExt[j]
                                     + kronDelta(i, j) * bigxExt[k] * bigxExt[l]
                                     + kronDelta(j, k) * bigxExt[i] * bigxExt[l]
                                     + kronDelta(i, k) * bigxExt[j] * bigxExt[l]
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

/*
  Real extcCij[3][3], extcCijk[3][3][3], extcCijkl[3][3][3][3];
  extcCij[0][0] = G_CGS * Mbh * ( rExt * rExt - 3. * bigxExt[0] * bigxExt[0] ) / r5Ext;
  extcCij[1][1] = G_CGS * Mbh * ( rExt * rExt - 3. * bigxExt[1] * bigxExt[1] ) / r5Ext;
  extcCij[2][2] = G_CGS * Mbh * ( rExt * rExt - 3. * bigxExt[2] * bigxExt[2] ) / r5Ext;

  extcCij[0][1] = - 3 * G_CGS * Mbh * bigxExt[0] * bigxExt[1] / r5Ext;
  chprintf("extcCij[0][1] = %.10e\n", extcCij[0][1] );
  extcCij[1][0] = extcCij[0][1];
  chprintf("extcCij[1][0] = %.10e\n", extcCij[1][0] );

  chprintf("extcCij[1][0] = %.10e\n", extcCij[1][0] );
  extcCij[0][2] = - 3 * G_CGS * Mbh * bigxExt[0] * bigxExt[2] / r5Ext;
  chprintf("extcCij[1][0] = %.10e\n", extcCij[1][0] );
  extcCij[2][0] = extcCij[0][2];
  chprintf("extcCij[1][0] = %.10e\n", extcCij[1][0] );

  extcCij[1][2] = - 3 * G_CGS * Mbh * bigxExt[1] * bigxExt[2] / r5Ext;
  chprintf("extcCij[1][0] = %.10e\n", extcCij[1][0] );
  extcCij[2][1] = extcCij[1][2];
  chprintf("extcCij[1][0] = %.10e\n", extcCij[1][0] );

  Real ratio;
  for ( int i = 0; i < 3; i++ ){
    for ( int j = 0; j < 3; j++ ){
      ratio = extCij[i][j] / extcCij[i][j];
      chprintf("C[%i][%i] / exact = %.10e/%.10e = %.10e\n", i, j, extCij[i][j], extcCij[i][j], ratio);
    }
  }
*/
/*
  extcCijk[0][0][0] = G_CGS * Mbh * ( - 9. * rExt * rExt * posFrameExt[0] + 15 * pow(posFrameExt[0], 3.) ) / pow(rExt, 7.);
  chprintf("C[1, 1, 1] / correct = %.10e\n", extCijk[0][0][0] / extcCijk[0][0][0]);
*/
//  chprintf("C[1, 1, 1, 1] / correct = %.10e\n", Cijk[0][0][0][0] / cCijkl[0][0][0][0]);

}

#endif
