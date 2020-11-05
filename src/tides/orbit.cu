#ifdef TIDES

#include "../global.h"
#include "../grid3D.h"
#include "tides.h"
#include "../io.h"
#include <math.h>

#define COMTPB (1024)

#ifdef MPI_CHOLLA
#include "../mpi_routines.h"
#endif

__global__ void comKernel(Real *rho, Real *momentum_x, Real *momentum_y, Real *momentum_z, Real *bounds, Real *dx, int *n, int n_ghost, Real *partialxstar, Real *partialvstar){

  __shared__ Real xstar[COMTPB * 3];
  __shared__ Real vstar[COMTPB * 3];

  for ( int i = 3 * threadIdx.x; i < 3 * ( threadIdx.x + 1 ); i++ ){
    xstar[i] = 0.;
    vstar[i] = 0.;
  }

  int nreal[3];
  for ( int i = 0; i < 3; i++ ) nreal[i] = n[i] - 2 * n_ghost;
  int nrealcells = nreal[0] * nreal[1] * nreal[2];

  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int tid_z = tid / ( nreal[0] * nreal[1] );
  int tid_y = ( tid - tid_z * nreal[0] * nreal[1] ) / nreal[0];
  int tid_x = tid - tid_z * nreal[0] * nreal[1] - tid_y * nreal[0];

  Real x[3];
  int fakeid;

  while ( tid < nrealcells ){

    tid_z = tid / ( nreal[0] * nreal[1] );
    tid_y = ( tid - tid_z * nreal[0] * nreal[1] ) / nreal[0];
    tid_x = tid - tid_z * nreal[0] * nreal[1] - tid_y * nreal[0];
    fakeid = ( tid_z + n_ghost ) * n[0] * n[1] + ( tid_y + n_ghost ) * n[0] + ( tid_x + n_ghost );

    x[0] = bounds[0] + dx[0] * ( tid_x + 0.5);
    x[1] = bounds[1] + dx[1] * ( tid_y + 0.5);
    x[2] = bounds[2] + dx[2] * ( tid_z + 0.5);

//  Position of the center of mass
    for ( int ii = 0; ii < 3; ii++ ) xstar[threadIdx.x * 3 + ii] += x[ii] * rho[fakeid];

//  Velocity of the center of mass
    vstar[threadIdx.x * 3    ] += momentum_x[fakeid];
    vstar[threadIdx.x * 3 + 1] += momentum_y[fakeid];
    vstar[threadIdx.x * 3 + 2] += momentum_z[fakeid];

    tid += blockDim.x * gridDim.x;
  }

  __syncthreads();

  int i = blockDim.x / 2;
  while ( i > 0 ){
    if ( threadIdx.x < i){
      for ( int ii = 0; ii < 3; ii++ ){
        xstar[threadIdx.x * 3 + ii] += xstar[( threadIdx.x + i) * 3 + ii];
        vstar[threadIdx.x * 3 + ii] += vstar[( threadIdx.x + i) * 3 + ii];
      }
    }
    __syncthreads();
    i /= 2;
  }

  if ( threadIdx.x == 0 ){
    for ( int ii = 0; ii < 3; ii++){
      partialxstar[3 * blockIdx.x + ii] = xstar[ii];
      partialvstar[3 * blockIdx.x + ii] = vstar[ii];
    }
  }

}

void Grid3D::updateCOM(){
  S.Mbox = Grav.ReQ[0] * sqrt( 4 * M_PI );
  Real totrho = S.Mbox / H.dx / H.dy / H.dz;
/*
  for ( int i = 0; i < 3; i++ ) S.xstar[i] = 0.;
  for ( int i = 0; i < 3; i++ ) S.vstar[i] = 0.;

  Real rho, x[3];
  int id;

  for (int k = H.n_ghost; k<H.nz-H.n_ghost; k++) {
    for (int j = H.n_ghost; j<H.ny-H.n_ghost; j++) {
      for (int i = H.n_ghost; i<H.nx-H.n_ghost; i++) {
        id = i + j*H.nx + k*H.nx*H.ny;

//      Get the centered cell positions at (i,j,k)
        Get_Position(i, j, k, &x[0], &x[1], &x[2]);
        rho = C.density[id];
//      Position of the center of mass
        for ( int ii = 0; ii < 3; ii++ ) S.xstar[ii] += x[ii] * rho;

//      Velocity of the center of mass
        S.vstar[0] += C.momentum_x[id];
        S.vstar[1] += C.momentum_y[id];
        S.vstar[2] += C.momentum_z[id];

      }
    }
  }

  #ifdef MPI_CHOLLA
  MPI_Allreduce(MPI_IN_PLACE, S.xstar, 3, MPI_CHREAL, MPI_SUM, world);
  MPI_Allreduce(MPI_IN_PLACE, S.vstar, 3, MPI_CHREAL, MPI_SUM, world);
  #endif

  Real vstarslow[3], xstarslow[3];
  for ( int i = 0; i < 3; i++ ){
    vstarslow[i] = S.vstar[i] / totrho;
    xstarslow[i] = S.xstar[i] / totrho;
  }
*/
  Real dx[3], bounds[3];
  int n[3];
  dx[0] = H.dx; dx[1] = H.dy; dx[2] = H.dz;
  bounds[0] = H.xblocal; bounds[1] = H.yblocal; bounds[2] = H.zblocal;
  n[0] = H.nx; n[1] = H.ny; n[2] = H.nz;
  Real *dev_rho, *dev_momentum_x, *dev_momentum_y, *dev_momentum_z, *dev_bounds, *dev_dx, *dev_partialxstar, *dev_partialvstar;
  int *dev_n;

//Allocate memory in GPU
  cudaMalloc( (void**)&dev_rho          , n[0] * n[1] * n[2] *sizeof(Real) );
  cudaMalloc( (void**)&dev_momentum_x   , n[0] * n[1] * n[2] *sizeof(Real) );
  cudaMalloc( (void**)&dev_momentum_y   , n[0] * n[1] * n[2] *sizeof(Real) );
  cudaMalloc( (void**)&dev_momentum_z   , n[0] * n[1] * n[2] *sizeof(Real) );
  cudaMalloc( (void**)&dev_bounds       , 3 * sizeof(Real)                 );
  cudaMalloc( (void**)&dev_n            , 3 * sizeof(int)                  );
  cudaMalloc( (void**)&dev_dx           , 3 * sizeof(Real)                 );
  cudaMalloc( (void**)&dev_partialxstar , S.comBlocks * 3 * sizeof(Real)   );
  cudaMalloc( (void**)&dev_partialvstar , S.comBlocks * 3 * sizeof(Real)   );

//Copy inputs to GPU
  cudaMemcpy( dev_rho       , C.density   , n[0] * n[1] * n[2] * sizeof(Real), cudaMemcpyHostToDevice);
  cudaMemcpy( dev_momentum_x, C.momentum_x, n[0] * n[1] * n[2] * sizeof(Real), cudaMemcpyHostToDevice);
  cudaMemcpy( dev_momentum_y, C.momentum_y, n[0] * n[1] * n[2] * sizeof(Real), cudaMemcpyHostToDevice);
  cudaMemcpy( dev_momentum_z, C.momentum_z, n[0] * n[1] * n[2] * sizeof(Real), cudaMemcpyHostToDevice);
  cudaMemcpy( dev_bounds    , bounds      , 3*sizeof(Real), cudaMemcpyHostToDevice);
  cudaMemcpy( dev_n         , n           , 3*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy( dev_dx        , dx          , 3*sizeof(Real), cudaMemcpyHostToDevice);

//Call Kernel
  cudaDeviceSynchronize();
  comKernel<<<S.comBlocks,COMTPB>>>(dev_rho, dev_momentum_x, dev_momentum_y, dev_momentum_z, dev_bounds, dev_dx, dev_n, H.n_ghost, dev_partialxstar, dev_partialvstar);

//Copy result to CPU
  cudaMemcpy(S.bufferxstar, dev_partialxstar, sizeof(Real) * S.comBlocks * 3, cudaMemcpyDeviceToHost);
  cudaMemcpy(S.buffervstar, dev_partialvstar, sizeof(Real) * S.comBlocks * 3, cudaMemcpyDeviceToHost);

//Free GPU
  cudaFree(dev_rho);
  cudaFree(dev_momentum_x);
  cudaFree(dev_momentum_y);
  cudaFree(dev_momentum_z);
  cudaFree(dev_n);
  cudaFree(dev_bounds);
  cudaFree(dev_dx);
  cudaFree(dev_partialxstar);
  cudaFree(dev_partialvstar);

  for ( int ii = 0; ii < 3; ii++ ){
    S.vstar[ii] = 0.;
    S.xstar[ii] = 0.;
    for ( int i = 0; i < S.comBlocks; i++ ){
      S.xstar[ii] += S.bufferxstar[3*i + ii];
      S.vstar[ii] += S.buffervstar[3*i + ii];
    }
  }

  for ( int i = 0; i < 3; i++ ){
    S.xstar[i] /= totrho;
    S.vstar[i] /= totrho;
  }

  #ifdef MPI_CHOLLA
  MPI_Allreduce(MPI_IN_PLACE, S.xstar, 3, MPI_CHREAL, MPI_SUM, world);
  MPI_Allreduce(MPI_IN_PLACE, S.vstar, 3, MPI_CHREAL, MPI_SUM, world);
  #endif
/*
  Real xeps[3];
  Real veps[3];
  for ( int i = 0; i < 3; i++ ){
    xeps[i] = S.xstar[i] / xstarslow[i] - 1.;
    veps[i] = S.vstar[i] / vstarslow[i] - 1.;
  }

  chprintf("xstar new: %.10e, %.10e, %.10e\n", S.xstar[0], S.xstar[1], S.xstar[2]);
  chprintf("xstar old: %.10e, %.10e, %.10e\n", xstarslow[0], xstarslow[1], xstarslow[2]);
  chprintf("vstar new: %.10e, %.10e, %.10e\n", S.vstar[0], S.vstar[1], S.vstar[2]);
  chprintf("vstar old: %.10e, %.10e, %.10e\n", vstarslow[0], vstarslow[1], vstarslow[2]);
  chprintf("xstar relative error: %.16e, %.16e, %.16e\n", xeps[0], xeps[1], xeps[2]);
  chprintf("vstar relative error: %.16e, %.16e, %.16e\n", veps[0], veps[1], veps[2]);
*/
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

#endif
