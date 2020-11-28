#ifdef TIDES

#include "../global.h"
#include "../global_cuda.h"
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

__global__ void potBHKernel(Real *bounds, Real *dx, Real *xFrame, Real *xBH, Real Mbh, int *n, int n_ghost, Real *potBH){

  int nreal[3], tid[3], fid;
  for ( int i = 0; i < 3; i++ ) nreal[i] = n[i] - 2 * n_ghost;
  int nrealcells = nreal[0] * nreal[1] * nreal[2];

  int tid1d = threadIdx.x + blockIdx.x * blockDim.x;
  tid[2] = tid1d / ( nreal[0] * nreal[1] );
  tid[1] = ( tid1d - tid[2] * nreal[0] * nreal[1] ) / nreal[0];
  tid[0] = tid1d - tid[2] * nreal[0] * nreal[1] - tid[1] * nreal[0];

  Real x[3];

  while ( tid1d < nrealcells ){

    tid[2] = tid1d / ( nreal[0] * nreal[1] );
    tid[1] = ( tid1d - tid[2] * nreal[0] * nreal[1] ) / nreal[0];
    tid[0] = tid1d - tid[2] * nreal[0] * nreal[1] - tid[1] * nreal[0];
    fid = ( tid[2] + n_ghost ) * n[0] * n[1] + ( tid[1] + n_ghost ) * n[0] + ( tid[0] + n_ghost );

    for ( int i = 0; i < 3; i++ ) x[i] = bounds[i] + dx[i] * ( tid[i] + 0.5) + xFrame[i] - xBH[i];
    potBH[fid] = - G_CGS * Mbh / sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );

    tid1d += blockDim.x * gridDim.x;
  }

}

#ifdef TIDES_OUTPUT_POTENTIAL_BH
void Grid3D::updatePotBH(){
  
  Real dx[3], bounds[3];
  int n[3];
  dx[0] = H.dx; dx[1] = H.dy; dx[2] = H.dz;
  bounds[0] = H.xblocal; bounds[1] = H.yblocal; bounds[2] = H.zblocal;
  n[0] = H.nx; n[1] = H.ny; n[2] = H.nz;

  Real *dev_potBH, *dev_bounds, *dev_dx, *dev_xFrame, *dev_xBH;
  int *dev_n;

  CudaSafeCall( cudaMalloc( (void**)&dev_bounds, 3                  * sizeof(Real) ) );
  CudaSafeCall( cudaMalloc( (void**)&dev_n     , 3                  * sizeof(int ) ) );
  CudaSafeCall( cudaMalloc( (void**)&dev_dx    , 3                  * sizeof(Real) ) );
  CudaSafeCall( cudaMalloc( (void**)&dev_xFrame, 3                  * sizeof(Real) ) );
  CudaSafeCall( cudaMalloc( (void**)&dev_xBH   , 3                  * sizeof(Real) ) );
  CudaSafeCall( cudaMalloc( (void**)&dev_potBH , n[0] * n[1] * n[2] * sizeof(Real) ) );

  CudaSafeCall( cudaMemcpy( dev_bounds, bounds    , 3*sizeof(Real), cudaMemcpyHostToDevice) );
  CudaSafeCall( cudaMemcpy( dev_n     , n         , 3*sizeof(int ), cudaMemcpyHostToDevice) );
  CudaSafeCall( cudaMemcpy( dev_dx    , dx        , 3*sizeof(Real), cudaMemcpyHostToDevice) );
  CudaSafeCall( cudaMemcpy( dev_xFrame, S.posFrame, 3*sizeof(Real), cudaMemcpyHostToDevice) );
  CudaSafeCall( cudaMemcpy( dev_xBH   , S.posBh   , 3*sizeof(Real), cudaMemcpyHostToDevice) );
  
  potBHKernel<<<ceil(n[0]*n[1]*n[2]/1024),1024>>>(dev_bounds, dev_dx, dev_xFrame, dev_xBH, S.Mbh, dev_n, H.n_ghost, dev_potBH);
  CudaCheckError();

  CudaSafeCall( cudaMemcpy(C.Grav_potential_BH, dev_potBH, n[0] * n[1] * n[2] * sizeof(Real), cudaMemcpyDeviceToHost) );
  
  cudaFree(dev_n);
  cudaFree(dev_bounds);
  cudaFree(dev_dx);
  cudaFree(dev_xFrame);
  cudaFree(dev_xBH);
  cudaFree(dev_potBH);

}
#endif

void Grid3D::updateCOM(){

  S.Mbox = Grav.ReQ[0] * sqrt( 4 * M_PI );
  Real totrho = S.Mbox / H.dx / H.dy / H.dz;

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

  #ifdef OUTPUT_ALWAYS_COM
//Write the COM, etc to the logfile
  char *message = (char*)malloc(500 * sizeof(char));
// Column headers:  <--t--> <--------xstar--------> <--------vstar--------> <---------xbh---------> <---------vbh---------> <-------xFrame -------> <-------vFrame -------> <-------aFrame ------->
  sprintf(message, "%17.10e %17.10e %17.10e %17.10e %17.10e %17.10e %17.10e %17.10e %17.10e %17.10e %17.10e %17.10e %17.10e %17.10e %17.10e %17.10e %17.10e %17.10e %17.10e %17.10e %17.10e %17.10e", H.t, S.xstar[0], S.xstar[1], S.xstar[2], S.vstar[0], S.vstar[1], S.vstar[2], S.posBh[0], S.posBh[1], S.posBh[2], S.velBh[0], S.velBh[1], S.velBh[2], S.posFrame[0], S.posFrame[1], S.posFrame[2], S.velFrame[0], S.velFrame[1], S.velFrame[2], S.accFrame[0], S.accFrame[1], S.accFrame[2]);
  Write_Message_To_Log_File("orbit_evolution.log", message);
  free(message);
  #endif

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
