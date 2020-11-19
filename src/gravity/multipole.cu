#if ( defined TIDES || defined POISSON_TEST ) && defined CUDA && defined GRAVITY

#include "../global.h"
#include "../grid3D.h"
#include "grav3D.h"
#include "../io.h"
#include "../global_cuda.h"

#ifdef MPI_CHOLLA
#include "../mpi_routines.h"
#endif

#ifdef TIDES
void Grav3D::AllocateMemoryBoundaries_GPU(int nx, int ny, int nz){


  chprintf("Allocating GPU memory for boundaries\n");
  chprintf("n alloc = %i %i %i\n", nx, ny, nz);
  CudaSafeCall( cudaMalloc( (void**)&dev_rho            , nx * ny * nz                                  * sizeof(Real) ) );
  CudaSafeCall( cudaMalloc( (void**)&dev_center         , 3                                             * sizeof(Real) ) );
  CudaSafeCall( cudaMalloc( (void**)&dev_bounds         , 3                                             * sizeof(Real) ) );
  CudaSafeCall( cudaMalloc( (void**)&dev_dx             , 3                                             * sizeof(Real) ) );
  CudaSafeCall( cudaMalloc( (void**)&dev_partialReQ     , Qblocks * ( ( 1 + LMAX ) * ( 2 + LMAX ) / 2 ) * sizeof(Real) ) );
  CudaSafeCall( cudaMalloc( (void**)&dev_partialImQ     , Qblocks * ( ( 1 + LMAX ) * ( 2 + LMAX ) / 2 ) * sizeof(Real) ) );
  CudaSafeCall( cudaMalloc( (void**)&dev_partialCenter  , 3 * centerBlocks                              * sizeof(Real) ) );
  CudaSafeCall( cudaMalloc( (void**)&dev_partialTotrhosq, 3 * centerBlocks                              * sizeof(Real) ) );
  CudaSafeCall( cudaMalloc( (void**)&dev_n              , 3                                             * sizeof(int)  ) );

}

void Grav3D::FreeMemoryBoundaries_GPU(){

  cudaFree(dev_rho            );
  cudaFree(dev_center         );
  cudaFree(dev_bounds         );
  cudaFree(dev_dx             );
  cudaFree(dev_partialReQ     );
  cudaFree(dev_partialImQ     );
  cudaFree(dev_partialCenter  );
  cudaFree(dev_partialTotrhosq);
  cudaFree(dev_n              );

}
#endif

//The arrays we use for Legendre polynomials are 1D, so we need to do some index juggling to turn the tuple (thread number, l, m) into a single number
int Grav3D::Qidx(int cidx, int l, int m){

  if ( m > l || l > LMAX || m < 0 ){
    printf("Wrong parameters!\n");
    return -1;
  }

  int stride = ( 1 + LMAX ) * ( 2 + LMAX ) / 2;
  int substride = l * ( l + 1 ) / 2;
  return cidx * stride + substride + m;

}

//Same but for the device
__device__ int dQidx(int cidx, int l, int m){

  if ( m > l || l > LMAX || m < 0 ){
    printf("Wrong parameters!\n");
    return -1;
  }

  int stride = ( 1 + LMAX ) * ( 2 + LMAX ) / 2;
  int substride = l * ( l + 1 ) / 2;
  return cidx * stride + substride + m;

}

//Recursively computes Legendre polynomials up to order LMAX. Returns 1D array
__device__ void fillLegP(Real* legP, Real x)
{

  for ( int l = 0; l <= LMAX; l++ ){
    for ( int m = 0; m <= l; m++ ){
      legP[dQidx(0,l,m)] = 0.;
    }
  }

//Initial polynomial for recursion relations
  legP[dQidx(0,0,0)] = 1./sqrt(4.*M_PI);

//Diagonal
  for( int m = 1; m <= LMAX; m++)
  {
    legP[dQidx(0,m,m)] = - sqrt( 1. + 1. / 2. / m ) * sqrt( 1. - x * x ) * legP[dQidx(0,m-1,m-1)];
  }

  for( int m = 0; m < LMAX; m++)
  {
    legP[dQidx(0,m+1,m)] = sqrt( 2. * m + 3. ) * x * legP[dQidx(0,m,m)];
  }

  for( int m = 0; m <= LMAX; m++){
    for( int l = m + 2; l <= LMAX; l++){
      Real c1 = sqrt( ((2.0*l+1)*(2.0*l-1)) / ((l+m)*(l-m)));
      Real c2 = sqrt( (2.0*l+1)*(l-m-1.0)*(l+m-1.0) / ((2.0*l-3)*(l-m)*(l+m)));
        legP[dQidx(0,l,m)] = c1 * x * legP[dQidx(0,l-1,m)] - c2 * legP[dQidx(0,l-2,m)];
    }
  }

}

void Grav3D::fillLegP(Real* legP, Real x){

  for ( int l = 0; l <= LMAX; l++ ){
    for ( int m = 0; m <= l; m++ ){
      legP[Qidx(0,l,m)] = 0.;
    }
  }

//Initial polynomial for recursion relations
  legP[Qidx(0,0,0)] = 1./sqrt(4.*M_PI);

//Diagonal 
  for( int m = 1; m <= LMAX; m++)
  {
    legP[Qidx(0,m,m)] = - sqrt( 1. + 1. / 2. / m ) * sqrt( 1. - x * x ) * legP[Qidx(0,m-1,m-1)];
  }

  for( int m = 0; m < LMAX; m++)
  {
    legP[Qidx(0,m+1,m)] = sqrt( 2. * m + 3. ) * x * legP[Qidx(0,m,m)];
  }

  for( int m = 0; m <= LMAX; m++){ 
    for( int l = m + 2; l <= LMAX; l++){
      Real c1 = sqrt( ((2.0*l+1)*(2.0*l-1)) / ((l+m)*(l-m)));
      Real c2 = sqrt( (2.0*l+1)*(l-m-1.0)*(l+m-1.0) / ((2.0*l-3)*(l-m)*(l+m)));
        legP[Qidx(0,l,m)] = c1 * x * legP[Qidx(0,l-1,m)] - c2 * legP[Qidx(0,l-2,m)];
    }
  }

}

__device__ int tidFake(int tid_x, int tid_y, int tid_z, int n_ghost, int *n){

  int tid = ( tid_z + n_ghost ) * n[0] * n[1] + ( tid_y + n_ghost ) * n[0] + ( tid_x + n_ghost );
  #ifdef POISSON_TEST
  int tid_z_fake = tid / ( n[0] * n[1]);
  int tid_y_fake = ( tid - tid_z_fake * n[0] * n[1] ) / n[0];
  int tid_x_fake = tid - tid_z_fake * n[0] * n[1] - tid_y_fake * n[0];
  if ( tid_x_fake < n_ghost || tid_y_fake < n_ghost || tid_z_fake < n_ghost || tid_z_fake > n[2] - n_ghost || tid_y_fake > n[1] - n_ghost || tid_x_fake > n[0] - n_ghost || tid >= n[0] * n[1] * n[2] || tid_z_fake != tid_z + n_ghost || tid_y_fake != tid_y + n_ghost || tid_x_fake != tid_x + n_ghost){
    printf("Something wrong in cell mapping.\n");
  }
  #endif
  return tid;
}

__global__ void QlmKernel(Real *rho, Real *center, Real *bounds, Real *dx, Real rmpole, int *n, int n_ghost, Real *partialReQ, Real *partialImQ){

  __shared__ Real ReQ[QTPB * (1 + LMAX ) * (2 + LMAX ) / 2];
  __shared__ Real ImQ[QTPB * (1 + LMAX ) * (2 + LMAX ) / 2];

  for ( int i = threadIdx.x * ( 1 + LMAX ) * ( 2 + LMAX ) / 2; i < ( threadIdx.x + 1 ) * ( 1 + LMAX ) * ( 2 + LMAX ) / 2; i++ ){
    ReQ[i] = 0.;
    ImQ[i] = 0.;
  }

  int nreal[3];
  for ( int i = 0; i < 3; i++ ) nreal[i] = n[i] - 2 * n_ghost;
  int nrealcells = nreal[0] * nreal[1] * nreal[2];

  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int tid_z = tid / ( nreal[0] * nreal[1] );
  int tid_y = ( tid - tid_z * nreal[0] * nreal[1] ) / nreal[0];
  int tid_x = tid - tid_z * nreal[0] * nreal[1] - tid_y * nreal[0];
  int cidx = threadIdx.x;
  int stride = blockDim.x * gridDim.x;

  Real r, phi, fac, pos[3], dev_legP[(1+LMAX)*(2+LMAX)/2];

  while ( tid < nrealcells ){

    tid_z = tid / ( nreal[0] * nreal[1] );
    tid_y = ( tid - tid_z * nreal[0] * nreal[1] ) / nreal[0];
    tid_x = tid - tid_z * nreal[0] * nreal[1] - tid_y * nreal[0];

    pos[0] = bounds[0] + dx[0] * ( tid_x + 0.5) - center[0];
    pos[1] = bounds[1] + dx[1] * ( tid_y + 0.5) - center[1];
    pos[2] = bounds[2] + dx[2] * ( tid_z + 0.5) - center[2];
    r = sqrt( pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2] );

    if ( r < rmpole ){
      phi = atan2(pos[1], pos[0]);

      fillLegP(dev_legP, pos[2] / r);

      for ( int l = 0; l <= LMAX; l++ ){
        fac = pow(r, l) * rho[tidFake(tid_x, tid_y, tid_z, n_ghost, n)];

        for ( int m = 0; m <= l; m++ ){
          ReQ[dQidx(cidx, l, m)] += dev_legP[dQidx(0,l,m)] * fac * cos(m * phi);
          ImQ[dQidx(cidx, l, m)] += dev_legP[dQidx(0,l,m)] * fac * sin(m * phi);
        }
      }
    }
    tid += stride;
  }

  __syncthreads();

  int i = blockDim.x / 2;
  while ( i > 0 ){
    if ( cidx < i){
      for ( int l = 0; l <= LMAX; l++ ){
        for ( int m = 0; m <= l; m++ ){
          ReQ[dQidx(cidx, l, m)] += ReQ[dQidx(cidx + i, l, m)];
          ImQ[dQidx(cidx, l, m)] += ImQ[dQidx(cidx + i, l, m)];
        }
      }
    }
    __syncthreads();
    i /= 2;
  }

  if ( cidx == 0 ){
    for ( int l = 0; l <= LMAX; l++ ){
      for ( int m = 0; m <= l; m++ ){
        partialReQ[dQidx(blockIdx.x, l, m)] = ReQ[dQidx(0, l, m)];
        partialImQ[dQidx(blockIdx.x, l, m)] = ImQ[dQidx(0, l, m)];
      }
    }
  }
}


__global__ void centerKernel(Real *rho, Real *bounds, Real *dx, int *n, int n_ghost, Real *partialCenter, Real *partialTotrhosq){

  __shared__ Real bCenter[3 * CENTERTPB];
  __shared__ Real bTotrhosq[CENTERTPB];

  for ( int i = threadIdx.x * 3; i < ( threadIdx.x + 1 ) * 3; i++ ) bCenter[i] = 0.;
  bTotrhosq[threadIdx.x] = 0.;

  int nreal[3], tid[3];
  for ( int i = 0; i < 3; i++ ) nreal[i] = n[i] - 2 * n_ghost;
  int nrealcells = nreal[0] * nreal[1] * nreal[2];
  int tid1d = threadIdx.x + blockIdx.x * blockDim.x;
  tid[2] = tid1d / ( nreal[0] * nreal[1] );
  tid[1] = ( tid1d - tid[2] * nreal[0] * nreal[1] ) / nreal[0];
  tid[0] = tid1d - tid[2] * nreal[0] * nreal[1] - tid[1] * nreal[0];

  Real x[3], rhosq;
  int fid;

  while ( tid1d < nrealcells ){
    tid[2] = tid1d / ( nreal[0] * nreal[1] );
    tid[1] = ( tid1d - tid[2] * nreal[0] * nreal[1] ) / nreal[0];
    tid[0] = tid1d - tid[2] * nreal[0] * nreal[1] - tid[1] * nreal[0];
    fid = tidFake(tid[0], tid[1], tid[2], n_ghost, n);

    rhosq = rho[fid] * rho[fid];

    bTotrhosq[threadIdx.x] += rhosq;
    for ( int i = 0; i < 3; i++ ){
      x[i] = bounds[i] + dx[i] * ( tid[i] + 0.5);
      bCenter[3 * threadIdx.x + i] += x[i] * rhosq;
    }

    tid1d += blockDim.x * gridDim.x;

  }

  __syncthreads();

  int i = blockDim.x / 2;
  while ( i > 0 ){
    if ( threadIdx.x < i ){
      for ( int ii = 0; ii < 3; ii++ ) bCenter[3 * threadIdx.x + ii] += bCenter[3 * ( threadIdx.x + i ) + ii];
      bTotrhosq[threadIdx.x] += bTotrhosq[threadIdx.x + i];
    }
    __syncthreads();
    i /= 2;
  }

  if ( threadIdx.x == 0 ){
    for ( int i = 0; i < 3; i++) partialCenter[3 * blockIdx.x + i] = bCenter[i];
    partialTotrhosq[blockIdx.x] = bTotrhosq[0];
  }

}

//TODO: rmpole should be the distance from the center of the expansion to the nearest boundary cell, not from the center of the domain to the nearest boundary cell
void Grid3D::setMoments(){

  Real dx[3], bounds[3];
  int n[3];
  dx[0] = H.dx; dx[1] = H.dy; dx[2] = H.dz;
  Real dV = dx[0] * dx[1] * dx[2];
  bounds[0] = H.xblocal; bounds[1] = H.yblocal; bounds[2] = H.zblocal;

  n[0] = H.nx; n[1] = H.ny; n[2] = H.nz;

  #ifdef DYNAMIC_GPU_ALLOC
  AllocateMemoryBoundaries_GPU(n[0], n[1], n[2];
  #endif

  chprintf("n copy: %i %i %i\n", n[0], n[1], n[2]);
//Find the center of the expansion according to Couch et al. 2013
  CudaSafeCall( cudaMemcpy( Grav.dev_rho   , C.density, n[0] * n[1] * n[2] * sizeof(Real), cudaMemcpyHostToDevice) );
  CudaSafeCall( cudaMemcpy( Grav.dev_bounds, bounds   , 3                  * sizeof(Real), cudaMemcpyHostToDevice) );
  CudaSafeCall( cudaMemcpy( Grav.dev_n     , n        , 3                  * sizeof(int ), cudaMemcpyHostToDevice) );
  CudaSafeCall( cudaMemcpy( Grav.dev_dx    , dx       , 3                  * sizeof(Real), cudaMemcpyHostToDevice) );

  centerKernel<<<Grav.centerBlocks,CENTERTPB>>>(Grav.dev_rho, Grav.dev_bounds, Grav.dev_dx, Grav.dev_n, H.n_ghost, Grav.dev_partialCenter, Grav.dev_partialTotrhosq);
  CudaCheckError();

  CudaSafeCall( cudaMemcpy(Grav.bufferCenter  , Grav.dev_partialCenter  , 3 * Grav.centerBlocks   * sizeof(Real), cudaMemcpyDeviceToHost) );
  CudaSafeCall( cudaMemcpy(Grav.bufferTotrhosq, Grav.dev_partialTotrhosq,     Grav.centerBlocks   * sizeof(Real), cudaMemcpyDeviceToHost) );

  Real totrhosq = 0.;
  for ( int i = 0; i < Grav.centerBlocks; i++ ){
    totrhosq += Grav.bufferTotrhosq[i];
  }
  #ifdef MPI_CHOLLA
  MPI_Allreduce(MPI_IN_PLACE, &totrhosq, 1, MPI_CHREAL, MPI_SUM, world);
  #endif

  for ( int i = 0; i < 3; i++ ) Grav.center[i] = 0.;
  for ( int i = 0; i < Grav.centerBlocks; i++ ){
     for ( int ii = 0; ii < 3; ii++ ) Grav.center[ii] += Grav.bufferCenter[3 * i + ii];
  }
  #ifdef MPI_CHOLLA
  MPI_Allreduce(MPI_IN_PLACE, Grav.center, 3, MPI_CHREAL, MPI_SUM, world);
  #endif

  for ( int i = 0; i < 3; i++ ) Grav.center[i] /= totrhosq;

  if ( H.n_step > 0) chprintf(" ");
  chprintf("Multipole center: %.10e, %.10e, %.10e\n", Grav.center[0], Grav.center[1], Grav.center[2]);

//Qnl
  cudaMemcpy( Grav.dev_rho, C.density, n[0] * n[1] * n[2]*sizeof(Real), cudaMemcpyHostToDevice);
  cudaMemcpy( Grav.dev_center, Grav.center, 3*sizeof(Real), cudaMemcpyHostToDevice);
  cudaMemcpy( Grav.dev_bounds, bounds, 3*sizeof(Real), cudaMemcpyHostToDevice);
  cudaMemcpy( Grav.dev_n, n, 3*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy( Grav.dev_dx, dx, 3*sizeof(Real), cudaMemcpyHostToDevice);

  QlmKernel<<<Grav.Qblocks,QTPB>>>(Grav.dev_rho, Grav.dev_center, Grav.dev_bounds, Grav.dev_dx, H.xdglobal / 2., Grav.dev_n, H.n_ghost, Grav.dev_partialReQ, Grav.dev_partialImQ);

  cudaMemcpy(Grav.bufferReQ, Grav.dev_partialReQ, sizeof(Real) * Grav.Qblocks * (1 + LMAX ) * (2 + LMAX ) / 2, cudaMemcpyDeviceToHost);
  cudaMemcpy(Grav.bufferImQ, Grav.dev_partialImQ, sizeof(Real) * Grav.Qblocks * (1 + LMAX ) * (2 + LMAX ) / 2, cudaMemcpyDeviceToHost);

  #ifdef DYNAMIC_GPU_ALLOC
  FreeMemoryBoundaries_GPU();
  #endif

  for ( int i = 0; i < ( 1 + LMAX ) * ( 2 + LMAX ) / 2; i++ ){
    Grav.ReQ[i] = 0.;
    Grav.ImQ[i] = 0.;
  }

  for ( int l = 0; l <= LMAX; l++ ){
    for ( int m = 0; m <= l; m++ ){
      for ( int b = 0; b < Grav.Qblocks; b++ ){

        Grav.ReQ[Grav.Qidx(0,l,m)] += Grav.bufferReQ[Grav.Qidx(b, l, m)];
        Grav.ImQ[Grav.Qidx(0,l,m)] += Grav.bufferImQ[Grav.Qidx(b, l, m)];

      }

      Grav.ReQ[Grav.Qidx(0,l,m)] *= dV;
      Grav.ImQ[Grav.Qidx(0,l,m)] *= dV;

    }
  }

  #ifdef MPI_CHOLLA
  MPI_Allreduce(MPI_IN_PLACE, Grav.ReQ, (1 + LMAX ) * (2 + LMAX ) / 2, MPI_CHREAL, MPI_SUM, world);
  MPI_Allreduce(MPI_IN_PLACE, Grav.ImQ, (1 + LMAX ) * (2 + LMAX ) / 2, MPI_CHREAL, MPI_SUM, world);
  #endif

  #ifdef POISSON_TEST
  int lmidx;
  for ( int l = 0; l <= LMAX; l++ ){
    for ( int m = 0; m <= l; m++ ){
      lmidx = Grav.Qidx(0,l,m);

      chprintf("ReQ[%i][%i]=%.20e\n", l, m, Grav.ReQ[lmidx]);
      chprintf("ImQ[%i][%i]=%.20e\n", l, m, Grav.ImQ[lmidx]);

    }
  }
  #endif

}

#endif
