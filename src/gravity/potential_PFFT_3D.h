#ifdef GRAVITY
#ifdef PFFT

#ifndef POTENTIAL_PFFT_3D_H
#define POTENTIAL_PFFT_3D_H

#include "../global.h"
#include <stdlib.h>
#include <cmath>
#include <time.h>

#include <pfft.h>

class Domain_PFFT_3D{
  
  public:
  ptrdiff_t nx_total;
  int ny_total;
  ptrdiff_t nz_total;
  
  ptrdiff_t nx_local;
  ptrdiff_t ny_local;
  ptrdiff_t nz_local;
  
  ptrdiff_t nx_local_start;
  ptrdiff_t ny_local_start;
  ptrdiff_t nz_local_start;

  ptrdiff_t nx_local_cholla;
  ptrdiff_t ny_local_cholla;
  ptrdiff_t nz_local_cholla;

  ptrdiff_t nx_local_start_cholla;
  ptrdiff_t ny_local_start_cholla;
  ptrdiff_t nz_local_start_cholla;
  
  Real domlen_x_cholla;
  Real domlen_y_cholla;
  Real domlen_z_cholla;
  
  Real xblocal_cholla;
  Real yblocal_cholla;
  Real zblocal_cholla;
  



  Real Lbox_x;
  Real Lbox_y;
  Real Lbox_z;
  
  Real dx;
  Real dy;
  Real dz;
  
  bool INITIALIZED;
  
  int procID_pfft;
  int nproc_pfft;
  MPI_Comm comm_pfft;
  
  ptrdiff_t alloc_local_fwd;
  ptrdiff_t alloc_local_bwd;
  int nprocs_grid_pfft[3];
  int pcoords_pfft[3];
  int poffset_pfft[3];
  ptrdiff_t n_pfft[3];
  ptrdiff_t local_n_in_pfft[3], local_in_start_pfft[3];
  ptrdiff_t local_n_out_pfft[3], local_out_start_pfft[3];
  ptrdiff_t local_n_transform_fwd_pfft[3], local_transform_fwd_start_pfft[3];
  ptrdiff_t local_n_transform_bwd_pfft[3], local_transform_bwd_start_pfft[3];
  

  Domain_PFFT_3D( void );
  void Initialize( struct parameters *P );
  
};

class Potential_PFFT_3D{
  
  public:
  
  Real Lbox_x;
  Real Lbox_y;
  Real Lbox_z;

  grav_int_t nx_total;
  grav_int_t ny_total;
  grav_int_t nz_total;

  int nx_local;
  int ny_local;
  int nz_local;

  Real dx;
  Real dy;
  Real dz;
  grav_int_t n_cells_total;
  grav_int_t n_cells_local;

  int procID_pfft;
  int nproc_pfft;
  MPI_Comm comm_pfft;

  int nprocs_grid_pfft[3];
  int pcoords_pfft[3];
  int poffset_pfft[3];
  ptrdiff_t n_pfft[3];
  ptrdiff_t alloc_local;
  ptrdiff_t local_ni_pfft[3], local_i_start_pfft[3];
  ptrdiff_t local_no_pfft[3], local_o_start_pfft[3];
  ptrdiff_t local_ntrans_pfft[3], local_trans_start_pfft[3];

  pfft_plan plan_fwd;
  pfft_plan plan_bwd;

  int index_0;

  Real xMin;
  Real yMin;
  Real zMin;
  
  Real *input_density;
  Real *output_potential;

  struct Fields
  {

    pfft_complex *transform;
    double *input;
    double *output;
    // Complex_fftw *transform;
    double *G;


  } F;

  Potential_PFFT_3D( void );

  void Initialize( Real Lx, Real Ly, Real Lz, Real x_min, Real y_min, Real z_min, int nx, int ny, int nz, int nx_real, int ny_real, int nz_real, Real dx, Real dy, Real dz, struct Header &H );
  
  void AllocateMemory_CPU( void );
  void Reset( void );
  
  void Copy_Input( Real *input_density, Real Grav_Constant, Real dens_avrg, Real current_a );
  void Copy_Output( Real *output_potential );
  void Get_K_for_Green_function( void );
  void Apply_G_Funtion( void );
  void Apply_K2_Funtion( void );
  void Get_Index_Global(int i, int j, int k, int *i_global, int *j_global, int *k_global);
  Real Get_Potential( Real *input_density,  Real *output_potential, Real Grav_Constant, Real dens_avrg, Real current_a );
  


};



#endif //POTENTIAL_PFFT_3D_H
#endif //PFFT
#endif //GRAVITY
