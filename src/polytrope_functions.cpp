#include "global.h"
#include "grid3D.h"
#include "io.h"
#include <math.h>  
#include <vector> 
     
using namespace std; 

class Real2 {
public:
  Real x = 0;
  Real y = 0;
  
  //Cosntructor
  Real2( Real x0=0, Real y0=0 ): x(x0), y(y0) {}
  
  //Sum operator
  Real2 operator+( const Real2 &v ){
    return Real2( x+v.x, y+v.y );
  }
  Real2 operator+( const Real &v ){
    return Real2( x+v, y+v );
  }
  Real2 operator*( const Real &v ){
    return Real2( x*v, y*v );
  }
};

Real2 Derivative( Real x, Real2 coords, Real n  ){
  Real theta, y;
  theta = coords.x;
  y = coords.y;

  Real2 deriv;
  deriv.x = y; 
  deriv.y = -1*pow( theta, n ) - 2*y/x;

  return deriv;  
}

Real2 RK4_Iteration( Real2 (*f)( Real x, Real2 coords, Real n ), Real t, Real2 x, Real dt, Real n ){
  Real2 k1, k2, k3, k4;
  k1 = f( t, x, n ) * dt;
  k2 = f( t + 0.5*dt, x+k1*0.5, n ) * dt;
  k3 = f( t + 0.5*dt, x+k2*0.5, n ) * dt;
  k4 = f( t + dt, x+k3, n ) * dt;
  return x + ( k1 + k2*2.0 + k3*2.0 + k4)*(1./6);  
}


int Binary_Search( int N, Real val, Real *data, int indx_l, int indx_r ){
  int n, indx;
  n = indx_r - indx_l;
  indx = indx_l + n/2;
  if ( val >= data[N-1] ) return indx_r;
  if ( val <= data[0]   ) return indx_l;
  if ( indx_r == indx_l + 1 ) return indx_l;
  if ( data[indx] <= val ) indx_l = indx;
  else indx_r = indx;
  return Binary_Search( N, val, data, indx_l, indx_r );
}

Real Interpolate( int n, Real x, Real *R_vals, Real *density_vals ){
  if ( x <= R_vals[0] )   return density_vals[0];
  if ( x >= R_vals[n-1] ) return density_vals[n-1];
  
  //Find the closest index which R_value is less than x;
  int indx = Binary_Search( n, x, R_vals, 0, n-1 );
  if ( x < R_vals[indx] ) printf(" interpolation Error\n" );
  if ( x > R_vals[indx+1] ) printf(" interpolation Error\n" );
  Real x_0, x_1, dx, y_0, y_1;
  x_0 = R_vals[indx];
  x_1 = R_vals[indx+1];
  dx = x_1 - x_0;
  y_0 = density_vals[indx];
  y_1 = density_vals[indx+1];
  return y_0 + ( y_1 - y_0 )/dx*(x-x_0);
}

void Grid3D::Polytropic_Star( struct parameters P ){
  
  Real n_poly = 3;
  
  chprintf( " Initializing Polytropic Star \n");
  chprintf( " Polytropic Index: %.1f  \n", n_poly );
  
  
  //Solve Lane–Emden equation for the polytrope 
  int n_points = 100000;
  Real *psi_vals =     new Real[n_points];
  Real *theta_vals =   new Real[n_points];
  Real *R_vals =       new Real[n_points];
  Real *density_vals = new Real[n_points];
  
  
  Real psi_min, psi_max, dpsi;
  psi_min = 1e-11;
  psi_max = 10;
  dpsi = ( psi_max - psi_min ) / n_points;
  for (int i=0; i<n_points; i++){
    psi_vals[i] = psi_min + (i+0.5)*dpsi;
  }
  vector<Real2> poly_coords;
  //Set initial condition for the polytrope solution
  poly_coords.push_back( Real2(1, 0) );
  theta_vals[0] = poly_coords[0].x;
  
  //Solve the polytrope equation usin the RK4 module
  Real2 coords, coords_new;
  Real dx, x;  
  for (int i=0; i<n_points-1; i++){
    dx = psi_vals[i+1] - psi_vals[i];
    x = psi_vals[i];
    coords = poly_coords[i];
    coords_new = RK4_Iteration( Derivative, x, coords, dx, n_poly );
    poly_coords.push_back( coords_new );
    theta_vals[i+1] = poly_coords[i+1].x;
  }
  
  //Convert tho physical values
  for ( int i=0; i<n_points; i++){
    R_vals[i] = psi_vals[i];
    density_vals[i] = theta_vals[i];
  }
  
  int i, j, k, id;
  Real x_pos, y_pos, z_pos, r, center_x, center_y, center_z;
  Real density, pressure,  energy, dens_min;
  Real vx, vy, vz, v2;
  center_x = 10;
  center_y = 10;
  center_z = 10;
  vx = 0;
  vy = 0;
  vz = 0;
  dens_min = 1e-5;
  // set the initial values of the conserved variables
  for (k=H.n_ghost; k<H.nz-H.n_ghost; k++) {
    for (j=H.n_ghost; j<H.ny-H.n_ghost; j++) {
      for (i=H.n_ghost; i<H.nx-H.n_ghost; i++) {
        id = i + j*H.nx + k*H.nx*H.ny;

        // // get the centered cell positions at (i,j,k)
        Get_Position(i, j, k, &x_pos, &y_pos, &z_pos);
        // density = 0.0005;
        pressure = 0.0005;

        r = sqrt( (x_pos-center_x)*(x_pos-center_x) + (y_pos-center_y)*(y_pos-center_y) + (z_pos-center_z)*(z_pos-center_z) );
        density = Interpolate( n_points, r, R_vals, density_vals );
        density = fmax( density, dens_min );
        
        v2 = vx*vx + vy*vy + vz*vz;
        energy = pressure/(gama-1) + 0.5*density*v2;
        C.density[id] = density;
        C.momentum_x[id] = density*vx;
        C.momentum_y[id] = density*vy;
        C.momentum_z[id] = density*vz;
        C.Energy[id] = energy;

        #ifdef DE
        C.GasEnergy[id] = pressure/(gama-1);
        #endif
      }
    }
  }
  
  //Free the Polytrope data
  delete[] psi_vals;
  delete[] theta_vals;
  delete[] R_vals;
  delete[] density_vals;
  poly_coords.clear();
  
}