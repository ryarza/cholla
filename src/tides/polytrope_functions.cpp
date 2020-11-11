#ifdef TIDES

#include "../global.h"
#include "../grid3D.h"
#include "../io.h"
#include <math.h>  
#include <vector>
//#include <cmath>
#include "../error_handling.h"

#ifdef MPI_CHOLLA
#include "mpi_routines.h"
#endif
     
using namespace std; 

class Real2 {
public:
  Real x;
  Real y;
  
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

Real Interpolate( int n, int rootIdx, Real xi, Real *xiVals, Real *thetaVals, Real *dthetaVals ){
//  if ( x <= xVals[0] )   return thetaVals[0];
//  if ( x >= xVals[n-1] ) return thetaVals[n-1];
  if ( xi < 0. ) chprintf("Error: radius requested < 0"    );
//    if ( r > P.Rstar   ) chprintf("Error: radius requested > Rstar");
  if ( xi > xiVals[rootIdx] ){
    Real val = thetaVals[rootIdx] + ( xi - xiVals[rootIdx] ) * dthetaVals[rootIdx];
    chprintf("This cell is after the last xi. Returning %.10e\n", val);
    return val;
  }
  if ( xi > xiVals[rootIdx+1]) chprintf("wtf\n");

// Find the closest index for which xVal is less than x;
  int indx = Binary_Search( n, xi, xiVals, 0, n-1 );
  if ( xi < xiVals[indx] ) printf(" interpolation error. req xi = %.10e, left xi = %.10e\n", xi, xiVals[indx]);
  if ( xi > xiVals[indx + 1] ) printf(" interpolation error. req xi = %.10e, right xi = %.10e\n", xi, xiVals[indx + 1]);
  Real x_0, x_1, dx, y_0, y_1;
  x_0 = xiVals[indx];
  x_1 = xiVals[indx+1];
  dx = x_1 - x_0;
  y_0 = thetaVals[indx    ];
  y_1 = thetaVals[indx + 1];
  return y_0 + ( y_1 - y_0 )/dx*(xi-x_0);
}

void Grid3D::Polytropic_Star( struct parameters &P ){

  chprintf(" Lane-Emden solver:\n");

  //Solve Lane–Emden equation for the polytrope 
//  int n_points = 500000000;
  int n_points = 10000000;
  Real *xi_vals      = new Real[n_points];
  Real *theta_vals   = new Real[n_points];
  Real *theta_deriv  = new Real[n_points];
  Real *R_vals       = new Real[n_points];
  Real *density_vals = new Real[n_points];
  
  Real xi_min, xi_max, dxi;

  xi_min = 1.e-10;

  if ( P.polyN == 0. ){
    xi_max = 2.46;
  }
  else if ( P.polyN == 1. ){
    xi_max = 3.15;
  }
  else if ( P.polyN == 1.5 ){
    xi_max = 3.65376;
  }
  else{
    xi_max = 7.;
  }

  dxi = xi_max / ( n_points - 1. );
  xi_vals[0] = 0.;
  xi_vals[1] = xi_min;
  for ( int i = 2; i < n_points; i++){
    xi_vals[i] = i * dxi;
  }

// Vector that will hold the solutions
  vector<Real2> poly_coords;

// The first elements of this vector will be the known boundary conditions
  poly_coords.push_back( Real2(1., 0.) );
  theta_vals[0] = 1.;
  theta_deriv[0] = 0.;

// We can't start the integration from the previous point because there'll be division by zero. Instead, integrate the first point from the known Taylor series solution to the equation
  Real thetaTaylor, dthetaTaylor;

  thetaTaylor = 1. - (1./6.) * xi_min * xi_min + P.polyN * pow(xi_min, 4.) / 120. - P.polyN * ( 8. * P.polyN - 5.) * pow(xi_min, 6.) / 15120.;
  dthetaTaylor = - xi_min / 3. + P.polyN * pow(xi_min, 3.) / 30. - pow(xi_min, 5.) * P.polyN * ( -5. + 8 * P.polyN ) / 2520.;

  poly_coords.push_back( Real2( thetaTaylor, dthetaTaylor ) );

//  chprintf("Taylor series: %.10e, %.10e", thetaTaylor, dthetaTaylor);

  theta_vals[1] = poly_coords[1].x;
  theta_deriv[1] = poly_coords[1].y;

  //Solve the polytrope equation using the RK4 module
  Real2 coords, coords_new;
  Real dx, x;
  for (int i = 1; i < n_points - 1; i++){
    dx = xi_vals[i+1] - xi_vals[i];
    x = xi_vals[i];
    coords = poly_coords[i];
    coords_new = RK4_Iteration( Derivative, x, coords, dx, P.polyN );
    poly_coords.push_back( coords_new );
    theta_vals[i+1] = poly_coords[i+1].x;
    theta_deriv[i+1] = poly_coords[i+1].y;
  }


  //Find the root of theta as a function of xi;
  int root_indx = 0;
  for (int i=0; i<n_points; i++ ){
    root_indx = i;
    if ( theta_vals[i]*theta_vals[i+1] < 0) break;
  }

/*
  for ( int i = 0; i < root_indx + 1; i++){
    chprintf("xi = %.10e, theta = %.10e, dtheta = %.10e\n", xi_vals[i], theta_vals[i], theta_deriv[i]);
  }
*/

//  Linear interpolation estimate of the root
  Real xi_root = ( xi_vals[root_indx + 1] * theta_vals[root_indx] - xi_vals[root_indx] * theta_vals[root_indx + 1] ) / ( theta_vals[root_indx] - theta_vals[root_indx + 1]  );
  chprintf( "  Root at xi = %.5e. Theta before and after: %.5e %.5e\n", xi_root, theta_vals[root_indx], theta_vals[root_indx+1] );
  
//  Linear extrapolation estimate of the derivative evaluated at the root
  Real theta_deriv_root = xi_vals[root_indx + 1] * theta_vals[root_indx] * ( theta_deriv[root_indx - 1] - theta_deriv[root_indx] );
  theta_deriv_root += xi_vals[root_indx - 1] * theta_deriv[root_indx] * ( theta_vals[root_indx] - theta_vals[root_indx + 1] );
  theta_deriv_root += xi_vals[root_indx] * ( theta_vals[root_indx + 1] * theta_deriv[root_indx] - theta_vals[root_indx] * theta_deriv[root_indx - 1] );
  theta_deriv_root /= ( xi_vals[root_indx - 1] - xi_vals[root_indx] ) * ( theta_vals[root_indx] - theta_vals[root_indx  + 1] );
  
  chprintf( "  d(theta)/d(xi) at the root: %.5e\n", theta_deriv_root );
  
  //Convert to physical values
  Real dens_avrg = ( 3 * P.Mstar ) / ( 4 * M_PI * pow( P.Rstar, 3) );
  Real beta = - ( xi_root / 3 / theta_deriv_root );
  Real dens_central = beta * dens_avrg;
  Real pressure_central = G_CGS * P.Mstar * P.Mstar / pow( P.Rstar, 4 ) / ( 4 * M_PI *( P.polyN+1 ) * theta_deriv_root * theta_deriv_root );
  Real K = pressure_central * pow( dens_central, -(P.polyN+1)/P.polyN );
  Real alpha = sqrt( (P.polyN + 1) * K / ( 4 * M_PI * G_CGS ) ) * pow( dens_central, (1.-P.polyN)/(2*P.polyN) );
  
  chprintf( "  rho_c / rho_av: %.5e g/cm^3\n", dens_central / dens_avrg );
  chprintf( "  p_c           : %.5e erg/cm^3\n", pressure_central );
  Real cs_center = sqrt(  pressure_central / dens_central * P.gamma );
  chprintf( "  t_cross       : %.5e s\n", P.Rstar / cs_center); 
  // chprintf( " K: %f \n", K );
  // chprintf( " alpha: %f \n", alpha );
  for ( int i=0; i<n_points; i++){
    R_vals[i] = alpha * xi_vals[i];
    density_vals[i] = dens_central * pow( theta_vals[i], P.polyN );
    if ( theta_vals[i] < 0 ) density_vals[i] = 0;
  }
  
  int i, j, k, id;
  Real x_pos, y_pos, z_pos, r, center_x, center_y, center_z;
  Real density, pressure, energy;
  #ifdef DE
  Real gasEnergy;
  #endif

  center_x = 0.;
  center_y = 0.;
  center_z = 0.;

  Real vx, vy, vz, v2;
  vx = 0.;
  vy = 0.;
  vz = 0.;
  v2 = vx*vx + vy*vy + vz*vz;

//We must assign the mean value of the polytropic solution to the cell. To do so, we divide the cell into subcells, evaluate the solution at every subcell, then assign to the cell the average of the subcells. We sample values using a "np.linspace" from x-dx/2 to x+dx/2 but not including the endpoints, so that we use only values inside the cell.

// Number of subcells in each dimension. The total number of subcells per cell is this number cubed
  int nSubcells = 5;
  Real totSubcells = pow(nSubcells, 3.);

//These variables hold the value of x-dx/2 (lo) and x+dx/2 (hi)
  Real xSubLo, xSubHi, ySubLo, ySubHi, zSubLo, zSubHi, dxSub, dySub, dzSub;

//  Subcell properties
  Real thetaSubcell, rhoSubcell, pressureSubcell, energySubcell;
  #ifdef DE
  Real gasEnergySubcell;
  #endif

  // set the initial values of the conserved variables
  for (k=H.n_ghost; k<H.nz-H.n_ghost; k++) {
    for (j=H.n_ghost; j<H.ny-H.n_ghost; j++) {
      for (i=H.n_ghost; i<H.nx-H.n_ghost; i++) {
        id = i + j*H.nx + k*H.nx*H.ny;

        // // get the centered cell positions at (i,j,k)
        Get_Position(i, j, k, &x_pos, &y_pos, &z_pos);
        r = sqrt( (x_pos-center_x)*(x_pos-center_x) + (y_pos-center_y)*(y_pos-center_y) + (z_pos-center_z)*(z_pos-center_z) );

        if ( r < P.Rstar + sqrt( H.dx * H.dx + H.dy * H.dy + H.dz * H.dz ) ){
//        This variable will hold the value at the cell center
          density = 0.;
          pressure = 0.;
          energy = 0.;
          #ifdef DE
          gasEnergy = 0.;
          #endif

//        Get the lower and upper limits of the linspace in every direction
          xSubLo = x_pos - H.dx / 2.;
          xSubHi = x_pos + H.dx / 2.;
          ySubLo = y_pos - H.dy / 2.;
          ySubHi = y_pos + H.dy / 2.;
          zSubLo = z_pos - H.dz / 2.;
          zSubHi = z_pos + H.dz / 2.;
          dxSub  = ( xSubHi - xSubLo ) / nSubcells;
          dySub  = ( ySubHi - ySubLo ) / nSubcells;
          dzSub  = ( zSubHi - zSubLo ) / nSubcells;

          for (int kk = 0; kk < nSubcells; kk++){
            z_pos = zSubLo + ( kk + 0.5 ) * dzSub;
            for (int jj = 0; jj < nSubcells; jj++){
              y_pos = ySubLo + ( jj + 0.5 ) * dySub;
              for (int ii = 0; ii < nSubcells; ii++){
                x_pos = xSubLo + ( ii + 0.5 ) * dxSub;

                r = sqrt( (x_pos-center_x)*(x_pos-center_x) + (y_pos-center_y)*(y_pos-center_y) + (z_pos-center_z)*(z_pos-center_z) );

                if ( r < P.Rstar ){
//                We interpolate linearly on theta. This interpolation is not equivalent to interpolating rho, since rho is a nonlinear function of theta!
                  thetaSubcell    = Interpolate( n_points, root_indx, r / alpha, xi_vals, theta_vals, theta_deriv );
                  rhoSubcell      = dens_central * pow( thetaSubcell, P.polyN );
                  pressureSubcell = K * pow(rhoSubcell, 1. + 1. / P.polyN );
                  energySubcell   = pressureSubcell / ( P.gamma - 1. ) + 0.5 * rhoSubcell * v2;
                }
                else{
                  rhoSubcell      = P.rhoAmb;
                  pressureSubcell = P.pAmb;
                  energySubcell   = pressure / ( P.gamma - 1. ) + 0.5 * density * v2;
                }

                density  += rhoSubcell / totSubcells;
                pressure += pressureSubcell / totSubcells;
                energy   += energySubcell / totSubcells;

                #ifdef DE
                gasEnergySubcell = pressureSubcell / ( P.gamma - 1. );
                gasEnergy += gasEnergySubcell / totSubcells;
                #endif

              }
            }
          }
        }

        else{
          density  = P.rhoAmb;
          pressure = P.pAmb;
          energy = pressure / ( P.gamma - 1. ) + 0.5 * density * v2;
          #ifdef DE
          gasEnergy = pressure / ( P.gamma - 1. );
          #endif
        }


        C.density[id] = density;
        C.momentum_x[id] = density*vx;
        C.momentum_y[id] = density*vy;
        C.momentum_z[id] = density*vz;
        C.Energy[id] = energy;

        #ifdef DE
        C.GasEnergy[id] = gasEnergy;
        #endif

      }
    }
  }
  
  //Free the polytrope data
  delete[] xi_vals;
  delete[] theta_vals;
  delete[] theta_deriv;
  delete[] R_vals;
  delete[] density_vals;
  poly_coords.clear();
  
}

void Grid3D::damp(){
  
  Real dens, vx, vy, vz, v2, E, U;
  Real dens_0, vx_0, vy_0, vz_0;
  
  Real dt = H.dt;
  Real t = H.t;

// Ohlmann+ 2018
/*
  Real tau;
  Real tau1 = tdyn / 10.;
  Real tau2 = tdyn      ;
*/

  Real relaxRate;
  if ( t <= S.tRelax ){
//  This relaxRate will apply only during the relaxation procedure to all cells
    relaxRate = ( t / S.tRelax ) * ( 1. - S.relaxRate0 ) + S.relaxRate0;
  }
  else{
//  This relaxRate will apply only after the relaxation procedure, and only to cells with low densities
    relaxRate = S.relaxRateBkgnd;
  }

  int i, j, k, id;
  for ( k = 0; k < Grav.nz_local; k++ ) {
    for ( j = 0; j < Grav.ny_local; j++ ) {
      for ( i = 0; i < Grav.nx_local; i++ ) {
        id = (i+H.n_ghost) + (j+H.n_ghost)*H.nx + (k+H.n_ghost)*H.nx*H.ny;

        dens = C.density[id];
        vx   = C.momentum_x[id] / dens;
        vy   = C.momentum_y[id] / dens;
        vz   = C.momentum_z[id] / dens;
        E    = C.Energy[id];
        v2 = vx*vx + vy*vy + vz*vz;
        U = E - 0.5*dens*v2;

//      Ohlmann+ 2018 relaxation
/*
        dens_0 = G.C.density_0[id];
        vx_0   = G.C.momentum_x_0[id] / dens_0;
        vy_0   = G.C.momentum_y_0[id] / dens_0;
        vz_0   = G.C.momentum_z_0[id] / dens_0;

        if (t < 2. * tdyn) {
          tau = tau1;
        }
        else {
          tau = tau1 * pow( tau2 / tau1, (t - 2. * tdyn) / ( 3. * tdyn ) );
        }

        vx -= vx_0 * dt / tau;
        vy -= vy_0 * dt / tau;
        vz -= vz_0 * dt / tau;
*/

//      Guillochon+ 2013 relaxation
//      The first criterion will apply to all cells when t < tRelax, and only to low density cells when t > tRelax
        if ( t < S.tRelax || dens < 1.e1 * 1.e-10 ){
          vx *= relaxRate;
          vy *= relaxRate;
          vy *= relaxRate;
        }

        v2 = vx*vx + vy*vy + vz*vz;
//        v = sqrt( v2 );

//      Compute the energy with the updated kinetic energy
        E = U + 0.5*dens*v2;
        
//      Save the updated values
        C.momentum_x[id] = dens*vx;
        C.momentum_y[id] = dens*vy;
        C.momentum_z[id] = dens*vz;
        C.Energy[id] = E;

      }
    }
  }
  
  //Get the global max speed
  /*
  Real max_speed_global = max_speed;
  #ifdef MPI_CHOLLA
  max_speed_global = ReduceRealMax( max_speed );
  #endif
  */
  
}


void Grid3D::Polytropic_Star_Relaxation(  struct parameters &P  ){
  
  chprintf("Polytropic star iterative relaxation\n");
  int nfile = 0;
  Real dti;
  
  int n_step = 0;
  
  chprintf(" Star dynamical time: %.5e\n", S.tdynStar);

  WriteData(*this, P, P.nfile);
  P.nfile++;

  while (H.t < S.tRelax ){
  
    chprintf(" Relaxation n_step: %d\n", n_step + 0 );

    S.update(H.t, H.dt);
    updateCOM();
    // calculate the timestep
    set_dt(dti);
    
    // Advance the grid by one timestep
    dti = Update_Hydro_Grid();
    
    // update the simulation time ( t += dt )
    Update_Time();
    
    // add one to the timestep count
    n_step++;
    
    #ifdef GRAVITY
    //Compute Gravitational potential for next step
    Compute_Gravitational_Potential( &P);
    #endif
    
    //Include the damping terms in momentum and energy
    damp();

//  TODO: Change from number of steps to time so that it's consistent with the rest of the code
    // Output
    if (n_step % int(P.outstep) == 0){
      WriteData(*this, P, P.nfile);
      P.nfile++;
    }

    // set boundary conditions for next time step 
    Set_Boundary_Conditions_Grid(P);
    
    chprintf("n_step: %d  sim time: %10.7f  sim timestep: %7.4e  \n\n", n_step, H.t, H.dt);

  }

  S.relaxed = 1;
  //Restart the simulation time
  H.t = 0;
  Grav.INITIAL = true;
  
  chprintf("Finished relaxation on step %d\n", n_step);
  
}

#endif
