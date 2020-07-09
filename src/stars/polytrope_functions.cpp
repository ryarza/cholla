#ifdef STARS

#include "../global.h"
#include "../grid3D.h"
#include "../io.h"
#include <math.h>  
#include <vector> 

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

void Grid3D::Polytropic_Star( struct parameters &P ){
  Star.Mstar = P.Mstar;
  Star.Rstar = P.Rstar;
  Star.polyN = P.polyN;

  Real n_poly = P.polyN;
  
  chprintf("Polytrope mass: %.3f \n", P.Mstar);
  chprintf("Polytrope radius: %.3f \n", P.Rstar);
  chprintf( " Initializing Polytropic Star \n");
  chprintf( " Polytropic Index: %.1f  \n", n_poly );
  
  
  //Solve Lane–Emden equation for the polytrope 
  int n_points = 100000000;
  Real *psi_vals =     new Real[n_points];
  Real *theta_vals =   new Real[n_points];
  Real *theta_deriv =  new Real[n_points];
  Real *R_vals =       new Real[n_points];
  Real *density_vals = new Real[n_points];
  
  
  Real psi_min, psi_max, dpsi;
  psi_min = 1e-15;
  psi_max = 3.65374;
  dpsi = ( psi_max - psi_min ) / n_points;
  for (int i=0; i<n_points; i++){
    psi_vals[i] = psi_min + (i+0.5)*dpsi;
  }
  vector<Real2> poly_coords;
  //Set initial condition for the polytrope solution
  poly_coords.push_back( Real2(1, 0) );
  theta_vals[0] = poly_coords[0].x;
  theta_deriv[0] = poly_coords[0].y;
  
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
    theta_deriv[i+1] = poly_coords[i+1].y;
  }
  
  //Get the derivative of Theta with respect to Psi
  // Real theta_l, theta_r, d_psi, dtheta_dpsi;
  // for ( int i=0; i<n_points; i++){
  //   dx = psi_vals[i+1] - psi_vals[i];
  //   if ( i == 0 ){
  //     theta_l = theta_vals[i];
  //     theta_r = theta_vals[i+1];
  //     d_psi = psi_vals[i+1] - psi_vals[i];
  //   }
  //   else if ( i == n_points-1 ){
  //     theta_l = theta_vals[i-1];
  //     theta_r = theta_vals[i];
  //     d_psi = psi_vals[i] - psi_vals[i-1];
  //   }
  //   else{
  //     theta_l = theta_vals[i-1];
  //     theta_r = theta_vals[i+1];
  //     d_psi = psi_vals[i+1] - psi_vals[i-1];
  //   }
  //   dtheta_dpsi = ( theta_r - theta_l ) / d_psi;
  //   theta_deriv[i] = dtheta_dpsi;  
  // }
  
  //Find the root of theta as a function of psi;
  int root_indx = 0;
  Real psi_root ;
  for (int i=0; i<n_points; i++ ){
    root_indx = i;
    if ( theta_vals[i]*theta_vals[i+1] < 0) break;
  }
  psi_root = psi_vals[root_indx];
  chprintf( " Theta Root:  %f   ->  theta values: %f %f \n", psi_root, theta_vals[root_indx], theta_vals[root_indx+1] );
  
  //Get the derivative of theta evaluated at the root of theta
  Real theta_deriv_root = theta_deriv[root_indx];
  
  chprintf( " Theta Deriv at Root:  %f   \n", theta_deriv_root );
    
  
  //Convert to physical values
  Real dens_avrg = ( 3 * P.Mstar ) / ( 4 * M_PI * pow( P.Rstar, 3) );
  Real beta = - ( psi_root / 3 / theta_deriv_root );
  Real dens_central = beta * dens_avrg;
  Real pressure_central = G_CGS * P.Mstar * P.Mstar / pow( P.Rstar, 4 ) / ( 4 * M_PI *( n_poly+1 ) * theta_deriv_root * theta_deriv_root );
  Real K = pressure_central * pow( dens_central, -(n_poly+1)/n_poly );
  Real alpha = sqrt( (n_poly + 1) * K / ( 4 * M_PI * G_CGS ) ) * pow( dens_central, (1-n_poly)/(2*n_poly) );
  
  chprintf( " dens_central / dens_avrg: %f \n", dens_central / dens_avrg );
  chprintf( " pressure_central: %e \n", pressure_central );
  Real cs_center = sqrt(  pressure_central / dens_central * P.gamma );
  chprintf ( " T cross: %f s\n ", P.Rstar / cs_center); 
  // chprintf( " K: %f \n", K );
  // chprintf( " alpha: %f \n", alpha );
  for ( int i=0; i<n_points; i++){
    R_vals[i] = alpha * psi_vals[i];
    density_vals[i] = dens_central * pow( theta_vals[i], n_poly );
    if ( theta_vals[i] < 0 ) density_vals[i] = 0;
  }
  
  int i, j, k, id;
  Real x_pos, y_pos, z_pos, r, center_x, center_y, center_z;
  Real density, pressure, energy;
	#ifdef DE
	Real gasEnergy;
	#endif
  Real vx, vy, vz, v2;
  Real bkgndRho = 1.e-17;
  center_x = 0.;
  center_y = 0.;
  center_z = 0.;
  vx = 0.;
  vy = 0.;
  vz = 0.;
	v2 = vx*vx + vy*vy + vz*vz;

//Mean dimensionless molecular weight
  Real mu = 1.;

//We must assign the mean value of the polytropic solution to the cell. To do so, we divide the cell into subcells, evaluate the solution at every subcell, then assign to the cell the average of the subcells. We sample values using a "np.linspace" from x-dx/2 to x+dx/2 but not including the endpoints, so that we use only values inside the cell.

// Number of subcells in each dimension. The total number of subcells per cell is this number cubed
	int nSubcells = 5;
	Real totSubcells = pow(nSubcells, 3.);

//These variables hold the value of x-dx/2 (lo) and x+dx/2 (hi)
  Real xSubLo, xSubHi, ySubLo, ySubHi, zSubLo, zSubHi, dxSub, dySub, dzSub;

//	Subcell properties
	Real rhoSubcell, pressureSubcell, energySubcell;
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

				if ( r < P.Rstar ){
//				This variable will hold the value at the cell center
					density = 0.;
					pressure = 0.;
					energy = 0.;
					#ifdef DE
					gasEnergy = 0.;
					#endif

//				Get the lower and upper limits of the linspace in every direction
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

								rhoSubcell      = Interpolate( n_points, r, R_vals, density_vals );
								pressureSubcell = K * pow(rhoSubcell, 1. + 1. / P.polyN );
								energySubcell   = pressureSubcell / ( P.gamma - 1. ) + 0.5 * rhoSubcell * v2;

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
					density  = bkgndRho;
					pressure = 1.e0 * bkgndRho * KB / mu / MP;
					energy = pressure / ( P.gamma - 1. ) + 0.5 * density * v2;
					gasEnergy = pressure / ( P.gamma - 1. );
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
  
  //Free the Polytrope data
  delete[] psi_vals;
  delete[] theta_vals;
  delete[] theta_deriv;
  delete[] R_vals;
  delete[] density_vals;
  poly_coords.clear();
  
}

void Apply_Damping_Step( Grid3D &G, Real tdyn){
  
  Real max_speed = 0;
  
  Real dens, vx, vy, vz, v2, v, E, U;
  Real dens_0, vx_0, vy_0, vz_0;
  
  Real dt = G.H.dt;
  Real t = G.H.t;

  Real relaxRate = ( t / 6.278 / tdyn ) * ( 1. - 0.9 ) + 0.9;
  // chprintf("Relax rate: %.5e", relaxRate);

  int i, j, k, id;
  for (k=0; k<G.Grav.nz_local; k++) {
    for (j=0; j<G.Grav.ny_local; j++) {
      for (i=0; i<G.Grav.nx_local; i++) {
        id = (i+G.H.n_ghost) + (j+G.H.n_ghost)*G.H.nx + (k+G.H.n_ghost)*G.H.nx*G.H.ny;
        
        dens_0 = G.C.density_0[id];
        vx_0   = G.C.momentum_x_0[id] / dens_0;
        vy_0   = G.C.momentum_y_0[id] / dens_0;
        vz_0   = G.C.momentum_z_0[id] / dens_0;
        
        dens = G.C.density[id];
        vx   = G.C.momentum_x[id] / dens;
        vy   = G.C.momentum_y[id] / dens;
        vz   = G.C.momentum_z[id] / dens;
        E    = G.C.Energy[id];
        v2 = vx*vx + vy*vy + vz*vz;
        U = E - 0.5*dens*v2;

        //Ohlmann+2018 relaxation scheme.
      	/*
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


      	// James relaxation procedure
      	vx *= relaxRate;
      	vy *= relaxRate;
      	vy *= relaxRate;

        v2 = vx*vx + vy*vy + vz*vz;
        v = sqrt( v2 );
        if ( v > max_speed && dens > 1.e-0) max_speed = v;
        
        //Compute the energy with the updated kinetic energy
        E = U + 0.5*dens*v2;
        
        //Save the updated values
        G.C.momentum_x[id] = dens*vx; 
        G.C.momentum_y[id] = dens*vy; 
        G.C.momentum_z[id] = dens*vz;
        G.C.Energy[id] = E; 
        
      }
    }
  }
  
  //Get the global max speed
  Real max_speed_global = max_speed;
  #ifdef MPI_CHOLLA
  max_speed_global = ReduceRealMax( max_speed );
  #endif
  
}


void Grid3D::Polytropic_Star_Relaxation(  struct parameters &P  ){
  
  chprintf( "Polytropic Star Iterative Relaxation \n");
  int nfile = 0;
  Real dti, max_speed, speed_threshold;
  
  int n_step = 0;
  
  bool converged = false;
  double tdyn = sqrt( pow(P.Rstar, 3.) / P.Mstar / G_CGS);//in seconds
  chprintf("Star dynamical time: %.5e", tdyn);

  WriteData(*this, P, P.nfile);
  P.nfile++;

  while( H.t < 6.278 * tdyn){
    
    chprintf("Relaxation n_step: %d \n", n_step + 1 );
    
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
    Apply_Damping_Step( *this , tdyn);

    // Output
    if (n_step % int(P.outstep) == 0){
      WriteData(*this, P, P.nfile);
      P.nfile++;
    }

    // set boundary conditions for next time step 
    Set_Boundary_Conditions_Grid(P);
    
    chprintf("n_step: %d   sim time: %10.7f   sim timestep: %7.4e  \n\n", n_step, H.t, H.dt);
    
    //Exit the iteartions if converged
//    if ( converged ) break;


  }


  //Restart the simulation time
  H.t = 0;
  Grav.INITIAL = true;
  
  chprintf( " Finished Relaxation Iteration: %d  \n", n_step);

  
}

#endif
