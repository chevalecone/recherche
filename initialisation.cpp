// C++ Standard includes
#include <iostream> 
#include <cmath>
#include <ctime>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <omp.h>
#include <cstring>

// Local includes
#include "domain.h"
#include "lattice.h"
#include "solutionExporterEulerian.h"
#include "function.h"
#include "geometry.h"
#include "boundaryc.h"
#include "relaxation_time.h"
#include "rarefied_models.h"
#include "compute_quantities.h"
#include "slip_velocity.h"
#include "wall_function.h"
#include "initialisation.h"


# define M_PI  3.14159265358979323846

void initialisation_domain (int N, int nx, int ny, Lattice lat, double rho_in, double rho_out, double** v_in, double** v_out, double* Uw, int D, int Q, double** xi, double sigma, double* omega_i, double** f_star, double tau, double** position)
{
	double dp_dx = (rho_out-rho_in)/nx; //Linéaire
	for (int j=0;j<N;j++)
	{	
		lat.rho_[j]=rho_in+dp_dx*position[j][0];
		//Conditions inlet
		if(j%nx==0)
		{
			lat.u_[j][0] = v_in[j/nx][0];
			lat.u_[j][1] = v_in[j/nx][1];	
		}
		//Conditions outlet
		else if((j%nx)==nx-1)
		{
			lat.u_[j][0] = v_out[j/nx][0];
			lat.u_[j][1] = v_out[j/nx][1];
		}
		//Conditions Nord Couette
		else if(j>=nx*ny-nx-1)
		{
			lat.u_[j][0] = Uw[0];
			lat.u_[j][1] = Uw[1];
		}
		velocity(j,D,Q,xi,lat,sigma);
		for (int k=0;k<Q;k++)
		{
			lat.f0_[j][k] = omega_i[k]*lat.rho_[j];
			lat.f_[j][k]  = lat.f0_[j][k]; //Initialement, on aura f = fi,eq (avec vitesse nulle)
			f_star[j][k] = lat.f_[j][k]-1/tau*(lat.f_[j][k] - lat.f0_[j][k]);
		}
	}	
}