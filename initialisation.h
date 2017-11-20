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
# define M_PI  3.14159265358979323846


void initialisation_domain (int N, int nx, int ny, Lattice lat, double rho_in, double rho_out, double** v_in, double** v_out, double* Uw, int D, int Q, double** xi, double sigma, double* omega_i, double** f_star, double tau, double** position);