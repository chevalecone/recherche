// C++ Standard includes
#include <iostream> 
#include <cmath>
#include <ctime>
#include <stdio.h>
#include <cstring>


// Local includes
#include "rarefied_models.h"
#include "domain.h"
#include "lattice.h"
#include "solutionExporterEulerian.h"
#include "function.h"
#include "relaxation_time.h"

# define M_PI  3.14159265358979323846

double* Wall_function(std::string wfunction, double Kn, double ymax, double dx, int N, double** rank, double rho_out, Lattice lat);