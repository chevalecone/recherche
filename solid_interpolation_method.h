#include <iostream> 
#include <cmath>
#include <ctime>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <omp.h>
#include <cstring>
#include "lattice.h"
#include "function.h"

//**FONCTIONS UTILES POUR L'INTERPOLATION DE FRONTIERES COMPLEXES**//
void solid_fraction_square(int N, int Q, double** solid_fraction_interpolation, int** conn, double abscisse, double ordonnee, double diametre,bool* typeLat, double* buffer, double** position);

void solid_fraction_circular(int N, int Q, double** solid_fraction_interpolation, int** conn, double abscisse, double ordonnee, double diametre, bool* typeLat, double** position);

//******************METHODES D'INTERPOLATION*******************//
void linear_interpolation_method( int j, int Q, Lattice lat, double** f_star, int**conn, bool*  typeLat,  int* bb,double** solid_fraction_interpolation, double* tab_marquage, int cas);
void quadratic_interpolation_method( int j, int Q, Lattice lat, double** f_star, int**conn, bool*  typeLat,  int* bb,double** solid_fraction_interpolation, double* tab_marquage, int cas);
void multireflection_interpolation_method( int j, int Q, Lattice lat, double** f_star, int**conn, bool*  typeLat,  int* bb,double** solid_fraction_interpolation, double* tab_marquage, int cas, double** C, double* t, double* teq, double Fpc, double mu, double rho, double** Si);
void central_interpolation_method( int j, int Q, Lattice lat, double** f_star, int**conn, bool*  typeLat,  int* bb,double** solid_fraction_interpolation, double* tab_marquage, int cas);