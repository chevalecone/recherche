#include "lattice.h"
#include <cstring>

//RAREFIED MODELS (WITH OR WITHOUT WALL FUNCTION)
double* MRT_Guo_2008_CBBSR_function(double Kn, double ymax, double r, double dx, int N, double** rank, double* PHI, double* A);
double* MRT_Guo_2008_CBBSR_function_v2(double Kn, double ymax, double r, double dx, int N, double** rank, double* PHI, double* A);
double* MRT_Guo_2008_CBBSR(double Kn, double ymax, double r, double dx, double* A);
double* MRT_Guo_2008_CBBSR_v2(double Kn, double ymax, double r, double dx, double* A);
double* MRT_Guo_2008_DBB(double Kn, double ymax, double beta, double dx, double* A);
double* MRT_Continuous(double cs, double dt, double nu);
double* MRT_Guo_2011_DBB(double Kn, double ymax, double beta, double dx, double* A);
double* MRT_Guo_2011_DBB_function(double Kn, double ymax, double beta, double dx, int N,double** rank,double** rank2, double* PHI, double* A);
double* MRT_Guo_2008_MR_function(double Kn, double ymax, double sigma, double dx, int N, double** rank, double* PHI, double*A);
double* MRT_Yudong_2016_DBB(double Kn, double ymax, double beta, double dx, double nu);
double* MRT_Yudong_2016_DBB_function(double Kn, double ymax, double beta, double dx, int N, double** rank, double** rank2,double* PHI);
double* MRT_Verhaeghe_2009_DBB(double Kn, double ymax, double beta, double dx, double* A);
double* MRT_Verhaeghe_2009_DBB_function(double Kn, double ymax, double beta, double dx, int N, double** rank, double** rank2,double* PHI);
double* MRT_Li_2011_CBBSR(double Kn, double ymax, double r, double dx, double* A);
double* MRT_Li_2011_CBBSR_function(double Kn, double ymax, double r, double dx, int N, double** rank, double** rank2,double* PHI, double* A);


//Tau-Kn relation
double Tau_Guo_2008(double Kn, double dx, double ymax);
double Tau_Tang_2005(double Kn, double dx, double ymax);

//FONCTIONS DE SIMPLIFICATION DU MAIN
double* rarefied_model(std::string rarefied_method, double cs, double dt, double nu, double Kn, double ymax, double r, double beta, double dx, double* A);
double* rarefied_model_function(std::string rarefied_method, double cs, double dt, double nu, double Kn, double ymax, double r, double sigma, double beta, double dx, int N, double** rank, double** rank2, double* PHI, double* A);