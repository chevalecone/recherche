#include "lattice.h"

double* MRT_Guo_2008_CBBSR_function(double Kn, double ymax, double tau_q, double r, double dx, int N, double** rank, double* PHI, double* A);
double* MRT_Guo_2008_CBBSR(double Kn, double ymax, double tau_s, double tau_q, double r, double dx, double* A);
double* MRT_Guo_2008_DBB(double Kn, double ymax, double tau_s, double tau_q, double beta, double dx, double* A);
double* MRT_Continuous(double cs, double dt, double nu);
double* MRT_Guo_2011_DBB(double Kn, double ymax, double beta, double dx, double* A);
double* MRT_Guo_2011_DBB_function(double Kn, double tau_q, double tau_s, double ymax, double beta, double dx, int N,double** rank,double** rank2, double* PHI, double* A);
double* MRT_Guo_2008_MR_function(double Kn, double ymax, double tau_q, double sigma1, double dx, int N, double** rank, double* PHI, double*A);
double* MRT_Yudong_2016_DBB(double Kn, double ymax, double tau_s, double tau_q, double beta, double dx, double nu);
double* MRT_Yudong_2016_DBB_function(double Kn, double ymax,double tau_s, double beta, double dx, int N, double** rank, double** rank2,double* PHI);
double* MRT_Verhaeghe_2009_DBB(double Kn, double ymax, double tau_s, double tau_q, double beta, double dx);
double* MRT_Verhaeghe_2009_DBB_function(double Kn, double ymax, double tau_s, double beta, double dx, int N, double** rank, double** rank2,double* PHI);

//Wall function
void PHI_Guo_2008(double Kn, double ymax, double dx, double* PHI, int N, double** rank);
void PHI_Zhang_2006(double Kn, double ymax, double dx, double* PHI, int N, double** rank);
void PHI_Dongari_2011(double Kn, double ymax, double dx, double* PHI, int N, double** rank);
void PHI_Arlemark_2010(double Kn, double ymax, double dx, double* PHI, int N, double** rank);
void PHI_Guo_Shu_2013(double Kn, double* PHI, int N);

//Tau-Kn relation
double Tau_Guo_2008(double Kn, double dx, double ymax);
double Tau_Tang_2005(double Kn, double dx, double ymax);

//Slip velocity
void slip_velocity_Guo_2008( double* A, double sigma);
void slip_velocity_Guo_2011( double* A, double sigma);
void slip_velocity_Hadjiconstantinou_2003(double* A, double sigma);
void slip_velocity_Wu_2008(double Kn, double* A, double sigma);
void slip_velocity_Wang_2017(double Kn, double* A, double sigma);
