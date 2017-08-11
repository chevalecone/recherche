// C++ Standard includes
#include <iostream> 
#include <cmath>
#include <ctime>
#include <stdio.h>

// Local includes
#include "rarefied_models.h"
#include "domain.h"
#include "lattice.h"
#include "solutionExporterEulerian.h"
#include "function.h"
#include "relaxation_time.h"

# define M_PI  3.14159265358979323846

//Cas CBBSR (pour l'instant sans wall function)
double* MRT_Guo_2008(double Kn, double ymax, double tau_s, double tau_q, double r, double dx)
{
	double* temp = new double[3];
	tau_s = 0.5+ sqrt(6/M_PI)*Kn*(ymax/dx)*0.16015557;
	double sigma = 1; //coefficient d'accommodation du moment égal à 1 dans le cas diffusif
	double A1 = (2-sigma)/sigma * (1-0.1817*sigma);
	double A2 = 1/(M_PI) + 0.5*A1*A1;
	double ki = sqrt(M_PI/6);
	double tau_s_tilde = tau_s-0.5;
	tau_q = 0.5 + (3+24*ki*ki*tau_s_tilde*tau_s_tilde*A2)/(16*tau_s_tilde);
	r = 1/(1+ki*A1);
	temp[0] = r;
	temp[1] = tau_s;
	temp[2] = tau_q;
	return temp;
}

void PHI_Guo_2008(double Kn, double ymax,int N, double dx, double* PHI, double** position, double mu, double rho, double cs)
{
	double lambda = Kn*ymax/dx;
	double* alpha = new double[N];
	double* alpha1 = new double[N];
	double* phi = new double[N];
	double* phi1 = new double[N];
	for (int j=0;j<N;j++)
	{
		alpha[j]  = position[j][1]/lambda; 
		alpha1[j] = (ymax-position[j][1])/lambda;
		phi[j]    = 1+(alpha[j]-1)*exp(-alpha[j])-alpha[j]*alpha[j]*Ei(alpha[j]);
		phi1[j]   = 1+(alpha1[j]-1)*exp(-alpha1[j])-alpha1[j]*alpha1[j]*Ei(alpha1[j]);
		PHI[j]    = 0.5*(phi[j] + phi1[j]);
		printf("alpha %f alpha1 %f phi %f phi1 %f PHI %f \n", alpha[j],alpha1[j],phi[j],phi1[j],PHI[j]);
	}
}

double* MRT_Continuous(double cs, double dt, double nu)
{
	double* temp = new double[3];
	double tau_q;
	double tau_s;
	tau_s = 0.5+ nu/(cs*cs*dt);
	tau_q = 1/1.2;
	temp[0] = 0;
	temp[1] = tau_s;
	temp[2] = tau_q;
	return temp;
}

double* MRT_Verhaeghe_2009(double Kn, double ymax, double tau_s, double tau_q, double beta, double mu, double dx, double dt, double cs)
{
	double* temp = new double[3];
	beta = (3*mu - Kn*ymax/dx)/(3*mu + Kn*ymax/dx);
	tau_s = 0.5 + sqrt(2/(3*M_PI))*Kn*ymax/(dx*cs*cs*dt);
	double s_nu = 1/tau_s;
	double s_q = 4* (2-s_nu)/(4+s_nu);
	tau_q = 1/s_q;
	temp[0] = beta;
	temp[1] = tau_s;
	temp[2] = tau_q;
	return temp;
}

