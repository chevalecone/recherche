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

//Cas CBBSR (toutes les formules sont tirés de l'article de Guo_2008)
double* MRT_Guo_2008(double Kn, double ymax, double tau_q, double r, double dx, int N, double** rank, double*PHI)
{
	double sigma = 1; //coefficient d'accommodation du moment égal à 1 dans le cas diffusif
	double lambda = Kn*ymax/dx;
	double tau_s0 = 0.5+sqrt(6/M_PI)*Kn*(ymax/dx)*(0.5*(1+(0-1)*exp(-0)+1+(ymax/lambda-1)*exp(-ymax/lambda)-ymax*ymax/(lambda*lambda)*Ei_big(1,ymax/lambda)));
	
	printf("tau_s0 : %f\n", tau_s0);
	double* temp = new double[3];
	double A1 = (2-sigma)/sigma * (1-0.1817*sigma);
	double A2 = 1/(M_PI) + 0.5*A1*A1;
	double ki = sqrt(M_PI/6);
	for (int i=0;i<N;i++)
	{
		rank[i][1] = 0.5+ sqrt(6/M_PI)*Kn*(ymax/dx)*PHI[i];
	}
	double tau_s_prime = 2*(rank[0][1]-tau_s0)/dx;
	printf("tau_s_prime : %f\n", tau_s_prime);
	double tau_s0_tilde = tau_s0-0.5;
	tau_q = (3+24*ki*ki*tau_s0_tilde*tau_s0_tilde*A2)/(16*tau_s0_tilde)+tau_s_prime*dx*(12+30*tau_s0_tilde*ki*A1)/(16*tau_s0_tilde*tau_s0_tilde);
	r = 1/(1+ki*A1+tau_s_prime*dx/(8*tau_s0_tilde*tau_s0_tilde));
	temp[0] = r;
	temp[1] = tau_s0;
	temp[2] = tau_q;
	return temp;
}

void PHI_Guo_2008(double Kn, double ymax, double dx, double* PHI, double** position, double mu, double rho, double cs, int N, double** rank)
{
	double lambda = Kn*ymax/dx;
	double* alpha = new double[N];
	double* alpha1 = new double[N];
	double* phi = new double[N];
	double* phi1 = new double[N];
	for (int j=0;j<N;j++)
	{
		alpha[j]  = rank[j][0]/lambda; 
		alpha1[j] = (ymax-rank[j][0])/lambda;
		phi[j]    = 1+(alpha[j]-1)*exp(-alpha[j])-alpha[j]*alpha[j]*Ei_big(1,alpha[j]);
		phi1[j]   = 1+(alpha1[j]-1)*exp(-alpha1[j])-alpha1[j]*alpha1[j]*Ei_big(1,alpha1[j]);
		PHI[j]    = 0.5*(phi[j] + phi1[j]);
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

