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


double* rarefied_model(std::string rarefied_method, double cs, double dt, double nu, double Kn, double ymax, double r, double beta, double dx, double* A)
{
	double* temp2 = new double[3];
	if(!rarefied_method.compare("_Continuous"))
	{
		printf("Modèle utilisé : Continuous\n");
		temp2 = MRT_Continuous(cs,dt,nu);
	}
	else if(!rarefied_method.compare("Verhaeghe-2009_DBB"))
	{
		printf("Modèle utilisé : Verhaeghe-2009_DBB\n");
		temp2 = MRT_Verhaeghe_2009_DBB(Kn,ymax,beta,dx,A);
	}
	else if(!rarefied_method.compare("Guo-2008_CBBSR"))
	{
		printf("Modèle utilisé : Guo-2008_CBBSR\n");
		temp2 = MRT_Guo_2008_CBBSR(Kn,ymax,r,dx,A);
	}
	else if(!rarefied_method.compare("Guo-2008_CBBSR_v2"))
	{
		printf("Modèle utilisé : Guo-2008_CBBSR_v2\n");
		temp2 = MRT_Guo_2008_CBBSR_v2(Kn,ymax,r,dx,A);
	}
	else if(!rarefied_method.compare("Guo-2008_DBB"))
	{
		printf("Modèle utilisé : Guo-2008_DBB\n");
		temp2 = MRT_Guo_2008_DBB(Kn,ymax,beta,dx,A);
	}
	else if(!rarefied_method.compare("Yudong-2016_DBB"))
	{
		printf("Modèle utilisé : Yudong-2016_DBB\n");
		temp2 = MRT_Yudong_2016_DBB(Kn,ymax,beta,dx,nu);
	}
	else if (!rarefied_method.compare("Li-2011_CBBSR"))
	{
		printf("Modèle utilisé : Li-2011_CBBSR\n");
		temp2 = MRT_Li_2011_CBBSR(Kn,ymax,r,dx,A);
	}
	else if (!rarefied_method.compare("Guo-2011_DBB"))
	{
		printf("Modèle utilisé : Guo-2011_DBB\n");
		temp2 = MRT_Guo_2011_DBB(Kn,ymax,beta,dx,A);
	}
	return temp2;
}

double* rarefied_model_function(std::string rarefied_method, double cs, double dt, double nu, double Kn, double ymax, double r, double sigma, double beta, double dx, int N, double** rank, double** rank2, double* PHI, double* A)
{
	double* temp2 = new double[3];
	if(!rarefied_method.compare("Guo-2008_CBBSR"))
	{
		printf("Modèle utilisé : Guo-2008_CBBSR_function\n");
		temp2 = MRT_Guo_2008_CBBSR_function(Kn,ymax,r,dx,N,rank,PHI,A);
	}
	else if(!rarefied_method.compare("Guo-2008_CBBSR_v2"))
	{
		printf("Modèle utilisé : Guo-2008_CBBSR_v2_function\n");
		temp2 = MRT_Guo_2008_CBBSR_function_v2(Kn,ymax,r,dx,N,rank,PHI,A);
	}
	else if(!rarefied_method.compare("Guo-2008_MR"))
	{
		printf("Modèle utilisé : Guo-2008_MR_function\n");
		temp2 = MRT_Guo_2008_MR_function(Kn,ymax,sigma,dx,N,rank,PHI,A);
	}
	else if(!rarefied_method.compare("Verhaeghe-2009_DBB"))
	{
		printf("Modèle utilisé : Verhaeghe-2009_DBB_function\n");
		temp2 =  MRT_Verhaeghe_2009_DBB_function(Kn,ymax,beta,dx,N,rank,rank2,PHI);
	}
	else if(!rarefied_method.compare("Li-2011_CBBSR"))
	{
		printf("Modèle utilisé : Li-2011_CBBSR_function\n");
		temp2 = MRT_Li_2011_CBBSR_function(Kn,ymax,r,dx,N,rank,rank2,PHI,A);
	}
	else if(!rarefied_method.compare("Yudong-2016_DBB"))
	{
		printf("Modèle utilisé : Yudong-2016_DBB_function\n");
		temp2 = MRT_Yudong_2016_DBB_function(Kn,ymax,beta,dx,N,rank,rank2,PHI);
	}
	else if(!rarefied_method.compare("Guo-2011_DBB"))
	{
		printf("Modèle utilisé : Guo-2011_DBB_function\n");
		temp2 = MRT_Guo_2011_DBB_function(Kn,ymax,beta,dx,N,rank,rank2,PHI,A);
	}
	return temp2;
}


//*****************************************RAREFIED METHODS WITHOUT WALL FUNCTION***********************************************************//
double* MRT_Continuous(double cs, double dt, double nu)
{
	double* temp = new double[3];
	double tau_q;
	double tau_s;
	tau_s = 0.5+ 3*nu;
	tau_q = 1/1.2;
	temp[0] = 0;
	temp[1] = tau_s;
	temp[2] = tau_q;
	return temp;
}


			      //***************************************DBB*****************************************//
double* MRT_Guo_2011_DBB(double Kn, double ymax, double beta, double dx, double* A) 
{
	double A1 = A[0];
	double A2 = A[1]; 
	double* temp = new double[3];
	beta = (sqrt(6/M_PI)-A[0])/(sqrt(6/M_PI)+A[0]);
	double tau_s = 0.5 + sqrt(6/M_PI)*Kn*ymax/dx;
	double tau_q =  (4*M_PI*A[1]*(tau_s-0.5)*(tau_s-0.5)+3)/(16*(tau_s-0.5));
	temp[0] = beta;
	temp[1] = tau_s;
	temp[2] = tau_q;
	return temp;
}

double* MRT_Yudong_2016_DBB(double Kn, double ymax, double beta, double dx, double nu)
{
	double lambda = Kn*ymax/dx;
	double* temp = new double[3];
	double c_bar = sqrt(8/(3*M_PI));
	double alpha_0 = 64/(15*M_PI);
	double beta_1 = 1.2;
	double alpha_1 = 6;
	double b = -1;
	double alpha = 2/M_PI*alpha_0*atan(alpha_1*pow(Kn,beta_1));
	nu =5*M_PI/32*c_bar/(1+alpha*Kn)*lambda; 
	//beta = (3*(1-b*Kn)*nu-Kn*ymax/dx)/(3*(1-b*Kn)*nu+Kn*ymax/dx);
	double tau_s  = 0.5 + sqrt(6/M_PI)*ymax/dx*Kn/(1+alpha*Kn);
	double s_nu = 1./tau_s;
	double s_q = 8*(2-s_nu)/(8-s_nu);
	double tau_q = 1/s_q;
	
	temp[0] = beta;
	temp[1] = tau_s;
	temp[2] = tau_q;
	return temp;
	
}
double* MRT_Verhaeghe_2009_DBB(double Kn, double ymax, double beta, double dx, double* A)
{
	double* temp = new double[3];
	double A1 = A[0];
	double tau_s = Tau_Guo_2008(Kn,dx,ymax);
	double s_nu = 1/tau_s;
	double s_q = 8*(2-s_nu)/(8-s_nu);
	double mu = sqrt(2/(3*M_PI))*Kn*ymax/dx;
	beta = (sqrt(6/M_PI)-A1)/(sqrt(6/M_PI)+A1);
	double tau_q = 1/s_q;
	temp[0] = beta;
	temp[1] = tau_s;
	temp[2] = tau_q;
	return temp;
}
double* MRT_Guo_2008_DBB(double Kn, double ymax, double beta, double dx, double* A)
{
	double* temp = new double[3];
	double A1 = A[0];
	double A2 = A[1];
	beta = 2/(1+sqrt(M_PI/6)*A1);
	double tau_s = 0.5 + sqrt(6/M_PI)*Kn*(ymax/dx);
	double tau_q = 0.5 + (M_PI*A2*(2*tau_s-1)*(2*tau_s-1)+3)/(8*(2*tau_s-1));
	temp[0] = beta;
	temp[1] = tau_s;
	temp[2] = tau_q;
	return temp;
}

				  
				  //***************************************CBBSR*****************************************//
double* MRT_Guo_2008_CBBSR(double Kn, double ymax, double r, double dx, double* A)
{
	double* temp = new double[3];
	double A1 = A[0];
	double A2 = A[1];

	double ki = sqrt(M_PI/6);
	double tau_s = Tau_Guo_2008(Kn,dx,ymax);
	double tau_s0_tilde = tau_s-0.5;
	double Ts = tau_s0_tilde;
	double tau_q = (3*Ts*(8*A2*Ts*Ts*Ts*ki*ki + Ts))/(4*(4*Ts*Ts*Ts));
	//tau_q = (3+24*ki*ki*tau_s0_tilde*tau_s0_tilde*A2)/(16*tau_s0_tilde);
	//r = (2*(4*Ts*Ts*Ts))/(8*Ts*Ts*Ts + 8*A1*Ts*Ts*Ts*ki);
	r = 1/(1+ki*A1);
	temp[0] = r;
	temp[1] = tau_s;
	temp[2] = tau_q;
	return temp;
}

double* MRT_Guo_2008_CBBSR_v2(double Kn, double ymax, double r, double dx, double* A)
{
	double* temp = new double[3];
	double A1 = A[0];
	double A2 = A[1];

	double ki = sqrt(M_PI/6);
	double tau_s = Tau_Guo_2008(Kn,dx,ymax);
	double Ts = tau_s-0.5;
	double tau_q = (3*Ts*(8*A2*Ts*Ts*Ts*ki*ki + Ts))/(4*(4*Ts*Ts*Ts));
	r = (2*(4*Ts*Ts*Ts))/(8*Ts*Ts*Ts + 8*A1*Ts*Ts*Ts*ki);
	temp[0] = r;
	temp[1] = tau_s;
	temp[2] = tau_q;
	return temp;
}

double* MRT_Li_2011_CBBSR(double Kn, double ymax, double r, double dx, double* A)
{
	double lambda = Kn*ymax/dx;
	double* temp = new double[3];
	double alpha = 2;
	double sigma = 1;
	double sigma_v = (2-sigma)/sigma;
	r = 1/(1+A[0]*sigma_v*sqrt(M_PI/6));
	double tau_s  = 0.5 + sqrt(6/M_PI)*ymax/dx*Kn/(1+alpha*Kn);
	double tau_s_tilde = tau_s-0.5;
	double tau_q = 0.5+ (3+4*M_PI*tau_s_tilde*tau_s_tilde*A[1])/(16*tau_s_tilde*tau_s_tilde);
	
	temp[0] = r;
	temp[1] = tau_s;
	temp[2] = tau_q;
	return temp;
	
}


//*****************************************RAREFIED METHODS WITH WALL FUNCTION***********************************************************//

			      //***************************************DBB*****************************************//
double* MRT_Verhaeghe_2009_DBB_function(double Kn, double ymax, double beta, double dx, int N, double** rank, double** rank2,double* PHI)
{
	double* temp = new double[3];
	double tau_s = Tau_Guo_2008(Kn,dx,ymax);
	double s_nu = 1/tau_s;
	double s_q = 8*(2-s_nu)/(8-s_nu);
	//beta = (sqrt(6/M_PI)-A1)/(sqrt(6/M_PI)+A1);
	for (int i=0;i<N;i++)
	{
		rank[i][1] = Tau_Guo_2008(Kn,dx,ymax)*PHI[i];
	}
	double tau_q = 1/s_q;
	temp[0] = 0;
	temp[1] = tau_s;
	temp[2] = tau_q;
	return temp;
}

double* MRT_Yudong_2016_DBB_function(double Kn, double ymax, double beta, double dx, int N, double** rank, double** rank2,double* PHI)
{
	double lambda = Kn*ymax/dx;
	double* temp = new double[3];
	double* nu_e = new double[N+1];
	double b = -1;
	for (int i=0;i<N;i++)
	{
		nu_e[i] = 0.5*sqrt(2./(3*M_PI))*lambda*PHI[i] ;
		rank[i][1] = Tau_Guo_2008(Kn,dx,ymax)*PHI[i];
		rank2[i][1] = (8-1/rank[i][1])/(8*(2-1/rank[i][1]));
		printf("Rang : %f  tau_s : %f tau_q : %f\n",rank[i][0], rank[i][1],rank2[i][1]);
	}
	nu_e[N] = 0.5*sqrt(2./(3*M_PI))*lambda*PHI[N]; 
	beta = (3*(1-b*Kn)*nu_e[N]-Kn*ymax/dx)/(3*(1-b*Kn)*nu_e[N]+Kn*ymax/dx);
	temp[0] = beta;
	temp[1] = rank[0][1];
	temp[2] = rank2[0][1];
	return temp;
	
	
}

double* MRT_Guo_2011_DBB_function(double Kn, double ymax, double beta, double dx, int N, double** rank, double** rank2, double* PHI, double* A)
{
	double A1 = A[0];
	double A2 = A[1];
	double ki = sqrt(M_PI/6);
	double* temp = new double[3];
	//beta = (sqrt(6/M_PI)-A1)/(sqrt(6/M_PI)+A1);
	for (int i=0;i<N;i++)
	{
		rank[i][1] = Tau_Guo_2008(Kn,dx,ymax)*PHI[i];
		//rank2[i][1] =  0.5 + (M_PI*A2*(2*rank[i][1]-1)*(2*rank[i][1]-1)+3)/(8*(2*rank[i][1]-1));
	}
	double tau_s0 = Tau_Guo_2008(Kn,dx,ymax)*PHI[N];
	printf("2*(rank[0][1]-tau_s0)/dx %f %f\n",rank[0][1],tau_s0);
	double tau_s_prime = 2*(rank[0][1]-tau_s0)/dx;
	printf("tau_s_prime : %f\n", tau_s_prime);
	double tau_s0_tilde = tau_s0-0.5;
	double Ts = tau_s0_tilde;
	double Ts_p = tau_s_prime;
	beta = (8*Ts*Ts*Ts  - Ts*Ts_p*dx - 8*A1*Ts*Ts*Ts*ki + 8*Ts*Ts*Ts_p*dx + 8*A1*Ts*Ts*Ts_p*ki*dx)/(8*Ts*Ts*Ts + Ts*Ts_p*dx + 8*A1*Ts*Ts*Ts*ki + 8*Ts*Ts*Ts_p*dx + 8*A1*Ts*Ts*Ts_p*ki*dx);
	double tau_q = 0.5 + (3*(8*A2*Ts*Ts*Ts*Ts*ki*ki + 8*A2*Ts*Ts*Ts*Ts_p*ki*ki*dx  + 10*A1*Ts*Ts*Ts_p*ki*dx + Ts*Ts + 5*Ts*Ts_p*dx))/(4*(4*Ts*Ts*Ts + 4*Ts*Ts*Ts_p*dx));
	//tau_q = 1/(16*(tau_s0-0.5))*(4*M_PI*A2*(tau_s0-0.5)*(tau_s0-0.5)+3);
	
	temp[0] = beta;
	temp[1] = rank[0][1];
	temp[2] = tau_q;
	return temp;
}

				  
			      //***************************************CBBSR*****************************************//				  
double* MRT_Guo_2008_CBBSR_function(double Kn, double ymax, double r, double dx, int N, double** rank, double* PHI, double* A)
{
	double* temp = new double[3];
	double A1 = A[0];
	double A2 = A[1];
	double ki = sqrt(M_PI/6);
	for (int i=0;i<N;i++)
	{
		rank[i][1] = Tau_Guo_2008(Kn,dx,ymax)*PHI[i];
	}
	double tau_s0 = Tau_Guo_2008(Kn,dx,ymax)*PHI[N];
	printf("2*(rank[0][1]-tau_s0)/dx %f %f\n",rank[0][1],tau_s0);
	double tau_s_prime = 2*(rank[0][1]-tau_s0)/dx;
	printf("tau_s_prime : %f\n", tau_s_prime);
	double tau_s0_tilde = tau_s0-0.5;
	double tau_q = (3+24*ki*ki*tau_s0_tilde*tau_s0_tilde*A2)/(16*tau_s0_tilde)+tau_s_prime*dx*(12+30*tau_s0_tilde*ki*A1)/(16*tau_s0_tilde*tau_s0_tilde);
	r = 1/(1+ki*A1+tau_s_prime*dx/(8*tau_s0_tilde*tau_s0_tilde));
	temp[0] = r;
	temp[1] = tau_s0;
	temp[2] = tau_q;
	printf("r : %f, Ts_p : %f, Ts : %f \n",r,tau_s_prime,tau_s0);
	return temp;
}

double* MRT_Guo_2008_CBBSR_function_v2(double Kn, double ymax, double r, double dx, int N, double** rank, double* PHI, double* A)
{
	double* temp = new double[3];
	double A1 = A[0];
	double A2 = A[1];
	double ki = sqrt(M_PI/6);
	for (int i=0;i<N;i++)
	{
		rank[i][1] = Tau_Guo_2008(Kn,dx,ymax)*PHI[i];
	}
	double tau_s0 = Tau_Guo_2008(Kn,dx,ymax)*PHI[N];
	printf("2*(rank[0][1]-tau_s0)/dx %f %f\n",rank[0][1],tau_s0);
	double Ts_p = 2*(rank[0][1]-tau_s0)/dx;
	printf("tau_s_prime : %f\n", Ts_p);
	double Ts = tau_s0-0.5;
	double tau_q = 0.5 + (3*Ts*(8*A2*Ts*Ts*Ts*ki*ki + 8*A2*Ts*Ts*Ts_p*ki*ki*dx + 10*A1*Ts*Ts_p*ki*dx + Ts + 5*Ts_p*dx))/(4*(4*Ts*Ts*Ts + 4*Ts*Ts*Ts_p*dx));
	r =   (2*(4*Ts*Ts*Ts + 4*Ts*Ts*Ts_p*dx))/(8*Ts*Ts*Ts + Ts*Ts_p*dx + 8*A1*Ts*Ts*Ts*ki + 8*Ts*Ts*Ts_p*dx + 8*A1*Ts*Ts*Ts_p*ki*dx);
	temp[0] = r;
	temp[1] = tau_s0;
	temp[2] = tau_q;
	printf("r : %f, Ts_p : %f, Ts : %f \n",r,Ts_p,tau_s0);
	return temp;
}

double* MRT_Li_2011_CBBSR_function(double Kn, double ymax, double r, double dx, int N, double** rank, double** rank2,double* PHI, double* A)
{
	double lambda = Kn*ymax/dx;
	double* temp = new double[3];
	double alpha = 2;
	double sigma = 1;
	double sigma_v = (2-sigma)/sigma;
	r = 1/(1+A[0]*sigma_v*sqrt(M_PI/6));
	for (int i=0;i<N;i++)
	{
		rank[i][1] = 0.5 + sqrt(6/M_PI)*ymax/dx*Kn*PHI[i]/(1+alpha*Kn);
		rank2[i][1] = 0.5+ (3+4*M_PI*(rank[i][1]-0.5)*(rank[i][1]-0.5)*A[1])/(16*(rank[i][1]-0.5)*(rank[i][1]-0.5));
	}
	
	temp[0] = r;
	temp[1] = rank[0][1];
	temp[2] = rank2[0][1];
	return temp;
	
}

				  //***************************************MR*****************************************//	
double* MRT_Guo_2008_MR_function(double Kn, double ymax, double sigma, double dx, int N, double** rank, double* PHI, double* A)
{
	double* temp = new double[3];
	double A1 = A[0];
	double A2 = A[1];
	double ki = sqrt(M_PI/6);
	for (int i=0;i<N;i++)
	{
		rank[i][1] = Tau_Guo_2008(Kn,dx,ymax)*PHI[i];
	}
	double tau_s0 = Tau_Guo_2008(Kn,dx,ymax)*PHI[N];
	printf("2*(rank[0][1]-tau_s0)/dx %f %f\n",rank[0][1],tau_s0);
	double tau_s_prime = 2*(rank[0][1]-tau_s0)/dx;
	printf("tau_s_prime : %f\n", tau_s_prime);
	double tau_s0_tilde = tau_s0-0.5;
	double Ts = tau_s0_tilde;
	double Ts_p = tau_s_prime;
	double tau_q = 0.5 + (3*(8*A2*Ts*Ts*Ts*Ts*ki*ki + 8*A2*Ts*Ts*Ts*Ts_p*ki*ki*dx + 10*A1*Ts*Ts*Ts_p*ki*dx + Ts*Ts + 5*Ts*Ts_p*dx))/(4*(4*Ts*Ts*Ts + 4*Ts*Ts*Ts_p*dx));
	sigma = (4*(4*Ts*Ts*Ts + 4*Ts*Ts*Ts_p*dx))/(8*Ts*Ts*Ts + Ts*Ts_p*dx + 8*A1*Ts*Ts*Ts*ki + 8*Ts*Ts*Ts_p*dx + 8*A1*Ts*Ts*Ts_p*ki*dx);
	temp[0] = sigma;
	temp[1] = tau_s0;
	temp[2] = tau_q;
	printf("Ts_p : %f Ts : %f \n",tau_s_prime,tau_s0);
	return temp;
}








double Tau_Guo_2008(double Kn, double dx, double ymax)
{
	double tau_s = 0.5+sqrt(6/M_PI)*Kn*(ymax/dx);
	return tau_s;
}

double Tau_Tang_2005(double Kn, double dx, double ymax)
{
	double tau_s = 0.5+sqrt((3*M_PI/8))	*Kn*(ymax/dx);
	return tau_s;
}

