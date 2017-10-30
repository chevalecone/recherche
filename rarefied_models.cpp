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

//Cas CBBSR (toutes les formules sont tirés de l'article de Guo_2008)
double* MRT_Guo_2008_CBBSR_function(double Kn, double ymax, double tau_q, double r, double dx, int N, double** rank, double* PHI, double* A)
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
	//tau_q = (3+24*ki*ki*tau_s0_tilde*tau_s0_tilde*A2)/(16*tau_s0_tilde)+tau_s_prime*dx*(12+30*tau_s0_tilde*ki*A1)/(16*tau_s0_tilde*tau_s0_tilde);
	//(48*A2*Ts^4*X^2 + 168*A2*Ts^3*Ts_p*X^2*dx - 48*A1*Ts^3*X - 48*A1*Ts^2*Ts_p*X*dx - 6*Ts^2 + 21*Ts*Ts_p*dx)/(64*(Ts^3 + 3*Ts^2*Ts_p*dx))
	//(3*Ts*(8*A2*Ts^3*X^2 + 8*A2*Ts^2*Ts_p*X^2*dx + 10*A1*Ts*Ts_p*X*dx + Ts + 5*Ts_p*dx))/(4*(4*Ts^3 + 4*Ts^2*Ts_p*dx))
	tau_q = 0.5 + (3*Ts*(8*A2*Ts*Ts*Ts*ki*ki + 8*A2*Ts*Ts*Ts_p*ki*ki*dx + 10*A1*Ts*Ts_p*ki*dx + Ts + 5*Ts_p*dx))/(4*(4*Ts*Ts*Ts + 4*Ts*Ts*Ts_p*dx));
	//r = 1/(1+ki*A1+tau_s_prime*dx/(8*tau_s0_tilde*tau_s0_tilde));
	//(2*(4*Ts^3 + 4*Ts^2*Ts_p*dx))/(8*Ts^3 + Ts*Ts_p*dx + 8*A1*Ts^3*X + 8*Ts^2*Ts_p*dx + 8*A1*Ts^2*Ts_p*X*dx)
	r =   (2*(4*Ts*Ts*Ts + 4*Ts*Ts*Ts_p*dx))/(8*Ts*Ts*Ts + Ts*Ts_p*dx + 8*A1*Ts*Ts*Ts*ki + 8*Ts*Ts*Ts_p*dx + 8*A1*Ts*Ts*Ts_p*ki*dx);
	temp[0] = r;
	temp[1] = tau_s0;
	temp[2] = tau_q;
	printf("Ts_p : %f Ts : %f \n",tau_s_prime,tau_s0);
	return temp;
}

double* MRT_Guo_2008_MR_function(double Kn, double ymax, double tau_q, double sigma1, double dx, int N, double** rank, double* PHI, double* A)
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
	tau_q = 0.5 + (3*(8*A2*Ts*Ts*Ts*Ts*ki*ki + 8*A2*Ts*Ts*Ts*Ts_p*ki*ki*dx + 10*A1*Ts*Ts*Ts_p*ki*dx + Ts*Ts + 5*Ts*Ts_p*dx))/(4*(4*Ts*Ts*Ts + 4*Ts*Ts*Ts_p*dx));
	sigma1 = (4*(4*Ts*Ts*Ts + 4*Ts*Ts*Ts_p*dx))/(8*Ts*Ts*Ts + Ts*Ts_p*dx + 8*A1*Ts*Ts*Ts*ki + 8*Ts*Ts*Ts_p*dx + 8*A1*Ts*Ts*Ts_p*ki*dx);
	temp[0] = sigma1;
	temp[1] = tau_s0;
	temp[2] = tau_q;
	printf("Ts_p : %f Ts : %f \n",tau_s_prime,tau_s0);
	return temp;
}

double* MRT_Verhaeghe_2009_DBB(double Kn, double ymax, double tau_s, double tau_q, double beta, double dx)
{
	double* temp = new double[3];
	double A1 = 1.146;
	tau_s = Tau_Guo_2008(Kn,dx,ymax);
	double s_nu = 1/tau_s;
	double s_q = 8*(2-s_nu)/(8-s_nu);
	double mu = sqrt(2/(3*M_PI))*Kn*ymax/dx;
	beta = (sqrt(6/M_PI)-A1)/(sqrt(6/M_PI)+A1);
	tau_q = 1/s_q;
	temp[0] = 0;
	temp[1] = tau_s;
	temp[2] = tau_q;
	return temp;
}

double* MRT_Verhaeghe_2009_DBB_function(double Kn, double ymax, double tau_s, double beta, double dx, int N, double** rank, double** rank2,double* PHI)
{
	double* temp = new double[3];
	tau_s = Tau_Guo_2008(Kn,dx,ymax);
	double s_nu = 1/tau_s;
	double s_q = 8*(2-s_nu)/(8-s_nu);
	//beta = (sqrt(6/M_PI)-A1)/(sqrt(6/M_PI)+A1);
	for (int i=0;i<N;i++)
	{
		rank[i][1] = Tau_Guo_2008(Kn,dx,ymax)*PHI[i];
		rank2[i][1] = (8-1/rank[i][1])/(8*(2-1/rank[i][1]));
	}
	double tau_q = 1/s_q;
	temp[0] = 0;
	temp[1] = tau_s;
	temp[2] = tau_q;
	return temp;
}

double* MRT_Guo_2008_CBBSR(double Kn, double ymax, double tau_s, double tau_q, double r, double dx, double* A)
{
	double* temp = new double[3];
	double A1 = A[0];
	double A2 = A[1];

	double ki = sqrt(M_PI/6);
	tau_s = Tau_Guo_2008(Kn,dx,ymax);
	double tau_s0_tilde = tau_s-0.5;
	double Ts = tau_s0_tilde;
	tau_q = (3*Ts*(8*A2*Ts*Ts*Ts*ki*ki + Ts))/(4*(4*Ts*Ts*Ts));
	//tau_q = (3+24*ki*ki*tau_s0_tilde*tau_s0_tilde*A2)/(16*tau_s0_tilde);
	r = (2*(4*Ts*Ts*Ts))/(8*Ts*Ts*Ts + 8*A1*Ts*Ts*Ts*ki);
	r = 1/(1+ki*A1);
	temp[0] = r;
	temp[1] = tau_s;
	temp[2] = tau_q;
	return temp;
}

double* MRT_Guo_2008_DBB(double Kn, double ymax, double tau_s, double tau_q, double beta, double dx, double* A)
{
	double* temp = new double[3];
	double A1 = A[0];
	double A2 = A[1];
	beta = 2/(1+sqrt(M_PI/6)*A1);
	tau_s = 0.5 + sqrt(6/M_PI)*Kn*(ymax/dx);
	tau_q = 0.5 + (M_PI*A2*(2*tau_s-1)*(2*tau_s-1)+3)/(8*(2*tau_s-1));
	temp[0] = beta;
	temp[1] = tau_s;
	temp[2] = tau_q;
	return temp;
}

double* MRT_Yudong_2016_DBB(double Kn, double ymax, double tau_s, double tau_q, double beta, double dx, double nu)
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
	tau_s  = 0.5 + sqrt(6/M_PI)*ymax/dx*Kn/(1+alpha*Kn);
	double s_nu = 1./tau_s;
	double s_q = 8*(2-s_nu)/(8-s_nu);
	tau_q = 1/s_q;
	
	temp[0] = beta;
	temp[1] = tau_s;
	temp[2] = tau_q;
	return temp;
	
}

double* MRT_Li_2011_CBBSR(double Kn, double ymax, double tau_s, double tau_q, double r, double dx, double* A)
{
	double lambda = Kn*ymax/dx;
	double* temp = new double[3];
	double alpha = 2;
	double sigma = 1;
	double sigma_v = (2-sigma)/sigma;
	r = 1/(1+A[0]*sigma_v*sqrt(M_PI/6));
	tau_s  = 0.5 + sqrt(6/M_PI)*ymax/dx*Kn/(1+alpha*Kn);
	double tau_s_tilde = tau_s-0.5;
	tau_q = 0.5+ (3+4*M_PI*tau_s_tilde*tau_s_tilde*A[1])/(16*tau_s_tilde*tau_s_tilde);
	
	temp[0] = r;
	temp[1] = tau_s;
	temp[2] = tau_q;
	return temp;
	
}

double* MRT_Li_2011_CBBSR_function(double Kn, double ymax, double tau_s, double r, double dx, int N, double** rank, double** rank2,double* PHI, double* A)
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
	temp[1] = tau_s;
	return temp;
	
}


double* MRT_Yudong_2016_DBB_function(double Kn, double ymax, double tau_s, double beta, double dx, int N, double** rank, double** rank2,double* PHI)
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




double* MRT_Guo_2011_DBB(double Kn, double ymax, double beta, double dx, double* A) //Pareil que Verhaeghe-2009 (mais body force driven)
{
	double L1 = A[0];
	double L2 = A[1]; 
	double* temp = new double[3];
	beta = (sqrt(6/M_PI)-L1)/(sqrt(6/M_PI)+L1);
	double tau_s = 0.5 + sqrt(6/M_PI)*Kn*ymax/dx;
	double tau_q = 1/(16*(tau_s-0.5))*(4*M_PI*L2*(tau_s-0.5)*(tau_s-0.5)+3);
	temp[0] = beta;
	temp[1] = tau_s;
	temp[2] = tau_q;
	return temp;
}


double* MRT_Guo_2011_DBB_function(double Kn, double tau_q, double tau_s, double ymax, double beta, double dx, int N, double** rank, double** rank2, double* PHI, double* A)
{
	double A1 = A[0];
	double A2 = A[1];
	double ki = sqrt(M_PI/6);
	double* temp = new double[3];
	//beta = (sqrt(6/M_PI)-A1)/(sqrt(6/M_PI)+A1);
	for (int i=0;i<N;i++)
	{
		rank[i][1] = Tau_Guo_2008(Kn,dx,ymax)*PHI[i];
		rank2[i][1] = (8-1/rank[i][1])/(8*(2-1/rank[i][1]));
	}
	double tau_s0 = Tau_Guo_2008(Kn,dx,ymax)*PHI[N];
	printf("2*(rank[0][1]-tau_s0)/dx %f %f\n",rank[0][1],tau_s0);
	double tau_s_prime = 2*(rank[0][1]-tau_s0)/dx;
	printf("tau_s_prime : %f\n", tau_s_prime);
	double tau_s0_tilde = tau_s0-0.5;
	double Ts = tau_s0_tilde;
	double Ts_p = tau_s_prime;
	beta = (8*Ts*Ts*Ts  - Ts*Ts_p*dx - 8*A1*Ts*Ts*Ts*ki + 8*Ts*Ts*Ts_p*dx + 8*A1*Ts*Ts*Ts_p*ki*dx)/(8*Ts*Ts*Ts + Ts*Ts_p*dx + 8*A1*Ts*Ts*Ts*ki + 8*Ts*Ts*Ts_p*dx + 8*A1*Ts*Ts*Ts_p*ki*dx);
	tau_q = 0.5 + (3*(8*A2*Ts*Ts*Ts*Ts*ki*ki + 8*A2*Ts*Ts*Ts*Ts_p*ki*ki*dx  + 10*A1*Ts*Ts*Ts_p*ki*dx + Ts*Ts + 5*Ts*Ts_p*dx))/(4*(4*Ts*Ts*Ts + 4*Ts*Ts*Ts_p*dx));
	
	temp[0] = beta;
	temp[1] = rank[0][1];
	temp[2] = tau_q;
	return temp;
}

void PHI_Guo_2008(double Kn, double ymax, double dx, double* PHI, int N, double** rank)
{
	double lambda = Kn*ymax/dx;
	double* alpha = new double[N+1];
	double* alpha1 = new double[N+1];
	double* phi = new double[N+1];
	double* phi1 = new double[N+1];
	for (int j=0;j<N;j++)
	{
		alpha[j]  = rank[j][0]/lambda; 
		alpha1[j] = (ymax-rank[j][0])/lambda;
		phi[j]    = 1+(alpha[j]-1)*exp(-alpha[j])-alpha[j]*alpha[j]*Ei_big(1,alpha[j]);
		phi1[j]   = 1+(alpha1[j]-1)*exp(-alpha1[j])-alpha1[j]*alpha1[j]*Ei_big(1,alpha1[j]);
		PHI[j]    = 0.5*(phi[j] + phi1[j]);
		printf("Rang %d, PHI %f\n",j,PHI[j]);
	}
	alpha[N] = 0;
	alpha1[N] = ymax/lambda;
	phi[N]    = 1+(alpha[N]-1)*exp(-alpha[N])-0;
	phi1[N]    = 1+(alpha1[N]-1)*exp(-alpha1[N])-alpha1[N]*alpha1[N]*Ei_big(1,alpha1[N]);
	PHI[N]    = 0.5*(phi[N] + phi1[N]);
	printf("Rang %d, PHI %f\n",N,PHI[N]);
}

void PHI_Zhang_2006(double Kn, double ymax, double dx, double* PHI, int N, double** rank)
{
	double C = 1;
	double lambda = Kn*ymax/dx;
	double* alpha = new double[N+1];
	double* alpha1 = new double[N+1];
	double* phi = new double[N+1];
	double* phi1 = new double[N+1];
	for (int j=0;j<N;j++)
	{
		alpha[j]  = rank[j][0]/lambda; 
		alpha1[j] = (ymax-rank[j][0])/lambda;
		phi[j]    = exp(-C*alpha[j]);
		phi1[j]   = exp(-C*alpha1[j]);
		PHI[j]    = 1/(1+0.7*(phi[j]+phi1[j]));
	}
		alpha[N]  = 0; 
		alpha1[N] = ymax/lambda;
		phi[N]    = exp(-C*alpha[N]);
		phi1[N]   = exp(-C*alpha1[N]);
		PHI[N]    = 1/(1+0.7*(phi[N]+phi1[N]));
}


void PHI_Dongari_2011(double Kn, double ymax, double dx, double* PHI, int N, double** rank)
{
	double lambda   = Kn*ymax/dx;
	double* alpha   = new double[N];
	double* alpha1  = new double[N];
	double* phi   = new double[N];
	double* phi1  = new double[N];
	double sum1,sum2,sum3,sum4;
	
	for (int j=0;j<N;j++)
	{
		alpha[j]  = rank[j][0]/lambda; 
		alpha1[j] = (ymax-rank[j][0])/lambda;
		for (int i = 1;i<9;i++)
		{
			sum1+=pow(1+alpha[j]/(cos((2*i-1)*M_PI/32)),-2);
			sum2+=pow(1+alpha1[j]/(cos((2*i-1)*M_PI/32)),-2);
		}
		for (int i =1;i<8;i++)
		{
			sum3+=pow(1+alpha[j]/cos(i*M_PI/16),-2);
			sum4+=pow(1+alpha1[j]/cos(i*M_PI/16),-2);
		}
		phi[j]    = pow(1+alpha[j],-2.)+4*sum1+2*sum3;
		phi1[j]   = pow(1+alpha1[j],-2.)+4*sum2+2*sum4;
		PHI[j]    = 1-1/96.*(phi[j]+phi1[j]);
		printf("Rang %d, PHI %f\n",j,PHI[j]);
		sum1=0;
		sum2=0;
		sum3=0;
		sum4=0;
	}
		alpha[N]  = 0; 
		alpha1[N] = ymax/lambda;
		for (int i = 1;i<9;i++)
		{
			sum1+=pow(1+alpha[N]/(cos((2*i-1)*M_PI/32)),-2);
			sum2+=pow(1+alpha1[N]/(cos((2*i-1)*M_PI/32)),-2);
		}
		for (int i =1;i<8;i++)
		{
			sum3+=pow(1+alpha[N]/cos(i*M_PI/16),-2);
			sum4+=pow(1+alpha1[N]/cos(i*M_PI/16),-2);
		}
		phi[N]    = pow(1+alpha[N],-2.)+4*sum1+2*sum3;
		phi1[N]   = pow(1+alpha1[N],-2.)+4*sum2+2*sum4;
		PHI[N]    = 1-1/96.*(phi[N]+phi1[N]);
		printf("Rang %d, PHI %f\n",N,PHI[N]);
}


void PHI_Arlemark_2010(double Kn, double ymax, double dx, double* PHI, int N, double** rank)
{
	double lambda   = Kn*ymax/dx;
	double* alpha   = new double[N];
	double* alpha1  = new double[N];
	double* phi   = new double[N];
	double* phi1  = new double[N];
	double sum1,sum2,sum3,sum4;
	
	for (int j=0;j<N;j++)
	{
		alpha[j]  =  (0.5*ymax+rank[j][0])/lambda; 
		if((0.5*ymax-rank[j][0])/lambda<0)
		{
			alpha1[j] =  -(0.5*ymax-rank[j][0])/lambda;
		}
		else
		{
			alpha1[j] = (0.5*ymax-rank[j][0])/lambda;
		}
		for (int i = 1;i<8;i++)
		{
			sum1+=exp(-alpha[j]/(cos((2*i-1)*M_PI/28)));
			sum2+=exp(-alpha1[j]/(cos((2*i-1)*M_PI/28)));
		}
		for (int i =1;i<7;i++)
		{
			sum3+=exp(-alpha[j]/cos(i*M_PI/14));
			sum4+=exp(-alpha1[j]/cos(i*M_PI/14));
		}
		phi[j]    = exp(-alpha[j])+4*sum1+2*sum3;
		phi1[j]   = exp(-alpha1[j])+4*sum2+2*sum4;
		PHI[j]    = 1-1/82.*(phi[j]+phi1[j]);
		printf("Rang %d, PHI %f\n",j,PHI[j]);
		sum1=0;
		sum2=0;
		sum3=0;
		sum4=0;
	}
	alpha[N]  =  0.5*ymax/lambda; 
	alpha1[N] =  0.5*ymax/lambda;
	for (int i = 1;i<8;i++)
		{
			sum1+=exp(-alpha[N]/(cos((2*i-1)*M_PI/28)));
			sum2+=exp(-alpha1[N]/(cos((2*i-1)*M_PI/28)));
		}
		for (int i =1;i<7;i++)
		{
			sum3+=exp(-alpha[N]/cos(i*M_PI/14));
			sum4+=exp(-alpha1[N]/cos(i*M_PI/14));
		}
		phi[N]    = exp(-alpha[N])+4*sum1+2*sum3;
		phi1[N]   = exp(-alpha1[N])+4*sum2+2*sum4;
		PHI[N]    = 1-1/82.*(phi[N]+phi1[N]);
		printf("Rang %d, PHI %f\n",N,PHI[N]);
}
void PHI_Guo_Shu_2013(double Kn, double* PHI, int N)
{
	for (int i =0;i<N;i++)
	{
		PHI[i] = 2/M_PI*atan(sqrt(2)*pow(Kn,-0.750));
		printf("Rang %d, PHI %f\n",i,PHI[i]);
	}
	PHI[N] = 2/M_PI*atan(sqrt(2)*pow(Kn,-0.750));
	printf("Rang %d, PHI %f\n",N,PHI[N]);
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

void slip_velocity_Wang_2017(double Kn, double* A, double sigma)
{
	double nA1 = 2-sigma+sigma*(-1+2*sigma)/(Kn*Kn*Kn)*(-Ei_big(1,1/Kn))+sigma*exp(-1/Kn)*((-1+2*sigma)/(Kn*Kn)+(1-2*sigma)/Kn+(-2+sigma))-(4-4*sigma)/(Kn*Kn*Kn)*(-Ei_big(1,2/Kn))-exp(-2/Kn)*((2-2*sigma)/(Kn*Kn)-(1-sigma)/Kn+(1-sigma));
	double dA1 = (2-sigma)*(1+sigma/(Kn*Kn)*(-Ei_big(1,1/Kn))+exp(-1/Kn)*(sigma/Kn-sigma)+(4-4*sigma)/(Kn*Kn)*(-Ei_big(1,2/Kn))+exp(-2/Kn)*((2-2*sigma)/Kn-(1-sigma)));
	double nA2 = 6+(-8+11*sigma)/(Kn*Kn*Kn*Kn)*(-Ei_big(1,1/Kn))+exp(-1/Kn)*((-8+11*sigma)/(Kn*Kn*Kn)+(8-11*sigma)/(Kn*Kn)+(-16+10*sigma)/Kn-6*sigma)+(16-16*sigma)/(Kn*Kn*Kn*Kn)*(-Ei_big(1,2/Kn))+exp(-2/Kn)*((8-8*sigma)/(Kn*Kn*Kn)-(4-4*sigma)/(Kn*Kn)+(4-4*sigma)/Kn-(6-6*sigma));
	double dA2 = 12+12*sigma/(Kn*Kn)*(-Ei_big(1,1/Kn))+exp(-1/Kn)*(12*sigma/Kn-12*sigma)+(48-48*sigma)/(Kn*Kn)*(-Ei_big(1,2/Kn))+exp(-2/Kn)*((24-24*sigma)/Kn-(12-12*sigma));
	printf("nA1 dA1 nA2 dA2 : %f %f %f %f\n",nA1,dA1,nA2,dA2);
	double A1 = 2./3.*nA1/dA1;
	double A2 = nA2/dA2;
	A[0] = A1;
	A[1] = A2;
}

void slip_velocity_Guo_2008(double* A, double sigma)
{
	A[0] = (2-sigma)/sigma * (1-0.1817*sigma);
	A[1] = 1/(M_PI) + 0.5*A[0]*A[0];
}
void slip_velocity_Guo_2011(double* A, double sigma)
{
	A[0] = 1.146;
	A[1] = 0.976;
}
void slip_velocity_Hadjiconstantinou_2003(double* A, double sigma)
{
	A[0] = 1.11*(2-sigma)/sigma;
	A[1] = 0.61;
}
void slip_velocity_Wu_2008(double Kn, double* A, double sigma)
{
	double f;
	if (Kn>1)
	{
		f = 1/Kn;
	}
	else
	{
		f = 1;
	}
	A[0] = 2./3.*(3-sigma*f*f*f/(sigma)-3*(1-f*f)/(2*Kn));
	A[1] = 1./4.*(f*f*f*f+2/(Kn*Kn)*(1-f*f));
}

void slip_velocity_Li_2011(double Kn, double* A, double sigma)
{
	A[0] = (1-0.1817*sigma);
	A[1] = 0.8;
}