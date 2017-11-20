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
#include "wall_function.h"

# define M_PI  3.14159265358979323846

double* Wall_function(std::string wfunction, double Kn, double ymax, double dx, int N, double** rank, double rho_out, Lattice lat)
{
	double* PHI = new double[N+1];
	if(!wfunction.compare("Guo-2008"))
	{
		printf("Wall function : Guo-2008\n");
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
		return PHI;
		delete[] alpha;
		delete[] alpha1;
		delete[] phi;
		delete[] phi1;
	}
	else if (!wfunction.compare("Zhang-2006"))
	{
		printf("Wall function : Zhang-2006\n");
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
		return PHI;	
		delete[] alpha;
		delete[] alpha1;
		delete[] phi;
		delete[] phi1;			
	}
	else if (!wfunction.compare("Dongari-2011"))
	{
		printf("Wall function : Dongari-2011\n");
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
		return PHI;
		delete[] alpha;
		delete[] alpha1;
		delete[] phi;
		delete[] phi1;			
	}
	else if (!wfunction.compare("Arlemark-2010"))
	{
		printf("Wall function : Arlemark-2010\n");
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
		return PHI;
		delete[] alpha;
		delete[] alpha1;
		delete[] phi;
		delete[] phi1;
	}
	else if (!wfunction.compare("Guo_Shu-2013"))
	{
		printf("Wall function : Guo_Shu-2013\n");
		for (int i =0;i<N;i++)
		{
			PHI[i] = 2/M_PI*atan(sqrt(2)*pow(Kn,-0.750));
			printf("Rang %d, PHI %f\n",i,PHI[i]);
		}
		PHI[N] = 2/M_PI*atan(sqrt(2)*pow(Kn,-0.750));
		printf("Rang %d, PHI %f\n",N,PHI[N]);
		return PHI;
	}
}