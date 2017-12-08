// C++ Standard includes
#include <iostream> 
#include <cmath>
#include <ctime>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <time.h>
#include <omp.h>
#include <cstring>

// Local includes
#include "domain.h"
#include "lattice.h"
#include "solutionExporterEulerian.h"
#include "function.h"
#include "solid_interpolation_method.h"
# define M_PI  3.14159265358979323846

#define EULER 0.5772156649 //Constante d'Euler
#define MAXIT 100 //Nombre maximum d'itérations
#define FPMIN 1.0e-30
#define EPS 1.0e-8 


void solid_fraction_square(int N, int Q, double** solid_fraction_interpolation, int** conn, double abscisse, double ordonnee, double diametre, bool* typeLat, double* buffer, double** position)
{
	
	double xC = abscisse, yC = ordonnee,xB,yB,xW,yW,xF,yF,lon,R1,R2, R = 0.5*diametre, q1,q2, buf;
	for (int j =0;j<N;j++)
	{
		for (int k=1;k<Q;k++)
		{
			if ( conn[j][k]!=-1 && typeLat[j] && (!typeLat[conn[j][k]]))
			{
				xB = position[j][0];
				yB = position[j][1];
				xF = position[conn[j][k]][0];
				yF = position[conn[j][k]][1];
				lon = sqrt((xF-xB)*(xF-xB)+(yF-yB)*(yF-yB));

				if(k==1)
				{
					xW = abscisse+R;
					solid_fraction_interpolation[j][k] = (xF-xW)/lon;
					//printf("lattice %d k = %d q : %f\n",j,k,solid_fraction_interpolation[j][k]);
				}
				else if(k==2)
				{
					yW = ordonnee+R;
					solid_fraction_interpolation[j][k] = (yF-yW)/lon;
					//printf("lattice %d k = %d q : %f\n",j,k,solid_fraction_interpolation[j][k]);
				}
				else if(k==3)
				{

					xW = abscisse-R;
					solid_fraction_interpolation[j][k] = (xW-xF)/lon;
					//printf("lattice %d k = %d q : %f\n",j,k,solid_fraction_interpolation[j][k]);
				}
				else if(k==4)
				{
					yW = ordonnee-R;
					solid_fraction_interpolation[j][k] = (yW-yF)/lon;
					//printf("lattice %d k = %d q : %f\n",j,k,solid_fraction_interpolation[j][k]);
				}
				else if(k==5)
				{

					xW = abscisse+R;
					buf = xW-xB;
					yW = yB+buf;
					q1 = sqrt((xF-xW)*(xF-xW)+(yF-yW)*(yF-yW))/lon;
					yW = ordonnee+R;
					buf = yW-yB;
					xW = xB+buf;
					q2 = sqrt((xF-xW)*(xF-xW)+(yF-yW)*(yF-yW))/lon;
					if(q1<1)
					{
						solid_fraction_interpolation[j][k] = q1;
					}
					else
					{
						solid_fraction_interpolation[j][k] = q2;
					}
					//printf("k = %d q1 : %f q2 = %f\n",k,q1,q2);
				}
				else if(k==6)
				{
					xW = abscisse-R;
					buf = xB-xW;
					yW = yB+buf;
					q1 = sqrt((xF-xW)*(xF-xW)+(yF-yW)*(yF-yW))/lon;
					yW = ordonnee+R;
					buf = yW-yB;
					xW = xB-buf;
					q2 = sqrt((xF-xW)*(xF-xW)+(yF-yW)*(yF-yW))/lon;
					if(q1<1)
					{
						solid_fraction_interpolation[j][k] = q1;
					}
					else
					{
						solid_fraction_interpolation[j][k] = q2;
					}
					//printf("k = %d q1 : %f q2 = %f\n",k,q1,q2);
				}
				else if(k==7)
				{
					xW = abscisse-R;
					buf = xB-xW;
					yW = yB-buf;
					q1 = sqrt((xF-xW)*(xF-xW)+(yF-yW)*(yF-yW))/lon;
					yW = ordonnee-R;
					buf = yB-yW;
					xW = xB-buf;
					q2 = sqrt((xF-xW)*(xF-xW)+(yF-yW)*(yF-yW))/lon;
					if(q1<1)
					{
						solid_fraction_interpolation[j][k] = q1;
					}
					else
					{
						solid_fraction_interpolation[j][k] = q2;
					}					
					//printf("k = %d q1 : %f q2 = %f\n",k,q1,q2);
				}
				else if(k==8)
				{
					xW = abscisse+R;
					buf = xW-xB;
					yW = yB-buf;
					q1 = sqrt((xF-xW)*(xF-xW)+(yF-yW)*(yF-yW))/lon;
					yW = ordonnee-R;
					buf = yB-yW;
					xW = xB+buf;
					q2 = sqrt((xF-xW)*(xF-xW)+(yF-yW)*(yF-yW))/lon;
					if(q1<1)
					{
						solid_fraction_interpolation[j][k] = q1;
					}
					else
					{
						solid_fraction_interpolation[j][k] = q2;
					}
					//printf("k = %d q1 : %f q2 = %f\n",k,q1,q2);
				}
			}
		}
	}
	
}

//Loi des sinus
void solid_fraction_circular(int N, int Q, double** solid_fraction_interpolation, int** conn, double abscisse, double ordonnee, double diametre, bool* typeLat, double** position)
{
	double xC = abscisse, yC = ordonnee,xB,yB,xW,yW,xF,yF,lon,R1,R2, R = 0.5*diametre,x, D1, gamma;
	double cT1, cT2, sT1, sT2, T1,T2;
	for (int j =0;j<N;j++)
	{
		for (int k=1;k<Q;k++)
		{
			if ( conn[j][k]!=-1 && typeLat[j] && (!typeLat[conn[j][k]]))
			{
					xB = position[j][0];
					yB = position[j][1];
					xF = position[conn[j][k]][0];
					yF = position[conn[j][k]][1];
					lon = sqrt((xF-xB)*(xF-xB)+(yF-yB)*(yF-yB));
					R1 = sqrt((xB-xC)*(xB-xC)+(yB-yC)*(yB-yC));
					R2 = sqrt((xF-xC)*(xF-xC)+(yF-yC)*(yF-yC));
					cT1 = (xB-xC)/R1;
					cT2 = (xF-xC)/R2;
					sT1 = (yB-yC)/R1;


					
					if(cT1>0 && sT1>0)
					{
						T1 = acos(cT1);
					}
					else if(cT1>0 && sT1<0)
					{
						T1 = asin(sT1);
					}
					else if(cT1<0 && sT1<0)
					{
						T1 = -M_PI-asin(sT1);
					}
					else if(cT1<0 && sT1>0)
					{
						T1 = acos(cT1);
					}
					else if(cT1 ==1 && sT1 ==0)
					{
						T1 = 0;
					}
					else if(cT1 ==-1 && sT1 ==0)
					{
						T1 = M_PI;
					}
					else if(cT1 ==0 && sT1 ==1)
					{
						T1 = M_PI/2;
					}
					else if(cT1 ==0 && sT1 ==-1)
					{
						T1 = 3*M_PI/2;
					}
					
					if(k==1)
					{
						T2 = 0;
					}
					else if(k==2)
					{
						T2 = M_PI/2;
					}
					else if(k==3)
					{
						T2 = M_PI;
					}
					else if(k==4)
					{
						T2 = 3*M_PI/2;
					}
					else if(k==5)
					{
						T2 = M_PI/4;
					}
					else if(k==6)
					{
						T2 =3*M_PI/4; 
					}
					else if(k==7)
					{
						T2 = -3*M_PI/4;
					}
					else if(k==8)
					{
						T2 = -M_PI/4;
					}
					gamma = cos(T1)*cos(T2)+sin(T1)*sin(T2);
					x = -gamma*R1+sqrt(gamma*gamma*R1*R1-(R1*R1-R*R));
				
					xW = R1*cos(T1)+x*cos(T2);
					yW = R1*sin(T1)+x*sin(T2);
					solid_fraction_interpolation[j][k] = (lon-x)/lon;
					
			}
		}
	}
	
}

void linear_interpolation_method( int j, int Q, Lattice lat, double** f_star, int**conn, bool*  typeLat,  int* bb,double** solid_fraction_interpolation, double* tab_marquage, int cas)
{
	double q;
	int xf1, xf2,xb;
	for (int k=0;k<Q;k++)
	{
		//Si la lattice j est fluide, et le voisin est solide
		if ( conn[j][k]!=-1 && !typeLat[j] && typeLat[conn[j][k]])
		{
			tab_marquage[j] = 1;
			xf1 = j;
			xf2 = conn[j][bb[k]];
			xb = conn[j][k];
			q = solid_fraction_interpolation[xb][bb[k]];
			if(q<0.5)
			{
				lat.f_[xf1][bb[k]] = (1-2*q)*f_star[xf2][k]+2*q*f_star[xf1][k];
			}
			else
			{
				lat.f_[xf1][bb[k]] = (1-1./2*q)*f_star[xf1][bb[k]]+1./2*q*f_star[xf1][k];
			}
			//printf("Lattice %d, direction %d, valeur de lat.f LIM : %f\n", j,bb[k],lat.f_[xf1][bb[k]]);	
		}
	}
}

void quadratic_interpolation_method( int j, int Q, Lattice lat, double** f_star, int**conn, bool*  typeLat,  int* bb,double** solid_fraction_interpolation, double* tab_marquage, int cas)
{
	double q;
	int xf1, xf2,xf3,xb;
	for (int k=0;k<Q;k++)
	{
		//Si la lattice j est fluide, et le voisin est solide
		if ( conn[j][k]!=-1 && !typeLat[j] && typeLat[conn[j][k]])
		{
			tab_marquage[j] = 1;
			xf1 = j;
			xf2 = conn[j][bb[k]];
			xf3 = conn[xf2][bb[k]];
			xb = conn[j][k];
			q = solid_fraction_interpolation[xb][bb[k]];
			if(q<0.5)
			{
				lat.f_[xf1][bb[k]] = q*(1+2*q)*f_star[xf1][k]+(1-4*q*q)*f_star[xf2][k]-q*(1-2*q)*f_star[xf3][k];
			}
			else
			{
				lat.f_[xf1][bb[k]] = (2*q-1)/q*f_star[xf1][bb[k]]+1/(q*(2*q+1))*f_star[xf1][k]+(1-2*q)/(2*q+1)*f_star[xf2][bb[k]];
			}
			//printf("Lattice %d, direction %d, valeur de lat.f QIM : %f\n", j,bb[k],lat.f_[xf1][bb[k]]);
		}
	}
}

void multireflection_interpolation_method( int j, int Q, Lattice lat, double** f_star, int**conn, bool*  typeLat,  int* bb,double** solid_fraction_interpolation, double* tab_marquage, int cas, double** C, double* t, double* teq, double Fpc, double mu, double rho, double** Si)
{
	double q;
	//t et teq sont uniquement les moments d'ordre 3 (liés au flux d'énergie). Les autres moments sont considérés nuls
	
	for (int k=0;k<Q;k++)
	{
		t[k] =0;
		teq[k]=0;
	}
	int xf1, xf2,xf3,xb,k1,k2,k3,k4,k5;
	if(cas>4)
	{
		for (int k=0;k<Q;k++)
		{
			//Si la lattice j est fluide, et le voisin est solide
			if ( conn[j][k]!=-1 && !typeLat[j] && typeLat[conn[j][k]])
			{
				q = solid_fraction_interpolation[xb][bb[k]];
				//Création de Fpc
				t[5] = lat.m_[j][5]; //qx
				t[7] = lat.m_[j][7]; //qy
				teq[5] = lat.m0_[j][5];//qx_eq
				teq[7] = lat.m0_[j][7];//qy_eq
				for (int l=0;l<Q;l++)
				{
					Fpc+=4*(Si[7][7]-0.5)*(Si[k][k]-0.5)/(3*mu/rho*(1+q)*(1+q))*C[k][l]*(t[l]-teq[l]);
				}
				tab_marquage[j] = 1;
				xf1 = j;
				xf2 = conn[j][bb[k]];
				xf3 = conn[xf2][bb[k]];
				xb = conn[j][k];
				q = solid_fraction_interpolation[xb][bb[k]];
				k1 = (1-2*q-2*q*q)/((1+q)*(1+q));
				k2 = q*q/((1+q)*(1+q));

				lat.f_[xf1][bb[k]] = k1*f_star[xf2][k]+k2*f_star[xf3][k]-k1*f_star[xf1][bb[k]]-k2*f_star[xf2][bb[k]]+f_star[xf1][k]+Fpc;
			}
		}
	}
}

void central_interpolation_method( int j, int Q, Lattice lat, double** f_star, int**conn, bool*  typeLat,  int* bb,double** solid_fraction_interpolation, double* tab_marquage, int cas)
{
	double q;
	int xf1, xf2,xb;
	for (int k=0;k<Q;k++)
	{
		//Si la lattice j est fluide, et le voisin est solide
		if ( conn[j][k]!=-1 && !typeLat[j] && typeLat[conn[j][k]])
		{
			tab_marquage[j] = 1;
			xf1 = j;
			xf2 = conn[j][bb[k]];
			xb = conn[j][k];
			q = solid_fraction_interpolation[xb][bb[k]];		
			lat.f_[xf1][bb[k]] = (1-2*q)/(1+2*q)*f_star[xf2][k]-(1-2*q)/(1+2*q)*f_star[xf1][bb[k]] +f_star[xf1][k]  ;
			//printf("Lattice %d, direction %d, valeur de lat.f LIM : %f\n", j,bb[k],lat.f_[xf1][bb[k]]);
		}
	}
}
