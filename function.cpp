// C++ Standard includes
#include <iostream> 
#include <cmath>
#include <ctime>
#include <stdio.h>

// Local includes
#include "domain.h"
#include "lattice.h"
#include "solutionExporterEulerian.h"
#include "function.h"

# define M_PI  3.14159265358979323846

#define EULER 0.57721566 //Constante d'Euler
#define MAXIT 100 //Nombre maximum d'itérations
#define FPMIN 1.0e-30
#define EPS 6.0e-8 

double pscal(double *a,double *b,  int D, double sigma) //fonction pour le produit scalaire en D dimensions
{
	sigma=0;
	for ( int i=0;i<D;i++)
	{
		sigma+=a[i]*b[i];
	}
return sigma;
}

void density ( int j, int Q, Lattice lat, double sigma) //fonction pour la somme des fi pour calculer la masse volumique
{
	sigma = 0;
	for ( int i=0;i<Q;i++)
	{
		sigma+=lat.f_[j][i];
	}
	lat.rho_[j] = sigma;
}
void velocity( int j, int D,  int Q, double** xi, Lattice lat, double sigma)
{
	for ( int k=0;k<D;k++)
	{
		sigma= 0;
		for ( int i=0;i<Q;i++)
		{
			sigma+=lat.f_[j][i]*xi[i][k];
		}
		lat.u_[j][k] = sigma/lat.rho_[j];
	}
}
void domainCondition(int nx, int ny,  int* cas)
{
	for ( int j=0;j<nx*ny;j++)
	{
		if ((j%nx)==0 && j!=0 && j!=nx*ny-nx) //Lattices situées sur la frontière ouest (sauf les coins)
		{
			cas[j]=1;
		}
			else if ((j%nx)==nx-1 && j!=nx-1 && j!=nx*ny-1)  //Lattices situées sur la frontière est (sauf les coins) 
		{
			cas[j]=2;
		}
		else if (j>nx*ny-nx && j!=nx*ny-1) //Lattices situées sur la frontière nord (sauf les coins)
		{
			cas[j]=3;
		}
		else if (j<nx && j!=0 && j!=nx-1) //Lattices situées sur la frontière sud (sauf les coins)
		{
			cas[j]=4;
		}
		else if(j==0) //Coin Sud-Ouest
		{
			cas[j]=5;
		}
		else if(j==nx-1) //Coin Sud-Est
		{
			cas[j]=6;
		}
		else if(j==nx*ny-nx) //Coin Nord-Ouest
		{
			cas[j]=7;
		}
		else if(j==nx*ny-1) //Coin Nord-Est
		{
			cas[j]=8;
		}
	}
}

void solidCondition(int Q, int** conn,int nx, int ny, int* cas, bool* typeLat)
{
	for (int j=0;j<nx*ny;j++)
	{
		for (int k=1;k<Q;k++)
		{
			if (conn[j][k]!=-1 && typeLat[j] && !typeLat[conn[j][k]])
			{
				cas[j]=9;
				cas[conn[j][k]] = 10;
			}
		}
	}
}

//Fonction donnant la localisation de tous les noeuds de lattices (pour des noeuds off-grid)
void localisation(int nx, int ny, double dx, double** position)
{
	for ( int i=0;i<nx*ny;i++)
	{
		position[i][0] = 0.5*dx+dx*(i%nx);
		position[i][1] = 0.5*dx+dx*(int)(i/nx);
	}
}

//Création d'un cylindre carré avec sa position du milieu et son diamètre
void SquareCylinder(double abscisse, double ordonnee, double diametre, double** coin)
{
	coin[0][0] = abscisse - 0.5*diametre;
	coin[0][1] = ordonnee - 0.5*diametre;
	coin[1][0] = abscisse + 0.5*diametre;
	coin[1][1] = ordonnee - 0.5*diametre;
	coin[2][0] = abscisse - 0.5*diametre;
	coin[2][1] = ordonnee + 0.5*diametre;
	coin[3][0] = abscisse + 0.5*diametre;
	coin[3][1] = ordonnee + 0.5*diametre;
}


//Fonction permettant de créer un tableau de booléens, où 0 représente un noeud fluide et 1 un noeud solide pour un cylindre carré
void typeSquare( int N, double** coin, double** position, bool* typeLat)
{
	for ( int i=0;i<N;i++)
	{
		if(position[i][0]<coin[1][0] && position[i][0]>coin[0][0] && position[i][1]<coin[2][1] && position[i][1]>coin[0][1])
		{
			typeLat[i]=true;
		}
	}
}



void connectivite(int nx,int ny,  int Q, int** conn)
{
	for ( int j=0;j<nx*ny;j++)
	{
		conn[j]= new int[Q];
		conn[j][0] = j;
		conn[j][1] = j+1;
		conn[j][2] = j+nx;
		conn[j][3] = j-1;
		conn[j][4] = j-nx;
		conn[j][5] = j+nx+1;
		conn[j][6] = j+nx-1;
		conn[j][7] = j-nx-1;
		conn[j][8] = j-nx+1;
	}
	for ( int j=0;j<nx*ny;j++)
	{
		if((j%nx)==0)
		{
			conn[j][6] = -1;
			conn[j][3] = -1;
			conn[j][7] = -1;
		}
		if((j%nx)==nx-1)
		{
			conn[j][1] = -1;
			conn[j][5] = -1;
			conn[j][8] = -1;
		}
		for ( int k=0;k<Q;k++)
		{
			if (conn[j][k]<0 || conn[j][k]>=nx*ny)
			{
				conn[j][k] = -1;
			}
		}
	}
}

void pos_solide (bool* typeLat,  int* pos, int nx, int ny) //Donne les positions min et max des lattices solides
{
	 int pos1 = 0; //position de la première lattice solide
	 int pos2 = nx*ny-1; //position de la dernière lattice solide
	while(!typeLat[pos1] && pos1<nx*ny)
	{
		pos1++;

	}
	pos[0] = pos1;
	//Recherche du dernier noeud solide
	while(!typeLat[pos2] &&  pos2>0)
	{
		pos2--;
	}
	pos[1] = pos2;
}

void bounceback_neighbour( int* bb,  int Q)
{
	bb[0] = 0;
	int Q2 = 0.5*Q;
	for ( int k=1;k<Q;k++)
	{
		if ((k%(Q2))==1 || (k%(Q2))==2)
		{
			bb[k] = k+2;
		}
		else
		{
			bb[k] = k-2;
		}
	}

}


//Explication
// fi,eq = rho * omega_i * (1+1/(cs*cs)*pscal(ci,u) + 1/(2*cs^4) * Qi : uu)
// fi,eq = rho * omega_i * (1+1/(cs*cs)*pscal(ci,u) + 1/(2*cs^4) * Qeq)
// Tout ça en D2Q9

//Calcul de la double contraction Qi : uu
void Qi_equilibre(int Q, double cs, double** xi, double*** Qi)
{
	for (int k=0;k<Q;k++)
	{
		Qi[k][0][0] = xi[k][0]*xi[k][0] - cs*cs;
		Qi[k][0][1] = xi[k][0]*xi[k][1];
		Qi[k][1][0] = xi[k][1]*xi[k][0];
		Qi[k][1][1] = xi[k][1]*xi[k][1] - cs*cs;
	}
}
void fi_equilibre(int Q, double*** Qi, Lattice lat, double* omega_i, double**xi, double sigma, int D, int j, double cs)
{
	for (int k =0;k<Q;k++)
	{
		double Qeq = Qi[k][0][0] * lat.u_[j][0]*lat.u_[j][0] + Qi[k][0][1] * lat.u_[j][0]*lat.u_[j][1] + Qi[k][1][0] * lat.u_[j][1]*lat.u_[j][0] + Qi[k][1][1] * lat.u_[j][1]*lat.u_[j][1];
		lat.f0_[j][k] = lat.rho_[j] * omega_i[k] * (1+ 1/(cs*cs)*pscal(xi[k],lat.u_[j],D,sigma) + 1/(2*cs*cs*cs*cs) * Qeq);
	}
	
}

double Ei(double x) //Donne la valeur de Ei(x)
{
	int k;
	double fact, prev, sum,term;

	if(x <= 0.0) printf("Bad argument in ei");
	if(x < FPMIN) return log(x) + EULER;
	if(x <= -log(EPS))
	{
		sum = 0.0;
		fact = 1.0;
		for (k=1;k<=MAXIT;k++)
		{
			fact*=x/k;
			term=fact/k;
			sum+=term;
			if(term < EPS*sum) break;
		}
		if (k > MAXIT) printf("Series failes in ei");
		return sum+log(x)+EULER;
	}
	else
	{
		sum = 0.0;
		term = 1.0;
		for(k=1;k<=MAXIT;k++)
		{
			prev = term;
			term *=k/x;
			if (term < EPS) break;
			if (term < prev) sum+=term;
			else
			{
				sum -=prev;
				break;
			}
		}
		return exp(x)*(1.0+sum)/x;
	}
}



