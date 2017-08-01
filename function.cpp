// C++ Standard includes
#include <iostream> 
#include <cmath>
#include <ctime>
#include <stdio.h>

// Library includes
#include "tinyxml2.h"

// Local includes
#include "parser.h"
#include "domain.h"
#include "lattice.h"
#include "solutionExporterEulerian.h"
#include "function.h"

double pscal(double *a,double *b,  int D) //fonction pour le produit scalaire en D dimensions
{
	double* produit =new double[1];
	for ( int i=0;i<D;i++)
	{
		produit[0]+=a[i]*b[i];
	}
return produit[0];
delete produit;
}

void density ( int j, int Q, Lattice lat) //fonction pour la somme des fi pour calculer la masse volumique
{
	double* sigma = new double[1];
	for ( int i=0;i<Q;i++)
	{
		sigma[0]+=lat.f_[j][i];
	}
	lat.rho_[j] = sigma[0];
	delete sigma;
}
void velocity( int j, int D,  int Q, double** xi, Lattice lat)
{
	double* sigma = new double[1];
	for ( int k=0;k<D;k++)
	{
		sigma[0]= 0;
		for ( int i=0;i<Q;i++)
		{
			sigma[0]+=lat.f_[j][i]*xi[i][k];
		}
		lat.u_[j][k] = sigma[0]/lat.rho_[j];
	}
	delete sigma;
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


