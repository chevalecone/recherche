// C++ Standard includes
#include <iostream> 
#include <cmath>
#include <ctime>
#include <stdio.h>
using namespace std;

// Local includes
#include "domain.h"
#include "lattice.h"
#include "solutionExporterEulerian.h"
#include "function.h" 
#include "compute_quantities.h"

//Fonction permettant le calcul de la force de traînée selon Mei et al. (2002) sur le principe du momentum exchange
double drag_force( int Q, bool* typeLat, int** conn, double** xi,  int* bb,  int* pos, Lattice lat)
{
	double dforce =0; //drag force

	for (int j=pos[0];j<pos[1]+1;j++)
	{
		for (int k=1;k<Q;k++)
		{
			if(!typeLat[conn[j][k]] && typeLat[j])//Si le voisin est un noeud fluide et qu'on est sur un noeud solide
												  //alors on peut appliquer le calcul
			{
				dforce-=xi[k][0]*(lat.f_[conn[j][k]][k]+lat.f_[conn[j][k]][bb[k]]);
			}
		}
	}
		return dforce;
}

//Fonction permettant le calcul de la force de portance selon Mei et al. (2002) sur le principe du momentum exchange
double lift_force( int Q, bool* typeLat, int** conn, double** xi,  int* bb,  int* pos, Lattice lat)
{
	double lforce =0; //drag force

	for (int j=pos[0];j<pos[1]+1;j++)
	{
		for (int k=1;k<Q;k++)
		{
			if(!typeLat[conn[j][k]] && typeLat[j])///Si le voisin est un noeud fluide et qu'on est sur un noeud solide
													   //alors on peut appliquer le calcul
			{
				lforce+=xi[k][1]*(lat.f_[conn[j][k]][k]+lat.f_[conn[j][k]][bb[k]]);
			}
		}
	}
		return lforce;
}

//Fonction permettant le calcul de la vorticité
//Rappel : Vorticité : V = 0.5 * (d(ux)/dy - d(uy)/dx) . z (vorticité sur l'axe z)
void vorticite(int nx,int ny, Lattice lat,double dx, bool* typeLat, int** conn)
{
	for (int j=0;j<nx*ny;j++)
	{
		if(conn[j][1]!=-1 && conn[j][2]!=-1 && conn[j][3]!=-1 && conn[j][4]!=-1 && !typeLat[conn[j][1]] && !typeLat[conn[j][2]] && !typeLat[conn[j][3]] && !typeLat[conn[j][4]])
		{
			lat.vorticity_[j][0] = 0.5*((lat.u_[j+nx][1]-lat.u_[j-nx][1])-(lat.u_[j+1][0]-lat.u_[j-1][0]))/dx;
		}
		else
		{
			lat.vorticity_[j][0] = 0;
		}
	}
}

//Fonction permettant le calcul de la vorticité
//T = somme(U)/somme(u) où U = sqrt(u²+v²) et u est la composante selon l'axe principal de la vitesse //(Provient de Nabovati et Sousa)
//Une autre formule est disponible : T = somme(U/u)
double tortuosite(int nx,int ny, Lattice lat,double dx, bool* typeLat, int** conn)
{
	double v_mag =0, v_x = 0;
	for (int j=0;j<nx*ny;j++)
	{
		if(typeLat[j]==false)
		{
			v_mag+=sqrt(lat.u_[j][0]*lat.u_[j][0]+lat.u_[j][1]*lat.u_[j][1]);
			v_x+=abs(lat.u_[j][0]);
		}
	}
	return v_mag/v_x;
}

