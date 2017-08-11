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


double drag_force( int Q, int nx,int ny, double** f_star, double xi_r, bool* typeLat, int** conn, double** xi,  int* bb,  int* pos, Lattice lat)
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

double lift_force( int Q, int nx,int ny, double** f_star, double xi_r, bool* typeLat, int** conn, double** xi,  int* bb,  int* pos)
{
	double lforce =0; //drag force

	for (int j=pos[0];j<pos[1]+1;j++)
	{
		for (int k=1;k<Q;k++)
		{
			if(!typeLat[conn[j][k]] && typeLat[j])///Si le voisin est un noeud fluide et qu'on est sur un noeud solide
													   //alors on peut appliquer le calcul
			{
				lforce+=xi[k][1]*(f_star[j][k]+f_star[j][bb[k]]);
			}
		}
	}
		return lforce;
}

void vorticite(int Q, int nx,int ny, Lattice lat, int* cas,double dx, bool* typeLat, int** conn)
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
