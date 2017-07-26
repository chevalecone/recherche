// C++ Standard includes
#include <iostream> 
#include <cmath>
#include <ctime>
#include <stdio.h>
#include <algorithm>

// Library includes
#include "tinyxml2.h"

// Local includes
#include "parser.h"
#include "domain.h"
#include "lattice.h"
#include "solutionExporterEulerian.h"
#include "boundaryc.h"
#include "function.h"

void solidPropagation( int j, int nx, int ny , int cas, Lattice lat, double** f_star)
{
	switch (cas)
	{
		case 10:
		lat.f_[j][0] = f_star[j][0]; //direction nulle
		lat.f_[j+1][1] = f_star[j][1];
		lat.f_[j+nx][2] = f_star[j][2];
		lat.f_[j-nx][4] = f_star[j][4]; //directions verticales + horizontales
		lat.f_[j+nx+1][5] = f_star[j][5];
		lat.f_[j-nx+1][8] = f_star[j][8]; //directions obliques
		break;
		case 11:
		lat.f_[j][0] = f_star[j][0]; //direction nulle
		lat.f_[j+nx][2] = f_star[j][2];
		lat.f_[j-1][3] = f_star[j][3];
		lat.f_[j-nx][4] = f_star[j][4]; //directions verticales + horizontales
		lat.f_[j+nx-1][6] = f_star[j][6];
		lat.f_[j-nx-1][7] = f_star[j][7];
		break;
		case 12:
		lat.f_[j][0] = f_star[j][0]; //direction nulle
		lat.f_[j+1][1] = f_star[j][1];
		lat.f_[j-1][3] = f_star[j][3];
		lat.f_[j-nx][4] = f_star[j][4]; //directions verticales + horizontales
		lat.f_[j-nx-1][7] = f_star[j][7];
		lat.f_[j-nx+1][8] = f_star[j][8]; //directions obliques
		break;
		case 13:
		lat.f_[j][0] = f_star[j][0]; //direction nulle
		lat.f_[j+1][1] = f_star[j][1];
		lat.f_[j+nx][2] = f_star[j][2];
		lat.f_[j-1][3] = f_star[j][3];
		lat.f_[j+nx+1][5] = f_star[j][5];
		lat.f_[j+nx-1][6] = f_star[j][6];
		break;
		case 14:
		lat.f_[j][0] = f_star[j][0]; //direction nulle
		lat.f_[j+1][1] = f_star[j][1];
		lat.f_[j+nx][2] = f_star[j][2];
		lat.f_[j-1][3] = f_star[j][3];
		lat.f_[j-nx][4] = f_star[j][4]; //directions verticales + horizontales
		lat.f_[j+nx+1][5] = f_star[j][5];
		lat.f_[j+nx-1][6] = f_star[j][6];
		lat.f_[j-nx+1][8] = f_star[j][8]; //directions obliques
		break;
		case 15:
		lat.f_[j][0] = f_star[j][0]; //direction nulle
		lat.f_[j+1][1] = f_star[j][1];
		lat.f_[j+nx][2] = f_star[j][2];
		lat.f_[j-1][3] = f_star[j][3];
		lat.f_[j-nx][4] = f_star[j][4]; //directions verticales + horizontales
		lat.f_[j+nx+1][5] = f_star[j][5];
		lat.f_[j+nx-1][6] = f_star[j][6];
		lat.f_[j-nx-1][7] = f_star[j][7];
		break;
		case 16:
		lat.f_[j][0] = f_star[j][0]; //direction nulle
		lat.f_[j+1][1] = f_star[j][1];
		lat.f_[j+nx][2] = f_star[j][2];
		lat.f_[j-1][3] = f_star[j][3];
		lat.f_[j-nx][4] = f_star[j][4]; //directions verticales + horizontales
		lat.f_[j+nx+1][5] = f_star[j][5];
		lat.f_[j-nx-1][7] = f_star[j][7];
		lat.f_[j-nx+1][8] = f_star[j][8]; //directions obliques
		break;
		case 17:
		lat.f_[j][0] = f_star[j][0]; //direction nulle
		lat.f_[j+1][1] = f_star[j][1];
		lat.f_[j+nx][2] = f_star[j][2];
		lat.f_[j-1][3] = f_star[j][3];
		lat.f_[j-nx][4] = f_star[j][4]; //directions verticales + horizontales
		lat.f_[j+nx-1][6] = f_star[j][6];
		lat.f_[j-nx-1][7] = f_star[j][7];
		lat.f_[j-nx+1][8] = f_star[j][8]; //directions obliques
		break;
		case 18:
		lat.f_[j][0] = f_star[j][0]; //direction nulle
		lat.f_[j+1][1] = f_star[j][1];
		lat.f_[j+nx][2] = f_star[j][2];
		lat.f_[j-1][3] = f_star[j][3];
		lat.f_[j+nx+1][5] = f_star[j][5];
		lat.f_[j+nx-1][6] = f_star[j][6];
		lat.f_[j-nx+1][8] = f_star[j][8]; //directions obliques
		break;
		case 19:
		lat.f_[j][0] = f_star[j][0]; //direction nulle
		lat.f_[j+1][1] = f_star[j][1];
		lat.f_[j+nx][2] = f_star[j][2];
		lat.f_[j-1][3] = f_star[j][3];
		lat.f_[j+nx+1][5] = f_star[j][5];
		lat.f_[j+nx-1][6] = f_star[j][6];
		lat.f_[j-nx-1][7] = f_star[j][7];
		break;
		case 20:
		lat.f_[j][0] = f_star[j][0]; //direction nulle
		lat.f_[j+nx][2] = f_star[j][2];
		lat.f_[j-1][3] = f_star[j][3];
		lat.f_[j-nx][4] = f_star[j][4]; //directions verticales + horizontales
		lat.f_[j+nx+1][5] = f_star[j][5];
		lat.f_[j+nx-1][6] = f_star[j][6];
		lat.f_[j-nx-1][7] = f_star[j][7];
		break;
		case 21:
		lat.f_[j][0] = f_star[j][0]; //direction nulle
		lat.f_[j+nx][2] = f_star[j][2];
		lat.f_[j-1][3] = f_star[j][3];
		lat.f_[j-nx][4] = f_star[j][4]; //directions verticales + horizontales
		lat.f_[j+nx-1][6] = f_star[j][6];
		lat.f_[j-nx-1][7] = f_star[j][7];
		lat.f_[j-nx+1][8] = f_star[j][8]; //directions obliques
		break;
		case 22:
		lat.f_[j][0] = f_star[j][0]; //direction nulle
		lat.f_[j+1][1] = f_star[j][1];
		lat.f_[j-1][3] = f_star[j][3];
		lat.f_[j-nx][4] = f_star[j][4]; //directions verticales + horizontales
		lat.f_[j+nx+1][5] = f_star[j][5];
		lat.f_[j+nx-1][6] = f_star[j][6];
		lat.f_[j-nx+1][8] = f_star[j][8]; //directions obliques
		break;
		case 23:
		lat.f_[j][0] = f_star[j][0]; //direction nulle
		lat.f_[j+1][1] = f_star[j][1];
		lat.f_[j-1][3] = f_star[j][3];
		lat.f_[j-nx][4] = f_star[j][4]; //directions verticales + horizontales
		lat.f_[j+nx+1][5] = f_star[j][5];
		lat.f_[j-nx-1][7] = f_star[j][7];
		lat.f_[j-nx+1][8] = f_star[j][8]; //directions obliques
		break;
		case 24:
		lat.f_[j][0] = f_star[j][0]; //direction nulle
		lat.f_[j+1][1] = f_star[j][1];
		lat.f_[j+nx][2] = f_star[j][2];
		lat.f_[j-nx][4] = f_star[j][4]; //directions verticales + horizontales
		lat.f_[j+nx+1][5] = f_star[j][5];
		lat.f_[j-nx-1][7] = f_star[j][7];
		lat.f_[j-nx+1][8] = f_star[j][8]; //directions obliques
		break;
		case 25:
		lat.f_[j][0] = f_star[j][0]; //direction nulle
		lat.f_[j+1][1] = f_star[j][1];
		lat.f_[j+nx][2] = f_star[j][2];
		lat.f_[j-nx][4] = f_star[j][4]; //directions verticales + horizontales
		lat.f_[j+nx+1][5] = f_star[j][5];
		lat.f_[j+nx-1][6] = f_star[j][6];
		lat.f_[j-nx+1][8] = f_star[j][8]; //directions obliques
		break;
	}
}
void propagation( int const& j, Lattice lat, double** f_star,int nx, int ny, int* cas, bool* typeLat, int** conn)
{
	for (int k=0;k<9;k++)
	{
		if(conn[j][k]!=-1 && !typeLat[j] && !typeLat[conn[j][k]])
		{
			lat.f_[conn[j][k]][k] = f_star[j][k];
		}
	}
}
void periodic_NS_BC( int j,int nx, int ny,  int cas, Lattice lat, double** f_star)
{
	switch (cas)
	{
		case 3:
		lat.f_[j][4] = f_star[j%nx][4]; //Condition périodique en entrée
		lat.f_[j][7] = f_star[j%nx+1][7];
		lat.f_[j][8] = f_star[j%nx-1][8];
		break;
		case 4:
		lat.f_[j][6] = f_star[nx*ny-nx+j+1][6]; //Condition périodique en sortie
		lat.f_[j][2] = f_star[nx*ny-nx+j][2];
		lat.f_[j][5] = f_star[nx*ny-nx+j-1][5];
		break;
		case 5:
		lat.f_[j][2] = f_star[nx*ny-nx][2]; //Conditions périodique en sortie
		lat.f_[j][6] = f_star[nx*ny-nx+1][6];
		lat.f_[j][5] = f_star[nx*ny-1][5];
		break;
		case 6:
		lat.f_[j][2] = f_star[nx*ny-1][2]; //Condition périodique en sortie
		lat.f_[j][5] = f_star[nx*ny-2][5];
		lat.f_[j][8] = f_star[nx*ny-nx][8];
		break;
		case 7:
		lat.f_[j][4] = f_star[0][4];//Condition périodique en entrée
		lat.f_[j][7] = f_star[1][7];
		lat.f_[j][8] = f_star[nx-1][8];
		break;
		case 8:
		lat.f_[j][4] = f_star[nx-1][4]; //Condition périodique en entrée
		lat.f_[j][8] = f_star[nx-2][8];
		lat.f_[j][7] = f_star[0][7];
		break;
	}
}
void periodic_WE_BC( int j,int nx, int ny,  int cas, Lattice lat, double** f_star)
{
	switch (cas)
	{
		case 1: //ouest
		lat.f_[j][1] = f_star[j+nx-1][1]; //Condition périodique en entrée
		lat.f_[j][5] = f_star[j-1][5];
		lat.f_[j][8] = f_star[j+2*nx-1][8];
		break;
		case 2: //est
		lat.f_[j][3] = f_star[j-nx+1][3]; //Condition périodique en sortie
		lat.f_[j][6] = f_star[j-2*nx+1][6];
		lat.f_[j][7] = f_star[j+1][7];
		break;
		case 5://S-O
		lat.f_[j][1] = f_star[nx-1][1]; //Conditions périodique en sortie
		lat.f_[j][8] = f_star[2*nx-1][8];
		lat.f_[j][5] = f_star[nx*ny-1][5];
		break;
		case 6://S-E
		lat.f_[j][3] = f_star[0][3]; //Condition périodique en sortie
		lat.f_[j][7] = f_star[nx][7];
		lat.f_[j][6] = f_star[nx*ny-nx][6];

		break;
		case 7://N-O
		lat.f_[j][1] = f_star[nx*ny-1][1];//Condition périodique en entrée
		lat.f_[j][5] = f_star[j-1][5];
		lat.f_[j][8] = f_star[nx-1][8];
		break;
		case 8://N-E
		lat.f_[j][3] = f_star[j-nx+1][3]; //Condition périodique en entrée
		lat.f_[j][6] = f_star[j-2*nx+1][6];
		lat.f_[j][7] = f_star[0][7];
		break;
	}
}


//HWBB pour le solide carré
void bounceback_solid_BC(int nx, int const& j, Lattice lat, double** f_star, int** const& conn, bool*  typeLat,  int* const& bb, double& nombre, int* pos) //Cas spéciaux pour le solide (coins et voisins)
{
	if(j>pos[0]-nx-2 && j<pos[1]+nx+2)
	{
		for (int k=0;k<9;k++)
		{	
			if ( conn[j][k]!=-1 && !typeLat[j] && typeLat[conn[j][k]])//Si la lattice est fluide et son voisin est un solide, alors il y a BB
			{
				lat.f_[j][bb[k]] = f_star[j][k];
				nombre++;
			}
		}
	}
}

/*void pression_in_BC( int j,  int cas, Lattice lat, double xi_r, double rho_in)
{
	double u_x;
	switch (cas)
	{
		case 1: //Ouest
		u_x = xi_r-xi_r/rho_in*(lat.f_[j][0]+lat.f_[j][2]+lat.f_[j][4]+2*(lat.f_[j][3]+lat.f_[j][6]+lat.f_[j][7]));
		lat.f_[j][1] = lat.f_[j][3]+2/(3*xi_r)*rho_in*u_x;
		lat.f_[j][5] = lat.f_[j][7]-0.5*(lat.f_[j][2]-lat.f_[j][4])+1/(6*xi_r)*rho_in*u_x ;//+ 1/(2*xi_r)*rho_in*lat.u_[j][1];
		lat.f_[j][8] = lat.f_[j][6]+0.5*(lat.f_[j][2]-lat.f_[j][4])+1/(6*xi_r)*rho_in*u_x ;//- 1/(2*xi_r)*rho_in*lat.u_[j][1];
		break;
		case 5: //Coin Sud-Ouest
		lat.f_[j][1] = lat.f_[j][3];
		lat.f_[j][8] = 0.5*(rho_in-(lat.f_[j][0]+lat.f_[j][1]+lat.f_[j][2]+lat.f_[j][3]+lat.f_[j][4]+lat.f_[j][5]+lat.f_[j][7]));
		lat.f_[j][6] = 0.5*(rho_in-(lat.f_[j][0]+lat.f_[j][1]+lat.f_[j][2]+lat.f_[j][3]+lat.f_[j][4]+lat.f_[j][5]+lat.f_[j][7]));
		/*u_x = xi_r-xi_r/rho_in*(lat.f_[j][0]+lat.f_[j][2]+lat.f_[j][4]+2*(lat.f_[j][3]+lat.f_[j][6]+lat.f_[j][7]));
		lat.f_[j][1] = lat.f_[j][3]+2/(3*xi_r)*rho_in*u_x;
		lat.f_[j][5] = lat.f_[j][7]-0.5*(lat.f_[j][2]-lat.f_[j][4])+1/(6*xi_r)*rho_in*u_x ;//+ 1/(2*xi_r)*rho_in*lat.u_[j][1];
		lat.f_[j][8] = lat.f_[j][6]+0.5*(lat.f_[j][2]-lat.f_[j][4])+1/(6*xi_r)*rho_in*u_x ;//- 1/(2*xi_r)*rho_in*lat.u_[j][1];
		break;*/

		/*case 7: //Coin Nord-Ouest
		lat.f_[j][1] = lat.f_[j][3];
		lat.f_[j][7] = 0.5*(rho_in-(lat.f_[j][0]+lat.f_[j][1]+lat.f_[j][2]+lat.f_[j][3]+lat.f_[j][4]+lat.f_[j][6]+lat.f_[j][8]));
		lat.f_[j][5] =0.5*(rho_in-(lat.f_[j][0]+lat.f_[j][1]+lat.f_[j][2]+lat.f_[j][3]+lat.f_[j][4]+lat.f_[j][6]+lat.f_[j][8]));
		/*u_x = xi_r-xi_r/rho_in*(lat.f_[j][0]+lat.f_[j][2]+lat.f_[j][4]+2*(lat.f_[j][3]+lat.f_[j][6]+lat.f_[j][7]));
		lat.f_[j][1] = lat.f_[j][3]+2/(3*xi_r)*rho_in*u_x;
		lat.f_[j][5] = lat.f_[j][7]-0.5*(lat.f_[j][2]-lat.f_[j][4])+1/(6*xi_r)*rho_in*u_x ;//+ 1/(2*xi_r)*rho_in*lat.u_[j][1];
		lat.f_[j][8] = lat.f_[j][6]+0.5*(lat.f_[j][2]-lat.f_[j][4])+1/(6*xi_r)*rho_in*u_x ;//- 1/(2*xi_r)*rho_in*lat.u_[j][1];*/
		/*break;
	}
}*/

	void pression_in_BC( int j,  int cas, Lattice lat, double xi_r, double rho_in)
{
	double u_x;
	switch (cas)
	{
		case 1: //Ouest
		u_x = xi_r-xi_r/rho_in*(lat.f_[j][0]+lat.f_[j][2]+lat.f_[j][4]+2*(lat.f_[j][3]+lat.f_[j][6]+lat.f_[j][7]));
		lat.f_[j][1] = lat.f_[j][3]+2/(3*xi_r)*rho_in*u_x;
		lat.f_[j][5] = lat.f_[j][7]-0.5*(lat.f_[j][2]-lat.f_[j][4])+1/(6*xi_r)*rho_in*u_x + 1/(2*xi_r)*rho_in*lat.u_[j][1];
		lat.f_[j][8] = lat.f_[j][6]+0.5*(lat.f_[j][2]-lat.f_[j][4])+1/(6*xi_r)*rho_in*u_x - 1/(2*xi_r)*rho_in*lat.u_[j][1];
		break;
		case 5: //Coin Sud-Ouest
		lat.f_[j][1] = lat.f_[j][3];
		lat.f_[j][5] = lat.f_[j][7];
		lat.f_[j][8] = lat.f_[j][6] - 0.5 * (lat.f_[j][2] - lat.f_[j][4]) + rho_in/(2*xi_r)*(lat.u_[j][0] - lat.u_[j][1]);
		break;

		case 7: //Coin Nord-Ouest
		lat.f_[j][1] = lat.f_[j][3];
		lat.f_[j][8] = lat.f_[j][6];
		lat.f_[j][5] = lat.f_[j][7] - 0.5 * (lat.f_[j][2] - lat.f_[j][4]) + rho_in/(2*xi_r)*(lat.u_[j][0] + lat.u_[j][1]);
		break;
	}
}
void pression_out_BC( int j,  int cas, Lattice lat, double xi_r,double rho_out)
{
	double u_x;
	switch (cas)
	{
		case 2: //Est
		u_x = xi_r/rho_out*(lat.f_[j][0]+lat.f_[j][2]+lat.f_[j][4]+2*(lat.f_[j][1]+lat.f_[j][5]+lat.f_[j][8]))-xi_r;
		lat.f_[j][3] = lat.f_[j][1] - 2/(3*xi_r)*rho_out*u_x;
		lat.f_[j][7] = lat.f_[j][5] + 0.5*(lat.f_[j][2]-lat.f_[j][4]) - 1/(6*xi_r)*rho_out*u_x ;//- 1/(2*xi_r)*rho_out*lat.u_[j][1];
		lat.f_[j][6] = lat.f_[j][8] - 0.5*(lat.f_[j][2]-lat.f_[j][4]) - 1/(6*xi_r)*rho_out*u_x ;//+ 1/(2*xi_r)*rho_out*lat.u_[j][1];
		break;
		case 6: //Coin Sud-Est
		lat.f_[j][3] = lat.f_[j][1];
		lat.f_[j][6] = lat.f_[j][8];
		lat.f_[j][7] = lat.f_[j][5] + 0.5 * (lat.f_[j][2] - lat.f_[j][4]) + rho_out/(2*xi_r)*(lat.u_[j][0] /*+ lat.u_[j][1]*/);
		break;
		case 8: //Coin Nord-Estj
		lat.f_[j][3] = lat.f_[j][1];
		lat.f_[j][7] = lat.f_[j][5];
		lat.f_[j][6] = lat.f_[j][8] + 0.5 * (lat.f_[j][2] - lat.f_[j][4]) + rho_out/(2*xi_r)*(lat.u_[j][0] /*- lat.u_[j][1]*/);
		break;
	}
}
void vitesse_in_BC( int j,int nx,int cas, Lattice lat, double xi_r, double** v_in)
{
	switch (cas)
	{
		case 1: //Ouest
		lat.rho_[j] = xi_r/(xi_r-v_in[j/nx][0])*(lat.f_[j][0]+lat.f_[j][2]+lat.f_[j][4]+2*(lat.f_[j][3]+lat.f_[j][6]+lat.f_[j][7]));
		lat.f_[j][1] = lat.f_[j][3] + 2/(3*xi_r)*lat.rho_[j]*v_in[j/nx][0];
		lat.f_[j][5] = lat.f_[j][7] - 0.5*(lat.f_[j][2]-lat.f_[j][4]) + 1/(6*xi_r)*lat.rho_[j]*v_in[j/nx][0] ;//+ 0.5*lat.rho_[j]/xi_r*lat.u_[j][1]; 
		lat.f_[j][8] = lat.f_[j][6] + 0.5*(lat.f_[j][2]-lat.f_[j][4]) + 1/(6*xi_r)*lat.rho_[j]*v_in[j/nx][0] ;//- 0.5*lat.rho_[j]/xi_r*lat.u_[j][1]; 
		break;
		case 5: //Coin Sud-Ouest
		/*lat.f_[j][1] = lat.f_[j][3];
		lat.f_[j][8] = 0.5*(lat.rho_[j+nx]-(lat.f_[j][0]+lat.f_[j][1]+lat.f_[j][2]+lat.f_[j][3]+lat.f_[j][4]+lat.f_[j][5]+lat.f_[j][7]));
		lat.f_[j][6] = 0.5*(lat.rho_[j+nx]-(lat.f_[j][0]+lat.f_[j][1]+lat.f_[j][2]+lat.f_[j][3]+lat.f_[j][4]+lat.f_[j][5]+lat.f_[j][7]));*/
		lat.rho_[j] = xi_r/(xi_r-v_in[j/nx][0])*(lat.f_[j][0]+lat.f_[j][2]+lat.f_[j][4]+2*(lat.f_[j][3]+lat.f_[j][6]+lat.f_[j][7]));
		lat.f_[j][1] = lat.f_[j][3] + 2/(3*xi_r)*lat.rho_[j]*v_in[j/nx][0];
		lat.f_[j][5] = lat.f_[j][7] - 0.5*(lat.f_[j][2]-lat.f_[j][4]) + 1/(6*xi_r)*lat.rho_[j]*v_in[j/nx][0] ;//+ 0.5*lat.rho_[j]/xi_r*lat.u_[j][1];
		lat.f_[j][8] = lat.f_[j][6] + 0.5*(lat.f_[j][2]-lat.f_[j][4]) + 1/(6*xi_r)*lat.rho_[j]*v_in[j/nx][0] ;//- 0.5*lat.rho_[j]/xi_r*lat.u_[j][1] ;
		break;
		case 7: //Coin Nord-Ouest
		/*lat.f_[j][1] = lat.f_[j][3];
		lat.f_[j][7] = 0.5*(lat.rho_[j-nx]-(lat.f_[j][0]+lat.f_[j][1]+lat.f_[j][2]+lat.f_[j][3]+lat.f_[j][4]+lat.f_[j][6]+lat.f_[j][8]));
		lat.f_[j][5] =0.5*(lat.rho_[j-nx]-(lat.f_[j][0]+lat.f_[j][1]+lat.f_[j][2]+lat.f_[j][3]+lat.f_[j][4]+lat.f_[j][6]+lat.f_[j][8]));*/
		lat.rho_[j] = xi_r/(xi_r-v_in[j/nx][0])*(lat.f_[j][0]+lat.f_[j][2]+lat.f_[j][4]+2*(lat.f_[j][3]+lat.f_[j][6]+lat.f_[j][7]));
		lat.f_[j][1] = lat.f_[j][3] + 2/(3*xi_r)*lat.rho_[j]*v_in[j/nx][0];
		lat.f_[j][5] = lat.f_[j][7] - 0.5*(lat.f_[j][2]-lat.f_[j][4]) + 1/(6*xi_r)*lat.rho_[j]*v_in[j/nx][0] ;//+ 0.5*lat.rho_[j]/xi_r*lat.u_[j][1];
		lat.f_[j][8] = lat.f_[j][6] + 0.5*(lat.f_[j][2]-lat.f_[j][4]) + 1/(6*xi_r)*lat.rho_[j]*v_in[j/nx][0] ;//- 0.5*lat.rho_[j]/xi_r*lat.u_[j][1];

		break;
	}
}

/*void pression_out_BC( int j,  int cas, Lattice lat, double xi_r,double rho_out)
{
	double u_x;
	switch (cas)
	{
		case 2: //Est
		u_x = xi_r/rho_out*(lat.f_[j][0]+lat.f_[j][2]+lat.f_[j][4]+2*(lat.f_[j][1]+lat.f_[j][5]+lat.f_[j][8]))-xi_r;
		lat.f_[j][3] = lat.f_[j][1] - 2/(3*xi_r)*rho_out*u_x;
		lat.f_[j][7] = lat.f_[j][5] + 0.5*(lat.f_[j][2]-lat.f_[j][4]) - 1/(6*xi_r)*rho_out*u_x ;//- 1/(2*xi_r)*rho_out*lat.u_[j][1];
		lat.f_[j][6] = lat.f_[j][8] - 0.5*(lat.f_[j][2]-lat.f_[j][4]) - 1/(6*xi_r)*rho_out*u_x ;//+ 1/(2*xi_r)*rho_out*lat.u_[j][1];
		break;
		case 6: //Coin Sud-Est
		lat.f_[j][3] = lat.f_[j][1];
		lat.f_[j][5] = 0.5*(rho_out-(lat.f_[j][0]+lat.f_[j][1]+lat.f_[j][2]+lat.f_[j][3]+lat.f_[j][4]+lat.f_[j][6]+lat.f_[j][8]));
		lat.f_[j][7] = 0.5*(rho_out-(lat.f_[j][0]+lat.f_[j][1]+lat.f_[j][2]+lat.f_[j][3]+lat.f_[j][4]+lat.f_[j][6]+lat.f_[j][8]));
		/*u_x = xi_r/rho_out*(lat.f_[j][0]+lat.f_[j][2]+lat.f_[j][4]+2*(lat.f_[j][1]+lat.f_[j][5]+lat.f_[j][8]))-xi_r;
		lat.f_[j][3] = lat.f_[j][1] - 2/(3*xi_r)*rho_out*u_x;
		lat.f_[j][7] = lat.f_[j][5] + 0.5*(lat.f_[j][2]-lat.f_[j][4]) - 1/(6*xi_r)*rho_out*u_x ;//- 1/(2*xi_r)*rho_out*lat.u_[j][1];
		lat.f_[j][6] = lat.f_[j][8] - 0.5*(lat.f_[j][2]-lat.f_[j][4]) - 1/(6*xi_r)*rho_out*u_x ;//+ 1/(2*xi_r)*rho_out*lat.u_[j][1];*/
		/*break;
		case 8: //Coin Nord-Est
		lat.f_[j][3] = lat.f_[j][1];
		lat.f_[j][6] = 0.5*(rho_out-(lat.f_[j][0]+lat.f_[j][1]+lat.f_[j][2]+lat.f_[j][3]+lat.f_[j][4]+lat.f_[j][5]+lat.f_[j][7]));
		lat.f_[j][8] = 0.5*(rho_out-(lat.f_[j][0]+lat.f_[j][1]+lat.f_[j][2]+lat.f_[j][3]+lat.f_[j][4]+lat.f_[j][5]+lat.f_[j][7]));
		/*u_x = xi_r/rho_out*(lat.f_[j][0]+lat.f_[j][2]+lat.f_[j][4]+2*(lat.f_[j][1]+lat.f_[j][5]+lat.f_[j][8]))-xi_r;
		lat.f_[j][3] = lat.f_[j][1] - 2/(3*xi_r)*rho_out*u_x;
		lat.f_[j][7] = lat.f_[j][5] + 0.5*(lat.f_[j][2]-lat.f_[j][4]) - 1/(6*xi_r)*rho_out*u_x ;//- 1/(2*xi_r)*rho_out*lat.u_[j][1];
		lat.f_[j][6] = lat.f_[j][8] - 0.5*(lat.f_[j][2]-lat.f_[j][4]) - 1/(6*xi_r)*rho_out*u_x ;//+ 1/(2*xi_r)*rho_out*lat.u_[j][1];*/
		/*break;
	}
}*/
void vitesse_out_BC( int j,int nx,int cas, Lattice lat, double xi_r,double** v_out)
{
	switch (cas)
	{
		case 2: //Est 
		lat.rho_[j] = xi_r/(xi_r+v_out[j/nx][0])*(lat.f_[j][0]+lat.f_[j][2]+lat.f_[j][4]+2*(lat.f_[j][1]+lat.f_[j][5]+lat.f_[j][8]));
		lat.f_[j][3] = lat.f_[j][1] - 2/(3*xi_r)*lat.rho_[j]*lat.u_[j][0];
		lat.f_[j][7] = lat.f_[j][5] + 0.5*(lat.f_[j][2]-lat.f_[j][4]) - 1/(6*xi_r)*lat.rho_[j]*v_out[j/nx][0] ;//- 1/(2*xi_r)*lat.rho_[j]*lat.u_[j][1];
		lat.f_[j][6] = lat.f_[j][8] - 0.5*(lat.f_[j][2]-lat.f_[j][4]) - 1/(6*xi_r)*lat.rho_[j]*v_out[j/nx][0] ;//+ 1/(2*xi_r)*lat.rho_[j]*lat.u_[j][1]; 
		break;
		case 6: //Coin Sud-Est
		lat.rho_[j] = xi_r/(xi_r+v_out[j/nx][0])*(lat.f_[j][0]+lat.f_[j][2]+lat.f_[j][4]+2*(lat.f_[j][1]+lat.f_[j][5]+lat.f_[j][8]));
		lat.f_[j][3] = lat.f_[j][1] - 2/(3*xi_r)*lat.rho_[j]*lat.u_[j][0];
		lat.f_[j][7] = lat.f_[j][5] + 0.5*(lat.f_[j][2]-lat.f_[j][4]) - 1/(6*xi_r)*lat.rho_[j]*v_out[j/nx][0] ;//- 1/(2*xi_r)*lat.rho_[j]*lat.u_[j][1];
		lat.f_[j][6] = lat.f_[j][8] - 0.5*(lat.f_[j][2]-lat.f_[j][4]) - 1/(6*xi_r)*lat.rho_[j]*v_out[j/nx][0] ;//+ 1/(2*xi_r)*lat.rho_[j]*lat.u_[j][1]; 
		break;
		case 8: //Coin Nord-Est
		lat.rho_[j] = xi_r/(xi_r+v_out[j/nx][0])*(lat.f_[j][0]+lat.f_[j][2]+lat.f_[j][4]+2*(lat.f_[j][1]+lat.f_[j][5]+lat.f_[j][8]));
		lat.f_[j][3] = lat.f_[j][1] - 2/(3*xi_r)*lat.rho_[j]*lat.u_[j][0];
		lat.f_[j][7] = lat.f_[j][5] + 0.5*(lat.f_[j][2]-lat.f_[j][4]) - 1/(6*xi_r)*lat.rho_[j]*v_out[j/nx][0] ;//- 1/(2*xi_r)*lat.rho_[j]*lat.u_[j][1];
		lat.f_[j][6] = lat.f_[j][8] - 0.5*(lat.f_[j][2]-lat.f_[j][4]) - 1/(6*xi_r)*lat.rho_[j]*v_out[j/nx][0] ;//+ 1/(2*xi_r)*lat.rho_[j]*lat.u_[j][1]; 
		break;
	}
}
	

void driven_cavity_nord( int j,  int cas, Lattice lat, double xi_r,double v_e)
{
	switch (cas)
	{
		case 3: //Nord
		lat.rho_[j] = 1/(1+lat.u_[j][1]/xi_r)*(lat.f_[j][0]+lat.f_[j][1]+lat.f_[j][3]+2*(lat.f_[j][2]+lat.f_[j][5]+lat.f_[j][6]));
		lat.f_[j][4] = lat.f_[j][2]-2/(3*xi_r)*lat.rho_[j]*lat.u_[j][1];
		lat.f_[j][7] = lat.f_[j][5]+0.5*(lat.f_[j][1]-lat.f_[j][3])-0.5*lat.rho_[j]/xi_r*v_e-lat.rho_[j]/(6*xi_r)*lat.u_[j][1];
		lat.f_[j][8] = lat.f_[j][6]-0.5*(lat.f_[j][1]-lat.f_[j][3])+0.5*lat.rho_[j]/xi_r*v_e-lat.rho_[j]/(6*xi_r)*lat.u_[j][1];
		break;
		case 7: //Coin Nord-Ouest
		lat.rho_[j] = 1/(1+lat.u_[j][1]/xi_r)*(lat.f_[j][0]+lat.f_[j][1]+lat.f_[j][3]+2*(lat.f_[j][2]+lat.f_[j][5]+lat.f_[j][6]));
		lat.f_[j][4] = lat.f_[j][2]-2/(3*xi_r)*lat.rho_[j]*lat.u_[j][1];
		lat.f_[j][7] = lat.f_[j][5]+0.5*(lat.f_[j][1]-lat.f_[j][3])-0.5*lat.rho_[j]/xi_r*v_e-lat.rho_[j]/(6*xi_r)*lat.u_[j][1];
		lat.f_[j][8] = lat.f_[j][6]-0.5*(lat.f_[j][1]-lat.f_[j][3])+0.5*lat.rho_[j]/xi_r*v_e-lat.rho_[j]/(6*xi_r)*lat.u_[j][1];
		break;
		case 8: //Coin Nord-Est
		lat.rho_[j] = 1/(1+lat.u_[j][1]/xi_r)*(lat.f_[j][0]+lat.f_[j][1]+lat.f_[j][3]+2*(lat.f_[j][2]+lat.f_[j][5]+lat.f_[j][6]));
		lat.f_[j][4] = lat.f_[j][2]-2/(3*xi_r)*lat.rho_[j]*lat.u_[j][1];
		lat.f_[j][7] = lat.f_[j][5]+0.5*(lat.f_[j][1]-lat.f_[j][3])-0.5*lat.rho_[j]/xi_r*v_e-lat.rho_[j]/(6*xi_r)*lat.u_[j][1];
		lat.f_[j][8] = lat.f_[j][6]-0.5*(lat.f_[j][1]-lat.f_[j][3])+0.5*lat.rho_[j]/xi_r*v_e-lat.rho_[j]/(6*xi_r)*lat.u_[j][1];
		break;
	}
}
void bounceback_W_BC( int j, int cas, Lattice lat, double** f_star)
{
	switch (cas)
	{
		case 1: //O
		lat.f_[j][1] = f_star[j][3];//bounce-back sur le mur ouest
		lat.f_[j][8] = f_star[j][6];
		lat.f_[j][5] = f_star[j][7];
		break;
		case 5://S-O
		lat.f_[j][1] = f_star[j][3];
		lat.f_[j][8] = f_star[j][6];//bounce-back sur la frontière ouest
		//lat.f_[j][5] = f_star[j][7];
		break;
		case 7://N-O
		lat.f_[j][1] = f_star[j][3];//Bounce-back sur la frontière ouest
		lat.f_[j][5] = f_star[j][7];
		//lat.f_[j][8] = f_star[j][6];
		break;
	}
}
void bounceback_E_BC( int j, int cas, Lattice lat, double** f_star)
{
	switch (cas)
	{
		case 2://E
		lat.f_[j][3] = f_star[j][1];//bounce-back sur le mur est
		lat.f_[j][6] = f_star[j][8];
		lat.f_[j][7] = f_star[j][5];
		break;
		case 6://S-E
		lat.f_[j][3] = f_star[j][1]; //bounce-back sur la frontière est
		//lat.f_[j][6] = f_star[j][8];
		lat.f_[j][7] = f_star[j][5];
		break;
		case 8://N-E
		lat.f_[j][3] = f_star[j][1]; //Bounce-back sur la frontière est
		lat.f_[j][6] = f_star[j][8];
		//lat.f_[j][7] = f_star[j][5];
		break;
	}
}
void bounceback_N_BC( int j, int cas, Lattice lat, double** f_star)
{
	switch (cas)
	{
		case 3: //N
		lat.f_[j][4] = f_star[j][2];
		lat.f_[j][7] = f_star[j][5];
		lat.f_[j][8] = f_star[j][6];
		break;
		case 7: //NO
		lat.f_[j][4] = f_star[j][2];
		//lat.f_[j][7] = f_star[j][5];
		//lat.f_[j][8] = f_star[j][6];
		break;
		case 8: //NE
		lat.f_[j][4] = f_star[j][2];
		//lat.f_[j][7] = f_star[j][5];
		//lat.f_[j][8] = f_star[j][6];
		break;
	}
}
void bounceback_S_BC( int j, int cas, Lattice lat, double** f_star)
{
	switch (cas)
	{
		case 4: //S
		lat.f_[j][2] = f_star[j][4];
		lat.f_[j][5] = f_star[j][7];
		lat.f_[j][6] = f_star[j][8];
		break;
		case 5: //SO
		lat.f_[j][2] = f_star[j][4];
		lat.f_[j][5] = f_star[j][7];
		lat.f_[j][6] = f_star[j][8];
		break;
		case 6: //SE
		lat.f_[j][2] = f_star[j][4];
		lat.f_[j][5] = f_star[j][7];
		lat.f_[j][6] = f_star[j][8];
		break;
	}
}





