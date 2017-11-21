// C++ Standard includes
#include <iostream> 
#include <cmath>
#include <ctime>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <time.h>


// Local includes
#include "domain.h"
#include "lattice.h"
#include "solutionExporterEulerian.h"
#include "function.h"

# define M_PI  3.14159265358979323846

#define EULER 0.5772156649 //Constante d'Euler
#define MAXIT 100 //Nombre maximum d'itérations
#define FPMIN 1.0e-30
#define EPS 1.0e-8 

//Fonction pour le produit scalaire en D dimensions
double pscal(double* a,double* b, int D, double sigma) 
{
	sigma=0;
	for ( int i=0;i<D;i++)
	{
		sigma+=a[i]*b[i];
	}
	return sigma;
}

//Fonction pour la somme des fi pour calculer la masse volumique
void density ( int j, int Q, Lattice lat, double sigma) 
{
	sigma = 0;
	for ( int k=0;k<Q;k++)
	{
		sigma+=lat.f_[j][k];
	}
	lat.rho_[j] = sigma;
}

//Fonction pour calculer la vitesse macroscopique en D dimensions
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

//Fonction donnant un numéro à la lattice, déterminant son cas (noté de 0 à 8 pour un domaine simple)
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


//Fonction donnant la localisation de tous les noeuds de lattices (pour des noeuds mid grid, le noeud est au milieu de la lattice)
void localisation(int nx, int ny, double dx, double** position)
{
	for ( int i=0;i<nx*ny;i++)
	{
		position[i][0] = 0.5*dx+dx*(i%nx);
		position[i][1] = 0.5*dx+dx*(int)(i/nx);
	}
}

//Localisation d'un cylindre carré en fonction de ses coordonnées de son milieu et de la longueur du côté
//Fonction donnant les 4 coins du carré
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

//Création d'un cylindre circulaire en fonction de ses coordonnées de son centre et de son diamètre
void typeCircular(double abscisse, double ordonnee, double diametre, int N, double** position, bool* typeLat)
{
	for (int i=0;i<N;i++)
	{
		// si (xm-xo)² + (ym-yo)² <= D²/4, alors la lattice fait partie du cylindre, il est solide
		if(((position[i][0]-abscisse)*(position[i][0]-abscisse)+(position[i][1]-ordonnee)*(position[i][1]-ordonnee))<=diametre*diametre*0.25)
		{
			typeLat[i] = true;
		}
	}
}

//Création d'une ellipse en fonction de son centre, les longueurs des demi-grands axes et de son orientation en degrés
void typeEllipse(double abscisse, double ordonnee, double a, double b, double orientation, int N, double** position, bool*typeLat)
{
	double* pos_f = new double[2]; //position du premier foyer
	double* pos_f_p = new double[2]; //position du second foyer
	double c  = sqrt(std::max(a,b)*std::max(a,b)-std::min(a,b)*std::min(a,b)); //c = sqrt(b²-a²)
	double d_f,d_fp ; //distance du point M aux foyer principal et secondaire
	pos_f[0] = abscisse + cos(orientation*M_PI/180)*c;
	pos_f[1] = ordonnee + sin(orientation*M_PI/180)*c;
	pos_f_p[0] = abscisse - cos(orientation*M_PI/180)*c;
	pos_f_p[1] = ordonnee - sin(orientation*M_PI/180)*c; 
	for (int i = 0; i<N;i++)
	{
		d_f = sqrt((position[i][0]-pos_f[0])*(position[i][0]-pos_f[0])+(position[i][1]-pos_f[1])*(position[i][1]-pos_f[1]));
		d_fp = sqrt((position[i][0]-pos_f_p[0])*(position[i][0]-pos_f_p[0])+(position[i][1]-pos_f_p[1])*(position[i][1]-pos_f_p[1]));
		//Condition pour qu'une lattice soit dans une ellipse : (d(MF) + d(MF'))<2a
		if ((d_f+d_fp)<=2*std::max(a,b)) 
		{
			typeLat[i] = true;
		}
	}
}

//Fonction générant au hasard des ellipses jusqu'à la porosité considérée, SANS encastrement
void randomEllipse(int nx, int ny, double xmin, double xmax, double ymin, double ymax, int N, double** position, bool* typeLat, double poro, double nombre)
{
	bool* typeLat_buf = new bool[N];
	double* tableau[50]; // Tableau qui enregistre les ellipses (abscisse, ordonnee, demi axes a et b, orientation), maximum de 50 cylindres ici
	for (int i=0;i<50;i++)
	{
		tableau[i] = new double[5];
		for (int j=0;j<5;j++)
		{
			tableau[i][j] = 0;
		}
	}
	int num = 0;
	double poro2;
	int buf = 0;
	srand(time(NULL));
	int rnd  = rand()%(nx-20)+10;
	int rnd2 = rand()%(ny-20)+10;
	int a = rand()%10 + 5; // demi grand axe variant entre 5 et 15 % de la largeur du domaine
	int b = rand()%8 + 2; //demi petit axe variant entre 2 et 10 % de la largeur du domaine
	int orientation = rand()%180; //Orientation entre 0 et 180°
	double abscisse =  position[rnd][0];  // abscisse entre 10 et 90 % de la longueur
	double ordonnee =  position[rnd2*100][1];  // ordonnee entre 10 et 90 % de la largeur
	srand(time(NULL));
	//On veut que la porosité générée se situe entre 99 et 101 % de la porosité idéale voulue
	while (poro2<=0.99*poro || poro2>=1.01*poro)
	{

		typeEllipse(abscisse,ordonnee,a,b,orientation, N,position,typeLat_buf);
		for(int j=0;j<N;j++)
		{
			if(typeLat_buf[j]==true && typeLat[j]==true)
			{
				buf++;
			}
			// Ici, on ne veut pas que le solide soit crée trop proche de l'entrée/sortie
			if(typeLat_buf[j]==true && (position[j][0]<0.03*nx || position[j][0]>0.97*nx)) 
			{
				buf++;
			}
		}
		//Si les conditions en haut sont toutes fausses (pas d'encastrement, et solide bien situé), alors on crée l'ellipse dans le domaine
		if(buf==0)
		{
			tableau[num][0] = abscisse;
			tableau[num][1] = ordonnee;
			tableau[num][2] = a;
			tableau[num][3] = b;
			tableau[num][4] = orientation;
			num++;
			for (int j=0;j<N;j++)
			{
				if(typeLat_buf[j]==true)
				{
					typeLat[j] = true;
				}
			}
		}
		rnd  = rand()%(nx-20)+10;
		rnd2 = rand()%(nx-20)+10;
		a = rand()%10 + 5; // demi grand axe variant entre 5 et 15 % de la largeur du domaine
	 	b = rand()%8 + 2; //demi petit axe variant entre 2 et 12 % de la largeur du domaine
	 	orientation = rand()%180; //Orientation entre 0 et 180°
		abscisse =  position[rnd][0];
		ordonnee =  position[rnd2*100][1];	
		buf = 0;
		poro2 = porosite (typeLat,nombre,N); //Calcul de la porosité
		printf("Porosité : %f\n",poro2);
		for (int j=0;j<N;j++)
		{
			typeLat_buf[j] = false;
		}
	
	}
	for (int i=0;i<50;i++)
	{
		if(tableau[i][0]!=0)
		{
			printf("Ellipse %d, abscisse %f ordonnee %f a : %f b : %f orientation : %f °\n",i,tableau[i][0],tableau[i][1],tableau[i][2],tableau[i][3],tableau[i][4]);
		}
	}
}

//Fonction générant des cylindres circulaires jusqu'à la porosité voulue, SANS encastrement
void randomCircular(int nx, int ny, double xmin,double xmax, double ymin, double ymax, int N, double** position, bool* typeLat, double poro, double nombre)
{
	bool* typeLat_buf = new bool[N];
	double* tableau[50]; // Tableau qui enregistre les cylindres circulaire (abscisse, ordonnee, diametre), maximum de 50 cylindres ici
	for (int i=0;i<50;i++)
	{
		tableau[i] = new double[3];
		for (int j=0;j<3;j++)
		{
			tableau[i][j] = 0;
		}
	}
	int num = 0;
	double poro2;
	int buf = 0;
	srand(time(NULL));
	int rnd  = rand()%(nx-20) + 10 ;
	int rnd2 = rand()%(ny-20) + 10 ;
	int ratio = rand()%30 + 10; //diametre entre 10 et 40 % de la largeur
	double abscisse =  position[rnd][0];  // abscisse entre 10 et 90 % de la longueur
	double ordonnee =  position[rnd2*100][1];  // ordonnée entre 10 et 90 % de la largeur
	srand(time(NULL));
    double diametre = (ratio)*0.01*ny;	

	while (poro2<=0.99*poro || poro2>=1.01*poro)
	{
		typeCircular(abscisse,ordonnee,diametre,N,position,typeLat_buf);
		for(int j=0;j<N;j++)
		{
			if(typeLat_buf[j]==true && typeLat[j]==true)
			{
				buf++;
			}
			if(typeLat_buf[j]==true && (position[j][0]<0.03*nx || position[j][0]>0.97*nx))
				{
					buf++;
				}
		}
		if(buf==0)
		{
			tableau[num][0] = abscisse;
			tableau[num][1]  = ordonnee;
			tableau[num][2]  = diametre;
			num++;
			for (int j=0;j<N;j++)
			{
				if(typeLat_buf[j]==true)
				{
					typeLat[j] = true;
				}
			}
		}
		rnd  = rand()%(nx-20)+10;
		rnd2 = rand()%(nx-20)+10;
		//ratio = 30;
		ratio = rand()%30 + 10;
		abscisse =  position[rnd][0];
		ordonnee =  position[rnd2*100][1];
		diametre = (ratio)*0.01*ny;	
		buf = 0;
		poro2 = porosite (typeLat,nombre,N);
		printf("Porosité : %f\n",poro2);
		for (int j=0;j<N;j++)
		{
			typeLat_buf[j] = false;
		}
	
	}
	for (int i=0;i<50;i++)
	{
		if(tableau[i][0]!=0)
		{
			printf("Cylindre %d, diametre %f\n",i,tableau[i][2]);
		}
	}
}

void randomSquare(int nx, int ny, double xmin,double xmax, double ymin, double ymax, int N, double** position, bool* typeLat, double poro, double nombre, double** cylinder)
{
	bool* typeLat_buf = new bool[N];
	double poro2 = 1;
	int buf = 0;
	srand(time(NULL));
	int rnd  = rand()%(nx-20) +10;
	int rnd2 = rand()%(ny-20) +10;
	int ratio = rand()%20 + 10; // Côté du carré entre 10 et 30 % de la largeur
	double abscisse =  position[rnd][0];  // abscisse in the range 1 to nx
	double ordonnee =  position[rnd2*100][1];  // ordonnee in the range 1 to ny
	srand(time(NULL));
    double diametre = (ratio)*0.01*ny;	
	printf("Abscisse : %f Ordonnée : %f Diametre : %f\n",abscisse, ordonnee, diametre);
	while (poro2<=0.98*poro || poro2>=1.02*poro)
	{
		SquareCylinder(abscisse,ordonnee,diametre,cylinder);
		typeSquare(N,cylinder,position,typeLat_buf);
		for(int j=0;j<N;j++)
		{
			if(typeLat_buf[j]==true && typeLat[j]==true)
			{
				buf++;
			}
			if(typeLat_buf[j]==true && (position[j][0]<0.03*nx || position[j][0]>0.97*nx || position[j][1]>0.97*ny || position[j][1]<0.03*ny))
			{
				buf++;
			}
		}
		if(buf==0)
		{
			printf("A PRENDRE\n");
			for (int j=0;j<N;j++)
			{
				if(typeLat_buf[j]==true)
				{
					typeLat[j] = true;
				}
			}
		}
		rnd  = rand()%(nx-20)+10;
		rnd2 = rand()%(nx-20)+10;
		ratio = rand()%20 + 10;
		abscisse =  position[rnd][0];
		ordonnee =  position[rnd2*100][1];
		diametre = (ratio)*0.01*ny;	
		printf("Abscisse : %f Ordonnée : %f Diametre : %f\n",abscisse, ordonnee, diametre);
		buf = 0;
		poro2 = porosite (typeLat,nombre,N);
		printf("Porosité : %f\n",poro2);
		for (int j=0;j<N;j++)
		{
			typeLat_buf[j] = false;
		}
	}
}

void nettoyage(bool* typeLat, int** conn, int N, int Q)
{
	int buf =0;
	for (int j=0;j<N;j++)
	{
		if (typeLat[j]==true)
		{
			for (int k=0;k<Q;k++)
			{
				if(conn[j][k]!=-1)
				{
					if (typeLat[conn[j][k]]==true)
					{
						buf++;
					}
				}
			}
			if(buf==0 || buf==1)
			{
				typeLat[j]=false;
			}
		}
		buf=0;
	}
}

double porosite (bool* typeLat, int nombre, int N)
{
	nombre=0;
	for (int i=0;i<N;i++)
	{
		if (typeLat[i]==true)
		{
			nombre++;
		}
	}
	return (double)(N-nombre)/N;
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


double Ei_big(int n, double x)//Donne la valeur de Ei(x) pour des grandes valeurs de x
{
	int i,ii,nm1;
	double a,b,c,d,del,fact,h,psi,ans;
	nm1=n-1;
	if (n < 0 || x < 0.0 || (x==0.0 && (n==0 || n==1)))
	printf("bad arguments in expint");
	else 
	{
		if (n == 0) ans=exp(-x)/x;
		else 
	    {
			if (x == 0.0) ans=1.0/nm1;
			else 
			{
				if (x > 1.0) 
				{
					b=x+n;
					c=1.0/FPMIN;
					d=1.0/b;
					h=d;
					for (i=1;i<=MAXIT;i++)
					{
						a = -i*(nm1+i);
						b += 2.0;
						d=1.0/(a*d+b);
						c=b+a/c;
						del=c*d;
						h *= del;
						if (fabs(del-1.0) < EPS) 
						{
							ans=h*exp(-x);
							return ans;
						}
					}
					printf("continued fraction failed in expint");
				} 
				else 
				{
					ans = (nm1!=0 ? 1.0/nm1 : -log(x)-EULER);
					fact=1.0;
					for (i=1;i<=MAXIT;i++)
					{
						fact *= -x/i;
						if (i != nm1) del = -fact/(i-nm1);
						else
						{
							psi = -EULER;
							for (ii=1;ii<=nm1;ii++) psi += 1.0/ii;
							del=fact*(-log(x)+psi);
						}
						ans += del;
						if (fabs(del) < fabs(ans)*EPS) return ans;
					}
					printf("series failed in expint");
				}
			}
		}
	}
	return ans;
}

//Rappel Qi = ci*ci - cs² * I
double*** Qi_function (int k, double cs, double** xi, int D, int Q)
{
	double*** Qi = new double**[Q];
	for (int i =0;i<Q;i++)
	{
		Qi[i] = new double*[D];
		for (int j=0;j<D;j++)
		{
			Qi[i][j] = new double[D];
			for (int k=0;k<D;k++)
			{
				if (k==j)
				{
					Qi[i][j][k] = xi[i][j]*xi[i][k]-cs*cs;
				}
				else
				{
					Qi[i][j][k] = xi[i][j]*xi[i][k];
				}
				
			}
		}
	}
	return Qi;
}

//Explication
// fi,eq = rho * omega_i * (1+1/(cs*cs)*pscal(ci,u) + 1/(2*cs^4) * Qi : uu)
// fi,eq = rho * omega_i * (1+1/(cs*cs)*pscal(ci,u) + 1/(2*cs^4) * Qeq)
// Tout ça en D2Q9
void fi_equilibre (int j, int k, double rho, double cs, Lattice lat, double* u, double** xi, int D, double*** Qi, double* buffer, double* omega_i, double sigma)
{
	for (int i =0;i<10;i++)
	{
		buffer[i] = 0;
	}
	double** sum3 = new double*[D];
	for (int i=0;i<D;i++)
	{
		sum3[i] = new double[D];
	}
	buffer[0] = pscal(xi[k],u,D,sigma);
	
	for (int i =0;i<D;i++)
	{
		for (int l =0;l<D;l++)
		{
			sum3[i][l] = u[i]*u[l];
		}
	}
	for (int i=0;i<D;i++)
	{
		for (int l=0;l<D;l++)
		{
			buffer[1]+=Qi[k][i][l]*sum3[l][i];
		}
	}
	lat.f0_[j][k] = omega_i[k]*rho*(1+1/(cs*cs)*buffer[0] + 1/(2*cs*cs*cs*cs)*buffer[1]);
	for (int i =0;i<D;i++)
	{
		delete[]sum3[i];
	}
	delete[] sum3;
}

void fi_equilibre_v2 (int j, int k, double rho, double cs, Lattice lat, double* u, double** xi, int D, double*** Qi, double* buffer, double* omega_i, double sigma, double* buffer2)
{
	for (int i =0;i<10;i++)
	{
		buffer[i] = 0;
	}
	double** sum3 = new double*[D];
	for (int i=0;i<D;i++)
	{
		sum3[i] = new double[D];
	}
	buffer[0] = pscal(xi[k],u,D,sigma);
	
	for (int i =0;i<D;i++)
	{
		for (int l =0;l<D;l++)
		{
			sum3[i][l] = u[i]*u[l];
		}
	}
	for (int i=0;i<D;i++)
	{
		for (int l=0;l<D;l++)
		{
			buffer[1]+=Qi[k][i][l]*sum3[l][i];
		}
	}
	buffer2[k] = omega_i[k]*rho*(1+1/(cs*cs)*buffer[0] + 1/(2*cs*cs*cs*cs)*buffer[1]);
	for (int i =0;i<D;i++)
	{
		delete[]sum3[i];
	}
	delete[] sum3;
	delete[] buffer;
}


	void simplified_fi_equilibre (int j, double rho, double cs, Lattice lat, double* u, double** xi, int D, double* omega_i, double* feq, int Q)
{
	for (int k=0;k<Q;k++)
	{
		feq[k] = omega_i[k]*rho*(1+1/(cs*cs)*(xi[k][0]*u[0]+xi[k][1]*u[1])+(xi[k][0]*u[0]+xi[k][1]*u[1])*(xi[k][0]*u[0]+xi[k][1]*u[1])/(2*cs*cs*cs*cs)-(u[0]*u[0]+u[1]*u[1])/(2*cs*cs));
	}
}

//0.0564 0.1128 0.1692 0.2257 0.3385 0.4514 0.6670 0.9027 1.1284 1.6926 2.2568 3.3851 4.5135 6.7703 9.0270 11.2838 16.9257
char FileName(double Kn)
{
	if (Kn==0.0564)
	{
		return 'M';
	}
	else if (Kn==0.1128)
	{
		return 'A';
	}
	else if (Kn==0.1692)
	{
		return 'N';
	}
	else if (Kn==0.2257)
	{
		return 'B';
	}
	else if (Kn==0.3385)
	{
		return 'O';
	}
	else if (Kn==0.4514)
	{
		return 'C';
	}
	else if (Kn==0.667)
	{
		return 'D';
	}
	else if (Kn==0.9027)
	{
		return 'E';
	}
	else if (Kn==1.1284)
	{
		return 'F';
	}
	else if (Kn==1.6926)
	{
		return 'P';
	}
	else if (Kn==2.2568)
	{
		return 'G';
	}
	else if (Kn==3.3851)
	{
		return 'R';
	}
	else if (Kn==4.5135)
	{
		return 'H';
	}
	else if (Kn==6.7703)
	{
		return 'I';
	}
	else if (Kn==9.0270)
	{
		return 'J';
	}
	else if (Kn==11.2838)
	{
		return 'K';
	}
	else if (Kn==16.9257)
	{
		return 'L';
	}
	else
	{
		return 'Z';
	}
}

void PI_neq_inlet(double** Pi_neq, int j, int Q, Lattice lat, double*** Qi, double* buffer, int* bb)
{
	for (int i =0;i<10;i++)
	{
		buffer[i] = 0;
	}
	for (int k=0;k<Q;k++)
	{ 
		if (k ==1 || k ==5 || k==8)
		{
			buffer[0]+=Qi[k][0][0]*(lat.f_[j][bb[k]]-lat.f0_[j][bb[k]]);
			buffer[1]+=Qi[k][0][1]*(lat.f_[j][bb[k]]-lat.f0_[j][bb[k]]);
			buffer[2]+=Qi[k][1][0]*(lat.f_[j][bb[k]]-lat.f0_[j][bb[k]]);
			buffer[3]+=Qi[k][1][1]*(lat.f_[j][bb[k]]-lat.f0_[j][bb[k]]);
		}
		else
		{
			buffer[0]+=Qi[k][0][0]*(lat.f_[j][k]-lat.f0_[j][k]);
			buffer[1]+=Qi[k][0][1]*(lat.f_[j][k]-lat.f0_[j][k]);
			buffer[2]+=Qi[k][1][0]*(lat.f_[j][k]-lat.f0_[j][k]);
			buffer[3]+=Qi[k][1][1]*(lat.f_[j][k]-lat.f0_[j][k]);
		}
	}
	Pi_neq[0][0] = buffer[0];
	Pi_neq[0][1] = buffer[1];
	Pi_neq[1][0] = buffer[2];
	Pi_neq[1][1] = buffer[3];
}

void PI_neq_outlet(double** Pi_neq, int j, int Q, Lattice lat, double*** Qi, double* buffer, int* bb)
{
	for (int i =0;i<10;i++)
	{
		buffer[i] = 0;
	}
	for (int k=0;k<Q;k++)
	{ 
		if (k ==3 || k ==6 || k==7)
		{
			buffer[0]+=Qi[k][0][0]*(lat.f_[j][bb[k]]-lat.f0_[j][bb[k]]);
			buffer[1]+=Qi[k][0][1]*(lat.f_[j][bb[k]]-lat.f0_[j][bb[k]]);
			buffer[2]+=Qi[k][1][0]*(lat.f_[j][bb[k]]-lat.f0_[j][bb[k]]);
			buffer[3]+=Qi[k][1][1]*(lat.f_[j][bb[k]]-lat.f0_[j][bb[k]]);
		}
		else
		{
			buffer[0]+=Qi[k][0][0]*(lat.f_[j][k]-lat.f0_[j][k]);
			buffer[1]+=Qi[k][0][1]*(lat.f_[j][k]-lat.f0_[j][k]);
			buffer[2]+=Qi[k][1][0]*(lat.f_[j][k]-lat.f0_[j][k]);
			buffer[3]+=Qi[k][1][1]*(lat.f_[j][k]-lat.f0_[j][k]);
		}
	}
	Pi_neq[0][0] = buffer[0];
	Pi_neq[0][1] = buffer[1];
	Pi_neq[1][0] = buffer[2];
	Pi_neq[1][1] = buffer[3];
}


void fi_bar(double* omega_i, double***Qi, double** Pi_neq, double cs, int j, int D, Lattice lat, double* buffer, int Q)
{
	for (int k=0;k<Q;k++)
	{
		for (int i =0;i<D;i++)
		{
			for (int l=0;l<D;l++)
			{
				buffer[0]+=Qi[k][i][l]*Pi_neq[l][i];
			}
		}
	lat.f_[j][k]=lat.f0_[j][k]+omega_i[k]/(2*cs*cs*cs*cs)*buffer[0];
	buffer[0]=0;
	}	
}

//Porosité : volume de solide / volume total
double porosity( bool* typeLat, int N)
{
	double poro = 0;
	int l = 0;
	for (int i=0;i<N;i++)
	{
		if(typeLat[i]==false) //Si la lattice est solide
		{
			l++;
		}
	}
	poro = l/N;
	return poro;
}

void tab_voisin(int N, int Q, bool* typeLat, int* tab_Voisin, int** conn)
{
	int buf;
	for (int j=0;j<N;j++)
	{
		buf = 0;
		for (int k=0;k<Q;k++)
		{
			if (conn[j][k]!=-1 && typeLat[conn[j][k]])
			{
				buf+=pow(2,k);
			}
		}
		tab_Voisin[j] = buf;
		if(typeLat[j])
		{
			tab_Voisin[j] = 511;
		}
	}
}



