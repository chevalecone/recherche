#include <stdlib.h>
#include <algorithm>
#include <time.h>
#include <iostream>
#include "lattice.h"


//*************************CALCULS MACROSCOPIQUES**************************//
double pscal(double *a,double *b,  int D, double sigma);
void density ( int j, int Q, Lattice lat, double sigma);
void velocity( int j, int D,  int Q, double** xi, Lattice lat, double sigma);

//***************************EQUILIBRE**************************************//
double*** Qi_function (int k, double cs, double** xi, int D, int Q);
void fi_equilibre (int j, int k, double rho, double cs, Lattice lat, double* u, double** xi, int D, double*** Qi, double* buffer, double* omega_i,double sigma);
void fi_equilibre_v2 (int j, int k, double rho, double cs, Lattice lat, double* u, double** xi, int D, double*** Qi, double* buffer, double* omega_i, double sigma, double* buffer2);
void simplified_fi_equilibre (int j, double rho, double cs, Lattice lat, double* u, double** xi, int D, double* omega_i, double* feq, int Q);
void PI_neq_inlet(double** Pi_neq, int j, int Q, Lattice lat, double*** Qi, double* buffer, int* bb);
void PI_neq_outlet(double** Pi_neq, int j, int Q, Lattice lat, double*** Qi, double* buffer, int* bb);
void fi_bar(double* omega_i, double***Qi, double** Pi_neq, double cs, int j, int D, Lattice lat, double* buffer, int Q);

//*************************GEOMETRIE DU DOMAINE**************************//

//Coordonnées des lattices
void localisation(int nx, int ny, double dx, double** position);
//Matrice de connectivité
void connectivite(int nx,int ny,  int Q, int** conn);
//Donne les directions opposées aux populations considérés (utile pour le HWBB)
void bounceback_neighbour( int* bb,  int Q);
//Cas pour les frontières du domaine et/ou solides
void domainCondition(int nx, int ny,  int* cas);
//Condition des lattices vis à vis des cylindres solides
void tab_voisin(int N, int Q, bool* typeLat, int* tab_Voisin, int** conn);

//********************************SOLIDES********************************//

//Donne la position min et max des lattices solides
void pos_solide (bool* typeLat,  int* pos, int nx, int ny);
//CYLINDRE CARRE
//Matrice des 4 coins de la section carrée du cylindre
void SquareCylinder(double abscisse, double ordonnee, double diametre, double** coin);
void typeCircular(double abscisse, double ordonnee, double diametre, int N, double** position, bool* typeLat);
void typeEllipse(double abscisse, double ordonnee, double a, double b, double orientation, int N, double** position, bool*typeLat);
//Remplit typeLat pour un cylindre carré
void typeSquare( int N, double** coin, double** position, bool* typeLat);
double porosite (bool* typeLat, int nombre, int N);
void randomCircular(int nx, int ny, double xmin,double xmax, double ymin, double ymax, int N, double** position, bool* typeLat, double poro, double nombre, int* cas);
void randomSquare(int nx, int ny, double xmin,double xmax, double ymin, double ymax, int N, double** position, bool* typeLat, double poro, double nombre, double** cylinder, int* cas);
void randomEllipse(int nx, int ny, double xmin, double xmax, double ymin, double ymax, int N, double** position, bool* typeLat, double poro, double nombre, int* cas, double a_ellipse, double b_ellipse);
void nettoyage(bool* typeLat, int** conn, int N, int Q);
//**************************WALL FUNCTION******************************//
//Donne la valeur de Ei(x)
double Ei_big(int n, double x);

//**********************UTILE POUR L'EXPORTATION****************//
char FileName(double Kn);

//**********************Caractéristiques du milieu************//
double porosity( bool* typeLat, int N);
