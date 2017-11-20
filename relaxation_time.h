#include <string>

//****************************MRT*******************************************//
void MRT_moment(int j,int k, Lattice lat, double** M, int Q); // calcul des moments (produit matriciel entre les populations et la matrice de passage)
void MRT_equilibre(int j,Lattice lat); //calcul des moments d'équilibre 
void MRT_equilibre_v2(int j, Lattice lat, double** 	M, int Q, double* temp);

void MRT_collision( int j, double** f_star, double** C, Lattice lat,  int Q, double* temp); // collision entre les populations en utilisant les moments
void MRT_collision_v2( int j, double** f_star, double** C, Lattice lat,  int Q, double* temp);
void MRT_forcing_collision(int j, double** f_star, double** C, Lattice lat, double dt, double** F, double** C3, double** F_bar, int Q, double* temp);
void MRT_forcing(int j, int Q, double cs, double* omega_i, double** xi, double* Fi, double** F_bar, Lattice lat, double* temp);

double** MRT_matrice_passage(int Q); //Matrice de passage des populations aux moments (Gram Schmidt)
double** MRT_S(int Q, double tau_s, double tau_q); //Matrice de relaxation
void MRT_S(int Q, double tau_s, double tau_q, double** Si);
void matrix_product(double **A, double **B,double **C,  int Q); //Produit matriciel


//Utiles pour calculer l'inverse d'une matrice
void MatrixInversion(double **A, int order, double **Y);
int GetMinor(double **src, double **dest, int row, int col, int order);
double CalcDeterminant( double **mat, int order);
//*****************************TRT************************************//
void TRT_creation(int j, Lattice lat, double** f_minus, double** f_plus, double** f0_plus, double** f0_minus, int* bb, int Q);

//Affichage matrice
void affichage_matrix(int Q, double** matrix);