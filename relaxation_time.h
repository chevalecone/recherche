

//******************************MRT*************************************//

void MRT_moment(int j,int k, Lattice lat, float** M, int Q); // calcul des moments (produit matriciel entre les populations et la matrice de passage)
void MRT_equilibre(int j, int k,Lattice lat, int Q, float** M); //calcul des moments d'équilibre 
void MRT_collision( int j, double** f_star, float** C, Lattice  lat,  int Q, double dt); // collision entre les populations en utilisant les moments
void MRT_forcing_collision(int j, double** f_star, float** C, Lattice lat,  int Q, double dt, float** F, float** C3, float** F_bar);
void MRT_forcing(int j, int k, double cs, double* omega_i, double** xi, double* Fi, float** F_bar, Lattice lat, int D);

float** MRT_matrice_passage(int Q); //Matrice de passage des populations aux moments (Gram Schmidt)
float** MRT_S(int Q, double nu, double cs, double dt); //Matrice de relaxation
void MRT_forcing(double** F, Lattice lat, double* Fi, int N); //Matrice de forçage
void matrix_product(float **A, float **B,float **C,  int Q); //Produit matriciel

//Affichage matrice
void affichage_matrix(int Q, float** M, float** invM, float** Si, float** C, float** F, float** F2);

//Utiles pour calculer l'inverse d'une matrice
void MatrixInversion(float **A, int order, float **Y);
int GetMinor(float **src, float **dest, int row, int col, int order);
double CalcDeterminant( float **mat, int order);

//*****************************TRT************************************//
void TRT_creation(int j, Lattice lat, double** f_minus, double** f_plus, double** f0_plus, double** f0_minus, int* bb, int Q);