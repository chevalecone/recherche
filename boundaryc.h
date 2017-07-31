
#include "lattice.h"
void propagation( int const& j, Lattice lat, double** f_star,int nx, int ny, int* cas,bool* typeLat, int** conn);
void solidPropagation( int j, int nx, int ny , int cas, Lattice lat, double** f_star); 

void periodic_WE_BC( int j,int nx, int ny,  int cas, Lattice lat, double** f_star);
void periodic_NS_BC( int j,int nx, int ny,  int cas, Lattice lat, double** f_star);


void bounceback_N_BC( int j, int cas, Lattice lat, double** f_star);
void bounceback_S_BC( int j, int cas, Lattice lat, double** f_star);
void bounceback_E_BC( int j, int cas, Lattice lat, double** f_star);
void bounceback_W_BC( int j, int cas, Lattice lat, double** f_star);

void bounceback_N_Collision(int j, int cas, Lattice lat, double** f_star);
void bounceback_S_Collision(int j, int cas, Lattice lat, double** f_star);
void bounceback_E_Collision(int j, int cas, Lattice lat, double** f_star);
void bounceback_W_Collision(int j, int cas, Lattice lat, double** f_star);

//SOLIDES
void bounceback_solid_BC( int nx,int const& j, Lattice lat, double** f_star, int** const& conn, bool* typeLat,  int* const& bb, double& nombre, int* pos);

//Special BC for lid-driven cavity
void driven_cavity_nord( int j,  int cas, Lattice lat, double xi_r,double v_e);

//ZOU-HE
void pression_in_BC(int j, int cas, Lattice lat, double xi_r,double rho_in);
void pression_out_BC(int j, int cas, Lattice lat, double xi_r,double rho_out);
void vitesse_in_BC( int j,int nx, int cas, Lattice lat, double xi_r,double** v_in);
void vitesse_out_BC( int j,int nx, int cas, Lattice lat, double xi_r,double** v_out);

