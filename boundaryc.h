
#include "lattice.h"
void propagation( int const& j, Lattice lat, double** f_star,int nx, int ny, int* cas,bool* typeLat, int** conn, int* bb);
 

void periodic_WE_BC( int j,int nx, int ny,  int cas, Lattice lat, double** f_star);
void periodic_NS_BC( int j,int nx, int ny,  int cas, Lattice lat, double** f_star);


void bounceback_N_BC( int j, int cas, Lattice lat, double** f_star);
void bounceback_S_BC( int j, int cas, Lattice lat, double** f_star);
void bounceback_E_BC( int j, int cas, Lattice lat, double** f_star);
void bounceback_W_BC( int j, int cas, Lattice lat, double** f_star);

//SOLIDES
void bounceback_solid_BC( int nx,int const& j, Lattice lat, double** f_star, int** const& conn, bool* typeLat,  int* const& bb, double& nombre, int* pos);

//Special BC for lid-driven cavity
void driven_cavity_nord( int j,  int cas, Lattice lat, double xi_r,double v_e);

//ZOU-HE
void pression_in_BC(int j, int cas, Lattice lat, double xi_r,double rho_in);
void pression_out_BC(int j, int cas, Lattice lat, double xi_r,double rho_out);
void vitesse_in_BC( int j,int nx, int cas, Lattice lat, double xi_r,double** v_in);
void vitesse_out_BC( int j,int nx, int cas, Lattice lat, double xi_r,double** v_out);

//COMBINED BOUNCE-BACK SPECULAR REFLECTION (CBBSR)
void CBBSR_N_BC(int j, int cas, Lattice lat, double r, double** f_star);
void CBBSR_S_BC(int j, int cas, Lattice lat, double r, double** f_star);

//EXTRAPOLATION METHOD
void extrapolation_inlet_BC(int j, int cas, Lattice lat);
void extrapolation_outlet_BC(int j, int cas, Lattice lat);

//EQUILIBRIUM METHOD
void equilibrium_inlet_BC(int j, int cas, Lattice lat, double rho_in, double cs, double* omega_i, double** xi, int Q, int nx);
void equilibrium_outlet_BC(int j, int cas, Lattice lat, double rho_out, double cs, double* omega_i, double** xi, int Q);


//DIFFUSE BOUNCE-BACK (DBB) non automatisé (ie. simplifié pour le D2Q9) avec une vitesse u_wall nulle
void DBB_N_BC(int j, int cas, Lattice lat, double beta, double** f_star);
void DBB_S_BC(int j, int cas, Lattice lat, double beta, double** f_star);

//COMBINED SPECULAR AND DIFFUSIVE REFLECTION (MR)
void MR_N_BC(int j, int cas, Lattice lat, double beta, double** f_star);
void MR_S_BC(int j, int cas, Lattice lat, double beta, double** f_star);