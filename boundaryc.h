
#include "lattice.h"
#include "function.h"

void propagation( int const& j, Lattice lat, double** f_star,bool* typeLat, int** conn, int* bb, int Q);
 

void periodic_WE_BC( int j,int nx, int ny,  int cas, Lattice lat, double** f_star);
void periodic_NS_BC( int j,int nx, int ny,  int cas, Lattice lat, double** f_star);
void periodic_pressure_WE_BC (int j, int nx, int ny, int cas, Lattice lat, double** f_star, double dRHO, double sigma, double cs);


void bounceback_N_BC( int j, int cas, Lattice lat, double** f_star);
void bounceback_S_BC( int j, int cas, Lattice lat, double** f_star);
void bounceback_E_BC( int j, int cas, Lattice lat, double** f_star);
void bounceback_W_BC( int j, int cas, Lattice lat, double** f_star);

//SOLIDES
void bounceback_solid_BC( int nx,int const& j, Lattice lat, double** f_star, int** const& conn, bool* typeLat,  int* const& bb, double& nombre, int* pos, int cas);
void CBBSR_solid_square_BC(int nx, int const& j, Lattice lat, double** f_star, int** const& conn, bool*  typeLat, double r, int* pos, int tabVoisin);

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
void CBBSR_N_BC_Couette(int j, int cas, Lattice lat, double r, double** f_star, double uw, double** xi, double cs, double* omega_i);


//DIFFUSE BOUNCE-BACK (DBB) non automatisé (ie. simplifié pour le D2Q9) avec une vitesse u_wall nulle
void DBB_N_BC(int j, int cas, Lattice lat, double beta, double** f_star);
void DBB_S_BC(int j, int cas, Lattice lat, double beta, double** f_star);
void DBB_N_BC_Couette(int j, int cas, Lattice lat, double beta, double** f_star,double cs, double* Uw, double* buffer, double* omega_i,double** xi, int D, int Q, double*** Qi, double sigma);

//COMBINED SPECULAR AND DIFFUSIVE REFLECTION (MR)
void MR_N_BC(int j, int cas, Lattice lat, double beta, double** f_star);
void MR_S_BC(int j, int cas, Lattice lat, double beta, double** f_star);

//REGULARIZED BOUNDARY CONDITION
void regularized_BC_v_inlet(int j,int k,Lattice lat,double cs,double** v_in, int nx,double** xi,int D,double*** Qi, double* buffer, double* omega_i,int cas,double** Pi_neq,int Q, double** f_neq, int* bb, double sigma);
void regularized_BC_p_inlet(int j,int k,Lattice lat,double cs,double rho_in,double** xi,int D,double*** Qi, double* buffer, double* omega_i,int cas,double** Pi_neq,int Q, double** f_neq, int* bb, double sigma);
void regularized_BC_v_outlet(int j,int k,Lattice lat,double cs,double** v_out, int nx,double** xi,int D,double*** Qi, double* buffer, double* omega_i,int cas,double** Pi_neq,int Q, double** f_neq, int* bb, double sigma);
void regularized_BC_p_outlet(int j,int k,Lattice lat,double cs,double rho_out,double** xi,int D,double*** Qi, double* buffer, double* omega_i,int cas,double** Pi_neq,int Q, double** f_neq, int* bb, double sigma);
