
#include "function.h"

double lift_force( int Q, int nx,int ny, double** f_star, double xi_r, bool* typeLat, int** conn, double** xi,  int* bb,  int* pos);
double drag_force( int Q, int nx,int ny, double** f_star, double xi_r, bool* typeLat, int** conn, double** xi,  int* bb,  int* pos, Lattice lat);
void vorticite(int Q, int nx,int ny, Lattice lat, int* cas,double dx, bool* typeLat, int** conn);

