
#include "function.h"

double lift_force( int Q, bool* typeLat, int** conn, double** xi,  int* bb,  int* pos, Lattice lat);
double drag_force( int Q, bool* typeLat, int** conn, double** xi,  int* bb,  int* pos, Lattice lat);
void vorticite(int nx,int ny, Lattice lat,double dx, bool* typeLat, int** conn);
double tortuosite(int nx,int ny, Lattice lat,double dx, bool* typeLat, int** conn);

