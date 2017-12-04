	#ifndef LATTICE_H
#define LATTICE_H


class Lattice
{
public:
    Lattice(int n,  int Q,  int D);
    double **f_;
    double **f0_;
    double **u_;
    double *rho_;
    double **vorticity_;
    double **m_;
    double **m0_;
	double** q_;
};

#endif // LATTICE_H
