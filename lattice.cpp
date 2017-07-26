#include "lattice.h"

Lattice::Lattice(int n,  int Q,  int D)
{
    f_ = new double*[n];
    f0_ = new double*[n];
    u_ = new double*[n];
    rho_ = new double[n];
    vorticity_ = new double*[n];
    m_ = new double*[n];
    m0_ = new double*[n];

    for (int i=0 ; i < n ; i ++)
    {
        f_[i] = new double[Q];
        f0_[i] = new double[Q];
        vorticity_[i] = new double[D];
        u_[i] = new double[D];
        m_[i] = new double[Q];
        m0_[i] = new double[Q];
    }
}
