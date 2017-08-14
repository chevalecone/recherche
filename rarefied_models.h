double* MRT_Guo_2008(double Kn, double ymax, double tau_q, double r, double dx, int N, double** rank, double* PHI);
void PHI_Guo_2008(double Kn, double ymax, double dx, double* PHI, double** position, double mu, double rho, double cs, int N, double** rank);
double* MRT_Continuous(double cs, double dt, double nu);
double* MRT_Verhaeghe_2009(double Kn, double ymax, double tau_s, double tau_q, double beta, double mu, double dx, double dt, double cs);
