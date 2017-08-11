double* MRT_Guo_2008(double Kn, double ymax, double tau_s, double tau_q, double r, double dx);
double* MRT_Continuous(double cs, double dt, double nu);
double* MRT_Verhaeghe_2009(double Kn, double ymax, double tau_s, double tau_q, double beta, double mu, double dx, double dt, double cs);
void PHI_Guo_2008(double Kn, double ymax,int N, double dx, double* PHI, double** position, double mu, double rho, double cs);