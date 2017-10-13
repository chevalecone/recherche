#include "geometry.h"

void D2Q9(double*omega_i,double**xi,double xi_r)
{
	//Initialisation des omega
	for (int k=0;k<9;k++)
	{
		if (k==0) 
		{omega_i[k] = (double) 4/9  ;}
		if (k<=4 && k>=1) 
		{omega_i[k] = (double) 1/9  ;}
		if (k>=5)
		{omega_i[k] = (double) 1/36 ;}
	}
	 xi[0][0]=0.0;
	 xi[0][1]=0.0;
	 xi[1][0]=xi_r;
	 xi[1][1]=0.0;
	 xi[2][0]=0.0;
	 xi[2][1]=xi_r;
	 xi[3][0]=-xi_r;
	 xi[3][1]=0.0;
	 xi[4][0]=0.0;
	 xi[4][1]=-xi_r;
	 xi[5][0]=xi_r;
	 xi[5][1]=xi_r;
	 xi[6][0]=-xi_r;
	 xi[6][1]=xi_r;
	 xi[7][0]=-xi_r;
	 xi[7][1]=-xi_r;
	 xi[8][0]=xi_r;
	 xi[8][1]=-xi_r;
}

void rang(int ny, double dx, int nx, double** rank)
{
	for (int i=0;i<ny;i++)
	{
		rank[i][0] = 0.5*dx+dx*i;
	}
}