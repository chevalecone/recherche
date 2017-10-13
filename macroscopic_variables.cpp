// C++ Standard includes
#include <iostream> 
#include <cmath>
#include <ctime>
#include <stdio.h>
#include <algorithm>

// Local includes
#include "lattice.h"
#include "macroscopic_variables.h"
#include "function.h"

void density()
{
	if(cas>8 && cas<27) //Tous cas confondu, les bornes représentent les lattices solides
	lat.rho_[j]=rho0;
	else
	{
		lat.rho_[j]=sum_fi(lat.f_[j],Q);
	}
}

void velocity()
{
	double vel = new double[D];
	if(cas>8 && cas<27)//Tous cas confondu, les bornes représentent les lattices solides
	{
		for (int k =0;k<D;k++)
		{
			lat.u_[j][k] = 0;
		}
	} 
	else
	{
		for (int k =0;k<D;k++)
		{
			lat.u_[j][k] = 0;
		}
	}
}