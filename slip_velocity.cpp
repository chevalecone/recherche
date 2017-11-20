// C++ Standard includes
#include <iostream> 
#include <cmath>
#include <ctime>
#include <stdio.h>
#include <cstring>


// Local includes
#include "function.h"
#include "slip_velocity.h"
#include "lattice.h"

# define M_PI  3.14159265358979323846

double* Slip_velocity(std::string slip, double sigma, double Kn)
{
	double* A = new double[2];
	if (!slip.compare("Guo-2008"))
	{
		A[0] = (2-sigma)/sigma * (1-0.1817*sigma);
		A[1] = 1/(M_PI) + 0.5*A[0]*A[0];
	}
	else if(!slip.compare("Guo-2011"))
	{
		A[0] = 1.146;
		A[1] = 0.976;
	}
	else if(!slip.compare("Wang-2017"))
	{
		double nA1 = 2-sigma+sigma*(-1+2*sigma)/(Kn*Kn*Kn)*(-Ei_big(1,1/Kn))+sigma*exp(-1/Kn)*((-1+2*sigma)/(Kn*Kn)+(1-2*sigma)/Kn+(-2+sigma))-(4-4*sigma)/(Kn*Kn*Kn)*(-Ei_big(1,2/Kn))-exp(-2/Kn)*((2-2*sigma)/(Kn*Kn)-(1-sigma)/Kn+(1-sigma));
		double dA1 = (2-sigma)*(1+sigma/(Kn*Kn)*(-Ei_big(1,1/Kn))+exp(-1/Kn)*(sigma/Kn-sigma)+(4-4*sigma)/(Kn*Kn)*(-Ei_big(1,2/Kn))+exp(-2/Kn)*((2-2*sigma)/Kn-(1-sigma)));
		double nA2 = 6+(-8+11*sigma)/(Kn*Kn*Kn*Kn)*(-Ei_big(1,1/Kn))+exp(-1/Kn)*((-8+11*sigma)/(Kn*Kn*Kn)+(8-11*sigma)/(Kn*Kn)+(-16+10*sigma)/Kn-6*sigma)+(16-16*sigma)/(Kn*Kn*Kn*Kn)*(-Ei_big(1,2/Kn))+exp(-2/Kn)*((8-8*sigma)/(Kn*Kn*Kn)-(4-4*sigma)/(Kn*Kn)+(4-4*sigma)/Kn-(6-6*sigma));
		double dA2 = 12+12*sigma/(Kn*Kn)*(-Ei_big(1,1/Kn))+exp(-1/Kn)*(12*sigma/Kn-12*sigma)+(48-48*sigma)/(Kn*Kn)*(-Ei_big(1,2/Kn))+exp(-2/Kn)*((24-24*sigma)/Kn-(12-12*sigma));
		printf("nA1 dA1 nA2 dA2 : %f %f %f %f\n",nA1,dA1,nA2,dA2);
		double A1 = 2./3.*nA1/dA1;
		double A2 = nA2/dA2;
		A[0] = A1;
		A[1] = A2;
	}
	else if (!slip.compare("Hadjiconstantinou-2003"))	
	{
		A[0] = 1.11*(2-sigma)/sigma;
		A[1] = 0.61;
	}
	else if (!slip.compare("Li-2011"))
	{
		A[0] = (1-0.1817*sigma);
		A[1] = 0.8;
	}
	else if (!slip.compare("Wu-2008"))
	{
		double f;
		if (Kn>1)
		{
			f = 1/Kn;
		}
		else
		{
			f = 1;
		}
		A[0] = 2./3.*(3-sigma*f*f*f/(sigma)-3*(1-f*f)/(2*Kn));
		A[1] = 1./4.*(f*f*f*f+2/(Kn*Kn)*(1-f*f));
	}
	else if (!slip.compare("de Izarra-2012"))
	{
		A[0] = 1;
		A[1] = 0.13;
	}
	return A;
}
	