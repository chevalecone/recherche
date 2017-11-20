// C++ Standard includes
#include <iostream> 
#include <cmath>
#include <ctime>
#include <stdio.h>

// Local includes
#include "domain.h"
#include "lattice.h"
#include "solutionExporterEulerian.h"
#include "function.h"
#include "relaxation_time.h"

# define M_PI  3.14159265358979323846

void TRT_creation(int j, Lattice lat, double** f_minus, double** f_plus, double** f0_plus, double** f0_minus, int* bb, int Q)
{
	for (int k =0;k<Q;k++)
	{
		f_minus[j][k] = 0.5*(lat.f_[j][k]-lat.f_[j][bb[k]]);
		f_plus[j][k] = 0.5*(lat.f_[j][k]+lat.f_[j][bb[k]]);
		f0_plus[j][k] = 0.5*(lat.f0_[j][k]+lat.f0_[j][bb[k]]);
		f0_minus[j][k] = 0.5*(lat.f0_[j][k]-lat.f0_[j][bb[k]]);
	}
}

void MRT_moment(int j,int k, Lattice lat, double** M, int Q)
{
	double sum = 0;
	for (int i=0;i<Q;i++)
		{
			sum+=M[k][i]*lat.f_[j][i];
		}
	lat.m_[j][k] = sum;
	//printf("Lattice %d, valeur du moment n° %d : %0.3f\n",j,k,lat.m_[j][k]);
}



double** MRT_matrice_passage(int Q) //Matrice de passage entre les fonctions de distribution et les moments
{
	double** matrix = new double*[Q];
	for ( int j=0;j<Q;j++)
	{
		matrix[j]=new double[Q];		
	}
	for ( int k =0;k<Q;k++)
	{
		matrix[0][k]=1;	
	}
	matrix[1][0] = -4;
	matrix[1][1] = -1;
	matrix[1][2] = -1;
	matrix[1][3] = -1;
	matrix[1][4] = -1;
	matrix[1][5] = 2;
	matrix[1][6] = 2;
	matrix[1][7] = 2;
	matrix[1][8] = 2;

	matrix[2][0] = 4;
	matrix[2][1] = -2;
	matrix[2][2] = -2;
	matrix[2][3] = -2;
	matrix[2][4] = -2;
	matrix[2][5] = 1;
	matrix[2][6] = 1;
	matrix[2][7] = 1;
	matrix[2][8] = 1;

	matrix[3][0] = 0;
	matrix[3][1] = 1;
	matrix[3][2] = 0;
	matrix[3][3] = -1;
	matrix[3][4] = 0;
	matrix[3][5] = 1;
	matrix[3][6] = -1;
	matrix[3][7] = -1;
	matrix[3][8] = 1;

	matrix[4][0] = 0;
	matrix[4][1] = -2;
	matrix[4][2] = 0;
	matrix[4][3] = 2;
	matrix[4][4] = 0;
	matrix[4][5] = 1;
	matrix[4][6] = -1;
	matrix[4][7] = -1;
	matrix[4][8] = 1;

	matrix[5][0] = 0;
	matrix[5][1] = 0;
	matrix[5][2] = 1;
	matrix[5][3] = 0;
	matrix[5][4] = -1;
	matrix[5][5] = 1;
	matrix[5][6] = 1;
	matrix[5][7] = -1;
	matrix[5][8] = -1;

	matrix[6][0] = 0;
	matrix[6][1] = 0;
	matrix[6][2] = -2;
	matrix[6][3] = 0;
	matrix[6][4] = 2;
	matrix[6][5] = 1;
	matrix[6][6] = 1;
	matrix[6][7] = -1;
	matrix[6][8] = -1;

	matrix[7][0] = 0;
	matrix[7][1] = 1;
	matrix[7][2] = -1;
	matrix[7][3] = 1;
	matrix[7][4] = -1;
	matrix[7][5] = 0;
	matrix[7][6] = 0;
	matrix[7][7] = 0;
	matrix[7][8] = 0;

	matrix[8][0] = 0;
	matrix[8][1] = 0;
	matrix[8][2] = 0;
	matrix[8][3] = 0;
	matrix[8][4] = 0;
	matrix[8][5] = 1;
	matrix[8][6] = -1;
	matrix[8][7] = 1;
	matrix[8][8] = -1;
	
	return matrix;

	//release memory
	for (int j=0;j<Q;j++)
	{
		delete [] matrix[j];
	}
	delete [] matrix;	
}


double** MRT_S(int Q, double tau_s, double tau_q) // Matrice S des temps de relaxation relatif aux différents moments
{
	double** matrix = new double*[Q];
	for ( int j=0;j<Q;j++)
	{
		matrix[j]=new double[Q];	
		for ( int k=0;k<Q;k++)
		{
			matrix[j][k]=0;	
		}	
	}
	
	matrix[0][0] = 1;
	matrix[1][1] = 1.4;
	matrix[2][2] = 1.4;
	matrix[3][3] = 1;
	matrix[5][5] = 1;
	

	matrix[4][4] = 1/tau_q;
	matrix[6][6] = 1/tau_q;
	matrix[7][7] = 1/tau_s;
	matrix[8][8] = 1/tau_s;

	
	return matrix;

	//release memory
	for(int j=0;j<Q;j++)
	{
		delete [] matrix[j];
	}
	delete [] matrix;
}


void MRT_S(int Q, double tau_s, double tau_q, double** Si) // Matrice S des temps de relaxation relatif aux différents moments
{
	
	Si[0][0] = 1;
	Si[1][1] = 1.4;
	Si[2][2] = 1.4;
	Si[3][3] = 1;
	Si[5][5] = 1;
	
	Si[4][4] = 1/tau_q;
	Si[6][6] = 1/tau_q;
	Si[7][7] = 1/tau_s;
	Si[8][8] = 1/tau_s;
}


//La collision avec m et m0
void MRT_collision( int j, double** f_star, double** C, Lattice lat,  int Q, double* temp)
{
	for (int k=0;k<Q;k++)
	{
		temp[0]=0;
		for ( int i=0;i<Q;i++)
		{
		temp[0]+=C[k][i]*(lat.m_[j][i]-lat.m0_[j][i]);
		}
		f_star[j][k] = lat.f_[j][k]-temp[0];
	}
}

//La collision avec f et f0 (on a donc une différence dans la matrice de passage)
void MRT_collision_v2( int j, double** f_star, double** C, Lattice lat,  int Q, double* temp)
{
	for (int k=0;k<Q;k++)
	{
		temp[0]=0;
		for ( int i=0;i<Q;i++)
		{
		temp[0]+=C[k][i]*(lat.f_[j][i]-lat.f0_[j][i]);
		}
		f_star[j][k] = lat.f_[j][k]-temp[0];
	}
}

void MRT_forcing_collision(int j, double** f_star, double** C, Lattice lat, double dt, double** F, double** C3, double** F_bar, int Q, double* temp)
{
	for (int k=0;k<Q;k++)
	{
		temp[0]=0;
		temp[1]=0;
		for (int i=0;i<Q;i++)
		{
			temp[0] += C3[k][i]*F_bar[j][i];
			temp[1] += C[k][i]*(lat.m_[j][i]-lat.m0_[j][i]);
		}
		F[j][k] = temp[0];
		f_star[j][k] = lat.f_[j][k] -temp[1] + dt*F[j][k];
	}
}


void MRT_equilibre(int j,Lattice lat)
{
	lat.m0_[j][0] = lat.rho_[j] * 1; //rho
	lat.m0_[j][1] = lat.rho_[j] * (-2 + 3 * (lat.u_[j][0]*lat.u_[j][0] + lat.u_[j][1] * lat.u_[j][1])); //e = rho * (-2 + 3 * |u|^2)
	lat.m0_[j][2] = lat.rho_[j] * (1 - 3 * (lat.u_[j][0]*lat.u_[j][0] + lat.u_[j][1] * lat.u_[j][1])); //epsilon = rho * (1- 3 * |u|^2)
	lat.m0_[j][3] = lat.rho_[j] * lat.u_[j][0]; //jx = rho * u
	lat.m0_[j][4] = lat.rho_[j] * (-lat.u_[j][0]); //qx = - rho * u
	lat.m0_[j][5] = lat.rho_[j] * lat.u_[j][1] ; //jy= rho * v 
	lat.m0_[j][6] = lat.rho_[j] * (-lat.u_[j][1]); //qy = - rho * v
	lat.m0_[j][7] = lat.rho_[j] * (lat.u_[j][0]-lat.u_[j][0] - lat.u_[j][1]*lat.u_[j][1]); //pxx = rho * (u*u - v*v)
	lat.m0_[j][8] = lat.rho_[j] * (lat.u_[j][0]*lat.u_[j][1]);//pxy = rho * u * v
}

void MRT_equilibre_v2(int j, Lattice lat, double** M, int Q, double* temp)
{
	for (int k=0;k<Q;k++)
	{
		temp[0]=0;
		for ( int i=0;i<Q;i++)
		{
			temp[0]+=M[k][i]*lat.f0_[j][i];
		}
		lat.m0_[j][k]= temp[0];
	}
}

// matrix inversion
// the result is put in Y
void MatrixInversion(double **A, int order, double **Y)
{
    // get the determinant of a
    double det = 1.0/CalcDeterminant(A,order);
 
    // memory allocation
    double *temp = new double[(order-1)*(order-1)];
    double **minor = new double*[order-1];
    for(int i=0;i<order-1;i++)
        minor[i] = temp+(i*(order-1));
 
    for(int j=0;j<order;j++)
    {
        for(int i=0;i<order;i++)
        {
            // get the co-factor (matrix) of A(j,i)
            GetMinor(A,minor,j,i,order);
            Y[i][j] = det*CalcDeterminant(minor,order-1);
            if( (i+j)%2 == 1)
                Y[i][j] = -Y[i][j];
        }
    }
 
    // release memory
    //delete [] minor[0];
    delete [] temp;
    delete [] minor;
}
 
// calculate the cofactor of element (row,col)
int GetMinor(double **src, double **dest, int row, int col, int order)
{
    // indicate which col and row is being copied to dest
    int colCount=0,rowCount=0;
    for(int i = 0; i < order; i++ )
    {
        if( i != row )
        {
            colCount = 0;
            for(int j = 0; j < order; j++ )
            {
                // when j is not the element
                if( j != col )
                {
                    dest[rowCount][colCount] = src[i][j];
                    colCount++;
                }
            }
            rowCount++;
        }
    }
 
    return 1;
}
// Calculate the determinant recursively.
double CalcDeterminant( double **mat, int order)
{
    // order must be >= 0
    // stop the recursion when matrix is a single element
    if( order == 1 )
        return mat[0][0];
 
    // the determinant value
    double det = 0;
 
    // allocate the cofactor matrix
    double **minor;
    minor = new double*[order-1];
    for(int i=0;i<order-1;i++)
        minor[i] = new double[order-1];
 
    for(int i = 0; i < order; i++ )
    {
        // get minor of element (0,i)
        GetMinor( mat, minor, 0, i , order);
        // the recusion is here!
 
        det += (i%2==1?-1.0:1.0) * mat[0][i] * CalcDeterminant(minor,order-1);
        //det += pow( -1.0, i ) * mat[0][i] * CalcDeterminant( minor,order-1 );
    }
    // release memory
    for(int i=0;i<order-1;i++)
    delete [] minor[i];

    delete [] minor;
    return det;
}

void matrix_product(double **A, double **B,double **C,  int Q) // produit matriciel de A et B dans C, C = A * B
{
	for ( int i=0;i<Q;i++)
	{
		for ( int j=0;j<Q;j++)
		{
			C[i][j] = 0;
		}
	}
	for ( int i=0;i<Q;i++)
	{
		for ( int j=0;j<Q;j++)
		{
			C[i][j] = 0;
			for ( int k=0;k<Q;k++)
			{
				C[i][j]+= A[i][k]*B[k][j];
			}
		}
	}
}

void MRT_forcing(int j, int Q, double cs, double* omega_i, double** xi, double* Fi, double** F_bar, Lattice lat, double* temp)
{
	//Rappel : F_bar = omega_i * [ci*F/cs² + uF : (cici-cs²I)/cs^4]
	//Rappel : F = invM*(I- 1/2 * S)*M*F_bar

	temp[0]=0; //produit scalaire ci * F
	temp[1]=0; //somme de uF:(cici-cs²I)
	for (int k=0;k<Q;k++)
	{
		temp[0] = (xi[k][0]*Fi[0] + xi[k][1]*Fi[1]);
		temp[1] = (lat.u_[j][0]*Fi[0]*(xi[k][0]*xi[k][0]-cs*cs) + lat.u_[j][0]*Fi[1]*xi[k][1]*xi[k][0] + lat.u_[j][1]*Fi[0]*xi[k][0]*xi[k][1] + lat.u_[j][1]*Fi[1]*(xi[k][1]*xi[k][1]-cs*cs));
		F_bar[j][k] = omega_i[k]*( 1/(cs*cs) * temp[0] + 1/(cs*cs*cs*cs) * temp[1]);
	}

	/*F_bar[j][0] = omega_i[0]*(1/(cs*cs)*(xi[0][0]*Fi[0]+xi[0][1]*Fi[1])+1/(cs*cs*cs*cs)*(lat.u_[j][0]*Fi[0]*(xi[0][0]*xi[0][0]-cs*cs) + lat.u_[j][0]*Fi[1]*xi[0][1]*xi[0][0] + lat.u_[j][1]*Fi[0]*xi[0][0]*xi[0][1] + lat.u_[j][1]*Fi[1]*(xi[0][1]*xi[0][1]-cs*cs)));
	F_bar[j][1] = omega_i[1]*(1/(cs*cs)*(xi[1][0]*Fi[0]+xi[1][1]*Fi[1])+1/(cs*cs*cs*cs)*(lat.u_[j][0]*Fi[0]*(xi[1][0]*xi[1][0]-cs*cs) + lat.u_[j][0]*Fi[1]*xi[1][1]*xi[1][0] + lat.u_[j][1]*Fi[0]*xi[1][0]*xi[1][1] + lat.u_[j][1]*Fi[1]*(xi[1][1]*xi[1][1]-cs*cs)));
	F_bar[j][2] = omega_i[2]*(1/(cs*cs)*(xi[2][0]*Fi[0]+xi[2][1]*Fi[1])+1/(cs*cs*cs*cs)*(lat.u_[j][0]*Fi[0]*(xi[2][0]*xi[2][0]-cs*cs) + lat.u_[j][0]*Fi[1]*xi[2][1]*xi[2][0] + lat.u_[j][1]*Fi[0]*xi[2][0]*xi[2][1] + lat.u_[j][1]*Fi[1]*(xi[2][1]*xi[2][1]-cs*cs)));
	F_bar[j][3] = omega_i[3]*(1/(cs*cs)*(xi[3][0]*Fi[0]+xi[3][1]*Fi[1])+1/(cs*cs*cs*cs)*(lat.u_[j][0]*Fi[0]*(xi[3][0]*xi[3][0]-cs*cs) + lat.u_[j][0]*Fi[1]*xi[3][1]*xi[3][0] + lat.u_[j][1]*Fi[0]*xi[3][0]*xi[3][1] + lat.u_[j][1]*Fi[1]*(xi[3][1]*xi[3][1]-cs*cs)));
	F_bar[j][4] = omega_i[4]*(1/(cs*cs)*(xi[4][0]*Fi[0]+xi[4][1]*Fi[1])+1/(cs*cs*cs*cs)*(lat.u_[j][0]*Fi[0]*(xi[4][0]*xi[4][0]-cs*cs) + lat.u_[j][0]*Fi[1]*xi[4][1]*xi[4][0] + lat.u_[j][1]*Fi[0]*xi[4][0]*xi[4][1] + lat.u_[j][1]*Fi[1]*(xi[4][1]*xi[4][1]-cs*cs)));
	F_bar[j][5] = omega_i[5]*(1/(cs*cs)*(xi[5][0]*Fi[0]+xi[5][1]*Fi[1])+1/(cs*cs*cs*cs)*(lat.u_[j][0]*Fi[0]*(xi[5][0]*xi[5][0]-cs*cs) + lat.u_[j][0]*Fi[1]*xi[5][1]*xi[5][0] + lat.u_[j][1]*Fi[0]*xi[5][0]*xi[5][1] + lat.u_[j][1]*Fi[1]*(xi[5][1]*xi[5][1]-cs*cs)));
	F_bar[j][6] = omega_i[6]*(1/(cs*cs)*(xi[6][0]*Fi[0]+xi[6][1]*Fi[1])+1/(cs*cs*cs*cs)*(lat.u_[j][0]*Fi[0]*(xi[6][0]*xi[6][0]-cs*cs) + lat.u_[j][0]*Fi[1]*xi[6][1]*xi[6][0] + lat.u_[j][1]*Fi[0]*xi[6][0]*xi[6][1] + lat.u_[j][1]*Fi[1]*(xi[6][1]*xi[6][1]-cs*cs)));
	F_bar[j][7] = omega_i[7]*(1/(cs*cs)*(xi[7][0]*Fi[0]+xi[7][1]*Fi[1])+1/(cs*cs*cs*cs)*(lat.u_[j][0]*Fi[0]*(xi[7][0]*xi[7][0]-cs*cs) + lat.u_[j][0]*Fi[1]*xi[7][1]*xi[7][0] + lat.u_[j][1]*Fi[0]*xi[7][0]*xi[7][1] + lat.u_[j][1]*Fi[1]*(xi[7][1]*xi[7][1]-cs*cs)));
	F_bar[j][8] = omega_i[8]*(1/(cs*cs)*(xi[8][0]*Fi[0]+xi[8][1]*Fi[1])+1/(cs*cs*cs*cs)*(lat.u_[j][0]*Fi[0]*(xi[8][0]*xi[8][0]-cs*cs) + lat.u_[j][0]*Fi[1]*xi[8][1]*xi[8][0] + lat.u_[j][1]*Fi[0]*xi[8][0]*xi[8][1] + lat.u_[j][1]*Fi[1]*(xi[8][1]*xi[8][1]-cs*cs)));
*/
}

void affichage_matrix(int Q, double** matrix)
{
	for (int i =0;i<Q;i++)
	{
		for (int j=0;j<Q;j++)
		{
			printf("%.3f \t",matrix[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}
