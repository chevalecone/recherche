// C++ Standard includes
#include <iostream> 
#include <cmath>
#include <ctime>
#include <stdio.h>

// Library includes
#include "tinyxml2.h"

// Local includes
#include "parser.h"
#include "domain.h"
#include "lattice.h"
#include "solutionExporterEulerian.h"
#include "function.h"
#include "relaxation_time.h"

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

void MRT_forcing(double** F, Lattice lat, double* Fi, int N)
{
	for (int j=0;j<N;j++)
	{
	F[j][0] = 0;
	F[j][1] = 6*(lat.u_[j][0]*Fi[0]+lat.u_[j][1]*Fi[1]);
	F[j][2] = -6*(lat.u_[j][0]*Fi[0]+lat.u_[j][1]*Fi[1]);
	F[j][3] = Fi[0];
	F[j][4] = -Fi[0];
	F[j][5] = Fi[1];
	F[j][6] = -Fi[1];
	F[j][7] = 2*(lat.u_[j][0]*Fi[0]-lat.u_[j][1]*Fi[1]);
	F[j][8] = lat.u_[j][0]*Fi[0]+lat.u_[j][1]*Fi[1];
	}
}

void MRT_moment(int j,int k, Lattice lat, float** M, int Q)
{
	double sum = 0;
	for (int i=0;i<Q;i++)
		{
			sum+=M[k][i]*lat.f_[j][i];
		}
	lat.m_[j][k] = sum;
	//printf("Lattice %d, valeur du moment n° %d : %0.3f\n",j,k,lat.m_[j][k]);
}



float** MRT_matrice_passage(int Q) //Matrice de passage entre les fonctions de distribution et les moments
{
	float** matrix = new float*[Q];
	for ( int j=0;j<Q;j++)
	{
		matrix[j]=new float[Q];		
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
}



float** MRT_S(int Q, double nu, double cs, double dt) // Matrice S des temps de relaxation relatif aux différents moments
{
	float** matrix = new float*[Q];
	for ( int j=0;j<Q;j++)
	{
		matrix[j]=new float[Q];	
		for ( int k=0;k<Q;k++)
		{
			matrix[j][k]=0;	
		}	
	}
	
	float s7 = 2/(1+6*nu);
	float s8 = 2/(1+6*nu);
	float tau_s =1/(0.5+ nu/(dt*cs*cs));
	matrix[0][0] = 1;
	matrix[1][1] = 1.4;
	matrix[2][2] = 1.4;
	matrix[3][3] = 1;
	matrix[4][4] = 1.2;
	matrix[5][5] = 1;
	matrix[6][6] = 1.2;
	matrix[7][7] = tau_s;
	matrix[8][8] = tau_s;
	return matrix;
}

void MRT_collision( int j, double** f_star, float** C, Lattice lat,  int Q, double dt)
{
	double sum=0;
	for (int k=0;k<Q;k++)
	{
		for ( int i=0;i<Q;i++)
		{
		sum+=C[k][i]*(lat.m_[j][i]-lat.m0_[j][i]);
		}
		f_star[j][k] = lat.f_[j][k]-sum;
		sum=0;
	}
	
}

void MRT_forcing_collision(int j, double** f_star, float** C, Lattice lat,  int Q, double dt, float** F, float** C3, float** F_bar)
{
	double sum=0;
	double sum2=0;
	for (int k=0;k<Q;k++)
	{
		for (int i=0;i<Q;i++)
		{
		sum+=C3[k][i]*F_bar[j][i];
		}
		F[j][k] = sum;
		sum=0;
	}

	for (int k=0;k<Q;k++)
	{
		for (int i=0;i<Q;i++)
		{
		sum2+=C[k][i]*(lat.m_[j][i]-lat.m0_[j][i]);
		}
		f_star[j][k] = lat.f_[j][k]-sum+dt*F[j][k];
		sum2=0;
	}
}


void MRT_equilibre(int j, int k,Lattice lat, int Q, float** M)
{
	/*double sum = 0;
	for (int i=0;i<Q;i++)
	{
		sum+=M[k][i]*lat.f0_[j][i];
	}
	lat.m0_[j][k] = sum;*/
	lat.m0_[j][0] = lat.rho_[j] * 1; //rho
	lat.m0_[j][1] = lat.rho_[j] * (-2 + 3 * (lat.u_[j][0]*lat.u_[j][0] + lat.u_[j][1] * lat.u_[j][1])); //e = rho * (-2 + 3 * |u|^2)
	lat.m0_[j][2] = lat.rho_[j] * (1 - 3 * (lat.u_[j][0]*lat.u_[j][0] + lat.u_[j][1] * lat.u_[j][1])); //epsilon = rho * (1- 3 * |u|^2)
	lat.m0_[j][3] = lat.rho_[j] * lat.u_[j][0]; //jx = rho * u
	lat.m0_[j][4] = lat.rho_[j] * (-lat.u_[j][0]); //qx = - rho * u
	lat.m0_[j][5] = lat.rho_[j] * lat.u_[j][1] ; //jy= rho * v 
	lat.m0_[j][6] = lat.rho_[j] * (-lat.u_[j][1]); //qy = - rho * v
	lat.m0_[j][7] = lat.rho_[j] * (lat.u_[j][0]-lat.u_[j][0] - lat.u_[j][1]*lat.u_[j][1]); //pxx = rho * (u*u - v*v)
	lat.m0_[j][8] = lat.rho_[j] * (lat.u_[j][0]*lat.u_[j][1]);//pxy = rho * u * v
	/*for (int i =0;i<9;i++)
	{
		printf("Lattice %d, valeur du moment d'équilibre n° %d : %0.3f\n",j,i,lat.m0_[j][i]);
	}*/
}

// matrix inversion
// the result is put in Y
void MatrixInversion(float **A, int order, float **Y)
{
    // get the determinant of a
    double det = 1.0/CalcDeterminant(A,order);
 
    // memory allocation
    float *temp = new float[(order-1)*(order-1)];
    float **minor = new float*[order-1];
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
int GetMinor(float **src, float **dest, int row, int col, int order)
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
double CalcDeterminant( float **mat, int order)
{
    // order must be >= 0
    // stop the recursion when matrix is a single element
    if( order == 1 )
        return mat[0][0];
 
    // the determinant value
    float det = 0;
 
    // allocate the cofactor matrix
    float **minor;
    minor = new float*[order-1];
    for(int i=0;i<order-1;i++)
        minor[i] = new float[order-1];
 
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

void matrix_product(float **A, float **B,float **C,  int Q) // produit matriciel de A et B dans C, C = A * B
{
	for ( int i=0;i<Q;i++)
	{
		for ( int j=0;j<Q;j++)
		{
			for ( int k=0;k<Q;k++)
			{
				C[i][j]+= A[i][k]*B[k][j];
			}
		}
	}
}

void affichage_matrix(int Q, float** M, float** invM, float** Si, float** C, float** F, float** F2)
{
	printf("Matrice de passage : \n");
	for (int i =0;i<Q;i++)
	{
		for (int j=0;j<Q;j++)
		{
			printf("%0.0f \t",M[i][j]);
		}
		printf("\n");
	}
	
	printf("Matrice inverse de passage : \n");
	for (int i =0;i<Q;i++)
	{
		for (int j=0;j<Q;j++)
		{
			printf("%0.3f \t",invM[i][j]);
		}
		printf("\n");
	}
	
	printf("Matrice de relaxation : \n");
	for (int i =0;i<Q;i++)
	{
		for (int j=0;j<Q;j++)
		{
			printf("%0.3f \t",Si[i][j]);
		}
		printf("\n");
	}
	
	printf("Matrice du produit M-1 * S : \n");
	for (int i =0;i<Q;i++)
	{
		for (int j=0;j<Q;j++)
		{
			printf("%0.3f \t",C[i][j]);
		}
		printf("\n");
	}

	printf("Matrice de forçage F: \n");
	for (int i=0;i<Q;i++)
	{
		printf("%0.3f \n",F[i]);
	}

	printf("Matrice du produit (I-0.5*S)*F: \n");
		for (int i =0;i<Q;i++)
	{
		for (int j=0;j<Q;j++)
		{
			printf("%0.3f \t",F2[i][j]);
		}
		printf("\n");
	}
}



void MRT_forcing(int j, int k, double cs, double* omega_i, double** xi, double* Fi, float** F_bar, Lattice lat, int D)
{
	//Rappel : F_bar = omega_i * [ci*F/cs² + uF : (cici-cs²I)/cs^4]
	//Rappel : F = invM*(I- 1/2 * S)*M*F_bar

	double sum = 0; //produit scalaire ci * F
	double sum2 = 0; //somme de uF:(cici-cs²I)
	double** matrice1 = new double*[D]; //matrice uF
	double** matrice2 = new double*[D]; //matrice cici
	double** matrice3 = new double*[D]; //matrice cs²I
	double** matrice4 = new double*[D]; //matrice cici-cs²I

	
	for (int l =0;l<D;l++)
	{
		matrice1[l] = new double[D];
		matrice2[l] = new double[D];
		matrice3[l] = new double[D];
		matrice4[l] = new double[D];
		for (int m=0;m<D;m++)
		{
			matrice1[l][m] = lat.u_[j][l]*Fi[m];
			matrice2[l][m] = xi[k][l]*xi[k][m];
			if(l==m)
			{
				matrice3[l][m] = cs*cs;
			}
			else
			{
				matrice3[l][m] = 0;	
			}
		}
		for (int m=0;m<D;m++)
		{
			matrice4[l][m] = matrice2[l][m]-matrice3[l][m];
		}
	}	
	for (int l =0;l<D;l++)
	{
		sum+=xi[k][l]*Fi[l];
		for (int m=0;m<D;m++)
		{
			sum2+=matrice1[l][m]*matrice4[m][l];
		}
	}
	F_bar[j][k] = omega_i[k]*( 1/(cs*cs) * sum + 1/(cs*cs*cs*cs) * sum2);
	free(matrice1);
	free(matrice2);
	free(matrice3);
	free(matrice4);
}

