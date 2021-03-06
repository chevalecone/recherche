/*******************************************************************
*
*
********************************************************************/
// C++ Standard includes
#include <iostream> 
#include <cmath>
#include <ctime>
#include <stdio.h>
#include <omp.h>
#include <string>

// Local includes
#include "domain.h"
#include "lattice.h"
#include "solutionExporterEulerian.h"
#include "function.h"
#include "geometry.h"
#include "boundaryc.h"
#include "relaxation_time.h"
#include "rarefied_models.h"

# define M_PI  3.14159265358979323846

//Variables de la simulation
# define RHO 1
# define MU 0.14
# define DX 1
# define dRHO 1.02
# define dFi 0.000001
# define U_IN_x 0.0
# define U_IN_y 0.
# define U_OUT_x 0.
# define U_OUT_y 0.
# define R 1
# define BETA 1
//0.0564 0.1128 0.1692 0.2257 0.3385 0.4514 0.6670 0.9027 1.1284 1.6926 2.2568 3.3851 4.5135 6.7703 9.0270 11.2838 16.9257 
//0.0564 0.1692 0.3385 1.6926 3.3851
# define KNU  0.4514
# define XMIN 0
# define XMAX 100		
# define YMIN 0
# define YMAX 100
# define OUTPUT 200
# define PRECISION 0.000000001
int main()
{
	Domain domain;

    // Maillage
    double dx= DX;
    domain.setDx(dx);

    // Set the size of the physical domain
    double xmin = XMIN;
    double xmax = XMAX;
    double ymin = YMIN;
    double ymax = YMAX;
    domain.setXMin(xmin);
    domain.setXMax(xmax);
    domain.setYMin(ymin);
    domain.setYMax(ymax);

    // Propriétés physiques
    double tau = 1; //Temps de relaxation
    double mu = MU; //Viscosité dynamique du fluide
    double rho0 = RHO; //Masse volumique du fluide	

    //Temps
	double error = PRECISION;
    int outputFrequency = OUTPUT; //Fréquence d'enregistrement des valeurs

    // Finalize creation of the domain
    domain.setFinalize();

    //Nombre de points de discrétisation en espace
    int nx = domain.getNx();
    int ny = domain.getNy();
    int N = domain.getNTot(); //Nombre de lattices dans le milieu

    //mu = 5/16*sqrt(2*M_PI/3)*rho0*KNU*YMAX;
    double nu = mu/rho0; //Viscosité cinématique du fluide
    
     tau = 0.5+3*nu/dx;
  

	// Information sur la boucle temporelle
	int it = 0;
    double dt = dx;

   
    int NbSave = outputFrequency+1; //Nombre de fichiers sauvés pour l'évolution temporelle des variables


    // Données générales sur le modèle géométrique utilisé
    int Q = 9; //Nombre de discrétisation de vitesse dans une lattice (ex : Q=9 pour D2Q9)
	int D = 2; //Nombre de dimensions (ex : D =2 pour 2D)
    Lattice lat(domain.getNTot(),Q,D);
	double xi_r = dx/dt; 

	time_t timer1 = time(NULL); //Temps départ
	//********************************ALLOCATION DE LA MEMOIRE*************************************//

	//Célérité de la lattice, Re, Vitesse max, variables pour la convergence en temps, force de drag et lift	
	double cs = 1/sqrt(3)*xi_r,Re,Umax,valeur1=0,valeur2=0,erreur=1,df,lf,sigma,w_time=0,Mach, poro;
	int i,j,k,l,m;
	double* buffer = new double[10];
	double* PHI = new double[N+1]; //Les wall function pour chaque lattice
	double** rank = new double*[ny];
	double** rank2 = new double*[ny];
	double* rho_bar = new double[N];
	double* coeff_r = new double[N];
	double* coeff_beta = new double[N];
	double* coeff_sigma = new double[N];
	double* temp = new double[2];
	double* temp2 = new double[3];
	double** f_star = new double*[N]; //Matrice temporaire f* entre t et t+dt (collisions)
	double** f_neq = new double*[N];
	int *cas = new int[N]; //Matrice des cas associées aux lattices (pour les BC du domaine)
	double **position = new double*[N]; //Matrice des coordonnées de chaque milieu de lattice
	int **conn = new int*[N]; //Matrice de connectivité des lattices
	bool* typeLat = new bool[N]; //booléen pour déterminer si le noeud est solide (true) ou fluide (false)
	double** S = new double*[N];
	double* A = new double[2]; //Coefficient du slip velocity
	
	double **xi = new double*[Q];//Tableau des vecteurs des vitesses du modèle
	double** Pi_neq = new double*[D];
	double omega_i[Q];//Vecteur des poids en fonction de la direction de la vitesse dans la lattice
	
	int* bb = new  int[Q]; //array des populations de bounceback
	bounceback_neighbour(bb,Q); //remplissage de la matrice bb
	
	//Cas avec terme de forçage
	double *Fi = new double[D]; //Terme de forçage
	 for (i=0;i<D;i++)
	{
		Pi_neq[i] = new double[D];
	}
	 for (i=0;i<Q;i++)
	{
		xi[i]=new double[D];
	}
	for (j=0;j<ny;j++)
	{
		rank[j] = new double[2];
		rank2[j] = new double[2];
	}
	//Matrice des f_star,S
	for (j=0;j<N;j++)
	{
		f_star[j] = new double[Q];
		f_neq[j] = new double[Q];
		typeLat[j] = false; //Initialement, tous les noeuds sont fluides
		position[j] = new double[D];
		conn[j] = new int[Q];
		cas[j] =0; // Initialisation des cas : tous des points intérieurs
		S[j] = new double[Q];
	}

//***********************INITIALISATION DE LA GEOMETRIE DE LATTICE***************************//
	D2Q9(omega_i,xi,xi_r); //Remplissage des poids omega_i et des vecteurs vitesses en fonction de la géométrie considérée
	connectivite(nx,ny,Q,conn);// Matrice de connectivité
//***********************************POSITION DES LATTICES***********************************//
	localisation(nx,ny,dx,position); //Coordonnées de chaque noeud
	domainCondition(nx,ny,cas);	//Cas pour les BC aux frontières du domaine
	rang(ny,dx,nx,rank);
	rang(ny,dx,nx,rank2);
	double*** Qi = Qi_function (k,cs,xi,D,Q);
//******************* PHENOMENES PHYSIQUES MOTEURS DE L'ECOULEMENT **************************//
	//Initialisation du terme de forçage (dans le cas de conditions périodiques)
	Fi[0]=dFi;
	Fi[1]=0; //Gravité uniquement verticale selon -y

	//Gradient de pression (tout en gardant l'hypothèse incompressible)
	double rho_in = dRHO*rho0;
	double rho_out = rho0;

	//Vitesse imposée (tout en gardant l'hypothèse incompressible)
	double* u_in = new double[D];
	double* u_out = new double[D];
	double** v_in = new double*[ny];
	double** v_out = new double*[ny];
	double Umaxx = U_IN_x, Umaxy = U_IN_y, Umax2x = U_OUT_x, Umax2y = U_OUT_y;
	double ratio2 = 0.01;
	for (j=0;j<ny;j++)
	{
		v_in[j] = new double[2];
		v_out[j] = new double[2];
		
		Umaxx= ratio2*cs;
		v_in[j][0] = Umaxx*(4*position[j*nx][1]/ymax-4*position[j*nx][1]*position[j*nx][1]/(ymax*ymax)); 
		v_in[j][1] = 0;
		v_out[j][0] = Umax2x*(4*position[j*nx][1]/ymax-4*position[j*nx][1]*position[j*nx][1]/(ymax*ymax)); 
		v_out[j][1] = 0;	
	}

//******************************CREATION DES SOLIDES******************************************//
	//CYLINDRE A SECTION CARREE
	double ratio = 0.45, nombre =0;
	double** cylinder1  = new double*[4];
	double** cylinder2  = new double*[4];
	double** cylinder3  = new double*[4];
	double** cylinder4  = new double*[4];
	double** cylinder5  = new double*[4];	

	for (j=0;j<4;j++)
	{
		cylinder1[j] = new double[D];
		cylinder2[j] = new double[D];
		cylinder3[j] = new double[D];
		cylinder4[j] = new double[D];
		cylinder5[j] = new double[D];
	}
	SquareCylinder(0,0,ratio*ymax,cylinder1);
	SquareCylinder(xmax,0,ratio*ymax,cylinder2);
	SquareCylinder(0,ymax,ratio*ymax,cylinder3);
	SquareCylinder(xmax,ymax,ratio*ymax,cylinder4);
	//SquareCylinder(xmax/2,ymax/2,ratio*ymax,cylinder5);
	typeSquare(N,cylinder1,position,typeLat); //Tableau de booléens pour les noeuds fluide/solide
	typeSquare(N,cylinder2,position,typeLat);
	typeSquare(N,cylinder3,position,typeLat);
	typeSquare(N,cylinder4,position,typeLat);
	//typeSquare(N,cylinder5,position,typeLat);
	//typeCircular(0,0,ratio*ymax,cylinder1,N,position,typeLat);
	//typeCircular(xmax,0,ratio*ymax,cylinder1,N,position,typeLat);
	//typeCircular(0,ymax,ratio*ymax,cylinder1,N,position,typeLat);
	//typeCircular(xmax,ymax,ratio*ymax,cylinder1,N,position,typeLat);
	//typeCircular(0.5*xmax,ymax/2,ratio*ymax,cylinder1,N,position,typeLat);
	nombre=0;
	for (int i=0;i<N;i++)
	{
		if (typeLat[i]==true)
		{
			nombre++;
		}
	}
	printf("poro : %f\n",(N-nombre)/N);
	printf("Diametre : %f\n",100*ratio);
	int *pos = new int[2];
	pos_solide(typeLat, pos,nx,ny);


//****************************ECOULEMENTS RAREFIES**********************************//
	double Kn = KNU;
	double tau_s = temp2[1];
	double tau_q = temp2[2];
	double r = temp2[0];
	double beta = temp2[0];
	double sigma1 = temp2[0];
	
		//*************************slip velocity***************************//
		sigma = 1;
		//slip_velocity_Wang_2017(Kn,A,sigma);
		//slip_velocity_Guo_2008(A,sigma);
		//slip_velocity_Guo_2011(Kn,A,sigma);
		//slip_velocity_Hadjiconstantinou_2003(Kn,A,sigma);
		//slip_velocity_Wu_2008(Kn,A,sigma);
		//printf("Coefficients de glissement : A1 = %f , A2 = %f\n",A[0],A[1]);
		//************************wall function***************************//
		//PHI_Guo_Shu_2013(Kn,PHI,N);
		//PHI_Guo_2008(Kn,ymax,dx,PHI,ny,rank); //Equivalent à Stops
		//PHI_Zhang_2006(Kn,ymax,dx,PHI,ny,rank);
		//PHI_Dongari_2011(Kn,ymax,dx,PHI,ny,rank);
		//PHI_Arlemark_2010(Kn,ymax,dx,PHI,ny,rank);
		
		//*************************method*********************************//
		temp2 = MRT_Continuous(cs,dt,nu);
	
		/*temp2 = MRT_Guo_2008_CBBSR(Kn,ymax,tau_s,tau_q,r,dx,A);
		r = temp2[0];
		printf("r : %f\n",r);*/
		
		//temp2 = MRT_Verhaeghe_2009_DBB_function(Kn,ymax,tau_s,beta,dx,ny,rank,rank2,PHI);
		//temp2 =  MRT_Verhaeghe_2009_DBB(Kn,ymax,tau_s,tau_q,beta,dx);
		
		//temp2 = MRT_Guo_2008_CBBSR_function(Kn,ymax,tau_q,r,dx,ny,rank, PHI,A);
		//r = temp2[0];
	
		/*temp2 =MRT_Guo_2011_DBB(Kn,ymax,tau_s,tau_q,beta,dx,A);
		beta = temp2[0];
		printf("beta : %f\n",beta);*/
		
		/*MRT_Yudong_2016_DBB(Kn,ymax,tau_s,tau_q,beta,dx,nu, A);
		beta = temp2[0];
		printf("beta : %f\n",beta);*/
		/*temp2 = MRT_Guo_2011_DBB_function(Kn,tau_q,tau_s,ymax,beta,dx,ny,rank,rank2,PHI,A);
		beta = temp2[0];
		printf("beta : %f\n",beta);*/
	
		/*temp2 = MRT_Guo_2008_MR_function(Kn,ymax,tau_q,sigma1,dx,ny,rank,PHI,A);
		sigma1 = temp2[0];
		printf("sigma : %f\n",sigma1);*/	
		
		/*temp2 =MRT_Yudong_2016_DBB_function(Kn,ymax,tau_s,beta,dx,ny,rank,rank2,PHI);
		beta = temp2[0];
		printf("beta : %f\n",beta);*/
		
		tau_s = temp2[1];
		tau_q = temp2[2];
		printf("tau_s : %f\n",tau_s);
		printf("tau_q : %f\n",tau_q);
		printf("nu = %f\n",nu);
		printf("Kn = %f\n",Kn);
	//Fi[0] = (rho_in-rho_out)*cs*cs/(xmax);
	//Umax = Fi[0]*sqrt(6)*ymax;
	//printf("Umax  = Fi * sqrt(6)*ymax %.8f\n", Umax);
	//printf("Re %f\n", Umax * ymax/nu);
	//printf("Ma %f\n",Umax/cs);	
	//Umax = Fi[0]*ymax*ymax/(8*nu);//Avec body force
	//mu = sqrt(2/(3*M_PI))*Kn*ymax/dx;
	Umax = ymax*ymax*(rho_in-rho_out)*cs*cs/(8*mu*xmax);//Vitesse max pour gradient de pression
	printf("Umax : %.10f\n",Umax);
	Re = Umax * ymax*ratio/nu;
	printf("Re %f\n", Re);
	printf("Ma %f\n",Umax/cs);
	printf("Domaine : [ %.0f %.0f ] x [ %.0f %.0f ]\n",xmin,xmax,ymin,ymax);
	printf("Erreur d'arrêt : %f\n", error);
	printf("Nombre de lattices : %d\n",N);  
	printf("s_nu (entre 0 et 2) : %f\n",1/tau_s);
	printf("s_q (entre 0 et 2): %f\n",1/tau_q);
	/* std:: cout << "  Number of processors available = " << omp_get_num_procs ( ) << "\n"<< std::endl;
	  std ::cout << "  Number of threads =              " << omp_get_max_threads ( ) << "\n" <<std::endl;*/
//************************INITIALISATION MRT*****************************************//
double**  M = MRT_matrice_passage(Q); //Matrice de passage des populations aux moments
double** invM = new double*[Q]; //Matrice de l'inverse de M		
for (int i=0;i<Q;i++)
	{
		invM[i] = new double[Q];
	}
MatrixInversion(M,Q,invM);  
//Matrices avec body force
double** F_bar = new double*[N];
double** F = new double*[N];
for (int i=0;i<N;i++)
	{
		F[i] = new double[Q];
		F_bar[i] = new double[Q];
	}	
				//********SANS WALL FUNCTION*********//
double** Si = MRT_S(Q,tau_s, tau_q); //Matrice de relaxation			
double** C = new double*[Q]; //Matrice du produit matriciel M-1 * Si
double** C1 = new double*[Q]; // Matrice I-0.5*S
double** C2 = new double*[Q]; //Matrice M-1 * (I-0.5*S) = invM * C1
double** C3 = new double*[Q]; //Matrice M-1 * (I-0.5*S)* M = C2 * M		
double** C4	= new double*[Q]; //Matrice du produit matriciel M-1 * Si * M = C * M
for (int i=0;i<Q;i++)
{
	C[i] = new double[Q];
	C1[i] = new double[Q];
	C2[i] = new double[Q];
	C3[i] = new double[Q];
	C4[i] = new double[Q];
}
matrix_product(invM,Si,C,Q);
for (int i=0;i<Q;i++)
{
	for (int j=0;j<Q;j++)
	{
		C1[i][j] = 0;
	}
	C1[i][i] = 1-0.5*Si[i][i];
}
matrix_product(invM,C1,C2,Q);
matrix_product(C2,M,C3,Q);
matrix_product(C,M,C4,Q);
//AFFICHAGE MATRICES
//affichage_matrix(Q,M);	
//affichage_matrix(Q,invM);	
//affichage_matrix(Q,Si);	
//affichage_matrix(Q,C);	
//affichage_matrix(Q,C3);
//affichage_matrix(Q,C4);
				
				//*******AVEC WALL FUNCTION*********//
//Différence sans/avec : tau_s dépend de la lattice avec wall function, chaque lattice 
//a sa propre matrice de relaxation
/*double*** Si_wall = new double**[ny];
double*** C_wall = new double**[ny]; //Matrice du produit matriciel M-1 * Si
double*** C1_wall = new double**[ny]; // Matrice I-0.5*S
double*** C2_wall = new double**[ny]; //Matrice M-1 * (I-0.5*S) = invM * C1
double*** C3_wall = new double**[ny]; //Matrice M-1 * (I-0.5*S)* M = C2 * M	
for (int i=0;i<ny;i++)
{
	Si_wall[i] = MRT_S(Q,rank[i][1], tau_q); // avec tau_q constant
	//Si_wall[i] = MRT_S(Q,rank[i][1], rank2[i][1]); //avec tau_q modifié
	C_wall[i] = new double*[Q];
	C1_wall[i] = new double*[Q];
	C2_wall[i] = new double*[Q];
	C3_wall[i] = new double*[Q];
	for (int j=0;j<Q;j++)
	{
		C_wall[i][j] = new double[Q];
		C1_wall[i][j] = new double[Q];
		C2_wall[i][j] = new double[Q];
		C3_wall[i][j] = new double[Q];
	}
	matrix_product(invM,Si_wall[i],C_wall[i],Q);
	for (int j=0;j<Q;j++)
	{
		C1_wall[i][j] = new double[Q];
		for (int k=0;k<Q;k++)
		{
			C1_wall[i][j][k] = 0;
		}
		C1_wall[i][j][j] = 1-0.5*Si_wall[i][j][j];
	}
	matrix_product(invM,C1_wall[i],C2_wall[i],Q);
	matrix_product(C2_wall[i],M,C3_wall[i],Q);
}*/
/*for (int i=0;i<ny;i++)
{
	printf("Lattice %d, ordonnée %f, PHI %f, tau_s %f, tau_q %f\n",i,rank[i][0],PHI[i],1/Si_wall[i][7][7],1/Si_wall[i][4][4]);
}*/
//*************************INITIALISATION DU DOMAINE******************************************//
for (j=0;j<N;j++)
{	
	lat.rho_[j]=rho0;
	//Initialisation des conditions liées aux positions des lattices
	//Conditions inlet
	if(j%nx==0)
	{
		lat.rho_[j] = rho_in;
		//lat.u_[j][0] = v_in[j/nx][0];
		//lat.u_[j][1] = v_in[j/nx][1];	
	}
	else if((j%nx)==nx-1)
	{
		lat.rho_[j] = rho_out;
		//lat.u_[j][0] = v_out[j/nx][0];
		//lat.u_[j][1] = v_out[j/nx][1];
	}
	velocity(j,D,Q,xi,lat,sigma);
	for (k=0;k<Q;k++)
	{
		S[j][k] = 0;
		lat.f0_[j][k] = omega_i[k]*lat.rho_[j];
		lat.f_[j][k]  = lat.f0_[j][k]; //Initialement, on aura f = fi,eq (avec vitesse nulle)
		f_star[j][k]  = lat.f_[j][k]-1/tau*(lat.f_[j][k]-lat.f0_[j][k]); //SRT
	}
}		

	
//******************************BOUCLE TEMPORELLE******************************************//
while((erreur>error || erreur<-error))	
{		
	for (j=0 ; j<N ; j++) //Pour chaque lattice
    {
		for (k=0;k<Q;k++)
		{
			fi_equilibre (j,k,lat.rho_[j],cs,lat,lat.u_[j],xi,D,Qi,buffer,omega_i,sigma);
		}
		//regularized_BC_v_inlet(j,k,lat,cs,v_in,nx,xi,D,Qi,buffer,omega_i,cas[j],Pi_neq,Q,f_neq, bb,sigma);
		regularized_BC_p_inlet(j,k,lat,cs,rho_in,xi,D,Qi,buffer,omega_i,cas[j],Pi_neq,Q,f_neq,bb,sigma);
		regularized_BC_p_outlet(j,k,lat,cs,rho_out,xi,D,Qi,buffer,omega_i,cas[j],Pi_neq,Q,f_neq,bb,sigma);
		MRT_equilibre_v2(j,lat,M,Q,buffer);
		/*mu = sqrt(2/(3*M_PI))*Kn*ymax/dx;
		tau_s = 0.5+3*mu/lat.rho_[j];
		tau_q = (8-1/tau_s)/(8*(2-1/tau_s));
		Si[4][4] = 1/tau_q;
		Si[6][6] = 1/tau_q;
		Si[7][7] = 1/tau_s;
		Si[8][8] = 1/tau_s;		
		matrix_product(invM,Si,C,Q);*/
		for (k=0;k<Q;k++)
		{
			MRT_moment(j,k,lat,M,Q);
		}
		
		//MRT_forcing(j,Q,cs,omega_i,xi,Fi,F_bar,lat,temp);
		MRT_collision(j,f_star,C,lat,Q,temp);
		//MRT_collision_v2(j,f_star,C4,lat,Q,temp);
		//MRT_forcing_collision(j,f_star,C,lat,dt,F,C3,F_bar,Q,temp);			
		//MRT_forcing(j,k,cs,omega_i,xi,Fi,F_bar,lat,temp);
		//MRT_collision(j,f_star,C_wall[(int)j/nx],lat,Q,dt,temp);
		//MRT_forcing_collision(j,f_star,C_wall[(int)j/nx],lat,dt,F,C3_wall[(int)j/nx],F_bar,Q,temp);	
		//MRT_collision(j,f_star,C_wall[(int)j/nx],lat,Q,temp);
		//DBB_N_BC(j,cas[j],lat,beta,f_star);
		//DBB_S_BC(j,cas[j],lat,beta,f_star);
		//CBBSR_N_BC(j,cas[j],lat,r,f_star);
		//CBBSR_S_BC(j,cas[j],lat,r,f_star);  
		bounceback_solid_BC(nx,j,lat,f_star,conn,typeLat,bb,nombre,pos,cas[j]);
		//pression_in_BC( j,cas[j],lat,xi_r,rho_in);
       // pression_out_BC( j,cas[j],lat,xi_r,rho_out);		
		periodic_NS_BC(j,nx,ny,cas[j],lat,f_star); 		
	    //bounceback_N_BC(j,cas[j],lat,f_star);
        //bounceback_S_BC(j,cas[j],lat,f_star);
			
		propagation(j,lat,f_star, typeLat, conn, bb,Q);
		if(!typeLat[j])
		{
			density(j,Q,lat,sigma); 
			velocity(j,D,Q,xi,lat,sigma);
			//lat.u_[j][0]+=0.5*dt*Fi[0]/lat.rho_[j];
		}
		else
		{
			lat.rho_[j] = 1;
			lat.u_[j][0] = 0;
			lat.u_[j][1] = 0;
		}
		
	}


        //ETAPE DE PROPAGATION
		//#pragma omp parallel for
		/* for (j=0 ; j<N ; j++) //Pour chaque lattice
        {
			
			propagation(j,lat,f_star, typeLat, conn, bb,Q);
        	pression_in_BC( j,cas[j],lat,xi_r,rho_in);
        	pression_out_BC( j,cas[j],lat,xi_r,rho_out);
			//periodic_WE_BC(j,nx,ny,cas[j],lat,f_star);        	
			//DBB_N_BC(j,cas[j],lat,beta,f_star);
        	//DBB_S_BC(j,cas[j],lat,beta,f_star);
        	bounceback_N_BC(j,cas[j],lat,f_star);
        	bounceback_S_BC(j,cas[j],lat,f_star);
			//CBBSR_N_BC(j,cas[j],lat,r,f_star);
			//CBBSR_S_BC(j,cas[j],lat,r,f_star);	
			//MR_N_BC(j,cas[j],lat,sigma1,f_star);
			//MR_S_BC(j,cas[j],lat,sigma1,f_star);
        	
		}
		for (j=0 ; j<N ; j++)
        {
			density(j,Q,lat,sigma); 
			velocity(j,D,Q,xi,lat,sigma);
			//lat.u_[j][0]+=0.5*dt*Fi[0]/lat.rho_[j];
			
			for (k=0;k<Q;k++)
			{
				fi_equilibre (j,k,lat.rho_[j],cs,lat,lat.u_[j],xi,D,Qi,buffer,omega_i,sigma);
			}
			//MRT_equilibre_v2(j,lat,M,Q,temp);
			MRT_moment(j,k,lat,M,Q); //Calcul des moments	
			//MRT_equilibre(j,lat);
			
			//MRT_forcing(j,Q,cs,omega_i,xi,Fi,F_bar,lat,temp);
			//MRT_collision(j,f_star,C,lat,Q,temp);
			MRT_collision_v2(j,f_star,C4,lat,Q,dt,temp);
			//MRT_forcing_collision(j,f_star,C,lat,dt,F,C3,F_bar,Q,temp);		
			
			//MRT_forcing(j,k,cs,omega_i,xi,Fi,F_bar,lat,temp);
			//MRT_collision(j,f_star,C_wall[(int)j/nx],lat,Q,dt,temp);
			//MRT_forcing_collision(j,f_star,C_wall[(int)j/nx],lat,dt,F,C3_wall[(int)j/nx],F_bar,Q,temp);

        }*/
	
	
               // Ecriture des résultats	
        if (it%outputFrequency==0 && it!=0) 
		{
			char name = FileName(Kn);
			//vorticite(Q,nx,ny,lat,cas,dx, typeLat, conn);
			//writeLattice(domain,"LBM",Re,name,it,lat);
			//CALCUL DE LA CONVERGENCE TEMPORELLE
			valeur2 = 0;
			m=0;
			for (j=0;j<N;j++)
			{
				if (sqrt(lat.u_[j][0]*lat.u_[j][0]+lat.u_[j][1]*lat.u_[j][1])>=0.1*sqrt(3))
				{
					m++;
				}
			}
			
			printf("Pourcentage de la zone fluide : %.2f et  dRHO = %.2f et mu = %.2f \n",m/nombre*100,dRHO,mu);	
			for (int j=0;j<N;j++)
			{
				valeur2+=sqrt(lat.u_[j][0]*lat.u_[j][0]+lat.u_[j][1]*lat.u_[j][1]);
			}
			erreur = 1/sqrt(N)*(valeur2-valeur1)/valeur2;
			valeur1 = valeur2;
			printf("Itération : %d\t Erreur  %.10f, Re = %.2f\n", it,erreur,Re);
		}
		it++;
}

	char name = FileName(Kn);
	writeLattice(domain,"LBM",dRHO,name,it,lat);
	time_t timer2 = time(NULL);
	printf("Temps d'exécution : %d s\n",(int)(timer2-timer1));
	printf("Temps d'exécution pour 1000 itérations: %.2f s\n",(double)((int)(timer2-timer1))/it*1000);

    return 0; 
}


