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

// Local includes
#include "domain.h"
#include "lattice.h"
#include "solutionExporterEulerian.h"
#include "function.h"
#include "geometry.h"
#include "boundaryc.h"
#include "compute_quantities.h"
#include "relaxation_time.h"
#include "rarefied_models.h"

# define M_PI  3.14159265358979323846

//Variables de la simulation
# define RHO 1
# define MU 0.025
# define DX 1
# define dRHO 1.00001
# define dFi 0.0000001
# define R 1
# define BETA 1
# define KNU 0.388
# define XMIN 0
# define XMAX 5
# define YMIN 0
# define YMAX 100
# define OUTPUT 10000
# define PRECISION 0.0000000001

int main(int argc, char *argv[])
{
	//*******************************CREATION DU DOMAINE A L'AIDE DU FICHIER XML******************************//
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

    
    double nu = mu/rho0; //Viscosité cinématique du fluide
    
     tau = 0.5+3*nu/dx;
  

	// Information sur la boucle temporelle
	int it = 0;
    //double dt=dx; //Calcul du pas de temps;
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
	double cs = 1/sqrt(3)*xi_r,Re,Umax,Umax2,valeur1=0,valeur2=0,erreur=1,df,lf,sigma,nombre,w_time=0,Mach;
	double* PHI = new double[N];
	double* temp = new double[2];
	double* temp2 = new double[3];
	double **f_star = new double*[N]; //Matrice temporaire f* entre t et t+dt (collisions)
	int *cas = new int[N]; //Matrice des cas associées aux lattices (pour les BC du domaine)
	double **position = new double*[N]; //Matrice des coordonnées de chaque milieu de lattice
	int **conn = new int*[N]; //Matrice de connectivité des lattices
	bool* typeLat = new bool[N]; //booléen pour déterminer si le noeud est solide (true) ou fluide (false)
	double*** Qi = new double**[Q];
	double** S = new double*[N];
	
	double **xi = new double*[Q];//Tableau des vecteurs des vitesses du modèle
	double omega_i[Q];//Vecteur des poids en fonction de la direction de la vitesse dans la lattice
	
	int* bb = new  int[Q]; //array des populations de bounceback
	bounceback_neighbour(bb,Q); //remplissage de la matrice bb
	
	int i,j,k,l,m;


	//Cas avec terme de forçage
	double *Fi = new double[D]; //Terme de forçage

	//Matrice des xi
	 for (i=0;i<Q;i++)
	{
		xi[i]=new double[D];
		Qi[i] = new double*[D];
		for (int j=0;j<D;j++)
		{
			Qi[i][j] = new double[D];
		}
	}
	//Matrice des f_star,S
	for (j=0;j<N;j++)
	{
		f_star[j] = new double[Q];
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
	Qi_equilibre(Q,cs,xi,Qi);


//******************* PHENOMENES PHYSIQUES MOTEURS DE L'ECOULEMENT **************************//
	//Initialisation du terme de forçage (dans le cas de conditions périodiques)
	Fi[0]=0.0000001;
	Fi[1]=0; //Gravité uniquement verticale selon -y

	//Gradient de pression (tout en gardant l'hypothèse incompressible)
	double rho_in = dRHO*rho0;
	double rho_out = rho0;

	//Vitesse imposée (tout en gardant l'hypothèse incompressible)
	Umax = 0.1*cs;
	Umax2 = 0.1*cs;

	double** v_in = new double*[ny];
	double** v_out = new double*[ny];

	//Inlet velocity
	for (j=0;j<ny;j++)
	{
		v_in[j] = new double[2];
		v_out[j] = new double[2];
		v_out[j][0] = Umax2*(4*position[j*nx][1]/ymax-4*position[j*nx][1]*position[j*nx][1]/(ymax*ymax)); 
		v_out[j][1] = 0;
		v_in[j][0] = Umax*(4*position[j*nx][1]/ymax-4*position[j*nx][1]*position[j*nx][1]/(ymax*ymax)); 
		v_in[j][1] = 0;
	}

//****************************ECOULEMENTS RAREFIES**********************************//
	double Kn = KNU;
	double tau_s = temp2[1];
	double tau_q = temp2[2];
	double r = temp2[0];
	double beta = temp2[0];

	//temp2 = MRT_Continuous(cs,dt,nu);


	/*temp2 = MRT_Guo_2008(Kn,ymax,tau_s,tau_q,r,dx);
	r = temp2[0];
	printf("r : %f\n",r);
	PHI_Guo_2008(Kn,ymax,N,dx,PHI,position, mu, rho0, cs);
	for (int i =0;i<N;i++)
	{
		printf("Lattice %d, ordonnée %f,  PHI : %.3f\n",i,position[i][1],PHI[i]);
	}*/
	
	//mu = sqrt(2/(3*M_PI))*Kn*ymax/dx;
	nu = mu/rho0;
	temp2 = MRT_Verhaeghe_2009(Kn,ymax,tau_s,tau_q,beta,mu,dx,dt,cs);
	beta = temp2[0];
	printf("beta : %f\n",beta);

	tau_s = temp2[1];
	tau_q = temp2[2];
	
	printf("tau_s : %f\n",tau_s);
	printf("tau_q : %f\n",tau_q);
	

	Umax = Fi[0]*ymax*ymax/(8*nu);//Avec body force
	//Umax = ymax*ymax*(rho_in-rho_out)*cs*cs/(8*nu*xmax); //Vitesse max pour gradient de pression
	Re = Umax * ymax/nu;
	Mach = Umax/cs;
	printf("Kn = %f\n",Kn);
	printf("Domaine : [ %.0f %.0f ] x [ %.0f %.0f ]\n",xmin,xmax,ymin,ymax);
	printf("Erreur d'arrêt : %f\n", error);
	printf("nu = %f\n",nu);
	printf("Re %f\n", Re);
	printf("Umax %f\n", Umax);
	printf("Ma %f\n",Mach);
	printf("Nombre de lattices : %d\n",N);  
	printf("s_nu (entre 0 et 2) : %f\n",1/tau_s);
	printf("s_q (entre 0 et 2): %f\n",1/tau_q);
	writeLattice(domain,"LBM_ini",0,lat);
	/* std:: cout << "  Number of processors available = " << omp_get_num_procs ( ) << "\n"<< std::endl;
	  std ::cout << "  Number of threads =              " << omp_get_max_threads ( ) << "\n" <<std::endl;*/
//********************INITIALISATION MRT********************************************//
	double**  M = MRT_matrice_passage(Q); //Matrice de passage des populations aux moments
	double** Si = MRT_S(Q,nu,cs,dt,tau_s, tau_q); //Matrice de relaxation
	double** invM = new double*[Q]; //Matrice de l'inverse de M
	double** C = new double*[Q]; //Matrice du produit matriciel M-1 * Si
	double** C1 = new double*[Q]; // Matrice I-0.5*S
	double** C2 = new double*[Q]; //Matrice M-1 * (I-0.5*S) = invM * C1
	double** C3 = new double*[Q]; //Matrice M-1 * (I-0.5*S)* M = C2 * M
	//Matrices avec body force
	double** F_bar = new double*[N];
	double** F = new double*[N];
	for (int i=0;i<Q;i++)
	{
		invM[i] = new double[Q];
		C[i] = new double[Q];
		C1[i] = new double[Q];
		C2[i] = new double[Q];
		C3[i] = new double[Q];
	}
	for (int i=0;i<N;i++)
	{
		F[i] = new double[Q];
		F_bar[i] = new double[Q];
	}
	MatrixInversion(M,Q,invM);
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

	//AFFICHAGE MATRICES
	//affichage_matrix(Q,M,invM,Si,C,C3);


//*************************INITIALISATION DU DOMAINE******************************************//
for (j=0;j<N;j++)
{	
	lat.rho_[j]=rho0;
	//Initialisation des conditions liées aux positions des lattices
	//Conditions inlet
	if(j%nx==0)
	{
		lat.rho_[j] = rho0;
		//lat.u_[j][0] = v_in[j/nx][0];
		//lat.u_[j][1] = v_in[j/nx][1];	
	}
	else if(j%nx==nx-1)
	{
		lat.rho_[j] = rho_out;
	}
	//velocity(j,D,Q,xi,lat,sigma);
	for (k=0;k<Q;k++)
	{
		S[j][k] = 0;
		//lat.f0_[j][k] = omega_i[k]*lat.rho_[j];
		lat.f0_[j][k] = omega_i[k]*lat.rho_[j]*(1+1/(cs*cs)*pscal(xi[k],lat.u_[j],D,sigma)+1/(2*cs*cs*cs*cs)*pscal(xi[k],lat.u_[j],D,sigma)*pscal(xi[k],lat.u_[j],D,sigma)-1/(2*cs*cs)*pscal(lat.u_[j],lat.u_[j],D,sigma));
		lat.f_[j][k]  = lat.f0_[j][k]; //Initialement, on aura f = fi,eq (avec vitesse nulle)
		f_star[j][k]  = lat.f_[j][k]-1/tau*(lat.f_[j][k]-lat.f0_[j][k]); //SRT
	}
	MRT_equilibre(j,lat,M);
	MRT_moment(j,lat,M,Q,temp);
	MRT_forcing(j,k,cs,omega_i,xi,Fi,F_bar,lat,temp);
	//MRT_collision(j,f_star,C,lat,Q,dt,temp);
	MRT_forcing_collision(j,f_star,C,lat,dt,F,C3,F_bar,Q,temp);
}		




	
//******************************BOUCLE TEMPORELLE******************************************//
while((erreur>error || erreur<-error))	
{		
        //ETAPE DE PROPAGATION
		//#pragma omp parallel for
		 for (j=0 ; j<N ; j++) //Pour chaque lattice
        {
        	propagation(j,lat,f_star,nx,ny,cas, typeLat, conn, bb);
        	//pression_in_BC( j,cas[j],lat,xi_r,rho_in);
        	//pression_out_BC( j,cas[j],lat,xi_r,rho_out);
        	//equilibrium_inlet_BC(j,cas[j],lat,rho_in,cs,omega_i,xi,Q,nx);
        	//equilibrium_outlet_BC(j,cas[j],lat,rho_out,cs,omega_i,xi,Q);
        	periodic_WE_BC(j,nx,ny,cas[j],lat,f_star);
        	DBB_N_BC(j,cas[j],lat,r,f_star);
        	DBB_S_BC(j,cas[j],lat,r,f_star);

        	//bounceback_N_BC(j,cas[j],lat,f_star);
        	//bounceback_S_BC(j,cas[j],lat,f_star);
		/*}	
       	for (j=0 ; j<N ; j++)
        {
        	CBBSR_N_BC(j,cas[j],lat,r,f_star);
			CBBSR_S_BC(j,cas[j],lat,r,f_star);	
		}
		for (j=0 ; j<N ; j++)
        {*/
			
			density(j,Q,lat,sigma); 
			velocity(j,D,Q,xi,lat,sigma);
			MRT_equilibre(j,lat,M); //Calcul des moments d'équilibre
			MRT_moment(j,lat,M,Q,temp); //Calcul des moments
			MRT_forcing(j,Q,cs,omega_i,xi,Fi,F_bar,lat,temp);
			//ETAPE DE COLLISION
			//MRT_collision(j,f_star,C,lat,Q,dt,temp);

			MRT_forcing_collision(j,f_star,C,lat,dt,F,C3,F_bar,Q,temp);
        }
	
	
               // Ecriture des résultats	
        if (it%outputFrequency==0 /*&& it!=0*/) 
		{
			//vorticite(Q,nx,ny,lat,cas,dx, typeLat, conn);
			writeLattice(domain,"LBM",it,lat);
			//CALCUL DE LA CONVERGENCE TEMPORELLE
			valeur2 = 0;
			for (int j=0;j<N;j++)
			{
				valeur2+=sqrt(lat.u_[j][0]*lat.u_[j][0]+lat.u_[j][1]*lat.u_[j][1]);
			}
			erreur = 1/sqrt(N)*(valeur2-valeur1)/valeur2;
			valeur1 = valeur2;
			printf("Itération : %d\t Erreur  %.10f\n", it,erreur);
		}
		it++;
}
	time_t timer2 = time(NULL);
   // vorticite(Q,nx,ny,lat,cas,dx, typeLat, conn);
    writeLattice(domain,"LBM_nu_0.025_r_0.95",it,lat);
	printf("Temps d'exécution : %d s\n",(int)(timer2-timer1));
	printf("Temps d'exécution pour 1000 itérations: %.2f s\n",(double)((int)(timer2-timer1))/it*1000);

    return 0; 
}


