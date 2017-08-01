/*******************************************************************
*
*
********************************************************************/
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
#include "geometry.h"
#include "boundaryc.h"
#include "compute_quantities.h"
#include "relaxation_time.h"

int main()
{
	//*******************************CREATION DU DOMAINE A L'AIDE DU FICHIER XML******************************//
    Parser parser;
    parser.parse();

    Domain domain;

    // Maillage
    double dx= parser.getDx();
    domain.setDx(dx);

    // Set the size of the physical domain
    double xmin = parser.getXmin();
    double xmax = parser.getXmax();
    double ymin = parser.getYmin();
    double ymax = parser.getYmax();
    domain.setXMin(xmin);
    domain.setXMax(xmax);
    domain.setYMin(ymin);
    domain.setYMax(ymax);

    // Finalize creation of the domain
    domain.setFinalize();

    //Nombre de points de discrétisation en espace
    int nx = domain.getNx();
    int ny = domain.getNy();
    int N = domain.getNTot(); //Nombre de lattices dans le milieu

    // Propriétés physiques
    double tau = parser.getTau(); //Temps de relaxation

    double mu = parser.getMu(); //Viscosité dynamique du fluide
    double rho0 = parser.getRho0(); //Masse volumique du fluide	
    double nu = mu/rho0; //Viscosité cinématique du fluide
    printf("nu = %f\n",nu);
     tau = 0.5+3*nu/dx;
      printf("tau = %f\n",tau);

	// Information sur la boucle temporelle
	int it = 0;
    //double dt=dx; //Calcul du pas de temps;
    double dt = dx;
    double timeEnd = parser.getTimeEnd(); //Temps total de simulation
	double error = parser.getError();
    int outputFrequency = parser.getOutputFrequency(); //Fréquence d'enregistrement des valeurs
    int ndt = int (timeEnd / dt);// Nombre de points de discrétisation en temps
    int NbSave = ndt/outputFrequency+1; //Nombre de fichiers sauvés pour l'évolution temporelle des variables


    // Données générales sur le modèle géométrique utilisé
    int Q = 9; //Nombre de discrétisation de vitesse dans une lattice (ex : Q=9 pour D2Q9)
	int D = 2; //Nombre de dimensions (ex : D =2 pour 2D)
    Lattice lat(domain.getNTot(),Q,D);
	double xi_r = dx/dt; 

	time_t timer1 = time(NULL); //Temps départ
	
	//********************************ALLOCATION DE LA MEMOIRE*************************************//

	//Célérité de la lattice, Re, Vitesse max, variables pour la convergence en temps, force de drag et lift	
	double cs = 1/sqrt(3)*xi_r,Re,Umax,Umax2,valeur1=0,valeur2=0,erreur=1,df,lf,sigma,nombre,w_time=0,Mach;
	double **f_star = new double*[N]; //Matrice temporaire f* entre t et t+dt (collisions)
	int *cas = new int[N]; //Matrice des cas associées aux lattices (pour les BC du domaine)
	double **position = new double*[N]; //Matrice des coordonnées de chaque milieu de lattice
	int **conn = new int*[N]; //Matrice de connectivité des lattices
	bool* typeLat = new bool[N]; //booléen pour déterminer si le noeud est solide (true) ou fluide (false)
	
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
	}
	//Matrice des f_star,S
	for (j=0;j<N;j++)
	{
		f_star[j] = new double[Q];
		typeLat[j] = false; //Initialement, tous les noeuds sont fluides
		position[j] = new double[D];
		conn[j] = new int[Q];
		cas[j] =0; // Initialisation des cas : tous des points intérieurs
	}
//***********************INITIALISATION DE LA GEOMETRIE DE LATTICE***************************//
	D2Q9(omega_i,xi,xi_r); //Remplissage des poids omega_i et des vecteurs vitesses en fonction de la géométrie considérée
	connectivite(nx,ny,Q,conn);// Matrice de connectivité
//***********************************POSITION DES LATTICES***********************************//
	localisation(nx,ny,dx,position); //Coordonnées de chaque noeud
	domainCondition(nx,ny,cas);	//Cas pour les BC aux frontières du domaine
//******************* PHENOMENES PHYSIQUES MOTEURS DE L'ECOULEMENT **************************//
	//Initialisation du terme de forçage (dans le cas de conditions périodiques)
	Fi[0]=0.000001*cs;
	Fi[1]=0; //Gravité uniquement verticale selon -y

	//Gradient de pression (tout en gardant l'hypothèse incompressible)
	double rho_in = 1.00005*rho0;
	double rho_out = rho0;

	//Vitesse imposée (tout en gardant l'hypothèse incompressible)
	Umax = 0.1*cs;
	Umax2 = 0.1*cs;

	double** v_in = new double*[ny];
	double** v_out = new double*[ny];

	//Création profil de vitesse inlet de Poiseuille selon l'axe x
	for (j=0;j<ny;j++)
	{
		v_in[j] = new double[2];
		v_out[j] = new double[2];
		v_out[j][0] = Umax2*(4*position[j*nx][1]/ymax-4*position[j*nx][1]*position[j*nx][1]/(ymax*ymax)); 
		v_out[j][1] = 0;
		v_in[j][0] = Umax*(4*position[j*nx][1]/ymax-4*position[j*nx][1]*position[j*nx][1]/(ymax*ymax)); 
		v_in[j][1] = 0;
	}
//********************INITIALISATION MRT********************************************//
	double**  M = MRT_matrice_passage(Q); //Matrice de passage des populations aux moments
	double** Si = MRT_S(Q,nu,cs,dt); //Matrice de relaxation
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
	affichage_matrix(Q,M,invM,Si,C,C3);


//*************************INITIALISATION DU DOMAINE******************************************//
for (j=0;j<N;j++)
{	
	//Initialisation des conditions liées aux positions des lattices
	//Conditions inlet
	lat.rho_[j]=rho0;
	velocity(j,D,Q,xi,lat);
	for (k=0;k<Q;k++)
	{
		lat.f0_[j][k] = omega_i[k]*lat.rho_[j];
		lat.f_[j][k]= lat.f0_[j][k]; //Initialement, on aura f = fi,eq (avec vitesse nulle)
		f_star[j][k] = lat.f_[j][k]-1/tau*(lat.f_[j][k]-lat.f0_[j][k]); //SRT
	}
	MRT_equilibre(j,lat,M);
	MRT_moment(j,lat,M);
	/*for (k=0;k<Q;k++)
	{	
		//printf("Lattice %d, valeur du moment n° %d : %0.3f\n",j,k,lat.m_[j][k]);
		MRT_forcing(j,k,cs,omega_i,xi,Fi,F_bar,lat);
	}
	//MRT_collision(j,f_star,C,lat,Q,dt);
	MRT_forcing_collision(j,f_star,C,lat,dt,F,C3,F_bar);*/
}		

Umax = Fi[0]*ymax*ymax/(8*nu);//Avec body force
Re = Umax * ymax/nu;

Mach = Umax/cs;


printf("dt %.8f\n", dt);
printf("cs %.8f\n", cs);
printf("Re %.5f\n", Re);
printf("Umax %.5f\n", Umax);
printf("Ma %.5f\n",Mach);
printf("Nombre de lattices : %d\n",N);  
writeLattice(domain,"LBM_ini",0,lat);


	
//******************************BOUCLE TEMPORELLE******************************************//
while((erreur>error || erreur<-error))	
{
        //ETAPE DE PROPAGATION
        for (j=0 ; j<N ; j++) //Pour chaque lattice
        {
   			propagation(j,lat,f_star,nx,ny,cas, typeLat, conn);
			bounceback_N_BC(j,cas[j],lat,f_star);
			bounceback_S_BC(j,cas[j],lat,f_star);
			periodic_WE_BC(j,nx,ny,cas[j],lat,f_star);
			density(j,Q,lat); 
			velocity(j,D,Q,xi,lat);
			MRT_equilibre(j,lat,M); //Calcul des moments d'équilibre
			MRT_moment(j,lat,M); //Calcul des moments
			MRT_forcing(j,Q,cs,omega_i,xi,Fi,F_bar,lat);
			//ETAPE DE COLLISION
			//MRT_collision(j,f_star,C,lat,Q,dt);
			MRT_forcing_collision(j,f_star,C,lat,dt,F,C3,F_bar);
		}

        // Ecriture des résultats
        if (it%outputFrequency==0 && it!=0) 
		{
			time_t timer3 = time(NULL);
			vorticite(Q,nx,ny,lat,cas,dx, typeLat, conn);
			writeLattice(domain,"LBM",it,lat);
		    time_t timer4 = time(NULL);

			//CALCUL DE LA CONVERGENCE TEMPORELLE
			valeur2 = 0;
			for (int j=0;j<N;j++)
			{
				valeur2+=sqrt(lat.u_[j][0]*lat.u_[j][0]+lat.u_[j][1]*lat.u_[j][1]);
			}
			erreur = 1/sqrt(N)*(valeur2-valeur1)/valeur2;
			if(w_time==0 && it!=0)
			{
				w_time = NbSave*(timer4-timer3);
			}
			valeur1 = valeur2;
			printf("Itération : %d\t Erreur  %.10f\n", it,erreur);
		}
		it++;
}
	time_t timer2 = time(NULL);
    vorticite(Q,nx,ny,lat,cas,dx, typeLat, conn);
    writeLattice(domain,"LBM",it,lat);
	printf("Temps d'exécution : %d s\n",(int)(timer2-timer1));
	printf("Temps d'écriture du fichier : %d s\n",(int)w_time/(NbSave));
	printf("Temps réel de calcul : %d s\n",(int)(timer2-timer1-w_time));
	printf("Temps d'exécution pour 1000 itérations: %.2f s\n",(double)((int)(timer2-timer1-w_time))/ndt*1000);

    return 0; 
}


