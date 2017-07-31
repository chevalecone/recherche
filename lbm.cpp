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
	double **f_minus = new double*[N]; //Matrice des populations et des populations d'équilibre symétriques et antisymétriques (TRT)
	double **f_plus = new double*[N];
	double **f0_minus = new double*[N];
	double **f0_plus = new double*[N];
	int *cas = new int[N]; //Matrice des cas associées aux lattices (pour les BC du domaine)
	double **position = new double*[N]; //Matrice des coordonnées de chaque milieu de lattice
	int **conn = new int*[N]; //Matrice de connectivité des lattices
	bool* typeLat = new bool[N]; //booléen pour déterminer si le noeud est solide (true) ou fluide (false)
	
	double **xi = new double*[Q];//Tableau des vecteurs des vitesses du modèle
	double omega_i[Q];//Vecteur des poids en fonction de la direction de la vitesse dans la lattice
	
	int* bb = new  int[Q]; //array des populations de bounceback
	bounceback_neighbour(bb,Q); //remplissage de la matrice bb
	
	double **drag_array = new double*[NbSave]; //Tableau pour stocker la traînée ou la portance en fct du temps
	int i,j,k,l,m;

	//Vitesse du mur (pour cavité entraînée)
	double u_w = 0;
	double v_w = 0.012*cs;
	
	//Ratio diamètre cylindre/ diamètre 
	double ratio = 0.125;

	//Cas avec terme de forçage
	double *Fi = new double[D]; //Terme de forçage
	double **S = new double*[N]; //Matrice des Si effets de force gravitationnelle sur chaque direction de la lattice

	for (i=0;i<NbSave;i++)
	{
		drag_array[i] = new double[2];
	 }
	//Matrice des xi
	 for (i=0;i<Q;i++)
	{
		xi[i]=new double[D];
	}
	//Matrice des f_star,S
	for (j=0;j<N;j++)
	{
		f_star[j] = new double[Q];
		f_minus[j] = new double[Q]; 
		f_plus[j] = new double[Q];
		f0_minus[j] = new double[Q];
		f0_plus[j] = new double[Q];
		S[j] = new double[Q]; //Matrice pour l'apport de la force de gravité
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
	Fi[0]=0.0001*cs;
	Fi[1]=0; //Gravité uniquement verticale selon -y

	//Gradient de pression (tout en gardant l'hypothèse incompressible)
	double rho_in = 1.0001*rho0;
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

//******************************CREATION DES SOLIDES******************************************//
/*	//CYLINDRE A SECTION CARREE
	double** cylinder1  = new double*[4];
	for (j=0;j<4;j++)
	{
		cylinder1[j] = new double[D];
	}
	SquareCylinder(0.25*xmax,0.5*ymax,ratio*ymax,cylinder1);
	typeSquare(N,cylinder1,position,typeLat); //Tableau de booléens pour les noeuds fluide/solide

	int *pos = new int[2];
	pos_solide(typeLat, pos,nx,ny);
	printf("Pos min et max des lattices solides : %d et %d\n",pos[0],pos[1]);
*/
//********************************INITIALISATION TRT*********************************//
/*
	double Lambda = 0.25; //Paramètre "magique"
	double tau_plus = nu/(dt*cs*cs)+0.5;
	double tau_minus = 0.5+Lambda/(tau_plus-0.5);
*/
//********************INITIALISATION MRT********************************************//
/*	float**  M = MRT_matrice_passage(Q); //Matrice de passage des populations aux moments
	float** Si = MRT_S(Q,nu,cs,dt); //Matrice des temps de relaxation
	float** invM = new float*[Q]; //Matrice de l'inverse de M
	float** C = new float*[Q]; //Matrice du produit matriciel M-1 * Si
	float** C1 = new float*[Q]; // Matrice I-0.5*S
	float** C2 = new float*[Q]; //Matrice M-1 * (I-0.5*S) = invM * C1
	float** C3 = new float*[Q]; //Matrice M-1 * (I-0.5*S)* M = C2 * M
	//Matrices avec body force
	float** F_bar = new float*[N];
	float** F = new float*[N];
	for (int i=0;i<Q;i++)
	{
		invM[i] = new float[Q];
		C[i] = new float[Q];
		C1[i] = new float[Q];
		C2[i] = new float[Q];
		C3[i] = new float[Q];
	}
	for (int i=0;i<N;i++)
	{
		F[i] = new float[Q];
		F_bar[i] = new float[Q];
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
	affichage_matrix(Q,M,invM,Si,C,F,F2);

	
*/
//*************************INITIALISATION DU DOMAINE******************************************//
for (j=0;j<N;j++)
{	
	//Initialisation des conditions liées aux positions des lattices
	//Conditions inlet
	if((j%nx)==0)
	{
		lat.rho_[j]  = rho0;
		//lat.u_[j][0] = v_in[j/nx][0];
		//lat.u_[j][1] = v_in[j/nx][1];		
	}
	//Conditions outlet
	else if ((j%nx) == nx-1)
	{
		lat.rho_[j]  = rho_out;
	}
	//Vitesse imposée Nord pour cavité entraînée + Couette
	/*else if(j>=nx*ny-nx)
	{
		lat.rho_[j]  = rho0;
		lat.u_[j][0] = v_w;
		lat.u_[j][1] = 0;
	}*/
	
	//Conditions pour les solides
	else if (typeLat[j]) //Sur les lattices solides, on initialise avec une vitesse nulle
	{
		lat.rho_[j]=rho0;
		lat.u_[j][0] = 0;
		lat.u_[j][1] = 0;
	}
	else //Sur les lattices fluides restantes du domaine (autres que inlet/outlet)
	{
		lat.rho_[j]=rho0;
		velocity(j,D,Q,xi,lat);
	}		
	//TRT_creation(j,lat,f_minus,f_plus,f0_plus,f0_minus,bb,Q);		
	for (k=0;k<Q;k++)
	{
		lat.f0_[j][k] = omega_i[k]*lat.rho_[j];
		lat.f_[j][k]= lat.f0_[j][k]; //Initialement, on aura f = fi,eq (avec vitesse nulle)
		f_star[j][k] = lat.f_[j][k]-1/tau*(lat.f_[j][k]-lat.f0_[j][k]); //SRT
		//f_star[j][k] = lat.f_[j][k]-1/tau_plus*(f_plus[j][k]-f0_plus[j][k])-1/tau_minus*(f_minus[j][k]-f0_minus[j][k]); //TRT
		S[j][k] = 0;	
	}
	velocity(j,D,Q,xi,lat);
	//printf("Lattice %d, valeur de la population n° %d : %0.3f\n",j,k,lat.f_[j][i]);
	/*for (k=0;k<Q;k++)
	{
		MRT_moment(j,k,lat,M,Q);
		//printf("Lattice %d, valeur du moment n° %d : %0.3f\n",j,k,lat.m_[j][k]);
		MRT_equilibre(j,k,lat,Q,M);
	}*/
	
}
//Re = v_w*ymax/nu; //Lid driven cavity
//Umax = ymax*ymax*(rho_in-rho_out)*cs*cs/(8*nu*xmax); //Vitesse max pour gradient de pression
Umax = Fi[0]*ymax*ymax/(8*nu);//Avec body force
Re = Umax * ymax/nu;

Mach = Umax/cs;


printf("dt %.8f\n", dt);
printf("cs %.8f\n", cs);
printf("Re %.5f\n", Re);
printf("Umax %.5f\n", Umax);
printf("Nombre de Mach : Ma = %.5f\n",Mach);
printf("Nombre de lattices : %d\n",N);  
//printf("tau_plus = %f\n",tau_plus);
//printf("tau_minus = %f\n",tau_minus);
writeLattice(domain,"LBM_ini",0,lat);


	
//******************************BOUCLE TEMPORELLE******************************************//
while(erreur>error || erreur<-error)
{
	it++;
	//for (it=0; it< ndt; it++)
    //{
        //ETAPE DE PROPAGATION
        for (j=0 ; j<N ; j++) //Pour chaque lattice
        {
   			propagation(j,lat,f_star,nx,ny,cas, typeLat, conn);
        	//pression_in_BC(j,cas[j],lat,xi_r,rho_in);
			//pression_out_BC(j,cas[j],lat,xi_r,rho_out);			
				
			bounceback_N_BC(j,cas[j],lat,f_star);
			bounceback_S_BC(j,cas[j],lat,f_star);
			periodic_WE_BC(j,nx,ny,cas[j],lat,f_star);	
			//vitesse_in_BC(j,nx,cas[j],lat,xi_r,v_in);
			//vitesse_out_BC(j,	nx,cas[j],lat,xi_r,v_out);
			//bounceback_E_BC(j,cas[j],lat,f_star);
			//bounceback_W_BC(j,cas[j],lat,f_star);
			//driven_cavity_nord(j,cas[j],lat,xi_r,v_w);
			//periodic_NS_BC(j,nx,ny,cas[j],lat,f_star);
			//bounceback_solid_BC(nx,j, lat, f_star,conn, typeLat, bb, nombre, pos);
		}
		//ETAPE DE COLLISION
    	for (j=0;j<N;j++)
    	{
    		//ETAPE DE CALCUL DES FORCES GRAVITATIONNELLES, VITESSE ET DENSITE
			//Calcul de u dans le cas où il n'y a pas de force à distance
			if (!typeLat[j])
			{
				density(j,Q,lat,sigma);
				velocity(j,D,Q,xi,lat);
				lat.u_[j][0]+=0.5*Fi[0]*dt/lat.rho_[j]; //SRT + body force
				lat.u_[j][1]+=0.5*Fi[1]*dt/lat.rho_[j];

			}
			else //Pour les lattices solides, on impose les variables macroscopiques pour visualiser
			{
				lat.rho_[j]=rho0;
				lat.u_[j][0]=0;
				lat.u_[j][1]=0;
			}
    		//**********************SRT***************************//
			for (k=0;k<Q;k++)
			{	
				S[j][k]=(1-1/(2*tau))*omega_i[k]*((xi[k][0]-lat.u_[j][0])/(cs*cs)+pscal(xi[k],lat.u_[j],D)/(cs*cs*cs*cs)*xi[k][0])*Fi[0];
				lat.f0_[j][k]=omega_i[k]*lat.rho_[j]*(1+1/(cs*cs)*pscal(xi[k],lat.u_[j],D)+1/(2*cs*cs*cs*cs)*pscal(xi[k],lat.u_[j],D)*pscal(xi[k],lat.u_[j],D)-1/(2*cs*cs)*pscal(lat.u_[j],lat.u_[j],D));
				f_star[j][k]=lat.f_[j][k]-1/tau*(lat.f_[j][k]-lat.f0_[j][k])+S[j][k]*dt; //SRT + body force
				//f_star[j][k]=lat.f_[j][k]-1/tau*(lat.f_[j][k]-lat.f0_[j][k]); //SRT
			}
			//***************************TRT************************//
			/*TRT_creation(j,lat,f_minus,f_plus,f0_plus,f0_minus,bb,Q);	
			for (k=0;k<Q;k++)
			{
				f_star[j][k] = lat.f_[j][k]-1/tau_plus*(f_plus[j][k]-f0_plus[j][k])-1/tau_minus*dt*(f_minus[j][k]-f0_minus[j][k]); //TRT
				//f_star[j][k] = lat.f_[j][k]-1/tau_plus*(f_plus[j][k]-f0_plus[j][k])-1/tau_minus*dt*(f_minus[j][k]-f0_minus[j][k]) + S[j][k]*dt; //TRT + body force
			}*/
			//**********************MRT***************************//
			/*MRT_equilibre(j,k,lat,Q,M);
			for (k=0;k<Q;k++)
			{
				MRT_moment(j,k,lat,M,Q);
				MRT_forcing(j,k,cs,omega_i,xi,Fi,F_bar,lat,D);	
			}
				MRT_collision(j,f_star,C,lat,Q,dt);
				MRT_forcing_collision(j,f_star,C,lat,Q,dt,F,C3,F_bar);
*/
    	}
        // When you want write the results
        if (it%outputFrequency==0) 
		{
			time_t timer3 = time(NULL);
			vorticite(Q,nx,ny,lat,cas,dx, typeLat, conn);
			writeLattice(domain,"LBM",it,lat);
		    time_t timer4 = time(NULL);
		    //df = drag_force(Q,nx,ny,f_star,xi_r,typeLat,conn,xi,bb,pos, lat);

		   //drag_array[it/outputFrequency][0] = it*dt;
		   //drag_array[it/outputFrequency][1] = df;

			//Calcul de la convergence temporelle
			valeur2 = 0;
			for (j =0;j<N;j++)
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
			//printf("Force de traînée %.8f\n",df);
			//printf("Coefficient de traînée : %.8f\n",df*4/(0.5*rho0*ymax*ratio*Umax*Umax*9));
		}
	}
	time_t timer2 = time(NULL);
	//drag_array[NbSave-1][0] = it*dt;
	//drag_array[NbSave-1][1] = drag_force(Q,nx,ny,f_star,xi_r,typeLat,conn,xi,bb,pos,lat);
	//writeTimeScalar("drag_force",NbSave,drag_array);
    vorticite(Q,nx,ny,lat,cas,dx, typeLat, conn);
    writeLattice(domain,"LBM",it,lat);
	printf("Temps d'exécution : %d s\n",(int)(timer2-timer1));
	printf("Temps d'écriture du fichier : %d s\n",(int)w_time/(NbSave));
	printf("Temps réel de calcul : %d s\n",(int)(timer2-timer1-w_time));
	printf("Temps d'exécution pour 1000 itérations: %.2f s\n",(double)((int)(timer2-timer1-w_time))/ndt*1000);

    return 0; 
}


