/*******************************************************************
*
*
********************************************************************/
// C++ Standard includes
#include <iostream> 
#include <cmath>
#include <ctime>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <omp.h>
#include <cstring>

// Local includes
#include "domain.h"
#include "lattice.h"
#include "solutionExporterEulerian.h"
#include "function.h"
#include "geometry.h"
#include "boundaryc.h"
#include "relaxation_time.h"
#include "rarefied_models.h"
#include "compute_quantities.h"
#include "slip_velocity.h"
#include "wall_function.h"
#include "initialisation.h"
#include "solid_interpolation_method.h"
# define M_PI  3.14159265358979323846

//Variables de la simulation
# define RHO 1
# define MU 0.5
# define DX 1

# define dRHO 1.2
# define dFi 0.00000001
# define U_IN_x 0.0
# define U_IN_y 0.
# define U_OUT_x 0.
# define U_OUT_y 0.
# define UW_X 0
# define UW_Y 0
# define R 1
# define BETA 1
# define PORO 0.6
# define ASPECT_RATIO 0.65
// 0.0564 0.1128 0.1692 0.2257 0.3385 0.4514     0.6670 0.9027 1.1284 1.6926 2.2568 3.3851    4.5135 6.7703 9.0270 11.2838 16.9257 
//0.0564 0.1692 0.3385 1.6926 3.3851
# define KNU  0.1128
# define XMIN 0
# define XMAX 250
# define YMIN 0
# define YMAX 40
# define OUTPUT 100
# define PRECISION 0.00000005


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
	double cs = 1/sqrt(3)*xi_r,Re,Umax,valeur1=0,valeur2=0,erreur=1,df,lf,sigma,w_time=0,Mach, poro, ttsity, Kne, Knee, lambdae, Fpc;
	int i,j,k,l,m;
	double* buffer = new double[10];
	double* t = new double[Q];
	double* teq = new double[Q];
	double* buffer2 = new double[10];
	double* tab_marquage = new double[N];
	double* test = new double[5];
	double** rank = new double*[ny];
	int* tab_Voisin = new int[N];
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
	bool* typeLat2 = new bool[N];
	double** S = new double*[N];
	double **xi = new double*[Q];//Tableau des vecteurs des vitesses du modèle
	double** Pi_neq = new double*[D];
	double omega_i[Q];//Vecteur des poids en fonction de la direction de la vitesse dans la lattice
	double* Uw = new double[2]; //Vecteur de vitesse de la frontière NORD (pour écoulement de Couette)
	 Uw[0] = UW_X;
	 Uw[1] = UW_Y;
	int* bb = new  int[Q]; //array des populations de bounceback
	bounceback_neighbour(bb,Q); //remplissage de la matrice bb
	
	double *Fi = new double[D]; //Body force à N dimensions
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
	for (j=0;j<N;j++)
	{
		tab_marquage[j] = 0;
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
	double ratio2 = 0.001;
	for (j=0;j<ny;j++)
	{
		v_in[j] = new double[2]; //vecteur vitesse d'entrée
		v_out[j] = new double[2];// vecteur vitesse de sortie
		
		//Créations profils de vitesse type Poiseuille
		Umaxx= ratio2*cs;
		v_in[j][0] = Umaxx*(4*position[j*nx][1]/ymax-4*position[j*nx][1]*position[j*nx][1]/(ymax*ymax)); 
		v_in[j][1] = 0;
		v_out[j][0] = Umax2x*(4*position[j*nx][1]/ymax-4*position[j*nx][1]*position[j*nx][1]/(ymax*ymax)); 
		v_out[j][1] = 0;	
	}

//******************************CREATION DES SOLIDES******************************************//
	double ratio = 0.6,nombre =0;
	double** solid_fraction_interpolation = new double*[N];
	for (j=0;j<N;j++)
	{
		solid_fraction_interpolation[j] = new double[Q-1];
	}
	
	/*double perimeter_ellipse,area_ellipse,aspect_ratio_ellipse,hydraulic_diameter_ellipse,h_ellipse;
	aspect_ratio_ellipse = ASPECT_RATIO;
	double b_ellipse = 0.075*nx;
	double a_ellipse = b_ellipse/aspect_ratio_ellipse;
	h_ellipse = (a_ellipse-b_ellipse)*(a_ellipse-b_ellipse)/((a_ellipse+b_ellipse)*(a_ellipse+b_ellipse));
	area_ellipse = M_PI*a_ellipse*b_ellipse;
	perimeter_ellipse = M_PI*(a_ellipse+b_ellipse)*(1+3*h_ellipse/(10+sqrt(4-3*h_ellipse)));
	hydraulic_diameter_ellipse = 4*area_ellipse/perimeter_ellipse;
	printf("Paramètres du milieu poreux aléatoire de cylindres à section elliptique : \n");
	printf("Demi grand axe : %f\n",a_ellipse);
	printf("Demi petit axe : %f\n",b_ellipse);
	printf("Aspect ratio : %f\n",aspect_ratio_ellipse);
	printf("Aire : %f\n",area_ellipse);
	printf("Périmètre : %f\n",perimeter_ellipse);
	printf("Diamètre hydraulique : %f\n",hydraulic_diameter_ellipse);*/
	
	//poro = PORO;
	//Milieu poreux aléatoire avec des cylindres à section cylindriques
	//randomCircular(nx,ny,xmin,xmax,ymin,ymax,N,position,typeLat,poro,nombre, cas);
	//Milieu poreux aléatoire avec des cylindres à section elliptique
	//randomEllipse(nx,ny,xmin,xmax,ymin,ymax,N,position,typeLat,poro,nombre,cas,a_ellipse,b_ellipse);

	
	
	//Milieu poreux aléatoire avec des cylindres à section carré
	//randomSquare(nx,ny,xmin,xmax,ymin,ymax,N,position,typeLat,poro, nombre, cylinder1);
	//typeSquare (0.48*xmax, 0.48*ymax, 0.188*ymax, 0.188*ymax, N, position,typeLat);
	//solid_fraction_square(N, Q, solid_fraction_interpolation,conn, 0.48*xmax, 0.48*ymax, 0.188*ymax, typeLat, buffer, position);
	//typeEllipse(0.5*xmax, 0.5*ymax, 0.1*ymax, 0.2*ymax, 45, N, position, typeLat);
	//typeCircular(0,0,ratio*ymax,N,position,typeLat, typeLat2);	
	//solid_fraction_circular(N, Q, solid_fraction_interpolation,conn, 0,0, ratio*ymax, typeLat2, position);
	//typeCircular(xmax,0,ratio*ymax,N,position,typeLat,typeLat2);
	//solid_fraction_circular(N, Q, solid_fraction_interpolation,conn, xmax, 0, ratio*ymax, typeLat2, position);
	//typeCircular(0,ymax,ratio*ymax,N,position,typeLat,typeLat2);
	//solid_fraction_circular(N, Q, solid_fraction_interpolation,conn, 0, ymax, ratio*ymax, typeLat2, position);
	//typeCircular(xmax,ymax,ratio*ymax,N,position,typeLat,typeLat2);
	//solid_fraction_circular(N, Q, solid_fraction_interpolation,conn, xmax,ymax, ratio*ymax, typeLat2, position);
	typeCircular(0.5*xmax,0.5*ymax,ratio*ymax,N,position,typeLat,typeLat2);
	solid_fraction_circular(N, Q, solid_fraction_interpolation,conn, 0.5*xmax, 0.5*ymax, ratio*ymax, typeLat2, position);


	for (j=0;j<N;j++)
	{
		for (k=0;k<Q;k++)
		{
				lat.q_[j][k] =solid_fraction_interpolation[j][k];
		} 
	}
	poro = porosite(typeLat,nombre,N);
	//printf("poro : %f\n",poro);
	//printf("Diametre : %f\n",100*ratio);
	
	int *pos = new int[2];
	pos_solide(typeLat, pos,nx,ny);

	//tab_voisin(N,Q,typeLat,tab_Voisin,conn);
	/*for (j=0;j<N;j++)
	{
		printf("Tab_Voisin %d : %d\n",j,tab_Voisin[j]);
	}*/
	nettoyage(typeLat,conn,N,Q);
	poro = porosite (typeLat,nombre,N);
	printf("Porosité : %f\n",poro);
	
//****************************ECOULEMENTS RAREFIES**********************************//

	double Kn = KNU;
	double tau_s = temp2[1];
	double tau_q = temp2[2];
	double r = temp2[0];
	double beta = temp2[0];
	double sigma1 = temp2[0];
	char name = FileName(Kn);

//*************************SLIP VELOCITY***************************//

	//Coefficients des vitesses de glissements considérés
	sigma = 1;
	// 'Guo-2008' , 'Guo-2011' , 'Wang-2017' , 'Hadjiconstantinou-2003' , 'Li-2011' , 'Wu-2008' , 'de Izarra-2012'
	std::string slip ("Guo-2008");
	double* A = Slip_velocity(slip,sigma,Kn);
	printf("Coefficients de glissement : A1 = %f , A2 = %f\n",A[0],A[1]);
	
//************************WALL FUNCTION***************************//

	//Définition de la wall function considérée
	// 'Guo-2008' , 'Zhang-2006' , 'Dongari-2011' , 'Arlemark-2010' , 'Guo_Shu-2013'
	bool boolean_wfunction = false;	
	std::string wfunction ("Guo-2008");
	double* PHI = Wall_function(wfunction,Kn,ymax,dx,ny,rank,rho_out,lat);
	
//*************************RAREFIED METHOD*********************************//

	//Choix de la méthode de raréfaction utilisée (dépend de la condition limite, présence ou non de wall function)
	std::string rarefied_method;
	// CBBSR DBB MR Continuous
	std::string rarefied_BC ("Continuous");
	// Guo-2008 Guo-2011 Yudong-2016 Verhaeghe-2009 Li-2011
	std::string rarefied_article ("");
	// Avec ou sans wall function
	rarefied_method = rarefied_article +"_"+rarefied_BC;
	if(!boolean_wfunction){temp2 = rarefied_model(rarefied_method,cs,dt,nu,Kn,ymax,r,beta,dx,A);}
	else{temp2 = rarefied_model_function(rarefied_method,cs,dt,nu,Kn,ymax,r,sigma,beta,dx,ny,rank,rank2,PHI,A);}
	r = temp2[0];
	beta = temp2[0];	
	tau_s = temp2[1];
	tau_q = temp2[2];
			
		
//*******************************AFFICHAGE INITIAL******************************//	

	printf("nu = %f\n",nu);
	printf("Kn = %f\n",Kn);
	printf("r : %f, tau_s : %f, tau_q : %f\n",r,tau_s,tau_q);
	Fi[0] = (rho_in-rho_out)*cs*cs/(xmax);
	//Umax = Fi[0]*sqrt(6)*ymax;
	//printf("Umax  = Fi * sqrt(6)*ymax %.8f\n", Umax);
	//printf("Re %f\n", Umax * ymax/nu);
	//printf("Ma %f\n",Umax/cs);	
	Umax = Fi[0]*ymax*ymax/(8*nu);//Avec body force
	//mu = sqrt(2/(3*M_PI))*Kn*ymax/dx;
	//Umax = ymax*ymax*(rho_in-rho_out)*cs*cs/(8*mu*xmax);//Vitesse max pour gradient de pression
	printf("Umax : %.10f\n",Umax);
	//Re = Umax * ymax*ratio/nu;
	printf("Re %f\n", Re);
	printf("Ma %f\n",Umax/cs);
	printf("Domaine : [ %.0f %.0f ] x [ %.0f %.0f ]\n",xmin,xmax,ymin,ymax);
	printf("Erreur d'arrêt : %f\n", error);
	printf("Nombre de lattices : %d\n",N);  
	printf("Porosité : %f\n",poro);
	//printf("Diamètre hydraulique : %f\n",hydraulic_diameter_ellipse);
	//printf("Aspect ratio : %f\n",aspect_ratio_ellipse);
	/*omp_set_num_threads(4);
	 std:: cout << "  Number of processors available = " << omp_get_num_procs ( ) << "\n"<< std::endl;
	  std ::cout << "  Number of threads =              " << omp_get_max_threads ( ) << "\n" <<std::endl;*/
	  
//*************************************INITIALISATION RELAXATION MRT*****************************************//

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

						
//*************************INITIALISATION DU DOMAINE******************************************//
initialisation_domain (N,nx,ny,lat,rho_in,rho_out,v_in,v_out,Uw, D,Q,xi,sigma, omega_i, f_star, tau, position);
for (j=0;j<N;j++)
{
	if(typeLat[j])
	{
		lat.rho_[j] = 0.5;
	}
}
//writeLattice(domain,"LBM",Kn,poro,name,0,lat);
//******************************BOUCLE TEMPORELLE******************************************//
while((erreur>error || erreur<-error))	
{		
	
	for (j=0 ; j<N ; j++) //Pour chaque lattice
    {
		
		/*mu = lat.rho_[j]* sqrt(2/(3*M_PI))*Kn*rho_out/lat.rho_[j]*ymax/dx;
		tau_s = 0.5 + 3*mu/lat.rho_[j];
		tau_q = (8-1/tau_s)/(8*(2-1/tau_s));
		beta = (3*mu-Kn*ny)/(3*mu+Kn*ny);
		//printf("mu, tau_s, tau_q, Kn, beta : %f %f %f %f %f \n",mu ,tau_s, tau_q, Kn, beta);		
		//Kne = Kn*rho_out/lat.rho_[j];
		MRT_S(Q,tau_s, tau_q,Si);
		matrix_product(invM,Si,C,Q);
		matrix_product(C,M,C4,Q);*/

		/*Kne = Kn*rho_out/lat.rho_[j];
		lambdae = Kne*ymax/dx;
		test[0] =  rank[j/nx][0]/lambdae; //alpha
		test[1] = (ymax-rank[j/nx][0])/lambdae;//alpha1
		test[2] = 1+(test[0]-1)*exp(-test[0])-test[0]*test[0]*Ei_big(1,test[0]);//phi
		test[3] = 1+(test[1]-1)*exp(-test[1])-test[1]*test[1]*Ei_big(1,test[1]);//phi1
		test[4] = 0.5*(test[2] + test[3]); //PHI
		Knee = Kne*test[4];
		tau = 0.5 + sqrt(6/M_PI)*ymax/dx*Knee;
		r = 1/(1+sqrt(M_PI/6)*(dx/ymax*dx/ymax/(4*Knee)+A[0]+(2*A[1]-8/M_PI)*Knee));*/
		//printf("Lattice %d, tau : %f, r : %f\n",j,tau,r);		
		for (k=0;k<Q;k++)
		{
			fi_equilibre (j,k,lat.rho_[j],cs,lat,lat.u_[j],xi,D,Qi,buffer,omega_i,sigma);	
		}
		//regularized_BC_v_inlet(j,k,lat,cs,v_in,nx,xi,D,Qi,buffer,omega_i,cas[j],Pi_neq,Q,f_neq, bb,sigma);
		regularized_BC_p_inlet(j,k,lat,cs,rho_in,xi,D,Qi,buffer,omega_i,cas[j],Pi_neq,Q,f_neq,bb,sigma);
		regularized_BC_p_outlet(j,k,lat,cs,rho_out,xi,D,Qi,buffer,omega_i,cas[j],Pi_neq,Q,f_neq,bb,sigma);
		for (k=0;k<Q;k++)
		{
			fi_equilibre (j,k,lat.rho_[j],cs,lat,lat.u_[j],xi,D,Qi,buffer,omega_i,sigma);	
		}
		/*for (k=0;k<Q;k++)
		{
			f_star[j][k] = lat.f_[j][k]-1/tau*(lat.f_[j][k]-lat.f0_[j][k]); //Collision en SRT
		}*/
		
		MRT_equilibre_v2(j,lat,M,Q,buffer);
		for (k=0;k<Q;k++)
		{
			MRT_moment(j,k,lat,M,Q);
		}
		
		//MRT_forcing(j,Q,cs,omega_i,xi,Fi,F_bar,lat,temp);
		//MRT_collision(j,f_star,C,lat,Q,temp);
		MRT_collision_v2(j,f_star,C4,lat,Q,temp);
		//MRT_forcing_collision(j,f_star,C,lat,dt,F,C3,F_bar,Q,temp);
		
		//MRT_forcing(j,k,cs,omega_i,xi,Fi,F_bar,lat,temp);
		//MRT_forcing_collision(j,f_star,C_wall[(int)j/nx],lat,dt,F,C3_wall[(int)j/nx],F_bar,Q,temp);	
		//MRT_collision(j,f_star,C_wall[(int)j/nx],lat,Q,temp);
		

		//periodic_WE_BC(j,nx,ny,cas[j],lat,f_star); 	
		//DBB_N_BC(j,cas[j],lat,beta,f_star);		
		//DBB_N_BC_Couette(j,cas[j],lat,beta,f_star,cs,Uw,buffer,omega_i,xi,D,Q,Qi,sigma);
		//DBB_S_BC(j,cas[j],lat,beta,f_star);
		//CBBSR_N_BC_Couette(j,cas[j],lat,r,f_star,uw,xi,cs, omega_i);
		//CBBSR_N_BC(j,cas[j],lat,r,f_star);
		//CBBSR_S_BC(j,cas[j],lat,r,f_star);  
		//bounceback_solid_BC(nx,j,lat,f_star,conn,typeLat,bb,nombre,pos,cas[j]);
		linear_interpolation_method(j,Q,lat,f_star,conn,typeLat,bb,solid_fraction_interpolation, tab_marquage, cas[j]);
		//quadratic_interpolation_method(j,Q,lat,f_star,conn,typeLat,bb,solid_fraction_interpolation, tab_marquage, cas[j]);
		multireflection_interpolation_method(j,Q,lat,f_star,conn,typeLat,bb,solid_fraction_interpolation,tab_marquage,cas[j],C,t,teq,Fpc,mu,lat.rho_[j],Si);
		//pression_in_BC( j,cas[j],lat,xi_r,rho_in);
        //pression_out_BC( j,cas[j],lat,xi_r,rho_out);		
		periodic_NS_BC(j,nx,ny,cas[j],lat,f_star); 
		//periodic_pressure_WE_BC (j, nx, ny, cas[j], lat, f_star, dRHO, sigma,cs);		
	   // bounceback_N_BC(j,cas[j],lat,f_star);
       // bounceback_S_BC(j,cas[j],lat,f_star);
	   
		propagation(j,lat,f_star, typeLat, conn, bb,Q);
		if(!typeLat[j])
		{
			density(j,Q,lat,sigma); 
			velocity(j,D,Q,xi,lat,sigma);
			//lat.u_[j][0]+=0.5*dt*Fi[0]/(lat.rho_[j]);
		}
		else
		{
			lat.rho_[j] = 1;
			lat.u_[j][0] = 0;
			lat.u_[j][1] = 0;
		}
	}
	
               // Ecriture des résultats	
        if (it%outputFrequency==0 && it!=0) 
		{
			char name = FileName(Kn);
			vorticite(nx,ny,lat,dx,typeLat,conn);
			ttsity = tortuosite(nx,ny,lat,dx,typeLat,conn);
			//writeLattice(domain,"LBM",Kn,poro,name,it,lat);
			valeur2 = 0;
			for (int j=0;j<N;j++)
			{
				valeur2+=sqrt(lat.u_[j][0]*lat.u_[j][0]+lat.u_[j][1]*lat.u_[j][1]);
			}
			erreur = 1/sqrt(N)*(valeur2-valeur1)/valeur2;
			valeur1 = valeur2;
			printf("Itération : %d\t Erreur  %.10f, tortuosité : %f\n", it,erreur,ttsity);
		}
		it++;
}
	//char name = FileName(Kn);
	writeScalar(domain, "marquage", it, tab_marquage);
	
	writeLattice(domain,"LBM",Kn,poro,name,it,lat);
	time_t timer2 = time(NULL);
	printf("Temps d'exécution : %d s\n",(int)(timer2-timer1));
	printf("Temps d'exécution pour 1000 itérations: %.2f s\n",(double)((int)(timer2-timer1))/it*1000);

    return 0; 
}


