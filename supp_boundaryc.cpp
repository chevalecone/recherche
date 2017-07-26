void curved_bounceback_collision_BC( int j,  int cas, Lattice lat, double** f_star, double tau, double** xi, double**cylinder,double* omega_i,int** conn, double cs,  int D)
{
	double* uw = new double[2];
	uw[0] = 0;
	uw[1] = 0;

	double* ubf = new double[2];
	double* uff = new double[2];
	double* uf = new double[2];
	double* uInt = new double[2];
	double* mem = new double[4];
	//mem[0] = feq, mem[1] = chi, mem[2] = rho, mem[3] = delta

	switch(cas)
	{
		case 26: //Calcul de la population post-collision sur la lattice solide
		{
			if((int)cylinder[j][0]!=0)
			{
				for (int k=1;k<9;k++)
				{
					mem[3] = cylinder[j][k];
					if(mem[3]>pow(10,-5))
					{
						if (mem[3]<0.5)
						{
							mem[1] = 1/(tau-2)*(2*mem[3]-1);
							uff[0] = lat.u_[conn[conn[j][k]][k]][0];
							uff[1] = lat.u_[conn[conn[j][k]][k]][1];
							ubf[0] = uff[0];
							ubf[1] = uff[1];
						}
						else
						{
							mem[1] = 1/(tau+0.5)*(2*mem[3]-1);
							ubf[0] = 0.5/(mem[3])*(2*mem[3]-3)*lat.u_[conn[j][k]][0];
							ubf[1] = 0.5/(mem[3])*(2*mem[3]-3)*lat.u_[conn[j][k]][1];
						}
						mem[0] = lat.f0_[conn[j][k]][k];
						mem[2] = lat.rho_[conn[j][k]];
						uf[0] = lat.u_[conn[j][k]][0];
						uf[1] = lat.u_[conn[j][k]][1];
						uInt[0] = ubf[0]-uf[0]-2*uw[0];
						uInt[1] = ubf[1]-uf[1]-2*uw[1];
						f_star[j][k] = f_star[conn[j][k]][k]-mem[1]*(f_star[conn[j][k]][k]-mem[0])+omega_i[k]*mem[2]/(cs*cs)*pscal(xi[k],uInt,D);
					}
				}
			}
		}
		break;
		default:
		{
			for (int k=0;k<9;k++)
			{
				f_star[j][k]=lat.f_[j][k]-1/tau*(lat.f_[j][k]-lat.f0_[j][k]);
			}
		}
		
	}

	//Free memory
	delete[] uw;
	delete[] ubf;
	delete[] uff;
	delete[] uf;
	delete[] uInt;
	delete[] mem;
}

void non_slip_Inamuro_BC( int j,int nx, int ny,  int cas, Lattice lat, double xi_r,double u_w,double v_w)
{
	double rho_p;
	double u_p;
	double rho_w;
	switch(cas)
	{
		case 4: //SUD
		rho_w = 1/(1-v_w)*(lat.f_[j][0]+lat.f_[j][1]+lat.f_[j][3]+2*(lat.f_[j][4]+lat.f_[j][7]+lat.f_[j][8]));
		rho_p = 6/(1+3*v_w+3*v_w*v_w)*(rho_w*v_w+lat.f_[j][4]+lat.f_[j][7]+lat.f_[j][8]);
		u_p   = 6*xi_r*xi_r/(rho_p*(xi_r+3*v_w))*(rho_w*u_w/xi_r-lat.f_[j][1]+lat.f_[j][3]+lat.f_[j][7]-lat.f_[j][8])-u_w; 
		lat.f_[j][2] = 1/9*rho_p  * (1 + 3 * v_w            + 9/2 * v_w*v_w                       - 3/2 * ((u_w+u_p)*(u_w+u_p)+v_w*v_w));
		lat.f_[j][5] = 1/36*rho_p * (1 + 3 * (u_w+u_p+v_w)  + 9/2 * (u_w+u_p+v_w)*(u_w+u_p+v_w)   - 3/2 * ((u_w+u_p)*(u_w+u_p)+v_w*v_w));
		lat.f_[j][6] = 1/36*rho_p * (1 + 3 * (-u_w-u_p+v_w) + 9/2 * (-u_w-u_p+v_w)*(-u_w-u_p+v_w) - 3/2 * ((u_w+u_p)*(u_w+u_p)+v_w*v_w));
		break;
		case 5: //SUD-OUEST
		rho_w = 1/(1-v_w)*(lat.f_[j][0]+lat.f_[j][1]+lat.f_[j][3]+2*(lat.f_[j][4]+lat.f_[j][7]+lat.f_[j][8]));
		rho_p = 6/(1+3*v_w+3*v_w*v_w)*(rho_w*v_w+lat.f_[j][4]+lat.f_[j][7]+lat.f_[j][8]);
		u_p   = 6*xi_r*xi_r/(rho_p*(xi_r+3*v_w))*(rho_w*u_w/xi_r-lat.f_[j][1]+lat.f_[j][3]+lat.f_[j][7]-lat.f_[j][8])-u_w; 
		lat.f_[j][2] = 1/9*rho_p  * (1 + 3 * v_w            + 9/2 * v_w*v_w                       - 3/2 * ((u_w+u_p)*(u_w+u_p)+v_w*v_w));
		lat.f_[j][5] = 1/36*rho_p * (1 + 3 * (u_w+u_p+v_w)  + 9/2 * (u_w+u_p+v_w)*(u_w+u_p+v_w)   - 3/2 * ((u_w+u_p)*(u_w+u_p)+v_w*v_w));
		lat.f_[j][6] = 1/36*rho_p * (1 + 3 * (-u_w-u_p+v_w) + 9/2 * (-u_w-u_p+v_w)*(-u_w-u_p+v_w) - 3/2 * ((u_w+u_p)*(u_w+u_p)+v_w*v_w));
		break;
		case 6: //SUD-EST
		rho_w = 1/(1-v_w)*(lat.f_[j][0]+lat.f_[j][1]+lat.f_[j][3]+2*(lat.f_[j][4]+lat.f_[j][7]+lat.f_[j][8]));
		rho_p = 6/(1+3*v_w+3*v_w*v_w)*(rho_w*v_w+lat.f_[j][4]+lat.f_[j][7]+lat.f_[j][8]);
		u_p   = 6*xi_r*xi_r/(rho_p*(xi_r+3*v_w))*(rho_w*u_w/xi_r-lat.f_[j][1]+lat.f_[j][3]+lat.f_[j][7]-lat.f_[j][8])-u_w; 
		lat.f_[j][2] = 1/9*rho_p  * (1 + 3 * v_w            + 9/2 * v_w*v_w                       - 3/2 * ((u_w+u_p)*(u_w+u_p)+v_w*v_w));
		lat.f_[j][5] = 1/36*rho_p * (1 + 3 * (u_w+u_p+v_w)  + 9/2 * (u_w+u_p+v_w)*(u_w+u_p+v_w)   - 3/2 * ((u_w+u_p)*(u_w+u_p)+v_w*v_w));
		lat.f_[j][6] = 1/36*rho_p * (1 + 3 * (-u_w-u_p+v_w) + 9/2 * (-u_w-u_p+v_w)*(-u_w-u_p+v_w) - 3/2 * ((u_w+u_p)*(u_w+u_p)+v_w*v_w));
		break;
		case 3: //NORD
		rho_w = 1/(1+v_w)*(lat.f_[j][0]+lat.f_[j][1]+lat.f_[j][3]+2*(lat.f_[j][2]+lat.f_[j][5]+lat.f_[j][6]));
		rho_p = 6/(1-3*v_w+3*v_w*v_w)*(lat.f_[j][2]+lat.f_[j][5]+lat.f_[j][6]-rho_w*v_w);
		u_p = 6/(rho_p*(1-3*v_w))*(rho_w*u_w-lat.f_[j][1]+lat.f_[j][3]-lat.f_[j][5]+lat.f_[j][6]);
		lat.f_[j][4] = 1/9*rho_p  * (1 - 3 * v_w            + 9/2 * v_w*v_w                     - 3/2 * ((u_w+u_p)*(u_w+u_p)+v_w*v_w));
		lat.f_[j][7] = 1/36*rho_p * (1 + 3 * (-u_w-u_p-v_w) + 9/2 * (u_w+u_p+v_w)*(u_w+u_p+v_w) - 3/2 * ((u_w+u_p)*(u_w+u_p)+v_w*v_w));
		lat.f_[j][8] = 1/36*rho_p * (1 + 3 * (u_w+u_p-v_w)  + 9/2 * (u_w+u_p-v_w)*(u_w+u_p-v_w) - 3/2 * ((u_w+u_p)*(u_w+u_p)+v_w*v_w));
		break;
		case 7: //NORD-OUEST
		rho_w = 1/(1+v_w)*(lat.f_[j][0]+lat.f_[j][1]+lat.f_[j][3]+2*(lat.f_[j][2]+lat.f_[j][5]+lat.f_[j][6]));
		rho_p = 6/(1-3*v_w+3*v_w*v_w)*(lat.f_[j][2]+lat.f_[j][5]+lat.f_[j][6]-rho_w*v_w);
		u_p = 6/(rho_p*(1-3*v_w))*(rho_w*u_w-lat.f_[j][1]+lat.f_[j][3]-lat.f_[j][5]+lat.f_[j][6]);
		lat.f_[j][4] = 1/9*rho_p  * (1 - 3 * v_w            + 9/2 * v_w*v_w                     - 3/2 * ((u_w+u_p)*(u_w+u_p)+v_w*v_w));
		lat.f_[j][7] = 1/36*rho_p * (1 + 3 * (-u_w-u_p-v_w) + 9/2 * (u_w+u_p+v_w)*(u_w+u_p+v_w) - 3/2 * ((u_w+u_p)*(u_w+u_p)+v_w*v_w));
		lat.f_[j][8] = 1/36*rho_p * (1 + 3 * (u_w+u_p-v_w)  + 9/2 * (u_w+u_p-v_w)*(u_w+u_p-v_w) - 3/2 * ((u_w+u_p)*(u_w+u_p)+v_w*v_w));
		break;
		case 8: //NORD-EST
		rho_w = 1/(1+v_w)*(lat.f_[j][0]+lat.f_[j][1]+lat.f_[j][3]+2*(lat.f_[j][2]+lat.f_[j][5]+lat.f_[j][6]));
		rho_p = 6/(1-3*v_w+3*v_w*v_w)*(lat.f_[j][2]+lat.f_[j][5]+lat.f_[j][6]-rho_w*v_w);
		u_p = 6/(rho_p*(1-3*v_w))*(rho_w*u_w-lat.f_[j][1]+lat.f_[j][3]-lat.f_[j][5]+lat.f_[j][6]);
		lat.f_[j][4] = 1/9*rho_p  * (1 - 3 * v_w            + 9/2 * v_w*v_w                     - 3/2 * ((u_w+u_p)*(u_w+u_p)+v_w*v_w));
		lat.f_[j][7] = 1/36*rho_p * (1 + 3 * (-u_w-u_p-v_w) + 9/2 * (u_w+u_p+v_w)*(u_w+u_p+v_w) - 3/2 * ((u_w+u_p)*(u_w+u_p)+v_w*v_w));
		lat.f_[j][8] = 1/36*rho_p * (1 + 3 * (u_w+u_p-v_w)  + 9/2 * (u_w+u_p-v_w)*(u_w+u_p-v_w) - 3/2 * ((u_w+u_p)*(u_w+u_p)+v_w*v_w));
		break;

	}

/*		rho_w = 1/(1+v_w/xi_r)*(lat.f_[j][0]+lat.f_[j][1]+lat.f_[j][3]+2*(lat.f_[j][2]+lat.f_[j][5]+lat.f_[j][6]));
		rho_p = 6*xi_r*xi_r/(xi_r*xi_r-3*xi_r*v_w+3*v_w*v_w)*(lat.f_[j][2]+lat.f_[j][5]+lat.f_[j][6]-rho_w/xi_r*v_w);
		u_p = 6*xi_r*xi_r/(rho_p*(xi_r-3*v_w))*(rho_w*u_w/xi_r-lat.f_[j][1]+lat.f_[j][3]+lat.f_[j][5]-lat.f_[j][6]);
		lat.f_[j][4] = 1/9*rho_p  * (1 - 3/xi_r * v_w            + 9/(2*xi_r*xi_r) * v_w*v_w                     - 3/(2*xi_r*xi_r) * ((u_w+u_p)*(u_w+u_p)+v_w*v_w));
		lat.f_[j][7] = 1/36*rho_p * (1 + 3/xi_r * (-u_w-u_p-v_w) + 9/(2*xi_r*xi_r) * (u_w+u_p+v_w)*(u_w+u_p+v_w) - 3/(2*xi_r*xi_r) * ((u_w+u_p)*(u_w+u_p)+v_w*v_w));
		lat.f_[j][8] = 1/36*rho_p * (1 + 3/xi_r * (u_w+u_p-v_w)  + 9/(2*xi_r*xi_r) * (u_w+u_p-v_w)*(u_w+u_p-v_w) - 3/(2*xi_r*xi_r) * ((u_w+u_p)*(u_w+u_p)+v_w*v_w));
*/

}

//Fonction changeant la matrice de cas pour le cas d'un cylindre à section circulaire
//Si la lattice solide n'intervient pas dans les conditions limites, le cas est 9
//Dans le cas où la lattice solide est BC, le cas est 26
//Si pour la lattice i considérée, le voisin k (k allant de 1 à 8 représentant les lattices voisines) a un delta non négligeable
// (ie. intervient dans le BC), alors le cas de cette lattice est 27
void solidCurvedCondition(int nx, int ny,double** cylinder, int* cas)
{
	 int N = nx*ny;
	for ( int i=0;i<N;i++)
	{
		if (cylinder[i][0]!=0) // si le noeud est solide
		{   
			double sum =0;
			for ( int j=1;j<9;j++)//on regarde parmi ses voisins le delta entre lui et ses voisins
			{
				sum+=cylinder[i][j];
			}
			if(sum<pow(10,-5)) //s'il n'y a pas de delta, c'est un noeud solide INTERIEUR
			{
				cas[(int)cylinder[i][0]] = 9;
			}
			else
			{
				cas[(int)cylinder[i][0]] = 26; //Sinon c'est un noeud solide proche des frontières
			}
		}
	}
}

double** curvedCylinder (double abscisse, double ordonnee, double diametre, int nx, int ny, double** position,int** conn)
{

	 int N = nx*ny;
	double norme = 0;
	double R = 0.5*diametre;
	 int pos_max = N-1;
	 int** tabInt = new  int*[N]; //Tableau intermédiaire avec position de la lattice solide
					  			//et booléens pour déterminer si noeud solide ou fluide
	double** cylinder = new double*[N];
	int pos_min = 0;
	double x0 = abscisse;
	double y0 = ordonnee;
	double xf;
	double yf;
	double xw;
	double yw;
	double xs;
	double ys;
	double L;
	double delta;

	for ( int i=0;i<N;i++)
	{
		tabInt[i] = new  int[9];
		cylinder[i] = new double[9];
		for ( int j=0;j<9;j++)
		{
			tabInt[i][j] = 0;
			cylinder[i][j] = 0;
		}
	}

	//Recherche des positions des lattices dans lequel est inscrit le cylindre
	while(position[pos_min][0]<x0-R || position[pos_min][1]<y0-R)
	{
		pos_min++;
	}
	pos_min--; // Position de la première lattice fluide avant le solide
	while(position[pos_max][0]>x0+R || position[pos_max][1]>y0+R)
	{
		pos_max--;
	}
	pos_max++; //Position de la première lattice fluide après le solide

	printf("Pos min : %d\n",pos_min);
	printf("Pos max : %d\n", pos_max);

	for ( int i=pos_min-1;i<pos_max+1;i++) //Sur le carré restreint, on regarde si la lattice en question est effectivement solide
	{
		norme = (position[i][0]-x0)*(position[i][0]-x0)+(position[i][1]-y0)*(position[i][1]-y0);
		if (norme<=R*R)
		{
			tabInt[i][0] = i; //La lattice i est solide
			for (int j=1;j<9;j++) //Si la lattice i est solide, on regarde parmi ses voisins lesquels sont les noeuds fluides
								  //les noeuds de la frontière solide x_w vont être situés entre les deux
			{
				norme = (position[conn[i][j]][0]-abscisse)*(position[conn[i][j]][0]-abscisse)+(position[conn[i][j]][1]-ordonnee)*(position[conn[i][j]][1]-ordonnee);
				
				if(norme<=R*R)
				{
					tabInt[i][j+1] = 1;
				}
			}
		}
	}

	//Il s'agit maintenant de déterminer la position des xw (coordonnées de la frontière entre xf et xb), et du delta conséquent.
	//On met la valeur de delta dans cylinder

	for ( int i=0;i<N;i++)
	{
		if(tabInt[i][0]!=0)
		{
			cylinder[i][0] = tabInt[i][0]; //Indice de la lattice solide

			//Récupération des coordonnées de la lattice solide
			xs = position[tabInt[i][0]][0];
			ys = position[tabInt[i][0]][1];

			//Voisin 1
			xf = position[tabInt[i][0]+1][0];
			yf = position[tabInt[i][0]+1][1];
			L = x0-xs+sqrt(R*R-(ys-y0)*(ys-y0));
			xw = xs+L;
			yw = ys;
			delta = sqrt(((xf-xw)*(xf-xw)+(yf-yw)*(yf-yw))/((xf-xs)*(xf-xs)+(yf-ys)*(yf-ys)));
			cylinder[i][1] = (1-tabInt[i][1])*delta;

			//Voisin 2
			xf = position[tabInt[i][0]+nx][0];
			yf = position[tabInt[i][0]+nx][1];	
			L = y0-ys+sqrt(R*R-(xs-x0)*(xs-x0));	
			xw = xs;
			yw = ys+L;
			delta = sqrt(((xf-xw)*(xf-xw)+(yf-yw)*(yf-yw))/((xf-xs)*(xf-xs)+(yf-ys)*(yf-ys)));
			cylinder[i][2] = (1-tabInt[i][2])*delta;

			//Voisin 3 
			xf = position[tabInt[i][0]-1][0];
			yf = position[tabInt[i][0]-1][1];
			L = xs-x0+sqrt(R*R-(ys-y0)*(ys-y0));	
			xw = xs-L;
			yw = ys;
			delta = sqrt(((xf-xw)*(xf-xw)+(yf-yw)*(yf-yw))/((xf-xs)*(xf-xs)+(yf-ys)*(yf-ys)));
			cylinder[i][3] = (1-tabInt[i][3])*delta;

			//Voisin 4
			xf = position[tabInt[i][0]-nx][0];
			yf = position[tabInt[i][0]-nx][1];	
			L = ys-y0+sqrt(R*R-(xs-x0)*(xs-x0));	
			xw = xs;
			yw = ys-L;
			delta = sqrt(((xf-xw)*(xf-xw)+(yf-yw)*(yf-yw))/((xf-xs)*(xf-xs)+(yf-ys)*(yf-ys)));
			cylinder[i][4] = (1-tabInt[i][4])*delta;

			//Voisin 5
			xf = position[tabInt[i][0]+nx+1][0];
			yf = position[tabInt[i][0]+nx+1][1];
			L = 0.5*(-(xs-x0+ys-y0)+sqrt(2*R*R-(xs-x0-ys+y0)*(xs-x0-ys+y0)));
			xw = xs+L;
			yw = ys+L;
			delta = sqrt(((xf-xw)*(xf-xw)+(yf-yw)*(yf-yw))/((xf-xs)*(xf-xs)+(yf-ys)*(yf-ys)));
			cylinder[i][5] = (1-tabInt[i][5])*delta;

			//Voisin 6
			xf = position[tabInt[i][0]+nx-1][0];
			yf = position[tabInt[i][0]+nx-1][1];
			L = 0.5*(-(-xs+x0+ys-y0)+sqrt(2*R*R-(xs-x0+ys-y0)*(xs-x0+ys-y0)));
			xw = xs-L;
			yw = ys+L;
			delta = sqrt(((xf-xw)*(xf-xw)+(yf-yw)*(yf-yw))/((xf-xs)*(xf-xs)+(yf-ys)*(yf-ys)));
			cylinder[i][6] = (1-tabInt[i][6])*delta;

			//Voisin 7
			xf = position[tabInt[i][0]-nx-1][0];
			yf = position[tabInt[i][0]-nx-1][1];
			L = 0.5*(xs-x0+ys-y0+sqrt(2*R*R-(xs-x0-ys+y0)*(xs-x0-ys+y0)));
			xw = xs-L;
			yw = ys-L;
			delta = sqrt(((xf-xw)*(xf-xw)+(yf-yw)*(yf-yw))/((xf-xs)*(xf-xs)+(yf-ys)*(yf-ys)));
			cylinder[i][7] = (1-tabInt[i][7])*delta;

			//Voisin 8
			xf = position[tabInt[i][0]-nx+1][0];
			yf = position[tabInt[i][0]-nx+1][1];
			L = 0.5*(-(xs-x0-ys+y0)+sqrt(2*R*R-(xs-x0+ys-y0)*(xs-x0+ys-y0)));
			xw = xs+L;
			yw = ys-L;
			delta = sqrt(((xf-xw)*(xf-xw)+(yf-yw)*(yf-yw))/((xf-xs)*(xf-xs)+(yf-ys)*(yf-ys)));
			cylinder[i][8] = (1-tabInt[i][8])*delta;
		}
			
	}

	//cylinder[i][0] donne l'indice du noeud (s'il est solide), ou 0
	//cylinder[i][1] à cylinder[i][8] donne le delta entre ce noeud et son voisin (de 1 à 8), ou bien 0 s'il n'y en a pas
	return cylinder;
	delete[] cylinder;
	delete[] tabInt;
}
//Fonction permettant de créer un tableau de booléens, où 0 représente un noeud fluide et 1 un noeud solide pour un cylindre circulaire
void typeCurved( int N, double** cylinder, bool* typeLat)
{
	for ( int i=0;i<N;i++)
	{
		if(cylinder[i][0]!=0)
		{
			typeLat[i]=true;
		}
	}
}