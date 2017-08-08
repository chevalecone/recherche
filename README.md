# Recherche
Projet Recherche POLYMTL

01/08/17
MRT - Périodique
Construction de la matrice de passage par Gram-Schmidt

Simulation sur 200 * 50 lattices

Pas de flexibilité : modèle D2Q9 codé en dur pour accélérer le calcul

CBBSR codé pour les frontières N/S avec coefficient d'accommodation arbitraire

CBBSR s'inspire de Succi-2002


02/08/17
Essai de lien coefficient d'accomodation - Kn

Extrapolation method pour inlet/outlet

Extrapolation method abandonnée : modified density ne semble pas intervenir dans la collision en MRT

On teste equilibrium distribution pour inlet/outlet : les densités sont imposées, et les vitesses sont extrapolées.

Problème avec equilibrium distribution : il faut extrapoler les vitesses, mais cela nécessite des noeuds solides qui ne sont pas des noeuds de calcul (bizarre).

Implémentation de Tau = 0.5 + sqrt(6/Pi)*KnH

Implémentation de r = 3*mu/ (KnH*rho_out + 3*mu)

Code MRT : enlèvement f_star dans le bounceback (n'intervient pas)

Calcul de Tau en utilisant d'abord la méthode IPL puis ensuite VHS (Esfahani 2014)

Essai du calcul de DBB (pas encore fait)

03/08/17
Continuer la recherche pour l'equilibrium BC

ZH ne marche vraiment pas

Fonction pour calculer les fi,eq en se basant sur les tenseurs (réussi)

Autres BC à part periodic conditions ne semblent pas marcher
Equilibrium BC : calcul de fi,eq pour les populations manquantes avec une densité imposée, et des vélocités extrapolées <-- comment le sont-elles ?

Parallélisation du code avec OpenMP : librairie omp.h manquante
   - Parallélisation non réussie : compilée, mais ralentit le code ??
  
Reprise des résultats de Arabjamoloei (2016) pour vérifier CBBSR et l'importance d'une viscosité fonction du temps de relaxation.
   
04/08/17
CBBSR : dépend de f_star et pas de lat.f_
   - Note : pas encore déterminé, on obtient pas les mêmes résultats que Arabjamoloei car inlet/outlet est basé sur velocity_inlet et pressure_outlet mep par He & Doolen (2002), alors qu'on utilise des conditions périodiques --> la vitesse max augmente quand on augmente la réflexion spéculaire

Changement pour la fonction propagation : on effectue tous les propagations sur la lattice j (comme ça une seule boucle for)

Changer la méthode pour equilibrium ? calculer sur (nx+1)*(ny+1)

08/08/17
Abandon du calcul sur (nx+1)*(ny+1) : la modification de la méthode de calcul est trop conséquente

Utilisation de l'equilibrium inlet/outlet sur un maillage nx*ny avec une extrapolation UNIQUEMENT pour la vitesse de la lattice voisine.
Autres extrapolations possibles : u(j) = 2*u(j+1)-u(j+2) (pour inlet), u(j) = 2*u(j-1)-u(j-2) (pour outlet)

PS : pour des Ma plus élevés, on retrouve les mêmes résultats que Z-H, à savoir qu'il y a un gonflement de la vitesse (une sorte de vorticité à l'inlet puis retour à un Poiseuille classique)
Retour aux conditions périodiques, calcul d'un Poiseuille classique, puis MEP de la CBBSR (retrouver la référence)
Lire les autres références qui peuvent être exploitables
