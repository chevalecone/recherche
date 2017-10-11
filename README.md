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

Guo-2008_MRT : wall function codée, profils de Poiseuille tracés et comparés avec Ohwada, mais pb avec homogénéisation par u_r : aH*sqrt(2/RT) il faut diviser par 2...

A FAIRE : comparaison avec N-S 1°slip order, 2°slip order et A1/A2

Regularized BC à regarder
Changement pour la fonction propagation : on effectue tous les propagations sur la lattice j (comme ça une seule boucle for)

Changer la méthode pour equilibrium ? calculer sur (nx+1)*(ny+1)

08/08/17
Abandon du calcul sur (nx+1)*(ny+1) : la modification de la méthode de calcul est trop conséquente

Utilisation de l'equilibrium inlet/outlet sur un maillage nx*ny avec une extrapolation UNIQUEMENT pour la vitesse de la lattice voisine.
Autres extrapolations possibles : u(j) = 2*u(j+1)-u(j+2) (pour inlet), u(j) = 2*u(j-1)-u(j-2) (pour outlet)

PS : pour des Ma plus élevés, on retrouve les mêmes résultats que Z-H, à savoir qu'il y a un gonflement de la vitesse (une sorte de vorticité à l'inlet puis retour à un Poiseuille classique) : le milieu du conduit ressemble à un Poiseuille "classique" ssi le conduit est assez allongé.
Retour aux conditions périodiques, calcul d'un Poiseuille classique, puis MEP de la CBBSR (retrouver la référence)
Lire les autres références qui peuvent être exploitables

09/08/17
ZH ou extrapolation inlet/outlet (Dirichlet BC) marche lorsque le conduit est assez long ie. le milieu de l'écoulement n'est pas influencé par les BC et les singularités.

Calcul en écoulement raréfié avec la MRT_Guo_2008 sans vitesse du mur (ie pas de Couette, on fait Poiseuille) et sans Wall function pour l'instant.
Kn = 0.001, 0.1, 0.2, 0.3, 0.4, 0.5

Création du MATLAB pour post- pour écoulement raréfié avec Guo_2008

10/08/17
Code DBB fait

Code MR (specular reflection + diffuse bounceback) réalisé

Problème avec CBBSR et autres BC avec du glissement : le traitement en utilisant f* donne des vitesses de glissement avec un coefficient de réflexion pure (r = 1) --> considérer de travailler avec lat.f_
Wall function pour Guo_2008 à réaliser (via Numerical Recipes) car sans wall function, les valeurs ne sont pas pertinentes
Problème dans le calcul de la wall function, l'exponentiel intégrale diverge : regarder comment est défini le MFP

11/08/17
Relire bien les inlet/outlet BC (equilibrium et extrapolation notamment)

DBB et CBBSR : trouver les références + écrire le cheminement (théorie, implémentation en fonction des grandeurs d'intérêt etc...)

14/08/17
Retour à Guo-2008_MRT : implantation de l'intégrale exponentielle de "Numerical Recipes" fructueux avec la première expression, les valeurs de PHI sont entre 0 et 1.

Revue rapide des différentes solutions analytique des vitesses de glissement, et des profils de Poiseuille (premier et second ordre)

Création d'une structure donnant grâce à l'ordonnée de la lattice, le rang dans lequel elle est et donc la valeur du wall function associée, et donc le temps de relaxation, les matrices de relaxation etc...

15/08/17
Guo_MRT fini : étude de profil de vitesse axiale, débit massique réalisé

10/11/2017
Etude des écoulements raréfiés terminée : Guo_DBB ou Guo_CBBSR équivalent, les coefficients proposés respectivement par Guo_2008 et Guo_2011 donne des résultats proches des résultats de comparaison d'Ohwada (1989, Boltzmann linéarisé).
La wall function de Stops permet d'obtenir un profil de vitesse (et donc un débit) proche des résultats d'Ohwada.
Néanmoins, les paramètres (coefficient de CBBSR/DBB et temps de relaxation tau_s et tau_q) dépendent des coefficients des vitesses de glissement provenant de la résolution de l'équation de Boltzmann (premier et second ordre), et ces coefficients de la vitesse de glissement ont été déterminés avec une fitting curve. Ceci donne des résultats différents (pour le profil de vitesse) en comparant avec la DSMC (par exemple) 
Condition limite régularisée implémentée en MRT pour les écoulements continus (comparaison avec un simple écoulement de Poiseuille)
Ecriture d'un artcicle à considérer par rapport à la comparaison entre DBB et CBBSR en utilisant wall function, tau_s et tau_q dépendant de Kn et en MRT --> Revue de littérature, voir ce qui existe ou non, dans quels journaux publier, puis faire un squelette d'article.
Perméabilité : article de référence : Jeong-2006 pour les résultats simplement en régime continu, calcul de la perméabilité en fonction des différentes structures (square/cylinder et in-line/staggered) et en fonction du Re, de Kn et de la porosité.
Utilisation des résultats de Mostafavi-2016 également
