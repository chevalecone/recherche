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

Essai du calcul de DBB
