function [erreur] = Lab5_normeL2(mg,NPOINTS,X,T) 
% Fonction appliquant la norme d'erreur L2([A,B]) 
% =============================================== 
[w,t] = Int_Gauss(mg);
normeL2 = 0; 
h=X(2)-X(1); 
for i = 1:NPOINTS-1, 
somme = 0; 
% Intégration de gauss 
% ==================== 
for j = 1:mg, 
xh =h/2 *t(j)+ (X(i)+X(i+1))/2; 
L0 = (X(i+1)-xh) / h; 
L1 = (xh-X(i) ) / h; 
Th = L0 * T(i) + L1 * T(i+1); 
Tex = ; 
% INSCRIRE ICI LA SOLUTION EXACTE 
somme = somme + (Tex-Th)^2 * w(j) * h/2; 
end 
normeL2 = normeL2 + somme; 
end 
% Sortie de la valeur de normeL2 
% ============================== 
normeL2 = sqrt(normeL2); 
erreur = normeL2;