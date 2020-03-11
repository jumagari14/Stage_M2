function [ rfit ] = fitRNA( i,nom,DPA,R,fit )
% Lissage des données RNA par un polynome ou par estimateur à noyau
global ptrans3 ptrans6 St6 mt6 St3 mt3 degre h

% fit=0,3,6
if fit==3
    % : regression polynomiale d'ordre 3 de (DPA_i, R_i)
    degre=3;
    % valeurs du poly à DPA (temps experimentaux)
    [rfit,deltat3] = polyval(ptrans3,DPA,St3,mt3);
elseif fit==6
    degre=6;
    % valeurs du poly à DPA (temps experimentaux)
    [rfit,deltat6] = polyval(ptrans6,DPA,St6,mt6);
    % affichage regression
elseif fit==0; %estimateur à noyau
    h=8;%estimateur à noyau largeur de bande
    rfit=ksr2(DPA,R,h);
else %estimateur à noyau
    h=6;%estimateur à noyau largeur de bande
    rfit=ksr2(DPA,R,h);
end

col='m';
[hAx1,hAx5]=affreg(i,nom,DPA,R,rfit,col)
title('Regression RNA ')
    xlabel('DPA ')
    ylabel('RNA')
