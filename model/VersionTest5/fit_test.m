function fit_test( i,nom,DPA,R)
% différents fits de R: polynomial, estimation par noyau
global ptrans3 ptrans6 St6 mt6 St3 mt3 degre h fitR

% Lissage des données RNA: regression polynomiale d'ordre 3 de (DPA_i, R_i)

% paramètres polynomes (DPA centré, réduit)
[ptrans3,St3,mt3] = polyfit(DPA,R,3);
% valeurs du poly à DPA (temps experimentaux)
[rpoly,deltat3] = polyval(ptrans3,DPA,St3,mt3);
col='c';
[hAx1,hAx2]=affreg(i,nom,DPA,R,rpoly,col)

[ptrans6,St6,mt6] = polyfit(DPA,R,6);
% valeurs du poly à DPA (temps experimentaux)
[rpoly6,deltat6] = polyval(ptrans6,DPA,St6,mt6);
% affichage regression
col='b';
[hAx1,hAx3]=affreg(i,nom,DPA,R,rpoly6,col)

rfit=ksr2(DPA,R,h);
col='r';
[hAx1,hAx4]=affreg(i,nom,DPA,R,rfit,col)

rfit2=ksr2(DPA,R,h-2);
col='m';
[hAx1,hAx5]=affreg(i,nom,DPA,R,rfit2,col)
erreur=[norm(R-rpoly,2),norm(R-rpoly6,2),norm(R-rfit,2),norm(R-rfit2,2)]./norm(R,2);
% fit par défaut: estimateur à noyau
%fitR=3;
%if choix==0
%      choix = input('fitR 0 3 6 ', 's');
%end
limaxis=axis;
xPos=limaxis(1);
yl = limaxis(4);
yPos = yl*.7;
texte=text(xPos,yPos,sprintf('er rel_{pol3,pol6,noyau1,noyau2}=%8.4f %8.4f %8.4f %8.4f',erreur),'FontSize',8);
legend([hAx1,hAx2,hAx3,hAx4,hAx5],{'points exp.','poly deg3','poly deg6', strcat('noyau h=',num2str(h)),strcat('noyau h=',num2str(h-2))},'FontSize',5,'Location','northwest');

end

