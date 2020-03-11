function [err]=fit_testRNA( i,nom,DPA_R,R)
% différents fits de R: polynomial, estimation par noyau
global ptrans3 ptrans6 St6 mt6 St3 mt3 plogtrans3 Slogt3 mlogt3 degre h fitR

% Lissage des données RNA: regression polynomiale d'ordre 3 de (DPA_i, R_i)
% paramètres polynomes (DPA de R centré, réduit)
[ptrans3,St3,mt3] = polyfit(DPA_R,R,3);
[ptrans6,St6,mt6] = polyfit(DPA_R,R,6);
% transformation log des données=>fit par exp(poly deg3)
[plogtrans3,Slogt3,mlogt3] = polyfit(DPA_R,log(R),3);
normR=norm(R,2);h=8;%valeur par défaut bande estimateur noyau
sauvefitR=fitR;

% valeurs du poly deg3 à DPA (temps experimentaux de mesures de P)
fitR=3;%poly degré3
rpoly=mRNA_deg(DPA_R);
fitR=6;%poly degré 6
rpoly6=mRNA_deg(DPA_R);
fitR=0;%estimateur à noyau h=6
rfit=mRNA_deg(DPA_R);
fitR=1;%estimateur noyau h=8
rfit2=mRNA_deg(DPA_R);% calcul des paramètres des polynomes de regression de R et visualisation
fitR=2;%exp(pol deg3)
rexppoly=mRNA_deg(DPA_R);% calcul des paramètres des polynomes de regression de R et visualisation
% calcul erreur
erreur=[norm(R-rpoly,2),norm(R-rpoly6,2),norm(R-rfit,2),norm(R-rfit2,2),norm(R-rexppoly,2)]./normR;

%affichage regression
col='c';[hAx1,hAx2]=affreg(i,nom,DPA_R,R,DPA_R,rpoly,col);%poly deg3
col='b';[hAx1,hAx3]=affreg(i,nom,DPA_R,R,DPA_R,rpoly6,col);%poly deg6
col='r';[hAx1,hAx4]=affreg(i,nom,DPA_R,R,DPA_R,rfit,col);% noyau h=6
col='m';[hAx1,hAx5]=affreg(i,nom,DPA_R,R,DPA_R,rfit2,col);% noyau h=8
col='g';[hAx1,hAx6]=affreg(i,nom,DPA_R,R,DPA_R,rexppoly,col);%exp(pol deg3)

limaxis=axis;
xPos=limaxis(1);
yl = limaxis(4);
yPos = yl*.7;
texte=text(xPos,yPos,sprintf('er rel_{pol3,pol6,noyau1,noyau2,exppol}=%8.4f %8.4f %8.4f %8.4f %8.4f',erreur),'FontSize',8);
legend([hAx1,hAx2,hAx3,hAx4,hAx5,hAx6],{'points exp.','poly deg3','poly deg6', strcat('noyau h=',num2str(h)),strcat('noyau h=',num2str(h-2)),'exp pol'},'FontSize',5,'Location','northwest');
title(['Regression RNA ' nom]);
xlabel('DPA ');
ylabel('RNA');

% courbe pour publi
figure(5001);
col='r';[hAx1,hAx6]=affreg(5001,nom,DPA_R,R,DPA_R,rexppoly,col);%exp(pol deg3)
title(['\fontsize{16}' nom ': RNA profile and fit']);
xlabel('\fontsize{14} Time (Days Post Anthesis)');
ylabel('\fontsize{14} RNA normalized by the mean');

% choix du fit des RNA
fitR=sauvefitR;
% on renvoit l'erreur correspondante
switch fitR
        case 3
            err=erreur(1);%poly degré3
        case 6
            err=erreur(2);%poly degré 6
        case 0
            err=erreur(3);%estimateur à noyau h=8
        case 1
            err=erreur(4);%estimateur à noyau h=6
        otherwise % exp(poly)
            err=erreur(5);
    end;


