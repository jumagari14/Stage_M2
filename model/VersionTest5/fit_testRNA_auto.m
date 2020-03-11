function [err]=fit_testRNA( i,nom,DPA_R,R)
% diff�rents fits de R: polynomial, estimation par noyau
global ptrans3 ptrans6 St6 mt6 St3 mt3 plogtrans3 Slogt3 mlogt3 degre h fitR

% Lissage des donn�es RNA: regression polynomiale d'ordre 3 de (DPA_i, R_i)
% param�tres polynomes (DPA de R centr�, r�duit)
[ptrans3,St3,mt3] = polyfit(DPA_R,R,3);
[ptrans6,St6,mt6] = polyfit(DPA_R,R,6);
% transformation log des donn�es=>fit par exp(poly deg3)
[plogtrans3,Slogt3,mlogt3] = polyfit(DPA_R,log(R),3);
normR=norm(R,2);h=8;%valeur par d�faut bande estimateur noyau
sauvefitR=fitR;

% valeurs du poly deg3 � DPA (temps experimentaux de mesures de P)
fitR=3;%poly degr�3
rpoly=mRNA_deg(DPA_R);
fitR=6;%poly degr� 6
rpoly6=mRNA_deg(DPA_R);
fitR=0;%estimateur � noyau h=6
rfit=mRNA_deg(DPA_R);
fitR=1;%estimateur noyau h=8
rfit2=mRNA_deg(DPA_R);% calcul des param�tres des polynomes de regression de R et visualisation
fitR=2;%exp(pol deg3)
rexppoly=mRNA_deg(DPA_R);% calcul des param�tres des polynomes de regression de R et visualisation
% calcul erreur
erreur=[norm(R-rpoly,2),norm(R-rpoly6,2),norm(R-rfit,2),norm(R-rfit2,2),norm(R-rexppoly,2)]./normR;



% choix du fit des RNA
fitR=sauvefitR;
% on renvoit l'erreur correspondante
switch fitR
        case 3
            err=erreur(1);%poly degr�3
        case 6
            err=erreur(2);%poly degr� 6
        case 0
            err=erreur(3);%estimateur � noyau h=8
        case 1
            err=erreur(4);%estimateur � noyau h=6
        otherwise % exp(poly)
            err=erreur(5);
    end;


