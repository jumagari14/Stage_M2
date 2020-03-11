function [sol,par]=traitement(i,ks,kd,DPA_R,DPA,nom,R,P)
% traitement d'un couple (RNA,PROT)=(r,p)
% i=numéro, DPA=temps de mesures,mesures RNA et Proteines
%
global ptrans3 ptrans6 St6 mt6 St3 mt3 degre
global xi yi h fitR

% calcul des paramètres des polynomes de regression de R et visualisation
h=8;fitR=0;
choix=fit_testRNA( i,nom,DPA_R,R);
switch choix
    case '1'
fitR=3;%poly degré3
    case '2'
fitR=6;%poly degré 6
    case '3'
fitR=0;%estimateur à noyau
    case '4'
fitR=1;
end;
rpoly=mRNA_deg(DPA);


% choix du fit ????
% en fonction du type de fit choisi, calcul des valeurs théoriques R(t)
%[ rpoly ] = fitRNA( i,nom,DPA_R,R,fitR );


%condition initiale: moyenne des proteines aux premiers DPA
P0=initiale(DPA,P);
lb=[P0*0.1;ks*0.8;0];
ub=[P0*2;ks*1.2;inf];
% intialisation paramètres
%kd=0.01;
affine=1;
parI=[P0,ks,kd];
[parMu,resnorm1,residual1,exitflag1,output1]  = lsqcurvefit(@resol_mu,parI,DPA,P,lb,ub);
parMu
go=1;
if go==1
% fit polynomiale de R 
    % calcul explicite de la solution possible
% calcul solution avec paramètres initiaux
 [solI]=solexplicite(parI,DPA);

% recherche des paramètres (ks,kd) minimisant les écarts
% avec solution explicite ED
[parK2,resnorm2,residual2,exitflag2,output2]  = lsqcurvefit(@solexplicite,parI,DPA,P,lb,ub);
if affine==1
% avec solution numérique ED
[parK,resnorm1,residual1,exitflag1,output1]  = lsqcurvefit(@resol_all,parK2,DPA,P,lb,ub);

% avec solution numérique ED + taux de croissance
[parMu,resnorm1,residual1,exitflag1,output1]  = lsqcurvefit(@resol_mu,parK2,DPA,P,lb,ub);
parMu
else
    parK=parK2.*0;
    
end;
% calcul exact de la solution de l'ED
 [sol2]=solexplicite( parK2,DPA);
  sol1=resol_all(parK,DPA);
  solMu=resol_mu(parMu,DPA);

 % sans Bornes paramètres
lb=[P0*0.1;0;0];
ub=[P0*2;inf;inf];% avec solution explicite ED
[parKL,resnormL,residualL,exitflagL,outputL]  = lsqcurvefit(@solexplicite,parK2,DPA,P,lb,ub);
% calcul de la solution explicite de l'ED
 [solL]=solexplicite( parKL,DPA);
par=[parI;parK;parK2;parKL;parMu];
sol=[solI,sol1,sol2,solL,solMu];
end%go
end

