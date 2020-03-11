function [kparg,sol,errg,score]=solg(NOM,DPA,P,ks_nmin,score);
% détermination de (ks,kd) à partir des données en gFW
% normalisées: (R,P) avec prise en compte du taux de croissance
% % dp/dt=ks*r-(kd+mu)*p (données en gFW)

% starting value condition initiale EDO: moyenne des valeurs initiales mesurees
[P0,Pmin,Pmax]=initiale(DPA,P);

% solution EDO en fonction des (ks,kd) détermines par les MC
sol=zeros(length(DPA),4);

% 3 starting values de ks 1)milieu=3*ksmin 2)ksmin 3)5ksmin
parI=[P0;ks_nmin*3;ks_nmin*0.3];
parInit=[P0;ks_nmin;ks_nmin*0.1];
parInit2=[P0;ks_nmin*5;ks_nmin*5];
% avec solution numérique ED + taux de croissance sans bornes (sauf
% positivité)
lb=[0;0;0];ub=[inf;inf;inf];
[parMu00,resnorm00,residual0,exitflag0,output0]  = lsqcurvefit(@resol_mu,parI,DPA,P,lb,ub);
[parMu01,resnorm01,residual01,exitflag01,output01]  = lsqcurvefit(@resol_mu,parInit,DPA,P,lb,ub);
[parMu02,resnorm02,residual02,exitflag02,output02]  = lsqcurvefit(@resol_mu,parInit2,DPA,P,lb,ub);
parMu00,parMu01,parMu02
resnorm00,resnorm01,resnorm02
%resnorm00,resnorm01,resnorm02

% Evaluation de l'optimisation (sans bornes)et choix du meilleur parametre
[ score,par,resnorm0 ] = EvaluationOptim( score,parMu00,parMu01,parMu02,resnorm00,resnorm01,resnorm02 );
parMu0=par(:,1);
% recherche avec des bornes sur les paramètres (ks,kd) avec ks dans [ksmin;5*ksmin]
% bounds for parameters search
lb=[Pmin;ks_nmin;0];
ub=[Pmax;ks_nmin*5;inf];
% avec solution numérique ED + taux de croissance 
[parMu1,resnorm1,residual1,exitflag1,output1]  = lsqcurvefit(@resol_mu,parI,DPA,P,lb,ub);
% avec solution numérique ED + taux de croissance
% autres starting points
[parMu2,resnorm2,residual2,exitflag2,output2]  = lsqcurvefit(@resol_mu,parMu0,DPA,P,lb,ub);
parInit=parI+randn(size(parI))*0.1.*parI;
[parMu3,resnorm3,residual3,exitflag3,output3]  = lsqcurvefit(@resol_mu,parInit,DPA,P,lb,ub);


% Evaluation de l'optimisation et choix du meilleur parametre (classement)
[ score,par,resnorm] = EvaluationOptim( score,parMu1,parMu2,parMu3,resnorm1,resnorm2,resnorm3);
% evaluation erreur fit proteine
[ score ] = ScoreErreur( sqrt(resnorm(1))/norm(P,2),score );
% paramètres déterminés en résolvant l'ED
% numeriquement avec bornes sur (ks,kd), num avec bornes, num avec bornes, ss bornes
kparg=[par(:,1:end) parMu0];
errg=[sqrt(resnorm) sqrt(resnorm0(1))]./norm(P,2);
parMu1=par(:,1);parMu2=par(:,2);parMu3=par(:,3);

% solutions t->P(t)
sol(:,1)=resol_mu(parMu1,DPA);
sol(:,2)=resol_mu(parMu2,DPA);
sol(:,3)=resol_mu(parMu3,DPA);
% solution sans borne sur (ks,kd)
sol(:,4)=resol_mu(parMu0,DPA);

end

