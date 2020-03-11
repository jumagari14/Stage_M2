function [kparg,sol,errg,score,message]=solg_ssBorne(NOM,DPA,P,ks_nmin, );
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
parInit3=[P0;ks_nmin*10;ks_nmin*10];
% avec solution numérique ED + taux de croissance sans bornes (sauf
% positivité)
lb=[0;0;0];ub=[inf;inf;inf];
[parMu00,resnorm00,residual0,exitflag0,output0]  = lsqcurvefit(@resol_mu,parI,DPA,P,lb,ub);
[parMu1,resnorm1,residual1,exitflag1,output1]  = lsqcurvefit(@resol_mu,parInit,DPA,P,lb,ub);
[parMu2,resnorm2,residual2,exitflag2,output2]  = lsqcurvefit(@resol_mu,parInit2,DPA,P,lb,ub);
% changement d'algorithm: levenberg-marquardt
options = optimset('Algorithm','levenberg-marquardt');
% attention pas de borne avec cet algo
lb1 = [];
ub1 = [];
[parMu3,resnorm3,residual3,exitflag3,output3]  = lsqcurvefit(@resol_mu,parInit3,DPA,P,lb1,ub1,options);
if find(parMu3<0) % si l'un des paramètres est négatifs, on reprend l'algo par défaut
    % avec la borne lb>0
    disp('negatif avec LevenM')
    parMu3
options = optimset('Algorithm','trust-region-reflective');
[parMu3,resnorm3,residual3,exitflag3,output3]  = lsqcurvefit(@resol_mu,parInit3,DPA,P,lb,ub,options);
end    
parMu00,parMu1,parMu2,parMu3
resnorm00,resnorm1,resnorm2,resnorm3
%resnorm00,resnorm01,resnorm02
resnorm=[resnorm00,resnorm1,resnorm2,resnorm3];
par=[ parMu00,parMu1,parMu2,parMu3];
flag=[exitflag0,exitflag1,exitflag2,exitflag3];
% Evaluation de l'optimisation (sans bornes)et choix du meilleur parametre
[ scoreO,par,resnorm,messageO] = EvaluationOptim( par,resnorm,flag );

parMu0=par(:,1);
parInit=parI+randn(size(parI))*0.1.*parI;
%[parMu3,resnorm3,residual3,exitflag3,output3]  = lsqcurvefit(@resol_mu,parInit,DPA,P,lb,ub);


% evaluation erreur fit proteine
[ scoreF ,messageF] = ScoreErreur( sqrt(resnorm(1))/norm(P,2),score );
% paramètres déterminés en résolvant l'ED
% numeriquement avec bornes sur (ks,kd), num avec bornes, num avec bornes, ss bornes
kparg=[par(:,1:end)];
errg=[sqrt(resnorm)]./norm(P,2);

parMu00=par(:,1);parMu1=par(:,2);parMu2=par(:,3);parMu3=par(:,4);

% % solution sans borne sur (ks,kd): t->P(t)
sol(:,1)=resol_mu(parMu00,DPA);
sol(:,2)=resol_mu(parMu1,DPA);
sol(:,3)=resol_mu(parMu2,DPA);
sol(:,4)=resol_mu(parMu3,DPA);

score=[scoreO;scoreF];
message={messageO;messageF};
end

