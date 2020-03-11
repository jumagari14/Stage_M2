function [X,Y,Z] = matrice_sensibilite(ti,param)
% calcul de la matrice de sensibilité 
% pour 3 paramétrages différents de l'ED
% X=(dp/dparam) en chaque ti pour chaque paramètre
% param=[p0,ks,kd] pour dP/dt=ks*r-(kd+mu)*p
% Y=(dp/dparam) en chaque ti pour chaque paramètre
% param=[p0,K=ks/kd,kd] pour dP/dt=kd*(K*r-p)-mu*p
% Z=(dp/dparam) en chaque ti pour chaque paramètre
% param=[p0,K=kd/ks,ks] pour dP/dt=ks*(r-K*p)-mu*p

% conditions initiales
s0=param(1);
CI=[s0;1;0;0];
[ti,ia,ic]=unique(ti);
[t,sol]= ode45(@(t,z)derive(t,z,param),ti, CI); %parametres après fit
X=sol(:,2:4);
[t,sol1]= ode45(@(t,z)derive_reparam1(t,z,param),ti, CI); %parametres après fit
Y=sol1(:,2:4);
[t,sol2]= ode45(@(t,z)derive_reparam2(t,z,param),ti, CI); %parametres après fit
Z=sol2(:,2:4);