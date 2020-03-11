function F = resol_Reparam(param,time)
% new parametrisation
% réécriture de l'ED: reparamétrisation
 % dy/dt=param(2)*(mRNA(t)-K*y)-mu*y
 % param(2)=ks et K=kd/ks=param(3)/param(2)
%paramètres  (condition initiale, kd et ks/kd)
% Utime: returns the same data as in time, but with no repetitions.
%  Utime = time(ia) and time = Utime(ic).
% calcul des solutions y de l'ED en t=Utime (sans répitition d'un ti)
[Utime,ia,ic]=unique(time);
[t y] = ode45(@(t,y)var_prot_mu_reparam2(t,y,param), Utime, param(1)); %# Obtain solutions at specific times
% F=valeur de y solution de l'ED en t=time (vecteur temps avec répétitions
% éventuelles)
F=y(ic);
