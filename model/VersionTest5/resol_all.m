function F = resol_all(param,time)
% r�solution ED en fonction des param�tres param au temps=time
% param�tres  (condition initiale ED, ks et kd)
% ode45(fun,[t0 tf],y0)
% Utime: returns the same data as in time, but with no repetitions.
%  Utime = time(ia) and time = Utime(ic).
% calcul des solutions y de l'ED en t=Utime (sans r�pitition d'un ti)
[Utime,ia,ic]=unique(time);
[t y] = ode45(@(t,y)var_prot_all(t,y,param), Utime, param(1)); %# Obtain solutions at specific times
% F=valeur de y solution de l'ED en t=time (vecteur temps avec r�p�titions
% �ventuelles)
F=y(ic);
%F=y;
