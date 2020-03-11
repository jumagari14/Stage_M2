function [ yv2 ] = resol_verhulst( par,xi)
%% résolution ED en fonction des paramètres param au temps=time
% paramètres  (condition initiale ED, ks et kd)
% ode45(fun,[t0 tf],y0)
% Utime: returns the same data as in time, but with no repetitions.
%  Utime = time(ia) and time = Utime(ic).
% calcul des solutions y de l'ED en t=Utime (sans répitition d'un ti)
r=par(1);K=par(2);y0=par(3);dverhulst=@(t,y,r,K)(r.*y.*(1-y./K));
[Utime,ia,ic]=unique(xi);
[t yv2] = ode45(@(t,y)dverhulst(t,y,r,K), Utime, y0); %# Obtain solutions at specific times
% F=valeur de y solution de l'ED en t=time (vecteur temps avec répétitions
% éventuelles)
yv2=yv2(ic);
end

