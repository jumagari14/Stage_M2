function [ yv2 ] = resol_contois( par,t)
%% r�solution ED en fonction des param�tres param au temps=time
% param�tres  (r,K,R,y0)
% ode45(fun,[t0 tf],y0)
% Utime: returns the same data as in time, but with no repetitions.
%  Utime = time(ia) and time = Utime(ic).
% calcul des solutions y de l'ED en t=Utime (sans r�pitition d'un ti)
r=par(1);K=par(2);R=par(3);y0=par(4);
dcontois=@(t,y,r,R,K)(r*y.*(1-y/K)./(K+(R-1)*y));
[Utime,ia,ic]=unique(t);
[t yv2] = ode45(@(t,y)dcontois(t,y,r,R,K), Utime, y0); %# Obtain solutions at specific times
% F=valeur de y solution de l'ED en t=time (vecteur temps avec r�p�titions
% �ventuelles)
yv2=yv2(ic);
end

