function [ yl ] = resol_logisgen( par,xi)
% r�solution logistique g�n�ralis�e (Pav� p164)
% 
% Utime: returns the same data as in time, but with no repetitions.
%  Utime = time(ia) and time = Utime(ic).
% calcul des solutions y de l'ED en t=Utime (sans r�pitition d'un ti)
p1=par(1);p2=par(2);p3=par(3);p4=par(4);
logisgen=@(t,p1,p2,p3,p4)(p1.*(1+exp((p2-t)/p3)).^(-1/p4));
yl=logisgen(xi,p1,p2,p3,p4);

end

