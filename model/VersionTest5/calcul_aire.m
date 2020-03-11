function [ Aire,message ] = calcul_aire(x1,x2)
% calcul de l'aire de la surface d�limit�e par (x1,x2)
% d�tection des points fronti�re "haute" 
% Aire1=aire en ces points et l'axe des abscisses
% Aire2=aire en les points fronti�re "basses" et axe des abs.
% Aire=Aire1-Aire2
if x1(1)==x1(end) &  x2(1)==x2(end)
    disp('courbe fermee')
    message='courbe fermee';
else
    disp('courbe non fermee')
    message='courbe non fermee';
end
% abscisse min et max du contour
idmax=find(x1==max(x1));idmin=find(x1==min(x1));
[x1sort,idsort]=sort(x1(idmin:idmax));
x2sort=x2(idsort);
% trac� des points du contour entre xmin et xmax
scatter(x1,x2,'red')
Aire1=trapz(x2sort,x1sort)
% suppression des points du contour entre xmin et xmax
x2(idmin:idmax)=[];x1(idmin:idmax)=[];
% ordonner les points restants du plus petit x1 au plus grand
[x1sort,idsort]=sort(x1);
x2sort=x2(idsort);
% trac� des points restants
scatter(x1,x2,'green')
if length(x2sort)>1
    Aire2=trapz(x2sort,x1sort)
else
    Aire2=0;
end;
Aire=abs(abs(Aire1)-abs(Aire2))

end

