function [Aire_cal,message]=detection_aire(c1);
% calcul de l'aire déterminée par le cotour c1
% si la courbe n'est pas une unique courbe fermée
% calcul de l'aire correspondant au plus grand nombre de points contenus
% dans cle contour c1
% c1: résultat de la fonction contour
% (niveau,k1=nbre_de points) coordonnées (xi,yi) des k1 points
% (niveau,k2=nbre_de points) coordonnées (xi,yi) des k2 points

% lecture du contour
% (niveau,nb_poits) (xi,yi) (niveau,nb_poits)(xi,yi)...
x1=c1(1,:);
x2=c1(2,:);

lg=length(x1);
pt_lecture=1;
i=1;
niv=[];rang=[];
% test s'il y a plusieurs courbes dans contour
if lg==x2(1)+1
    % calcul de l'aire du contour
    [Aire_cal,message]=calcul_aire(x1(2:end),x2(2:end));
else
    % détection dans (x1,x2) des entêtes (niveau,nb_poits) de chaque courbe
    while pt_lecture<lg
        % chaque courbe commence par (niveau,nb_poits)
        niv(i)=x1(pt_lecture);
        nb(i)=x2(pt_lecture);
        rang(i)=pt_lecture;
        pt_lecture=pt_lecture+nb(i)+1; % pointe sur courbe suivante
        i=i+1;
    end
    % recherche de la courbe ayant le plus gand nombre de points
    % a priori la courbe d'aire la plus grande
    j=find(nb==max(nb)); % j=indice de nb où nb(j) est le plus grand
    % vérification que cette courbe est fermée
    idebut=rang(j)+1;ifin=rang(j)+nb(j);
    
    % calcul de l'aire
    [Aire_cal,message]=calcul_aire(x1(idebut:ifin),x2(idebut:ifin));
end;

end

