function [P0,Pmin,Pmax] = initiale( DPA,P)
% intialisation de la condition initiale de l'EDO
% calcul de la valeur moyenne aux premiers DPA
DPA_min=min(DPA);
indice=find(DPA==DPA_min);
if numel(indice)==1
    % un seul point à DPA_min
    P0= P(indice);
    Pmax=5*P0;
    Pmin=0.1*P0;
else %plusieurs points
    P0= mean(P(indice));
    Pmin=min(P(indice))-std(P)*2;
    Pmax=max(P(indice))+std(P)*2;
    
end
end

