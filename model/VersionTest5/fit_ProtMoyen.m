function [sol]=fit_ProtMoyen(DPA,kparbm)
% calcul Prot à partir de la valeur moyenne des bootstrap
sol=zeros(length(DPA),4);

sol(:,1)=resol_mu(kparbm(1:3),DPA);
sol(:,2)=resol_mu(kparbm(4:6),DPA);
%sol(:,3)=resol_mu(kparbm(7:9),DPA);
%sol(:,4)=resol_mu(kparbm(10:12),DPA);
end
