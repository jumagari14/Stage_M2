function [ R,P,ks_n,r0,p0 ] = normalisation(i,nom,DPA_R, r,DPA_P,p,ks )
% normalisation et affichage des donn�es

%  Normalisation des donn�es par valeur moyenne
r0=mean(r);p0=mean(p);
R=r/r0;P=p/p0;
% normalisation ks
ks_n=ks*r0/(p0);
% affichage donn�es normalis�es
affdata(i,nom,DPA_R,R,DPA_P,P)

end

