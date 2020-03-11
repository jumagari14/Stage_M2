function  courbe_niveau(X,kparg,Pg,cpt_fen)
%affichage des courbes de niveau des écarts en fonction de (ks,kd)
% autour du (kparg(11),kparg(12)) et de (kparg(2),kparg(3))
% valeur min de (ks,kd) et valeur bornée de (ks,kd)
global DPA_P
parI=kparg(10:12);par=kparg(1:3);
dist=abs(parI-par)
pas1=parI(2)/10;pas2=parI(3)/10;
% nombre de points de la grille nbgrid^2
nbgrid=50;edo=zeros(length(DPA_P),nbgrid,nbgrid);
ecart=zeros(nbgrid,nbgrid);
ksseq=linspace(parI(2)*0.1,parI(2)*2+dist(2),nbgrid);
kdseq=linspace(parI(3)*0.1,parI(3)*2+dist(3),nbgrid);
nPg=norm(Pg);
for in=1:nbgrid
    for jn=1:nbgrid
        xseq=[parI(1);ksseq(in);kdseq(jn)];
        edo(:,in,jn)=resol_mu(xseq,DPA_P);
        ecart(in,jn)=norm(edo(:,in,jn)-Pg)/nPg;
    end
end
figure(cpt_fen);
emin=min(min(ecart));pas=emin/20;
z=[emin;emin+pas;emin+2*pas;emin+4*pas;emin+8*pas;emin+16*pas];
contour(ksseq,kdseq,transpose(ecart),10,'LevelList',z,'LevelListMode','manual','ShowText','on')
xlabel('k_s ')
ylabel('k_d') %
hold on;
plot(kparg(11),kparg(12),'r*');
plot(kparg(5),kparg(6),'c*');
plot(kparg(2),kparg(3),'y*');
plot(kparg(8),kparg(9),'y*');
% ellipse de confiance
alpha=0.05;nl=length(Pg);nb_par=2;emin*nPg
seuil=emin*nPg*(1+nb_par/(nl-nb_par)*qf(1-alpha,nb_par,nl-nb_par))
contour(ksseq,kdseq,transpose(ecart),10,'LevelList',seuil,'LevelListMode','manual','ShowText','on')
% aire de l'ellipse de confiance
U=transpose(X)*X;det(U)
aire=emin*nPg*nb_par/(nl-nb_par)*qf(1-alpha,nb_par,nl-nb_par)*4*pi/sqrt(det(U))
end

