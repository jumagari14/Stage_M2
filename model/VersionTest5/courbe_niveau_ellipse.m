function  [c1,seuil2,Aire,score,message]=courbe_niveau_ellipse(nom,X,kparg,errg,Pg,cpt_fen)
%affichage des courbes de niveau des écarts en fonction de (ks,kd)
% autour du (kparg(11),kparg(12)) et de (kparg(2),kparg(3))
% valeur min de (ks,kd) et valeur bornée de (ks,kd)
global DPA_P
parI=kparg(:,1);
%dist=abs(parI-par)
%pas1=parI(2)/10;pas2=parI(3)/10;
% nombre de points de la grille nbgrid^2
nbgrid=50;edo=zeros(length(DPA_P),nbgrid,nbgrid);
ecart=zeros(nbgrid,nbgrid);
%ksseq=linspace(parI(2)*0.1,parI(2)*2+dist(2),nbgrid);
ksseq=linspace(parI(2)*0.1,parI(2)*5,nbgrid);
%kdseq=linspace(parI(3)*0.1,parI(3)*2+dist(3),nbgrid);
kdseq=linspace(parI(3)*0.1,max(parI(3)*5,0.1),nbgrid);
%nPg=norm(Pg);
for in=1:nbgrid
    for jn=1:nbgrid
        xseq=[parI(1);ksseq(in);kdseq(jn)];
        edo(:,in,jn)=resol_mu(xseq,DPA_P);
        ecart(in,jn)=norm(edo(:,in,jn)-Pg);
    end
end
figure(cpt_fen);
emin=min(min(ecart));pas=emin/20;
%z=[emin;emin+pas;emin+2*pas;emin+4*pas;emin+8*pas;emin+16*pas];
%contour(ksseq,kdseq,transpose(ecart),10,'LevelList',z,'LevelListMode','manual','ShowText','on')
xlabel('k_s ')
ylabel('k_d') %

hold on;
plot(kparg(2,1),kparg(3,1),'r*');
%plot(kparg(5),kparg(6),'c*');
%plot(kparg(2),kparg(3),'y*');
%plot(kparg(8),kparg(9),'y*');
% ellipse de confiance
alpha=0.1;nl=length(Pg);nb_par=2;emin
seuil1=emin*(1+nb_par/(nl-nb_par)*qf(1-alpha,nb_par,nl-nb_par))
seuil2=emin*(1+nb_par/(nl-nb_par)*qf(1-2.5*alpha,nb_par,nl-nb_par))
seuil3=emin*(1+nb_par/(nl-nb_par)*qf(1-5*alpha,nb_par,nl-nb_par))
z=[seuil1;seuil2;seuil3]

title(['Domaine de Confiance ' nom],'FontWeight','bold');
[c2,hc2]=contour(ksseq,kdseq,transpose(ecart),10,'LevelList',z,'LevelListMode','manual','ShowText','on');
hcont = get(hc2,'children');
%legend('\alpha=10%','\alpha=25%','\alpha=50%','FontSize',5,'Location','northwest');

% aire de l'ellipse de confiance
U=transpose(X(:,2:3))*X(:,2:3);
d2=emin*(nb_par/(nl-nb_par)*qf(1-2.5*alpha,nb_par,nl-nb_par));
%[Vect,Diag]=eig(U)
lg_el=sqrt(d2./eig(U));
aire=d2*pi/sqrt(det(U));
score=0.1/aire;
limaxis=axis;
xPos=limaxis(1);
yl = limaxis(4);
yPos = yl*.7;
texte=text(xPos,yPos,sprintf('ks=%8.4f kd=%8.4f',kparg(2,1),kparg(3,1)),'FontSize',8);
yPos = yl*.6;
texte=text(xPos,yPos,sprintf('aire_{th}=%8.4f a=%8.4f b=%8.4f score_{th}=%2d',aire,lg_el,floor(score)),'FontSize',8);

% calcul approx de l'aire
% si rho n'est pas trop proche de 1
rho=matrice_correlation(X,numel(Pg),(errg(1)*norm(Pg,2))^2)
% test sur rho + si kd<>0 =>courbe fermée
c1=[];
%if rho(2,3)>0.99 % très grande correlation entre ks et kd
% le domaine de confiance n'est pas fermé
%    score=0;Aire=999;
%else %le domaine de confiance est fermé

[c1,h]=contour(ksseq,kdseq,transpose(ecart),10,'LevelList',[seuil2 seuil2],'LevelListMode','manual','ShowText','on');
% calcul de l'aire
[Aire,message]=detection_aire(c1);

if Aire>1e-2
    score=0.1/Aire;
    
else
    score=20;
end
yPos = yl*.5;
texte=text(xPos,yPos,sprintf('aire calculée=%8.4f Score=%2d',Aire,floor(score)),'FontSize',8);
yPos = yl*.4;
texte=text(xPos,yPos,sprintf(message),'FontSize',8);
%end
end

