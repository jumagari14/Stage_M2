function fit_testMu( i,nom,xi,yi)
% différents fits du Poids et de mu (taux de croissance) : polynomial, estimation par noyau
global Wp6 SW6 mW6 Wp3 SW3 mW3 degreW h fitMu

% Lissage des données Poids: regression polynomiale d'ordre 3 de (DPA_i, R_i)
% paramètres polynomes (DPA centré, réduit)
[Wp3,SW3,mW3] = polyfit(xi,yi,3);
% valeurs du poly à xi (temps experimentaux)
[Wval3,dW3] = polyval(Wp3,xi,SW3,mW3);
[dWp3]=polyder(Wp3);
    [dWval3]=polyval(dWp3,xi,SW3,mW3)./mW3(2);
    mu1=dWval3./Wval3;
col='c';
[hAx1,hAx2]=affreg(i,nom,xi,yi,Wval3,col);

% fit de log(y): pour avoir des y>0
yi(1)=0.1;yi(2)=0.1;
[Wp6,SW6,mW6] = polyfit(xi,log(yi),3);
% valeurs du poly à yi (temps experimentaux)
[Wval6,dW6] = polyval(Wp6,xi,SW6,mW6);
[dWp6]=polyder(Wp6);
    [dWval6]=polyval(dWp6,xi,SW6,mW6)./mW6(2);
    mu2=dWval6
% affichage regression
col='b';
[hAx1,hAx3]=affreg(i,nom,xi,yi,exp(Wval6),col);

Wfit=ksr2(xi,yi,h);dWfit=dksr2(xi,yi,h);
mufit=dWfit./Wfit;
col='r';
[hAx1,hAx4]=affreg(i,nom,xi,yi,Wfit,col);

Wfit2=ksr2(xi,yi,h-2);dWfit2=dksr2(xi,yi,h-2);
mufit2=dWfit2./Wfit2;
col='m';
[hAx1,hAx5]=affreg(i,nom,xi,yi,Wfit2,col);

erreur=[norm(yi-Wval3,2),norm(yi-exp(Wval6),2),norm(yi-Wfit,2),norm(yi-Wfit2,2)]./norm(yi,2);
% fit par défaut: estimateur à noyau

%if choix==0
%      choix = input('fitR 0 3 6 ', 's');
%end
limaxis=axis;
xPos=limaxis(1);
yl = limaxis(4);
yPos = yl*.7;
texte=text(xPos,yPos,sprintf('er rel_{pol3,exp(pol3),noyau1,noyau2}=%8.4f %8.4f %8.4f %8.4f',erreur),'FontSize',8);
legend([hAx1,hAx2,hAx3,hAx4,hAx5],{'points exp.','poly deg3','expo poly deg3', strcat('noyau h=',num2str(h)),strcat('noyau h=',num2str(h-2))},'FontSize',5,'Location','northwest');
title('Regression Poids du fruit')
    xlabel('DPA ')
    ylabel('Poids')
    hold off;
    
    
%
[hAx21]=affMu(i,nom,xi,dWval3,'c');
[hAx22]=affMu(i,nom,xi,dWval6,'b');
[hAx23]=affMu(i,nom,xi,dWfit,'r');
[hAx24]=affMu(i,nom,xi,dWfit2,'m');
     hold off;
% Taux de croissance calculé à partir de la regression du poids
[hAx11]=affMu(i+50,nom,xi,mu1./10,'c');
 [hAx12]=affMu(i+50,nom,xi,mu2,'b');
   
    [hAx13]=affMu(i+50,nom,xi,mufit,'r');
    [hAx14]=affMu(i+50,nom,xi,mufit2,'m');
    
    limaxis=axis;
xPos=limaxis(1);
yl = limaxis(4);
yPos = yl*.7;
legend([hAx11,hAx12,hAx13,hAx14],{'poly deg3','exp(poly)', strcat('noyau h=',num2str(h)),strcat('noyau h=',num2str(h-2))},'FontSize',5,'Location','northeast');
title('Taux de croissance du fruit')
    xlabel('DPA ')
    ylabel('Mu')
    hold off;
    
end

