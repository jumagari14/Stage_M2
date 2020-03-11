function  [err,message]=fit_testpoids( xi,yi )
% essai de differents fits pour le poids du fruit
% choix de l'utilisateur

global fitPoids parv parg parc parl pare pard Wp3 Sp3 mp3 h
% logistique
% dy/dx=ry(1-y/K) + y(0)=y0
% y(x)=Ky0/(y0+(K-y0)exp(-rx)
% gompertz
% dy/dx=ry ln(K/y)
% y(x)=K exp(b exp(-ax)) avec b=ln(y0/K)
% logistique généralisée
% dy/dx=ry(1-[y/K]^a)
% contois
% dy/dx=ry(1-y/K)/(K+(RK'-1)y) si Rk'<1 inflexion <K/2

figure(1);
% points expérimentaux (seuls)
hAx=scatter(xi,yi);

figure(2);
% points expérimentaux et fit
hAx=scatter(xi,yi);normyi=norm(yi);
% croissance exponentielle: modèle de Malthus
% malthus=@(t,r,y0)(y0 * exp(r.*t));

% croissance logistique: modèle de Verhulst
verhulst=@(par,t)(par(2)*par(3)./(par(3)+(par(2)-par(3))*exp(-par(1).*t)));
parv=lsqcurvefit(verhulst,[0.1,100,1],xi,yi);
yv=verhulst(parv,xi);% y0=CI à t=0 DPA
err_ver=norm(yv-yi)/normyi;
hold on; hAx2=plot(xi,yv,'r-');

% test ODE (attention CI à 7 DPA)
r=0.11;K=101.8;y1=7.08;
dverhulst=@(t,y,r,K)(r.*y.*(1-y./K));
% Verhulst détermination paramètre 
% avec solution numérique ED verhulst y0=CI ici à t=7DPA
r0=0.1;K0=100;par=[0.1;100;7];
[parv2,resnormv,residualv,exitflagv,outputv]  = lsqcurvefit(@resol_verhulst,par,xi,yi);
yv2=resol_verhulst(parv2,xi);% y0=CI à t=7DPA


% Croissance : modèle de Gompertz
gompertz=@(par,t)(par(2)*exp(log(par(3)/par(2)).*exp(-par(1).*t)));
parg=lsqcurvefit(gompertz,[0.065,114.39,0.52],xi,yi);
yg=gompertz(parg,xi);
err_gom=norm(yg-yi)/normyi;
hold on; hAx3=plot(xi,yg,'c-');

% Modèle de Contois
r=0.5;K=110;R=2;y0=0.5;par=[5;2;100;1];lb=[0;0;90;0.05];ub=[inf;inf;130;2];
dcontois=@(t,y,r,R,K)(r.*y*(1-y/K)/(K+(R-1)*y)) ;
[parc,resnormc,residualc,exitflagc,outputc]  = lsqcurvefit(@resol_contois,par,xi,yi,lb,ub);
yc=resol_contois(parc,xi);
err_con=norm(yc-yi)/normyi;
hold on; hAx4=plot(xi,yc,'y-');

% logistique généralisée
%logisgen=@(t,r,K,y0)(K*...);
%yl=logisgen(xi,0.065,114.39,0.52)
par=[110.8;0;13.7;0.144];
[parl,resnorml,residuall,exitflagl,outputl]  = lsqcurvefit(@resol_logisgen,par,xi,yi);
yl=resol_logisgen(parl,xi);
err_lg=norm(yl-yi)/normyi;
hold on; hAx5=plot(xi,yl,'m:');



% Modèle empirique: log(y(t))=V(t-7)/(k+t-7)
empirique=@(par,t)(exp(par(1).*(t-par(3))./(par(2)+t-par(3))));
pare=lsqcurvefit(empirique,[5.38,8,7],xi,yi);
ye=empirique(pare,xi);
hold on; hAx6=plot(xi,ye,'k-');
err_emp=norm(ye-yi)/normyi;

% estimateur à noyau h=6
h=6;
yn=ksr2(xi,xi,yi,h);       
err_noy=norm(yn-yi)/normyi;
hold on; hAx7=plot(xi,yn,'m:');
% Chapman-Richards
% K(1-exp(-2*r*t/K))^2;

% Modèle  log(y(t))=poly degre 3
[Wp3,Sp3,mp3] = polyfit(xi,log(yi),3);
% valeurs du poly à yi (temps experimentaux)
[ylp,dW3] = polyval(Wp3,xi,Sp3,mp3);ylp=exp(ylp);
hold on; hAx8=plot(xi,ylp,'k:');
err_lp=norm(ylp-yi)/normyi;

xlabel('DPA ')
ylabel('Poids gFW') % y-axis
title('Test fit poids tomate',...
        'FontWeight','bold')
legend([hAx,hAx2,hAx3,hAx4,hAx5,hAx6,hAx7,hAx8],{'exp.','logistic','Gompertz','Contois','logis. gen','empirique','noyau h=6','log poly'},'FontSize',5,'Location','northwest');
limaxis=axis;
xPos=limaxis(1);
yl = limaxis(4);
yPos = yl*.7;
texte=text(xPos,yPos,sprintf('err_{Ver,Gom,Con,LG,emp,noyau,logpoly}=%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f',err_ver,err_gom,err_con,err_lg,err_emp,err_noy,err_lp),'FontSize',8);


figure(3);
% points expérimentaux et autres fits
hAx=scatter(xi,yi);

% Croissance : simple sigmoide
sigmoide=@(par,t)(par(4)+par(1)./(1+exp(-par(2).*(t-par(3)) )));
lb=[90;0;25;0.1];ub=[120;inf;35;2];
pars=lsqcurvefit(sigmoide,[110;0.1;30;1],xi,yi,lb,ub);
ys=sigmoide(pars,xi);
err_sig=norm(ys-yi)/normyi;
hold on; hAx9=plot(xi,ys,'c--');

% Croissance : double sigmoide
sigmoide=@(par,t)(par(4)+par(1)./(1+exp(-par(2).*(t-par(3)) ))...
    +par(5)./(1+exp(-par(6).*(t-par(7)) )));
lb=[20;0;15;0.1;60;0;40];ub=[70;inf;25;2;110;inf;55];
pard=lsqcurvefit(sigmoide,[40;1;20;1;100;2;45],xi,yi,lb,ub);
yd=sigmoide(pard,xi);
err_sid=norm(yd-yi)/normyi;
hold on; hAx10=plot(xi,yd,'r--');


% Chapman-Richards
% K(1-exp(-2*r*t/K))^2;
% fonction de Hill
% V*t^n/(k+t^n)

xlabel('DPA ')
ylabel('Poids gFW') % y-axis
title('Test fit poids tomate',...
        'FontWeight','bold')
legend([hAx,hAx9,hAx10],{'exp.','sigmoide','double sigmoide'},'FontSize',5,'Location','northwest');
limaxis=axis;
xPos=limaxis(1);
yl = limaxis(4);
yPos = yl*.7;
texte=text(xPos,yPos,sprintf('err_{sig,doublesig}= %8.4f %8.4f',err_sig,err_sid),'FontSize',8);

%%%%%%%%
% courbe pour publication: double sigmoôide
figure(5000);
% points expérimentaux
hAx=scatter(xi,yi);xlabel('\fontsize{14} Time (Days Post Anthesis)');
ylabel('\fontsize{14} Weight (gFW)');title('\fontsize{16} Tomato weight curve fit',...
        'FontWeight','bold') %
% Croissance : double sigmoide
hold on; hAx10=plot(xi,yd,'r--');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fit des données avec un point 1/sqrt(yi)
figure(2);
hAx=scatter(xi,yi);
% ajout d'un point inversement proportionnel à la racine carrée du poids
weights=1./sqrt(yi);

% croissance logistique avec poids
verw=@(par,t)weights.*(par(2)*par(3)./(par(3)+(par(2)-par(3))*exp(-par(1).*t)));
parvw=lsqcurvefit(verw,[0.1,100,1],xi,yi.*weights);
yvw=verw(parvw,xi)./weights;% 
err_verp=norm(yi-yvw)/normyi;
hold on; hAx2=plot(xi,yvw,'r-');


% Gompertz 
%weights=1./yi;lb=[0,0,0];ub=[10,150,10];
%gompw=@(par,t)weights.*(par(2)*exp(log(par(3)/par(2)).*exp(-par(1).*t)));
%pargw=lsqcurvefit(gompw,[0.065,114.39,0.52],xi,yi.*weights,lb,ub);
%ygw=gompw(pargw,xi)./weights;

% Gompertz avec poids 
lb=[0,0,0];ub=[10,150,10];
gompw=@(par,t)weights.*(par(2)*exp(log(par(3)/par(2)).*exp(-par(1).*t)));
pargw=lsqcurvefit(gompw,[0.065,114.39,0.52],xi,yi.*weights,lb,ub);
ygw=gompw(pargw,xi)./weights;
err_gomp=norm(yi-ygw)/normyi;
hold on; hAx3=plot(xi,ygw,'c-');


% à faire............

% Contois avec poids
r=0.5;K=110;R=2;y0=0.5;par=[5;2;100;1];lb=[0;0;90;0.05];ub=[inf;inf;130;2];
%dcontois=@(t,y,r,R,K)(r.*y*(1-y/K)/(K+(R-1)*y)) 
%[parc,resnormc,residualc,exitflagc,outputc]  = lsqcurvefit(@resol_contois,par,xi,yi,lb,ub);
%yc=resol_contois(parc,xi);
%hold on; hAx4=plot(xi,yc,'y-');

% logistique généralisée
%logisgen=@(t,r,K,y0)(K*...);
%yl=logisgen(xi,0.065,114.39,0.52)
%par=[110.8;0;13.7;0.144];
%[parl,resnorml,residuall,exitflagl,outputl]  = lsqcurvefit(@resol_logisgen,par,xi,yi);
%yl=resol_logisgen(parl,xi);
%hold on; hAx5=plot(xi,yl,'c-');


xlabel('DPA ')
ylabel('Poids gFW') % y-axis
title('Test fit poids tomate Pondéré',...
        'FontWeight','bold')
legend([hAx,hAx2,hAx3],{'exp.','logistic Pondéré','Gompertz Pondéré'},'FontSize',5,'Location','northwest');
limaxis=axis;
xPos=limaxis(1);
yl = limaxis(4);
yPos = yl*.7;
texte=text(xPos,yPos,sprintf('err_{VerP,GomP}=%8.4f %8.4f ',err_verp,err_gomp),'FontSize',8);

% calcul du taux de croissance
figure(4);
fitPoids=1;
hold on; hAx1=plot(xi,mu(xi),'r-');
fitPoids=2;
hold on; hAx2=plot(xi,mu(xi),'c-');
fitPoids=3;
hold on; hAx3=plot(xi,mu(xi),'y-');
fitPoids=5;
hold on; hAx5=plot(xi,mu(xi),'k--');
fitPoids=6;
hold on; hAx6=plot(xi,mu(xi),'m:');
fitPoids=7;
hold on; hAx7=plot(xi,mu(xi),'b:');
fitPoids=8;
hold on; hAx8=plot(xi,mu(xi),'b--');
xlabel('DPA ')
ylabel('taux') % y-axis
title('Test fit taux de croissance',...
        'FontWeight','bold')
legend([hAx1,hAx2,hAx3,hAx5,hAx6,hAx7,hAx8],{'logistic','Gompertz','Contois','empirique','noyau','log poly','double sig'},'FontSize',5,'Location','northwest');

% choix du fit par l'utilisateur
fitPoids=0;
while fitPoids==0
    prompt = 'Fit 1-logistic 2-Gompertz 3-Contois  5-empirique 6-noyau 7-log poly 8-double sigmoide';
    dlg_title = 'Choix courbe croissance';
    num_lines = 1;
    defaultans = {'8'};
    choix = inputdlg(prompt,dlg_title,num_lines,defaultans);
    
    
    switch choix{1}
        case '1'
            fitPoids=1;err=err_ver;message='logistic';%logistic
        case '2'
            fitPoids=2;err=err_gom;message='Gompertz';%Gompertz
        case '3'
            fitPoids=3;err=err_con;message='Contois';%Contois
        %case '4'
         %   fitPoids=4;err=err_lg;message='log généralisé';% log gene
        case '5'
            fitPoids=5;err=err_emp;message='empirique';%empirique
        case '6'
            fitPoids=6;err=err_noy;message='estimateur à noyau';%noyau
        case '7'
            fitPoids=7;err=err_lp;message='exp(poly deg 3';%log poly
        case '8'
            fitPoids=8;err=err_sid;message='double sigmoide';%double sigmoide
        otherwise
            fitPoids=0;%
            disp('indisponible');
            
    end;
end;

end

