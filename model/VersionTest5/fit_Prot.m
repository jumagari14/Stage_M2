function fit_Prot(i,nom,DPA,P,kpar,DPA_fit,sol,resnorm,R0,P0,ksn,score)
% affichage fit de Proteines (quantité en gFw ou en fruit)/
% comparaison (ks,kd)

% coleurs successives courbes
 col=['c','b','r','m','g','k','y'];
% nombre de fit des proteines
nb_fit=numel(kpar)/3;


% pour chaque fit
switch nb_fit
    case 2
    %affichage regression
    % les courbes
    [hAx1,hAx2]=affreg(i,nom,DPA,P,DPA_fit,sol(:,1),col(1));
    [hAx1,hAx3]=affreg(i,nom,DPA,P,DPA_fit,sol(:,2),col(2));
    movegui('northeast')
    % texte: ks,kd et erreur relative
    limaxis=axis;
    xPos=limaxis(1);
    yl = limaxis(4);
    yPos = yl*.7;
    % calcul (ks,kd) sans normalisation
    ksnorm=[ksn, kpar(2),kpar(5)];
    ks=[ksn, kpar(2),kpar(5)].*P0/R0;
    kd=[0.1*ksn, kpar(3),kpar(6)];
    texte=text(xPos,yPos,sprintf('ks_{min, borne,sans}=%8.4f %8.4f %8.4f',ks),'FontSize',8);
    texte=text(xPos,yPos-yPos*0.1,sprintf('ks_{normalisé min, borne,sans}=%8.4f %8.4f %8.4f',ksnorm),'FontSize',8);
    texte=text(xPos,yPos-yPos*0.2,sprintf('kd_{initial, borne,sans}=%8.4f %8.4f %8.4f',kd),'FontSize',8);
    texte=text(xPos,yPos-yPos*0.3,sprintf('er rel_{borne,ssborne}=%8.4f %8.4f ',resnorm),'FontSize',8);
    texte=text(xPos,yPos-yPos*0.4,sprintf('RNA moy P moy score=%8.4f %8.4f %2d',R0,P0,score),'FontSize',8);
    legend([hAx1,hAx2,hAx3],{'points exp.','(ks,kd) bornes','(ks,kd) non bornés'},'FontSize',5,'Location','northwest');
    
    case 3
    [hAx1,hAx2]=affreg(i,nom,DPA,P,DPA_fit,sol(:,1),col(1));
    [hAx1,hAx3]=affreg(i,nom,DPA,P,DPA_fit,sol(:,2),col(2));
    [hAx1,hAx4]=affreg(i,nom,DPA,P,DPA_fit,sol(:,3),col(3));
    movegui('northwest')
    limaxis=axis;
    xPos=limaxis(1);
    yl = limaxis(4);
    yPos = yl*.7;
    ksnorm=[ksn, kpar(2),kpar(5),kpar(8)];
    ks=[ksn, kpar(2),kpar(5),kpar(8)].*P0/R0;
    kd=[ksn*0.1, kpar(3),kpar(6),kpar(9)];
    texte=text(xPos,yPos,sprintf('ks_{min, expl borne,num borne, sans}=%8.4f %8.4f %8.4f %8.4f',ks),'FontSize',8);
    texte=text(xPos,yPos-yPos*0.1,sprintf('ks_{normalisé min, expl borne,num borne,sans}=%8.4f %8.4f %8.4f %8.4f',ksnorm),'FontSize',8);
    texte=text(xPos,yPos-yPos*0.2,sprintf('kd_{initial,expl borne,num borne,sans}=%8.4f %8.4f %8.4f %8.4f',kd),'FontSize',8);
    texte=text(xPos,yPos-yPos*0.3,sprintf('er rel_{expl borne,num borne,ssborne}=%8.4f %8.4f %8.4f ',resnorm),'FontSize',8);
    texte=text(xPos,yPos-yPos*0.4,sprintf('RNA moy P moy score=%8.4f %8.4f %2d',R0,P0,score),'FontSize',8);
    legend([hAx1,hAx2,hAx3,hAx4],{'points exp.','explicite borne','num borne', 'sans borne'},'FontSize',5,'Location','northwest');
    
    case 4
    [hAx1,hAx2]=affreg(i,nom,DPA,P,DPA_fit,sol(:,1),col(1));
    [hAx1,hAx3]=affreg(i,nom,DPA,P,DPA_fit,sol(:,2),col(2));
    [hAx1,hAx4]=affreg(i,nom,DPA,P,DPA_fit,sol(:,3),col(3));
    [hAx1,hAx5]=affreg(i,nom,DPA,P,DPA_fit,sol(:,4),col(4));
    movegui('northwest')
    limaxis=axis;
    xPos=limaxis(1);
    yl = limaxis(4);
    yPos = yl*.7;
    ksnorm=[ksn, kpar(2),kpar(5),kpar(8),kpar(11)];
    ks=[ksn, kpar(2),kpar(5),kpar(8),kpar(11)].*P0/R0;
    kd=[ksn*0.1, kpar(3),kpar(6),kpar(9),kpar(12)];
    texte=text(xPos,yPos,sprintf('ks_{min, sans borne}=%8.4f %8.4f %8.4f %8.4f %8.4f',ks),'FontSize',8);
    texte=text(xPos,yPos-yPos*0.1,sprintf('ks_{normalisé min, sans}=%8.4f %8.4f %8.4f %8.4f %8.4f',ksnorm),'FontSize',8);
    texte=text(xPos,yPos-yPos*0.2,sprintf('kd_{initial,sans}=%8.4f %8.4f %8.4f %8.4f %8.4f',kd),'FontSize',8);
    texte=text(xPos,yPos-yPos*0.3,sprintf('er rel=%8.4f %8.4f %8.4f %8.4f',resnorm),'FontSize',8);
    texte=text(xPos,yPos-yPos*0.4,sprintf('RNA_{moy}=%8.4f P_{moy}=  %8.4f ',R0,P0),'FontSize',8);
    texte=text(xPos,yPos-yPos*0.5,sprintf('scoreO=%2d scoreF=%2d ',score(1),score(2)),'FontSize',8);
    legend([hAx1,hAx2,hAx3,hAx4,hAx5],{'points exp.','optim1','optim2','optim3', 'optim4'},'FontSize',5,'Location','northwest');
    
end
title(['Regression PROT ' nom])
xlabel('DPA ')
ylabel('PROT')

% courbe pour publi
figure(5002);
col='r';
 [hAx1,hAx2]=affreg(5002,nom,DPA,P,DPA_fit,sol(:,1),col);
title(['\fontsize{16}' nom ': Protein profile and resolution']);
xlabel('\fontsize{14} Time (Days Post Anthesis)');
ylabel('\fontsize{14} Protein normalized by the mean');



end

