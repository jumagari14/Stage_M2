function affsol(i,ks_n,kd,par,nom,DPA,P,sol)

global ptrans3 ptrans6 St6 mt6 St3 mt3 degre

parI=par(1,:);par1=par(2,:);par2=par(3,:);parL=par(4,:);parMu=par(5,:);
solI=sol(:,1);sol1=sol(:,2);sol2=sol(:,3);solL=sol(:,4);solMu=sol(:,5);
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
h=figure(i+100);

[hAx] = [scatter(DPA,P,'r*')];hold on;
% ks/kd*R 
Trans=par2(2)/par2(3)*mRNA_deg(DPA);
[hAx2] = plot(DPA,[sol1,sol2,solI,solL,Trans,solMu]);
    
    set(hAx2(1),'LineStyle','-','Color','r');
    set(hAx2(2),'LineStyle','--','Color','c');
    set(hAx2(3),'LineStyle','--','Color','g');
     set(hAx2(4),'LineStyle','-','Color','m');
     set(hAx2(5),'LineStyle','-','Color','b');
     set(hAx2(6),'LineStyle','--','Color','k');
     legend([hAx,hAx2(1),hAx2(2),hAx2(3),hAx2(4),hAx2(5),hAx2(6)],{'points exp.','sol ED num', 'sol ED expl','sol SC','sol ED (ss bornes)','trans*ks/kd','avec mu'},'FontSize',5,'Location','northwest');
    % set(hAx(3),'LineStyle',':','Color','b');
     title(nom,...
        'FontWeight','bold')
    xlabel('DPA ')
    ylabel('Proteines')
    limaxis=axis;
    %xl = limaxis(2);
    %xPos = xl + (xl-limaxis(1)) / 2;
    xPos=limaxis(1);
    yl = limaxis(4);
    yPos = yl*.7;
    texte=text(xPos,yPos,sprintf('ks_{sc,num,expl,ssb,mu}=%8.4f %8.4f %8.4f %8.4f %8.4f',ks_n,par1(2),par2(2),parL(2),parMu(2)),'FontSize',8);
    texte=text(xPos,yPos-yPos*0.1,sprintf('kd_{sc,num,expl,ssb,mu}=%8.4f %8.4f %8.4f %8.4f %8.4f',kd,par1(3),par2(3),parL(3),parMu(3)),'FontSize',8);
    texte=text(xPos,yPos-yPos*0.15,sprintf('ks/kd_{sc,num,expl,ssb,mu}=%8.4f %8.4f %8.4f %8.4f %8.4f',ks_n/kd,par1(2)/par1(3),par2(2)/par2(3),parL(2)/parL(3),parMu(2)/parMu(3)),'FontSize',8);
    %set(texte, 'HorizontalAlignment','BackgroundColor',[.7 .9 .7],'FontSize',8);
    %yPos2 = yl*.7;
    erreur=[norm(P-solI,2),norm(P-sol1,2),norm(P-sol2,2),norm(P-solL,2),norm(P-solMu,2)]./norm(P,2);
     texte=text(xPos,yPos-yPos*0.2,sprintf('er rel_{sc,num,expl,ssb,mu}=%8.4f %8.4f %8.4f %8.4f %8.4f',erreur),'FontSize',8);
  %texte=text(xPos,yPos2,sprintf('err=%8.4f ou %8.4f (%8.4f) ','FontSize',8);
    
end

