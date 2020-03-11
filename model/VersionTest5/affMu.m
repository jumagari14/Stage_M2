function [hAx]=affMu(i,nom,x,y,col)
% affichage regression du taux de croissance
%   
figure(60+i);
[hAx] = plot(x,y);hold on;
    
set(hAx,'LineStyle','-','Color',col);
    %legend([hAx1,hAx2],{'points exp.',leg},'FontSize',5,'Location','northwest');
    


end

