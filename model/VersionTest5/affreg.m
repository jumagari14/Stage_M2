function [hAx1,hAx2]=affreg(i,nom,DPA_R,R,DPA_fit,rpoly,col)
% affichage regression
%   
figure(i);
[hAx1] = [scatter(DPA_R,R,'b*')];hold on;
[hAx2] = [plot(DPA_fit,[rpoly])];
 
set(hAx2,'LineStyle','-','Color',col);
    %legend([hAx1,hAx2],{'points exp.',leg},'FontSize',5,'Location','northwest');
    


end

