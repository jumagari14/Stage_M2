function estimateur_noyau( nom,x,y,h)
%estimatuer à noyau
% largeur de bande 
 
rf=ksr2(x,y,h);
[hAx1] = [plot(x,[rf]),'c-'];hold on;
[hAx2] = [scatter(x,y,'b*')];
    
    %set(hAx(2),'LineStyle','-','Color','c');
    legend([hAx1,hAx2],{'points exp.','reg '},'FontSize',5,'Location','northwest');
    title('Regression par estimateur noyau')
    xlabel('x ')
    ylabel('f(x)')
       