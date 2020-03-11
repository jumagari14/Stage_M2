function affdata(i, nom, DPA_R,R,DPA_P,P )
% affichage des donn�es
h=figure(i);


hAx1=scatter(DPA_P,P,'r+');
hold on;hAx2=scatter(DPA_R,R,'b*');
title(nom,'FontWeight','bold')
xlabel('DPA ')
ylabel('Conc P/R fmol normalis�e par moyenne') % y-axis

legend([hAx1,hAx2],{'P','R'},'FontSize',5,'Location','northeast');

    