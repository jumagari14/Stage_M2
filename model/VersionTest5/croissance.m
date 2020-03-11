function [RNAf,PROTf]=croissance(RNA,DPA_R,PROT,DPA_P)
% prise en compte de la croissance: quantit� dans le fruit
% A partir des donn�es brutes (RNA,PROT) quantit�s par gFW, et de la fonction poids du fruit
% renvoi les quantit�s (RNAf,PROTf) par fruit.

RNAf=RNA.*Poids(DPA_R);
PROTf=PROT.*Poids(DPA_P);

end

