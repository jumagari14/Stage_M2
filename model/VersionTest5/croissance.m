function [RNAf,PROTf]=croissance(RNA,DPA_R,PROT,DPA_P)
% prise en compte de la croissance: quantité dans le fruit
% A partir des données brutes (RNA,PROT) quantités par gFW, et de la fonction poids du fruit
% renvoi les quantités (RNAf,PROTf) par fruit.

RNAf=RNA.*Poids(DPA_R);
PROTf=PROT.*Poids(DPA_P);

end

