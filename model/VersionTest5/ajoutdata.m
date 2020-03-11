function [ PROT,RNA,poids] = ajoutdata(PROT,RNA,poids);
% ajout manuel de données
% proteines : ajout de (0,0) (0,0) (0,0) 
nb_Prot=size(PROT.data,2);
data0=zeros(3,nb_Prot)
PROT.data=[data0;PROT.data];
% RNA : ajout de (0,0) (0,0) (0,0) 
nb_RNA=size(RNA.data,2);
data0=zeros(3,nb_RNA)
RNA.data=[data0;RNA.data];
% poids ajout de 10 (0,0)
data0=zeros(10,2);
poids=[data0;poids];


