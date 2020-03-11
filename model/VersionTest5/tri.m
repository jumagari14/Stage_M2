function [ PROT ] = tri( RNA,PROT )
% tri de la structure PROT pour l'ordonner comme RNA à partir du nom
%  
% Nbre de Proteines (1ere colonne=DPA)
nb_Prot=size(PROT.data,2)-1;
if size(RNA.data)==size(PROT.data)
% 1ère colonne: DPA des mesures
    PROTORD.data(:,1)=PROT.data(:,1);
for i=1:nb_Prot 
    index=find(strcmp(PROT.NOM,RNA.NOM(i)));
    PROTORD.data(:,i+1)=PROT.data(:,index);
    PROTORD.NOM(:,i)=PROT.NOM(:,index);
end
PROT=PROTORD;

end

