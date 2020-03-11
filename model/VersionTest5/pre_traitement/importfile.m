function [CONC]=importfile(fileToRead1)
% sous programme de pre_RNA22PROT
% importer les donn�es depuis le fichier txt fileToRead1
% le fichier txt ne doit pas contenir de NA
% les donnees seront class�s par ordre alphab�tique selon leur nom (vars{2})
% les valeurs NaN 999999 seront supprim�es


% Import the file
newData1 = importdata(fileToRead1);

% Noms des champs de la structure
vars = fieldnames(newData1);

% s'il y a deux champs, on continue
if length(vars)>=2
    % si les nombres d'�l�m�nts sont �gaux
    % nb colonnes de data=nb de textdata
    %if numel(newData1.(vars{1})(1,:))== numel(newData1.(vars{2}))
        % tri selon le champ vars{2}
        % (ordre alphab�tique du nom des concentrations)
        [pasbesoin,idx]=sort(newData1.(vars{2}));
        newData1.(vars{2})=newData1.(vars{2})(idx);
        newData1.(vars{1})=newData1.(vars{1})(:,idx);
        % initialisation structure: trois champs NOM, DPA et data
        CONC=struct('NOM',[],'DPA',[],'data',[]);
        % nb_prot+1 (1�re colonne=DPA)
        nb_prot=length(newData1.(vars{1}));
        
        % remplir la structure CONC � partir du fichier
        % compteur CONC sauv�es
        cpt_c=1;
        for i = 2:nb_prot
            CONC(cpt_c).NOM=newData1.(vars{2})(i);
            CONC(cpt_c).DPA=newData1.(vars{1})(:,1);
            CONC(cpt_c).data=newData1.(vars{1})(:,i);
            % d�tection des NaN ou 999999 ou 0 et suppression de ces donn�es
            index=find(isnan(CONC(cpt_c).data)|CONC(cpt_c).data==999999|CONC(cpt_c).data==0);
            %if length(index)>7 % moins de 20 donn�es, on supprime ce couple
            %    disp('suppression: nombre donn�es insuffisantes');
            %    CONC(cpt_c).NOM
            %    CONC(cpt_c).NOM=[];
            %    CONC(cpt_c).DPA=[];
            %    CONC(cpt_c).data=[];
            %else % on ne supprime que les nan 999999 ou 0
                CONC(cpt_c).data(index)=[];
                CONC(cpt_c).DPA(index)=[];
                cpt_c=cpt_c+1;% incrementation du compteur nb de CONC sauv�s
            %end
         end
    %else
    %    error('pb nombre entete <> data')
    %end
else
    error('pb de fichier');
end


