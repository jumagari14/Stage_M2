function [ indice ] = choix_data(nb_indice)
% choix de l'utilisateur des (RNA,PROT) � �tudier
% ie recherche de leur (ks,kd)
% au choix: toutes donn�es, un tirage al�atoire, debut-fin, un choix
% d'indices

indice=-1;
while indice==-1
    prompt = '1-toutes donn�es 2-tirage al�atoire 3-debut-fin 4-saisie d indices ';
    dlg_title = 'Choix des donn�es';
    num_lines = 1;
    defaultans = {'1'};
    choix = inputdlg(prompt,dlg_title,num_lines,defaultans);
    
    switch choix{1}
        case '1'
            indice=1:nb_indice;
        case '2'
            nb_tir=0;
            % tirage al�atoire
            while nb_tir<1|nb_tir>nb_indice
                prompt = 'nbre de tirages: ';
                dlg_title = 'Tirage al�atoire ';
                num_lines = 1;
                defaultans = {num2str(nb_indice)};
                choix = inputdlg(prompt,dlg_title,num_lines,defaultans);
                nb_tir=floor(str2num(choix{1}));
            end;
            V = 1:nb_indice;
            nV=numel(V);
            for n=1:nb_tir
                idx = randperm(nV);
                indice(n) = V(idx(1));
            end
            indice=unique(sort(indice));
        case '3'
            % saisie d'indices
            indice_debut=0;indice_fin=1
            while indice_debut<1|indice_fin>nb_indice|indice_debut>indice_fin
                prompt = {'d�but indice','fin indice:'};
                dlg_title = 'Input';
                num_lines = 1;
                defaultans = {'1',num2str(nb_indice)};
                answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
                indice_debut=str2num(answer{1});indice_fin=str2num(answer{2})
            end
            indice=indice_debut:indice_fin;
        case '4'
            indice_debut=0;indice_fin=1;
            while indice_debut<1|indice_fin>nb_indice|indice_debut>indice_fin
                x = inputdlg('entrez des indices s�par�s par un espace:',...
                    'Sample', [1 50]);
                indice = str2num(x{:});
                indice=unique(sort(floor(indice)));
                indice_debut=indice(1);indice_fin=indice(end);
            end;
    end;
end;



