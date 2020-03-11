% Pré_RNA2PROT Programme principal: pré-traitement avant RNA2PROT
% lecture données dans fichiers texte sauvegarde dans un fichier matlab

% En entrée: 3 fichier de données poids.mat, nomfich1.txt ET nomfich2.txt
% poids.mat se trouve dans ../../../data/Isma/donnees_a_modeliser
% par défaut repertoire_brutes=../../../data/Isma/donnees_brutes/
% nomfich1.txt=PROT.txt 
% nomfich2.txt=RNA.txt
% (ou autre nom fourni par l'utilisateur)
% les 2 fichiers texte contenant les données en colonne
% nomfich1= concentration proteines en fmol/gFW
% nomfich2= concentration transcrits en fmol/gFW
% données en colonnes : 1ère ligne entête temps/DPA nom proteines(solyc)
% lignes suivantes: DPAi, conc. solyc1 à DPAi, solyc2 à DPA_i, ...


% appel: importfile(nomfichx)
% lecture du fichier , classement des solyc par ordre alpha, 
% suppression des données NaN ou 999 999 ou 0

% vérification adéquation nomfich1/nomfich2 (nom solyc similaire)
% affichage des nuages de points (normalisés par la moyenne)
% sauvegarde des données dans un fichier matlab dans le répertoire
% ../../../data/Isma/donnees_a_modeliser
% par défaut nomfich1.m et nomfich2.m prennent les noms d'entrée 
% nomfich1.txt et nomfich2.txt (.txt remplacé par .m)

% suppression/ajout affichage graphique
affichage=1; % mettre 1 si affichage souhaité

close all;
% Répertoires et nom fichiers par défaut
repertoire_brutes = '../../../data/Isma/donnees_brutes/';
repertoire_a_modeliser ='../../../data/Isma/donnees_a_modeliser/';
%nomfich1 = 'PROT_34var_filtrees';
%nomfich2 = 'RNA_34var_filtrees';
nomfich1 = '37P';
nomfich2 = '37R';

% Saisie répertoire et nom de fichier utilisateur (0=existe 1=n'existe pas )
[File1] = LectureNomFichier(repertoire_brutes,[nomfich1 '.txt'],0);
[File2] = LectureNomFichier(repertoire_brutes,[nomfich2 '.txt'],0);

% Lecture des données et nettoyage
PROT=importfile(File1);
RNA=importfile(File2);

% fichier des poids
load ([repertoire_a_modeliser 'POIDS.mat'])
 
% nombre de PROTEINES 
nb_prot=numel(PROT);
stop=0;idx=[];
% Vérification nombre proteines=nombre transcrits
if numel(RNA)~=nb_prot
    disp('pb fichier: nb RNA et nb PROT ne correspondent pas')
    stop=1;
else
    for i=1:nb_prot
        if (strcmp(PROT(i).NOM,RNA(i).NOM)~=1)
            disp('pb fichier: RNA/PROT ne correspondent pas')
            stop=1;
            %break
        else % si le nombre de données est <20
            if (numel(PROT(i).data)<20|numel(RNA(i).data)<20)
                %détection des couples
                idx=[idx i]
            end;
            
        end
    end
    % suppression des couples
    disp('suppression')
    idx
    PROT(idx).NOM
    PROT(idx)=[];RNA(idx)=[];
end
nb_prot=numel(PROT);
if stop==0
    if affichage==1
        for i=1:nb_prot
            %i=1;
            p0=mean(PROT(i).data);
            figure(i);
            hAx=scatter(PROT(i).DPA,PROT(i).data/p0,'r+');
            q0=mean(RNA(i).data);
            hold on;hAx2=scatter(RNA(i).DPA,RNA(i).data/q0,'b*');
            title(PROT(i).NOM,...
                'FontWeight','bold')
            xlabel('DPA ')
            ylabel('Conc P/R fmol/gFW normalisée par moyenne') % y-axis
            
            legend([hAx,hAx2],{'P','R'},'FontSize',5,'Location','northeast');
            figure(i+nb_prot);
            if length(PROT(i).DPA)==length(RNA(i).DPA)
                scatter(PROT(i).data/p0,RNA(i).data/q0,'r+');title(PROT(i).NOM,...
                    'FontWeight','bold')
                xlabel('P fmol/gFW normalisée')
                ylabel('R fmol/gFW normalisée') % y-axis
            end         
        end;      
    end;
     % sauvegarde dans un fichier matlab (qui ne doit pas déjà exister)
     %(0=existe 1=n'existe pas )
        [File3] = LectureNomFichier(repertoire_a_modeliser,[nomfich1 '.mat'],1);
        [File4] = LectureNomFichier(repertoire_a_modeliser,[nomfich2 '.mat'],1);
        save(File3,'PROT');
        save(File4,'RNA');
end;
