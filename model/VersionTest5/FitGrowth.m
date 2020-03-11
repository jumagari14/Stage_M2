% Modèles simples pour fitter la croissance des fruits
% Lecture des données: nom, (DPA_i,r_i,p_i) i=1..n

%provisoire
%global ptrans3 ptrans6 St6 mt6 St3 mt3 plogtrans3 Slogt3 mlogt3 degre h fitR
global xi yi
%global DPA_R DPA_P R
global fitPoids parv parg parc parl pare pard Wp3 Sp3 mp3

% repertoire courant;
chemin=pwd();
% ajout au path du répertoire stixbox
addpath([chemin,'\stixbox']);

% Lecture des données dans fichier(s)
% DPA matrice colonne
% data matrice colonne

close all;
repertoire_a_modeliser ='../../data/Isma/donnees_a_modeliser/';
% nom des fichiers par défaut
nomfich= 'POIDS';
% Saisie nom répertoire et fichier des données (0=existence-mode lecture)
[FileP] = LectureNomFichier(repertoire_a_modeliser,[nomfich '.mat'],0);

% lecture des fichiers POIDS
load(FileP);

xi=poids(:,1);yi=poids(:,2);

% choix de la courbe de croissance (fitPoids variable globale)
[errPoids,messagePoids]=fit_testpoids(xi,yi);
cpt_fen=4;%compteur fenetre

% sauvegarde figures poids et fermeture fenetres

    for i=1:cpt_fen
        nomfich=['graphes/Poids'  num2str(i),'.pdf']
        saveas(figure(i),nomfich)
    end
    
% sauvegarde fichiers
        
        RESULT.Poids{1,1}={fitPoids;messagePoids;errPoids;'erreur relative fit'};
        RESULT.Poids{1,2}=[{'verhulst' parv 'gompertz' parg 'contois' parc 'logistic' parl 'empirique' pare 'double sigmoide' pard 'poly' Wp3 Sp3 mp3 'noyau' h}];
        RESULT.Poids{1,3}=[xi,yi];
                
        save('graphes/RESULT.mat');
