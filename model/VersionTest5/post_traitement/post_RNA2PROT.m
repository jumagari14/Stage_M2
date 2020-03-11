% post traitement de RNA2PROT
% lecture des fichiers dans le répertoire resultats

% sauvegarde de tous les résultats dans un fichier 'table_compilation.mat'
clear all
% détection des fichiers .mat contenu dans le repertoire resultats
rep = '../resultats/';
ext = '*.mat';
chemin = fullfile(rep,ext);
list = dir(chemin);

%sauvegarde files.mat in folder 'resultats'
%  RESULT.NOM=NOM{1};%         
%         RESULT.Poids{1,1}={fitPoids;messagePoids;errPoids;'erreur relative fit'}; 
%         RESULT.Poids{1,2}=[{'verhulst' parv 'gompertz' parg 'contois' parc 'logistic' parl 'empirique' pare 'double sigmoide' pard 'poly' Wp3 Sp3 mp3 'noyau' h}];
%         RESULT.Poids{1,3}=[xi,yi];
%         
%         RESULT.Rna{1,1}={fitR;messageR;errR;'erreur relative fit'}; 
%         RESULT.Rna{1,2}=[{'poly deg3' ptrans3 St3 mt3 'poly deg 6' ptrans6 St6 mt6  'log poly deg3' plogtrans3 Slogt3 mlogt3 degre h}];
%         RESULT.Rna{1,3}=[DPA_R,Rg];
%         RESULT.Rna{1,4}=[{'moyenne',R0g}];

%         RESULT.Prot{1,1}=[{'P0' kparg(1,1) 'ks' kparg(2,1)*P0g/R0g 'Ks' kparg(2,1) 'kd' kparg(3,1)}]
%         RESULT.Prot{1,2}=[{{'erreur relative fit',errg(1)}; {messageOF(1) scoreOF(1);messageOF(2) scoreOF(2)}}];
%         RESULT.Prot{1,3}=[DPA_P,Pg];
%         RESULT.Prot{1,4}=[{'moyenne',P0g}];
  
%fitPoids=8; errPoids=0; fitR=2; errR=0; score=0

% pour chacun des fichiers
ResultatMatrix=[];
for n=1 : numel(list)
      nom_n = list(n).name;
 % nom=[nom;nom_n]
     load(fullfile(rep,list(n).name))
   
     % stockage de tous les résultats dans une table/matrice
 
     ligne_n=[RESULT.NOM, RESULT.Poids{1,1}{2,1}, RESULT.Poids{1,1}{3,1},...
         RESULT.Rna{1,1}{2,1}, RESULT.Rna{1,1}{3,1}, RESULT.Rna{1,4}{1,2},...
         RESULT.Prot{1,4}{1,2}, RESULT.Prot{1,1}{1,2}, RESULT.Prot{1,1}{1,4}, RESULT.Prot{1,1}{1,6}, RESULT.Prot{1,1}{1,8},...
        RESULT.Prot{1,2}{1,1}{1,2}, RESULT.Prot{1,2}{2,1}{1,1}, RESULT.Prot{1,2}{2,1}{1,2}, RESULT.Prot{1,2}{2,1}{2,1}, RESULT.Prot{1,2}{2,1}{2,2}];
     % Solyc fitPoids(8=messagePoids(doublesigmoid, errPoids) fitR(2=messageR(exp polydeg3) errR moyenne=R0g...
     ...moyenne=P0g P0=kparg(1,1) ks=kparg(2,1)*P0g/R0g Ks = kparg(2,1) kd = kparg(3,1)...
         ...errg(1)=erreur relative fit  messageOF(1) scoreOF(1) messageOF(2) scoreOF(2)
         
     ResultatMatrix=[ResultatMatrix; ligne_n];   
end

line1={'Solyc', 'fitPoids', 'errFitPoids', 'fitR', 'errFitRNA', 'R0g','P0g','P0', 'ks', 'Ks', 'kd', 'errg', 'messageOF', 'scoreOF', 'messageOF', 'scoreOF'};
ResultatMatrix=[line1;ResultatMatrix]; 
    
% sauvegarde des résultats dans une matrice
% [File1] = LectureNomFichier(rep, ['nomfich' '.mat'], 1);
% je n'arrive pas à sauvegarder le fichier avec le nom saisi 'nomfich' ??
 save ../resultats/table_compilation ResultatMatrix
 