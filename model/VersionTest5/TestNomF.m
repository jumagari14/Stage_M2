function [Test] = TestNomF(varargin)
% test l'existence d'un fichier varagin{1}
% Renvoi -1 aucun argument, 0 existe, 1 n'existe pas dans le répertoire, 2
% le répertoire n'existe pas

if ~nargin 
    %Pas d'argument
    %disp('aucun nom de fichier');
    Test=-1;
else
    File1=varargin{1};
    if exist(File1,'file')
        % le fichier existe
        Test=0;
    else
        %le fichier n'existe pas => on va chercher à savoir pourquoi
        [rep,nom,ext] = fileparts(File1);
        if ~exist(rep,'dir')
            Test=2;%le repertoire n'existe pas
        else
            Test=1;%le fichier n'existe pas dans ce répertoire
         end
        
    end
end
        