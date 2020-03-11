function [FileToRead] = LectureNomFichier(repertoire,nomfich,existence)
% lecture du nom du fichier et vérification de son existence ou non
% existence

Test=2;
while Test~=existence;
prompt = {'Répertoire','Nom du fichier'};
                dlg_title = 'Saisie Nom du Fichier';
                num_lines = 1;
                defaultans = {repertoire,nomfich};
                answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
                FileToRead=[answer{1} answer{2}];
                
                [Test] = TestNomF(FileToRead);
end
    