function [ score,message ] = ScoreErreur( erreur,score )
% modification d'un score selon erreur relative 
% entre points exp�rimentaux et valeurs th�oriques
% 
if erreur<0.10
    message='fit excellent';
elseif erreur<0.15
    message='tr�s bon fit';
elseif erreur<0.20
    message='bon fit';
elseif erreur <0.30
    message='assez bon fit';
elseif erreur<0.4
    message='fit m�diocre';
else
    message='fit mauvais';
end
score=floor((0.5-erreur)*20);
end