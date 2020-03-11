function [ score,message ] = ScoreErreur( erreur,score )
% modification d'un score selon erreur relative 
% entre points expérimentaux et valeurs théoriques
% 
if erreur<0.10
    message='fit excellent';
elseif erreur<0.15
    message='très bon fit';
elseif erreur<0.20
    message='bon fit';
elseif erreur <0.30
    message='assez bon fit';
elseif erreur<0.4
    message='fit médiocre';
else
    message='fit mauvais';
end
score=floor((0.5-erreur)*20);
end