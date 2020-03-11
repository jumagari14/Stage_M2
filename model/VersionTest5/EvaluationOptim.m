function [ score,par,resnorm,message] = EvaluationOptim( par,resnorm,flag)
% évaluation qualité du fit par attribution d'un score


% classement des résultats en fonction de resnom (du plus petit au plus
% grand)
% on renvoie les paramètres classés selon resnorm du plus petit au plus grand
[resnorm indice]=sort(resnorm);
par=par(:,indice);
parMu00=par(:,1);
parMu1=par(:,2);
parMu2=par(:,3);
parMu3=par(:,4);
score=0;message='';

% l'optimisation s'est correctement dérouler 
if flag>0
    % On garde les 3 meilleures optimisations
    if  (max(parMu00,parMu1))>1e-4
        % ecart relatif
        disp('ecart rel1')
        EcartP01=abs(parMu00-parMu1)./max(parMu00,parMu1)
    else % l'un des max est petit
        % on regarde les ecarts absolus
        EcartP01=abs(parMu00-parMu1)
        disp('ecart 1')
    end
    if (max(parMu1,parMu2))>1e-4
        disp('ecart rel2')
        EcartP12=abs(parMu1-parMu2)./max(parMu1,parMu2)
    else
        disp('ecart 1')
        EcartP12=abs(parMu1-parMu2)
    end
    % si l'un des max des paramètres est nul ou très petit
    
    % si les 3 fits sont à peu près égaaux
    if (EcartP01<5e-2&EcartP12<5e-2)
        % % si l'optimiseur a convergé 3 fois vers la même valeur en partant
        % d'initialisations différentes
        message='Optimisation: converge vers 3 valeurs identiques'
        score=10;
    elseif  (EcartP12<1e-1&EcartP01<1e-1)%
        score=8;
        message='Optimisation: converge vers 3 valeurs identiques à peu près 3 fois'
    elseif EcartP01<5e-2 %2 fits égaux
        score=6;
        message='Optimisation: converge vers 2 valeurs identiques'
    elseif EcartP12<1e-1
        score=4;
        message='Optimisation: converge vers 2 valeurs identiques à peu près'
    else % optimisation a convergé vers des valeurs différentes
        score=1;
        message='Optimisation: convergence vers différentes valeurs'
    end;
else
    message='Optimisation: aucune convergence'
    score=0;
end;

end

