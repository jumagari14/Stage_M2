function [ score,par,resnorm,message] = EvaluationOptim( par,resnorm,flag)
% �valuation qualit� du fit par attribution d'un score


% classement des r�sultats en fonction de resnom (du plus petit au plus
% grand)
% on renvoie les param�tres class�s selon resnorm du plus petit au plus grand
[resnorm indice]=sort(resnorm);
par=par(:,indice);
parMu00=par(:,1);
parMu1=par(:,2);
parMu2=par(:,3);
parMu3=par(:,4);
score=0;message='';

% l'optimisation s'est correctement d�rouler 
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
    % si l'un des max des param�tres est nul ou tr�s petit
    
    % si les 3 fits sont � peu pr�s �gaaux
    if (EcartP01<5e-2&EcartP12<5e-2)
        % % si l'optimiseur a converg� 3 fois vers la m�me valeur en partant
        % d'initialisations diff�rentes
        message='Optimisation: converge vers 3 valeurs identiques'
        score=10;
    elseif  (EcartP12<1e-1&EcartP01<1e-1)%
        score=8;
        message='Optimisation: converge vers 3 valeurs identiques � peu pr�s 3 fois'
    elseif EcartP01<5e-2 %2 fits �gaux
        score=6;
        message='Optimisation: converge vers 2 valeurs identiques'
    elseif EcartP12<1e-1
        score=4;
        message='Optimisation: converge vers 2 valeurs identiques � peu pr�s'
    else % optimisation a converg� vers des valeurs diff�rentes
        score=1;
        message='Optimisation: convergence vers diff�rentes valeurs'
    end;
else
    message='Optimisation: aucune convergence'
    score=0;
end;

end

