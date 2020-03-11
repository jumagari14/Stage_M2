function [ val ] = Poids( t)
% fonction t->Poids(t) poids du fruit

global fitPoids parv parg parc parl pare pard Wp3 Sp3 mp3 h xi yi

switch fitPoids
    case 1 % croissance logistique: modèle de Verhulst
        verhulst=@(par,t)(par(2)*par(3)./(par(3)+(par(2)-par(3))*exp(-par(1).*t)));
        val=verhulst(parv,t);
    case 2 % modèle Gompertz
        gompertz=@(par,t)(par(2)*exp(log(par(3)/par(2)).*exp(-par(1).*t)));
        val=gompertz(parg,t);
    case 3 % modèle Contois
        val=resol_contois(parc,t);
        
    case 4 % log gen
        val=resol_logisgen(parl,t);
    case 5 % empirique log(y(t))=V(t-a)/(k+t-a)
        empirique=@(par,t)(exp(par(1).*(t-par(3))./(par(2)+t-par(3))));
        val=empirique(pare,t);
    case 6 % noyau
        val=ksr2(t,xi,yi,h);
    case 7 % log poly
        % valeurs du poly à t 
        [val,dW3] = polyval(Wp3,t,Sp3,mp3);
        val=exp(val);%(transfo log des donnees)
    case 8 % double sigmoide
        sigmoide=@(par,t)(par(4)+par(1)./(1+exp(-par(2).*(t-par(3)) ))...
    +par(5)./(1+exp(-par(6).*(t-par(7)) )));
        val=sigmoide(pard,t);
    otherwise
        val=[];
        disp('erreur');
              
end





