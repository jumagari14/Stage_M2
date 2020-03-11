function [ val ] = mu( t)
% fonction t->mu(t) taux de croissance

global fitPoids parv parg parc parl pare pard Wp3 Sp3 mp3 h xi yi

switch fitPoids
    case 1 % croissance logistique: modèle de Verhulst
        verhulst=@(par,t)(par(2)*par(3)./(par(3)+(par(2)-par(3))*exp(-par(1).*t)));
        y=verhulst(parv,t);
        mu=@(t,parv)(parv(1).*(1-y./parv(2)));
        val=mu(t,parv);
    case 2 % modèle Gompertz
        gompertz=@(par,t)(par(2)*exp(log(par(3)/par(2)).*exp(-par(1).*t)));
        y=gompertz(parg,t);
        mu=@(t,parg)(parg(1).*log(parg(2)./y));
        val=mu(t,parg);
    case 3 % modèle Contois
        yc=resol_contois(parc,t);
        val=parc(1).*(1-yc./parc(2))./(parc(2)+(parc(3)-1).*yc);
    case 4 % log gen
    case 5 % empirique log(y(t))=V(t-a)/(k+t-a)
        val=pare(1)*pare(2)./(pare(2)+t-pare(3)).^2;
    case 6 % noyau
        Wfit=ksr2(t,xi,yi,h);dWfit=dksr2(t,xi,yi,h);
        val=dWfit./Wfit;
    case 7 % log poly
        % valeurs du poly à yi (temps experimentaux)
        %[ylp,dW3] = polyval(Wp3,t,Sp3,mp3);
        [dWp3]=polyder(Wp3);
        [dWval6]=polyval(dWp3,t,Sp3,mp3)./mp3(2);
        val=dWval6;
    case 8 %double sigmoide       
         val=pard(2).*exp(-pard(2).*(t-pard(3)) )./( 1+exp(-pard(2).*(t-pard(3))) )...
           +pard(6).*exp( -pard(6).*(t-pard(7)) )./(1+exp(-pard(6).*(t-pard(7)) ));
end





