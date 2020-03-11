function [ val ] = mRNA_deg( t)
% fonction t->mRNA(t)=transcrit(t) 
% selon la valeur de la variable globale fitR
% soit par un polynome, soit par exp(poly)
% soit par estimateur à noyau
% degré du polynome de regression des transcrits
global DPA_R R h fitR
global ptrans3 ptrans6 St6 mt6 St3 mt3 degre
global plogtrans3 Slogt3 mlogt3

switch fitR
    case 6% polynome degré 6
       [pval, deltat6]=polyval(ptrans6,t,St6,mt6);
    case 3% polynome degré 3
       [pval,deltat3] = polyval(ptrans3,t,St3,mt3);
    case 2%exp(poly deg 3)
        % valeur du polynome fittant le log(donnéees)
        [pval,deltat3] = polyval(plogtrans3,t,Slogt3,mlogt3);
        pval=exp(pval);
    otherwise % estimateur à noyau
        if fitR~=0 % h=6 (sinon h=8)
            h=6;
        end;
        kerf=@(z)exp(-z.*z/2)/sqrt(2*pi);
        N=length(DPA_R);Nx=length(t);
        z=zeros(N,1);
        for k=1:Nx
            z=kerf((t(k)-DPA_R)/h);
            pval(k)=sum(z.*R)/sum(z);
        end
        pval=pval';
end
val=pval;




