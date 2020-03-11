function [ c,ceq,dc,dceq ] = contraintes( par0, param )
%contraintes sur les param�tres (kd,K)
% parameters=[enz0,kd,K] avec K=ks/kd*q0/p0*1e5
% KS=kd*K=parameters(2)*parameters(3) (o� KS=ks*q0/p0*1e5)
% KS doit �tre proche des valeurs calcul�es par la
% formule de V. Mengin 
% KS in[80%*KS_calcul�, 120%*KS_calcul�]
% donc contraintes G(kd,K)=(g1(kd,K),g2(kd,K))
% avec g1(kd,K)=80%*KS_calcul�-KS
% g2(kd,K)=KS-120%*KS_calcul�
ceq=[];
dceq=[];
% contraintes G<0
c(1)=(par0(2)*par0(3))*0.8-param(2)*param(3);
c(2)=param(2)*param(3)-(par0(2)*par0(3))*1.2;
% gradient de G
dc=[-param(3),-param(2);param(3),param(2)];
end

