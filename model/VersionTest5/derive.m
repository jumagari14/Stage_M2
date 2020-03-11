function spoint = derive(t,s,param)
% dérivée de p par rapport à t et aux paramètres
%
% dérivée de P par rapport à t, à param(1)=P0, à param(2)=ks et param(3)=kd
% dP/dt=Ks RNA(t)-(kd+mu)P(t)

% s1=P ds1/dt=ks*mRNA(t)-(kd+mu)s1
ds1dt=param(2)*mRNA_deg(t)-(param(3)+mu(t))*s(1);
% en posant s2=dP/dy0=-(kd+mu)s2
ds2dt=-(param(3)+mu(t))*s(2);
% en posant s3=dP/dks=-(kd+mu)s3+mRNA
ds3dt=-(param(3)+mu(t))*s(3)+mRNA_deg(t);
% en posant s4=dP/dkd=-(kd+mu)*s4 -s1
ds4dt=-(param(3)+mu(t))*s(4)-s(1);

spoint=[ds1dt;ds2dt;ds3dt;ds4dt];