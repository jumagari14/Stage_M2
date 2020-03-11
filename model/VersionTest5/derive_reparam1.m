function spoint = derive_reparam1(t,s,param)
% reparammetrisation dérivée de p par rapport à t et aux nouveaux paramètres
kd=param(3);
K=param(2)/param(3);

% dérivée de P par rapport à t, à param(1)=P0, à K et kd
% dp/dt=ks(K*r-p)-mu*p

% s1=P ds1/dt=ks*(K*mRNA(t)-s1)-mu s1
ds1dt=kd*(K*mRNA_deg(t)-s(1))-mu(t)*s(1);
% en posant s2=dP/dy0=-(kd+mu)s2
ds2dt=-(kd+mu(t))*s(2);
% en posant s3=dP/dks=-(kd+mu)s3+kd*mRNA
ds3dt=-(kd+mu(t))*s(3)+kd*mRNA_deg(t);
% en posant s4=dP/dkd=-(kd+mu)*s4 -s1+K*mRNA
ds4dt=-(kd+mu(t))*s(4)-s(1)+K*mRNA_deg(t);

spoint=[ds1dt;ds2dt;ds3dt;ds4dt];