function spoint = derive_reparam2(t,s,param)
% reparammetrisation dérivée de p par rapport à t et aux nouveaux paramètres
ks=param(2);
K=param(3)/param(2);

% dérivée de P par rapport à t, à param(1)=P0, à K et kd
% dp/dt=ks(r-K*p)-mu*p

% s1=P ds1/dt=ks*(mRNA(t)-K*s1)-mu s1
ds1dt=ks*(mRNA_deg(t)-K*s(1))-mu(t)*s(1);
% en posant s2=dP/dy0=-(ks*K+mu)s2
ds2dt=-(ks*K+mu(t))*s(2);
% en posant s3=dP/dks=-(ks*K+mu)s3-ks*s1
ds3dt=-(ks*K+mu(t))*s(3)-ks*s(1);
% en posant s4=dP/dkd=-(ks*K+mu)*s4 -K*s1+mRNA
ds4dt=-(ks*K+mu(t))*s(4)-K*s(1)+mRNA_deg(t);

spoint=[ds1dt;ds2dt;ds3dt;ds4dt];