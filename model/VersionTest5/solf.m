function [kparf,sol,errf]=solf(NOM,DPA,Pf,ks_nf);
% d�termination de (ks,kd) � partir des donn�es FRUIT
% normalis�es: (Rf,Pf) dp/dt=ks*r-kd*p
% 
% starting value condition initiale EDO: moyenne des valeurs initiales mesurees
[P0,Pmin,Pmax]=initiale(DPA,Pf);

% on cherche ks dans [ks_nf;5*ks_nf]
% starting value de ks=milieu=3*ks_nf
% starting value parameters
parI=[P0;ks_nf;ks_nf*0.1];

% bounds for parameters search
lb=[Pmin;ks_nf;0];
ub=[Pmax;ks_nf*5;inf];

sol=zeros(length(DPA),3);

% avec solution num�rique ED sans bornes (sauf positivit�) sur les param�tres
lb=[0;0;0];ub=[inf;inf;inf];
[parK3,resnorm3,residual3,exitflag3,output3]  = lsqcurvefit(@resol_all,parI,DPA,Pf,lb,ub);
sol(:,3)=resol_all(parK3,DPA);
% avec solution num�rique ED
[parK1,resnorm1,residual1,exitflag1,output1]  = lsqcurvefit(@resol_all,parI,DPA,Pf,lb,ub);
sol(:,1)=resol_all(parK1,DPA);
% avec solution explicite ED si R=poly deg3
[parK2,resnorm2,residual2,exitflag2,output2]  = lsqcurvefit(@solexplicite,parI,DPA,Pf,lb,ub);
sol(:,2)=resol_all(parK2,DPA);

kparf=[parK1;parK2;parK3];
 errf=[sqrt(resnorm1);sqrt(resnorm2); sqrt(resnorm3)]./norm(Pf,2);

 
end
