function obj = fonction_ecart(param,time,data)
% minimisation des écarts au carré
 F = resol_Reparam(param,time);
  obj=norm(F-data,2)^2;
