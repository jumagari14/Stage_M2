function obj = fonction_ecart(param,time,data)
% minimisation des �carts au carr�
 F = resol_Reparam(param,time);
  obj=norm(F-data,2)^2;
