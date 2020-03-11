function ypoint = var_prot_mu_reparam2(t,y,param)
% EDO repraramétrée dy/dt=param(2)*(mRNA(t)-param(3)/param(2)*y)-mu(t)*y
% 
% appel function mRNA_deg
global xi yi h

ks=param(2);
K=param(3)/param(2);
ypoint = ks*(mRNA_deg(t)-K*y)-mu(t)*y;