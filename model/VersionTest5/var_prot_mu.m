function ypoint = var_prot_mu(t,y,param)
% calcul de dy/dt=param(2)*mRNA(t)-(param(3)+mu(t))*y
% dy/dt=ks*mRNA(t)-(ks+mu)*y
% appel function mRNA_deg
global xi yi h

ypoint = param(2)*mRNA_deg(t)-(param(3)+mu(t))*y;
 
 