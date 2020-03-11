function ypoint = var_prot_all(t, y,param)
 % dy/dt=param(2)*mRNA(t)-param(3)*y
 % appel function mRNA_deg 
 %ypoint = ksi*mRNA(t)...-kdi*y;
     ypoint = param(2)*mRNA_deg(t)-param(3)*y;


 
 