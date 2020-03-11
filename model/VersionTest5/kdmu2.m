function [ mu ] = kdmu2( x)
% fonction x->mu(x)=taux de croissance du fruit
% soit par un polynome soit par estimateur à noyau
global xi yi h
global Wp6 SW6 mW6 Wp3 SW3 mW3  fitW


if fitW==6
    [Wval6, dW6]=polyval(Wp6,x,SW6,mW6);
    [dWp6]=polyder(Wp6);
    [dWval6]=polyval(dWp6,x);
    mu=dWval6./Wval6;
elseif fitW==3
    [Wval3, dW3]== polyval(Wp3,x,SW3,mW3);
    [dWp3]=polyder(Wp3);
    [dWval3]=polyval(dWp3,x);
    mu=dWval3./Wval3;
else % estimateur à noyau
    if fitW~=0 % h=6 (sinon h=8)
        h=6;
    end;
    
    % KSR   Kernel smoothing regression
    
    % Gaussian kernel function
    kerf=@(z)exp(-z.*z/2)/sqrt(2*pi);
    dkerf=@(z)-exp(-z.*z/2).*z./sqrt(2*pi);
    
    %r.f=zeros(1,N);
    N=length(xi);Nx=length(x);
    z=zeros(N,1);dz=zeros(N,1);
    for k=1:Nx
        z=kerf((x(k)-xi)/h);
        dz=dkerf((x(k)-xi)/h)/h;
        rf(k)=sum(z.*yi)/sum(z);
        drf(k)=(sum(dz.*yi)*sum(z)-sum(z.*yi)*sum(dz))/(sum(z))^2;
        mu(k)=drf(k)/rf(k);
    end
    mu=mu';
end







