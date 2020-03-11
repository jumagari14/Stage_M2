function rf=ksr2(t,x,y,h)
% KSR2   Kernel smoothing regression
%
% (x,y) points expérimentaux
% rf(t)=estimateur à noyau de gauss à partir des points (x,y)à l'instant t
% rf(t) = sum(kerf((t-x)/h).*y)/sum(kerf((t-x)/h))

% Gaussian kernel function
kerf=@(z)exp(-z.*z/2)/sqrt(2*pi);

N=length(x);Nx=length(t);
z=zeros(N,1);rf=zeros(Nx,1);
for k=1:Nx
    z=kerf((t(k)-x)/h);
    rf(k)=sum(z.*y)/sum(z);
end