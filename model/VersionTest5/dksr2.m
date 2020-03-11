function drf=dksr2(t,x,y,h)
% derivative KSR   Kernel smoothing regression
% (x,y) points expérimentaux
% dr=dksr(t,x,y,h) returns the derivative of Gaussian kernel regression such that
% rf(t) = sum(kerf((t-x)/h).*y)/sum(kerf((t-x)/h))
%
% Gaussian kernel function
kerf=@(z)exp(-z.*z/2)/sqrt(2*pi);
dkerf=@(z)-exp(-z.*z/2).*z./sqrt(2*pi);

N=length(x);Nx=length(t);
z=zeros(N,1);
drf=zeros(N,1);
for k=1:Nx
    z=kerf((t(k)-x)/h);
    dz=dkerf((t(k)-x)/h)/h;
    drf(k)=(sum(dz.*y)*sum(z)-sum(z.*y)*sum(dz))/(sum(z))^2;
end