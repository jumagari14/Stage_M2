function rho=matrice_correlation(X,nobs,ecart)
% matrice de corrélation entre les paramètres
% Correlation= ecart/n-p * (tXX)^[-1)

U=transpose(X)*X;
% matrice de covariance
cov=ecart/(nobs-3)*inv(U);
% matrice de coorelation
for i=1:3
    for j=1:3
        rho(i,j)=cov(i,j)/sqrt(cov(i,i)*cov(j,j));
    end
end

