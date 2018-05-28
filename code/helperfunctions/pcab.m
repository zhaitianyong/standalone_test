% This function performs PCA
%
% Babak Alipanahi
% University of Waterloo
% July 29, 2010

function [lX sv V] = pcab(X)


[dim num] = size(X);
mX = mean(X,2);
cX = bsxfun(@minus,X,mX);

if dim < num   
    C = cov(cX',1);
    [V S] = eigb(0.5*(C+C'));
    
    sv = diag(S);
    nz = sv/max(sv) > 10^-10;
    sv = sv(nz);
    
    V = V(:,nz);
    lX = V'*cX;    
else
    K = cX'*cX;
    [V S] = eigb(0.5*(K+K'));
    sv = diag(S);
    nz = sv/max(sv) > 10^-10;
    sv = sv(nz);
    
    V = V(:,nz);
    dS = diag(S(nz,nz));
    dS = sqrt(dS);
    lX = bsxfun(@times,V',dS);
end

