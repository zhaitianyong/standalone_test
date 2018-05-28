% This function performs procrustes
%
% Babak Alipanahi
% University of Waterloo
% 
% June 16, 2011

function [R t err] = procrustesb(X,Y,reflection)

if ~exist('reflection','var')
    reflection = 1;
end

[dx nx] = size(X);
[dy ny] = size(Y);

if dx ~= dy || nx ~= ny
    error('X and Y must be of the same size')
end

mX = sum(X,2)/nx;
mY = sum(Y,2)/ny;

cX = bsxfun(@minus,X,mX);
cY = bsxfun(@minus,Y,mY);

normX = 1; sqrt(sum(cX(:).^2));
normY = 1; sqrt(sum(cY(:).^2));

cX = cX/normX;
cY = cY/normY;

M = cY*cX';
[U S V] = svd(M);


% outputs
R = V*U';

if ~reflection && det(R) < 0
    V(:,end) = -V(:,end);
    R = V*U';
end


t = mX - R*mY;
err = norm(cX - R*cY,'fro');