% Kdagger function K = kdagger(D)
%
% Babak Alipanahi
% University of Waterloo
% Jan 12, 2010


function K = kdagger(D)

num = size(D,1);
J = eye(num) - (1/num)*ones(num,1)*ones(num,1)';
K = -0.5*J*(D - diag(diag(D)))*J;
