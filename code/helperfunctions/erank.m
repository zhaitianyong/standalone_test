% Computes the effective rank of a matrix
%
% Babak Alipanahi
% University of Waterloo
% Jan. 12, 2010


function out = erank(A)
sv = svd(A);
out = sum(sv/max(sv) > 1e-6);

