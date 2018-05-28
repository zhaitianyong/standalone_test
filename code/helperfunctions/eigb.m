% Modified eigedecomposition function
% (gives sorted eigenvalues)
%
% Babak Alipanahi
% University of Waterloo
% June 25, 2010

function [V D] = eigb(A)

if nargout < 2
   V = eig(A);
   V = sort(V,'descend');
else
   [V D] = eig(A);
   dD = diag(D);
   [dD index_sort] = sort(dD,'descend');
   V = V(:,index_sort);
   D = diag(dD);
end