% This function determines the coordinates of
% side chain atoms 
%
% Babak Alipanahi
% University of Waterloo
% June 18, 2010

function outX = SidechainFinder(X,index_align,A)


tX = A.X(:,A.first_meta:A.last_meta);
talign = tX(:,index_align);

% num_atoms = size(tX,2);

% [err, ~, R] = procrustes(X',talign','scaling',false,'reflection',false);
% if err > 1e-2
%     disp('Warning: error is too large.');
% end
% tvec = repmat(R.c(1,:),num_atoms,1);
% outX = (tX'*R.T + tvec)';        

reflection = 0;
[R t err]  = procrustesb(X,talign,reflection);

if err > 5e-1
    disp('Warning: error is too large.');
end

outX = bsxfun(@plus,R*tX,t);    
