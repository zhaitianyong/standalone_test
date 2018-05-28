% This function determines the coordinates of
% atoms in the peptide plane
%
% Babak Alipanahi
% University of Waterloo
% May 10, 2010

function outX = PlaneFinder(X,cur,next)

% atoms: CA O C N HN CA 
% (last three are from the next residue)
% X contains coordinates of CA, C, and N
[bond_len bond_ang] = parsebonds(cur,next);
% bond lens
% 1:C-O, 2:C-N, 3:C-CA, 4:CA-N, 5:N-H 
% bond angs
% 1: CA-C-O, 2: N-C-O, 3: CA-C-N
% 4: C-N-CA, 5: C-N-H, 6: H-N-CA

if cur == 15
    bond_len(5) = 1.473;
    bond_ang(5) = 125;
end

% C is at origin
tX = zeros(3,6);
% N is above C
tX(2,4) = bond_len(2);
% O 
tX(1,2) = bond_len(1)*cosd(bond_ang(2) - 90);
tX(2,2) = -bond_len(1)*sind(bond_ang(2) - 90);
% CA (i)
tX(1,1) = -bond_len(3)*cosd(bond_ang(3) - 90);
tX(2,1) = -bond_len(3)*sind(bond_ang(3) - 90);
% CA (i+1)
tX(1,6) = bond_len(4)*cosd(bond_ang(4) - 90);
tX(2,6) = tX(2,4) + bond_len(4)*sind(bond_ang(4) - 90);
% HN / CD for PRO
tX(1,5) = -bond_len(5)*cosd(bond_ang(5) - 90);
tX(2,5) = tX(2,4) + bond_len(5)*sind(bond_ang(5) - 90);


% CA, C, and N are given
talign = tX(:,[1 3 4]);      

[R t err] = procrustesb(X,talign);
if err > 1e-2
    disp('Warning: error is too large.');
end

outX = bsxfun(@plus,R*tX,t);

% [err, ~, R] = procrustes(X',talign','scaling',false,'reflection',false);
% if err > 1e-2
%     disp('Warning: error is too large.');
% end
% tvec = repmat(R.c(1,:),6,1);
% outX = (tX'*R.T + tvec)';        

