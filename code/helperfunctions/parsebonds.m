% bond lens
% 1:C-O, 2:C-N, 3:C-CA, 4:CA-N, 5:N-H 
% bond angs
% 1: CA-C-O, 2: N-C-O, 3: CA-C-N
% 4: C-N-CA, 5: C-N-H, 6: H-N-CA
% 7: C-CA-N

function [bond_len bond_ang] = parsebonds(cur,next)

bond_len = zeros(5,1);

% -Lengths
% 
% C-O bond
% C-O      : 1.231

bond_len(1) = 1.231;

% C-N bond
% C-N      : 1.341 (Pro)
% C-NH1    : 1.329 (N with one H)

if next == 15
    bond_len(2) = 1.341;
else
    bond_len(2) = 1.329;
end

% C-CA bond
% C-CH2G   : 1.516 (Gly)
% C-CH1E   : 1.525 (All but Gly)

if cur == 8
    bond_len(3) = 1.516;
else
    bond_len(3) = 1.525;
end

% CA-N bond    
% CH2G-NH1 : 1.451 (Gly)    
% CH1E-N   : 1.466 (Pro)
% CH1E-NH1 : 1.458 (All but Pro and Gly)

if cur == 8 
    bond_len(4) = 1.451;
elseif cur == 15
    bond_len(4) = 1.466;
else
    bond_len(4) = 1.458;
end
% N-H bond
bond_len(5) = 0.98;

bond_ang = zeros(7,1);

%=========================================
% -Angles
% CA-C-O
% CH1E-C-O    : 120.8
% CH2G-C-O    : 120.8

if cur == 8
    bond_ang(1) = 120.8;
else
    bond_ang(1) = 120.8;
end

% N-C-O
% NH-C-O  : 123
% N-C-O   : 122

if next == 15
    bond_ang(2) = 122;
else
    bond_ang(2) = 123;
end

% CA-C-N
% CH1E-C-N    : 116.9
% CH2G-C-N    : 118.2
% CH1E-C-NH1  : 116.2
% CH2G-C-NH1  : 116.4

bond_ang(3) = 360 - sum(bond_ang(1:2));

% 
% C-N-CA
% C-N-CH1E    : 122.6 (Pro)
% C-NH1-CH1E  : 121.7
% C-NH1-CH2G  : 120.6

if next == 8
    bond_ang(4) = 120.6;
elseif cur == 15
    bond_ang(4) = 122.6;
else
    bond_ang(4) = 121.7;
end

% C-NH1-H
bond_ang(5) = 120;

bond_ang(6) = 360 - sum(bond_ang(4:5));

% C-CA-N
% C-CH1E-N  : 111.8
% C-CH1E-NH1: 111.2
if cur == 15
    bond_ang(7) = 111.8;
else
    bond_ang(7) = 111.2;
end