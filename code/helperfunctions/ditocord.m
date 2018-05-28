% This function computes the fourth atom given
% the distances between atoms and phi or psi angle
%
%
% Babak Alipanahi
% University of Waterloo
% May 10, 2010

function point = ditocord(X,res,nextres,ang,type)

% bond lens
% 1:C-O, 2:C-N, 3:C-CA, 4:CA-N, 5:N-H
% bond angs
% 1: CA-C-O, 2: N-C-O, 3: CA-C-N
% 4: C-N-CA, 5: C-N-H, 6: H-N-CA

[bond_len bond_ang] = parsebonds(res,nextres);

v12 = X(:,2) - X(:,1); v12 = v12/norm(v12);
v23 = X(:,3) - X(:,2); v23 = v23/norm(v23);


ang_12_23 = angvec(v12,v23);

na = cross(v12,v23);
na = na/norm(na);

switch type
    case 'PHI'
        % CN--CAC (i-1,i)
        % C-CA bond length
        d34 = bond_len(3);
        % C-CA-N angle
        ang_23_34 = pi - pi*bond_ang(7)/180;

    case 'PSI'
        % NCA--CN (i,i+1)
        % C-N bond length
        d34 = bond_len(2);
        % CA-C-N angle
        ang_23_34 = pi - pi*bond_ang(3)/180;
       
end

p_na    = d34*sin(ang_23_34)*cos(ang-pi/2);
p_v12   = d34*sin(ang_23_34)*sin(ang-pi/2)/sin(ang_12_23);
p_v23   = d34*cos(ang_23_34)- d34*sin(ang_23_34)*sin(ang-pi/2)/tan(ang_12_23);

v34 = p_na*na + p_v12*v12 + p_v23*v23;
point = X(:,3) + v34;
% nX = [X(:,1:3) point];



function ang = angvec(a,b)

ang = acos(dot(a,b)/norm(a)/norm(b));

% function dist = distplane(n,X0,point)
% 
% d = n'*X0;
% dist = abs(n'*point + d);

