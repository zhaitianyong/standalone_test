% This function computes dihedral angle based
% on the four points given
% 
% Babak Alipanahi
% University of Waterloo
% May 10, 2010

function angle = dicalc(X)

v1 = X(:,2) - X(:,1);
v2 = X(:,3) - X(:,2);
v3 = X(:,4) - X(:,3);

arg1 = dot(norm(v2)*v1,cross(v2,v3));
arg2 = dot(cross(v1,v2),cross(v2,v3));

angle = atan2(arg1,arg2);

u1 = cross(v1,v2);
u2 = cross(v2,v3);

ang = acos(dot(u1,u2)/norm(u1)/norm(u2));

if (abs(angle) - ang) > 1e-3
    error('dihedral error');
end

