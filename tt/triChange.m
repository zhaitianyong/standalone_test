function  [Dret,Dretsanmple ] = triChange( Dtri, d,n )
%TRICHANGE Summary of this function goes here
%   Detailed explanation goes here
Dret = Dtri;
for ii = 1:n
    for jj = 1:n
         if d(ii,jj)~=0
             Dret(ii,jj) = d(ii,jj);
         end
    end
end

Dretsanmple = Dret-d;





end

