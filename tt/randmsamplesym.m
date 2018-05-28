function B = randmsamplesym(A,ratio)
%RANDMSAMPLESYM Summary of this function goes here
%   Detailed explanation goes here
m = size(A,1);
n = size(A,2);
if m~=n
    error {'A must be symmetry matrix. '};
end
B=zeros(n,n);
N=floor((n^2+n)/2);
b = zeros(N,1);
d = round(N*ratio);
vv = randperm(N);
vv=vv(1:d);
b(vv)=ones(d,1);
k=1;
for ii=1:n
    B(ii,ii)=b(ii)/2;
end
for jj=1:n-1
    for ii=jj+1:n
      B(ii,jj)=b(k+n);
      k=k+1;
    end
end
B=B+B';
B = B.*A;

end

