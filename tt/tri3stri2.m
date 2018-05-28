function [m,num_seg] = tri3stri( DH0 )
%TRISTRI Summary of this function goes here
%   Detailed explanation goes here
m = DH0;
n = size(DH0,1)
num_seg=zeros(n);
% n=length(m);
%tic

%for i=1:n
% if (mod(i,100)==0)
%   display(i)
% end
%    for j=1:n
%        if m(i,j)==0
%            m(i,j)=1000000;
%        end
%    end
%end
  m(find(m(:)==0))=1000000;
%m
m=m-1000000*eye(n);    
% %
% it depends on the size of the real protein;
inf=1000000;
% %
for i=1:n
%tic
 if (mod(i,100)==0)
 i
 end
    for j=1:n
        %if m(j,i)<inf
	 if m(n*(i-1)+j)<inf
            %if i~=j && num_seg(i,j)==0
            %    num_seg(i,j)=1;  
if i~=j && num_seg((j-1)*n+i)==0
  num_seg((j-1)*n+1)=1;
end
            % %
            for k=1:n
%                if m(i,k)<inf
		    if m((k-1)*n+i)<inf
                    % %
%                    if i~=k && num_seg(i,k)==0
		    if i~=k && num_seg((k-1)*n+i)==0
                      %  num_seg(i,k)=1;
num_seg((k-1)*n+i)=1;
                    end
                    % %
                 %    s=(m(j,i)+m(i,k));
s=(m((i-1)*n+j)+m((k-1)*n+i));

                    %if s<m(j,k)
		    if s<m((k-1)*n+j)
%                        m(j,k)=s;
m((k-1)*n+j)=s;
                        % %
           %             num_seg(j,k)=(num_seg(i,j)+num_seg(i,k));
num_seg((k-1)*n+j)=(num_seg((j-1)*n+i)+num_seg((k-1)*n+i));
                    end
                end
            end
        end
    end
   % toc
end
display('triangle inequality estimation finished!')
% figure;
% mesh(m);

end

