function [m,num_seg] = tri3stri( DH0 )
%TRISTRI Summary of this function goes here
%   Detailed explanation goes here
m = DH0;
n = size(DH0,1);
num_seg=zeros(n);
% n=length(m);
for i=1:n
    for j=1:n
        if m(i,j)==0
            m(i,j)=1000000; % ��û����ֵ�ĵط���1000000����
        end
    end
end
m=m-1000000*eye(n); % ȥ���Խ��ߣ����Խ�������
      
%%
% it depends on the size of the real protein;
inf=1000000;
%%
for i=1:n
    for j=1:n
        if m(j,i)<inf
            %% ��������if����и�num_seg(i,j)��ֵ1��˵����Ϊʲô��1��������������1
              %��ֵ:�����÷�֤������m(i,j)<inf�� i~=j��num_seg(i,j)==0��num_seg(i,j)=1
              %��num_seg(i,j)>1���Ƴ�num_seg(i,j)�ѱ���ֵ���ȼ���m(i,j)�Ѿ�������������Ƴ�num_seg(i,j)~=0�����������ì�ܡ�
              %�±�ͬ��
            if i~=j && num_seg(i,j)==0
                num_seg(i,j)=1; %  ����ֵ��λ�ó��˷ǶԽ���Ԫ�� ��1
            end
            %%
            for k=1:n
                if m(i,k)<inf
                    %%
                    if i~=k && num_seg(i,k)==0
                        num_seg(i,k)=1;
                    end
                    %%
                     s=(m(j,i)+m(i,k));

                    if s<m(j,k)
                        m(j,k)=s;
                        %%
                        num_seg(j,k)=(num_seg(i,j)+num_seg(i,k));%����õ�����ʼ�ն�Ӧ��С·����number of segments;1ʼ�ն�Ӧ�ⶨ�ľ��롣
  
%                         if num_seg(i,j)>=num_seg(i,k)
%                            num_seg(j,k)=num_seg(i,j)+1+(num_seg(i,k)-1)*rand;
%                         else
%                            num_seg(j,k)=num_seg(i,k)+1+(num_seg(i,j)-1)*rand; 
%                         end
                        %%
                    end
                end
            end
        end
    end
end
% figure;
% mesh(m);

end

